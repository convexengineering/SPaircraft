" gp model of GA aircraft "
import numpy as np
from gpkit import Model, Variable, Vectorize
from gpkit.tools import te_exp_minus1
from gpkit.constraints.tight import Tight as TCS

# pylint: disable=attribute-defined-outside-init

class Aircraft(Model):
    "Aircraft class"
    def setup(self, wing):

        self.wing = wing

        Wpay = Variable('W_{payload}', 'N', 'Aircraft Payload Weight')
        Wstruct = Variable("W_{struct}", "N", "structural weight")
        fstruct = Variable("f_{struct}", "-", "fractional structural weight")
        Wdry = Variable("W_{dry}", "N", "aircraft dry weight")

        self.components = [self.wing]

        constraints = [Wdry >= Wstruct + Wpay + self.wing["W"]]

        return constraints, self.components

    def dynamic(self, state):
        "aircraft flight model"
        return AircraftP(self, state)

class AircraftP(Model):
    " aircraft performance "
    def  setup(self, aircraft, state):

        self.aircraft = aircraft
        self.wing = aircraft.wing.dynamic(state)
        self.dynamic_models = [self.wing]

        Vstall = Variable('V_{stall}', 'knots', 'Aircraft Stall Speed')
        D = Variable('D', 'N', 'Total Aircraft Drag')
        Wavg = Variable('W_{avg}', 'N', 'current average aircraft weight')
        CD = Variable('C_{D}', '-', 'Overall Drag Coefficient')
        cda = Variable("CDA", "-", "non-wing drag")
        tsfc = Variable('TSFC', '1/hr', 'Thrust Specific Fuel Consumption')

        constraints = [
            state['V'] >= Vstall,
            D >= (self.wing['D_{wing}'] + cda*0.5*state["\\rho"]
                  * state["V"]**2*self.aircraft.wing["S"]),
            CD == D/(0.5*state['\\rho']*state['V']**2
                     * self.aircraft.wing['S']),
            Wavg == (0.5*self.wing['C_{L}']*self.aircraft.wing['S']
                     * state['\\rho']*state['V']**2),
            ]

        return constraints, self.dynamic_models

class FlightSegment(Model):
    "flight segment that include cruise and climb"
    def setup(self, aircraft, N=3):

        rho = Variable('\\rho', 'kg/m^3', 'Density of air')

        with Vectorize(N):
            self.state = FlightState(rho)
            self.perf = aircraft.dynamic(self.state)
            Wstart = Variable('W_{start}', 'N', 'Segment Start Weight')
            Wend = Variable('W_{end}', 'N', 'Segment End Weight')
            Wburn = Variable("W_{burn}", "N", "fuel burned weight")
            d = Variable("d", "nautical_miles", "flight distance")
            t = Variable("t", "hrs", "flight time")
            zbre = Variable('z_{bre}', '-', 'Breguet Parameter')

        Wfuel = Variable("W_{fuel}", "N", "flight segment fuel weight")
        R = Variable("R", "nautical_miles", "flight segment range")

        constraints = [
            self.perf["W_{avg}"] == (Wstart*Wend)**0.5,
            Wfuel >= sum(Wburn),
            zbre >= self.perf["TSFC"]*t*self.perf["D"]/self.perf["W_{avg}"],
            Wburn/Wend >= te_exp_minus1(zbre, nterm=3),
            Wstart >= Wend + Wburn,
            d == t*self.state["V"],
            R == d*N]

        if N > 1:
            constraints.extend([Wend[:-1] >= Wstart[1:]])

        return constraints, self.state, self.perf

class FlightState(Model):
    " flight state, air and velocity "
    def setup(self, rho):

        V = Variable('V', 'kts', 'Aircraft Flight Speed')
        a = Variable('a', 'm/s', 'Speed of Sound')
        Rspec = Variable('R_{spec}', 287, 'J/kg/K', 'Air Specific Heat')
        gamma = Variable('\\gamma', 1.4, '-', 'Air Specific Heat Ratio')
        m = Variable('M', '-', 'Mach Number')
        Tatm = Variable("T_{atm}", "K", "air temperature")
        Tsl = Variable("T_{sl}", 288.15, "K", "sea level temperature")
        rhosl = Variable('\\rho_{sl}', 1.225, 'kg/m^3', 'Density of air')

        constraints = [a == (gamma*Rspec*Tatm)**.5,
                       V == m*a,
                       Tatm/Tsl == (rho/rhosl)**0.23490721
                      ]

        return constraints

class Wing(Model):
    " generate lift thing "
    def setup(self):

        S = Variable('S', 'm^2', 'Wing Planform Area')
        AR = Variable('AR', '-', 'Aspect Ratio')
        b = Variable('b', 'm', 'Wing Span')
        K = Variable('K', '-', 'K for Parametric Drag Model')
        e = Variable('e', '-', 'Oswald Span Efficiency Factor')
        W = Variable('W', 1, 'N', 'Wing Weight')

        constraints = [
            AR == b**2/S,
            K == (np.pi*e*AR)**-1,
            ]

        return constraints

    def dynamic(self, state):
        " wing aero "
        return WingPerf(self, state)

class WingRaymer(Model):
    """
    place holder wing model
    """
    def setup(self):

        W = Variable('W', 'N', 'Wing Weight')
        S = Variable('S', 'm^2', 'Wing Planform Area')
        AR = Variable('AR', '-', 'Aspect Ratio')
        b = Variable('b', 'm', 'Wing Span')
        K = Variable('K', '-', 'K for Parametric Drag Model')
        e = Variable('e', '-', 'Oswald Span Efficiency Factor')
        dum1 = Variable('dum1', 124.58, 'm^2')
        dum2 = Variable('dum2', 105384.1524, 'N')

        constraints = [
            (S/(dum1))**.65 * (AR/10.1)**.5 == W/(dum2),
            AR == b**2/S,
            K == (np.pi*e*AR)**-1,
            ]

        return constraints

    def dynamic(self, state):
        " wing aero "
        return WingPerf(self, state)

class WingPerf(Model):
    " wing aero "
    def setup(self, wing, state):

        CL = Variable('C_{L}', '-', 'Lift Coefficient')
        Cdw = Variable('C_{d_w}', '-', 'Cd for a NC130 Airfoil at Re=2e7')
        Dwing = Variable('D_{wing}', 'N', 'Total Wing Drag')
        Lwing = Variable('L_{wing}', 'N', 'Wing Lift')

        constraints = [
            Lwing == (0.5*wing['S']*state['\\rho']*state['V']**2)*CL,
            Cdw**6.5 >= (
                1.02458748e10 * CL**15.587947404823325
                * state['M']**156.86410659495155
                + 2.85612227e-13 * CL**1.2774976672501526
                * state['M']**6.2534328002723703
                + 2.08095341e-14 * CL**0.8825277088649582
                * state['M']**0.0273667615730107
                + 1.94411925e+06 * CL**5.6547413360261691
                * state['M']**146.51920742858428),
            Dwing >= ((0.5*wing['S']*state['\\rho']*state['V']**2)
                      * (Cdw + wing['K']*CL**2)),
            ]

        return constraints

class Mission(Model):
    " fly a mission "
    def setup(self, aircraft, Ncruise=3):

        cruise = FlightSegment(aircraft, N=Ncruise)

        mission = [cruise]

        Wtot = Variable("W_{total}", "N", "total aircraft weight")
        Wfueltot = Variable("W_{fuel-tot}", "N", "total fuel weight")

        constraints = [
            Wtot == mission[0]["W_{start}"][0],
            Wtot >= Wfueltot + aircraft["W_{dry}"],
            Wfueltot >= sum([fs["W_{fuel}"] for fs in mission]),
            mission[-1]["W_{end}"][-1] >= aircraft["W_{dry}"],
            aircraft["W_{struct}"] >= Wtot*aircraft["f_{struct}"]
            ]

        return aircraft, mission, constraints

def subbing(model, substitutions):
    " sub in dict "
    for s in substitutions:
        if hasattr(model[s], "__len__"):
            for vk in model.varkeys[s]:
                model.substitutions.update({vk: substitutions[s]})
        else:
            model.substitutions.update({s: substitutions[s]})

def init_subs(model):
    " initialize substitutions "
    subs = {"W_{payload}": 1000, "f_{struct}": 0.5, "V_{stall}": 120,
            "TSFC": 0.7, "R": 1000, "\\rho": 0.4135,
            "AR": 15, "e": 0.9, "CDA": 0.1}
    subbing(model, subs)

if __name__ == "__main__":
    wing = Wing()
    aircraft = Aircraft(wing)
    M = Mission(aircraft)
    init_subs(M)
    M.cost = M["W_{fuel-tot}"]
    sol = M.solve("mosek")
