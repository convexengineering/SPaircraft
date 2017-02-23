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

        with Vectorize(N):
            self.state = FlightState()
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
    def setup(self):

        rho = Variable('\\rho', 'kg/m^3', 'Density of air')
        V = Variable('V', 'kts', 'Aircraft Flight Speed')
        a = Variable('a', 'm/s', 'Speed of Sound')
        Rspec = Variable('R_{spec}', 287, 'J/kg/K', 'Air Specific Heat')
        gamma = Variable('\\gamma', 1.4, '-', 'Air Specific Heat Ratio')
        m = Variable('M', '-', 'Mach Number')
        Tatm = Variable("T_{atm}", "K", "air temperature")
        Tsl = Variable("T_{sl}", 288.15, "K", "sea level temperature")
        rhosl = Variable('\\rho_{sl}', 1.225, 'kg/m^3', 'Density of air')
        mu = Variable("\\mu", 1.45e-5, "kg/m/s", "dynamic viscosity")

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
        return WingPerfMach(self, state)

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
        return WingPerfMach(self, state)

class WingBox(Model):
    """
    Structural model for a wing
    source: Hoburg, "Geometric Programming for Aircraft Design Optimization"
    Note - does not have a performance model
    """
    def setup(self):

        Icap = Variable('I_{cap}', '-',
                        'Non-dim spar cap area moment of inertia')
        Mr = Variable('M_r', 'N', 'Root moment per root chord')
        nu = Variable('\\nu', '-', '$(\\lambda^2 + \\lambda + 1)/(\\lambda+1)$')
        Wcap = Variable('W_{cap}', 'N', 'Weight of spar caps')
        Wweb = Variable('W_{web}', 'N', 'Weight of shear web')
        Wstruct = Variable('W', 'N', 'Structural weight')
        lam = Variable('\\lambda', '-', 'Taper ratio')
        fwadd = Variable('f_{w,add}', '-', 'Wing added weight fraction')
        g = Variable('g', 9.81, 'm/s^2', 'Gravitational acceleration')
        Nlift = Variable('N_{lift}', '-', 'Wing loading multiplier')
        rh = Variable('r_h', 0.75, '-',
                      'Fractional wing thickness at spar web')
        rhocap = Variable('\\rho_{cap}', 2700, 'kg/m^3',
                          'Density of spar cap material')
        rhoweb = Variable('\\rho_{web}', 2700, 'kg/m^3',
                          'Density of shear web material')
        sigmax = Variable('\\sigma_{max}', 250e6, 'Pa',
                          'Allowable tensile stress')
        sigmaxshear = Variable('\\sigma_{max,shear}', 167e6, 'Pa',
                               'Allowable shear stress')
        wwb = Variable('wwb', 0.5, '-', 'Wingbox-width-to-chord ratio')
        tcap = Variable('t_{cap}', '-', 'Non-dim. spar cap thickness')
        tweb = Variable('t_{web}', '-', 'Non-dim. shear web thickness')
        p = Variable('p', '-', 'Substituted variable = 1 + 2*lam')
        q = Variable('q', '-', 'Substituted variable = 1 + lam')
        ymac = Variable('y_{mac}', 'm',
                        'Spanwise location of mean aerodynamic chord')
        rho0 = Variable('\\rho_0', 1.225, 'kg/m^3', 'Air density (0 ft)')
        AR = Variable('AR', '-', 'Wing aspect ratio')
        Lmax = Variable('L_{max}', 'N', 'Maximum load')
        Sw = Variable('S', 'm^2', 'Wing area')
        Vne = Variable('V_{ne}', 60, 'm/s', 'Never exceed velocity')
        b = Variable('b', 'm', 'Wing span')
        CLwmax = Variable('C_{L_{wmax}}', 1.2, '-', 'Max lift coefficient wing')
        tau = Variable('\\tau', '-', 'Wing thickness/chord ratio')
        mac = Variable('mac', 'm', 'Mean aerodynamic chord (wing)')
        e = Variable('e', '-', 'Oswald efficiency factor')

        constraints = [
            AR == b**2/Sw,
            mac == Sw/b,
            p >= 1 + 2*lam,
            2*q >= 1 + p,
            ymac == (b/3)*q/p,
            Lmax == 0.5*rho0*Vne**2*Sw*CLwmax,
            Mr >= Lmax*AR*p/24,
            0.92*wwb*tau*tcap**2 + Icap <= 0.92**2/2*wwb*tau**2*tcap,
            8 >= Nlift*Mr*AR*q**2*tau/(Sw*Icap*sigmax),
            12 >= AR*Lmax*Nlift*q**2/(tau*Sw*tweb*sigmaxshear),
            nu**3.94 >= 0.86*p**(-2.38)+ 0.14*p**0.56,
            Wcap >= 8*rhocap*g*wwb*tcap*Sw**1.5*nu/(3*AR**0.5),
            TCS([Wweb >= 8*rhoweb*g*rh*tau*tweb*Sw**1.5*nu/(3*AR**0.5)]),
            Wstruct >= (1 + fwadd)*(Wweb + Wcap),
            ]

        return constraints

    def dynamic(self, state):
        " wing aero "
        return WingPerfRe(self, state)


class WingPerfMach(Model):
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

class WingPerfRe(Model):
    """
    Wing performance model
    """
    def setup(self, wing, state):
        self.wing = wing

        Re = Variable('Re_w', '-', 'Reynolds number (wing)')
        CDp = Variable('C_{D_{p_w}}', '-',
                       'Wing parasitic drag coefficient')
        CDw = Variable('C_{d_w}', '-', 'Drag coefficient, wing')
        CLw = Variable('C_{L}', '-', 'Lift coefficient, wing')
        D = Variable('D_{wing}', 'N', 'Wing drag')
        Lw = Variable('L_w', 'N', 'Wing lift')

        constraints = [
            Lw == 0.5*state['\\rho']*state['V']**2*self.wing['S']*CLw,
            D == 0.5*state['\\rho']*state['V']**2*self.wing['S']*CDw,
            CDw >= CDp + CLw**2/(np.pi*self.wing['e']*self.wing['AR']),
            Re == state['\\rho']*state['V']*self.wing['mac']/state['\\mu'],
            1 >= (
                2.56*CLw**5.88/(Re**1.54*self.wing['\\tau']**3.32*CDp**2.62)
                +3.8e-9*self.wing['\\tau']**6.23/(CLw**0.92*Re**1.38*CDp**9.57)
                +2.2e-3*Re**0.14*self.wing['\\tau']**0.033/(CLw**0.01*CDp**0.73)
                +6.14e-6*CLw**6.53/(Re**0.99*self.wing['\\tau']**0.52*CDp**5.19)
                +1.19e4*CLw**9.78*self.wing['\\tau']**1.76/(Re*CDp**0.91)),
            ]

        return constraints

class Mission(Model):
    " fly a mission "
    def setup(self, aircraft, Ncruise=3):

        cruise = FlightSegment(aircraft, N=Ncruise)

        mission = [cruise]

        Wtot = Variable("W_{total}", "N", "total aircraft weight")
        Wfueltot = Variable("W_{fuel-tot}", "N", "total fuel weight")
        WLmax = Variable("W_{Load_max}", 6664, "N/m^2", "max wing loading")

        constraints = [
            Wtot == mission[0]["W_{start}"][0],
            Wtot >= Wfueltot + aircraft["W_{dry}"],
            Wfueltot >= sum([fs["W_{fuel}"] for fs in mission]),
            mission[-1]["W_{end}"][-1] >= aircraft["W_{dry}"],
            aircraft["W_{struct}"] >= Wtot*aircraft["f_{struct}"],
            WLmax >= Wtot/aircraft["S"]
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
            "TSFC": 0.7, "R": 1000, "AR": 15, "e": 0.9, "CDA": 0.1,
            "\\rho": 0.4, "\\lambda": 0.5, "f_{w,add}": 0.4, "N_{lift}": 2.0,
            "\\tau": 0.15, "e": 0.9}
    subbing(model, subs)

if __name__ == "__main__":
    # wing = Wing()
    wing = WingBox()
    aircraft = Aircraft(wing)
    M = Mission(aircraft)
    init_subs(M)
    M.cost = M["W_{fuel-tot}"]
    sol = M.solve("mosek")
