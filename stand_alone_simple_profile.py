"""Simple commercial aircraft flight profile and aircraft model"""
from numpy import pi
import numpy as np
from gpkit import Variable, Model, units, SignomialsEnabled, Vectorize
from gpkit.constraints.sigeq import SignomialEquality as SignomialEquality
from gpkit.tools import te_exp_minus1
from gpkit.constraints.tight import Tight as TCS
import matplotlib.pyplot as plt

"""
Models required to minimize the aircraft total fuel weight. Rate of climb equation taken from John
Anderson's Aircraft Performance and Design (eqn 5.85).
Inputs
-----
- Number of passengers
- Passenger weight [N]
- Fuselage area per passenger (recommended to use 1 m^2 based on research) [m^2]
- Engine weight [N]
- Number of engines
- Required mission range [nm]
- Oswald efficiency factor
- Max allowed wing span [m]
- Cruise altitude [ft]
"""

class Aircraft(Model):
    "Aircraft class"
    def  setup(self, **kwargs):
        #create submodels
        self.fuse = Fuselage()
        self.wing = Wing()
        self.engine = Engine()

        #variable definitions
        numeng = Variable('numeng', '-', 'Number of Engines')

        self.components = [self.fuse, self.wing, self.engine]

        return self.components

    def climb_dynamic(self, state):
        """
        creates an aircraft climb performance model, given a state
        """
        return ClimbP(self, state)

    def cruise_dynamic(self, state):
        """
        creates an aircraft cruise performance model, given a state
        """
        return CruiseP(self, state)

class AircraftP(Model):
    """
    aircraft performance models superclass, contains constraints true for
    all flight segments
    """
    def  setup(self, aircraft, state, **kwargs):
        #make submodels
        self.aircraft = aircraft
        self.wingP = aircraft.wing.dynamic(state)
        self.fuseP = aircraft.fuse.dynamic(state)
        self.engineP = aircraft.engine.dynamic(state)

        self.Pmodels = [self.wingP, self.fuseP, self.engineP]

        #variable definitions
        Vstall = Variable('V_{stall}', 'knots', 'Aircraft Stall Speed')
        D = Variable('D', 'N', 'Total Aircraft Drag')
        W_avg = Variable('W_{avg}', 'N', 'Geometric Average of Segment Start and End Weight')
        W_start = Variable('W_{start}', 'N', 'Segment Start Weight')
        W_end = Variable('W_{end}', 'N', 'Segment End Weight')
        W_burn = Variable('W_{burn}', 'N', 'Segment Fuel Burn Weight')
        WLoadmax = Variable('W_{Load_{max}}', 'N/m^2', 'Max Wing Loading')
        WLoad = Variable('W_{Load}', 'N/m^2', 'Wing Loading')
        t = Variable('tmin', 'min', 'Segment Flight Time in Minutes')
        thours = Variable('thr', 'hour', 'Segment Flight Time in Hours')

        constraints = []

        constraints.extend([
            #speed must be greater than stall speed
            state['V'] >= Vstall,


            #Figure out how to delete
            Vstall == 120*units('kts'),
            WLoadmax == 6664 * units('N/m^2'),

            #compute the drag
            TCS([D >= self.wingP['D_{wing}'] + self.fuseP['D_{fuse}']]),

            #constraint CL and compute the wing loading
            W_avg == .5*self.wingP['C_{L}']*self.aircraft['S']*state.atm['\\rho']*state['V']**2,
            WLoad == .5*self.wingP['C_{L}']*self.aircraft['S']*state.atm['\\rho']*state['V']**2/self.aircraft.wing['S'],

            #set average weight equal to the geometric avg of start and end weight
            W_avg == (W_start * W_end)**.5,

            #constrain the max wing loading
            WLoad <= WLoadmax,

            #compute fuel burn from TSFC
            W_burn == aircraft['numeng']*self.engineP['TSFC'] * thours * self.engineP['thrust'],

            #time unit conversion
            t == thours,

            #make lift equal weight --> small angle approx in climb
            self.wingP['L_{wing}'] == W_avg,
            ])

        return constraints, self.wingP, self.engineP, self.fuseP

class ClimbP(Model):
    """
    Climb constraints
    """
    def setup(self, aircraft, state, **kwargs):
        #submodels
        self.aircraft = aircraft
        self.aircraftP = AircraftP(aircraft, state)
        self.wingP = self.aircraftP.wingP
        self.fuseP = self.aircraftP.fuseP
        self.engineP = self.aircraftP.engineP

        #variable definitions
        theta = Variable('\\theta', '-', 'Aircraft Climb Angle')
        excessP = Variable('P_{excess}', 'W', 'Excess Power During Climb')
        RC = Variable('RC', 'feet/min', 'Rate of Climb/Decent')
        dhft = Variable('dhft', 'feet', 'Change in Altitude Per Climb Segment [feet]')
        RngClimb = Variable('R_{climb}', 'nautical_miles', 'Down Range Covered in Each Climb Segment')

        #constraints
        constraints = []

        constraints.extend([
           #constraint on drag and thrust
            self.aircraft['numeng']*self.engineP['thrust'] >= self.aircraftP['D'] + self.aircraftP['W_{avg}'] * theta,

            #climb rate constraints
            TCS([excessP + state['V'] * self.aircraftP['D'] <=  state['V'] * aircraft['numeng'] * self.engineP['thrust']]),

            RC == excessP/self.aircraftP['W_{avg}'],
            RC >= 500*units('ft/min'),

            #make the small angle approximation and compute theta
            theta * state['V']  == RC,

            dhft == self.aircraftP['tmin'] * RC,

            #makes a small angle assumption during climb
            RngClimb == self.aircraftP['thr']*state['V'],
            ])

        return constraints, self.aircraftP

class CruiseP(Model):
    """
    Cruise constraints
    """
    def setup(self, aircraft, state, **kwargs):
        self.aircraft = aircraft
        self.aircraftP = AircraftP(aircraft, state)
        self.wingP = self.aircraftP.wingP
        self.fuseP = self.aircraftP.fuseP
        self.engineP = self.aircraftP.engineP

        #variable definitions
        z_bre = Variable('z_{bre}', '-', 'Breguet Parameter')
        Rng = Variable('Rng', 'nautical_miles', 'Cruise Segment Range')

        constraints = []

        constraints.extend([
             #steady level flight constraint on D
             self.aircraftP['D'] == aircraft['numeng'] * self.engineP['thrust'],

             #taylor series expansion to get the weight term
             TCS([self.aircraftP['W_{burn}']/self.aircraftP['W_{end}'] >=
                  te_exp_minus1(z_bre, nterm=3)]),

             #breguet range eqn
             # old version -- possibly unneeded numeng
 #            TCS([z_bre >= (self.aircraft['numeng'] * self.engineP['TSFC'] * self.aircraftP['thr']*
 #                           self.aircraftP['D']) / self.aircraftP['W_{avg}']]),

            # new version -- needs to be thought through carefully
             # seems correct to me - I switched T to D below (steady level flight) but fogot
             #about the Negn term
             TCS([z_bre >= (self.engineP['TSFC'] * self.aircraftP['thr']*
                            self.aircraftP['D']) / self.aircraftP['W_{avg}']]),

             #time
             self.aircraftP['thr'] * state['V'] == Rng,
             ])

        return constraints, self.aircraftP

class CruiseSegment(Model):
    """
    Combines a flight state and aircrat to form a cruise flight segment
    """
    def setup(self, aircraft, **kwargs):
        self.state = FlightState()
        self.cruiseP = aircraft.cruise_dynamic(self.state)

        return self.state, self.cruiseP

class ClimbSegment(Model):
    """
    Combines a flight state and aircrat to form a cruise flight segment
    """
    def setup(self, aircraft, **kwargs):
        self.state = FlightState()
        self.climbP = aircraft.climb_dynamic(self.state)

        return self.state, self.climbP

class FlightState(Model):
    """
    creates atm model for each flight segment, has variables
    such as veloicty and altitude

    SKIP VERIFICATION
    """
    def setup(self,**kwargs):
        #make an atmosphere model
        self.alt = Altitude()
        self.atm = Atmosphere(self.alt)

        #declare variables
        V = Variable('V', 'kts', 'Aircraft Flight Speed')
        a = Variable('a', 'm/s', 'Speed of Sound')

        R = Variable('R', 287, 'J/kg/K', 'Air Specific Heat')
        gamma = Variable('\\gamma', 1.4, '-', 'Air Specific Heat Ratio')
        M = Variable('M', '-', 'Mach Number')

        #make new constraints
        constraints = []

        constraints.extend([
            V == V, #required so velocity variable enters the model

            #compute the speed of sound with the state
            a  == (gamma * R * self.atm['T_{atm}'])**.5,

            #compute the mach number
            V == M * a,
            ])

        #build the model
        return constraints, self.alt, self.atm

class Altitude(Model):
    """
    holds the altitdue variable

    Upper Unbounded
    ---------------
    h

    Lower Unbounded
    ---------------
    h
    """
    def setup(self, **kwargs):
        #define altitude variables
        h = self.h = Variable('h', 'm', 'Segment Altitude [meters]')
        hft = Variable('hft', 'feet', 'Segment Altitude [feet]')

        constraints = []

        constraints.extend([
            h == hft, #convert the units on altitude
            ])

        return constraints


class Atmosphere(Model):
    """

    Upper Unbounded
    ---------------  # TOOD: h is considered unbounded because SignomialEqualities are not detected
    rho, T_atm, h

    Lower Unbounded
    ---------------
    rho
    """
    def setup(self, alt, **kwargs):
        p_sl = Variable("p_{sl}", 101325, "Pa", "Pressure at sea level")
        T_sl = Variable("T_{sl}", 288.15, "K", "Temperature at sea level")
        L_atm = Variable("L_{atm}", .0065, "K/m", "Temperature lapse rate")
        M_atm = Variable("M_{atm}", .0289644, "kg/mol",
                         "Molar mass of dry air")
        p_atm = Variable("P_{atm}", "Pa", "air pressure")
        R_atm = Variable("R_{atm}", 8.31447, "J/mol/K", "air specific heating value")
        TH = 5.257386998354459 #(g*M_atm/R_atm/L_atm).value
        rho = self.rho = Variable('\\rho', 'kg/m^3', 'Density of air')
        T_atm = self.T_atm = Variable("T_{atm}", "K", "air temperature")

        h = self.h = alt['h']

        """
        Dynamic viscosity (mu) as a function of temperature
        References:
        http://www-mdp.eng.cam.ac.uk/web/library/enginfo/aerothermal_dvd_only/aero/
            atmos/atmos.html
        http://www.cfd-online.com/Wiki/Sutherland's_law
        """
        mu  = Variable('\\mu', 'kg/(m*s)', 'Dynamic viscosity')

        T_s = Variable('T_s', 110.4, "K", "Sutherland Temperature")
        C_1 = Variable('C_1', 1.458E-6, "kg/(m*s*K^0.5)",
                       'Sutherland coefficient')

        with SignomialsEnabled():
            constraints = [
                # Pressure-altitude relation
                (p_atm/p_sl)**(1/5.257) == T_atm/T_sl,

                # Ideal gas law
                rho == p_atm/(R_atm/M_atm*T_atm),

                #temperature equation
                SignomialEquality(T_sl, T_atm + L_atm*h),

                #constraint on mu
                mu == C_1 * T_atm**1.5 /(6.64*units('K^.28')*T_s**0.72),
                ]

        #like to use a local subs here in the future
        subs = None

        return constraints

class Engine(Model):
    """
    place holder engine model
    """
    def setup(self, **kwargs):
        #new variables
        W_engine = Variable('W_{engine}', 1000, 'N', 'Weight of a Single Turbofan Engine')

        constraints = []

        constraints.extend([
            W_engine == W_engine
            ])

        return constraints

    def dynamic(self, state):
            """
            returns an engine performance model
            """
            return EnginePerformance(self, state)

class EnginePerformance(Model):
    """
    place holder engine perofrmacne model
    """
    def setup(self, engine, state, **kwargs):
        #new variables
        TSFC = Variable('TSFC', '1/hr', 'Thrust Specific Fuel Consumption')
        thrust = Variable('thrust', 'N', 'Thrust')

        #constraints
        constraints = []

        constraints.extend([
            TSFC == TSFC,

            thrust == thrust, #want thrust to enter the model
            ])

        return constraints


class Wing(Model):
    """
    place holder wing model
    """
    def setup(self, ** kwargs):
        #new variables
        W_wing = Variable('W_{wing}', 'N', 'Wing Weight')

        #aircraft geometry
        S = Variable('S', 'm^2', 'Wing Planform Area')
        AR = Variable('AR', '-', 'Aspect Ratio')
        span = Variable('b', 'm', 'Wing Span')
        span_max = Variable('b_{max}', 'm', 'Max Wing Span')

        K = Variable('K', '-', 'K for Parametric Drag Model')
        e = Variable('e', '-', 'Oswald Span Efficiency Factor')

        dum1 = Variable('dum1', 124.58, 'm^2')
        dum2 = Variable('dum2', 105384.1524, 'N')

        constraints = []

        constraints.extend([
            #wing weight constraint
            #based off of a raymer weight and 737 data from TASOPT output file
            (S/(dum1))**.65 == W_wing/(dum2),

            #compute wing span and aspect ratio, subject to a span constraint
            AR == (span**2)/S,
            #AR == 9,

            span <= span_max,

            #compute K for the aircraft
            K == (pi * e * AR)**-1,
            ])

        return constraints

    def dynamic(self, state):
        """
        creates an instance of the wing's performance model
        """
        return WingPerformance(self, state)


class WingPerformance(Model):
    """
    wing aero modeling
    """
    def setup(self, wing, state, **kwargs):
        #new variables
        CL= Variable('C_{L}', '-', 'Lift Coefficient')
        Cdw = Variable('C_{d_w}', '-', 'Cd for a NC130 Airfoil at Re=2e7')
        Dwing = Variable('D_{wing}', 'N', 'Total Wing Drag')
        Lwing = Variable('L_{wing}', 'N', 'Wing Lift')

        #constraints
        constraints = []

        constraints.extend([
            #airfoil drag constraint
            Lwing == (.5*wing['S']*state.atm['\\rho']*state['V']**2)*CL,
            TCS([Cdw**6.5 >= (1.02458748e10 * CL**15.587947404823325 * state['M']**156.86410659495155 +
                         2.85612227e-13 * CL**1.2774976672501526 * state['M']**6.2534328002723703 +
                         2.08095341e-14 * CL**0.8825277088649582 * state['M']**0.0273667615730107 +
                         1.94411925e+06 * CL**5.6547413360261691 * state['M']**146.51920742858428)]),
            TCS([Dwing >= (.5*wing['S']*state.atm['\\rho']*state['V']**2)*(Cdw + wing['K']*CL**2)]),
            ])

        return constraints

class Fuselage(Model):
    """
    place holder fuselage model
    """
    def setup(self, **kwargs):
        #new variables
        n_pax = Variable('n_{pass}', '-', 'Number of Passengers to Carry')

        #weight variables
        W_payload = Variable('W_{payload}', 'N', 'Aircraft Payload Weight')
        W_e = Variable('W_{e}', 'N', 'Empty Weight of Aircraft')
        W_pax = Variable('W_{pass}', 'N', 'Estimated Average Passenger Weight, Includes Baggage')

        A_fuse = Variable('A_{fuse}', 'm^2', 'Estimated Fuselage Area')
        pax_area = Variable('pax_{area}', 'm^2', 'Estimated Fuselage Area per Passenger')

        constraints = []

        constraints.extend([
            #compute fuselage area for drag approximation
            A_fuse == pax_area * n_pax,

            #constraints on the various weights
            W_payload == n_pax * W_pax,

            #estimate based on TASOPT 737 model
            W_e == .75*W_payload,
            ])

        return constraints

    def dynamic(self, state):
        """
        returns a fuselage performance model
        """
        return FuselagePerformance(self, state)

class FuselagePerformance(Model):
    """
    Fuselage performance model
    """
    def setup(self, fuse, state, **kwargs):
        #new variables
        Cdfuse = Variable('C_{D_{fuse}}', '-', 'Fuselage Drag Coefficient')
        Dfuse = Variable('D_{fuse}', 'N', 'Total Fuselage Drag')

        #constraints
        constraints = []

        constraints.extend([
            Dfuse == Cdfuse * (.5 * fuse['A_{fuse}'] * state.atm['\\rho'] * state['V']**2),

            Cdfuse == .005,
            ])

        return constraints


class Mission(Model):
    """
    mission class, links together all subclasses
    """
    def setup(self, ac, substitutions = None, **kwargs):
        #define the number of each flight segment
        Nclimb = 2
        Ncruise = 2

        #Vectorize
        with Vectorize(Nclimb):
            climb = ClimbSegment(ac)

        with Vectorize(Ncruise):
            cruise = CruiseSegment(ac)

        #declare new variables
        W_ftotal = Variable('W_{f_{total}}', 'N', 'Total Fuel Weight')
        W_fclimb = Variable('W_{f_{climb}}', 'N', 'Fuel Weight Burned in Climb')
        W_fcruise = Variable('W_{f_{cruise}}', 'N', 'Fuel Weight Burned in Cruise')
        W_total = Variable('W_{total}', 'N', 'Total Aircraft Weight')
        CruiseAlt = Variable('CruiseAlt', 'ft', 'Cruise Altitude [feet]')
        ReqRng = Variable('R_{req}', 'nautical_miles', 'Required Cruise Range')

        h = climb['h']
        hftClimb = climb['hft']
        dhft = climb['dhft']
        hftCruise = cruise['hft']

        #make overall constraints
        constraints = []

        constraints.extend([
            #weight constraints
            TCS([ac['W_{e}'] + ac['W_{payload}'] + W_ftotal + ac['numeng'] * ac['W_{engine}'] + ac['W_{wing}'] <= W_total]),

            climb['W_{start}'][0] == W_total,
            climb['W_{end}'][-1] == cruise['W_{start}'][0],

            # similar constraint 1
            TCS([climb['W_{start}'] >= climb['W_{end}'] + climb['W_{burn}']]),
            # similar constraint 2
            TCS([cruise['W_{start}'] >= cruise['W_{end}'] + cruise['W_{burn}']]),

            climb['W_{start}'][1:] == climb['W_{end}'][:-1],
            cruise['W_{start}'][1:] == cruise['W_{end}'][:-1],

            TCS([ac['W_{e}'] + ac['W_{payload}'] + ac['numeng'] * ac['W_{engine}'] + ac['W_{wing}'] <= cruise['W_{end}'][-1]]),

            TCS([W_ftotal >=  W_fclimb + W_fcruise]),
            TCS([W_fclimb >= sum(climb['W_{burn}'])]),
            TCS([W_fcruise >= sum(cruise['W_{burn}'])]),

            #altitude constraints
            hftCruise == CruiseAlt,
            TCS([hftClimb[1:Ncruise] >= hftClimb[:Ncruise-1] + dhft]),
            TCS([hftClimb[0] >= dhft[0]]),
            hftClimb[-1] <= hftCruise,

            #compute the dh
            dhft == hftCruise/Nclimb,

            #constrain the thrust
            climb['thrust'] <= 2 * max(cruise['thrust']),

            #set the range for each cruise segment, doesn't take credit for climb
            #down range disatnce covered
            cruise.cruiseP['Rng'] == ReqRng/(Ncruise),

            #set the TSFC
            climb['TSFC'] == .7*units('1/hr'),
            cruise['TSFC'] == .5*units('1/hr'),
            ])

        # Model.setup(self, W_ftotal + s*units('N'), constraints + ac + climb + cruise, subs)
        return constraints + ac + climb + cruise


if __name__ == '__main__':
    #build required submodels
    ac = Aircraft()

    substitutions = {
##            'V_{stall}': 120,
            'R_{req}': 500, #('sweep', np.linspace(500,2000,4)),
            'CruiseAlt': 30000, #('sweep', np.linspace(20000,40000,4)),
            'numeng': 1,
##            'W_{Load_{max}}': 6664,
            'W_{pass}': 91 * 9.81,
            'n_{pass}': 150,
            'pax_{area}': 1,
##            'C_{D_{fuse}}': .005, #assumes flat plate turbulent flow, from wikipedia
            'e': .9,
            'b_{max}': 35,
            }

    mission = Mission(ac)
    m = Model(mission['W_{f_{total}}'], mission, substitutions)
    sol = m.localsolve(solver='mosek', verbosity = 4)

    substitutions = {
##            'V_{stall}': 120,
            'R_{req}': ('sweep', np.linspace(500,3000,8)),
            'CruiseAlt': 30000, #('sweep', np.linspace(20000,40000,4)),
            'numeng': 1,
##            'W_{Load_{max}}': 6664,
            'W_{pass}': 91 * 9.81,
            'n_{pass}': 150,
            'pax_{area}': 1,
##            'C_{D_{fuse}}': .005, #assumes flat plate turbulent flow, from wikipedia
            'e': .9,
            'b_{max}': 35,
            }
    mission = Mission(ac)
    m = Model(mission['W_{f_{total}}'], mission, substitutions)
    solRsweep = m.localsolve(solver='mosek', verbosity = 4)

##    plt.plot(solRsweep('R_{req}'), solRsweep('W_{f_{total}}'), '-r')
##    plt.xlabel('Mission Range [nm]')
##    plt.ylabel('Total Fuel Burn [N]')
##    plt.title('Fuel Burn vs Range')
##    plt.show()

    substitutions = {
##            'V_{stall}': 120,
            'R_{req}': 2000,#('sweep', np.linspace(500,2000,4)),
            'CruiseAlt': ('sweep', np.linspace(20000,40000,8)),
            'numeng': 1,
##            'W_{Load_{max}}': 6664,
            'W_{pass}': 91 * 9.81,
            'n_{pass}': 150,
            'pax_{area}': 1,
##            'C_{D_{fuse}}': .005, #assumes flat plate turbulent flow, from wikipedia
            'e': .9,
            'b_{max}': 35,
            }

    mmission = Mission(ac)
    m = Model(mission['W_{f_{total}}'], mission, substitutions)
    solAltsweep = m.localsolve(solver='mosek', verbosity = 4)

##    plt.plot(solAltsweep('CruiseAlt'), solAltsweep('W_{f_{total}}'), '-r')
##    plt.xlabel('Cruise Alt [ft]')
##    plt.ylabel('Total Fuel Burn [N]')
##    plt.title('Fuel Burn vs Range')
##    plt.show()
