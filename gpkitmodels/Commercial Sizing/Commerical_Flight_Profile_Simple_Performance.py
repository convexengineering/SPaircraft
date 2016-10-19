"""Simple commercial aircraft flight profile and aircraft model"""
from numpy import pi
import gpkit
import numpy as np
from gpkit import VectorVariable, Variable, Model, units, ConstraintSet, LinkedConstraintSet, SignomialsEnabled, SignomialEquality, vectorize
from gpkit.tools import te_exp_minus1
from gpkit.constraints.tight import TightConstraintSet as TCS
import matplotlib.pyplot as plt
from gpkit.small_scripts import mag
from atm_test import Atmosphere

#packages just needed for plotting since this is for sweeps
import matplotlib.pyplot as plt

#only needed for the local bounded debugging tool
from collections import defaultdict

"""
Minimizes the aircraft total fuel weight. Rate of climb equation taken from John
Anderson's Aircraft Performance and Design (eqn 5.85).

Inputs
-----

- Number of passtengers,
- Fusealge area per passenger (recommended to use 1 m^2 based on research)
- Engine weight
"""

class Aircraft(Model):
    "Aircraft"
    def __init__(self, **kwargs):
        #create submodels
        self.fuse = Fuselage()
        self.wing = Wing()
        self.engine = Engine()

        
        #variable definitions
        #number of engines
        numeng = Variable('numeng', '-', 'Number of Engines')

        constraints = []

        constraints.extend([
            numeng == numeng, #need numeng in the model
            ])

        self.components = [self.fuse, self.wing, self.engine]

        Model.__init__(self, None, [self.components + constraints], **kwargs)
        
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
    def  __init__(self, aircraft, state, **kwargs):
        #make submodels
        self.aircraft = aircraft
        self.wingP = wingP
        self.fuseP = fuseP
        self.engineP = engineP

        self.Pmodels = [self.wingP, self.fuseP, self.engineP]

        #variable definitions
        Vstall = Variable('V_{stall}', 'knots', 'Aircraft Stall Speed')
        D = Variable('D', 'N', 'Total Aircraft Drag')
        W_avg = Variable('W_{avg}', 'N', 'Geometric Average of Segment Start and End Weight')
        W_start = Variable('W_{start}', 'N', 'Segment Start Weight')
        W_end = Variable('W_{end}', 'N', 'Segment End Weight')
        W_burn = Variable('W_{burn}', 'N', 'Segment Fuel Burn Weight')
        WLoadmax = Variable('W_{Load_max}', 'N/m^2', 'Max Wing Loading')
        WLoad = Variable('W_{Load}', 'N/m^2', 'Wing Loading')
        t = Variable('tmin', 'min', 'Segment Flight Time in Minutes')
        thours = Variable('thrC', 'hour', 'Segment Flight Time in Hours')

        constraints = []

        constraints.extend([
            #speed must be greater than stall speed
            state['V'] >= Vstall,

            #compute the drag
            TCS([D >= self.wingP['D_{wing}'] + self.fuseP['D_{fuse}']]),

            #constraint CL and compute the wing loading
            W_avg == .5*self.wingP['C_{L}']*self.aircraft['S']*state['\\rho']*state['V']**2,      
            WLoad == .5*self.wingP['C_{L}']*self.aircraft['S']*state['\\rho']*state['V']**2/self.aircraft.wing['S'],

            #constraint segment start weight
            W_start >= W_end + W_burn,

            #set average weight equal to the geometric avg of start and end weight
            W_avg == (W_start * W_end)**.5,

            #constrain the max wing loading
            WLoad <= WLoadmax,

            #compute fuel burn from TSFC
            W_burn == aircraft['numeng']*self.engineP['TSFC'] * self.aircraft['thours'] * self.engineP['thrust'],
               
            #time unit conversion
            t == thours,
            ])

        Model.__init__(self, None, [self.Pmodels + constraints], **kwargs)

class ClimbP(Model):
    """
    Climb constraints
    """
    def __init__(self, aircraft, state, **kwargs):
        #submodels
        self.aircraft = aircraft
        self.aircraftP = AircraftP(aircraft, state)
        self.wingP = self.aircraftP.wingP
        self.fuseP = self.aircraftP.fuseP
        self.engineP = self.aircraftP.engineP
        
                           
        #variable definitions
        theta = Variable('\\theta', '-', 'Aircraft Climb Angle')
        excessP = Variable('excessP', 'W', 'Excess Power During Climb')
        RC = Variable('RC', 'feet/min', 'Rate of Climb/Decent')
        dhft = Variable('dhft', 'feet', 'Change in Altitude Per Climb Segment [feet]')
        RngClimb = Variable('RngClimb', 'mi', 'Down Range Covered in Each Climb Segment')

        #constraints
        constraints = []
        
        constraints.extend([
           #constraint on drag and thrust
            self.aircraft['numeng']*self.engineP['thrust'] >= D + W_avg * theta,
            
            #climb rate constraints
            TCS([excessP + state['V'] * D <=  state['V'] * aircraft['numeng'] *self.engineP['thrust']]),
            
            RC == excessP/W_avg,
            RC >= 500*units('ft/min'),
            
            #make the small angle approximation and compute theta
            theta * state['V']  == RC,
           
            dhft == self.aircraftP['tmin'] * RC,
        
            #makes a small angle assumption during climb
            RngClimb == self.aircraft['thr']*state['V'],
            ])

        Model.__init__(self, )

class CruiseP(Model):
    """
    Cruise constraints
    """
    def __init__(self, aircraft, state, **kwargs):
        self.aircraft = aircraft
        self.aircraftP = AircraftP(aircraft, aircraft)
        self.wingP = self.aircraftP.wingP
        self.fuseP = self.aircraftP.fuseP
        self.engineP = self.aircraftP.engineP
                           
        #variable definitions
        z_bre = Variable('z_{bre}', '-', 'Breguet Parameter')
        Rng = Variable('Rng', 'mi', 'Cruise Segment Range')

        constraints = []
        
        constraints.extend([
             #steady level flight constraint on D 
             D == aircraft['numeng'] * engineP['thrust'],

             #taylor series expansion to get the weight term
             TCS([W_burn/W_end >= te_exp_minus1(z_bre, nterm=3)]),

             #breguet range eqn
             TCS([z_bre >= (self.aircraft['numeng'] * self.engineP['TSFC'] * thours * D) / W_avg]),

             #time
             thours * state['V'] == aircraft['Rng'],
             ])

        Model.__init__(self, )

class CruiseSegment(Model):
    """
    Combines a flight state and aircrat to form a cruise flight segment
    """
    def __init__(self, aircraft, **kwargs):
        self.state = FlightState()
        self.cruiseP = aircraft.cruise_dynamic(self.state)

        Model.__init__(self, None, [self.state, self.cruiseP], **kwargs)

class ClimbSegment(Model):
    """
    Combines a flight state and aircrat to form a cruise flight segment
    """
    def __init__(self, aircraft, **kwargs):
        self.state = FlightState()
        self.climbP = aircraft.climb_dynamic(self.state)

        Model.__init__(self, None, [self.state, self.climbP], **kwargs)

class FlightState(Model):
    """
    creates atm model for each flight segment, has variables
    such as veloicty and altitude
    """
    def __init__(self, atmconstraints = [], **kwargs):
        #declare variables
        V = Variable('V', 'kts', 'Aircraft Flight Speed')
        a = Variable('a', 'm/s', 'Speed of Sound')
        h = Variable('h', 'm', 'Segment Altitude [meters]')
        hft = Variable('hft', 'feet', 'Segment Altitude [feet]')
        R = Variable('R', 287, 'J/kg/K', 'Air Specific Heat')
        T_atm = Variable("T_{atm}", "K", "air temperature")
        gamma = Variable('\\gamma', 1.4, '-', 'Air Specific Heat Ratio')
        M = Variable('M', '-', 'Mach Number')
        self.x = 1
        #make new constraints
        constraints = []

        constraints.extend([
            V == V, #required so velocity variable enters the model

            h == hft, #convert the units on altitude

            #compute the speed of sound with the state
            a  == (gamma * R * T_atm)**.5,

            #compute the mach number
            V == M * a,
            ])

        #build the model
        Model.__init__(self, None, constraints + atmconstraints, **kwargs)
        
class Atmosphere(FlightState):
    def __init__(self, **kwargs):
        g = Variable('g', 'm/s^2', 'Gravitational acceleration')
        p_sl = Variable("p_{sl}", "Pa", "Pressure at sea level")
        T_sl = Variable("T_{sl}", "K", "Temperature at sea level")
        L_atm = Variable("L_{atm}", "K/m", "Temperature lapse rate")
        M_atm = Variable("M_{atm}", "kg/mol",
                         "Molar mass of dry air")
        p_atm = Variable("P_{atm}", "Pa", "air pressure")
        R_atm = Variable("R_{atm}", "J/mol/K", "air specific heating value")
        TH = (g*M_atm/R_atm/L_atm).value
        rho = Variable('\\rho', 'kg/m^3', 'Density of air')
        T_atm = Variable("T_{atm}", "K", "air temperature")



        h = Variable("h", "ft", "Altitude")

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
        
        t_plus_ts_approx = (T_atm + T_s).mono_approximation({T_atm: 288.15,
                                                         T_s: T_s.value})

        with SignomialsEnabled():
            constraints = [
                # Pressure-altitude relation
                (p_atm/p_sl)**(1/5.257) == T_atm/T_sl,

                # Ideal gas law
                rho == p_atm/(R_atm/M_atm*T_atm),

                #temperature equation
                SignomialEquality(T_sl, T_atm + L_atm*h),

                #constraint on mu
##                SignomialEquality((T_atm + T_s) * mu, C_1 * T_atm**1.5),
                TCS([(T_atm + T_s) * mu >= C_1 * T_atm**1.5])
                ]


        super(Atmosphere, self).__init__(constraints)

class Engine(Model):
    """
    place holder engine model
    """
    def __init(self, **kwargs):
        #new variables
        W_engine = Variable('W_{engine}', 'N', 'Weight of a Single Turbofan Engine')
        
        constraints = []

        constraints.extend([
            W_engine == 8000 * units('N')
            ])

        Model.__init__(self, None, constraints)

    def dynamic(self, state):
            """
            returns an engine performance model
            """
            return EnginePerformance(self, state)

class EnginePerformance(Model):
    """
    place holder engine perofrmacne model
    """
    def __init__(self, engine, state, **kwargs):
            #new variables
            TSFC = Variable('TSFC', '1/hr', 'Thrust Specific Fuel Consumption')
            thrust = Variable('thrust', 'N', 'Thrust')
            #constraints
            constraints = []

            constraints.extend([
                TSFC ==  .5*units('1/hr'),

                thrust == thrust, #want thrust to enter the model
                ])

            Model.__init__(self, None, constraints)


class Wing(Model):
    """
    place holder wing model
    """
    def __init__(self, ** kwargs):
        #new variables
        W_wing = Variable('W_{wing}', 'N', 'Wing Weight')
                           
        #aircraft geometry
        S = Variable('S', 'm^2', 'Wing Planform Area')
        AR = Variable('AR', '-', 'Aspect Ratio')
        span = Variable('b', 'm', 'Wing Span')
        span_max = Variable('b_{max}', 'm', 'Max Wing Span')

        K = Variable('K', '-', 'K for Parametric Drag Model')
        e = Variable('e', '-', 'Oswald Span Efficiency Factor')
        
        constraints = []

        constraints.extend([
            #wing weight constraint
            #based off of a raymer weight and 737 data from TASOPT output file
            (S/(124.58*units('m^2')))**.65 == W_wing/(105384.1524*units('N')),

            #compute wing span and aspect ratio, subject to a span constraint
            AR == (span**2)/S,
            span <= span_max,

            #compute K for the aircraft
            K == (pi * e * AR)**-1,
            ])

        Model.__init__(self, None, constraints)

    def dynamic(self, state):
        """
        creates an instance of the wing's performance model
        """
        return WingPerformance(self, state)
        

class WingPerformance(Model):
    """
    wing aero modeling
    """
    def __init__(self, wing, state, **kwargs):
        #new variables
        CL= Variable('C_{L}', '-', 'Lift Coefficient')
        Cdw = Variable('C_{d_w}', '-', 'Cd for a NC130 Airfoil at Re=2e7')
        Dwing = Variable('D_{wing}', 'N', 'Total Wing Drag')

        #constraints
        constraints = []

        constraints.extend([
            #airfoild drag constraint
            Cdw**6.5 >= (1.02458748e10 * CL**15.587947404823325 * state['M']**156.86410659495155 +
                         2.85612227e-13 * CL**1.2774976672501526 * state['M']**6.2534328002723703 +
                         2.08095341e-14 * CL**0.8825277088649582 * state['M']**0.0273667615730107 +
                         1.94411925e+06 * CL**5.6547413360261691 * state['M']**146.51920742858428),


            Dwing >= (.5*wing['S']*state['\\rho']*state['V']**2)*(Cdw + wing['K']*CL**2)
            ])

        Model.__init__(self, None, constraints)

class Fuselage(Model):
    """
    place holder fuselage model
    """
    def __init__(self, **kwargs):
        #new variables
        n_pax = Variable('n_{pax}', '-', 'Number of Passengers to Carry')
                           
        #weight variables
        W_payload = Variable('W_{payload}', 'N', 'Aircraft Payload Weight')
        W_e = Variable('W_{e}', 'N', 'Empty Weight of Aircraft')
        W_pax = Variable('W_{pax}', 'N', 'Estimated Average Passenger Weight, Includes Baggage')

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

        Model.__init__(self, None, constraints)

    def dynamic(self, state):
        """
        returns a fuselage performance model
        """
        return FuselagePerformance(self, state)

class FuselagePerformance(Model):
    """
    Fuselage performance model
    """
    def __init__(self, fuse, state, **kwargs):
        #new variables
        Cdfuse = Variable('C_{D_{fuse}}', '-', 'Fuselage Drag Coefficient')
        Dfuse = Variable('D_{fuse}', 'N', 'Total Fuselage Drag')
        
        #constraints
        constraints = []

        constraints.extend([
            Dfuse == Cdfuse * (.5 * fuse['A_{fuse}'] * state['\\rho'] * state['V']**2),
            ])

        Model.__init__(self, None, constraints)
    

class Mission(Model):
    """
    mission class, links together all subclasses
    """
    def __init__(self, subs = None, **kwargs):
        #define the number of each flight segment
        Nclimb = 2
        Ncruise = 2

        #build required submodels
        ac = Aircraft()

        #vectorize
        with vectorize(Nclimb):
            cls = ClimbSegment(ac)

        with vectorize(Ncruise):
            crs = CruiseSegment(ac)

        #declare new variables
        W_ftotal = Variable('W_{f_{total}}', 'N', 'Total Fuel Weight')
        W_total = Variable('W_{total}', 'N', 'Total Aircraft Weight')
        hftCruise = Variable('hftCruise', 'ft', 'Cruise Altitude [feet]')

        h = cls.state['h']
        hft = cls.state['hft']
        dhft = cls.aircraftP['dhft']

        #make overall constraints
        constraints = []

        constraints.extend([
            #weight constraints
            TCS([ac['W_{e}'] + ac['W_{payload}'] + W_ftotal + ac['numeng'] * ac['W_{engine}'] + ac['W_{wing}'] <= W_total]),

            acp['W_{start}'][0] == W_total,
            cls.AircraftP['W_{end}'][-1] == crs.AircraftP['W_{start}'][0],

            TCS([ac['W_{e}'] + ac['W_{payload}'] + ac['numeng'] * engine['W_{engine }'] + ac['W_{wing}'] <= W_end[-1]]),
            TCS([W_ftotal >= sum(cls.ClimbP['W_{burn}']) + sum(crs.CruiseP['W_{burn}'])]),

            #altitude constraints
            hft[Nclimb:] == hftCruise,
            TCS([hft[1:Ncruise] >= hftClimb[:Ncruise-1] + dhftClimb]),
            TCS([hft[0] >= dhft[0]]),
            hft[Nclimb] <= hftCruise,
            
            #compute the dh
            dhftClimb == hftCruise/Nclimb,

            #constrain the thrust
            thrust[1:Nclimb] <= 2 * max(thrust[Ncruise:]),

            #set the range for each cruise segment, doesn't take credit for climb
            #down range disatnce covered
            crs.CruiseP['Rng'] == ReqRng/(Ncruise),
            ])

        vecconstraints = LinkedConstraintSet([cls, crs])

        Model.__init__(self, W_ftotal, constraints + vecconstraints, subs)

    def bound_all_variables(self, model, eps=1e-30, lower=None, upper=None):
        "Returns model with additional constraints bounding all free variables"
        lb = lower if lower else eps
        ub = upper if upper else 1/eps
        constraints = []
        freevks = tuple(vk for vk in model.varkeys if "value" not in vk.descr)
        for varkey in freevks:
            units = varkey.descr.get("units", 1)
            constraints.append([ub*units >= Variable(**varkey.descr),
                                Variable(**varkey.descr) >= lb*units])
        m = Model(model.cost, [constraints, model], model.substitutions)
        m.bound_all = {"lb": lb, "ub": ub, "varkeys": freevks}
        return m

    # pylint: disable=too-many-locals
    def determine_unbounded_variables(self, model, solver=None, verbosity=0,
                                      eps=1e-30, lower=None, upper=None, **kwargs):
        "Returns labeled dictionary of unbounded variables."
        m = self.bound_all_variables(model, eps, lower, upper)
        sol = m.localsolve(solver, verbosity, **kwargs)
        solhold = sol
        lam = sol["sensitivities"]["la"][1:]
        out = defaultdict(list)
        for i, varkey in enumerate(m.bound_all["varkeys"]):
            lam_gt, lam_lt = lam[2*i], lam[2*i+1]
            if abs(lam_gt) >= 1e-7:  # arbitrary threshold
                out["sensitive to upper bound"].append(varkey)
            if abs(lam_lt) >= 1e-7:  # arbitrary threshold
                out["sensitive to lower bound"].append(varkey)
            value = mag(sol["variables"][varkey])
            distance_below = np.log(value/m.bound_all["lb"])
            distance_above = np.log(m.bound_all["ub"]/value)
            if distance_below <= 3:  # arbitrary threshold
                out["value near lower bound"].append(varkey)
            elif distance_above <= 3:  # arbitrary threshold
                out["value near upper bound"].append(varkey)
        return out, solhold

if __name__ == '__main__':
    substitutions = {      
            'V_{stall}': 120,
            'ReqRng': 1000, #('sweep', np.linspace(500,2000,4)),
            'hftCruise': 30000, #('sweep', np.linspace(20000,40000,4)),
            'numeng': 2,
            'W_{Load_max}': 6664,
            'W_{engine}': 1000,
            'W_{pax}': 91 * 9.81,
            'n_{pax}': 150,
            'pax_{area}': 1,
            'C_{d_fuse}': .005, #assumes flat plate turbulent flow, from wikipedia
            'e': .9,
            'span_{max}': 35,

            #atm subs
            "p_{sl}": 101325,
            "T_{sl}": 288.15,
            "L_{atm}": .0065,
            "M_{atm}":.0289644,
            "R_{atm}": 8.31447,
            }
           
    m = Mission(substitutions)
    m.localsolve(solver='mosek')
