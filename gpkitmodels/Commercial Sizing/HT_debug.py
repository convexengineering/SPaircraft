"""Simple commercial aircraft flight profile and aircraft model"""
from numpy import pi, cos, tan
import numpy as np
from gpkit import Variable, Model, units, SignomialsEnabled, vectorize, VectorVariable
from gpkit.constraints.sigeq import SignomialEqualityConstraint as SignomialEquality
from gpkit.tools import te_exp_minus1
from gpkit.constraints.tight import TightConstraintSet as TCS
import matplotlib.pyplot as plt
#only needed for the local bounded debugging tool
from collections import defaultdict
from gpkit.small_scripts import mag

"""
Models requird to minimze the aircraft total fuel weight. Rate of climb equation taken from John
Anderson's Aircraft Performance and Design (eqn 5.85).
Inputs
-----
- Number of passtengers
- Passegner weight [N]
- Fusealge area per passenger (recommended to use 1 m^2 based on research) [m^2]
- Engine weight [N]
- Number of engines
- Required mission range [nm]
- Oswald efficiency factor
- Max allowed wing span [m]
- Cruise altitude [ft]
"""

class Aircraft(Model):
    "Aircraft class"
    def __init__(self, **kwargs):
        #create submodels
        self.fuse = Fuselage()
        self.wing = Wing()
        self.engine = Engine()
        self.ht = HorizontalTail(self.fuse, self.wing)

        #variable definitions
        numeng = Variable('numeng', '-', 'Number of Engines')

        constraints = []

        constraints.extend([
            numeng == numeng, #need numeng in the model
            ])

        self.components = [self.fuse, self.wing, self.engine, self.ht]

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
        self.wingP = aircraft.wing.dynamic(state)
        self.fuseP = aircraft.fuse.dynamic(state)
        self.engineP = aircraft.engine.dynamic(state)
        self.htP = aircraft.ht.dynamic(aircraft.fuse, aircraft.wing, self.fuseP, self.wingP, state)

        self.Pmodels = [self.wingP, self.fuseP, self.engineP, self.htP]

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
        thours = Variable('thr', 'hour', 'Segment Flight Time in Hours')

        constraints = []
        with SignomialsEnabled():
            constraints.extend([
                #speed must be greater than stall speed
                state['V'] >= Vstall,

                #Figure out how to delete
                Vstall == 120*units('kts'),
                WLoadmax == 6664 * units('N/m^2'),

                #compute the drag
                TCS([D >= self.wingP['D_{wing}'] + self.fuseP['D_{fuse}'] + self.htP['D_{ht}']]),

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
                self.wingP['L_{wing}'] >= W_avg + self.htP['L_h'],
##                SignomialEquality(self.wingP['L_{wing}'] , W_avg + self.htP['L_h'])
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
        RngClimb = Variable('RngClimb', 'nautical_miles', 'Down Range Covered in Each Climb Segment')

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

        Model.__init__(self, None, constraints + self.aircraftP)

class CruiseP(Model):
    """
    Cruise constraints
    """
    def __init__(self, aircraft, state, **kwargs):
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

            # new version -- needs to be thought through carefully
             # seems correct to me - I switched T to D below (steady level flight) but fogot
             #about the Negn term
             TCS([z_bre >= (self.engineP['TSFC'] * self.aircraftP['thr']*
                            self.aircraftP['D']) / self.aircraftP['W_{avg}']]),

             #time
             self.aircraftP['thr'] * state['V'] == Rng,
             ])

        Model.__init__(self, None, constraints + self.aircraftP)

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
    def __init__(self,**kwargs):
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

            M == .8,
            ])

        #build the model
        Model.__init__(self, None, constraints + self.atm + self.alt, **kwargs)

class Altitude(Model):
    """
    holds the altitdue variable
    """
    def __init__(self, **kwargs):
        #define altitude variables
        h = Variable('h', 'm', 'Segment Altitude [meters]')
        hft = Variable('hft', 'feet', 'Segment Altitude [feet]')

        constraints = []

        constraints.extend([
            h == hft, #convert the units on altitude
            ])

        Model.__init__(self, None, constraints, **kwargs)
            

class Atmosphere(Model):
    def __init__(self, alt, **kwargs):
        g = Variable('g', 'm/s^2', 'Gravitational acceleration')
        p_sl = Variable("p_{sl}", 101325, "Pa", "Pressure at sea level")
        T_sl = Variable("T_{sl}", 288.15, "K", "Temperature at sea level")
        L_atm = Variable("L_{atm}", .0065, "K/m", "Temperature lapse rate")
        M_atm = Variable("M_{atm}", .0289644, "kg/mol",
                         "Molar mass of dry air")
        p_atm = Variable("P_{atm}", "Pa", "air pressure")
        R_atm = Variable("R_{atm}", 8.31447, "J/mol/K", "air specific heating value")
        TH = 5.257386998354459 #(g*M_atm/R_atm/L_atm).value
        rho = Variable('\\rho', 'kg/m^3', 'Density of air')
        T_atm = Variable("T_{atm}", "K", "air temperature")
        h = Variable("h", "m", "Altitude")

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
        
##        t_plus_ts_approx = (T_atm + T_s).mono_approximation({T_atm: 288.15,
##                                                         T_s: T_s.value})

        with SignomialsEnabled():
            constraints = [
                # Pressure-altitude relation
                (p_atm/p_sl)**(1/5.257) == T_atm/T_sl,

                # Ideal gas law
                rho == p_atm/(R_atm/M_atm*T_atm),

                #temperature equation
##                SignomialEquality(T_sl, T_atm + L_atm*alt['h']),
                T_atm == 218*units('K'),

                #constraint on mu
##                SignomialEquality((T_atm + T_s) * mu, C_1 * T_atm**1.5),
                mu == 1.4e-5*units('kg/(m*s)'),
##                TCS([(T_atm + T_s) * mu >= C_1 * T_atm**1.5])
                ]

        #like to use a local subs here in the future
        subs = None

        Model.__init__(self, None, constraints, subs)

class Engine(Model):
    """
    place holder engine model
    """
    def __init__(self, **kwargs):
        #new variables
        W_engine = Variable('W_{engine}', 'N', 'Weight of a Single Turbofan Engine')
        
        constraints = []

        constraints.extend([
            W_engine == 1000 * units('N')
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
            TSFC == TSFC,

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

        xw     = Variable('x_w', 'm', 'Position of wing aerodynamic center')
        cwma    = Variable('\\bar{c}_w', 'm',
                          'Mean aerodynamic chord (wing)')
        Vne     = Variable('V_{ne}', 'm/s', 'Never exceed velocity')
        
        constraints = []

        constraints.extend([
            #wing weight constraint
            #based off of a raymer weight and 737 data from TASOPT output file
            (S/(124.58*units('m^2')))**.65 == W_wing/(105384.1524*units('N')),

            #compute wing span and aspect ratio, subject to a span constraint
            AR == (span**2)/S,
            #AR == 9,

            span <= span_max,

            #compute K for the aircraft
            K == (pi * e * AR)**-1,

##            xw == xw,
##            
##            Vne == Vne,
##            cwma == cwma,
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
        Lwing = Variable('L_{wing}', 'N', 'Wing Lift')
        CLaw    = Variable('C_{L_{aw}}', '-', 'Lift curve slope, wing')

        #constraints
        constraints = []

        constraints.extend([
            
            Lwing == (.5*wing['S']*state.atm['\\rho']*state['V']**2)*CL,

            #airfoil drag constraint
            TCS([Cdw**6.5 >= (1.02458748e10 * CL**15.587947404823325 * state['M']**156.86410659495155 +
                         2.85612227e-13 * CL**1.2774976672501526 * state['M']**6.2534328002723703 +
                         2.08095341e-14 * CL**0.8825277088649582 * state['M']**0.0273667615730107 +
                         1.94411925e+06 * CL**5.6547413360261691 * state['M']**146.51920742858428)]),
            TCS([Dwing >= (.5*wing['S']*state.atm['\\rho']*state['V']**2)*(Cdw + wing['K']*CL**2)]),

            CLaw == 5,
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

        lfuse   = Variable('l_{fuse}', 'm', 'Fuselage length')
        wfuse   = Variable('w_{fuse}', 'm', 'Fuselage width')

        constraints = []
        
        constraints.extend([
            #compute fuselage area for drag approximation
            A_fuse == pax_area * n_pax,

            A_fuse == lfuse * wfuse,

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
        Cmfu    = Variable('C_{m_{fuse}}', '-', 'Moment coefficient (fuselage)')
        
        #constraints
        constraints = []

        constraints.extend([
            Dfuse == Cdfuse * (.5 * fuse['A_{fuse}'] * state.atm['\\rho'] * state['V']**2),

            Cdfuse == .005,

            Cmfu == .05,
            ])

        Model.__init__(self, None, constraints)
    

class Mission(Model):
    """
    mission class, links together all subclasses
    """
    def __init__(self, subs = None, **kwargs):
        #define the number of each flight segment
        Nclimb = 2
        Ncruise = 1

        #build required submodels
        ac = Aircraft()

        #vectorize
##        with vectorize(Nclimb):
##            climb = ClimbSegment(ac)

        with vectorize(Ncruise):
            cruise = CruiseSegment(ac)

        #declare new variables
        W_ftotal = Variable('W_{f_{total}}', 'N', 'Total Fuel Weight')
        W_fclimb = Variable('W_{f_{climb}}', 'N', 'Fuel Weight Burned in Climb')
        W_fcruise = Variable('W_{f_{cruise}}', 'N', 'Fuel Weight Burned in Cruise')
        W_total = Variable('W_{total}', 'N', 'Total Aircraft Weight')
        CruiseAlt = Variable('CruiseAlt', 'ft', 'Cruise Altitude [feet]')
        ReqRng = Variable('ReqRng', 'nautical_miles', 'Required Cruise Range')

##        h = climb.state.alt['h']
##        hftClimb = climb.state.alt['hft']
##        dhft = climb.climbP['dhft']
        hftCruise = cruise.state.alt['hft']

        #make overall constraints
        constraints = []

        if Ncruise == 1:
            constraints.extend([
                #weight constraints
                TCS([ac['W_{e}'] + ac['W_{payload}'] + W_ftotal + ac['numeng'] * ac['W_{engine}'] + ac['W_{wing}'] <= W_total]),

 
                cruise.cruiseP.aircraftP['W_{start}'] == W_total,

                # similar constraint 2
                TCS([cruise.cruiseP.aircraftP['W_{start}'] >= cruise.cruiseP.aircraftP['W_{end}'] + cruise.cruiseP.aircraftP['W_{burn}']]),

                TCS([ac['W_{e}'] + ac['W_{payload}'] + ac['numeng'] * ac['W_{engine}'] + ac['W_{wing}']  <= cruise.cruiseP.aircraftP['W_{end}']]),


                TCS([W_fcruise >= sum(cruise.cruiseP['W_{burn}'])]),
                TCS([W_ftotal >=  W_fcruise]),

                #altitude constraints
                hftCruise == CruiseAlt,

                #set the range for each cruise segment, doesn't take credit for climb
                #down range disatnce covered
                cruise.cruiseP['Rng'] == ReqRng/(Ncruise),

                #set the TSFC
                cruise.cruiseP.engineP['TSFC'] == .5*units('1/hr'),
                ])
        
        if Ncruise > 1:
            constraints.extend([
                #weight constraints
                TCS([ac['W_{e}'] + ac['W_{payload}'] + W_ftotal + ac['numeng'] * ac['W_{engine}'] + ac['W_{wing}'] + ac.ht['W_{struct}'] <= W_total]),

    ##            climb.climbP.aircraftP['W_{start}'][0] == W_total,
    ##            climb.climbP.aircraftP['W_{end}'][-1] == cruise.cruiseP.aircraftP['W_{start}'][0],
                cruise.cruiseP.aircraftP['W_{start}'][0] == W_total,

                # similar constraint 1
    ##            TCS([climb.climbP.aircraftP['W_{start}'] >= climb.climbP.aircraftP['W_{end}'] + climb.climbP.aircraftP['W_{burn}']]),
                # similar constraint 2
                TCS([cruise.cruiseP.aircraftP['W_{start}'] >= cruise.cruiseP.aircraftP['W_{end}'] + cruise.cruiseP.aircraftP['W_{burn}']]),

    ##            climb.climbP.aircraftP['W_{start}'][1:] == climb.climbP.aircraftP['W_{end}'][:-1],
                cruise.cruiseP.aircraftP['W_{start}'][1:] == cruise.cruiseP.aircraftP['W_{end}'][:-1],

                TCS([ac['W_{e}'] + ac['W_{payload}'] + ac['numeng'] * ac['W_{engine}'] + ac['W_{wing}'] + ac.ht['W_{struct}'] <= cruise.cruiseP.aircraftP['W_{end}'][-1]]),

    ##            TCS([W_ftotal >=  W_fclimb + W_fcruise]),
    ##            TCS([W_fclimb >= sum(climb.climbP['W_{burn}'])]),
                TCS([W_fcruise >= sum(cruise.cruiseP['W_{burn}'])]),
                TCS([W_ftotal >=  W_fclimb]),

                #altitude constraints
                hftCruise == CruiseAlt,
    ##            TCS([hftClimb[1:Ncruise] >= hftClimb[:Ncruise-1] + dhft]),
    ##            TCS([hftClimb[0] >= dhft[0]]),
    ##            hftClimb[-1] <= hftCruise,

                #compute the dh
    ##            dhft == hftCruise/Nclimb,

                #constrain the thrust
    ##            climb.climbP.engineP['thrust'] <= 2 * max(cruise.cruiseP.engineP['thrust']),

                #set the range for each cruise segment, doesn't take credit for climb
                #down range disatnce covered
                cruise.cruiseP['Rng'] == ReqRng/(Ncruise),

                #set the TSFC
    ##            climb.climbP.engineP['TSFC'] == .7*units('1/hr'),
                cruise.cruiseP.engineP['TSFC'] == .5*units('1/hr'),
                ])
        
        # Model.__init__(self, W_ftotal + s*units('N'), constraints + ac + climb + cruise, subs)
##        Model.__init__(self, W_ftotal, constraints + ac + climb + cruise, subs)
        Model.__init__(self, W_ftotal, constraints + ac + cruise, subs)

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
        sol = m.solve(solver, verbosity, **kwargs)
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

class HorizontalTailNoStruct(Model):
    """
    horiziontal tail model from Philippe's thesis
    as a performance model without the wing box

    References:
    [1] TASOPT code
    [2] http://adg.stanford.edu/aa241/stability/staticstability.html
    """
    def __init__(self, fuse, wing,**kwargs):
        self.fuse = fuse
        self.wing = wing
        
        #variables
        p       = Variable('p_{ht}', '-', 'Substituted variable = 1 + 2*taper')
        q       = Variable('q_{ht}', '-', 'Substituted variable = 1 + taper')

        
        rho0    = Variable('\\rho_0', 1.225, 'kg/m^3', 'Air density (0 ft)')
        tanLh   = Variable('\\tan(\\Lambda_{ht})', '-',
                           'tangent of horizontal tail sweep')
        taper   = Variable('\lambda_h', '-', 'Horizontal tail taper ratio')
        tau     = Variable('\\tau_h', '-',
                           'Horizontal tail thickness/chord ratio')
        wf      = Variable('w_{fuse}', 'm', 'Fuselage width')
        xcght   = Variable('x_{CG_{ht}}', 'm', 'Horizontal tail CG location')
        ymac    = Variable('y_{\\bar{c}_{ht}}', 'm',
                           'Spanwise location of mean aerodynamic chord')
        lht     = Variable('l_{ht}', 'm', 'Horizontal tail moment arm')
        ARh     = Variable('AR_h', '-', 'Horizontal tail aspect ratio')
        Vne     = Variable('V_{ne}', 'm/s', 'Never exceed velocity')
        amax    = Variable('\\alpha_{max,h}', '-', 'Max angle of attack, htail')
        e       = Variable('e_h', '-', 'Oswald efficiency factor')
        Sw      = Variable('S_w', 'm^2', 'Wing area')
        Sh      = Variable('S_h', 'm^2', 'Horizontal tail area')
        bht     = Variable('b_{ht}', 'm', 'Horizontal tail span')
        chma    = Variable('\\bar{c}_{ht}', 'm', 'Mean aerodynamic chord (ht)')
        croot   = Variable('c_{root_h}', 'm', 'Horizontal tail root chord')
        ctip    = Variable('c_{tip_h}', 'm', 'Horizontal tail tip chord')
        cwma    = Variable('\\bar{c}_w', 'm',
                           'Mean aerodynamic chord (wing)')
        Lmax    = Variable('L_{{max}_h}', 'N', 'Maximum load')
        Kf      = Variable('K_f', '-',
                           'Empirical factor for fuselage-wing interference')
        fl      = Variable(r"f(\lambda_h)d", '-',
                           'Empirical efficiency function of taper')
        SMmin   = Variable('S.M._{min}', '-', 'Minimum stability margin')
        CLhmax  = Variable('C_{L_{hmax}}', '-', 'Max lift coefficient')

        #constraints
        constraints = []

        with SignomialsEnabled():

            constraints.extend([
                # Moment arm and geometry -- same as for vtail
##                p >= 1 + 2*taper,
##                2*q >= 1 + p,
##                ymac == (bht/3)*q/p,
##                ARh == bht**2/Sh,
##                TCS([(2./3)*(1 + taper + taper**2)*croot/q >=
##                     chma]), # [SP]
##                taper == ctip/croot,
##                TCS([Sh <= bht*(croot + ctip)/2]), # [SP]

                # K_f as f(wing position) -- (fitted posynomial)
                # from from UMich AE-481 course notes Table 9.1
##                Kf >= (1.5012*(self.wing['x_w']/self.fuse['l_{fuse}'])**2 +
##                       0.538*(self.wing['x_w']/self.fuse['l_{fuse}']) +
##                       0.0331),

##                self.wing['x_w'] == .76*self.fuse['l_{fuse}'],
                # Oswald efficiency
                # Nita, Scholz,
                # "Estimating the Oswald factor from basic
                # aircraft geometrical parameters"
##                TCS([fl >= (0.0524*taper**4 - 0.15*taper**3
##                            + 0.1659*taper**2
##                            - 0.0706*taper + 0.0119)], reltol=0.2),
                # NOTE: slightly slack
##                TCS([e*(1 + fl*ARh) <= 1]),

##                taper >= 0.2, # TODO: make less arbitrary

##                Lmax == 0.5*rho0*self.wing['V_{ne}']**2*Sh*CLhmax,

##                tanLh == tanLh,
##
##                tau == tau,
##
##                lht == lht,
##
##                SMmin == SMmin,
##
##                amax == amax,
##
##                xcght == xcght,

                Sh == 10*units('m^2'),
                ])

        Model.__init__(self, None, constraints)


class HorizontalTailPerformance(Model):
    """
    Horizontal tail performance model
    """
    def __init__(self, ht, fuse, wing, fuseP, wingP, state, **kwargs):
        self.ht = ht
        self.fuse = fuse
        self.wing = wing
        self.fuseP = fuseP
        self.wingP = wingP
        
        #variables
        alpha   = Variable('\\alpha', '-', 'Horizontal tail angle of attack')
        D       = Variable('D_{ht}', 'N', 'Horizontal tail drag')
        Lh      = Variable('L_h', 'N', 'Horizontal tail downforce')
        SM      = VectorVariable(2, 'S.M.', '-', 'Stability margin')
        Rec     = Variable('Re_{c_h}', '-',
                           'Cruise Reynolds number (Horizontal tail)')
        CLah    = Variable('C_{L_{ah}}', '-', 'Lift curve slope (htail)')
        CLah0   = Variable('C_{L_{ah_0}}', '-',
                           'Isolated lift curve slope (htail)')
        CLaw    = Variable('C_{L_{aw}}', '-', 'Lift curve slope (wing)')

        
        CLh     = Variable('C_{L_h}', '-', 'Lift coefficient (htail)')
        dxlead  = Variable('\\Delta x_{{lead}_h}', 'm',
                           'Distance from CG to horizontal tail leading edge')
        dxtrail = Variable('\\Delta x_{{trail}_h}', 'm',
                           'Distance from CG to horizontal tail trailing edge')
        dxw     = VectorVariable(2, '\\Delta x_w', 'm',
                           'Distance from aerodynamic centre to CG')
        
        etaht   = Variable('\\eta_{ht}', '-', 'Tail efficiency')
        eta     = Variable('\\eta_h', '-',
                           ("Lift efficiency (diff between sectional and "
                            "actual lift)"))
        xcg     = VectorVariable(2, 'x_{CG}', 'm', 'CG location')

        Cmac    = Variable('|C_{m_{ac}}|', '-', # Absolute value of CMwing
                           'Moment coefficient about aerodynamic centre (wing)')
        CDh     = Variable('C_{D_h}', '-', 'Horizontal tail drag coefficient')
        CD0h    = Variable('C_{D_{0_h}}', '-',
                           'Horizontal tail parasitic drag coefficient')

        left = Variable('left', '-', 'debug')
        right = Variable('right', '-', 'debug')

        #cosntraints
        constraints = []

        with SignomialsEnabled():

            constraints.extend([
##                TCS([SM + dxw/self.wing['\\bar{c}_w'] + self.ht['K_f']*self.fuse['w_{fuse}']**2*self.fuse['l_{fuse}']/(self.wingP['C_{L_{aw}}']*self.wing['S']*self.wing['\\bar{c}_w'])
##                                    <= CLah*self.ht['S_h']*self.ht['l_{ht}']/(self.wingP['C_{L_{aw}}']*self.wing['S']*self.wing['\\bar{c}_w'])]),
##
##                SM >= self.ht['S.M._{min}'],

                # Trim from UMich AE-481 course notes
##                CLh*self.ht['S_h']*self.ht['l_{ht}']/(self.wing['S']*self.wing['\\bar{c}_w']) + Cmac >= self.wingP['C_{L}']*dxw/self.wing['\\bar{c}_w'] + self.fuseP['C_{m_{fuse}}'],
##                SignomialEquality(self.wingP['C_{L}']*dxw/self.wing['\\bar{c}_w'] + self.fuseP['C_{m_{fuse}}'], right),
##                Lh >= .11*units('N'),
##                Lh == 0.5*state['\\rho']*state['V']**2*self.ht['S_h']*CLh,
                Lh == 100*units('N'),

                # Moment arm and geometry -- same as for vtail
##                TCS([dxlead + self.ht['c_{root_h}'] <= dxtrail]),
##                TCS([xcg + dxtrail <= self.fuse['l_{fuse}']], reltol=0.002),
##                TCS([dxlead + self.ht['y_{\\bar{c}_{ht}}']*self.ht['\\tan(\\Lambda_{ht})'] + 0.25*self.ht['\\bar{c}_{ht}'] >= self.ht['l_{ht}']],
##                    reltol=1e-2), # [SP]

                # DATCOM formula (Mach number makes it SP)
##                TCS([(self.ht['AR_h']/eta)**2*(1+self.ht['\\tan(\\Lambda_{ht})']**2-state['M']**2) + 8*pi*self.ht['AR_h']/CLah0
##                     <= (2*pi*self.ht['AR_h']/CLah0)**2]),

                # K_f as f(wing position) -- (fitted posynomial)
                # from from UMich AE-481 course notes Table 9.1
##                TCS([self.wing['x_w'] >= xcg + dxw]),
##                self.wing['x_w'] <= xcg + dxlead,

                # Loss of tail effectiveness due to wing downwash
##                CLah + (2*self.wingP['C_{L_{aw}}']/(pi*self.wing['AR']))*etaht*CLah0 <= CLah0*etaht,
##                CLah == CLah0,
##                CLh == CLah*alpha,
##                alpha <= self.ht['\\alpha_{max,h}'],

##                alpha >= .01,
##                dxlead >= 1 * units('m'),

                # Drag
                D == 0.5*state['\\rho']*state['V']**2*self.ht['S_h']*CDh,
                CDh == .05,
                
##                TCS([CDh >= CD0h + CLh**2/(pi*self.ht['e_h']*self.ht['AR_h'])]),

                # same drag model as vtail
##                CD0h**0.125 >= 0.19*(self.ht['\\tau_h'])**0.0075 *(Rec)**0.0017
##                            + 1.83e+04*(self.ht['\\tau_h'])**3.54*(Rec)**-0.494
##                            + 0.118*(self.ht['\\tau_h'])**0.0082 *(Rec)**0.00165
##                            + 0.198*(self.ht['\\tau_h'])**0.00774*(Rec)**0.00168,
##                Rec == state['\\rho']*state['V']*self.ht['\\bar{c}_{ht}']/state['\\mu'],

##                self.ht['x_{CG_{ht}}'] >= xcg+(dxlead+dxtrail)/2,
##                self.ht['x_{CG_{ht}}'] <= self.fuse['l_{fuse}'],

                #fix later
##                dxw == [2] * units('m'),
##                Cmac == .1,
##                etaht == .9,
##                eta == .97,
##                xcg == [17] * units('m'),
                ])

        Model.__init__(self, None, constraints)

class HorizontalTail(Model):
    """
    horiziontal tail model from Philippe's thesis
    as a performance model without the wing box

    References:
    [1] TASOPT code
    [2] http://adg.stanford.edu/aa241/stability/staticstability.html
    """
    def __init__(self, fuse, wing,**kwargs):
        self.fuse = fuse
        self.wing = wing
        
        self.htns = HorizontalTailNoStruct(self.fuse, self.wing)
        self.wb = WingBox(self.htns)

        Model.__init__(self, None, self.htns + self.wb)


    def dynamic(self, fuse, wing, fuseP, wingP, state):
        """"
        creates a horizontal tail performance model
        """
        return HorizontalTailPerformance(self, fuse, wing, fuseP, wingP, state)

class WingBox(Model):
    """
    Structural model for a wing
    source: Hoburg, "Geometric Programming for Aircraft Design Optimization"

    Note - does not have a performance model
    """

    def __init__(self, surface, **kwargs):
        # Variables
        Icap    = Variable('I_{cap}', '-',
                           'Non-dim spar cap area moment of inertia')
        Mr      = Variable('M_r', 'N', 'Root moment per root chord')
        nu      = Variable('\\nu', '-',
                           'Dummy variable = $(t^2 + t + 1)/(t+1)$')
        Wcap    = Variable('W_{cap}', 'N', 'Weight of spar caps')
        Wweb    = Variable('W_{web}', 'N', 'Weight of shear web')
        Wstruct = Variable('W_{struct}', 'N', 'Structural weight')

        # Constants
        taper = Variable('taper', 0.2, '-', 'Taper ratio')
        fwadd  = Variable('f_{w,add}', 0.4, '-',
                          'Wing added weight fraction') # TASOPT code (737.tas)
        g      = Variable('g', 9.81, 'm/s^2', 'Gravitational acceleration')
        Nlift  = Variable('N_{lift}', 2.0, '-', 'Wing loading multiplier')
        rh     = Variable('r_h', 0.75, '-',
                          'Fractional wing thickness at spar web')
        rhocap = Variable('\\rho_{cap}', 2700, 'kg/m^3',
                          'Density of spar cap material')
        rhoweb = Variable('\\rho_{web}', 2700, 'kg/m^3',
                          'Density of shear web material')
        sigmax = Variable('\\sigma_{max}', 250e6, 'Pa',
                          'Allowable tensile stress')
        sigmaxshear = Variable('\\sigma_{max,shear}', 167e6, 'Pa',
                               'Allowable shear stress')
        w      = Variable('w', 0.5, '-', 'Wingbox-width-to-chord ratio')
        tcap    = Variable('t_{cap}' ,'-', 'Non-dim. spar cap thickness')
        tweb    = Variable('t_{web}', '-', 'Non-dim. shear web thickness')
        
        objective = Wstruct

##        if isinstance(surface, HorizontalTailNoStruct):
##            AR = surface['AR_h']
##            b = surface['b_{ht}']
##            S = surface['S_h']
##            p = surface['p_{ht}']
##            q = surface['q_{ht}']
##            tau = surface['\\tau_h']
##            Lmax = surface['L_{{max}_h}']

        constraints = [
                       # Upper bound on maximum thickness
##                       tau <= 0.15,

                       # Root moment calculation (see Hoburg 2014)
                       # Depends on a given load the wing must support, Lmax
                       # Assumes lift per unit span proportional to local chord
##                       Mr >= Lmax*AR*p/24,

                       # Root stiffness (see Hoburg 2014)
                       # Assumes rh = 0.75, so that rms box height = ~0.92*tmax
##                       0.92*w*tau*tcap**2 + Icap <= 0.92**2/2*w*tau**2*tcap,

                       # Stress limit
                       # Assumes bending stress carried by caps (Icap >> Iweb)
##                       8 >= Nlift*Mr*AR*q**2*tau/(S*Icap*sigmax),

                       # Shear web sizing
                       # Assumes all shear loads are carried by web and rh=0.75
##                       12 >= AR*Lmax*Nlift*q**2/(tau*S*tweb*sigmaxshear),

                       # Posynomial approximation of nu=(1+lam+lam^2)/(1+lam^2)
##                       nu**3.94 >= 0.86*p**(-2.38)+ 0.14*p**0.56,

                       # Weight of spar caps and shear webs
##                       Wcap >= 8*rhocap*g*w*tcap*S**1.5*nu/(3*AR**0.5),
##                       Wweb >= 8*rhoweb*g*rh*tau*tweb*S**1.5*nu/(3*AR**0.5),

                       # Total wing weight using an additional weight fraction
                       Wstruct == 1000*units('N'),#(1 + fwadd)*(Wweb + Wcap),
                       ]
        
        Model.__init__(self, objective, constraints, **kwargs)

if __name__ == '__main__':
    substitutions = {      
##            'V_{stall}': 120,
            'ReqRng': 500, #('sweep', np.linspace(500,2000,4)),
            'CruiseAlt': 30000, #('sweep', np.linspace(20000,40000,4)),
            'numeng': 2,
##            'W_{Load_max}': 6664,
            'W_{pax}': 91 * 9.81,
            'n_{pax}': 150,
            'pax_{area}': 1,
##            'C_{D_{fuse}}': .005, #assumes flat plate turbulent flow, from wikipedia
            'e': .9,
            'b_{max}': 35,

##             'S.M._{min}': 0.05,
##             'V_{ne}': 144,
##             'C_{L_{hmax}}': 2.5,

##             '\\alpha_{max,h}': .3,#0.1, # (6 deg)
##             '\\bar{c}_w': 5,
##             '\\rho_0': 1.225,
##             '\\tan(\\Lambda_{ht})': tan(30*pi/180),
             'w_{fuse}': 6,
            }
           
    m = Mission(substitutions)
##    sol = m.localsolve(solver='mosek', verbosity = 4)
    bounds, sol = m.determine_unbounded_variables(m, solver="mosek",verbosity=4, iteration_limit=100)
 
    m = Mission(substitutions)
##    solRsweep = m.localsolve(solver='mosek', verbosity = 4)

##    plt.plot(solRsweep('ReqRng'), solRsweep('W_{f_{total}}'), '-r')
##    plt.xlabel('Mission Range [nm]')
##    plt.ylabel('Total Fuel Burn [N]')
##    plt.title('Fuel Burn vs Range')
##    plt.show()


    m = Mission(substitutions)
##    solAltsweep = m.localsolve(solver='mosek', verbosity = 4)

##    plt.plot(solAltsweep('CruiseAlt'), solAltsweep('W_{f_{total}}'), '-r')
##    plt.xlabel('Cruise Alt [ft]')
##    plt.ylabel('Total Fuel Burn [N]')
##    plt.title('Fuel Burn vs Range')
##    plt.show()

