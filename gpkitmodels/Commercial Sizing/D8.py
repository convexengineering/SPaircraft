"""Fully integrated D8 model"""
from numpy import pi, tan, cos
import numpy as np
from gpkit import Variable, Model, units, SignomialsEnabled, vectorize
from gpkit.constraints.sigeq import SignomialEqualityConstraint as SignomialEquality
from gpkit.tools import te_exp_minus1
from gpkit.constraints.tight import TightConstraintSet as TCS
import matplotlib.pyplot as plt

#only needed for the local bounded debugging tool
from collections import defaultdict
from gpkit.small_scripts import mag

#import aircraft subsystems
from stand_alone_simple_profile import FlightState, Altitude, Atmosphere
from VT_simple_profile import VerticalTail, VerticalTailPerformance
from Wing_simple_performance import Wing, WingPerformance
from D8_integration import Fuselage, Engine, EnginePerformance
from CFP_Fuselage_Performance import FuselagePerformance
#from CFP_Fuselage_Performance import Fuselage
#from D8_integration import Engine, EnginePerformance
#from CFP_Fuselage_Performance import Fuselage, FuselagePerformance
#from HT_Simple_Profile_Performance import HorizontalTail, HorizontalTailPerformance

#set up models
class D8(Model):
    "Aircraft class"
    def __init__(self, **kwargs):
        #create submodels
        self.fuse = Fuselage()
        self.wing = Wing()
        self.engine = Engine()
        self.VT = VerticalTail(self.fuse, self.engine)
        #self.HT = HorizontalTail(self.fuse, self.wing)

        #variable definitions
        numeng = Variable('numeng', '-', 'Number of Engines')

        constraints = []

        constraints.extend([
            numeng == numeng, #need numeng in the model
            ])

        self.components = [self.fuse, self.wing, self.engine, self.VT]#, self.HT]

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

class D8P(Model):
    """
    aircraft performance models superclass, contains constraints true for
    all flight segments
    """
    def  __init__(self, aircraft, state, **kwargs):
        #make submodels
        self.aircraft = D8()
        self.wingP = aircraft.wing.dynamic(state)
        self.fuseP = aircraft.fuse.dynamic(state)
        self.engineP = aircraft.engine.dynamic(state)
        self.VTP = aircraft.VT.dynamic(aircraft.fuse, self.fuseP, state)
        #self.HTP = aircraft.HT.dynamic(self.aircraft.fuse, self.aircraft.wing, self.fuseP, self.wingP, state)

        self.Pmodels = [self.wingP, self.fuseP, self.engineP, self.VTP]#, self.HTP]

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

        constraints.extend([
            #speed must be greater than stall speed
            state['V'] >= Vstall,


            #Figure out how to delete
            Vstall == 120*units('kts'),
            WLoadmax == 6664 * units('N/m^2'),

            #compute the drag
            TCS([D >= self.wingP['D_{wing}'] + self.fuseP['D_{fuse}'] + self.VTP['D_{vt}']]), #+ self.HTP['D_{ht}']

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
            self.wingP['L_w'] == W_avg,
            ])

        Model.__init__(self, None, [self.Pmodels + constraints], **kwargs)

class ClimbP(Model):
    """
    Climb constraints
    """
    def __init__(self, aircraft, state, **kwargs):
        #submodels
        self.aircraft = aircraft
        self.aircraftP = D8P(aircraft, state)
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
        self.aircraftP = D8P(aircraft, state)
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

class Mission(Model):
    """
    mission class, links together all subclasses
    """
    def __init__(self, ac, substitutions = None, **kwargs):
        #define the number of each flight segment
        Nclimb = 2
        Ncruise = 2

        #vectorize
        with vectorize(Nclimb):
            climb = ClimbSegment(ac)

        with vectorize(Ncruise):
            cruise= CruiseSegment(ac)

        #declare new variables
        W_ftotal = Variable('W_{f_{total}}', 'N', 'Total Fuel Weight')
        W_fclimb = Variable('W_{f_{climb}}', 'N', 'Fuel Weight Burned in Climb')
        W_fcruise = Variable('W_{f_{cruise}}', 'N', 'Fuel Weight Burned in Cruise')
        W_total = Variable('W_{total}', 'N', 'Total Aircraft Weight')
        W_dry = Variable('W_{dry}', 'N', 'Total Aircraft Weight with No Fuel')
        CruiseAlt = Variable('CruiseAlt', 'ft', 'Cruise Altitude [feet]')
        ReqRng = Variable('ReqRng', 'nautical_miles', 'Required Cruise Range')

        h = climb.state.alt['h']
        hftClimb = climb.state.alt['hft']
        dhft = climb.climbP['dhft']
        hftCruise = cruise.state.alt['hft']

        #make overall constraints
        constraints = []

        constraints.extend([
            #weight constraints
            TCS([ac['W_{e}'] + ac['W_{payload}'] + ac['numeng'] * ac['W_{engine}'] + ac.wing['W_{struct}'] + ac.VT['W_{struct}']<= W_dry]), # + ac.HT['W_{struct}']
            
            TCS([W_dry + W_ftotal <= W_total]),

            climb.climbP.aircraftP['W_{start}'][0] == W_total,
            climb.climbP.aircraftP['W_{end}'][-1] == cruise.cruiseP.aircraftP['W_{start}'][0],

            # similar constraint #1
            TCS([climb.climbP.aircraftP['W_{start}'] >= climb.climbP.aircraftP['W_{end}'] + climb.climbP.aircraftP['W_{burn}']]),
            # similar constraint #2
            TCS([cruise.cruiseP.aircraftP['W_{start}'] >= cruise.cruiseP.aircraftP['W_{end}'] + cruise.cruiseP.aircraftP['W_{burn}']]),

            climb.climbP.aircraftP['W_{start}'][1:] == climb.climbP.aircraftP['W_{end}'][:-1],
            cruise.cruiseP.aircraftP['W_{start}'][1:] == cruise.cruiseP.aircraftP['W_{end}'][:-1],

            TCS([W_dry <= cruise.cruiseP.aircraftP['W_{end}'][-1]]),

            TCS([W_ftotal >=  W_fclimb + W_fcruise]),
            TCS([W_fclimb >= sum(climb.climbP['W_{burn}'])]),
            TCS([W_fcruise >= sum(cruise.cruiseP['W_{burn}'])]),

            #altitude constraints
            hftCruise == CruiseAlt,
            TCS([hftClimb[1:Ncruise] >= hftClimb[:Ncruise-1] + dhft]),
            TCS([hftClimb[0] >= dhft[0]]),
            hftClimb[-1] <= hftCruise,

            #compute the dh
            dhft == hftCruise/Nclimb,

            #constrain the thrust
            climb.climbP.engineP['thrust'] <= 2 * max(cruise.cruiseP.engineP['thrust']),

            #set the range for each cruise segment, doesn't take credit for climb
            #down range disatnce covered
            cruise.cruiseP['Rng'] == ReqRng/(Ncruise),

            #set the TSFC
            climb.climbP.engineP['TSFC'] == .7*units('1/hr'),
            cruise.cruiseP.engineP['TSFC'] == .5*units('1/hr'),

            #VT constriants
            ac.VT['T_e'] == climb.climbP.engineP['thrust'][0],

            #wing constraints
            ac.wing['W_{fuel_{wing}}'] == W_ftotal,
            ])
        
        # Model.__init__(self, W_ftotal + s*units('N'), constraints + ac + climb + cruise, subs)
        Model.__init__(self, W_ftotal, constraints + ac + climb + cruise, substitutions)

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
            if distance_below <= 3:  # arbitrary thresholdprint 
                out["value near lower bound"].append(varkey)
            elif distance_above <= 3:  # arbitrary threshold
                out["value near upper bound"].append(varkey)
        return out, solhold


if __name__ == '__main__':
    #build a D8 model
    ac = D8()

    wing_sweep = 30 #[deg]
    
    substitutions = {      
##            'V_{stall}': 120,
            'ReqRng': 500, #('sweep', np.linspace(500,2000,4)),
            'CruiseAlt': 30000, #('sweep', np.linspace(20000,40000,4)),
            'numeng': 1,
##            'W_{Load_max}': 6664,
            'W_{pax}': 91 * 9.81,
            'n_{pax}': 150,
            'pax_{area}': 1,
##            'C_{D_{fuse}}': .005, #assumes flat plate turbulent flow, from wikipedia

            #VT subs
           'C_{D_{wm}}': 0.5, # [2]
           'C_{L_{vmax}}': 2.6, # [2]
           'V_1': 70,
##           'V_{ne}': 144, # [2]
           '\\rho_{TO}': 1.225,
           '\\tan(\\Lambda_{vt})': np.tan(40*np.pi/180),
##           'c_{l_{vt}}': 0.5, # [2]
           'c_{l_{vtEO}}': 0.5,
           'A_2': np.pi*(.5*1.75)**2, # [1]
           'e_v': 0.8,
##           'x_{CG}': 18,
           'y_{eng}': 4.83, # [3]

            #wing subs
            'C_{L_{wmax}}': 2.5,
            '\\alpha_{max,w}': 0.1, # (6 deg)
            '\\cos(\\Lambda)': cos(wing_sweep*pi/180),
            '\\eta': 0.97,
##            '\\rho_0': 1.225,
            '\\rho_{fuel}': 817, # Kerosene [TASOPT]
            '\\tan(\\Lambda)': tan(wing_sweep*pi/180),

#             #HT subs
#              'S.M._{min}': 0.05,
# ##             'V_{ne}': 144,
#              'C_{L_{hmax}}': 2.5,

             # '\\alpha_{max,h}': .3,#0.1, # (6 deg)
##             '\\bar{c}_w': 5,
##             '\\rho_0': 1.225,
             # '\\tan(\\Lambda_{ht})': tan(30*pi/180),
##             'w_{fuse}': 6,
            }
           
    m = Mission(ac, substitutions)
    sol = m.localsolve(solver='mosek', verbosity = 2)
    # bounds, sol = m.determine_unbounded_variables(m, solver="mosek",verbosity=2, iteration_limit=100)
