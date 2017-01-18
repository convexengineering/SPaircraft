"""D8 horizontal tail model linked with a simple flight profile"""
from numpy import pi, tan
import numpy as np
from gpkit import Variable, Model, units, SignomialsEnabled, Vectorize, SignomialEquality
from gpkit.tools import te_exp_minus1
from gpkit.constraints.tight import Tight as TCS
import matplotlib.pyplot as plt
#only needed for the local bounded debugging tool
from collections import defaultdict
from gpkit.small_scripts import mag

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
    def setup(self):
        #create submodels
        self.fuse = Fuselage()
        self.wing = Wing()
        self.engine = Engine()
        self.HT = HorizontalTail()

        #variable definitions
        numeng = Variable('numeng', '-', 'Number of Engines')
        Vne     = Variable('V_{ne}', 144, 'm/s', 'Never exceed velocity')
        rho0    = Variable('\\rho_0', 1.225, 'kg/m^3', 'Air density (0 ft)')

        SMmin = Variable('SM_{min}', '-', 'Minimum Static Margin')
        dxCG = Variable('\\Delta x_{CG}', 'm', 'Max CG Travel Range')


        constraints = []

        constraints.extend([
            numeng == numeng, #need numeng in the model
            ])

        self.components = [self.fuse, self.wing, self.engine, self.HT]

        return self.components, constraints
        
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
    def  setup(self, aircraft, state):
        #make submodels
        self.aircraft = aircraft
        self.wingP = aircraft.wing.dynamic(state)
        self.fuseP = aircraft.fuse.dynamic(state)
        self.engineP = aircraft.engine.dynamic(state)
        self.HTP = aircraft.HT.dynamic(aircraft.fuse, aircraft.wing, state)

        self.Pmodels = [self.wingP, self.fuseP, self.engineP, self.HTP]

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

        xAC = Variable('x_{AC}', 'm', 'Aerodynamic Center Location')
        xCG     = Variable('x_{CG}', 'm', 'CG location')


        constraints = []
        with SignomialsEnabled():
            constraints.extend([
                #speed must be greater than stall speed
                state['V'] >= Vstall,


                #Figure out how to delete
                Vstall == 120*units('kts'),
                WLoadmax == 6664 * units('N/m^2'),

                #compute the drag
                TCS([D >= self.wingP['D_{wing}'] + self.fuseP['D_{fuse}'] + self.HTP['D_{ht}']]),

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
                self.wingP['L_{wing}'] >= W_avg,
                 ])

        return self.Pmodels, constraints

class ClimbP(Model):
    """
    Climb constraints
    """
    def setup(self, aircraft, state):
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

        return constraints, self.aircraftP

class CruiseP(Model):
    """
    Cruise constraints
    """
    def setup(self, aircraft, state):
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
    def setup(self, aircraft):
        self.state = FlightState()
        self.cruiseP = aircraft.cruise_dynamic(self.state)

        return self.state, self.cruiseP

class ClimbSegment(Model):
    """
    Combines a flight state and aircrat to form a cruise flight segment
    """
    def setup(self, aircraft):
        self.state = FlightState()
        self.climbP = aircraft.climb_dynamic(self.state)

        return self.state, self.climbP

class FlightState(Model):
    """
    creates atm model for each flight segment, has variables
    such as veloicty and altitude
    """
    def setup(self):
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
        return constraints, self.atm, self.alt

class Altitude(Model):
    """
    holds the altitdue variable
    """
    def setup(self):
        #define altitude variables
        h = Variable('h', 'm', 'Segment Altitude [meters]')
        hft = Variable('hft', 'feet', 'Segment Altitude [feet]')

        constraints = []

        constraints.extend([
            h == hft, #convert the units on altitude
            ])

        return constraints
            

class Atmosphere(Model):
    def setup(self, alt):
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

        """
        Dynamic viscosity (mu) as a function of temperature
        References:
        http://www-mdp.eng.cam.aircraft.uk/web/library/enginfo/aerothermal_dvd_only/aero/
            atmos/atmos.HTml
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
                SignomialEquality(T_sl, T_atm + L_atm*alt['h']),

                #constraint on mu
                SignomialEquality((T_atm + T_s) * mu, C_1 * T_atm**1.5),
                ]

        return constraints

class Engine(Model):
    """
    place holder engine model
    """
    def setup(self):
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
    def setup(self, engine, state):
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
    def setup(self):
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

        mac    = Variable('mac', 'm',
                          'Mean aerodynamic chord (wing)')
        cmw = Variable('c_{m_{w}}', '-', 'Wing Pitching Moment Coefficient')

        CLmax = Variable('C_{L_{max}}', '-', 'Max Wing Lift Coefficient')

        xw     = Variable('x_w', 'm', 'Position of wing aerodynamic center')
        
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
            mac == mac,
            cmw == 1.0,
            CLmax == CLmax,

            xw == xw,
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
    def setup(self, wing, state):
        #new variables
        CL= Variable('C_{L}', '-', 'Lift Coefficient')
        Cdw = Variable('C_{d_w}', '-', 'Cd for a NC130 Airfoil at Re=2e7')
        Dwing = Variable('D_{wing}', 'N', 'Total Wing Drag')
        Lwing = Variable('L_{wing}', 'N', 'Wing Lift')

        CLaw    = Variable('C_{L_{aw}}', '-', 'Lift curve slope, wing')

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
            
            CLaw == 5,
            ])

        return constraints

class Fuselage(Model):
    """
    place holder fuselage model
    """
    def setup(self):
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
    def setup(self, fuse, state):
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

        return constraints
    

class Mission(Model):
    """
    mission class, links together all subclasses
    """
    def setup(self, aircraft):
        #define the number of each flight segment
        Nclimb = 2
        Ncruise = 2

        #Vectorize
        with Vectorize(Nclimb):
            climb = ClimbSegment(aircraft)

        with Vectorize(Ncruise):
            cruise = CruiseSegment(aircraft)

        #declare new variables
        W_ftotal = Variable('W_{f_{total}}', 'N', 'Total Fuel Weight')
        W_fclimb = Variable('W_{f_{climb}}', 'N', 'Fuel Weight Burned in Climb')
        W_fcruise = Variable('W_{f_{cruise}}', 'N', 'Fuel Weight Burned in Cruise')
        W_total = Variable('W_{total}', 'N', 'Total Aircraft Weight')
        CruiseAlt = Variable('CruiseAlt', 'ft', 'Cruise Altitude [feet]')
        ReqRng = Variable('ReqRng', 'nautical_miles', 'Required Cruise Range')

        h = climb['h']
        hftClimb = climb['hft']
        dhft = climb['dhft']
        hftCruise = cruise['hft']

        #make overall constraints
        constraints = []

        constraints.extend([
            #weight constraints
            TCS([aircraft['W_{e}'] + aircraft['W_{payload}'] + W_ftotal + aircraft['numeng'] * aircraft['W_{engine}'] + aircraft['W_{wing}'] + aircraft.HT['W_{struct}'] <= W_total]),

            climb['W_{start}'][0] == W_total,
            climb['W_{end}'][-1] == cruise['W_{start}'][0],

            # similar constraint 1
            TCS([climb['W_{start}'] >= climb['W_{end}'] + climb['W_{burn}']]),
            # similar constraint 2
            TCS([cruise['W_{start}'] >= cruise['W_{end}'] + cruise['W_{burn}']]),

            climb['W_{start}'][1:] == climb['W_{end}'][:-1],
            cruise['W_{start}'][1:] == cruise['W_{end}'][:-1],

            TCS([aircraft['W_{e}'] + aircraft['W_{payload}'] + aircraft['numeng'] * aircraft['W_{engine}'] + aircraft['W_{wing}'] + aircraft.HT['W_{struct}'] <= cruise['W_{end}'][-1]]),

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

            # climb['C_{L_h}'] == 2*3.14*climb['\\alpha_{ht}],
            # cruise['C_{L_h}'] == 2*3.14*cruise['\\alpha_{ht}],
            ])
        
        #Horizontal Tail Constraints
        with SignomialsEnabled():
            constraints.extend([

                # Trim condition for each flight segment
                TCS([cruise['x_{AC}']/aircraft.wing['mac'] <= aircraft.wing['c_{m_{w}}']/cruise['C_{L}'] + \
                     cruise['x_{CG}']/aircraft.wing['mac'] + aircraft.HT['V_{h}']*(cruise['C_{L_h}']/cruise['C_{L}'])]),
                TCS([climb['x_{AC}']/aircraft.wing['mac'] <= aircraft.wing['c_{m_{w}}']/climb['C_{L}'] + \
                     climb['x_{CG}']/aircraft.wing['mac'] + aircraft.HT['V_{h}']*(climb['C_{L_h}']/climb['C_{L}'])]),


                aircraft.HT['L_{{max}_h}'] == 0.5*aircraft['\\rho_0']*aircraft['V_{ne}']**2*aircraft.HT['S_h']*aircraft.HT['C_{L_{hmax}}'],
                #compute mrat, is a signomial equality
                SignomialEquality(aircraft.HT['m_{ratio}']*(1+2/aircraft.wing['AR']), 1 + 2/aircraft.HT['AR_h']),

                #tail volume coefficient
                aircraft.HT['V_{h}'] == aircraft.HT['S_h']*aircraft.HT['x_{CG_{ht}}']/(aircraft.wing['S']*aircraft.wing['mac']),

                #enforce max tail location is the end of the fuselage
                aircraft.HT['x_{CG_{ht}}'] <= aircraft.fuse['l_{fuse}'],
                aircraft.HT['l_{ht}'] >= aircraft.HT['x_{CG_{ht}}'] - cruise['x_{CG}'],
                aircraft.HT['l_{ht}'] >= aircraft.HT['x_{CG_{ht}}'] - climb['x_{CG}'],

                #Stability constraint, is a signomial
                TCS([aircraft['SM_{min}'] + aircraft['\\Delta x_{CG}']/aircraft.wing['mac'] <= aircraft.HT['V_{h}']*aircraft.HT['m_{ratio}'] + aircraft.wing['c_{m_{w}}']/aircraft.wing['C_{L_{max}}'] + aircraft.HT['V_{h}']*aircraft.HT['C_{L_{hmax}}']/aircraft.wing['C_{L_{max}}']]),

                TCS([aircraft.wing['x_w'] >= cruise['x_{CG}'] + cruise['\\Delta x_w']]),
                TCS([aircraft.wing['x_w'] >= climb['x_{CG}'] + climb['\\Delta x_w']]),


                TCS([cruise['x_{CG}'] + cruise['\\Delta x_{{trail}_h}'] <= aircraft.fuse['l_{fuse}']], reltol=0.002),
                TCS([climb['x_{CG}'] + climb['\\Delta x_{{trail}_h}'] <= aircraft.fuse['l_{fuse}']], reltol=0.002),

                #compute the aerodynamic center location
                TCS([climb['x_{AC}'] <= climb['x_{CG}'] + climb['\\Delta x_w'] ]),
                TCS([cruise['x_{AC}'] <= cruise['x_{CG}'] + cruise['\\Delta x_w'] ]),

##                SignomialEquality(cruise['x_{ac}'],xcg + cruise['\\Delta x_w'] ),
##                SignomialEquality(climb['x_{ac}'],xcg + climb['\\Delta x_w'] ),
                TCS([aircraft.HT['x_{CG_{ht}}'] >= climb['x_{CG}'] + (climb['\\Delta x_{{lead}_h}']+climb['\\Delta x_{{trail}_h}'])/2]),
                TCS([aircraft.HT['x_{CG_{ht}}'] >= cruise['x_{CG}'] + (cruise['\\Delta x_{{lead}_h}']+cruise['\\Delta x_{{trail}_h}'])/2]),
                #---------------------------------------------------------#

                # Substitutions for xCG and xAC
                cruise['x_{CG}'] == 15*units('m'),
                climb['x_{CG}'] == 15*units('m'),
                cruise['x_{AC}'] == aircraft.wing['x_w'],
                climb['x_{AC}'] == aircraft.wing['x_w'],

                #compute the HT chord at its attachment point to the VT
                (aircraft.HT['b_{ht}']/aircraft.fuse['w_{fuse}'])*aircraft.HT['\lambda_h']*aircraft.HT['c_{root_h}'] == aircraft.HT['c_{attach}']
                                              
                ])

        return climb, cruise, constraints

    def bound_all_variables(self, model, eps=1e-30, lower=None, upper=None):
        "Returns model with additional constraints bounding all free variables"
        lb = lower if lower else eps
        ub = upper if upper else 1/eps
        constraints = []
        freevks = tuple(vk for vk in model.varkeys if "value" not in vk.descr)
        for varkey in freevks:
            units = varkey.descr.get("units", 1)
            varub = Variable('varub', ub, units)
            varlb = Variable('varls', lb, units)
            constraints.append([varub >= Variable(**varkey.descr),
                                Variable(**varkey.descr) >= varlb])
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


class HorizontalTailNoStruct(Model):
    """
    horiziontal tail model from Philippe's thesis
    as a performance model without the wing box

    References:
    [1] TASOPT code
    [2] http://adg.stanford.edu/aa241/stability/staticstability.HTml

    This model does not include the effects of wing downwash on tail
    effectiveness.
    """
    def setup(self):
        
        #variables
        p       = Variable('p_{ht}', '-', 'Substituted variable = 1 + 2*taper')
        q       = Variable('q_{ht}', '-', 'Substituted variable = 1 + taper')

        
        tanLh   = Variable('\\tan(\\Lambda_{ht})', '-',
                           'tangent of horizontal tail sweep')
        taper   = Variable('\lambda_h', '-', 'Horizontal tail taper ratio')
        tau     = Variable('\\tau_h', '-',
                           'Horizontal tail thickness/chord ratio')
        xcght   = Variable('x_{CG_{ht}}', 'm', 'Horizontal tail CG location')
        ymac    = Variable('y_{\\bar{c}_{ht}}', 'm',
                           'Spanwise location of mean aerodynamic chord')
        lht     = Variable('l_{ht}', 'm', 'Horizontal tail moment arm')
        ARh     = Variable('AR_h', '-', 'Horizontal tail aspect ratio')
        amax    = Variable('\\alpha_{max,h}', '-', 'Max angle of attack, htail')
        e       = Variable('e_h', '-', 'Oswald efficiency factor')
        Sh      = Variable('S_h', 'm^2', 'Horizontal tail area')
        bht     = Variable('b_{ht}', 'm', 'Horizontal tail span')
        chma    = Variable('\\bar{c}_{ht}', 'm', 'Mean aerodynamic chord (ht)')
        croot   = Variable('c_{root_h}', 'm', 'Horizontal tail root chord')
        ctip    = Variable('c_{tip_h}', 'm', 'Horizontal tail tip chord')
        Lmax    = Variable('L_{{max}_h}', 'N', 'Maximum load')
        # Kf      = Variable('K_f', '-',
        #                    'Empirical factor for fuselage-wing interference')
        fl      = Variable(r"f(\lambda_h)", '-',
                           'Empirical efficiency function of taper')
        CLhmax  = Variable('C_{L_{hmax}}', '-', 'Max lift coefficient')

        #new variables
        Vh = Variable('V_{h}', '-', 'Horizontal Tail Volume Coefficient')
        mrat = Variable('m_{ratio}', '-', 'Wing to Tail Lift Slope Ratio')
        #variable just for the D8
        cattach = Variable('c_{attach}', 'm', 'HT Chord Where it is Mountded to the VT')

        #constraints
        constraints = []

        with SignomialsEnabled():

            constraints.extend([
                # Moment arm and geometry -- same as for vtail
                p >= 1 + 2*taper,
                2*q >= 1 + p,
                ymac == (bht/3)*q/p,
                TCS([(2./3)*(1 + taper + taper**2)*croot/q >=
                     chma]), # [SP]
                taper == ctip/croot,
                SignomialEquality(Sh, bht*(croot + ctip)/2),
                # TCS([Sh <= bht*(croot + ctip)/2]), # [SP]

                # Oswald efficiency
                # Nita, Scholz,
                # "Estimating the Oswald factor from basic
                # aircraft geometrical parameters"
                TCS([fl >= (0.0524*taper**4 - 0.15*taper**3
                            + 0.1659*taper**2
                            - 0.0706*taper + 0.0119)], reltol=0.2),
                # NOTE: slightly slack
                TCS([e*(1 + fl*ARh) <= 1]),

                ARh == bht**2/Sh,
                
                taper >= 0.2, # TODO: make less arbitrary
                taper <= 1,

##                Sh >= .1*units('m^2'),
                ])

        return constraints

class HorizontalTailPerformance(Model):
    """
    Horizontal tail performance model
    """
    def setup(self, ht, fuse, wing, state):
        self.HT = ht
        self.fuse = fuse
        self.wing = wing
        
        #variables

        D       = Variable('D_{ht}', 'N', 'Horizontal tail drag')
        Lh      = Variable('L_h', 'N', 'Horizontal tail downforce')
        Rec     = Variable('Re_{c_h}', '-',
                           'Cruise Reynolds number (Horizontal tail)')
        CLah    = Variable('C_{L_{ah}}', '-', 'Lift curve slope (htail)')
        # CLah0   = Variable('C_{L_{ah_0}}', '-',
        #                    'Isolated lift curve slope (htail)')
        # CLaw    = Variable('C_{L_{aw}}', '-', 'Lift curve slope (wing)')

        
        CLh     = Variable('C_{L_h}', '-', 'Lift coefficient (htail)')
        dxlead  = Variable('\\Delta x_{{lead}_h}', 'm',
                           'Distance from CG to horizontal tail leading edge')
        dxtrail = Variable('\\Delta x_{{trail}_h}', 'm',
                           'Distance from CG to horizontal tail trailing edge')
        dxw     = Variable('\\Delta x_w', 'm',
                           'Distance from aerodynamic centre to CG')
        
        # etaht   = Variable('\\eta_{ht}', '-', 'Tail efficiency')
        # eta     = Variable('\\eta_h', '-',
        #                    ("Lift efficiency (diff between sectional and "
        #                     "actual lift)"))
        CDh     = Variable('C_{D_h}', '-', 'Horizontal tail drag coefficient')
        CD0h    = Variable('C_{D_{0_h}}', '-',
                           'Horizontal tail parasitic drag coefficient')

        alphah   = Variable('\\alpha_{ht}', '-', 'Horizontal tail angle of attack')

        constraints = []

        with SignomialsEnabled():
           
            constraints.extend([
                Lh == 0.5*state['\\rho']*state['V']**2*self.HT['S_h']*CLh,

                # Angle of attack and lift slope constraints
                CLh <= 1.1*CLah*alphah,
                CLh >= 0.9*CLah*alphah,
                alphah <= self.HT['\\alpha_{max,h}'],

                # Moment arm and geometry -- same as for vtail
##                dxlead >= self.wing['x_w'] + 1.5*units('m'),
                TCS([dxlead + self.HT['c_{root_h}'] <= dxtrail]),
                TCS([dxlead + self.HT['y_{\\bar{c}_{ht}}']*self.HT['\\tan(\\Lambda_{ht})'] + 0.25*self.HT['\\bar{c}_{ht}'] >= self.HT['l_{ht}']],
                    reltol=1e-2), # [SP]               
                dxtrail <= self.fuse['l_{fuse}'],

                # Currently using TAT to approximate
                CLah == 2*3.14,

                # Drag
                D == 0.5*state['\\rho']*state['V']**2*self.HT['S_h']*CDh,
                CDh >= CD0h + CLh**2/(pi*self.HT['e_h']*self.HT['AR_h']),

                # Same drag model as VerticalTail
                CD0h**0.125 >= 0.19*(self.HT['\\tau_h'])**0.0075 *(Rec)**0.0017
                            + 1.83e+04*(self.HT['\\tau_h'])**3.54*(Rec)**-0.494
                            + 0.118*(self.HT['\\tau_h'])**0.0082 *(Rec)**0.00165
                            + 0.198*(self.HT['\\tau_h'])**0.00774*(Rec)**0.00168,
                Rec == state['\\rho']*state['V']*self.HT['\\bar{c}_{ht}']/state['\\mu'],
                ])

        return constraints

class HorizontalTail(Model):
    """
    horiziontal tail model from Philippe's thesis
    as a performance model without the wing box

    References:
    [1] TASOPT code
    [2] http://adg.stanford.edu/aa241/stability/staticstability.HTml
    """
    def setup(self):
        self.HTns = HorizontalTailNoStruct()
        self.wb = WingBox(self.HTns)

        return self.HTns, self.wb


    def dynamic(self, fuse, wing, state):
        """"
        creates a horizontal tail performance model
        """
        return HorizontalTailPerformance(self, fuse, wing, state)
    
class WingBox(Model):
    """
    Structural model for a wing
    source: Hoburg, "Geometric Programming for Aircraft Design Optimization"

    Note - does not have a performance model
    """

    def setup(self, surface):
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

        if isinstance(surface, HorizontalTailNoStruct):
            AR = surface['AR_h']
            b = surface['b_{ht}']
            S = surface['S_h']
            p = surface['p_{ht}']
            q = surface['q_{ht}']
            tau = surface['\\tau_h']
            Lmax = surface['L_{{max}_h}']

        constraints = [
                       # Upper bound on maximum thickness
                       tau <= 0.15,

                       # Root moment calculation (see Hoburg 2014)
                       # Depends on a given load the wing must support, Lmax
                       # Assumes lift per unit span proportional to local chord
                       Mr >= Lmax*AR*p/24,

                       # Root stiffness (see Hoburg 2014)
                       # Assumes rh = 0.75, so that rms box height = ~0.92*tmax
                       0.92*w*tau*tcap**2 + Icap <= 0.92**2/2*w*tau**2*tcap,

                       # Stress limit
                       # Assumes bending stress carried by caps (Icap >> Iweb)
                       8 >= Nlift*Mr*AR*q**2*tau/(S*Icap*sigmax),

                       # Shear web sizing
                       # Assumes all shear loads are carried by web and rh=0.75
                       12 >= AR*Lmax*Nlift*q**2/(tau*S*tweb*sigmaxshear),

                       # Posynomial approximation of nu=(1+lam+lam^2)/(1+lam^2)
                       nu**3.94 >= 0.86*p**(-2.38)+ 0.14*p**0.56,

                       # Weight of spar caps and shear webs
                       Wcap >= 8*rhocap*g*w*tcap*S**1.5*nu/(3*AR**0.5),
                       Wweb >= 8*rhoweb*g*rh*tau*tweb*S**1.5*nu/(3*AR**0.5),

                       # Total wing weight using an additional weight fraction
                       Wstruct >= (1 + fwadd)*(Wweb + Wcap),
                       ]
        
        return constraints

if __name__ == '__main__':
    plot = True
    
    #build required submodels
    aircraft = Aircraft()
        
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
            'e': .9,
            'b_{max}': 35,
             'C_{L_{hmax}}': 2.5,
             '\\tan(\\Lambda_{ht})': tan(30*pi/180),
             'w_{fuse}': 6,

##            'l_{fuse}': 30,
            # 'c_{m_{w}}': 1,
            'C_{L_{max}}': 2,

            '\\alpha_{max,h}': 2.5,

     
            'SM_{min}': 0.5,

            '\\Delta x_{CG}': 4,

##            'x_{CG}': [17, 18],
             #think about how to constrain this
             'x_w': 19,

            'mac': 2,
            }
    mission = Mission(aircraft)
    m = Model(mission['W_{f_{total}}'], [aircraft, mission], substitutions)
    sol = m.localsolve(solver='mosek', verbosity = 2)
##    bounds, sol = m.determine_unbounded_variables(m, solver="mosek",verbosity=4, iteration_limit=100)

    # if plot == True:

    #     substitutions = {      
    # ##            'V_{stall}': 120,
    #             'ReqRng': ('sweep', np.linspace(500,2000,4)),
    #             'CruiseAlt': 30000, #('sweep', np.linspace(20000,40000,4)),
    #             'numeng': 1,
    # ##            'W_{Load_max}': 6664,
    #             'W_{pax}': 91 * 9.81,
    #             'n_{pax}': 150,
    #             'pax_{area}': 1,
    # ##            'C_{D_{fuse}}': .005, #assumes flat plate turbulent flow, from wikipedia
    #             'e': .9,
    #             'b_{max}': 35,

    #              'V_{ne}': 144,
    #              'C_{L_{hmax}}': 2.5,

    #              '\\rho_0': 1.225,
    #              '\\tan(\\Lambda_{ht})': tan(30*pi/180),
    #              'w_{fuse}': 6,

    # ##            'l_{fuse}': 30,
    #             'c_{m_{w}}': 1,
    #             'C_{L_{max}}': 2,

    #             '\\alpha_{max,h}': 2.5,

         
    #             'SM_{min}': 0.5,

    #             '\\Delta x_{CG}': 4,

    # ##            'x_{CG}': [17, 18],
    #              #think about how to constrain this
    #              'x_w': 19,

    #             'mac': 2,
    #             }
               
    #     mission = Mission(aircraft)
    #     m = Model(mission['W_{f_{total}}'], [aircraft, mission], substitutions)
    #     solRsweep = m.localsolve(solver='mosek', verbosity = 4)
        
    #     plt.plot(solRsweep('ReqRng'), solRsweep('AR_h'), '-r')
    #     plt.xlabel('Mission Range [nm]')
    #     plt.ylabel('Horizontal Tail Aspect Ratio')
    #     plt.title('Horizontal Tail Aspect Ratio vs Range')
    # ##    plt.savefig('HT_Sweeps/VT_rng_RC.pdf')
    #     plt.show()

    #     plt.plot(solRsweep('ReqRng'), solRsweep('S_h'), '-r')
    #     plt.xlabel('Mission Range [nm]')
    #     plt.ylabel('Horizontal Tail Area [m$^2$]')
    #     plt.title('Horizontal Tail Area vs Range')
    # ##    plt.savefig('HT_Sweeps/VT_rng_RC.pdf')
    #     plt.show()

        
    #     substitutions = {      
    # ##            'V_{stall}': 120,
    #             'ReqRng': 500, #('sweep', np.linspace(500,2000,4)),
    #             'CruiseAlt': ('sweep', np.linspace(20000,40000,4)),
    #             'numeng': 1,
    # ##            'W_{Load_max}': 6664,
    #             'W_{pax}': 91 * 9.81,
    #             'n_{pax}': 150,
    #             'pax_{area}': 1,
    # ##            'C_{D_{fuse}}': .005, #assumes flat plate turbulent flow, from wikipedia
    #             'e': .9,
    #             'b_{max}': 35,

    #              'V_{ne}': 144,
    #              'C_{L_{hmax}}': 2.5,

    #              '\\rho_0': 1.225,
    #              '\\tan(\\Lambda_{ht})': tan(30*pi/180),
    #              'w_{fuse}': 6,

    # ##            'l_{fuse}': 30,
    #             'c_{m_{w}}': 1,
    #             'C_{L_{max}}': 2,

    #             '\\alpha_{max,h}': 2.5,

         
    #             'SM_{min}': 0.5,

    #             '\\Delta x_{CG}': 4,

    # ##            'x_{CG}': [17, 18],
    #              #think about how to constrain this
    #              'x_w': 19,

    #             'mac': 2,
    #             }
               
    #     mission = Mission(aircraft)
    #     m = Model(mission['W_{f_{total}}'], [aircraft, mission], substitutions)
    #     solAltsweep = m.localsolve(solver='mosek', verbosity = 4)
        
    #     plt.plot(solAltsweep('CruiseAlt'), solAltsweep('AR_h'), '-r')
    #     plt.xlabel('Cruise Altitude [ft]')
    #     plt.ylabel('Horizontal Tail Aspect Ratio')
    #     plt.title('Horizontal Tail Aspect Ratio vs Cruise Altitude')
    # ##    plt.savefig('HT_Sweeps/VT_rng_RC.pdf')
    #     plt.show()

    #     plt.plot(solAltsweep('CruiseAlt'), solAltsweep('S_h'), '-r')
    #     plt.xlabel('Cruise Altitude [ft]')
    #     plt.ylabel('Horizontal Tail Area [m$^2$]')
    #     plt.title('Horizontal Tail Area vs Cruise Altitude')
    # ##    plt.savefig('HT_Sweeps/VT_rng_RC.pdf')
    #     plt.show()
