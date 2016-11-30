"""Creates a D8 empenage static and performance model"""
from numpy import pi, tan
import numpy as np
from gpkit import Variable, Model, units, SignomialsEnabled, Vectorize, SignomialEquality
from gpkit.tools import te_exp_minus1
from gpkit.constraints.tight import TightConstraintSet as TCS
import matplotlib.pyplot as plt

from D8_HT_simple_profile import HorizontalTail
from D8_VT_yaw_rate_and_EO_simple_profile import VerticalTail

#only needed for the local bounded debugging tool
from collections import defaultdict
from gpkit.small_scripts import mag

class Aircraft(Model):
    "Aircraft class"
    def __init__(self, **kwargs):
        #create submodels
        self.fuse = Fuselage()
        self.wing = Wing()
        self.engine = Engine()
        self.ht = HorizontalTail()
        self.VT = VerticalTail()

        #variable definitions
        numeng = Variable('numeng', '-', 'Number of Engines')

        constraints = []

        constraints.extend([
            numeng == numeng, #need numeng in the model
            ])

        self.components = [self.fuse, self.wing, self.engine, self.ht, self.VT]
##        self.components = [self.fuse, self.wing, self.engine, self.VT]

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
        self.htP = aircraft.ht.dynamic(aircraft.fuse, aircraft.wing, state)
        self.VTP = aircraft.VT.dynamic(aircraft.fuse, state)

        self.Pmodels = [self.wingP, self.fuseP, self.engineP, self.htP, self.VTP]
##        self.Pmodels = [self.wingP, self.fuseP, self.engineP, self.VTP]

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
                TCS([D >= self.wingP['D_{wing}'] + self.fuseP['D_{fuse}'] + self.htP['D_{ht}'] + self.VTP['D_{vt}']]),
##                TCS([D >= self.wingP['D_{wing}'] + self.fuseP['D_{fuse}'] +  self.VTP['D_{vt}']]),

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
                SignomialEquality(T_sl, T_atm + L_atm*alt['h']),

                #constraint on mu
                SignomialEquality((T_atm + T_s) * mu, C_1 * T_atm**1.5),
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
        W_engine = Variable('W_{engine}', 1000, 'N', 'Weight of a Single Turbofan Engine')
        A2 = Variable('A_2', 'm^2', 'Fan Area')
        
        constraints = []

        constraints.extend([
            W_engine == W_engine,

            A2 == A2,
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

        dum1 = Variable('dum1', 124.58, 'm^2')
        dum2 = Variable('dum2', 105384.1524, 'N')

        cwma    = Variable('\\bar{c}_w', 'm',
                          'Mean aerodynamic chord (wing)')
        Vne     = Variable('V_{ne}', 'm/s', 'Never exceed velocity')
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

            span <= span_max,

            #compute K for the aircraft
            K == (pi * e * AR)**-1,

            Vne == 144*units('m/s'),
            cwma == cwma,
            cmw == cmw,
            CLmax == CLmax,

            xw == xw,
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

        alpha   = Variable('\\alpha', '-', 'Horizontal tail angle of attack')

        #constraints
        constraints = []

        constraints.extend([
            #airfoil drag constraint
            Lwing == (.5*wing['S']*state.atm['\\rho']*state['V']**2)*CL,
            CL == 2*3.14*alpha,
            
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

        xcg     = Variable('x_{CG}', 'm', 'CG location')
        
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
    def __init__(self, ac, substitutions = None, **kwargs):
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
        Wnofuel = Variable('W_{nofuel}', 'N', 'Aircraft Weight witn Zero Fuel')
        CruiseAlt = Variable('CruiseAlt', 'ft', 'Cruise Altitude [feet]')
        ReqRng = Variable('ReqRng', 'nautical_miles', 'Required Cruise Range')

        xcg     = Variable('x_{CG}', 'm', 'CG location')

        h = climb['h']
        hftClimb = climb['hft']
        dhft = climb['dhft']
        hftCruise = cruise['hft']

        #make overall constraints
        constraints = []

        constraints.extend([
            #weight constraints
            TCS([ac['W_{e}'] + ac['W_{payload}'] + ac['numeng'] * ac['W_{engine}'] + ac['W_{wing}'] + ac.ht['W_{struct}'] + 2*ac.VT['W_{struct}'] <= Wnofuel]),
##            TCS([ac['W_{e}'] + ac['W_{payload}'] + ac['numeng'] * ac['W_{engine}'] + ac['W_{wing}'] + 2*ac.VT['W_{struct}'] <= Wnofuel]),
            TCS([Wnofuel + W_ftotal <= W_total]),

            climb['W_{start}'][0] == W_total,
            climb['W_{end}'][-1] == cruise['W_{start}'][0],

            # similar constraint 1
            TCS([climb['W_{start}'] >= climb['W_{end}'] + climb['W_{burn}']]),
            # similar constraint 2
            TCS([cruise['W_{start}'] >= cruise['W_{end}'] + cruise['W_{burn}']]),

            climb['W_{start}'][1:] == climb['W_{end}'][:-1],
            cruise['W_{start}'][1:] == cruise['W_{end}'][:-1],

            TCS([Wnofuel <= cruise['W_{end}'][-1]]),

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


            xcg == 15*units('m'),
            ])
        
        #Horizontal Tail Constraints
        with SignomialsEnabled():
            constraints.extend([
                cruise['C_{L_h}'] <= 1.1*cruise['C_{L_{ah}}']*cruise['\\alpha'],
                cruise['C_{L_h}'] >= 0.9*cruise['C_{L_{ah}}']*cruise['\\alpha'],

                # Trim condidtion for each flight segment
                TCS([cruise['x_{ac}']/ac.wing['\\bar{c}_w'] <= ac.wing['c_{m_{w}}']/cruise['C_{L}'] + xcg/ac.wing['\\bar{c}_w'] + ac.ht['V_{h}']*(cruise['C_{L_h}']/cruise['C_{L}'])]),
                
                ac.ht['L_{{max}_h}'] == 0.5*ac.ht['\\rho_0']*ac.wing['V_{ne}']**2*ac.ht['S_h']*ac.ht['C_{L_{hmax}}'],
                #compute mrat, is a signomial equality
                SignomialEquality(ac.ht['m_{ratio}']*(1+2/ac.wing['AR']), 1 + 2/ac.ht['AR_h']),

                #tail volume coefficient
                
                ac.ht['V_{h}'] == ac.ht['S_h']*ac.ht['x_{CG_{ht}}']/(ac.wing['S']*ac.wing['\\bar{c}_w']),

                #enforce max tail location is the end of the fuselage
                ac.ht['x_{CG_{ht}}'] <= ac.fuse['l_{fuse}'],

                #Stability constraint, is a signomial
                TCS([ac.ht['SM_{min}'] + ac.ht['\\Delta x_{CG}']/ac.wing['\\bar{c}_w'] <= ac.ht['V_{h}']*ac.ht['m_{ratio}'] + ac.wing['c_{m_{w}}']/ac.wing['C_{L_{max}}'] + ac.ht['V_{h}']*ac.ht['C_{L_{hmax}}']/ac.wing['C_{L_{max}}']]),

                ac.ht['l_{ht}'] >= ac.ht['x_{CG_{ht}}'] - xcg,
                TCS([ac.wing['x_w'] >= xcg + cruise['\\Delta x_w']]),
                TCS([xcg + cruise['\\Delta x_{{trail}_h}'] <= ac.fuse['l_{fuse}']], reltol=0.002),

                TCS([ac.ht['x_{CG_{ht}}'] >= xcg+(cruise['\\Delta x_{{lead}_h}']+cruise['\\Delta x_{{trail}_h}'])/2]),

                climb['C_{L_h}'] <= 1.1*climb['C_{L_{ah}}']*climb['\\alpha'],
                climb['C_{L_h}'] >= 0.9*climb['C_{L_{ah}}']*climb['\\alpha'],

                # Trim condidtion for each flight segment
                TCS([climb['x_{ac}']/ac.wing['\\bar{c}_w'] <= ac.wing['c_{m_{w}}']/climb['C_{L}'] + xcg/ac.wing['\\bar{c}_w'] + ac.ht['V_{h}']*(climb['C_{L_h}']/climb['C_{L}'])]),
                



                ac.ht['l_{ht}'] >= ac.ht['x_{CG_{ht}}'] - xcg,
                TCS([ac.wing['x_w'] >= xcg + climb['\\Delta x_w']]),
                TCS([xcg + climb['\\Delta x_{{trail}_h}'] <= ac.fuse['l_{fuse}']], reltol=0.002),
                #compute the aerodynamic center location
                TCS([climb['x_{ac}'] <= xcg + climb['\\Delta x_w'] ]),

##                SignomialEquality(climb['x_{ac}'],xcg + climb['\\Delta x_w'] ),
                TCS([ac.ht['x_{CG_{ht}}'] >= xcg+(climb['\\Delta x_{{lead}_h}']+climb['\\Delta x_{{trail}_h}'])/2]),
                    
                TCS([cruise['x_{ac}'] <= xcg + cruise['\\Delta x_w'] ]),
##                SignomialEquality(cruise['x_{ac}'],xcg + cruise['\\Delta x_w'] ),
                TCS([ac.ht['x_{CG_{ht}}'] >= xcg+(cruise['\\Delta x_{{lead}_h}']+cruise['\\Delta x_{{trail}_h}'])/2]),
                #---------------------------------------------------------#
                #REMOVE FROM OVERALL MODEL
                cruise['\\alpha'] <= ac.ht['\\alpha_{max,h}'],
                climb['\\alpha'] <= ac.ht['\\alpha_{max,h}'],
                    
                

                #compute the HT chord at its attachment point to the VT
                (ac.ht['b_{ht}']/ac.fuse['w_{fuse}'])*ac.ht['\lambda_h']*ac.ht['c_{root_h}'] == ac.ht['c_{attach}'],
                ])
 
        #VT constriants
        constraints.extend([
            ac.VT['T_e'] == climb.climbP.engineP['thrust'][0],

            # Drag of a windmilling engine
            ac.VT['D_{wm}'] >= 0.5*ac.VT['\\rho_{TO}']*ac.VT['V_1']**2*ac.engine['A_2']*ac.VT['C_{D_{wm}}'],
            
            ac.VT['x_{CG_{vt}}'] <= ac.fuse['l_{fuse}'],
            
            #VTP constraints
            ac.fuse['l_{fuse}'] >= ac.VT['\\Delta x_{lead_v}'] + xcg,
            ac.VT['x_{CG_{vt}}'] >= xcg+(ac.VT['\\Delta x_{lead_v}']+ac.VT['\\Delta x_{trail_v}'])/2,

            ac.fuse['l_{fuse}'] >= ac.VT['\\Delta x_{lead_v}'] + xcg,
            ac.VT['x_{CG_{vt}}'] >= xcg+(ac.VT['\\Delta x_{lead_v}']+ac.VT['\\Delta x_{trail_v}'])/2,
            ])

        Model.__init__(self, W_ftotal, constraints + ac + climb + cruise, substitutions)

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


if __name__ == '__main__':
    plot = False
    
    #build required submodels
    ac = Aircraft()
        
    substitutions = {      
            #overall aircraft 
            'ReqRng': 500, #('sweep', np.linspace(500,2000,4)),
            'CruiseAlt': 30000, #('sweep', np.linspace(20000,40000,4)),
            'numeng': 2,
##            'V_{ne}': 144,
            
            #fuselage place holder model subs
            'W_{pax}': 91 * 9.81,
            'n_{pax}': 150,
            'pax_{area}': 1,
##            'C_{D_{fuse}}': .005, #assumes flat plate turbulent flow, from wikipedia
             'w_{fuse}': 6,

            #subs for the HT
            'C_{L_{max}}': 2,
            '\\alpha_{max,h}': 2.5,
            'SM_{min}': 0.5,
            '\\Delta x_{CG}': 4,
            'C_{L_{hmax}}': 2.5,
             '\\rho_0': 1.225,
             '\\tan(\\Lambda_{ht})': tan(30*pi/180),


             #think about how to constrain this
             'x_w': 19,

            #subs into the wing palceholder model
            '\\bar{c}_w': 2,
            'c_{m_{w}}': 1,
            'e': .9,
            'b_{max}': 35,

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
##           'l_{fuse}': 39,
##           'x_{CG}': 18,
           'y_{eng}': 4.83, # [3]
           'V_{land}': 72,
           'I_{z}': 12495000, #estimate for late model 737 at max takeoff weight (m l^2/12)
           '\\dot{r}_{req}': 0.174533, #10 deg/s yaw rate
            'N_{spar}': 2,
            }
           
    m = Mission(ac, substitutions)
    sol = m.localsolve(solver='mosek', verbosity = 4)
##    bounds, sol = m.determine_unbounded_variables(m, solver="mosek",verbosity=4, iteration_limit=100)
    
