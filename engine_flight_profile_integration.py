"""Simple commercial aircraft flight profile and aircraft model"""
from numpy import pi, cos, tan, sin
import numpy as np
from gpkit import Variable, Model, units, SignomialsEnabled, Vectorize
from gpkit.constraints.sigeq import SignomialEquality as SignomialEquality
from gpkit.tools import te_exp_minus1
from gpkit.constraints.tight import Tight as TCS
import matplotlib.pyplot as plt
from TASOPT_engine import Engine
from gpkit.small_scripts import mag
from Wing_simple_performance import Wing, WingPerformance
from D8_VT_yaw_rate_and_EO_simple_profile import VerticalTail, VerticalTailPerformance
from D8_HT_simple_profile import HorizontalTail, HorizontalTailPerformance
from gpkit.constraints.bounded import Bounded as BCS

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
    def setup(self, Nclimb, Ncruise, enginestate, **kwargs):
        #create submodels
        self.fuse = Fuselage()
        self.wing = Wing()
        self.engine = Engine(0, True, Nclimb + Ncruise, enginestate)
        self.VT = VerticalTail()
        self.HT = HorizontalTail()

        #variable definitions
        numeng = Variable('numeng', '-', 'Number of Engines')
        SMmin = Variable('SM_{min}','-', 'Minimum Static Margin')
        dxCG = Variable('\\Delta x_{CG}', 'm', 'Max CG Travel Range')
        xCGmin = Variable('x_{CG_{min}}','m','Maximum Forward CG')

        self.components = [self.fuse, self.wing, self.engine, self.VT, self.HT]

        constraints = []
        
        with SignomialsEnabled():
            constraints.extend([
                            # Lift curve slope ratio for HT and Wing
                            SignomialEquality(self.HT['m_{ratio}']*(1+2/self.wing['AR']), 1 + 2/self.HT['AR_h']),

                            # HT Location and Volume Coefficient
                            self.HT['x_{CG_{ht}}'] <= self.fuse['l_{fuse}'],
##                            self.fuse['x_{tail}'] == self.VT['x_{CG_{vt}}'],
                            TCS([self.HT['V_{h}'] == self.HT['S_h']*self.HT['l_{ht}']/(self.wing['S']*self.wing['mac'])]),
                            self.HT['V_{h}'] >= 0.1,

                            # HT Max Loading
                            TCS([self.HT['L_{{max}_h}'] >= 0.5*self.VT['\\rho_{TO}']*self.VT['V_{ne}']**2*self.HT['S_h']*self.HT['C_{L_{hmax}}']]),

                            # HT/VT joint constraint
                            self.HT['b_{ht}']/self.fuse['w_{fuse}']*self.HT['\lambda_h']*self.HT['c_{root_h}'] == self.HT['c_{attach}'],

                            # VT height constraint (4*engine radius)
##                            self.VT['b_{vt}']**2 >= 16.*self.engine['A_2']/np.pi,

                            # VT root chord constraint #TODO find better constraint
##                            self.VT['c_{root_{vt}}'] <= self.fuse['l_{cone}'],
                            ])

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
    def  setup(self, aircraft, state, **kwargs):
        #make submodels
        self.aircraft = aircraft
        self.wingP = aircraft.wing.dynamic(state)
        self.fuseP = aircraft.fuse.dynamic(state)
        self.VTP = aircraft.VT.dynamic(aircraft.fuse, state)
        self.HTP = aircraft.HT.dynamic(aircraft.fuse,aircraft.wing,state)

        self.Pmodels = [self.wingP, self.fuseP, self.VTP, self.HTP]

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

        xAC = Variable('x_{AC}','m','Aerodynamic Center of Aircraft')
        xCG = Variable('x_{CG}','m','Center of Gravity of Aircraft')

        constraints = []

        constraints.extend([
            #speed must be greater than stall speed
            state['V'] >= Vstall,


            #Figure out how to delete
            Vstall == 120*units('kts'),
            WLoadmax == 6664 * units('N/m^2'),

            #compute the drag
            TCS([D >= self.wingP['D_{wing}'] + self.fuseP['D_{fuse}'] + self.VTP['D_{vt}'] + self.HTP['D_{ht}']]),

            #constraint CL and compute the wing loading
            W_avg == .5*self.wingP['C_{L}']*self.aircraft['S']*state.atm['\\rho']*state['V']**2,      
            WLoad == .5*self.wingP['C_{L}']*self.aircraft['S']*state.atm['\\rho']*state['V']**2/self.aircraft.wing['S'],

            #set average weight equal to the geometric avg of start and end weight
            W_avg == (W_start * W_end)**.5,

            #constrain the max wing loading
            WLoad <= WLoadmax,

            #time unit conversion
            t == thours,

            #make lift equal weight --> small angle approx in climb
            self.wingP['L_w'] == W_avg,
           # Tail downforce penalty to wing lift
##            self.wingP['L_w'] >= W_avg + self.HTP['L_h'],

           # Center of gravity constraints #TODO Refine
            xCG <= 0.7*aircraft.fuse['l_{fuse}'],
            xCG >= aircraft['x_{CG_{min}}'],
            xAC >= xCG,

            # Tail aspect ratio and lift constraints
            aircraft.HT['AR_h'] >= 6, #TODO change to tip Re constraint
            self.HTP['C_{L_h}'] >= 0.01, #TODO remove

            # HT/VT moment arm constraints
            TCS([aircraft.HT['l_{ht}'] + xCG <= aircraft.HT['x_{CG_{ht}}'] ]),
            TCS([aircraft.VT['l_{vt}'] + xCG <= aircraft.VT['x_{CG_{vt}}']]),

             # Wing location and AC constraints
            TCS([xCG + self.HTP['\\Delta x_{{trail}_h}'] <= aircraft.fuse['l_{fuse}']]), #TODO tighten
            xAC == aircraft['x_{wing}'], #TODO improve, only works because cmw == 0.1

            aircraft.fuse['x_{wing}'] >= aircraft.fuse['l_{fuse}']*0.5, #TODO remove
            aircraft.fuse['x_{wing}'] <= aircraft.fuse['l_{fuse}']*0.6, #TODO remove
            aircraft.fuse['x_{wing}'] == aircraft.wing['x_w'],

            TCS([aircraft.HT['x_{CG_{ht}}'] >= xCG + 0.5*(self.HTP['\\Delta x_{{trail}_h}'] + self.HTP['\\Delta x_{{lead}_h}'])]), #TODO tighten
            ])

        with SignomialsEnabled():
            constraints.extend([
                # Aircraft trim conditions
                TCS([xAC/aircraft.wing['mac'] <= self.wingP['c_{m_{w}}']/self.wingP['C_{L}'] + xCG/aircraft.wing['mac'] + \
                                  aircraft.HT['V_{h}']*(self.HTP['C_{L_h}']/self.wingP['C_{L}'])]),
                
                # Wing location and AC constraints
                SignomialEquality(xAC,xCG + self.HTP['\\Delta x_w']),

                
                # Static margin constraint with and without dxCG #TODO validate if this works as intended
                TCS([aircraft['SM_{min}'] + aircraft['\\Delta x_{CG}']/aircraft.wing['mac'] <=
                                                aircraft.HT['V_{h}']*aircraft.HT['m_{ratio}'] \
                                              + self.wingP['c_{m_{w}}']/aircraft.wing['C_{L_{wmax}}'] + \
                                                aircraft.HT['V_{h}']*aircraft.HT['C_{L_{hmax}}']/aircraft.wing['C_{L_{wmax}}']]), # [SP]
                TCS([aircraft['SM_{min}'] <=
                                                aircraft.HT['V_{h}']*aircraft.HT['m_{ratio}'] \
                                              + self.wingP['c_{m_{w}}']/aircraft.wing['C_{L_{wmax}}'] + \
                                                aircraft.HT['V_{h}']*aircraft.HT['C_{L_{hmax}}']/aircraft.wing['C_{L_{wmax}}']]), # [SP]
                ])

        return constraints, self.Pmodels

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
                                  
        #variable definitions
        theta = Variable('\\theta', '-', 'Aircraft Climb Angle')
        excessP = Variable('excessP', 'W', 'Excess Power During Climb')
        RC = Variable('RC', 'feet/min', 'Rate of Climb/Decent')
        dhft = Variable('dhft', 'feet', 'Change in Altitude Per Climb Segment [feet]')
        RngClimb = Variable('RngClimb', 'nautical_miles', 'Down Range Covered in Each Climb Segment')

        #constraints
        constraints = []
        
        constraints.extend([ 
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
                        
        #variable definitions
        z_bre = Variable('z_{bre}', '-', 'Breguet Parameter')
        Rng = Variable('Rng', 'nautical_miles', 'Cruise Segment Range')

        constraints = []

        constraints.extend([
             #taylor series expansion to get the weight term
             TCS([self.aircraftP['W_{burn}']/self.aircraftP['W_{end}'] >=
                  te_exp_minus1(z_bre, nterm=3)]),

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

class StateLinking(Model):
    """
    link all the state model variables
    """
    def setup(self, climbstate, cruisestate, enginestate, Nclimb, Ncruise):
        statevarkeys = ['p_{sl}', 'T_{sl}', 'L_{atm}', 'M_{atm}', 'P_{atm}', 'R_{atm}',
                        '\\rho', 'T_{atm}', '\\mu', 'T_s', 'C_1', 'h', 'hft', 'V', 'a', 'R', '\\gamma', 'M']
        constraints = []
        for i in range(len(statevarkeys)):
            varkey = statevarkeys[i]
            for i in range(Nclimb):
                constraints.extend([
                    climbstate[varkey][i] == enginestate[varkey][i]
                    ])
            for i in range(Ncruise):
                constraints.extend([
                    cruisestate[varkey][i] == enginestate[varkey][i+Nclimb]
                    ])           
        
        return constraints

class FlightState(Model):
    """
    creates atm model for each flight segment, has variables
    such as veloicty and altitude
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
            #compute the speed of sound with the state
            a  == (gamma * R * self.atm['T_{atm}'])**.5,

            #compute the mach number
            V == M * a,
            ])

        #build the model
        return self.alt, self.atm, constraints

class Altitude(Model):
    """
    holds the altitdue variable
    """
    def setup(self, **kwargs):
        #define altitude variables
        h = Variable('h', 'm', 'Segment Altitude [meters]')
        hft = Variable('hft', 'feet', 'Segment Altitude [feet]')

        constraints = []

        constraints.extend([
            h == hft, #convert the units on altitude
            ])

        return constraints
            

class Atmosphere(Model):
    def setup(self, alt, **kwargs):
##        g = Variable('g', 'm/s^2', 'Gravitational acceleration')
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
                SignomialEquality(T_sl, T_atm + L_atm*alt['h']),

                #constraint on mu
                SignomialEquality((T_atm + T_s) * mu, C_1 * T_atm**1.5),
                ]

        #like to use a local subs here in the future
        subs = None

        return constraints

class Fuselage(Model):
    """
    place holder fuselage model
    """
    def setup(self, **kwargs):
        #new variables
        n_pax = Variable('n_{pax}', '-', 'Number of Passengers to Carry')
                           
        #weight variables
        W_payload = Variable('W_{payload}', 'N', 'Aircraft Payload Weight')
        W_e = Variable('W_{e}', 'N', 'Empty Weight of Aircraft')
        W_pax = Variable('W_{pax}', 'N', 'Estimated Average Passenger Weight, Includes Baggage')
        Wfuse = Variable('w_{fuse}', 'm', 'Fuselage Width')

        A_fuse = Variable('A_{fuse}', 'm^2', 'Estimated Fuselage Area')
        pax_area = Variable('pax_{area}', 'm^2', 'Estimated Fuselage Area per Passenger')

        lfuse   = Variable('l_{fuse}', 'm', 'Fuselage length')

        xwing = Variable('x_{wing}', 'm', 'Wing Location')

        constraints = []
        
        constraints.extend([
            #compute fuselage area for drag approximation
            A_fuse == pax_area * n_pax,

            #constraints on the various weights
            W_payload == n_pax * W_pax,
            
            #estimate based on TASOPT 737 model
            W_e == .75*W_payload,

            lfuse == lfuse,

            Wfuse == 5*units('m'),
            
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
##        xCG    = Variable('x_{CG}', 'm', 'x-location of CG')
        #constraints
        constraints = []

        constraints.extend([
            Dfuse == Cdfuse * (.5 * fuse['A_{fuse}'] * state.atm['\\rho'] * state['V']**2),

            Cdfuse == .005,

##            xCG == 18*units('m'),
            ])

        return constraints

class Mission(Model):
    """
    mission class, links together all subclasses
    """
    def setup(self, substitutions = None, **kwargs):
        #define the number of each flight segment
        Nclimb = 2
        Ncruise = 2

       # vectorize
        with Vectorize(Nclimb + Ncruise):
            enginestate = FlightState()

        ac = Aircraft(Nclimb, Ncruise, enginestate)

        #Vectorize
        with Vectorize(Nclimb):
            climb = ClimbSegment(ac)

        with Vectorize(Ncruise):
            cruise = CruiseSegment(ac)

        statelinking = StateLinking(climb.state, cruise.state, enginestate, Nclimb, Ncruise)

        #declare new variables
        W_ftotal = Variable('W_{f_{total}}', 'N', 'Total Fuel Weight')
        W_fclimb = Variable('W_{f_{climb}}', 'N', 'Fuel Weight Burned in Climb')
        W_fcruise = Variable('W_{f_{cruise}}', 'N', 'Fuel Weight Burned in Cruise')
        W_total = Variable('W_{total}', 'N', 'Total Aircraft Weight')
        W_dry = Variable('W_{dry}', 'N', 'Total Zero Fuel Aircraft Weight')
        CruiseAlt = Variable('CruiseAlt', 'ft', 'Cruise Altitude [feet]')
        ReqRng = Variable('ReqRng', 'nautical_miles', 'Required Cruise Range')

        RCmin = Variable('RC_{min}', 'ft/min', 'Minimum allowed climb rate')

        h = climb['h']
        hftClimb = climb['hft']
        dhft = climb['dhft']
        hftCruise = cruise['hft']

        #make overall constraints
        constraints = []

        constraints.extend([
            #weight constraints
            TCS([W_dry + W_ftotal <= W_total]),
            TCS([ac['W_{e}'] + ac['W_{payload}'] + ac['numeng'] * ac['W_{engine}'] + ac.wing['W_{struct}'] + 2*ac.VT['W_{struct}'] + ac.HT['W_{struct}'] <= W_dry]),

            climb['W_{start}'][0] == W_total,
            climb['W_{end}'][-1] == cruise['W_{start}'][0],

            # similar constraint 1
            TCS([climb['W_{start}'] >= climb['W_{end}'] + climb['W_{burn}']]),
            # similar constraint 2
            TCS([cruise['W_{start}'] >= cruise['W_{end}'] + cruise['W_{burn}']]),

            climb['W_{start}'][1:] == climb['W_{end}'][:-1],
            cruise['W_{start}'][1:] == cruise['W_{end}'][:-1],

            TCS([W_dry <= cruise['W_{end}'][-1]]),

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

            #set the range for each cruise segment, doesn't take credit for climb
            #down range disatnce covered
            cruise.cruiseP['Rng'] == ReqRng/(Ncruise),

            #compute fuel burn from TSFC
            cruise['W_{burn}'][0] == ac['numeng']*ac.engine['TSFC'][2] * cruise['thr'] * ac.engine['F'][2],
            cruise['W_{burn}'][1] == ac['numeng']*ac.engine['TSFC'][3] * cruise['thr'] * ac.engine['F'][3],              
            climb['W_{burn}'][0] == ac['numeng']*ac.engine['TSFC'][0] * climb['thr'] * ac.engine['F'][0],
            climb['W_{burn}'][1] == ac['numeng']*ac.engine['TSFC'][1] * climb['thr'] * ac.engine['F'][1],

            #min climb rate constraint
            climb['RC'][0] >= RCmin,
            ])

        M2 = .8
        M25 = .6
        M4a = .1025
        Mexit = 1
        M0 = .8

        engineclimb = [
            ac.engine.engineP['M_2'][0] == climb['M'][0],
            ac.engine.engineP['M_2'][1] == climb['M'][1],
            ac.engine.engineP['M_{2.5}'][0] == M25,
##            ac.engine.engineP['M_{2.5}'][1] == M25,
            ac.engine.compressor['hold_{2}'] == 1+.5*(1.398-1)*M2**2,
            ac.engine.compressor['hold_{2.5}'] == 1+.5*(1.354-1)*M25**2,
            ac.engine.compressor['c1'] == 1+.5*(.401)*M0**2,

            #constraint on drag and thrust
            ac['numeng']*ac.engine['F_{spec}'][0] >= climb['D'][0] + climb['W_{avg}'][0] * climb['\\theta'][0],
            ac['numeng']*ac.engine['F_{spec}'][1] >= climb['D'][1] + climb['W_{avg}'][1] * climb['\\theta'][1],

            #climb rate constraints
            TCS([climb['excessP'][0] + climb.state['V'][0] * climb['D'][0] <=  climb.state['V'][0] * ac['numeng'] * ac.engine['F_{spec}'][0]]),
            TCS([climb['excessP'][1] + climb.state['V'][1] * climb['D'][1] <=  climb.state['V'][1] * ac['numeng'] * ac.engine['F_{spec}'][1]]),
            ]

        M2 = .8
        M25 = .6
        M4a = .1025
        Mexit = 1
        M0 = .8

        enginecruise = [
            ac.engine.engineP['M_2'][2] == cruise['M'][0],
            ac.engine.engineP['M_2'][3] == cruise['M'][1],
            ac.engine.engineP['M_{2.5}'][2] == M25,
            ac.engine.engineP['M_{2.5}'][3] == M25,

            #steady level flight constraint on D 
            cruise['D'][0] == ac['numeng'] * ac.engine['F_{spec}'][2],
            cruise['D'][1] == ac['numeng'] * ac.engine['F_{spec}'][3],

            #breguet range eqn
            TCS([cruise['z_{bre}'][0] >= (ac.engine['TSFC'][2] * cruise['thr'][0]*
            cruise['D'][0]) / cruise['W_{avg}'][0]]),
            TCS([cruise['z_{bre}'][1] >= (ac.engine['TSFC'][3] * cruise['thr'][1]*
            cruise['D'][1]) / cruise['W_{avg}'][1]]),
            ]

        wing = [
            #wing constraints
            ac.wing['W_{fuel_{wing}}'] == W_ftotal,
            climb.climbP.wingP['L_w'] == climb.climbP.aircraftP['W_{avg}'],
            cruise.cruiseP.wingP['L_w'] == cruise.cruiseP.aircraftP['W_{avg}'],
            climb['c_{m_{w}}'] == .10, # for boundedness
            cruise['c_{m_{w}}'] == .10, # for boundedness
            ac.wing['V_{ne}'] == 144*units('m/s'),
            ]

        #VT constriants
        VT = [
            ac.VT['T_e'] == ac.engine['F'][0],

            # Drag of a windmilling engine
            ac.VT['D_{wm}'] >= 0.5*ac.VT['\\rho_{TO}']*ac.VT['V_1']**2*ac.engine['A_2']*ac.VT['C_{D_{wm}}'],
            
            ac.VT['x_{CG_{vt}}'] <= ac.fuse['l_{fuse}'],
            
            #VTP constraints
            ac.fuse['l_{fuse}'] >= ac.VT['\\Delta x_{lead_v}'] + climb.climbP.aircraftP['x_{CG}'],
            ac.VT['x_{CG_{vt}}'] >= climb.climbP.aircraftP['x_{CG}']+(ac.VT['\\Delta x_{lead_v}']+ac.VT['\\Delta x_{trail_v}'])/2,

            ac.fuse['l_{fuse}'] >= ac.VT['\\Delta x_{lead_v}'] + cruise.cruiseP.aircraftP['x_{CG}'],
            ac.VT['x_{CG_{vt}}'] >= cruise.cruiseP.aircraftP['x_{CG}']+(ac.VT['\\Delta x_{lead_v}']+ac.VT['\\Delta x_{trail_v}'])/2,

            ac.VT['V_{ne}'] == 144*units('m/s'),
            ac.VT['\\rho_{TO}'] == 1.225 * units('kg/m^3'),
            ]

##        HT = [
##            #compute the HT chord at its attachment point to the VT
##            (aircraft.HT['b_{ht}']/aircraft.fuse['w_{fuse}'])*aircraft.HT['\lambda_h']*aircraft.HT['c_{root_h}'] == aircraft.HT['c_{attach}']
##            ]
        
        return constraints + ac + climb + cruise + enginecruise + engineclimb + enginestate + statelinking + VT + wing

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
    def determine_unbounded_variables(self, model, verbosity=4,
                                      eps=1e-30, lower=None, upper=None, **kwargs):
        "Returns labeled dictionary of unbounded variables."
        m = self.bound_all_variables(model, eps, lower, upper)
        sol = m.localsolve(solver='mosek', verbosity=4, iteration_limit = 100, **kwargs)
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
    M4a = .1025
    fan = 1.685
    lpc  = 4.744
    hpc = 3.75
 
    sweep = 30 #[deg]
    
    substitutions = {      
            'ReqRng': 500, #('sweep', np.linspace(500,2000,4)),
            'numeng': 2,
            'W_{pax}': 91 * 9.81,
            'n_{pax}': 150,
            'pax_{area}': 1,
##            'e': .9,

            #wing subs
            'C_{L_{wmax}}': 2.5,
            '\\alpha_{max,w}': 0.1, # (6 deg)
            '\\cos(\\Lambda)': cos(sweep*pi/180),
            '\\eta': 0.97,
            '\\rho_0': 1.225,
            '\\rho_{fuel}': 817, # Kerosene [TASOPT]
            '\\tan(\\Lambda)': tan(sweep*pi/180),

            #VT subs
           'C_{D_{wm}}': 0.5, # [2]
           'C_{L_{vmax}}': 2.6, # [2]
           'V_1': 70,
##           '\\rho_{TO}': 1.225,
           '\\tan(\\Lambda_{vt})': np.tan(40*np.pi/180),
           'c_{l_{vtEO}}': 0.5,
           'e_v': 0.8,
           'l_{fuse}': 39,
           'y_{eng}': 4.83, # [3]
           'V_{land}': 72,
           'I_{z}': 12495000, #estimate for late model 737 at max takeoff weight (m l^2/12)
           '\\dot{r}_{req}': 0.174533, #10 deg/s yaw rate
            'N_{spar}': 2,

            # HT substitutions
            '\\alpha_{max,h}': 2.5,
            '\\tan(\\Lambda_{ht})': tan(30*pi/180),
            'C_{L_{hmax}}': 2.5,
            'SM_{min}': 0.05,
            '\\Delta x_{CG}': 2.0*units('m'),
            'x_{CG_{min}}' : 1,#13.0*units('m'),

            #engine subs
            '\\pi_{tn}': .98,
            '\pi_{b}': .94,
            '\pi_{d}': .98,
            '\pi_{fn}': .98,
            'T_{ref}': 288.15,
            'P_{ref}': 101.325,
            '\eta_{HPshaft}': .97,
            '\eta_{LPshaft}': .97,
            'eta_{B}': .9827,

            '\pi_{f_D}': fan,
            '\pi_{hc_D}': hpc,
            '\pi_{lc_D}': lpc,

            '\\alpha_{OD}': 5.105,

            'hold_{4a}': 1+.5*(1.313-1)*M4a**2,#sol('hold_{4a}'),
            'r_{uc}': .01,
            '\\alpha_c': .19036,
            'T_{t_f}': 435,

            'M_{takeoff}': .9556,

            'G_f': 1,

            'h_f': 43.03,

            'Cp_t1': 1280,
            'Cp_t2': 1184,
            'Cp_c': 1216,

            'RC_{min}': 1000,
            }
           
    mission = Mission()
    m = Model(mission['W_{f_{total}}'], mission, substitutions)
##    m = Model(mission['W_{f_{total}}'], BCS(mission), substitutions)
    sol = m.localsolve(solver='mosek', verbosity = 4)
##    bounds, sol = mission.determine_unbounded_variables(m)
