"""Simple commercial aircraft flight profile and aircraft model"""
from numpy import pi, cos, tan
import numpy as np
from gpkit import Variable, Model, units, SignomialsEnabled, vectorize
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
        self.wing = Wing()

        #variable definitions
        numeng = Variable('numeng', '-', 'Number of Engines')

        constraints = []

        constraints.extend([
            numeng == numeng, #need numeng in the model
            ])

        self.components = [self.fuse, self.wing, self.engine, self.wing]

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
        self.wingP = aircraft.wing.dynamic(state)

        self.Pmodels = [self.wingP, self.fuseP, self.engineP, self.wingP]

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
            Dfuse == Cdfuse * (.5 * fuse['A_{fuse}'] * state.atm['\\rho'] * state['V']**2),

            Cdfuse == .005,
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
        W_fclimb = Variable('W_{f_{climb}}', 'N', 'Fuel Weight Burned in Climb')
        W_fcruise = Variable('W_{f_{cruise}}', 'N', 'Fuel Weight Burned in Cruise')
        W_total = Variable('W_{total}', 'N', 'Total Aircraft Weight')
        CruiseAlt = Variable('CruiseAlt', 'ft', 'Cruise Altitude [feet]')
        ReqRng = Variable('ReqRng', 'nautical_miles', 'Required Cruise Range')

        h = cls.state['h']
        hftClimb = cls.state['hft']
        dhft = cls.climbP['dhft']
        hftCruise = crs.state['hft']

        #make overall constraints
        constraints = []

        constraints.extend([
            #weight constraints
            TCS([ac['W_{e}'] + ac['W_{payload}'] + W_ftotal + ac['numeng'] * ac['W_{engine}'] + ac.wing.wb['W_{struct}'] <= W_total]),

            cls.climbP.aircraftP['W_{start}'][0] == W_total,
            cls.climbP.aircraftP['W_{end}'][-1] == crs.cruiseP.aircraftP['W_{start}'][0],

            # similar constraint 1
            TCS([cls.climbP.aircraftP['W_{start}'] >= cls.climbP.aircraftP['W_{end}'] + cls.climbP.aircraftP['W_{burn}']]),
            # similar constraint 2
            TCS([crs.cruiseP.aircraftP['W_{start}'] >= crs.cruiseP.aircraftP['W_{end}'] + crs.cruiseP.aircraftP['W_{burn}']]),

            cls.climbP.aircraftP['W_{start}'][1:] == cls.climbP.aircraftP['W_{end}'][:-1],
            crs.cruiseP.aircraftP['W_{start}'][1:] == crs.cruiseP.aircraftP['W_{end}'][:-1],

            TCS([ac['W_{e}'] + ac['W_{payload}'] + ac['numeng'] * ac['W_{engine}'] + ac.wing.wb['W_{struct}'] <= crs.cruiseP.aircraftP['W_{end}'][-1]]),

            TCS([W_ftotal >=  W_fclimb + W_fcruise]),
            TCS([W_fclimb >= sum(cls.climbP['W_{burn}'])]),
            TCS([W_fcruise >= sum(crs.cruiseP['W_{burn}'])]),

            #altitude constraints
            hftCruise == CruiseAlt,
            TCS([hftClimb[1:Ncruise] >= hftClimb[:Ncruise-1] + dhft]),
            TCS([hftClimb[0] >= dhft[0]]),
            hftClimb[-1] <= hftCruise,

            #compute the dh
            dhft == hftCruise/Nclimb,

            #constrain the thrust
            cls.climbP.engineP['thrust'] <= 2 * max(crs.cruiseP.engineP['thrust']),

            #set the range for each cruise segment, doesn't take credit for climb
            #down range disatnce covered
            crs.cruiseP['Rng'] == ReqRng/(Ncruise),

            #set the TSFC
            cls.climbP.engineP['TSFC'] == .7*units('1/hr'),
            crs.cruiseP.engineP['TSFC'] == .5*units('1/hr'),

            #wing constraints
            ac.wing['W_{fuel_{wing}}'] == W_ftotal,
            cls.climbP.wingP['L_w'] == cls.climbP.aircraftP['W_{avg}'],
            crs.cruiseP.wingP['L_w'] == crs.cruiseP.aircraftP['W_{avg}'],
            ])
        
        # Model.__init__(self, W_ftotal + s*units('N'), constraints + ac + cls + crs, subs)
        Model.__init__(self, W_ftotal, constraints + ac + cls + crs, subs)

class Wing(Model):
    """
    Philippe's thesis wing model
    """
    def __init__(self, **kwargs):
        self.wns = WingNoStruct()
        self.wb = WingBox(self.wns)

        Model.__init__(self, None, self.wns + self.wb)
        
    def dynamic(self, state):
        """
        returns an instance of the wing perofrmance model
        """
        return WingPerformance(self, state)

class WingNoStruct(Model):
    """
    Philippe's wing model minus structure
    """
    def __init__(self, **kwargs):
        #declare variables
               #Variables
        Afuel   = Variable('\\bar{A}_{fuel, max}', '-', 'Non-dim. fuel area')
        
        CLwmax  = Variable('C_{L_{wmax}}', '-', 'Max lift coefficient, wing')
        
        
        Vfuel   = Variable('V_{fuel, max}', 'm^3', 'Available fuel volume')
        
        
        
        
        amax    = Variable('\\alpha_{max,w}', '-', 'Max angle of attack')
        
        cosL    = Variable('\\cos(\\Lambda)', '-',
                           'Cosine of quarter-chord sweep angle')
        croot   = Variable('c_{root}', 'm', 'Wing root chord')
        ctip    = Variable('c_{tip}', 'm', 'Wing tip chord')
        
        
        eta     = Variable('\\eta', '-', 'Lift efficiency (diff b/w sectional, actual lift)')
        fl      = Variable('f(\\lambda_w)', '-', 'Empirical efficiency function of taper')
        g       = Variable('g', 9.81, 'm/s^2', 'Gravitational acceleration')
        p       = Variable('p', '-', 'Substituted variable = 1 + 2*taper')
        q       = Variable('q', '-', 'Substituted variable = 1 + taper')
        rho0    = Variable('\\rho_0', 1.225, 'kg/m^3', 'Air density (0 ft)')
        rhofuel = Variable('\\rho_{fuel}', 'kg/m^3', 'Density of fuel')
        tanL    = Variable('\\tan(\\Lambda)', '-',
                           'Tangent of quarter-chord sweep angle')
        taper   = Variable('\\lambda', '-', 'Wing taper ratio')
        tau     = Variable('\\tau', '-', 'Wing thickness/chord ratio')
        tcap    = Variable('t_{cap}' ,'-', 'Non-dim. spar cap thickness')
        tweb    = Variable('t_{web}', '-', 'Non-dim. shear web thickness')
        w       = Variable('w', 0.5, '-', 'Wingbox-width-to-chord ratio')
        xw     = Variable('x_w', 'm', 'Position of wing aerodynamic center')
        ymac    = Variable('y_{\\bar{c}_w}', 'm',
                           'Spanwise location of mean aerodynamic chord')

        #Linked Variables
        AR      = Variable('AR', '-', 'Wing aspect ratio')
        Lmax    = Variable('L_{max}', 'N', 'Maximum load')
        Sw      = Variable('S', 'm^2', 'Wing area')
        Vne     = Variable('V_{ne}', 144, 'm/s', 'Never exceed velocity')
        WfuelWing   = Variable('W_{fuel_{wing}}', 'N', 'Fuel weight')
        b       = Variable('b', 'm', 'Wing span')
        #the following two variables have the same name in the flight profile and
        #will be automatically linked by the linked constraint set
        cwma    = Variable('\\bar{c}_w', 'm',
                          'Mean aerodynamic chord (wing)')
        e       = Variable('e', '-', 'Oswald efficiency factor')

        
        #make cosntraints
        constraints = []

        with SignomialsEnabled():

            constraints.extend([
                 p >= 1 + 2*taper,
                 2*q >= 1 + p,
                 ymac == (b/3)*q/p,
                 TCS([(2./3)*(1+taper+taper**2)*croot/q <= cwma],
                                   reltol=1E-2),

                 taper == ctip/croot,

                 SignomialEquality(Sw, b*(croot + ctip)/2),

                 # Oswald efficiency
                 # Nita, Scholz, "Estimating the Oswald factor from
                 # basic aircraft geometrical parameters"
                 TCS([fl >= 0.0524*taper**4 - 0.15*taper**3
                         + 0.1659*taper**2 - 0.0706*taper + 0.0119],
                    reltol=1E-2),
                TCS([e*(1 + fl*AR) <= 1]),
                taper >= 0.2, # TODO

                # Fuel volume [TASOPT doc]
                TCS([Afuel <= w*0.92*tau]),
                # GP approx of the signomial constraint:
                # Afuel <= (w - 2*tweb)*(0.92*tau - 2*tcap),
                Vfuel <= croot**2 * (b/6) * (1+taper+taper**2)*cosL,
                WfuelWing <= rhofuel*Afuel*Vfuel*g,
                  
                Lmax == 0.5*rho0*Vne**2*Sw*CLwmax,

                tanL == tanL,
                eta == eta,
                amax == amax,
                 xw == xw,
                ])

        Model.__init__(self, None, constraints)

class WingPerformance(Model):
    """
    Wing performance model
    """
    def __init__(self, wing, state, **kwargs):
        self.wing = wing
        
        #declare variables
        #Vector Variables
        alpha   = Variable('\\alpha_w', '-', 'Wing angle of attack')
        CLaw    = Variable('C_{L_{aw}}', '-', 'Lift curve slope, wing')
        
        M       = Variable('M', '-', 'Cruise Mach number')
        Re      = Variable('Re_w', '-', 'Reynolds number (wing)')
        CDp     = Variable('C_{D_{p_w}}', '-',
                           'Wing parasitic drag coefficient')
        CDw     = Variable('C_{d_w}', '-', 'Drag coefficient, wing')
        CLw     = Variable('C_{L}', '-', 'Lift coefficient, wing')
 
        D       = Variable('D_{wing}', 'N', 'Wing drag')
        Lw      = Variable('L_w', 'N', 'Wing lift')
        W0      = Variable('W_0', 'N', 'Weight excluding wing')
        #make constraints
        constraints = []

        with SignomialsEnabled():
            
            constraints.extend([
                Lw == 0.5*state['\\rho']*state['V']**2*self.wing['S']*CLw,

                # DATCOM formula (Mach number makes it SP)

                #THIS IS THE NOT TIGHT CONSTRAINT
                
##                TCS([(self.wing['AR']/self.wing['\\eta'])**2*(1 + self.wing['\\tan(\\Lambda)']**2 - state['M']**2) + 8*pi*self.wing['AR']/CLaw
##                      <= (2*pi*self.wing['AR']/CLaw)**2]),

                SignomialEquality((self.wing['AR']/self.wing['\\eta'])**2*(1 + self.wing['\\tan(\\Lambda)']**2 - state['M']**2) + 8*pi*self.wing['AR']/CLaw
                      , (2*pi*self.wing['AR']/CLaw)**2),
                
                CLw == CLaw*alpha,
                alpha <= self.wing['\\alpha_{max,w}'],

                # Drag
                D == 0.5*state['\\rho']*state['V']**2*self.wing['S']*CDw,
                TCS([CDw >= CDp + CLw**2/(pi*self.wing['e']*self.wing['AR'])]),
                Re == state['\\rho']*state['V']*self.wing['\\bar{c}_w']/state['\\mu'],
                1 >= (2.56*CLw**5.88/(Re**1.54*self.wing['\\tau']**3.32*CDp**2.62)
                   + 3.8e-9*self.wing['\\tau']**6.23/(CLw**0.92*Re**1.38*CDp**9.57)
                   + 2.2e-3*Re**0.14*self.wing['\\tau']**0.033/(CLw**0.01*CDp**0.73)
                   + 6.14e-6*CLw**6.53/(Re**0.99*self.wing['\\tau']**0.52*CDp**5.19)
                   + 1.19e4*CLw**9.78*self.wing['\\tau']**1.76/(Re*CDp**0.91)),
                ])

        Model.__init__(self, None, constraints)

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
        taper = Variable('taper', '-', 'Taper ratio')
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

        if isinstance(surface, WingNoStruct):
            AR = surface['AR']
            b = surface['b']
            S = surface['S']
            p = surface['p']
            q = surface['q']
            tau = surface['\\tau']
            Lmax = surface['L_{max}']

 

        constraints = [
                       # Aspect ratio definition
                       AR == b**2/S,

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
        
        Model.__init__(self, None, constraints, **kwargs)

if __name__ == '__main__':

    sweep = 30 #[deg]
    
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

            #wing subs
            'C_{L_{wmax}}': 2.5,
            'V_{ne}': 144,
            '\\alpha_{max,w}': 0.1, # (6 deg)
            '\\cos(\\Lambda)': cos(sweep*pi/180),
            '\\eta': 0.97,
            '\\rho_0': 1.225,
            '\\rho_{fuel}': 817, # Kerosene [TASOPT]
            '\\tan(\\Lambda)': tan(sweep*pi/180),
            }
           
    m = Mission(substitutions)
##    sol = m.localsolve(solver='mosek', verbosity = 4)

    substitutions = {      
##            'V_{stall}': 120,
            'ReqRng': ('sweep', np.linspace(500,3000,30)),
            'CruiseAlt': 30000,
            'numeng': 1,
##            'W_{Load_max}': 6664,
            'W_{pax}': 91 * 9.81,
            'n_{pax}': 150,
            'pax_{area}': 1,
##            'C_{D_{fuse}}': .005, #assumes flat plate turbulent flow, from wikipedia

            #wing subs
            'C_{L_{wmax}}': 2.5,
            'V_{ne}': 144,
            '\\alpha_{max,w}': 0.1, # (6 deg)
            '\\cos(\\Lambda)': cos(sweep*pi/180),
            '\\eta': 0.97,
            '\\rho_0': 1.225,
            '\\rho_{fuel}': 817, # Kerosene [TASOPT]
            '\\tan(\\Lambda)': tan(sweep*pi/180),
            }
           
    m = Mission(substitutions)
    solRsweep = m.localsolve(solver='mosek', verbosity = 4)

    plt.plot(solRsweep('ReqRng'), solRsweep('W_{struct}'), '-r')
    plt.xlabel('Mission Range [nm]')
    plt.ylabel('Wing Weight [N]')
    plt.title('Wing Weight vs Range')
    plt.savefig('Wing_Sweeps/wing_rng_Wstruct.pdf')
    plt.show()

    plt.plot(solRsweep('ReqRng'), solRsweep('AR'), '-r')
    plt.xlabel('Mission Range [nm]')
    plt.ylabel('Wing Aspect Ratio')
    plt.title('Wing Aspect Ratio vs Range')
    plt.savefig('Wing_Sweeps/wing_rng_AR.pdf')
    plt.show()

    plt.plot(solRsweep('ReqRng'), solRsweep('S'), '-r')
    plt.xlabel('Mission Range [nm]')
    plt.ylabel('VWing Area [m$^2$]')
    plt.title('Wing Area vs Range')
    plt.savefig('Wing_Sweeps/wing_rng_S.pdf')
    plt.show()

    plt.plot(solRsweep('ReqRng'), solRsweep('b'), '-r')
    plt.xlabel('Mission Range [nm]')
    plt.ylabel('VWing Span [m]')
    plt.title('Wing Span vs Range')
    plt.savefig('Wing_Sweeps/wing_rng_b.pdf')
    plt.show()

    plt.plot(solRsweep('ReqRng'), solRsweep['sensitivities']['constants']['V_{ne}'], '-r')
    plt.xlabel('Mission Range [nm]')
    plt.ylabel('Sensitivity of $V_{ne}$')
    plt.title('Sensitivity of $V_{ne}$ vs Range')
    plt.savefig('Wing_Sweeps/wing_rng_SensVne.pdf')
    plt.show()

    plt.plot(solRsweep('ReqRng'), solRsweep['sensitivities']['constants']['N_{lift}'], '-r')
    plt.xlabel('Mission Range [nm]')
    plt.ylabel('Sensitivity to Seciontal Lift Multiplier')
    plt.title('Sensitivity to Seciontal Lift Multiplier vs Range')
    plt.savefig('Wing_Sweeps/wing_rng_Nlift.pdf')
    plt.show()

    plt.plot(solRsweep('ReqRng'), solRsweep['sensitivities']['constants']['C_{L_{wmax}}'], '-r')
    plt.xlabel('Mission Range [nm]')
    plt.ylabel('Sensitivity to Wing $C_{L_{max}}$')
    plt.title('Sensitivity to Wing $C_{L_{max}}$ vs Range')
    plt.savefig('Wing_Sweeps/wing_rng_SensClMax.pdf')
    plt.show()

    plt.plot(solRsweep('ReqRng'), solRsweep['sensitivities']['constants']['ReqRng'], '-r')
    plt.xlabel('Mission Range [nm]')
    plt.ylabel('Sensitivity to Required Range')
    plt.title('Sensitivity to Required Range vs Range')
    plt.savefig('Wing_Sweeps/wing_rng_SensRng.pdf')
    plt.show()

    plt.plot(solRsweep('ReqRng'), solRsweep['sensitivities']['constants']['CruiseAlt'], '-r')
    plt.xlabel('Mission Range [nm]')
    plt.ylabel('Sensitivity to Cruise Altitude')
    plt.title('Sensitivity to Cruise Altitude vs Range')
    plt.savefig('Wing_Sweeps/wing_rng_SensAlt.pdf')
    plt.show()

    substitutions = {      
##            'V_{stall}': 120,
            'ReqRng': 500,
            'CruiseAlt': ('sweep', np.linspace(20000,40000,8)),
            'numeng': 1,
##            'W_{Load_max}': 6664,
            'W_{pax}': 91 * 9.81,
            'n_{pax}': 150,
            'pax_{area}': 1,
##            'C_{D_{fuse}}': .005, #assumes flat plate turbulent flow, from wikipedia

            #wing subs
            'C_{L_{wmax}}': 2.5,
            'V_{ne}': 144,
            '\\alpha_{max,w}': 0.1, # (6 deg)
            '\\cos(\\Lambda)': cos(sweep*pi/180),
            '\\eta': 0.97,
            '\\rho_0': 1.225,
            '\\rho_{fuel}': 817, # Kerosene [TASOPT]
            '\\tan(\\Lambda)': tan(sweep*pi/180),
            }
           
    m = Mission(substitutions)
    solAltsweep = m.localsolve(solver='mosek', verbosity = 4)

    plt.plot(solAltsweep('CruiseAlt'), solAltsweep('W_{struct}'), '-r')
    plt.xlabel('Cruise Altitude [ft]')
    plt.ylabel('Wing Weight [N]')
    plt.title('Wing Weight vs Cruise Altitude')
    plt.savefig('Wing_Sweeps/wing_alt_Wstruct.pdf')
    plt.show()

    plt.plot(solAltsweep('CruiseAlt'), solAltsweep('AR'), '-r')
    plt.xlabel('Cruise Altitude [ft]')
    plt.ylabel('Wing Aspect Ratio')
    plt.title('Wing Aspect Ratio vs Cruise Altitude')
    plt.savefig('Wing_Sweeps/wing_alt_AR.pdf')
    plt.show()

    plt.plot(solAltsweep('CruiseAlt'), solAltsweep('S'), '-r')
    plt.xlabel('Cruise Altitude [ft]')
    plt.ylabel('Wing Area [m$^2$]')
    plt.title('Wing Area vs Cruise Altitude')
    plt.savefig('Wing_Sweeps/wing_alt_S.pdf')
    plt.show()

    plt.plot(solAltsweep('CruiseAlt'), solAltsweep('b'), '-r')
    plt.xlabel('Cruise Altitude [ft]')
    plt.ylabel('Wing Span [m]')
    plt.title('Wing Span vs Cruise Altitude')
    plt.savefig('Wing_Sweeps/wing_alt_b.pdf')
    plt.show()

    plt.plot(solAltsweep('CruiseAlt'), solAltsweep['sensitivities']['constants']['V_{ne}'], '-r')
    plt.xlabel('Cruise Altitude [ft]')
    plt.ylabel('Sensitivity of $V_{ne}$')
    plt.title('Sensitivity of $V_{ne}$ vs Cruise Altitude')
    plt.savefig('Wing_Sweeps/wing_alt_SensNve.pdf')
    plt.show()

    plt.plot(solAltsweep('CruiseAlt'), solAltsweep['sensitivities']['constants']['N_{lift}'], '-r')
    plt.xlabel('Cruise Altitude [ft]')
    plt.ylabel('Sensitivity to Seciontal Lift Multiplier')
    plt.title('Sensitivity to Seciontal Lift Multiplier vs Cruise Altitude')
    plt.savefig('Wing_Sweeps/wing_alt_SensNLift.pdf')
    plt.show()

    plt.plot(solAltsweep('CruiseAlt'), solAltsweep['sensitivities']['constants']['C_{L_{wmax}}'], '-r')
    plt.xlabel('Cruise Altitude [ft]')
    plt.ylabel('Sensitivity to Wing $C_{L_{max}}$')
    plt.title('Sensitivity to Wing $C_{L_{max}}$ vs Cruise Altitude')
    plt.savefig('Wing_Sweeps/wing_alt_SensClMax.pdf')
    plt.show()

    plt.plot(solAltsweep('CruiseAlt'), solAltsweep['sensitivities']['constants']['ReqRng'], '-r')
    plt.xlabel('Cruise Altitude [ft]')
    plt.ylabel('Sensitivity to Required Range')
    plt.title('Sensitivity to Required Range vs Cruise Altitude')
    plt.savefig('Wing_Sweeps/wing_alt_SensRng.pdf')
    plt.show()

    plt.plot(solAltsweep('CruiseAlt'), solAltsweep['sensitivities']['constants']['CruiseAlt'], '-r')
    plt.xlabel('Cruise Altitude [ft]')
    plt.ylabel('Sensitivity to Cruise Altitude')
    plt.title('Sensitivity to Cruise Altitude vs Cruise Altitude')
    plt.savefig('Wing_Sweeps/wing_alt_SensAlt.pdf')
    plt.show()
