"""Simple commercial aircraft flight profile and aircraft model"""
from numpy import pi
import numpy as np
from gpkit import Variable, Model, units, SignomialsEnabled, Vectorize
from gpkit.constraints.sigeq import SignomialEquality
from gpkit.tools import te_exp_minus1
from gpkit.constraints.tight import Tight as TCS
import matplotlib.pyplot as plt
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
    def setup(self, **kwargs):
        #create submodels
        self.fuse = Fuselage()
        self.wing = Wing()
        self.engine = Engine()
        self.VT = VerticalTail()

        #variable definitions
        numeng = Variable('numeng', '-', 'Number of Engines')

        constraints = []

        constraints.extend([
            numeng == numeng, #need numeng in the model
            ])

        self.components = [self.fuse, self.wing, self.engine, self.VT]

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
        self.engineP = aircraft.engine.dynamic(state)
        self.VTP = aircraft.VT.dynamic(aircraft.fuse, state)

        self.Pmodels = [self.wingP, self.fuseP, self.engineP, self.VTP]

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
            TCS([D >= self.wingP['D_{wing}'] + self.fuseP['D_{fuse}'] + self.VTP['D_{vt}']]),

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

        return self.Pmodels, constraints

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
        return constraints, self.atm, self.alt

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
        substitutions = None

        return constraints

class Engine(Model):
    """
    place holder engine model
    """
    def setup(self, **kwargs):
        #new variables
        W_engine = Variable('W_{engine}', 'N', 'Weight of a Single Turbofan Engine')
        A2 = Variable('A_2', 'm^2', 'Fan Area')
        
        constraints = []

        constraints.extend([
            W_engine == 1000 * units('N'),
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
        n_pax = Variable('n_{pax}', '-', 'Number of Passengers to Carry')
                           
        #weight variables
        W_payload = Variable('W_{payload}', 'N', 'Aircraft Payload Weight')
        W_e = Variable('W_{e}', 'N', 'Empty Weight of Aircraft')
        W_pax = Variable('W_{pax}', 'N', 'Estimated Average Passenger Weight, Includes Baggage')

        A_fuse = Variable('A_{fuse}', 'm^2', 'Estimated Fuselage Area')
        pax_area = Variable('pax_{area}', 'm^2', 'Estimated Fuselage Area per Passenger')

        lfuse   = Variable('l_{fuse}', 'm', 'Fuselage length')

        constraints = []
        
        constraints.extend([
            #compute fuselage area for drag approximation
            A_fuse == pax_area * n_pax,

            #constraints on the various weights
            W_payload == n_pax * W_pax,
            
            #estimate based on TASOPT 737 model
            W_e == .75*W_payload,

            lfuse == lfuse,
            
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
        xCG    = Variable('x_{CG}', 'm', 'x-location of CG')
        #constraints
        constraints = []

        constraints.extend([
            Dfuse == Cdfuse * (.5 * fuse['A_{fuse}'] * state.atm['\\rho'] * state['V']**2),

            Cdfuse == .005,

            xCG == 18*units('m'),
            ])

        return constraints

class Mission(Model):
    """
    mission class, links together all subclasses
    """
    def setup(self):
        #define the number of each flight segment
        Nclimb = 2
        Ncruise = 2

        #build required submodels
        ac = Aircraft()

        #vectorize
        with Vectorize(Nclimb):
            cls = ClimbSegment(ac)

        with Vectorize(Ncruise):
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
            TCS([ac['W_{e}'] + ac['W_{payload}'] + W_ftotal + ac['numeng'] * ac['W_{engine}'] + ac['W_{wing}'] + 2*ac.VT['W_{struct}'] <= W_total]),

            cls.climbP.aircraftP['W_{start}'][0] == W_total,
            cls.climbP.aircraftP['W_{end}'][-1] == crs.cruiseP.aircraftP['W_{start}'][0],

            # similar constraint 1
            TCS([cls.climbP.aircraftP['W_{start}'] >= cls.climbP.aircraftP['W_{end}'] + cls.climbP.aircraftP['W_{burn}']]),
            # similar constraint 2
            TCS([crs.cruiseP.aircraftP['W_{start}'] >= crs.cruiseP.aircraftP['W_{end}'] + crs.cruiseP.aircraftP['W_{burn}']]),

            cls.climbP.aircraftP['W_{start}'][1:] == cls.climbP.aircraftP['W_{end}'][:-1],
            crs.cruiseP.aircraftP['W_{start}'][1:] == crs.cruiseP.aircraftP['W_{end}'][:-1],

            TCS([ac['W_{e}'] + ac['W_{payload}'] + ac['numeng'] * ac['W_{engine}'] + ac['W_{wing}'] + 2*ac.VT['W_{struct}'] <= crs.cruiseP.aircraftP['W_{end}'][-1]]),

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
            ])

        #VT constriants
        constraints.extend([
            ac.VT['T_e'] == cls.climbP.engineP['thrust'][0],

            # Drag of a windmilling engine
            ac.VT['D_{wm}'] >= 0.5*ac.VT['\\rho_{TO}']*ac.VT['V_1']**2*ac.engine['A_2']*ac.VT['C_{D_{wm}}'],
            
            ac.VT['x_{CG_{vt}}'] <= ac.fuse['l_{fuse}'],
            
            #VTP constraints
            ac.fuse['l_{fuse}'] >= ac.VT['\\Delta x_{lead_v}'] + cls.climbP.aircraftP['x_{CG}'],
            ac.VT['x_{CG_{vt}}'] >= cls.climbP.aircraftP['x_{CG}']+(ac.VT['\\Delta x_{lead_v}']+ac.VT['\\Delta x_{trail_v}'])/2,

            ac.fuse['l_{fuse}'] >= ac.VT['\\Delta x_{lead_v}'] + crs.cruiseP.aircraftP['x_{CG}'],
            ac.VT['x_{CG_{vt}}'] >= crs.cruiseP.aircraftP['x_{CG}']+(ac.VT['\\Delta x_{lead_v}']+ac.VT['\\Delta x_{trail_v}'])/2,
            ])
        
        # Model.__init__(self, W_ftotal + s*units('N'), constraints + ac + cls + crs, subs)
        self.cost = W_ftotal
        return constraints, ac, crs, cls

class VerticalTail(Model):
    """
    Vertical tail sizing

    References:
    1: LEAP engine specs (1B)
    2: TASOPT 737 code
    3: Boeing 737 airport doc
    4: Boeing 737 Max airport doc
    5: http://www.digitaldutch.com/atmoscalc
    6: Engineering toolbox
    7: Boeing.com
    """
    def setup(self, **kwargs):
        self.vtns = VerticalTailNoStruct()
        self.wb = WingBox(self.vtns)

        return self.vtns, self.wb

    def dynamic(self, fuse, state):
        """"
        creates a horizontal tail performance model
        """
        return VerticalTailPerformance(self, fuse, state)

class VerticalTailNoStruct(Model):
    """
    Vertical tail sizing w/no structural model

    References:
    1: LEAP engine specs (1B)
    2: TASOPT 737 code
    3: Boeing 737 airport doc
    4: Boeing 737 Max airport doc
    5: http://www.digitaldutch.com/atmoscalc
    6: Engineering toolbox
    7: Boeing.com
    """
    def setup(self, **kwargs):
        #define new variables
        Avt    = Variable('A_{vt}', '-', 'Vertical tail aspect ratio')
        CDwm   = Variable('C_{D_{wm}}', '-', 'Windmill drag coefficient')
        Dwm    = Variable('D_{wm}', 'N', 'Engine out windmill drag')
        Lvmax  = Variable('L_{v_{max}}', 'N',
                          'Maximum load for structural sizing')
        CLvmax = Variable('C_{L_{vmax}}', '-', 'Max lift coefficient')
        CLvtEO   = Variable('C_{L_{vtEO}}', '-', 'Vertical tail lift coefficient (Engine Out)')
        clvtEO   = Variable('c_{l_{vtEO}}', '-',
                            'Sectional lift force coefficient (engine out)')
        LvtEO    = Variable('L_{vtEO}', 'N', 'Vertical tail lift in engine out')
        Svt    = Variable('S_{vt}', 'm^2', 'Vertical tail reference area (half)')
        V1     = Variable('V_1', 'm/s', 'Minimum takeoff velocity')
        Vne    = Variable('V_{ne}', 144, 'm/s', 'Never exceed velocity')
        Wvt    = Variable('W_{vt}', 'N', 'Vertical tail weight')
        bvt    = Variable('b_{vt}', 'm', 'Vertical tail half span')
        cma    = Variable('\\bar{c}_{vt}', 'm', 'Vertical tail mean aero chord')
        croot  = Variable('c_{root_{vt}}', 'm', 'Vertical tail root chord')
        ctip   = Variable('c_{tip_{vt}}', 'm', 'Vertical tail tip chord')
        dxlead = Variable('\\Delta x_{lead_v}', 'm',
                          'Distance from CG to vertical tail leading edge')
        dxtrail= Variable('\\Delta x_{trail_v}', 'm',
                          'Distance from CG to vertical tail trailing edge')
        e      = Variable('e_v', '-', 'Span efficiency of vertical tail')
        lvt    = Variable('l_{vt}', 'm', 'Vertical tail moment arm')
        mu0    = Variable('\\mu_0', 1.8E-5, 'N*s/m^2', 'Dynamic viscosity (SL)')
        p      = Variable('p_{vt}', '-', 'Substituted variable = 1 + 2*taper')
        plamv  = Variable('p_{\\lambda_v}', '-',
                          'Dummy variable = 1 + 2\\lambda') # fuselage
        q      = Variable('q_{vt}', '-', 'Substituted variable = 1 + taper')
        rho0   = Variable('\\rho_{TO}', 'kg/m^3', 'Air density (SL))')
        rho_c  = Variable('\\rho_c', 'kg/m^3', 'Air density (35,000ft)')
        tanL   = Variable('\\tan(\\Lambda_{vt})', '-',
                          'Tangent of leading edge sweep (40 deg)')
        taper  = Variable('\\lambda_{vt}', '-', 'Vertical tail taper ratio')
        tau    = Variable('\\tau_{vt}', '-', 'Vertical tail thickness/chord ratio')
        xCGvt  = Variable('x_{CG_{vt}}', 'm', 'x-location of tail CG')
        y_eng  = Variable('y_{eng}', 'm', 'Engine moment arm')
        zmac   = Variable('z_{\\bar{c}_{vt}}', 'm',
                          'Vertical location of mean aerodynamic chord')

        #engine values
        Te     = Variable('T_e', 'N', 'Thrust per engine at takeoff')

        #variables specific to yaw rate sizing
        Vland = Variable('V_{land}', 'm/s', 'Aricraft Landing Speed')
        CLvyaw = Variable('C_{L_{vyaw}}', '-', 'VT CL at rotation')
        Iz = Variable('I_{z}', 'kg*m^2', 'Aircrat Z-axis Moment of Inertia')
        rreq = Variable('\\dot{r}_{req}', 's^-2', 'Required Yaw Rate at Landing')
        
        #constraints
        constraints = []

        with SignomialsEnabled():
            #non vectorized constraints
            constraints.extend([
               #constraint tail Cl at flare
                CLvyaw == .85*CLvmax,

                #meet yaw rate constraint at flare
                2*.5*rho0*Vland**2*Svt*lvt*CLvyaw >= rreq*Iz,
                
                2*LvtEO*lvt >= Te*y_eng + Dwm*y_eng,
                # Force moment balance for one engine out condition
                # TASOPT 2.0 p45

                TCS([dxlead + zmac*tanL + 0.25*cma >= lvt]), # [SP]
                # Tail moment arm

                LvtEO == 0.5*rho0*V1**2*Svt*CLvtEO,
                # Vertical tail force (y-direction) for engine out

                TCS([CLvtEO*(1 + clvtEO/(np.pi*e*Avt)) <= clvtEO]),
                #engine out CL computation

                Avt == bvt**2/Svt,
                          
                Svt <= bvt*(croot + ctip)/2, # [SP]
                # Tail geometry relationship

                TCS([dxtrail >= croot + dxlead]),
                # Tail geometry constraint

                
                # Fuselage length constrains the tail trailing edge

                TCS([p >= 1 + 2*taper]),
                TCS([2*q >= 1 + p]),
                zmac == (bvt/3)*q/p,
                TCS([(2./3)*(1 + taper + taper**2)*croot/q >= cma]), # [SP]
                taper == ctip/croot,
                # Define vertical tail geometry

                Lvmax == 0.5*rho0*Vne**2*Svt*CLvmax,
                #compute the max force

                taper >= 0.27,
                # TODO: Constrain taper by tip Reynolds number
                # source: b737.org.uk
                ])

        return constraints

class VerticalTailPerformance(Model):
    """
    Vertical tail perofrmance model
    """
    def setup(self, vt, fuse, state):
        self.fuse = fuse
        self.vt = vt

        #define new variables
 
        CLvt   = Variable('C_{L_{vt}}', '-', 'Vertical tail lift coefficient')
        clvt   = Variable('c_{l_{vt}}', '-',
                          'Sectional lift force coefficient')
        Dvt    = Variable('D_{vt}', 'N', 'Vertical tail viscous drag, cruise')
        Rec    = Variable('Re_{vt}', '-', 'Vertical tail reynolds number, cruise')
        Vinf   = Variable('V_{\\infty}', 'm/s', 'Cruise velocity')
        CDvis  = Variable('C_{D_{vis}}', '-', 'Viscous drag coefficient')

        #constraints
        constraints = []

         #vectorized constraints
        constraints.extend([
                    
##            TCS([CLvt*(1 + clvt/(np.pi*self.vt['e_v']*self.vt['A_{vt}'])) <= clvt]),

            # Finite wing theory
            # people.clarkson.edu/~pmarzocc/AE429/AE-429-4.pdf
            # Valid because tail is untwisted and uncambered
            # (lift curve slope passes through origin)

            .5*Dvt >= 0.5*state['\\rho']*state['V']**2*self.vt['S_{vt}']*CDvis,
            CDvis**0.125 >= 0.19*(self.vt['\\tau_{vt}'])**0.0075 *(Rec)**0.0017
                        + 1.83e+04*(self.vt['\\tau_{vt}'])**3.54*(Rec)**-0.494
                        + 0.118*(self.vt['\\tau_{vt}'])**0.0082 *(Rec)**0.00165
                        + 0.198*(self.vt['\\tau_{vt}'])**0.00774*(Rec)**0.00168,
            # Vertical tail viscous drag in cruise
            # Data fit from Xfoil

            Rec == state.atm['\\rho']*state['V']*self.vt['\\bar{c}_{vt}']/state.atm['\\mu'],
            # Cruise Reynolds number


            ])

        return constraints

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

        numspar = Variable('N_{spar}', '-', 'Number of Spars in Each VT Carrying Stress in 1 in 20 Case')

        objective = Wstruct

        if isinstance(surface, VerticalTailNoStruct):
            #factors of 2 required since the VT is only half span and the wing model
            #is for a full span wing
            AR = 2*surface['A_{vt}']
            b = 2*surface['b_{vt}']
            S = 2*surface['S_{vt}']
            p = surface['p_{vt}']
            q = surface['q_{vt}']
            tau = surface['\\tau_{vt}']
            Lmax = 2*surface['L_{v_{max}}']
            taper = surface['\\lambda_{vt}']

        constraints = [
                       # Aspect ratio definition
                       AR == b**2/S,

                       # Defining taper dummy variables
                       TCS([p >= 1 + 2*taper]),
                       TCS([2*q >= 1 + p]),

                       # Upper bound on maximum thickness
                       tau <= 0.15,

                       # Root moment calculation (see Hoburg 2014)
                       # Depends on a given load the wing must support, Lmax
                       # Assumes lift per unit span proportional to local chord
                       numspar*Mr >= Lmax*AR*p/24,

                       # Root stiffness (see Hoburg 2014)
                       # Assumes rh = 0.75, so that rms box height = ~0.92*tmax
                       0.92*w*tau*tcap**2 + Icap <= 0.92**2/2*w*tau**2*tcap,

                       # Stress limit
                       # Assumes bending stress carried by caps (Icap >> Iweb)
                       8 >= Nlift*Mr*AR*q**2*tau/(S*Icap*sigmax),

                       # Shear web sizing
                       # Assumes all shear loads are carried by web and rh=0.75
                       12 >= AR*Lmax*Nlift*q**2/(tau*S*tweb*sigmaxshear*numspar),

                       # Posynomial approximation of nu=(1+lam+lam^2)/(1+lam^2)
                       nu**3.94 >= 0.86*p**(-2.38)+ 0.14*p**0.56,

                       # Weight of spar caps and shear webs
                       Wcap >= 8*rhocap*g*w*tcap*S**1.5*nu/(3*AR**0.5),
                       Wweb >= 8*rhoweb*g*rh*tau*tweb*S**1.5*nu/(3*AR**0.5),

                       # Total wing weight using an additional weight fraction
                       Wstruct/3 >= (1 + fwadd)*(Wweb + Wcap),
                       ]
        
        return constraints

if __name__ == '__main__':
    plot = True
    
    substitutions = {      
            'ReqRng': 500,
            'CruiseAlt': 30000,
            'numeng': 2,
##            'W_{Load_max}': 6664,
            'W_{pax}': 91 * 9.81,
            'n_{pax}': 150,
            'pax_{area}': 1,
##            'C_{D_{fuse}}': .005, #assumes flat plate turbulent flow, from wikipedia
            'e': .9,
            'b_{max}': 35,

            #VT subs
           'C_{D_{wm}}': 0.5, # [2]
           'C_{L_{vmax}}': 2.6, # [2]
           'V_1': 70,
           'V_{ne}': 144, # [2]
           '\\rho_{TO}': 1.225,
           '\\tan(\\Lambda_{vt})': np.tan(40*np.pi/180),
##           'c_{l_{vt}}': 0.5, # [2]
           'c_{l_{vtEO}}': 0.5,
           'A_2': np.pi*(.5*1.75)**2, # [1]
           'e_v': 0.8,
           'l_{fuse}': 39,
##           'x_{CG}': 18,
           'y_{eng}': 4.83, # [3]

           'V_{land}': 72,
           'I_{z}': 12495000, #estimate for late model 737 at max takeoff weight (m l^2/12)
           '\\dot{r}_{req}': 0.174533, #10 deg/s yaw rate

            'N_{spar}': 2,
            }

    if plot == True:
        mission = Mission()
        m = Model(mission['W_{f_{total}}'], mission)
        m.substitutions.update(substitutions)
        sol = m.localsolve(solver='mosek', verbosity = 4)

    #     substitutions = {
    #            'ReqRng': ('sweep', np.linspace(500,3000,8)),
    #             'CruiseAlt': 30000,
    #             'numeng': 2,
    # ##            'W_{Load_max}': 6664,
    #             'W_{pax}': 91 * 9.81,
    #             'n_{pax}': 150,
    #             'pax_{area}': 1,
    # ##            'C_{D_{fuse}}': .005, #assumes flat plate turbulent flow, from wikipedia
    #             'e': .9,
    #             'b_{max}': 35,
    #
    #             #VT subs
    #            'C_{D_{wm}}': 0.5, # [2]
    #            'C_{L_{vmax}}': 2.6, # [2]
    #            'V_1': 70,
    #            'V_{ne}': 144, # [2]
    #            '\\rho_{TO}': 1.225,
    #            '\\tan(\\Lambda_{vt})': np.tan(40*np.pi/180),
    # ##           'c_{l_{vt}}': 0.5, # [2]
    #            'c_{l_{vtEO}}': 0.5,
    #            'A_2': np.pi*(.5*1.75)**2, # [1]
    #            'e_v': 0.8,
    #            'l_{fuse}': 39,
    # ##           'x_{CG}': 18,
    #            'y_{eng}': 4.83, # [3]
    #
    #            'V_{land}': 72,
    #            'I_{z}': 12495000, #estimate for late model 737 at max takeoff weight (m l^2/12)
    #            '\\dot{r}_{req}': 0.174533, #10 deg/s yaw rate
    #
    #            'N_{spar}': 2,
    #             }
    #
    #     mission = Mission()
    #     m = Model(mission['W_{f_{total}}'], mission)
    #     m.substitutions.update(substitutions)
    #     solRsweep = m.localsolve(solver='mosek', verbosity = 4)
    #
    #     RC = []
    #
    #     for i in range(len(solRsweep('ReqRng'))):
    #         RC.append(mag(solRsweep('RC')[i][0]))
    #
    #     plt.plot(solRsweep('ReqRng'), RC, '-r')
    #     plt.xlabel('Mission Range [nm]')
    #     plt.ylabel('Initial  Rate of Climb [ft/min]')
    #     plt.title('Initial Rate of Climb vs Range')
    # ##    plt.savefig('HT_Sweeps/VT_rng_RC.pdf')
    #     plt.show()
    #
    # ##    plt.plot(solRsweep('ReqRng'), solRsweep('L_{vtEO}'), '-r')
    # ##    plt.xlabel('Mission Range [nm]')
    # ##    plt.ylabel('VT Lift and Takeoff ENgine Out [N]')
    # ##    plt.title('Initial Climb Thrust vs Range')
    # ####    plt.savefig('HT_Sweeps/VT_rng_LVTTO.pdf')
    # ##    plt.show()
    #
    #     plt.plot(solRsweep('ReqRng'), solRsweep('L_{v_{max}}'), '-r')
    #     plt.xlabel('Mission Range [nm]')
    #     plt.ylabel('Max VT Lift Force [N]')
    #     plt.title('Initial Climb Thrust vs Range')
    # ##    plt.savefig('HT_Sweeps/VT_rng_LVmax.pdf')
    #     plt.show()
    #
    # ##    plt.plot(solRsweep('ReqRng'), solRsweep('T_e'), '-r')
    # ##    plt.xlabel('Mission Range [nm]')
    # ##    plt.ylabel('Initial  Climb Thrust [N]')
    # ##    plt.title('Initial Climb Thrust vs Range')
    # ####    plt.savefig('HT_Sweeps/VT_rng_TTO.pdf')
    # ##    plt.show()
    #
    #     plt.plot(solRsweep('ReqRng'), solRsweep('W_{struct}'), '-r')
    #     plt.xlabel('Mission Range [nm]')
    #     plt.ylabel('Vertical Tail Weight [N]')
    #     plt.title('Vertical Tail Weight vs Range')
    # ##    plt.savefig('HT_Sweeps/VT_rng_VTweight.pdf')
    #     plt.show()
    #
    #     plt.plot(solRsweep('ReqRng'), solRsweep('A_{vt}'), '-r')
    #     plt.xlabel('Mission Range [nm]')
    #     plt.ylabel('Vertical Tail Aspect Ratio')
    #     plt.title('Vertical Tail Aspect Ratio vs Range')
    # ##    plt.savefig('HT_Sweeps/VT_rng_VTAR.pdf')
    #     plt.show()
    #
    #     plt.plot(solRsweep('ReqRng'), solRsweep('S_{vt}'), '-r')
    #     plt.xlabel('Mission Range [nm]')
    #     plt.ylabel('Vertical Tail Area [m$^2$]')
    #     plt.title('Vertical Tail Area vs Range')
    # ##    plt.savefig('HT_Sweeps/VT_rng_SVT.pdf')
    #     plt.show()
    #
    #     plt.plot(solRsweep('ReqRng'), solRsweep['sensitivities']['constants']['ReqRng'], '-r')
    #     plt.xlabel('Mission Range [nm]')
    #     plt.ylabel('Sensitivity of Mission Range')
    #     plt.title('Sensitivity to Range vs Range')
    # ####    plt.savefig('HT_Sweeps/VT_rng_RngSens.pdf')
    #     plt.show()
    #
    #     plt.plot(solRsweep('ReqRng'), solRsweep['sensitivities']['constants']['CruiseAlt'], '-r')
    #     plt.xlabel('Mission Range [nm]')
    #     plt.ylabel('Sensitivity of Cruise Altitude')
    #     plt.title('Sensitivity to Cruise Altitude vs Range')
    # ##    plt.savefig('HT_Sweeps/VT_rng_AltSens.pdf')
    #     plt.show()
    #
    #     plt.plot(solRsweep('ReqRng'), solRsweep['sensitivities']['constants']['b_{max}'], '-r')
    #     plt.xlabel('Mission Range [nm]')
    #     plt.ylabel('Sensitivity of Max Wing Span')
    #     plt.title('Sensitivity of Max Wing Span vs Range')
    # ##    plt.savefig('HT_Sweeps/VT_rng_bMaxSens.pdf')
    #     plt.show()
    #
    #     plt.plot(solRsweep('ReqRng'), solRsweep['sensitivities']['constants']['W_{pax}'], '-r')
    #     plt.xlabel('Mission Range [nm]')
    #     plt.ylabel('Sensitivity of Passenger Weight')
    #     plt.title('Sensitivity of Passegner Weight vs Range')
    # ##    plt.savefig('HT_Sweeps/VT_rng_WPaxSens.pdf')
    #     plt.show()
    #
    #     substitutions = {
    #             'ReqRng': 2000,
    #             'CruiseAlt': ('sweep', np.linspace(20000,40000,10)),
    #             'numeng': 2,
    # ##            'W_{Load_max}': 6664,
    #             'W_{pax}': 91 * 9.81,
    #             'n_{pax}': 150,
    #             'pax_{area}': 1,
    # ##            'C_{D_{fuse}}': .005, #assumes flat plate turbulent flow, from wikipedia
    #             'e': .9,
    #             'b_{max}': 35,
    #
    #             #VT subs
    #            'C_{D_{wm}}': 0.5, # [2]
    #            'C_{L_{vmax}}': 2.6, # [2]
    #            'V_1': 70,
    #            'V_{ne}': 144, # [2]
    #            '\\rho_{TO}': 1.225,
    #            '\\tan(\\Lambda_{vt})': np.tan(40*np.pi/180),
    # ##           'c_{l_{vt}}': 0.5, # [2]
    #            'c_{l_{vtEO}}': 0.5,
    #            'A_2': np.pi*(.5*1.75)**2, # [1]
    #            'e_v': 0.8,
    #            'l_{fuse}': 39,
    # ##           'x_{CG}': 18,
    #            'y_{eng}': 4.83, # [3]
    #
    #            'V_{land}': 72,
    #            'I_{z}': 12495000, #estimate for late model 737 at max takeoff weight (m l^2/12)
    #            '\\dot{r}_{req}': 0.174533, #10 deg/s yaw rate
    #
    #             'N_{spar}': 2,
    #             }
    #
    #     mission = Mission()
    #     m = Model(mission['W_{f_{total}}'], mission, substitutions)
    #     solAltsweep = m.localsolve(solver='mosek', verbosity = 4)
    #
    #     RC = []
    #
    #     for i in range(len(solAltsweep('CruiseAlt'))):
    #         RC.append(mag(solAltsweep('RC')[i][0]))
    #
    #     plt.plot(solAltsweep('CruiseAlt'), RC, '-r')
    #     plt.xlabel('Mission Range [nm]')
    #     plt.ylabel('Initial  Rate of Climb [ft/min]')
    #     plt.title('Initial Rate of Climb vs Range')
    # ##    plt.savefig('HT_Sweeps/VT_alt_RC.pdf')
    #     plt.show()
    #
    # ##    plt.plot(solAltsweep('CruiseAlt'), solAltsweep('L_{vtEO}'), '-r')
    # ##    plt.xlabel('Cruise Alt [ft]')
    # ##    plt.ylabel('Initial  Climb Thrust [N]')
    # ##    plt.title('Initial Climb Thrust vs Range')
    # ####    plt.savefig('HT_Sweeps/VT_alt_LvEOmax.pdf')
    # ##    plt.show()
    #
    #     plt.plot(solAltsweep('CruiseAlt'), solAltsweep('L_{v_{max}}'), '-r')
    #     plt.xlabel('Cruise Alt [ft]')
    #     plt.ylabel('Initial  Climb Thrust [N]')
    #     plt.title('Initial Climb Thrust vs Range')
    # ##    plt.savefig('HT_Sweeps/VT_alt_Lvmax.pdf')
    #     plt.show()
    #
    # ##    plt.plot(solAltsweep('CruiseAlt'), solAltsweep('T_e'), '-r')
    # ##    plt.xlabel('Cruise Altitude [ft]')
    # ##    plt.ylabel('Initial  Climb Thrust [N]')
    # ##    plt.title('Initial Climb Thrust vs Cruise Altitude')
    # ####    plt.savefig('HT_Sweeps/VT_alt_TTO.pdf')
    # ##    plt.show()
    #
    #     plt.plot(solAltsweep('CruiseAlt'), solAltsweep('W_{struct}'), '-r')
    #     plt.xlabel('Cruise Altitude [ft]')
    #     plt.ylabel('Vertical Tail Weight [N]')
    #     plt.title('Vertical Tail Weight vs Cruise Altitude')
    # ##    plt.savefig('HT_Sweeps/VT_alt_VTweight.pdf')
    #     plt.show()
    #
    #     plt.plot(solAltsweep('CruiseAlt'), solAltsweep('A_{vt}'), '-r')
    #     plt.xlabel('Cruise Altitude [ft]')
    #     plt.ylabel('Vertical Tail Aspect Ratio')
    #     plt.title('Vertical Tail Aspect Ratio vs Cruise Altitude')
    # ##    plt.savefig('HT_Sweeps/VT_alt_VTAR.pdf')
    #     plt.show()
    #
    #     plt.plot(solAltsweep('CruiseAlt'), solAltsweep('S_{vt}'), '-r')
    #     plt.xlabel('Cruise Altitude [ft]')
    #     plt.ylabel('Vertical Tail Area [m$^2$]')
    #     plt.title('Vertical Tail Area vs Cruise Altitude')
    # ##    plt.savefig('HT_Sweeps/VT_alt_VTarea.pdf')
    #     plt.show()
    #
    #     plt.plot(solAltsweep('CruiseAlt'), solAltsweep['sensitivities']['constants']['ReqRng'], '-r')
    #     plt.xlabel('Cruise Altitude [ft]')
    #     plt.ylabel('Sensitivity of Mission Range')
    #     plt.title('Sensitivity to Range vs Cruise Altitude')
    # ##    plt.savefig('HT_Sweeps/VT_alt_RngSens.pdf')
    #     plt.show()
    #
    #     plt.plot(solAltsweep('CruiseAlt'), solAltsweep['sensitivities']['constants']['CruiseAlt'], '-r')
    #     plt.xlabel('Cruise Altitude [ft]')
    #     plt.ylabel('Sensitivity of Cruise Altitude')
    #     plt.title('Sensitivity to Cruise Alt vs Cruise Alt')
    # ##    plt.savefig('HT_Sweeps/VT_alt_AltSens.pdf')
    #     plt.show()
    #
    #     plt.plot(solAltsweep('CruiseAlt'), solAltsweep['sensitivities']['constants']['W_{pax}'], '-r')
    #     plt.xlabel('Cruise Altitude [ft]')
    #     plt.ylabel('Sensitivity of Mission Range')
    #     plt.title('Sensitivity to Passenger Weight vs Cruise Altitude')
    # ##    plt.savefig('HT_Sweeps/VT_alt_WPaxSens.pdf')
    #     plt.show()
    #
    #     plt.plot(solAltsweep('CruiseAlt'), solAltsweep['sensitivities']['constants']['b_{max}'], '-r')
    #     plt.xlabel('Cruise Altitude [ft]')
    #     plt.ylabel('Sensitivity of Max Wing Span')
    #     plt.title('Sensitivity to Max Wing Span vs Cruise Alt')
    # ##    plt.savefig('HT_Sweeps/VT_alt_bMaxSens.pdf')
    #     plt.show()
    #
    #     substitutions = {
    #             'ReqRng': 2000,
    #             'CruiseAlt': ('sweep', np.linspace(20000,40000,10)),
    #             'numeng': 2,
    # ##            'W_{Load_max}': 6664,
    #             'W_{pax}': 91 * 9.81,
    #             'n_{pax}': 150,
    #             'pax_{area}': 1,
    # ##            'C_{D_{fuse}}': .005, #assumes flat plate turbulent flow, from wikipedia
    #             'e': .9,
    #             'b_{max}': 35,
    #
    #             #VT subs
    #            'C_{D_{wm}}': 0.5, # [2]
    #            'C_{L_{vmax}}': 2.6, # [2]
    #            'V_1': 70,
    #            'V_{ne}': 144, # [2]
    #            '\\rho_{TO}': 1.225,
    #            '\\tan(\\Lambda_{vt})': np.tan(40*np.pi/180),
    # ##           'c_{l_{vt}}': 0.5, # [2]
    #            'c_{l_{vtEO}}': 0.5,
    #            'A_2': np.pi*(.5*1.75)**2, # [1]
    #            'e_v': 0.8,
    #            'l_{fuse}': 39,
    # ##           'x_{CG}': 18,
    #            'y_{eng}': 4.83, # [3]
    #
    #            'V_{land}': 72,
    #            'I_{z}': 12495000, #estimate for late model 737 at max takeoff weight (m l^2/12)
    #            '\\dot{r}_{req}': ('sweep', np.linspace(.01,.25,10)), #10 deg/s yaw rate
    #
    #             'N_{spar}': 2,
    #             }
    #
    #     mission = Mission()
    #     m = Model(mission['W_{f_{total}}'], mission, substitutions)
    #     solYawsweep = m.localsolve(solver='mosek', verbosity = 4)
    #
    #     plt.plot(solYawsweep('\\dot{r}_{req}'), solYawsweep('W_{struct}'), '-r')
    #     plt.xlabel('Yaw Rate Acceleration [rad/s&^2$]')
    #     plt.ylabel('Vertical Tail Weight [N]')
    #     plt.title('Vertical Tail Weight vs Yaw Rate Acceleration')
    #     plt.savefig('VT_Sweeps/VT_r_VTweight.pdf')
    #     plt.show()
    #
    #     plt.plot(solYawsweep('\\dot{r}_{req}'), solYawsweep('A_{vt}'), '-r')
    #     plt.xlabel('Yaw Rate Acceleration [rad/s&^2$]')
    #     plt.ylabel('Vertical Tail Aspect Ratio')
    #     plt.title('Vertical Tail Aspect Ratio vs Yaw Rate Acceleration')
    #     plt.savefig('VT_Sweeps/VT_r_VTAR.pdf')
    #     plt.show()
    #
    #     plt.plot(solYawsweep('\\dot{r}_{req}'), solYawsweep('S_{vt}'), '-r')
    #     plt.xlabel('Yaw Rate Acceleration [rad/s&^2$]')
    #     plt.ylabel('Vertical Tail Area [m$^2$]')
    #     plt.title('Vertical Tail Area vs Yaw Rate Acceleration')
    #     plt.savefig('VT_Sweeps/VT_r_VTarea.pdf')
    #     plt.show()
