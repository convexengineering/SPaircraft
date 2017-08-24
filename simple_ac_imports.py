"""
simple aircraft classes to import
"""
from gpkit import Model, Variable, units, SignomialsEnabled
from gpkit.constraints.sigeq import SignomialEquality as SignomialEquality
from gpkit.tools import te_exp_minus1
from gpkit.constraints.tight import Tight as TCS
from numpy import pi
import numpy as np
from turbofan.engine_validation import Engine
from Wing_simple_performance import Wing

class Aircraft(Model):
    "Aircraft class"
    def  setup(self, Nclimb, Ncruise, enginestate, eng, Nfleet=0, **kwargs):
        #create submodels
        self.fuse = Fuselage()
        self.wing = Wing()
        if Nfleet != 0:
            self.engine = Engine(0, True, Nclimb+Ncruise, enginestate, eng, Nfleet)
        else:
           self.engine = Engine(0, True, Nclimb+Ncruise, enginestate, eng)            

        #variable definitions
        numeng = Variable('numeng', '-', 'Number of Engines')
        Vne     = Variable('V_{ne}', 144, 'm/s', 'Never exceed velocity')
        rho0    = Variable('\\rho_0', 1.225, 'kg/m^3', 'Air density (0 ft)')

##        SMmin = Variable('SM_{min}', '-', 'Minimum Static Margin')
##        dxCG = Variable('\\Delta x_{CG}', 'm', 'Max CG Travel Range')

        
        constraints = [self.wing['x_w'] == self.fuse['l_{fuse}']*0.6,
                       ]

        self.components = [self.fuse, self.wing, self.engine]

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
    
    def cruise_climb_dynamic(self, state):
        """
        creates an aircraft cruise performance model, given a state
        """
        return CruiseClimbP(self, state)

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

        self.Pmodels = [self.wingP, self.fuseP]

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
        CD = Variable('C_{D}', '-', 'Overall Drag Coefficient')

##        xAC = Variable('x_{AC}', 'm', 'Aerodynamic Center Location')
##        xCG     = Variable('x_{CG}', 'm', 'CG location')

        constraints = []

        constraints.extend([
            #speed must be greater than stall speed
            state['V'] >= Vstall,

            #Figure out how to delete
            Vstall == 120*units('kts'),
            WLoadmax == 6664 * units('N/m^2'),

            #compute the drag
            TCS([D >= self.wingP['D_{wing}'] + self.fuseP['D_{fuse}']]),

            #compute the drag coefficient
            CD == D/(.5*state.atm['\\rho']*state['V']**2*self.aircraft['S']),

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
            self.wingP['L_w'] >= W_avg,
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

            self.aircraftP['W_{burn}'] == self.aircraft.engine['TSFC'][:2]*self.aircraft.engine['F'][:2]*self.aircraftP['thr']
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

             self.aircraftP['W_{burn}'] == self.aircraft.engine['TSFC'][2:]*self.aircraft.engine['F'][:2:]*self.aircraftP['thr']
             ])

        return constraints, self.aircraftP

class CruiseClimbP(Model):
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
        RngCruise = Variable('RngCruise', 'nautical_miles', 'Down Range Covered in Each Cruise Segment')

        #constraints
        constraints = []
        
        constraints.extend([ 
            RC == excessP/self.aircraftP['W_{avg}'],
            
            #make the small angle approximation and compute theta
            theta * state['V']  == RC,
           
            dhft == self.aircraftP['tmin'] * RC,
        
            #makes a small angle assumption during climb
            RngCruise == self.aircraftP['thr']*state['V'],
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

class CruiseClimbSegment(Model):
    """
    Combines a flight state and aircrat to form a cruise flight segment
    """
    def setup(self, aircraft, **kwargs):
        self.state = FlightState()
        self.cruiseP = aircraft.cruise_climb_dynamic(self.state)

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

##class Wing(Model):
##    """
##    place holder wing model
##    """
##    def setup(self, ** kwargs):
##        #new variables
##        W_wing = Variable('W_{wing}', 'N', 'Wing Weight')
##                           
##        #aircraft geometry
##        S = Variable('S', 'm^2', 'Wing Planform Area')
##        AR = Variable('AR', '-', 'Aspect Ratio')
##        span = Variable('b', 'm', 'Wing Span')
##        span_max = Variable('b_{max}', 'm', 'Max Wing Span')
##
##        K = Variable('K', '-', 'K for Parametric Drag Model')
##        e = Variable('e', '-', 'Oswald Span Efficiency Factor')
##
##        dum1 = Variable('dum1', 124.58, 'm^2')
##        dum2 = Variable('dum2', 105384.1524, 'N')
##
##        mac    = Variable('mac', 'm',
##                  'Mean aerodynamic chord (wing)')
##        
##        cmw = Variable('c_{m_{w}}', '-', 'Wing Pitching Moment Coefficient')
##
##        CLmax = Variable('C_{L_{max}}', '-', 'Max Wing Lift Coefficient')
##
##        xw     = Variable('x_w', 'm', 'Position of wing aerodynamic center')
##        
##        constraints = []
##
##        constraints.extend([
##            #wing weight constraint
##            #based off of a raymer weight and 737 data from TASOPT output file
##            (S/(dum1))**.65 * (AR/10.1)**.5 == W_wing/(dum2),
##
##            #compute wing span and aspect ratio, subject to a span constraint
##            AR == (span**2)/S,
##            AR <= 10,
##
##            #compute K for the aircraft
##            K == (pi * e * AR)**-1,
##
##            mac == mac,
##            CLmax == CLmax,
##            cmw == cmw,
##            xw == xw,
##            ])
##
##        return constraints
##
##    def dynamic(self, state):
##        """
##        creates an instance of the wing's performance model
##        """
##        return WingPerformance(self, state)
##        
##
##class WingPerformance(Model):
##    """
##    wing aero modeling
##    """
##    def setup(self, wing, state, **kwargs):
##        #new variables
##        CL= Variable('C_{L}', '-', 'Lift Coefficient')
##        Cdw = Variable('C_{d_w}', '-', 'Cd for a NC130 Airfoil at Re=2e7')
##        Dwing = Variable('D_{wing}', 'N', 'Total Wing Drag')
##        Lwing = Variable('L_{wing}', 'N', 'Wing Lift')
##
##        CLaw    = Variable('C_{L_{aw}}', '-', 'Lift curve slope, wing')
##
##        #constraints
##        constraints = []
##
##        constraints.extend([
##            #airfoil drag constraint
##            Lwing == (.5*wing['S']*state.atm['\\rho']*state['V']**2)*CL,
##            TCS([Cdw**6.5 >= (1.02458748e10 * CL**15.587947404823325 * state['M']**156.86410659495155 +
##                         2.85612227e-13 * CL**1.2774976672501526 * state['M']**6.2534328002723703 +
##                         2.08095341e-14 * CL**0.8825277088649582 * state['M']**0.0273667615730107 +
##                         1.94411925e+06 * CL**5.6547413360261691 * state['M']**146.51920742858428)]),
##            TCS([Dwing >= (.5*wing['S']*state.atm['\\rho']*state['V']**2)*(Cdw + wing['K']*CL**2)]),
##
##            CLaw == 5,
##            ])
##
##        return constraints

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
    def setup(self, fuse, state, **kwargs):
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
