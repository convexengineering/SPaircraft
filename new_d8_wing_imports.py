"""
simple aircraft classes to import
"""
from gpkit import Model, Variable, units, SignomialsEnabled
from gpkit.constraints.sigeq import SignomialEquality as SignomialEquality
from gpkit.tools import te_exp_minus1
from gpkit.constraints.tight import Tight as TCS
from numpy import pi
import numpy as np
from engine_validation import Engine

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
            self.wingP['L_w'] == W_avg,
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

        return constraints


class Wing(Model):
    """
    Philippe's thesis wing model
    """
    def setup(self, **kwargs):
        self.wns = WingNoStruct()
        self.wb = WingBox(self.wns)

        constraints = [self.wns['\\lambda'] == self.wb['taper']]

        return self.wns, self.wb, constraints
        
    def dynamic(self, state):
        """
        returns an instance of the wing perofrmance model
        """
        return WingPerformance(self, state)

class WingNoStruct(Model):
    """
    place holder wing model
    """
    def setup(self, ** kwargs):
        #new variables
##        W_wing = Variable('W_{wing}', 'N', 'Wing Weight')
                           
        #aircraft geometry
        xw     = Variable('x_w', 'm', 'Position of wing aerodynamic center')
        S = Variable('S', 'm^2', 'Wing Planform Area')
        AR = Variable('AR', '-', 'Aspect Ratio')
        span = Variable('b', 'm', 'Wing Span')

        K = Variable('K', '-', 'K for Parametric Drag Model')
        e = Variable('e', '-', 'Oswald Span Efficiency Factor')

        amax    = Variable('\\alpha_{max,w}', '-', 'Max angle of attack')
        eta     = Variable('\\eta', '-', 'Lift efficiency (diff b/w sectional, actual lift)')
        tanL    = Variable('\\tan(\\Lambda)', '-',
                           'Tangent of quarter-chord sweep angle')
        croot   = Variable('c_{root}', 'm', 'Wing root chord')
        ctip    = Variable('c_{tip}', 'm', 'Wing tip chord')
        taper   = Variable('\\lambda', '-', 'Wing taper ratio')
        fl      = Variable('f(\\lambda_w)', '-', 'Empirical efficiency function of taper')
        WfuelWing   = Variable('W_{fuel_{wing}}', 'N', 'Fuel weight')
        Vfuel   = Variable('V_{fuel, max}', 'm^3', 'Available fuel volume')
        tau     = Variable('\\tau', '-', 'Wing thickness/chord ratio')
        rhofuel = Variable('\\rho_{fuel}', 'kg/m^3', 'Density of fuel')
        wwn       = Variable('wwn', 0.5, '-', 'Wingbox-width-to-chord ratio')
        g       = Variable('g', 9.81, 'm/s^2', 'Gravitational acceleration')
        Lmax    = Variable('L_{max}', 'N', 'Maximum load')

        mac    = Variable('mac', 'm',
                  'Mean aerodynamic chord (wing)')
        ymac    = Variable('y_{mac}', 'm',
                   'Spanwise location of mean aerodynamic chord')
        p       = Variable('p', '-', 'Substituted variable = 1 + 2*taper')
        q       = Variable('q', '-', 'Substituted variable = 1 + taper')
        CLwmax  = Variable('C_{L_{wmax}}', '-', 'Max lift coefficient, wing')

        constraints = []

        constraints.extend([
            taper == ctip/croot,

            # Oswald efficiency
            # Nita, Scholz, "Estimating the Oswald factor from
            # basic aircraft geometrical parameters"
            TCS([e*(1 + fl*AR) <= 1]),
            taper >= 0.2, # TODO


            #compute wing span and aspect ratio, subject to a span constraint
##            AR == (span**2)/S,
##            AR <= 10,

            #compute K for the aircraft
            K == (pi * e * AR)**-1,

            Vfuel <= croot*span*taper*wwn*tau*croot,
            WfuelWing <= rhofuel*Vfuel*g,

             p >= 1 + 2*taper,
             2*q >= 1 + p,

            p<=10,
            q<=10,

            ymac == (span/3)*q/p,
             TCS([(2./3)*(1+taper+taper**2)*croot/q <= mac],
                               reltol=1E-2),

            xw==xw,
            ])

        with SignomialsEnabled():
            constraints.extend([
                SignomialEquality(S, span*(croot + ctip)/2),

                # Oswald efficiency
                # Nita, Scholz, "Estimating the Oswald factor from
                # basic aircraft geometrical parameters"
                TCS([fl >= 0.0524*taper**4 - 0.15*taper**3
                     + 0.1659*taper**2 - 0.0706*taper + 0.0119],
                reltol=1E-2),
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
        Lwing = Variable('L_w', 'N', 'Wing Lift')

        alpha   = Variable('\\alpha_w', '-', 'Wing angle of attack')
        CLaw    = Variable('C_{L_{aw}}', '-', 'Lift curve slope, wing')

        #wing moment variables -- need a good way to model this, currently using TAT
        cmw = Variable('c_{m_{w}}', '-', 'Wing Pitching Moment Coefficient')
        
        Re      = Variable('Re_w', '-', 'Reynolds number (wing)')

        #constraints
        constraints = []

        constraints.extend([
            #airfoil drag constraint
            Lwing == (.5*wing['S']*state.atm['\\rho']*state['V']**2)*CL,
##            TCS([Cdw**6.5 >= (1.02458748e10 * CL**15.587947404823325 * state['M']**156.86410659495155 +
##                         2.85612227e-13 * CL**1.2774976672501526 * state['M']**6.2534328002723703 +
##                         2.08095341e-14 * CL**0.8825277088649582 * state['M']**0.0273667615730107 +
##                         1.94411925e+06 * CL**5.6547413360261691 * state['M']**146.51920742858428)]),
            TCS([Dwing >= (.5*wing['S']*state.atm['\\rho']*state['V']**2)*(Cdw + wing['K']*CL**2)]),

                Re == state['\\rho']*state['V']*wing['mac']/state['\\mu'],
                1 >= (2.56*CL**5.88/(Re**1.54*wing['\\tau']**3.32*Cdw**2.62)
                   + 3.8e-9*wing['\\tau']**6.23/(CL**0.92*Re**1.38*Cdw*9.57)
                   + 2.2e-3*Re**0.14*wing['\\tau']**0.033/(CL**0.01*Cdw**0.73)
                   + 6.14e-6*CL**6.53/(Re**0.99*wing['\\tau']**0.52*Cdw*5.19)
                   + 1.19e4*CL**9.78*wing['\\tau']**1.76/(Re*Cdw**0.91)),

            CL == CLaw*alpha,
            alpha <= wing['\\alpha_{max,w}'],
            ])

        with SignomialsEnabled():
            constraints.extend([
               SignomialEquality((wing['AR']/wing['\\eta'])**2*(1 + wing['\\tan(\\Lambda)']**2 - state['M']**2) + 8*pi*wing['AR']/CLaw
                  , (2*pi*wing['AR']/CLaw)**2),
                ])

        return constraints

class WingBox(Model):
    """
    Structural model for a wing
    source: Hoburg, "Geometric Programming for Aircraft Design Optimization"

    Note - does not have a performance model
    """

    def setup(self, surface, **kwargs):
        # Variables
        g = 9.81*units('m*s^-2')
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
                          'Wing added weight fraction') # [TAS]
        Nlift  = Variable('N_{lift}', 3.0, '-', 'Wing loading multiplier') # [TAS]
        rh     = Variable('r_h', 0.75, '-',
                          'Fractional wing thickness at spar web') # [TAS]
        rhocap = Variable('\\rho_{cap}', 2700, 'kg/m^3',
                          'Density of spar cap material') # [TAS]
        rhoweb = Variable('\\rho_{web}', 2700, 'kg/m^3',
                          'Density of shear web material') # [TAS]
        sigmax = Variable('\\sigma_{max}', 250e6, 'Pa',
                          'Allowable tensile stress') # [TAS]
        sigmaxshear = Variable('\\sigma_{max,shear}', 167e6, 'Pa',
                               'Allowable shear stress')
        wwb      = Variable('wwb', 0.5, '-', 'Wingbox-width-to-chord ratio') # [TAS]
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

        dum1 = Variable('dum1', 124.58, 'm^2')
        dum2 = Variable('dum2', 105384.1524, 'N')

        constraints = [
                       # Aspect ratio definition
                       AR == b**2/S,

                       # Upper bound on maximum thickness
                       tau <= 0.15,

                       # Root moment calculation (see Hoburg 2014)
                       # Depends on a given load the wing must support, Lmax
                       # Assumes lift per unit span proportional to local chord
                       TCS([Mr >= Lmax*AR*p/24]),

                       # Root stiffness (see Hoburg 2014)
                       # Assumes rh = 0.75, so that rms box height = ~0.92*tmax
                       TCS([0.92*wwb*tau*tcap**2 + Icap <= 0.92**2/2*wwb*tau**2*tcap]),

                       # Stress limit
                       # Assumes bending stress carried by caps (Icap >> Iweb)
                       TCS([8 >= Mr*AR*q**2*tau/(S*Icap*sigmax)]),

                       # Shear web sizing
                       # Assumes all shear loads are carried by web and rh=0.75
                       TCS([12 >= AR*Lmax*q**2/(tau*S*tweb*sigmaxshear)]),

                       # Posynomial approximation of nu=(1+lam+lam^2)/(1+lam^2)
                       nu**3.94 >= 0.86*p**(-2.38)+ 0.14*p**0.56,

                       # Weight of spar caps and shear webs
                       Wcap >= 8*rhocap*g*wwb*tcap*S**1.5*nu/(3*AR**0.5),
                       TCS([Wweb >= 8*rhoweb*g*rh*tau*tweb*S**1.5*nu/(3*AR**0.5)]),

                       # Total wing weight using an additional weight fraction
                       Wstruct >= (1 + fwadd)*(Wweb + Wcap),
                       #based off of a raymer weight and 737 data from TASOPT output file
##                        (S/(dum1))**.65 * (AR/10.1)**.5 == Wstruct/(dum2),
                       ]

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
