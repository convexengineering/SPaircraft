"""Simple commercial aircraft flight profile and aircraft model"""
from numpy import pi, cos, tan
import numpy as np
from gpkit import Variable, Model, units, SignomialsEnabled, Vectorize
from gpkit.constraints.sigeq import SignomialEquality
from gpkit.tools import te_exp_minus1
from gpkit.constraints.tight import Tight as TCS
import matplotlib.pyplot as plt
from simple_ac_imports_no_engine import Fuselage, Engine, CruiseP, ClimbP, FlightState, CruiseSegment, ClimbSegment

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

D80 = True

class Aircraft(Model):
    "Aircraft class"
    def setup(self, **kwargs):
        #create submodels
        self.fuse = Fuselage()
        self.wing = Wing()
        self.engine = Engine()

        #variable definitions
        numeng = Variable('numeng', '-', 'Number of Engines')

        constraints = [self.wing['x_w'] == self.fuse['l_{fuse}']*0.6,
                       ]

        constraints.extend([
            numeng == numeng, #need numeng in the model
            ])

        self.components = [self.fuse, self.wing, self.engine]

        return self.components, constraints

    def dynamic(self, state):
        """
        creates an aircraft climb performance model, given a state
        """
        return AircraftP(self, state)
        
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
            W_burn == aircraft['numeng']*self.engineP['TSFC'] * thours * self.engineP['F'],
               
            #time unit conversion
            t == thours,
            ])

        return self.Pmodels, constraints
    
class Mission(Model):
    """
    mission class, links together all subclasses
    """
    def setup(self, subs = None, **kwargs):
        #define the number of each flight segment
        Nclimb = 2
        Ncruise = 2

        #build required submodels
        aircraft = Aircraft()

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

        h = climb.state['h']
        hftClimb = climb.state['hft']
        dhft = climb.climbP['dhft']
        hftCruise = cruise.state['hft']

        #make overall constraints
        constraints = []

        constraints.extend([
            # Wing max loading constraint
            aircraft['L_{max}'] >= aircraft['N_{lift}'] * W_total,

            #weight constraints
            TCS([aircraft['W_{e}'] + aircraft['W_{payload}'] + W_ftotal + aircraft['numeng'] * aircraft['W_{engine}'] + aircraft.wing.wb['W_{struct}'] <= W_total]),

            climb.climbP.aircraftP['W_{start}'][0] == W_total,
            climb.climbP.aircraftP['W_{end}'][-1] == cruise.cruiseP.aircraftP['W_{start}'][0],

            # similar constraint 1
            TCS([climb.climbP.aircraftP['W_{start}'] >= climb.climbP.aircraftP['W_{end}'] + climb.climbP.aircraftP['W_{burn}']]),
            # similar constraint 2
            TCS([cruise.cruiseP.aircraftP['W_{start}'] >= cruise.cruiseP.aircraftP['W_{end}'] + cruise.cruiseP.aircraftP['W_{burn}']]),

            climb.climbP.aircraftP['W_{start}'][1:] == climb.climbP.aircraftP['W_{end}'][:-1],
            cruise.cruiseP.aircraftP['W_{start}'][1:] == cruise.cruiseP.aircraftP['W_{end}'][:-1],

            TCS([aircraft['W_{e}'] + aircraft['W_{payload}'] + aircraft['numeng'] * aircraft['W_{engine}'] + aircraft.wing.wb['W_{struct}'] <= cruise.cruiseP.aircraftP['W_{end}'][-1]]),

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
            climb.climbP.engineP['F'] <= 2 * max(cruise.cruiseP.engineP['F']),

            #set the range for each cruise segment, doesn't take credit for climb
            #down range disatnce covered
            cruise.cruiseP['Rng'] == ReqRng/(Ncruise),

            #set the TSFC
            climb.climbP.engineP['TSFC'] == .7*units('1/hr'),
            cruise.cruiseP.engineP['TSFC'] == .5*units('1/hr'),

            #wing constraints
            aircraft.wing['W_{fuel_{wing}}'] == W_ftotal,
            climb.climbP.wingP['L_w'] == climb.climbP.aircraftP['W_{avg}'],
            cruise.cruiseP.wingP['L_w'] == cruise.cruiseP.aircraftP['W_{avg}'],
            climb['c_{m_{w}}'] == .10, # for boundedness
            cruise['c_{m_{w}}'] == .10, # for boundedness
            ])
        
        return aircraft, climb, cruise, constraints

class Wing(Model):
    """
    Philippe's thesis wing model
    """
    def setup(self, **kwargs):
        self.wns = WingNoStruct()
        self.wb = WingBox(self.wns)

        Wwing = Variable('W_{wing_system}', 'N', 'Total Wing Weight')

        dxACwing = Variable('\\Delta x_{AC_{wing}}','m','Wing Aerodynamic Center Shift')
            # w.r.t. the quarter chord of the root of the wing.

        constraints = []
        with SignomialsEnabled():
            constraints.extend([
            self.wns['\\lambda'] == self.wb['taper'],

            TCS([Wwing >= self.wb['W_{struct}'] + self.wb['W_{struct}']*(self.wns['f_{flap}'] + \
                    self.wns['f_{slat}'] + self.wns['f_{aileron}'] + self.wns['f_{lete}'] + self.wns['f_{ribs}'] + \
                    self.wns['f_{spoiler}'] + self.wns['f_{watt}'])]),
            SignomialEquality(dxACwing,(1/12.*self.wns['A_{tri}'] + 1/4.*self.wns['A_{rect}'])/self.wns['S'] \
                 *self.wns['b']*self.wns['\\tan(\\Lambda)']),

##            self.wns['\\bar{A}_{fuel, max}'] <= (self.wns['wwn'] - 2*self.wb['t_{web}'])*(0.92*self.wns['\\tau'] - 2*self.wb['t_{web}']),
            ])

        return self.wns, self.wb, constraints
        
    def dynamic(self, state):
        """
        returns an instance of the wing perofrmance model
        """
        return WingPerformance(self, state)

class WingNoStruct(Model):
    """
    Philippe's wing model minus structure
    """
    def setup(self, **kwargs):
        #declare variables
               #Variables
        Afuel   = Variable('\\bar{A}_{fuel, max}', '-', 'Non-dim. fuel area')
        
        CLwmax  = Variable('C_{L_{wmax}}', '-', 'Max lift coefficient, wing')
        
        
        Vfuel   = Variable('V_{fuel, max}', 'm^3', 'Available fuel volume')
        
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
        wwn       = Variable('wwn', 0.5, '-', 'Wingbox-width-to-chord ratio')
        xw     = Variable('x_w', 'm', 'Position of wing aerodynamic center')
        ymac    = Variable('y_{mac}', 'm',
                           'Spanwise location of mean aerodynamic chord')
        bmax = Variable('b_{max}', 'm', 'Max Wing Span')

        #Linked Variables
        AR      = Variable('AR', '-', 'Wing aspect ratio')
        Lmax    = Variable('L_{max}', 'N', 'Maximum load')
        Sw      = Variable('S', 'm^2', 'Wing area')
        Vne     = Variable('V_{ne}', 'm/s', 'Never exceed velocity')
        WfuelWing   = Variable('W_{fuel_{wing}}', 'N', 'Fuel weight')
        b       = Variable('b', 'm', 'Wing span')
        #the following two variables have the same name in the flight profile and
        #will be automatically linked by the linked constraint set
        mac    = Variable('mac', 'm',
                          'Mean aerodynamic chord (wing)')
        e       = Variable('e', '-', 'Oswald efficiency factor')
        FuelFrac = Variable('FuelFrac', '-', 'Usability Factor of Max Fuel Volume')

        #fractional componenet weights
        fflap = Variable('f_{flap}', '-', 'Flap Fractional Weight')
        fslat = Variable('f_{slat}', '-', 'Slat Fractional Weight')
        faile = Variable('f_{aileron}', '-', 'Aileron Fractional Weight')
        flete = Variable('f_{lete}', '-', 'Lete Fractional Weight')
        fribs = Variable('f_{ribs}', '-', 'Rib Fractional Weight')
        fspoi = Variable('f_{spoiler}', '-', 'Spoiler Fractional Weight')
        fwatt = Variable('f_{watt}', '-', 'Watt Fractional Weight')

        # Area fractions
        Atri = Variable('A_{tri}','m^2','Triangular Wing Area')
        Arect = Variable('A_{rect}','m^2','Rectangular Wing Area')
        
        #make constraints
        constraints = []

        with SignomialsEnabled():

            constraints.extend([
                 SignomialEquality(Arect, ctip*b),
                 SignomialEquality(Atri, 1./2.*(1-taper)*croot*b), #[SP]
                 p >= 1 + 2*taper,
                 2*q >= 1 + p,
                 ymac == (b/3)*q/p,
                 TCS([(2./3)*(1+taper+taper**2)*croot/q <= mac],
                                   reltol=1E-2),

                 taper == ctip/croot,
                 # ymac  <= b/2,

                 SignomialEquality(Sw, b*(croot + ctip)/2),

                 # Oswald efficiency
                 # Nita, Scholz, "Estimating the Oswald factor from
                 # basic aircraft geometrical parameters"
                 TCS([fl >= 0.0524*taper**4 - 0.15*taper**3
                         + 0.1659*taper**2 - 0.0706*taper + 0.0119],
                    reltol=1E-2),
                TCS([e*(1 + fl*AR) <= 1]),
                taper >= 0.15, # TODO

                # Fuel volume [TASOPT doc]
                # TCS([Afuel <= wwn*0.92*tau]),
                # # GP approx of the signomial constraint:
                # Vfuel <= mac**2*b*tau,
                # WfuelWing <= rhofuel*Vfuel*g,
                #
                # b <= bmax,

                # Fuel volume [TASOPT doc]
                TCS([Afuel <= wwn*0.92*tau]),
                # GP approx of the signomial constraint:
                # Afuel <= (w - 2*tweb)*(0.92*tau - 2*tcap),
                Vfuel <= croot**2 * (b/6) * (1+taper+taper**2)*cosL,
                WfuelWing <= rhofuel*Afuel*Vfuel*g,

                b <= bmax,
                ])

        return constraints

class WingPerformance(Model):
    """
    Wing performance model
    """
    def setup(self, wing, state, **kwargs):
        self.wing = wing
        
        #declare variables
        #Vector Variables
        alpha   = Variable('\\alpha_w', '-', 'Wing angle of attack')
        CLaw    = Variable('C_{L_{aw}}', '-', 'Lift curve slope, wing')
        
        Re      = Variable('Re_w', '-', 'Reynolds number (wing)')
        CDp     = Variable('C_{D_{p_w}}', '-',
                           'Wing parasitic drag coefficient')
        CDw     = Variable('C_{d_w}', '-', 'Drag coefficient, wing')
        CLw     = Variable('C_{L}', '-', 'Lift coefficient, wing')
 
        D       = Variable('D_{wing}', 'N', 'Wing drag')
        Lw      = Variable('L_w', 'N', 'Wing lift')

        # Center wing section lift reduction variables
        dLo     = Variable('\\Delta L_{o}','N','Center wing lift loss')
        etao    = Variable('\\eta_{o}','-','Center wing span coefficient')
        po      = Variable('p_{o}','N/m','Center section theoretical wing loading')
        fLo     = Variable('f_{L_{o}}',0.5,'-','Center wing lift reduction coefficient')

        # Wing tip lift reduction variables
        dLt = Variable('\\Delta L_{t}','N','Wing tip lift loss')
        fLt = Variable('f_{L_{t}}',0.05,'-','Wing tip lift reduction coefficient')

        #wing moment variables -- need a good way to model this, currently using TAT
        cmw = Variable('c_{m_{w}}', '-', 'Wing Pitching Moment Coefficient')

        amax    = Variable('\\alpha_{max,w}', '-', 'Max angle of attack')
        
        #make constraints
        constraints = []

        with SignomialsEnabled():
            constraints.extend([
                # Lw == 0.5*state['\\rho']*state['V']**2*self.wing['S']*CLw,
                0.5*state['\\rho']*state['V']**2*self.wing['S']*CLw >= Lw + dLo + 2.*dLt,
                dLo == etao*fLo*self.wing['b']/2*po,
                dLt == fLt*po*self.wing['c_{root}']*self.wing['taper']**2, # TODO improve approximations croot~co and taper~gammat

                # DATCOM formula (Mach number makes it SP)
                # Swept wing lift curve slope constraint
                SignomialEquality((self.wing['AR']/self.wing['\\eta'])**2*(1 + self.wing['\\tan(\\Lambda)']**2 - state['M']**2) + 8*pi*self.wing['AR']/CLaw
                      , (2*pi*self.wing['AR']/CLaw)**2),
                
                CLw == CLaw*alpha,
                alpha <= amax,

                # Drag
                D == 0.5*state['\\rho']*state['V']**2*self.wing['S']*CDw,
                TCS([CDw >= CDp + CLw**2/(pi*self.wing['e']*self.wing['AR'])]),
                Re == state['\\rho']*state['V']*self.wing['mac']/state['\\mu'],

                TCS([CDp**6.5 >= (1.02458748e10 * CLw**15.587947404823325 * (self.wing['\\cos(\\Lambda)']*state['M'])**156.86410659495155 +
                                 2.85612227e-13 * CLw**1.2774976672501526 * (self.wing['\\cos(\\Lambda)']*state['M'])**6.2534328002723703 +
                                 2.08095341e-14 * CLw**0.8825277088649582 * (self.wing['\\cos(\\Lambda)']*state['M'])**0.0273667615730107 +
                                 1.94411925e+06 * CLw**5.6547413360261691 * (self.wing['\\cos(\\Lambda)']*state['M'])**146.51920742858428)]),
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
        Icap    = Variable('I_{capprint s}', '-',
                           'Non-dim spar cap area moment of inertia')
        Mr      = Variable('M_r', 'N', 'Root moment per root chord')
        nu      = Variable('\\nu', '-',
                           'Dummy variable = $(t^2 + t + 1)/(t+1)$')
        Wcap    = Variable('W_{cap}', 'N', 'Weight of spar caps')
        Wweb    = Variable('W_{web}', 'N', 'Weight of shear web')
        Wstruct = Variable('W_{struct}', 'N', 'Structural weight')

        # Constants
        taper = Variable('taper', '-', 'Taper ratio')
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

                       # Posynomial approximation of nu=(1+lam+lam^2)/(1+lam)
                       # nu**3.94 >= 0.86*p**(-2.38)+ 0.14*p**0.56, # PHILIPPE'S FIT
                        (nu/1.09074074)**.166 >= 0.205*(p/1.7)**0.772 + 0.795*(p/1.7)**-0.125, # BERK'S FIT


                       # Weight of spar caps and shear webs
                       Wcap >= 8*rhocap*g*wwb*tcap*S**1.5*nu/(3*AR**0.5),
                       TCS([Wweb >= 8*rhoweb*g*rh*tau*tweb*S**1.5*nu/(3*AR**0.5)]),

                       # Total wing weight using an additional weight fraction
                       Wstruct >= (Wweb + Wcap),
                       ]
        
        return constraints

if __name__ == '__main__':
    plot = False

    sweep = 30 #[deg]
    
    substitutions = {      
            'ReqRng': 500, #('sweep', np.linspace(500,2000,4)),
            'CruiseAlt': 30000, #('sweep', np.linspace(20000,40000,4)),
            'numeng': 1,
            'W_{pax}': 91 * 9.81,
            'n_{pax}': 150,
            'pax_{area}': 1,

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

    if D80:
        sweep = 27.566
        substitutions.update({
            'ReqRng': 3000,
            'CruiseAlt':36632,
            'numeng':2,
            'W_{pax}': 91*9.81,
            'n_{pax}':180,
            'pax_{area}': 1,

            #wing subs
            'C_{L_{wmax}}': 2.5,
            'V_{ne}': 144,
            '\\alpha_{max,w}': 0.1, # (6 deg)
            '\\cos(\\Lambda)': cos(sweep*pi/180),
            '\\eta': 0.97,
            '\\rho_0': 1.225,
            '\\rho_{fuel}': 817, # Kerosene [TASOPT]
            '\\tan(\\Lambda)': tan(sweep*pi/180),

            'W_{engine}': 15100.3*0.454*9.81, #units('N')
            'AR':10.8730,
            'b':116.548*0.3048,#units('ft'),
            'c_{root}': 17.4*0.3048,#units('ft'),

                # Minimum Cruise Mach Number
                # 'M': 0.8,
        })
    mission = Mission()
    m = Model(mission['W_{f_{total}}'], mission, substitutions)
    sol = m.localsolve(solver='mosek', verbosity = 2)

    if plot == True:
         substitutions = {      
                 'ReqRng': 2000,#('sweep', np.linspace(500,3000,8)),
                 'CruiseAlt': 30000,
                 'numeng': 1,
                 'W_{pax}': 91 * 9.81,
                 'n_{pax}': 150,
                 'pax_{area}': 1,

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
               
         mission = Mission()
         m = Model(mission['W_{f_{total}}'], mission, substitutions)
         solRsweep = m.localsolve(solver='mosek', verbosity = 4)




#          plt.plot(solRsweep('ReqRng'), solRsweep('W_{struct}'), '-r')
#          plt.xlabel('Mission Range [nm]')
#          plt.ylabel('Wing Weight [N]')
#          plt.title('Wing Weight vs Range')
# ##         plt.savefig('Wing_Sweeps/wing_rng_Wstruct.pdf')
#          plt.show()

#          plt.plot(solRsweep('ReqRng'), solRsweep('AR'), '-r')
#          plt.xlabel('Mission Range [nm]')
#          plt.ylabel('Wing Aspect Ratio')
#          plt.title('Wing Aspect Ratio vs Range')
# ##         plt.savefig('Wing_Sweeps/wing_rng_AR.pdf')
#          plt.show()

#          plt.plot(solRsweep('ReqRng'), solRsweep('S'), '-r')
#          plt.xlabel('Mission Range [nm]')
#          plt.ylabel('VWing Area [m$^2$]')
#          plt.title('Wing Area vs Range')
# ##         plt.savefig('Wing_Sweeps/wing_rng_S.pdf')
#          plt.show()

#          plt.plot(solRsweep('ReqRng'), solRsweep('b'), '-r')
#          plt.xlabel('Mission Range [nm]')
#          plt.ylabel('VWing Span [m]')
#          plt.title('Wing Span vs Range')
# ##         plt.savefig('Wing_Sweeps/wing_rng_b.pdf')
#          plt.show()

#          plt.plot(solRsweep('ReqRng'), solRsweep['sensitivities']['constants']['V_{ne}'], '-r')
#          plt.xlabel('Mission Range [nm]')
#          plt.ylabel('Sensitivity of $V_{ne}$')
#          plt.title('Sensitivity of $V_{ne}$ vs Range')
# ##         plt.savefig('Wing_Sweeps/wing_rng_SensVne.pdf')
#          plt.show()

#          plt.plot(solRsweep('ReqRng'), solRsweep['sensitivities']['constants']['N_{lift}'], '-r')
#          plt.xlabel('Mission Range [nm]')
#          plt.ylabel('Sensitivity to Seciontal Lift Multiplier')
#          plt.title('Sensitivity to Seciontal Lift Multiplier vs Range')
# ##         plt.savefig('Wing_Sweeps/wing_rng_Nlift.pdf')
#          plt.show()

#          plt.plot(solRsweep('ReqRng'), solRsweep['sensitivities']['constants']['C_{L_{wmax}}'], '-r')
#          plt.xlabel('Mission Range [nm]')
#          plt.ylabel('Sensitivity to Wing $C_{L_{max}}$')
#          plt.title('Sensitivity to Wing $C_{L_{max}}$ vs Range')
# ##         plt.savefig('Wing_Sweeps/wing_rng_SensClMax.pdf')
#          plt.show()

#          plt.plot(solRsweep('ReqRng'), solRsweep['sensitivities']['constants']['ReqRng'], '-r')
#          plt.xlabel('Mission Range [nm]')
#          plt.ylabel('Sensitivity to Required Range')
#          plt.title('Sensitivity to Required Range vs Range')
# ##         plt.savefig('Wing_Sweeps/wing_rng_SensRng.pdf')
#          plt.show()

#          plt.plot(solRsweep('ReqRng'), solRsweep['sensitivities']['constants']['CruiseAlt'], '-r')
#          plt.xlabel('Mission Range [nm]')
#          plt.ylabel('Sensitivity to Cruise Altitude')
#          plt.title('Sensitivity to Cruise Altitude vs Range')
# ##         plt.savefig('Wing_Sweeps/wing_rng_SensAlt.pdf')
#          plt.show()

#          substitutions = {      
#      ##            'V_{stall}': 120,
#                  'ReqRng': 500,
#                  'CruiseAlt': ('sweep', np.linspace(20000,40000,8)),
#                  'numeng': 1,
#      ##            'W_{Load_max}': 6664,
#                  'W_{pax}': 91 * 9.81,
#                  'n_{pax}': 150,
#                  'pax_{area}': 1,
#      ##            'C_{D_{fuse}}': .005, #assumes flat plate turbulent flow, from wikipedia

#                  #wing subs
#                  'C_{L_{wmax}}': 2.5,
#                  'V_{ne}': 144,
#                  '\\alpha_{max,w}': 0.1, # (6 deg)
#                  '\\cos(\\Lambda)': cos(sweep*pi/180),
#                  '\\eta': 0.97,
#                  '\\rho_0': 1.225,
#                  '\\rho_{fuel}': 817, # Kerosene [TASOPT]
#                  '\\tan(\\Lambda)': tan(sweep*pi/180),
#                  }
               
#          mission = Mission()
#          m = Model(mission['W_{f_{total}}'], mission, substitutions)
#          solAltsweep = m.localsolve(solver='mosek', verbosity = 4)

#          plt.plot(solAltsweep('CruiseAlt'), solAltsweep('W_{struct}'), '-r')
#          plt.xlabel('Cruise Altitude [ft]')
#          plt.ylabel('Wing Weight [N]')
#          plt.title('Wing Weight vs Cruise Altitude')
# ##         plt.savefig('Wing_Sweeps/wing_alt_Wstruct.pdf')
#          plt.show()

#          plt.plot(solAltsweep('CruiseAlt'), solAltsweep('AR'), '-r')
#          plt.xlabel('Cruise Altitude [ft]')
#          plt.ylabel('Wing Aspect Ratio')
#          plt.title('Wing Aspect Ratio vs Cruise Altitude')
# ##         plt.savefig('Wing_Sweeps/wing_alt_AR.pdf')
#          plt.show()

#          plt.plot(solAltsweep('CruiseAlt'), solAltsweep('S'), '-r')
#          plt.xlabel('Cruise Altitude [ft]')
#          plt.ylabel('Wing Area [m$^2$]')
#          plt.title('Wing Area vs Cruise Altitude')
# ##         plt.savefig('Wing_Sweeps/wing_alt_S.pdf')
#          plt.show()

#          plt.plot(solAltsweep('CruiseAlt'), solAltsweep('b'), '-r')
#          plt.xlabel('Cruise Altitude [ft]')
#          plt.ylabel('Wing Span [m]')
#          plt.title('Wing Span vs Cruise Altitude')
# ##         plt.savefig('Wing_Sweeps/wing_alt_b.pdf')
#          plt.show()

#          plt.plot(solAltsweep('CruiseAlt'), solAltsweep['sensitivities']['constants']['V_{ne}'], '-r')
#          plt.xlabel('Cruise Altitude [ft]')
#          plt.ylabel('Sensitivity of $V_{ne}$')
#          plt.title('Sensitivity of $V_{ne}$ vs Cruise Altitude')
# ##         plt.savefig('Wing_Sweeps/wing_alt_SensNve.pdf')
#          plt.show()

#          plt.plot(solAltsweep('CruiseAlt'), solAltsweep['sensitivities']['constants']['N_{lift}'], '-r')
#          plt.xlabel('Cruise Altitude [ft]')
#          plt.ylabel('Sensitivity to Seciontal Lift Multiplier')
#          plt.title('Sensitivity to Seciontal Lift Multiplier vs Cruise Altitude')
# ##         plt.savefig('Wing_Sweeps/wing_alt_SensNLift.pdf')
#          plt.show()

#          plt.plot(solAltsweep('CruiseAlt'), solAltsweep['sensitivities']['constants']['C_{L_{wmax}}'], '-r')
#          plt.xlabel('Cruise Altitude [ft]')
#          plt.ylabel('Sensitivity to Wing $C_{L_{max}}$')
#          plt.title('Sensitivity to Wing $C_{L_{max}}$ vs Cruise Altitude')
# ##         plt.savefig('Wing_Sweeps/wing_alt_SensClMax.pdf')
#          plt.show()

#          plt.plot(solAltsweep('CruiseAlt'), solAltsweep['sensitivities']['constants']['ReqRng'], '-r')
#          plt.xlabel('Cruise Altitude [ft]')
#          plt.ylabel('Sensitivity to Required Range')
#          plt.title('Sensitivity to Required Range vs Cruise Altitude')
# ##         plt.savefig('Wing_Sweeps/wing_alt_SensRng.pdf')
#          plt.show()

#          plt.plot(solAltsweep('CruiseAlt'), solAltsweep['sensitivities']['constants']['CruiseAlt'], '-r')
#          plt.xlabel('Cruise Altitude [ft]')
#          plt.ylabel('Sensitivity to Cruise Altitude')
#          plt.title('Sensitivity to Cruise Altitude vs Cruise Altitude')
# ##         plt.savefig('Wing_Sweeps/wing_alt_SensAlt.pdf')
#          plt.show()

