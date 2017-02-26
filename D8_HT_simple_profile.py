"""D8 horizontal tail model linked with a simple flight profile"""
from numpy import pi, tan
import numpy as np
from gpkit import Variable, Model, units, SignomialsEnabled, Vectorize, SignomialEquality
from gpkit.tools import te_exp_minus1
from gpkit.constraints.tight import Tight as TCS
import matplotlib.pyplot as plt
from gpkit.small_scripts import mag
from simple_ac_imports_no_engine import Wing, Fuselage, Engine, CruiseP, ClimbP, FlightState, CruiseSegment, ClimbSegment

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
                W_burn == aircraft['numeng']*self.engineP['TSFC'] * thours * self.engineP['F'],
                   
                #time unit conversion
                t == thours,

                #make lift equal weight --> small angle approx in climb
                self.wingP['L_{wing}'] >= W_avg,
                 ])

        return self.Pmodels, constraints

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
        W_dry = Variable('W_{dry}', 'N', 'Aircraft Dry Weight')

        h = climb['h']
        hftClimb = climb['hft']
        dhft = climb['dhft']
        hftCruise = cruise['hft']

        #make overall constraints
        constraints = []

        constraints.extend([
            #weight constraints
            TCS([aircraft['W_{e}'] + aircraft['W_{payload}'] + aircraft['numeng'] * aircraft['W_{engine}'] + aircraft['W_{wing}'] + aircraft.HT['W_{struct}'] <= W_dry]),
            TCS([W_ftotal + W_dry <= W_total]),

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

            #constrain the thrust
            climb.climbP['F'] <= 2 * max(cruise.cruiseP['F']),

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
        CLfCG = Variable('C_{L_{hfcG}}', '-', 'HT CL During Max Forward CG')

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
                TCS([(2./3)*(1 + taper + taper**2)*croot/q <=
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

        #HT system weight variable
        WHT = Variable('W_{HT_system}', 'N', 'HT System Weight')
        fHT = Variable('f_{HT}' ,'-', 'Rudder etc. fractional weight')

        constraints = []
        with SignomialsEnabled():
            constraints.append([
                self.wb['L_{h_{rect}}'] >= self.wb['L_{{max}_h}']/2.*self.HTns['c_{tip_h}']*self.HTns['b_{ht}']/self.HTns['S_h'],
                self.wb['L_{h_{tri}}'] >= self.wb['L_{{max}_h}']/4.*(1-self.wb['taper'])*self.HTns['c_{root_h}']*self.HTns['b_{ht}']/self.HTns['S_h'], #[SP]

                WHT >= self.wb['W_{struct}'] + self.wb['W_{struct}']  * fHT,
            ])

        return self.HTns, self.wb, constraints


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
        taper = Variable('taper', 0.3, '-', 'Taper ratio')
        fwadd  = Variable('f_{w,add}', 0.3, '-',
                          'Wing added weight fraction') # [TAS]
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


        # Pi tail sizing variables
        # Splits the max lift into triangular and rectangular components
        # for root bending sizing.
        Lhtri = Variable('L_{h_{tri}}','N','Triangular HT load')
        Lhrect = Variable('L_{h_{rect}}','N','Rectangular HT load')
        
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
                       tau <= 0.14,

                       # Root moment calculation (see Hoburg 2014)
                       # Depends on a given load the wing must support, Lmax
                       # Assumes lift per unit span proportional to local chord
                       # Mr >= Lmax*AR*p/24,

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
            'numeng': 2,
            'W_{pax}': 91 * 9.81,
            'n_{pax}': 150,
            'pax_{area}': 1,
            'e': .9,
            'b_{max}': 60,

            #HT subs
            'C_{L_{hmax}}': 2.5,
            '\\tan(\\Lambda_{ht})': tan(30*pi/180),
            'w_{fuse}': 6,
            'c_{m_{w}}': 1,
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
    sol = m.localsolve(solver='mosek', verbosity = 4)

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
