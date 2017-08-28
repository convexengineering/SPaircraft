"""Simple commercial aircraft flight profile and aircraft model"""
from numpy import pi, cos, tan
import numpy as np
from gpkit import Variable, Model, units, SignomialsEnabled, Vectorize
from gpkit.constraints.sigeq import SignomialEquality
from gpkit.tools import te_exp_minus1
from gpkit.constraints.tight import Tight as TCS
from wingbox import WingBox

class Wing(Model):
    """
    Philippe's thesis wing model
    """
    def setup(self, **kwargs):
        self.wns = WingNoStruct()
        self.wb = WingBox(self.wns, "wing")

        Wwing = Variable('W_{wing}', 'N', 'Wing System Weight')

        Cwing = Variable('C_{wing}', 1, '-', 'Wing Weight Margin and Sensitivity Factor')

        # w.r.t. the quarter chord of the root of the wing.
        dxACwing = Variable('\\Delta x_{AC_{wing}}','m','Wing Aerodynamic Center Shift')

        #wing induced drag reduction due to wing tip devices
        TipReduct = Variable('TipReduct', '-', 'Induced Drag Reduction Factor from Wing Tip Devices')

        constraints = []
        with SignomialsEnabled():
            constraints.extend([
            self.wns['\\lambda'] == self.wb['taper'],

            TCS([Wwing >= Cwing * self.wb['W_{struct}'] + self.wb['W_{struct}']*(self.wns['f_{flap}'] + \
                    self.wns['f_{slat}'] + self.wns['f_{aileron}'] + self.wns['f_{lete}'] + self.wns['f_{ribs}'] + \
                    self.wns['f_{spoiler}'] + self.wns['f_{watt}'])]),
            TCS([dxACwing <= 1./24.*(self.wns['c_{root}'] + 5.*self.wns['c_{tip}'])/self.wns['S'] \
                            *self.wns['b']**2*self.wns['\\tan(\\Lambda)']]),
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
        
        CLwmax  = Variable('C_{L_{w,max}}', '-', 'Max lift coefficient, wing')
        
        
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
        tau_max = Variable('\\tau_{max_w}', '-', 'Max allowed wing thickness')
        wwn       = Variable('wwn', 0.5, '-', 'Wingbox-width-to-chord ratio')
        xw     = Variable('x_w', 'm', 'Position of wing aerodynamic center')
        ymac    = Variable('y_{mac}', 'm',
                           'Spanwise location of mean aerodynamic chord')
        bmax = Variable('b_{max}', 'm', 'Max Wing Span')

        #Linked Variables
        AR      = Variable('AR', '-', 'Wing aspect ratio')
        Lmax    = Variable('L_{max}', 'N', 'Maximum load')
        Sw      = Variable('S', 'm^2', 'Wing area')
        WfuelWing   = Variable('W_{fuel_{wing}}', 'N', 'Fuel weight')
        b       = Variable('b', 'm', 'Wing span')
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
                 Arect == ctip*b,
                 Atri >= 1./2.*(1-taper)*croot*b, #[SP]
                 p >= 1 + 2*taper,
                 2*q >= 1 + p,
                 ymac == (b/3)*q/p,
                 TCS([(2./3)*(1+taper+taper**2)*croot/q <= mac],
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

                taper >= 0.15, # TODO

                # Fuel volume [TASOPT doc]
                # GP approx of the signomial constraint:
                Vfuel <= .3026*mac**2*b*tau,
                 
                WfuelWing <= rhofuel*Vfuel*g,

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
        CLaw    = Variable('C_{L_{\\alpha,w}}', '-', 'Lift curve slope, wing')
        
        Re      = Variable('Re_w', '-', 'Reynolds number (wing)')
        CDp     = Variable('C_{D_{p_w}}', '-',
                           'Wing parasitic drag coefficient')
        CDi     = Variable('C_{D_{i_w}}', '-',
                           'Wing induced drag coefficient')
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
                TCS([CDw >= CDp + CDi]),
                TCS([CDi >= self.wing['TipReduct']*CLw**2/(pi*(self.wing['e'])*self.wing['AR'])]),
                Re == state['\\rho']*state['V']*self.wing['mac']/state['\\mu'],

                #original Philippe thesis fit
##                TCS([CDp**6.5 >= (1.02458748e10 * CLw**15.587947404823325 * (self.wing['\\cos(\\Lambda)']*state['M'])**156.86410659495155 +
##                                 2.85612227e-13 * CLw**1.2774976672501526 * (self.wing['\\cos(\\Lambda)']*state['M'])**6.2534328002723703 +
##                                 2.08095341e-14 * CLw**0.8825277088649582 * (self.wing['\\cos(\\Lambda)']*state['M'])**0.0273667615730107 +
##                                 1.94411925e+06 * CLw**5.6547413360261691 * (self.wing['\\cos(\\Lambda)']*state['M'])**146.51920742858428)]),
##                ])

            #Martin's TASOPT c series airfoil fit
            TCS([
                CDp**1.6515 >= 1.61418 * (Re/1000)**-0.550434 * (self.wing['\\tau'])**1.29151 * (self.wing['\\cos(\\Lambda)']*state['M'])**3.03609 * CLw**1.77743
                    + 0.0466407 * (Re/1000)**-0.389048 * (self.wing['\\tau'])**0.784123 * (self.wing['\\cos(\\Lambda)']*state['M'])**-0.340157 * CLw**0.950763
                    + 190.811 * (Re/1000)**-0.218621 * (self.wing['\\tau'])**3.94654 * (self.wing['\\cos(\\Lambda)']*state['M'])**19.2524 * CLw**1.15233
                    + 2.82283e-12 * (Re/1000)**1.18147 * (self.wing['\\tau'])**-1.75664 * (self.wing['\\cos(\\Lambda)']*state['M'])**0.10563 * CLw**-1.44114
                ]),
            ])

        return constraints

## Old code, not currenlty working. Can be used for subsystem testing
## in the future

##if __name__ == '__main__':
##    plot = False
##
##    sweep = 30 #[deg]
##
##    substitutions = {      
##            'R_{req}': 500, #('sweep', np.linspace(500,2000,4)),
##            'CruiseAlt': 30000, #('sweep', np.linspace(20000,40000,4)),
##            'numeng': 1,
##            'W_{pass}': 91 * 9.81,
##            'n_{pass}': 150,
##            'pax_{area}': 1,
##
##            #wing subs
##            'C_{L_{w,max}}': 2.5,
##            'V_{ne}': 144,
##            '\\alpha_{max,w}': 0.1, # (6 deg)
##            '\\cos(\\Lambda)': cos(sweep*pi/180),
##            '\\eta': 0.97,
##            '\\rho_0': 1.225,
##            '\\rho_{fuel}': 817, # Kerosene [TASOPT]
##            '\\tan(\\Lambda)': tan(sweep*pi/180),
##            }
##
##    if D80:
##        sweep = 27.566
##        substitutions.update({
##            'R_{req}': 3000,
##            'CruiseAlt':36632,
##            'numeng':2,
##            'W_{pass}': 91*9.81,
##            'n_{pass}':180,
##            'pax_{area}': 1,
##
##            #wing subs
##            'C_{L_{w,max}}': 2.5,
##            'V_{ne}': 144,
##            '\\alpha_{max,w}': 0.1, # (6 deg)
##            '\\cos(\\Lambda)': cos(sweep*pi/180),
##            '\\eta': 0.97,
##            '\\rho_0': 1.225,
##            '\\rho_{fuel}': 817, # Kerosene [TASOPT]
##            '\\tan(\\Lambda)': tan(sweep*pi/180),
##
##            'W_{engine}': 15100.3*0.454*9.81, #units('N')
##            'AR':10.8730,
##            'b':116.548*0.3048,#units('ft'),
##            'c_{root}': 17.4*0.3048,#units('ft'),
##
##                # Minimum Cruise Mach Number
##                # 'M': 0.8,
##        })
##    mission = Mission()
##    m = Model(mission['W_{f_{total}}'], mission, substitutions)
##    sol = m.localsolve(solver='mosek', verbosity = 2)
##
##    if plot == True:
##         substitutions = {      
##                 'R_{req}': 2000,#('sweep', np.linspace(500,3000,8)),
##                 'CruiseAlt': 30000,
##                 'numeng': 1,
##                 'W_{pass}': 91 * 9.81,
##                 'n_{pass}': 150,
##                 'pax_{area}': 1,
##
##                 #wing subs
##                 'C_{L_{w,max}}': 2.5,
##                 'V_{ne}': 144,
##                 '\\alpha_{max,w}': 0.1, # (6 deg)
##                 '\\cos(\\Lambda)': cos(sweep*pi/180),
##                 '\\eta': 0.97,
##                 '\\rho_0': 1.225,
##                 '\\rho_{fuel}': 817, # Kerosene [TASOPT]
##                 '\\tan(\\Lambda)': tan(sweep*pi/180),
##                 }
##               
##         mission = Mission()
##         m = Model(mission['W_{f_{total}}'], mission, substitutions)
##         solRsweep = m.localsolve(solver='mosek', verbosity = 4)
##
##
##
##from simple_ac_imports_no_engine import Fuselage, Engine, CruiseP, ClimbP, FlightState, CruiseSegment, ClimbSegment
##
##D80 = True
##
##class Aircraft(Model):
##    "Aircraft class"
##    def setup(self, **kwargs):
##        #create submodels
##        self.fuse = Fuselage()
##        self.wing = Wing()
##        self.engine = Engine()
##
##        #variable definitions
##        numeng = Variable('numeng', '-', 'Number of Engines')
##
##        constraints = [self.wing['x_w'] == self.fuse['l_{fuse}']*0.6,
##                       ]
##
##        constraints.extend([
##            numeng == numeng, #need numeng in the model
##            ])
##
##        self.components = [self.fuse, self.wing, self.engine]
##
##        return self.components, constraints
##
##    def dynamic(self, state):
##        """
##        creates an aircraft climb performance model, given a state
##        """
##        return AircraftP(self, state)
##        
##    def climb_dynamic(self, state):
##        """
##        creates an aircraft climb performance model, given a state
##        """
##        return ClimbP(self, state)
##
##    def cruise_dynamic(self, state):
##        """
##        creates an aircraft cruise performance model, given a state
##        """
##        return CruiseP(self, state)
##
##class AircraftP(Model):
##    """
##    aircraft performance models superclass, contains constraints true for
##    all flight segments
##    """
##    def  setup(self, aircraft, state, **kwargs):
##        #make submodels
##        self.aircraft = aircraft
##        self.wingP = aircraft.wing.dynamic(state)
##        self.fuseP = aircraft.fuse.dynamic(state)
##        self.engineP = aircraft.engine.dynamic(state)
##        self.wingP = aircraft.wing.dynamic(state)
##
##        self.Pmodels = [self.wingP, self.fuseP, self.engineP, self.wingP]
##
##        #variable definitions
##        Vstall = Variable('V_{stall}', 'knots', 'Aircraft Stall Speed')
##        D = Variable('D', 'N', 'Total Aircraft Drag')
##        W_avg = Variable('W_{avg}', 'N', 'Geometric Average of Segment Start and End Weight')
##        W_start = Variable('W_{start}', 'N', 'Segment Start Weight')
##        W_end = Variable('W_{end}', 'N', 'Segment End Weight')
##        W_burn = Variable('W_{burn}', 'N', 'Segment Fuel Burn Weight')
##        WLoadmax = Variable('W_{Load_{max}}', 'N/m^2', 'Max Wing Loading')
##        WLoad = Variable('W_{Load}', 'N/m^2', 'Wing Loading')
##        t = Variable('tmin', 'min', 'Segment Flight Time in Minutes')
##        thours = Variable('thr', 'hour', 'Segment Flight Time in Hours')
##
##        constraints = []
##
##        constraints.extend([
##            #speed must be greater than stall speed
##            state['V'] >= Vstall,
##
##
##            #Figure out how to delete
##            Vstall == 120*units('kts'),
##            WLoadmax == 6664 * units('N/m^2'),
##
##            #compute the drag
##            TCS([D >= self.wingP['D_{wing}'] + self.fuseP['D_{fuse}']]),
##
##            #constraint CL and compute the wing loading
##            W_avg == .5*self.wingP['C_{L}']*self.aircraft['S']*state.atm['\\rho']*state['V']**2,      
##            WLoad == .5*self.wingP['C_{L}']*self.aircraft['S']*state.atm['\\rho']*state['V']**2/self.aircraft.wing['S'],
##
##            #set average weight equal to the geometric avg of start and end weight
##            W_avg == (W_start * W_end)**.5,
##
##            #constrain the max wing loading
##            WLoad <= WLoadmax,
##
##            #compute fuel burn from TSFC
##            W_burn == aircraft['numeng']*self.engineP['TSFC'] * thours * self.engineP['F'],
##               
##            #time unit conversion
##            t == thours,
##            ])
##
##        return self.Pmodels, constraints
##    
##class Mission(Model):
##    """
##    mission class, links together all subclasses
##    """
##    def setup(self, subs = None, **kwargs):
##        #define the number of each flight segment
##        Nclimb = 2
##        Ncruise = 2
##
##        #build required submodels
##        aircraft = Aircraft()
##
##        #Vectorize
##
##        with Vectorize(Nclimb):
##            climb = ClimbSegment(aircraft)
##
##        with Vectorize(Ncruise):
##            cruise = CruiseSegment(aircraft)
##
##        #declare new variables
##        W_ftotal = Variable('W_{f_{total}}', 'N', 'Total Fuel Weight')
##        W_fclimb = Variable('W_{f_{climb}}', 'N', 'Fuel Weight Burned in Climb')
##        W_fcruise = Variable('W_{f_{cruise}}', 'N', 'Fuel Weight Burned in Cruise')
##        W_total = Variable('W_{total}', 'N', 'Total Aircraft Weight')
##        CruiseAlt = Variable('CruiseAlt', 'ft', 'Cruise Altitude [feet]')
##        ReqRng = Variable('R_{req}', 'nautical_miles', 'Required Cruise Range')
##
##        h = climb.state['h']
##        hftClimb = climb.state['hft']
##        dhft = climb.climbP['dhft']
##        hftCruise = cruise.state['hft']
##
##        #make overall constraints
##        constraints = []
##
##        constraints.extend([
##            # Wing max loading constraint
##            aircraft['L_{max}'] >= aircraft['N_{lift}'] * W_total,
##
##            #weight constraints
##            TCS([aircraft['W_{e}'] + aircraft['W_{payload}'] + W_ftotal + aircraft['numeng'] * aircraft['W_{engine}'] + aircraft.wing.wb['W_{struct}'] <= W_total]),
##
##            climb.climbP.aircraftP['W_{start}'][0] == W_total,
##            climb.climbP.aircraftP['W_{end}'][-1] == cruise.cruiseP.aircraftP['W_{start}'][0],
##
##            # similar constraint 1
##            TCS([climb.climbP.aircraftP['W_{start}'] >= climb.climbP.aircraftP['W_{end}'] + climb.climbP.aircraftP['W_{burn}']]),
##            # similar constraint 2
##            TCS([cruise.cruiseP.aircraftP['W_{start}'] >= cruise.cruiseP.aircraftP['W_{end}'] + cruise.cruiseP.aircraftP['W_{burn}']]),
##
##            climb.climbP.aircraftP['W_{start}'][1:] == climb.climbP.aircraftP['W_{end}'][:-1],
##            cruise.cruiseP.aircraftP['W_{start}'][1:] == cruise.cruiseP.aircraftP['W_{end}'][:-1],
##
##            TCS([aircraft['W_{e}'] + aircraft['W_{payload}'] + aircraft['numeng'] * aircraft['W_{engine}'] + aircraft.wing.wb['W_{struct}'] <= cruise.cruiseP.aircraftP['W_{end}'][-1]]),
##
##            TCS([W_ftotal >=  W_fclimb + W_fcruise]),
##            TCS([W_fclimb >= sum(climb.climbP['W_{burn}'])]),
##            TCS([W_fcruise >= sum(cruise.cruiseP['W_{burn}'])]),
##
##            #altitude constraints
##            hftCruise == CruiseAlt,
##            TCS([hftClimb[1:Ncruise] >= hftClimb[:Ncruise-1] + dhft]),
##            TCS([hftClimb[0] >= dhft[0]]),
##            hftClimb[-1] <= hftCruise,
##
##            #compute the dh
##            dhft == hftCruise/Nclimb,
##
##            #constrain the thrust
##            climb.climbP.engineP['F'] <= 2 * max(cruise.cruiseP.engineP['F']),
##
##            #set the range for each cruise segment, doesn't take credit for climb
##            #down range disatnce covered
##            cruise.cruiseP['Rng'] == ReqRng/(Ncruise),
##
##            #set the TSFC
##            climb.climbP.engineP['TSFC'] == .7*units('1/hr'),
##            cruise.cruiseP.engineP['TSFC'] == .5*units('1/hr'),
##
##            #wing constraints
##            aircraft.wing['W_{fuel_{wing}}'] == W_ftotal,
##            climb.climbP.wingP['L_w'] == climb.climbP.aircraftP['W_{avg}'],
##            cruise.cruiseP.wingP['L_w'] == cruise.cruiseP.aircraftP['W_{avg}'],
##            climb['c_{m_{w}}'] == .10, # for boundedness
##            cruise['c_{m_{w}}'] == .10, # for boundedness
##            ])
##        
##        return aircraft, climb, cruise, constraints
