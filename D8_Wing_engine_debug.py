"""Simple commercial aircraft flight profile and aircraft model"""
from numpy import pi, cos, tan
import numpy as np
from gpkit import Variable, Model, units, SignomialsEnabled, Vectorize
from gpkit.constraints.sigeq import SignomialEquality
from gpkit.tools import te_exp_minus1
from gpkit.constraints.tight import Tight as TCS
import matplotlib.pyplot as plt
from simple_ac_imports import Aircraft, Fuselage, CruiseP, ClimbP, FlightState, CruiseSegment, ClimbSegment
from engine_validation import Engine

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
    
class Mission(Model):
    """
    mission class, links together all subclasses
    """
    def setup(self, subs = None, **kwargs):
        #define the number of each flight segment
        Nclimb = 2
        Ncruise = 2

        # vectorize
        with Vectorize(Nclimb + Ncruise):
            enginestate = FlightState()

        #build required submodels
        aircraft = Aircraft(Nclimb, Ncruise, enginestate, 0)

        #Vectorize

        with Vectorize(Nclimb):
            climb = ClimbSegment(aircraft)

        with Vectorize(Ncruise):
            cruise = CruiseSegment(aircraft)

        statelinking = StateLinking(climb.state, cruise.state, enginestate, Nclimb, Ncruise)

        #declare new variables
        W_ftotal = Variable('W_{f_{total}}', 'N', 'Total Fuel Weight')
        W_fclimb = Variable('W_{f_{climb}}', 'N', 'Fuel Weight Burned in Climb')
        W_fcruise = Variable('W_{f_{cruise}}', 'N', 'Fuel Weight Burned in Cruise')
        W_total = Variable('W_{total}', 'N', 'Total Aircraft Weight')
        CruiseAlt = Variable('CruiseAlt', 'ft', 'Cruise Altitude [feet]')
        ReqRng = Variable('ReqRng', 'nautical_miles', 'Required Cruise Range')
        RngCruise = Variable('RngCruise', 'nautical_miles', 'Total Cruise Range')

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

            #set the range for each cruise segment, doesn't take credit for climb
            #down range disatnce covered
            cruise.cruiseP['Rng'] == RngCruise/(Ncruise),

            #wing constraints
            aircraft.wing['W_{fuel_{wing}}'] == W_ftotal,
            climb.climbP.wingP['L_w'] == climb.climbP.aircraftP['W_{avg}'],
            cruise.cruiseP.wingP['L_w'] == cruise.cruiseP.aircraftP['W_{avg}'],
            climb['c_{m_{w}}'] == .10, # for boundedness
            cruise['c_{m_{w}}'] == .10, # for boundedness

            cruise['hft'] >= 1000*units('ft'),
            ])

        with SignomialsEnabled():
            constraints.extend([
                sum(climb['RngClimb']) + RngCruise >= ReqRng
                ])

        M2 = .8
        M25 = .6
        M4a = .1025
        M0 = .8

        engineclimb = [
            aircraft.engine.engineP['M_2'][:Nclimb] == climb['M'],
            aircraft.engine.engineP['M_{2.5}'][:Nclimb] == M25,
            aircraft.engine.engineP['hold_{2}'] == 1+.5*(1.398-1)*M2**2,
            aircraft.engine.engineP['hold_{2.5}'] == 1+.5*(1.354-1)*M25**2,
            aircraft.engine.engineP['c1'] == 1+.5*(.401)*M0**2,

            #constraint on drag and thrust
            aircraft['numeng']*aircraft.engine['F_{spec}'][:Nclimb] >= climb['D'] + climb['W_{avg}'] * climb['\\theta'],

            #climb rate constraints
            TCS([climb['excessP'] + climb.state['V'] * climb['D'] <=  climb.state['V'] * aircraft['numeng'] * aircraft.engine['F_{spec}'][:Nclimb]]),
            ]

        M25 = .6

        enginecruise = [
            aircraft.engine.engineP['M_2'][Nclimb:] == cruise['M'],
            aircraft.engine.engineP['M_{2.5}'][Nclimb:] == M25,
            
            #steady level flight constraint on D 
            cruise['D'] == aircraft['numeng'] * aircraft.engine['F_{spec}'][Nclimb:],

            #breguet range eqn
            TCS([cruise['z_{bre}'] >= (aircraft.engine['TSFC'][Nclimb:] * cruise['thr']*
            cruise['D']) / cruise['W_{avg}']]),
            ]
        
        return aircraft, climb, cruise, constraints, statelinking, enginestate, engineclimb, enginecruise

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
    Philippe's wing model minus structure
    """
    def setup(self, **kwargs):
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
        wwn       = Variable('wwn', 0.5, '-', 'Wingbox-width-to-chord ratio')
        xw     = Variable('x_w', 'm', 'Position of wing aerodynamic center')
        ymac    = Variable('y_{mac}', 'm',
                           'Spanwise location of mean aerodynamic chord')

        #Linked Variables
        AR      = Variable('AR', '-', 'Wing aspect ratio')
        Lmax    = Variable('L_{max}', 'N', 'Maximum load')
        Sw      = Variable('S', 'm^2', 'Wing area')
##        Vne     = Variable('V_{ne}', 144, 'm/s', 'Never exceed velocity')
        WfuelWing   = Variable('W_{fuel_{wing}}', 'N', 'Fuel weight')
        b       = Variable('b', 'm', 'Wing span')
        #the following two variables have the same name in the flight profile and
        #will be automatically linked by the linked constraint set
        mac    = Variable('mac', 'm',
                          'Mean aerodynamic chord (wing)')
        e       = Variable('e', '-', 'Oswald efficiency factor')

        
        #make cosntraints
        constraints = []

        with SignomialsEnabled():

            constraints.extend([
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
                taper >= 0.2, # TODO

                # Fuel volume [TASOPT doc]
##                TCS([Afuel <= wwn*0.92*tau]),
                # GP approx of the signomial constraint:
                # Afuel <= (w - 2*tweb)*(0.92*tau - 2*tcap),
##                Vfuel <= croot**2 * (b/6) * (1+taper+taper**2)*cosL,
                Vfuel <= croot*b*taper*wwn*tau*croot,
                WfuelWing <= rhofuel*Afuel*Vfuel*g,
                  
                # Lmax == 0.5*rho0*Vne**2*Sw*CLwmax,
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

        #wing moment variables -- need a good way to model this, currently using TAT
        cmw = Variable('c_{m_{w}}', '-', 'Wing Pitching Moment Coefficient')
        
        #make constraints
        constraints = []

        with SignomialsEnabled():
            constraints.extend([
                Lw == 0.5*state['\\rho']*state['V']**2*self.wing['S']*CLw,

                # DATCOM formula (Mach number makes it SP)
                # Swept wing lift curve slope constraint
                SignomialEquality((self.wing['AR']/self.wing['\\eta'])**2*(1 + self.wing['\\tan(\\Lambda)']**2 - state['M']**2) + 8*pi*self.wing['AR']/CLaw
                      , (2*pi*self.wing['AR']/CLaw)**2),
                
                CLw == CLaw*alpha,
                alpha <= self.wing['\\alpha_{max,w}'],

                # Drag
                D == 0.5*state['\\rho']*state['V']**2*self.wing['S']*CDw,
                TCS([CDw >= CDp + CLw**2/(pi*self.wing['e']*self.wing['AR'])]),
                Re == state['\\rho']*state['V']*self.wing['mac']/state['\\mu'],
                1 >= (2.56*CLw**5.88/(Re**1.54*self.wing['\\tau']**3.32*CDp**2.62)
                   + 3.8e-9*self.wing['\\tau']**6.23/(CLw**0.92*Re**1.38*CDp**9.57)
                   + 2.2e-3*Re**0.14*self.wing['\\tau']**0.033/(CLw**0.01*CDp**0.73)
                   + 6.14e-6*CLw**6.53/(Re**0.99*self.wing['\\tau']**0.52*CDp**5.19)
                   + 1.19e4*CLw**9.78*self.wing['\\tau']**1.76/(Re*CDp**0.91)),
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
                       ]
        
        return constraints

if __name__ == '__main__':
    plot = False

    sweep = 30 #[deg]

    M4a = .1025
    fan = 1.685
    lpc  = 1.935
    hpc = 9.369

    D80 = False
    
    substitutions = {      
            'ReqRng': 2656, #('sweep', np.linspace(500,2000,4)),
##            'CruiseAlt': 30000, #('sweep', np.linspace(20000,40000,4)),
            'numeng': 2,
            'W_{pax}': 91 * 9.81,
            'n_{pax}': 150,
            'pax_{area}': 1,

            #wing subs
            'C_{L_{wmax}}': 2.5,
##            'V_{ne}': 144,
            '\\alpha_{max,w}': 0.1, # (6 deg)
            '\\cos(\\Lambda)': cos(sweep*pi/180),
            '\\eta': 0.97,
##            '\\rho_0': 1.225,
            '\\rho_{fuel}': 817, # Kerosene [TASOPT]
            '\\tan(\\Lambda)': tan(sweep*pi/180),

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
            '\\alpha_{max}': 5.105,

            'hold_{4a}': 1+.5*(1.313-1)*M4a**2,
            'r_{uc}': .01,
            '\\alpha_c': .19036,
            'T_{t_f}': 435,

            'M_{takeoff}': .9556,

            'G_f': 1,

            'h_f': 43.003,

            'Cp_t1': 1280,
            'Cp_t2': 1184,
            'Cp_c': 1216,
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
##            'V_{ne}': 144,
            '\\alpha_{max,w}': 0.1, # (6 deg)
            '\\cos(\\Lambda)': cos(sweep*pi/180),
            '\\eta': 0.97,
##            '\\rho_0': 1.225,
            '\\rho_{fuel}': 817, # Kerosene [TASOPT]
            '\\tan(\\Lambda)': tan(sweep*pi/180),

            'AR':10.8730,
            'b':116.548*0.3048,#units('ft'),
            'c_{root}': 17.4*0.3048,#units('ft'),

            # Minimum Cruise Mach Number
            # 'M': 0.8,
        })
    mission = Mission()
    m = Model(mission['W_{f_{total}}'], mission, substitutions)
    sol = m.localsolve(solver='mosek', verbosity = 4)

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

