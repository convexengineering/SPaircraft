from gpkit import Variable, Model, SignomialsEnabled, LinkedConstraintSet, VectorVariable, units
from gpkit.constraints.costed import CostedConstraintSet
from gpkit.constraints.tight import TightConstraintSet as TCS
from numpy import cos, pi, tan
from wingbox import WingBox
from collections import defaultdict
from gpkit.small_scripts import mag
import numpy as np
# pylint:disable=bad-whitespace

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
        #xw     = Variable('x_w', 'm', 'Position of wing aerodynamic center')
        ymac    = Variable('y_{\\bar{c}_w}', 'm',
                           'Spanwise location of mean aerodynamic chord')

        #Linked Variables
        AR      = Variable('AR', '-', 'Wing aspect ratio')
        Lmax    = Variable('L_{max}', 'N', 'Maximum load')
        Sw      = Variable('S', 'm^2', 'Wing area')
        Vne     = Variable('V_{ne}', 'm/s', 'Never exceed velocity')
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

                 TCS([Sw <= b*(croot + ctip)/2], reltol=1E-2), # [SP]

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
                TCS([(self.wing['AR']/self.wing['\\eta'])**2*(1 + self.wing['\\tan(\\Lambda)']**2 - state['M']**2) + 8*pi*self.wing['AR']/CLaw
                      <= (2*pi*self.wing['AR']/CLaw)**2]),
                CLw == CLaw*alpha,
                alpha <= self.wing['\\alpha_{max,w}'],

                # Drag
                D == 0.5*state['\\rho']*state['V']**2*self.wing['S']*CDw,
                CDw >= CDp + CLw**2/(pi*self.wing['e']*self.wing['AR']),
                Re == state['\\rho']*state['V']*self.wing['\\bar{c}_w']/state['\\mu'],
                1 >= (2.56*CLw**5.88/(Re**1.54*self.wing['\\tau']**3.32*CDp**2.62)
                   + 3.8e-9*self.wing['\\tau']**6.23/(CLw**0.92*Re**1.38*CDp**9.57)
                   + 2.2e-3*Re**0.14*self.wing['\\tau']**0.033/(CLw**0.01*CDp**0.73)
                   + 6.14e-6*CLw**6.53/(Re**0.99*self.wing['\\tau']**0.52*CDp**5.19)
                   + 1.19e4*CLw**9.78*self.wing['\\tau']**1.76/(Re*CDp**0.91)),
                ])

        Model.__init__(self, CDw, constraints)

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

class StandAlone(Model):
    """
    constraints needed to make the model run alone
    """
    def __init__(self, wb, wingP):
        constraints = []

        constraints.append([
            wingP['L_w'] >= wb['W_{struct}'] + (5E5 + 1E5) * units('N')
            ])

        Model.__init__(self, None, constraints)

class TestState(Model):
    """
    state class only to be used for testing purposes
    """
    def __init__(self):
        #define variables
        p_atm = Variable("P_{atm}", "Pa", "air pressure")
        R_atm = Variable("R_{atm}", "J/mol/K", "air specific heating value")
        TH = 5.257386998354459 #(g*M_atm/R_atm/L_atm).value
        rho = Variable('\\rho', 'kg/m^3', 'Density of air')
        T_atm = Variable("T_{atm}", "K", "air temperature")
  

        mu  = Variable('\\mu', 'kg/(m*s)', 'Dynamic viscosity')

        V = Variable('V', 'kts', 'Aircraft Flight Speed')
        a = Variable('a', 'm/s', 'Speed of Sound')
        h = Variable('h', 'm', 'Segment Altitude [meters]')
        hft = Variable('hft', 'feet', 'Segment Altitude [feet]')
        R = Variable('R', 287, 'J/kg/K', 'Air Specific Heat')
        gamma = Variable('\\gamma', 1.4, '-', 'Air Specific Heat Ratio')
        M = Variable('M', '-', 'Mach Number')

        #make constraints
        constraints = []

        constraints.extend([
            rho == rho,
            V == V,
            V == M * a,
            a  == (gamma * R * T_atm)**.5,
            a == 297 * units('m/s'),
            mu == mu,
            ])

        Model.__init__(self, None, constraints)


if __name__ == "__main__":
    wing = Wing()
    state= TestState()

    wingP = wing.dynamic(state)


    standalone = StandAlone(wing.wb, wingP)

    sweep = 30 #[deg]

    subs = {
                 'C_{L_{wmax}}': 2.5,
                 'M': 0.78,
                 'V_{ne}': 144,
##                 'W_0': 5E5,
                 'W_{fuel_{wing}}': 1E5,
                 '\\alpha_{max,w}': 0.1, # (6 deg)
                 '\\cos(\\Lambda)': cos(sweep*pi/180),
                 '\\eta': 0.97,
                 '\\mu': 1.4E-5,
                 '\\rho': 0.38,
                 '\\rho_0': 1.225,
                 '\\rho_{fuel}': 817, # Kerosene [TASOPT]
                 '\\tan(\\Lambda)': tan(sweep*pi/180),
                 'a': 297,
##                 'g': 9.81,
##                 'C_{L}': 0.5
##                 'L_w': 1.296e+06,
##                 'I_{cap}': .0001,
                }
    m = Model(wingP['D_{wing}'], [wing, state, wingP, standalone])
    m.substitutions.update(subs)
    sol = m.localsolve(solver='mosek', verbosity = 4)
##    bounds, sol = state.determine_unbounded_variables(m, solver="mosek",verbosity=4, iteration_limit=100)
