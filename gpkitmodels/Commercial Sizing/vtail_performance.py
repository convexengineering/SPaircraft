# coding=utf-8
"Implements a Vertical Tail model"
import numpy as np
from gpkit import Model, Variable, SignomialsEnabled, LinkedConstraintSet, units
from gpkit.constraints.costed import CostedConstraintSet
from gpkit.constraints.tight import TightConstraintSet as TCS
from wingbox import WingBox
from collections import defaultdict
from gpkit.small_scripts import mag
import numpy as np

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
    def __init__(self, fuse, engine, **kwargs):
        self.fuse = fuse
        self.engine = engine
        
        self.vtns = VerticalTailNoStruct(self.fuse, self.engine)
        self.wb = WingBox(self.vtns)

        Model.__init__(self, None, self.vtns + self.wb)

    def dynamic(self, state):
        """"
        creates a horizontal tail performance model
        """
        return VerticalTailPerformance(self, state)

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
    def __init__(self, fuse, engine, **kwargs):
        self.fuse = fuse
        self.engine = engine
        
        #define new variables
        Avt    = Variable('A_{vt}', '-', 'Vertical tail aspect ratio')
        CDwm   = Variable('C_{D_{wm}}', '-', 'Windmill drag coefficient')
        Dwm    = Variable('D_{wm}', 'N', 'Engine out windmill drag')
        Lmax   = Variable('L_{max_{vt}}', 'N',
                          'Maximum load for structural sizing') # fuselage
        Lvmax  = Variable('L_{v_{max}}', 'N',
                          'Maximum load for structural sizing')
        CLvmax = Variable('C_{L_{vmax}}', '-', 'Max lift coefficient')
        CLvtEO   = Variable('C_{L_{vtEO}}', '-', 'Vertical tail lift coefficient (Engine Out)')
        clvtEO   = Variable('c_{l_{vtEO}}', '-',
                            'Sectional lift force coefficient (engine out)')
        LvtEO    = Variable('L_{vtEO}', 'N', 'Vertical tail lift in engine out')
        S      = Variable('S', 'm^2', 'Vertical tail reference area (full)')
        Svt    = Variable('S_{vt}', 'm^2', 'Vertical tail reference area (half)')
        V1     = Variable('V_1', 'm/s', 'Minimum takeoff velocity')
        Vne    = Variable('V_{ne}', 'm/s', 'Never exceed velocity')
        Wstruct= Variable('W_{struct}', 'N', 'Full span weight')
        Wvt    = Variable('W_{vt}', 'N', 'Vertical tail weight')
        b      = Variable('b', 'm', 'Vertical tail full span')
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

        #constraints
        constraints = []

        with SignomialsEnabled():
            #non vectorized constraints
            constraints.extend([
                LvtEO*lvt >= Te*y_eng + Dwm*y_eng,
                # Force moment balance for one engine out condition
                # TASOPT 2.0 p45

                TCS([dxlead + zmac*tanL + 0.25*cma >= lvt]), # [SP]
                # Tail moment arm

                LvtEO == 0.5*rho0*V1**2*Svt*CLvtEO,
                # Vertical tail force (y-direction) for engine out

                TCS([CLvtEO*(1 + clvtEO/(np.pi*e*Avt)) <= clvtEO]),
                #engine out CL computation

                Avt == bvt**2/Svt,

                Dwm >= 0.5*rho0*V1**2*self.engine['A_2']*CDwm,
                # Drag of a windmilling engine
                          
                Svt <= bvt*(croot + ctip)/2, # [SP]
                # Tail geometry relationship

                TCS([dxtrail >= croot + dxlead]),
                # Tail geometry constraint

                self.fuse['l_{fuse}'] >= dxtrail + self.fuse['x_{CG}'],
                # Fuselage length constrains the tail trailing edge

                TCS([p >= 1 + 2*taper]),
                TCS([2*q >= 1 + p]),
                zmac == (bvt/3)*q/p,
                TCS([(2./3)*(1 + taper + taper**2)*croot/q >= cma]), # [SP]
                taper == ctip/croot,
                # Define vertical tail geometry

                Lvmax == 0.5*rho0*Vne**2*Svt*CLvmax,
                #compute the max force

                S == Svt*2,
                b == bvt*2,
                Lmax == 2*Lvmax,
                # Relate vertical tail geometry/weight to generic
                # wing used in structural model
     
                taper >= 0.27,
                # TODO: Constrain taper by tip Reynolds number
                # source: b737.org.uk

                TCS([xCGvt >= self.fuse['x_{CG}']+(dxlead+dxtrail)/2],
                         raiseerror=False),
                xCGvt <= self.fuse['l_{fuse}'],

                tau == tau,
                ])

        Model.__init__(self, None, constraints)

        

class VerticalTailPerformance(Model):
    """
    Vertical tail perofrmance model
    """
    def __init__(self, vt, state):
        self.vt = vt

        #define new variables
 
        CLvt   = Variable('C_{L_{vt}}', '-', 'Vertical tail lift coefficient')
        clvt   = Variable('c_{l_{vt}}', '-',
                          'Sectional lift force coefficient')
        Dvt    = Variable('D_{vt}', 'N', 'Vertical tail viscous drag, cruise')
        Rec    = Variable('Re_{vt}', '-', 'Vertical tail reynolds number, cruise')
        Vinf   = Variable('V_{\\infty}', 'm/s', 'Cruise velocity')
        xCG    = Variable('x_{CG}', 'm', 'x-location of CG')
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

            Dvt >= 0.5*state['\\rho']*state['V']**2*self.vt['S_{vt}']*CDvis,
            CDvis**0.125 >= 0.19*(self.vt['\\tau_{vt}'])**0.0075 *(Rec)**0.0017
                        + 1.83e+04*(self.vt['\\tau_{vt}'])**3.54*(Rec)**-0.494
                        + 0.118*(self.vt['\\tau_{vt}'])**0.0082 *(Rec)**0.00165
                        + 0.198*(self.vt['\\tau_{vt}'])**0.00774*(Rec)**0.00168,
            # Vertical tail viscous drag in cruise
            # Data fit from Xfoil

            Rec == state['\\rho']*state['V']*self.vt['\\bar{c}_{vt}']/state['\\mu'],
            # Cruise Reynolds number
            ])

        Model.__init__(self, None, constraints)


class TestFuse(Model):
    """
    fuselage class used for testing purposes only
    """
    def __init__(self, **kwargs):
        #variables
        lfuse   = Variable('l_{fuse}', 'm', 'Fuselage length')
        xCG    = Variable('x_{CG}', 'm', 'x-location of CG')
        #constraints
        constraints = []

        constraints.extend([
            lfuse == lfuse,
            xCG == xCG,
            ])

        Model.__init__(self, None, constraints)

class TestEngine(Model):
    """
    fuselage class used for testing purposes only
    """
    def __init__(self, **kwargs):
        #variables
        A2 = Variable('A_2', 'm^2', 'Fan Area')
        
        #constraints
        constraints = []

        constraints.extend([
            A2 == A2,
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
            rho == .38*units('kg/m^3'),
            V == V,
##            V == M * a,
##            a  == (gamma * R * T_atm)**.5,
##            a == 297 * units('m/s'),
            mu == mu,
            ])

        Model.__init__(self, None, constraints)

    def bound_all_variables(self, model, eps=1e-30, lower=None, upper=None):
        "Returns model with additional constraints bounding all free variables"
        lb = lower if lower else eps
        ub = upper if upper else 1/eps
        constraints = []
        freevks = tuple(vk for vk in model.varkeys if "value" not in vk.descr)
        for varkey in freevks:
            units = varkey.descr.get("units", 1)
            constraints.append([ub*units >= Variable(**varkey.descr),
                                Variable(**varkey.descr) >= lb*units])
        m = Model(model.cost, [constraints, model], model.substitutions)
        m.bound_all = {"lb": lb, "ub": ub, "varkeys": freevks}
        return m

    # pylint: disable=too-many-locals
    def determine_unbounded_variables(self, model, solver=None, verbosity=0,
                                      eps=1e-30, lower=None, upper=None, **kwargs):
        "Returns labeled dictionary of unbounded variables."
        m = self.bound_all_variables(model, eps, lower, upper)
        sol = m.localsolve(solver, verbosity, **kwargs)
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

        if isinstance(surface, VerticalTailNoStruct):
            AR = 2*surface['A_{vt}']
            b = surface['b']
            S = surface['S']
            p = surface['p_{vt}']
            q = surface['q_{vt}']
            tau = surface['\\tau_{vt}']
            Lmax = surface['L_{max_{vt}}']
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
        
        Model.__init__(self, objective, constraints, **kwargs)


if __name__ == "__main__":
    state = TestState()
    engine = TestEngine()
    fuse = TestFuse()
    vt = VerticalTail(fuse, engine)
    vtP = vt.dynamic(state)

    subs = {
           'C_{D_{wm}}': 0.5, # [2]
           'C_{L_{vmax}}': 2.6, # [2]
           'T_e': 1.29e5, # [4]
           'V_1': 70,
           'V': 454.86, # [7]
           'V_{ne}': 144, # [2]
           '\\mu': 1.4e-5, # [5]
           '\\rho': 0.38, # [6]
           '\\rho_{TO}': 1.225,
           '\\tan(\\Lambda_{vt})': np.tan(40*np.pi/180),
##           'c_{l_{vt}}': 0.5, # [2]
           'c_{l_{vtEO}}': 0.5,
           'A_2': np.pi*(.5*1.75)**2, # [1]
           'e_v': 0.8,
           'l_{fuse}': 39,
           'x_{CG}': 18,
           'y_{eng}': 4.83, # [3]
           }

    m = Model( vtP['D_{vt}'] + 0.05*.5*vt.wb['W_{struct}'], state + engine + fuse + vt + vtP, subs)
##    sol = m.localsolve(solver="mosek", verbosity = 4)
    bounds, sol = state.determine_unbounded_variables(m, solver="mosek",verbosity=4, iteration_limit=100)
