"Implement HorizontalTail model"
from numpy import pi, tan
from gpkit import Variable, Model, SignomialsEnabled, LinkedConstraintSet, units, VectorVariable
from gpkit.constraints.costed import CostedConstraintSet
from gpkit.constraints.tight import TightConstraintSet as TCS
from wingbox import WingBox
from collections import defaultdict
from gpkit.small_scripts import mag
import numpy as np
# pylint:disable=bad-whitespace

class HorizontalTailNoStruct(Model):
    """
    horiziontal tail model from Philippe's thesis
    as a performance model without the wing box

    References:
    [1] TASOPT code
    [2] http://adg.stanford.edu/aa241/stability/staticstability.html
    """
    def __init__(self, fuse, wing,**kwargs):
        self.fuse = fuse
        self.wing = wing
        
        #variables
        p       = Variable('p_{ht}', '-', 'Substituted variable = 1 + 2*taper')
        q       = Variable('q_{ht}', '-', 'Substituted variable = 1 + taper')

        
        rho0    = Variable('\\rho_0', 1.225, 'kg/m^3', 'Air density (0 ft)')
        tanLh   = Variable('\\tan(\\Lambda_{ht})', '-',
                           'tangent of horizontal tail sweep')
        taper   = Variable('\lambda_h', '-', 'Horizontal tail taper ratio')
        tau     = Variable('\\tau_h', '-',
                           'Horizontal tail thickness/chord ratio')
        wf      = Variable('w_{fuse}', 'm', 'Fuselage width')
        xcght   = Variable('x_{CG_{ht}}', 'm', 'Horizontal tail CG location')
        ymac    = Variable('y_{\\bar{c}_{ht}}', 'm',
                           'Spanwise location of mean aerodynamic chord')
        lht     = Variable('l_{ht}', 'm', 'Horizontal tail moment arm')
        ARh     = Variable('AR_h', '-', 'Horizontal tail aspect ratio')
        Vne     = Variable('V_{ne}', 'm/s', 'Never exceed velocity')
        amax    = Variable('\\alpha_{max,h}', '-', 'Max angle of attack, htail')
        e       = Variable('e_h', '-', 'Oswald efficiency factor')
        Sw      = Variable('S_w', 'm^2', 'Wing area')
        Sh      = Variable('S_h', 'm^2', 'Horizontal tail area')
        bht     = Variable('b_{ht}', 'm', 'Horizontal tail span')
        chma    = Variable('\\bar{c}_{ht}', 'm', 'Mean aerodynamic chord (ht)')
        croot   = Variable('c_{root_h}', 'm', 'Horizontal tail root chord')
        ctip    = Variable('c_{tip_h}', 'm', 'Horizontal tail tip chord')
        cwma    = Variable('\\bar{c}_w', 'm',
                           'Mean aerodynamic chord (wing)')
        Lmax    = Variable('L_{{max}_h}', 'N', 'Maximum load')
        Kf      = Variable('K_f', '-',
                           'Empirical factor for fuselage-wing interference')
        fl      = Variable(r"f(\lambda_h)", '-',
                           'Empirical efficiency function of taper')
        SMmin   = Variable('S.M._{min}', '-', 'Minimum stability margin')
        CLhmax  = Variable('C_{L_{hmax}}', '-', 'Max lift coefficient')

        #constraints
        constraints = []

        with SignomialsEnabled():

            constraints.extend([
                # Moment arm and geometry -- same as for vtail
                p >= 1 + 2*taper,
                2*q >= 1 + p,
                ymac == (bht/3)*q/p,
                TCS([(2./3)*(1 + taper + taper**2)*croot/q >=
                     chma]), # [SP]
                taper == ctip/croot,
                TCS([Sh <= bht*(croot + ctip)/2]), # [SP]

                # K_f as f(wing position) -- (fitted posynomial)
                # from from UMich AE-481 course notes Table 9.1
                Kf >= (1.5012*(self.wing['x_w']/self.fuse['l_{fuse}'])**2 +
                       0.538*(self.wing['x_w']/self.fuse['l_{fuse}']) +
                       0.0331),

                # Oswald efficiency
                # Nita, Scholz,
                # "Estimating the Oswald factor from basic
                # aircraft geometrical parameters"
                TCS([fl >= (0.0524*taper**4 - 0.15*taper**3
                            + 0.1659*taper**2
                            - 0.0706*taper + 0.0119)], reltol=0.2),
                # NOTE: slightly slack
                TCS([e*(1 + fl*ARh) <= 1]),
                taper >= 0.2, # TODO: make less arbitrary

                Lmax == 0.5*rho0*self.wing['V_{ne}']**2*Sh*CLhmax,

                tanLh == tanLh,

                tau == tau,

                lht == lht,

                SMmin == SMmin,

                amax == amax,

                xcght == xcght,
                ])

        Model.__init__(self, None, constraints)


class HorizontalTailPerformance(Model):
    """
    Horizontal tail performance model
    """
    def __init__(self, ht, fuse, wing, state, **kwargs):
        self.ht = ht
        self.fuse = fuse
        self.wing = wing
        
        #variables
        alpha   = Variable('\\alpha', '-', 'Horizontal tail angle of attack')
        D       = Variable('D_{ht}', 'N', 'Horizontal tail drag')
        Lh      = Variable('L_h', 'N', 'Horizontal tail downforce')
        SM      = VectorVariable(2, 'S.M.', '-', 'Stability margin')
        Rec     = Variable('Re_{c_h}', '-',
                           'Cruise Reynolds number (Horizontal tail)')
        CLah    = Variable('C_{L_{ah}}', '-', 'Lift curve slope (htail)')
        CLah0   = Variable('C_{L_{ah_0}}', '-',
                           'Isolated lift curve slope (htail)')
        CLaw    = Variable('C_{L_{aw}}', '-', 'Lift curve slope (wing)')

        
        CLh     = Variable('C_{L_h}', '-', 'Lift coefficient (htail)')
        dxlead  = Variable('\\Delta x_{{lead}_h}', 'm',
                           'Distance from CG to horizontal tail leading edge')
        dxtrail = Variable('\\Delta x_{{trail}_h}', 'm',
                           'Distance from CG to horizontal tail trailing edge')
        dxw     = VectorVariable(2, '\\Delta x_w', 'm',
                           'Distance from aerodynamic centre to CG')
        
        etaht   = Variable('\\eta_{ht}', '-', 'Tail efficiency')
        eta     = Variable('\\eta_h', '-',
                           ("Lift efficiency (diff between sectional and "
                            "actual lift)"))
        xcg     = VectorVariable(2, 'x_{CG}', 'm', 'CG location')

        Cmac    = Variable('|C_{m_{ac}}|', '-', # Absolute value of CMwing
                           'Moment coefficient about aerodynamic centre (wing)')
        CDh     = Variable('C_{D_h}', '-', 'Horizontal tail drag coefficient')
        CD0h    = Variable('C_{D_{0_h}}', '-',
                           'Horizontal tail parasitic drag coefficient')
        
        #cosntraints
        constraints = []

        with SignomialsEnabled():
           
            constraints.extend([
                TCS([SM + dxw/self.wing['\\bar{c}_w'] + self.ht['K_f']*self.fuse['w_{fuse}']**2*self.fuse['l_{fuse}']/(self.wing['C_{L_{aw}}']*self.wing['S']*self.wing['\\bar{c}_w'])
                                    <= CLah*self.ht['S_h']*self.ht['l_{ht}']/(self.wing['C_{L_{aw}}']*self.wing['S']*self.wing['\\bar{c}_w'])]),
                
                SM >= self.ht['S.M._{min}'],

                # Trim from UMich AE-481 course notes
                TCS([CLh*self.ht['S_h']*self.ht['l_{ht}']/(self.wing['S']*self.wing['\\bar{c}_w']) + Cmac >=
                self.wing['C_{L}']*dxw/self.wing['\\bar{c}_w'] + self.fuse['C_{m_{fuse}}']], reltol=0.02),
                Lh == 0.5*state['\\rho']*state['V']**2*self.ht['S_h']*CLh,

                # Moment arm and geometry -- same as for vtail
                TCS([dxlead + self.ht['c_{root_h}'] <= dxtrail]),
                TCS([xcg + dxtrail <= self.fuse['l_{fuse}']], reltol=0.002),
                TCS([dxlead + self.ht['y_{\\bar{c}_{ht}}']*self.ht['\\tan(\\Lambda_{ht})'] + 0.25*self.ht['\\bar{c}_{ht}'] >= self.ht['l_{ht}']],
                    reltol=1e-2), # [SP]

                # DATCOM formula (Mach number makes it SP)
                TCS([(self.ht['AR_h']/eta)**2*(1+self.ht['\\tan(\\Lambda_{ht})']**2-state['M']**2) + 8*pi*self.ht['AR_h']/CLah0
                     <= (2*pi*self.ht['AR_h']/CLah0)**2]),

                # K_f as f(wing position) -- (fitted posynomial)
                # from from UMich AE-481 course notes Table 9.1
                self.wing['x_w'] >= xcg + dxw,

                # Loss of tail effectiveness due to wing downwash
                CLah + (2*self.wing['C_{L_{aw}}']/(pi*self.wing['AR']))*etaht*CLah0 <= CLah0*etaht,
                CLh == CLah*alpha,
                alpha <= self.ht['\\alpha_{max,h}'],

                # Drag
                D == 0.5*state['\\rho']*state['V']**2*self.ht['S_h']*CDh,
                CDh >= CD0h + CLh**2/(pi*self.ht['e_h']*self.ht['AR_h']),

                # same drag model as vtail
                CD0h**0.125 >= 0.19*(self.ht['\\tau_h'])**0.0075 *(Rec)**0.0017
                            + 1.83e+04*(self.ht['\\tau_h'])**3.54*(Rec)**-0.494
                            + 0.118*(self.ht['\\tau_h'])**0.0082 *(Rec)**0.00165
                            + 0.198*(self.ht['\\tau_h'])**0.00774*(Rec)**0.00168,
                Rec == state['\\rho']*state['V']*self.ht['\\bar{c}_{ht}']/state['\\mu'],

                self.ht['x_{CG_{ht}}'] >= xcg+(dxlead+dxtrail)/2,
                self.ht['x_{CG_{ht}}'] <= self.fuse['l_{fuse}'],

                dxw == [4, 2]*units('m'),
                xcg == [15, 17]*units('m'),
                ])

        Model.__init__(self, None, constraints)

class HorizontalTail(Model):
    """
    horiziontal tail model from Philippe's thesis
    as a performance model without the wing box

    References:
    [1] TASOPT code
    [2] http://adg.stanford.edu/aa241/stability/staticstability.html
    """
    def __init__(self, fuse, wing,**kwargs):
        self.fuse = fuse
        self.wing = wing
        
        self.htns = HorizontalTailNoStruct(self.fuse, self.wing)
        self.wb = WingBox(self.htns)

        Model.__init__(self, None, self.htns + self.wb)


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
        taper = Variable('taper', 0.2, '-', 'Taper ratio')
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

        if isinstance(surface, HorizontalTailNoStruct):
            AR = surface['AR_h']
            b = surface['b_{ht}']
            S = surface['S_h']
            p = surface['p_{ht}']
            q = surface['q_{ht}']
            tau = surface['\\tau_h']
            Lmax = surface['L_{{max}_h}']

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

class TestFuse(Model):
    """
    fuselage class used for testing purposes only
    """
    def __init__(self, **kwargs):
        #variables
        lfuse   = Variable('l_{fuse}', 'm', 'Fuselage length')
        Cmfu    = Variable('C_{m_{fuse}}', '-', 'Moment coefficient (fuselage)')
        wf      = Variable('w_{fuse}', 'm', 'Fuselage width')
        
        #constraints
        constraints = []

        constraints.extend([
            wf == wf,
            lfuse == lfuse,
            Cmfu == Cmfu,
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
            V == M * a,
            a  == (gamma * R * T_atm)**.5,
            a == 297 * units('m/s'),
            mu == mu,
            ])

        Model.__init__(self, None, constraints)

class TestWing(Model):
    """
    wing class for testing purposes
    """
    def __init__(self, **kwargs):
        #variables
        AR      = Variable('AR', '-', 'Wing aspect ratio')
        Vne     = Variable('V_{ne}', 'm/s', 'Never exceed velocity')
        CLw     = Variable('C_{L}', '-', 'Lift coefficient, wing')
        CLaw    = Variable('C_{L_{aw}}', '-', 'Lift curve slope, wing')
        cwma    = Variable('\\bar{c}_w', 'm',
                          'Mean aerodynamic chord (wing)')
        S       = Variable('S', 'm^2', 'Reference area')
        xw     = Variable('x_w', 'm', 'Position of wing aerodynamic center')
        
        #cosntraints
        constraints = []

        constraints.extend([
            AR == AR,
            Vne == Vne,
            CLw == CLw,
            CLaw == CLaw,
            cwma == cwma,
            S == S,
            xw==xw,
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

if __name__ == "__main__":
     #create the required submodels
     fuse = TestFuse()
     state = TestState()
     wing = TestWing()
     ht = HorizontalTail(fuse, wing)

     htP = ht.dynamic(fuse, wing, state)

     subs = {
             'AR': 9,
             'C_{L}': 0.5,
             'C_{L_{aw}}': 5,
             'C_{m_{fuse}}': 0.05, # [1]
             'M': 0.78,
             'S.M._{min}': 0.05,
             'S': 125,
             'V_{ne}': 144,
             'C_{L_{hmax}}': 2.5,

##             '\\Delta x_w': [4, 2],
             '\\alpha_{max,h}': 0.1, # (6 deg)
             '\\bar{c}_w': 5,
             '\\eta_h': 0.97, # [2]
             '\\eta_{ht}': 0.9, # [2]
             '\\mu': 1.4E-5,
             '\\rho': 0.38,
             '\\rho_0': 1.225,
             '\\tan(\\Lambda_{ht})': tan(30*pi/180),
             'a': 297,
             'l_{fuse}': 40,
             'w_{fuse}': 6,
##             'x_{CG}': [15, 17],
             '|C_{m_{ac}}|': 0.1,
            }

     m = Model(htP['D_{ht}'] + 0.1*ht.wb['W_{struct}'], [htP, ht, wing, state, fuse], subs)
     sol = m.localsolve(solver='mosek', verbosity=4)
##     bounds, sol = wing.determine_unbounded_variables(m, solver="mosek",verbosity=4, iteration_limit=100)

