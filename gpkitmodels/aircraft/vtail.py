# coding=utf-8
"Implements a Vertical Tail model"
import numpy as np
from gpkit import Model, Variable, SignomialsEnabled, LinkedConstraintSet
from gpkit.constraints.costed import CostedConstraintSet
from gpkit.constraints.tight import TightConstraintSet as TCS
from wingbox import WingBox

class VerticalTail(CostedConstraintSet):
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
    def __init__(self, **kwargs):

        Afan   = Variable('A_{fan}', 'm^2', 'Engine reference area')
        Avt    = Variable('A_{vt}', '-', 'Vertical tail aspect ratio')
        CDvis  = Variable('C_{D_{vis}}', '-', 'Viscous drag coefficient')
        CDwm   = Variable('C_{D_{wm}}', '-', 'Windmill drag coefficient')
        CLvmax = Variable('C_{L_{vmax}}', '-', 'Max lift coefficient')
        CLvt   = Variable('C_{L_{vt}}', '-', 'Vertical tail lift coefficient')
        Dvt    = Variable('D_{vt}', 'N', 'Vertical tail viscous drag, cruise')
        Dwm    = Variable('D_{wm}', 'N', 'Engine out windmill drag')
        Lmax   = Variable('L_{max_{vt}}', 'N',
                          'Maximum load for structural sizing') # fuselage
        Lvmax  = Variable('L_{v_{max}}', 'N',
                          'Maximum load for structural sizing')
        Lvt    = Variable('L_{vt}', 'N', 'Vertical tail lift in engine out')
        Rec    = Variable('Re_{vt}', '-', 'Vertical tail reynolds number, cruise')
        S      = Variable('S', 'm^2', 'Vertical tail reference area (full)')
        Svt    = Variable('S_{vt}', 'm^2', 'Vertical tail ref. area (half)')
        Te     = Variable('T_e', 'N', 'Thrust per engine at takeoff')
        V1     = Variable('V_1', 'm/s', 'Minimum takeoff velocity')
        Vc     = Variable('V_c', 'm/s', 'Cruise velocity')
        Vne    = Variable('V_{ne}', 'm/s', 'Never exceed velocity')
        Wstruct= Variable('W_{struct}', 'N', 'Full span weight')
        Wvt    = Variable('W_{vt}', 'N', 'Vertical tail weight')
        b      = Variable('b', 'm', 'Vertical tail full span')
        bvt    = Variable('b_{vt}', 'm', 'Vertical tail half span')
        clvt   = Variable('c_{l_{vt}}', '-',
                          'Sectional lift force coefficient (engine out)')
        cma    = Variable('\\bar{c}_{vt}', 'm', 'Vertical tail mean aero chord')
        croot  = Variable('c_{root_{vt}}', 'm', 'Vertical tail root chord')
        ctip   = Variable('c_{tip_{vt}}', 'm', 'Vertical tail tip chord')
        cvt    = Variable('c_{vt}', 'm', 'Vertical tail root chord') # fuselage
        dfan   = Variable('d_{fan}', 'm', 'Fan diameter')
        dxlead = Variable('\\Delta x_{lead}', 'm',
                          'Distance from CG to vertical tail leading edge')
        dxtrail= Variable('\\Delta x_{trail}', 'm',
                          'Distance from CG to vertical tail trailing edge')
        e      = Variable('e', '-', 'Span efficiency of vertical tail')
        lfuse  = Variable('l_{fuse}', 'm', 'Length of fuselage')
        lvt    = Variable('l_{vt}', 'm', 'Vertical tail moment arm')
        mu     = Variable('\\mu', 'N*s/m^2', 'Dynamic viscosity (35,000ft)')
        mu0    = Variable('\\mu_0', 1.8E-5, 'N*s/m^2', 'Dynamic viscosity (SL)')
        p      = Variable('p_{vt}', '-', 'Substituted variable = 1 + 2*taper')
        plamv  = Variable('p_{\\lambda_v}', '-',
                          'Dummy variable = 1 + 2\\lambda') # fuselage
        q      = Variable('q_{vt}', '-', 'Substituted variable = 1 + taper')
        rho0   = Variable('\\rho_{TO}', 'kg/m^3', 'Air density (SL))')
        rho_c  = Variable('\\rho_c', 'kg/m^3', 'Air density (35,000ft)')
        tanL   = Variable('\\tan(\\Lambda_{LE})', '-',
                          'Tangent of leading edge sweep (40 deg)')
        taper  = Variable('\\lambda_{vt}', '-', 'Vertical tail taper ratio')
        tau    = Variable('\\tau_{vt}', '-', 'Vertical tail thickness/chord ratio')
        xCG    = Variable('x_{CG}', 'm', 'x-location of CG')
        xCGvt  = Variable('x_{CG_{vt}}', 'm', 'x-location of tail CG')
        y_eng  = Variable('y_{eng}', 'm', 'Engine moment arm')
        zmac   = Variable('z_{\\bar{c}_{vt}}', 'm',
                          'Vertical location of mean aerodynamic chord')


        with SignomialsEnabled():
            objective = Dvt + 0.5*Wvt
            constraints = [
                           Lvt*lvt >= Te*y_eng + Dwm*y_eng,
                           # Force moment balance for one engine out condition
                           # TASOPT 2.0 p45

                           TCS([dxlead + zmac*tanL + 0.25*cma >= lvt]), # [SP]
                           # Tail moment arm

                           Lvt == 0.5*rho0*V1**2*Svt*CLvt,
                           # Vertical tail force (y-direction) for engine out 

                           Avt == bvt**2/Svt,
                           TCS([CLvt*(1 + clvt/(np.pi*e*Avt)) <= clvt]),
                           # Finite wing theory
                           # people.clarkson.edu/~pmarzocc/AE429/AE-429-4.pdf
                           # Valid because tail is untwisted and uncambered
                           # (lift curve slope passes through origin)

                           Dwm >= 0.5*rho0*V1**2*Afan*CDwm,
                           Afan >= np.pi*(dfan/2)**2,
                           # Drag of a windmilling engine

                           Svt <= bvt*(croot + ctip)/2, # [SP]
                           # Tail geometry relationship

                           TCS([dxtrail >= croot + dxlead]),
                           # Tail geometry constraint

                           lfuse >= dxtrail + xCG,
                           # Fuselage length constrains the tail trailing edge

                           p >= 1 + 2*taper,
                           2*q >= 1 + p,
                           zmac == (bvt/3)*q/p,
                           TCS([(2./3)*(1 + taper + taper**2)*croot/q >= cma]), # [SP]
                           taper == ctip/croot,
                           # Define vertical tail geometry

                           taper >= 0.27,
                           # TODO: Constrain taper by tip Reynolds number
                           # source: b737.org.uk

                           Dvt >= 0.5*rho_c*Vc**2*Svt*CDvis,
                           CDvis**0.125 >= 0.19*(tau)**0.0075 *(Rec)**0.0017
                                        + 1.83e+04*(tau)**3.54*(Rec)**-0.494
                                        + 0.118*(tau)**0.0082 *(Rec)**0.00165
                                        + 0.198*(tau)**0.00774*(Rec)**0.00168,
                           # Vertical tail viscous drag in cruise
                           # Data fit from Xfoil

                           Rec == rho_c*Vc*cma/mu,
                           # Cruise Reynolds number

                           S == Svt*2,
                           b == bvt*2,
                           Wvt == Wstruct/2,
                           Lmax == 2*Lvmax,
                           # Relate vertical tail geometry/weight to generic
                           # wing used in structural model

                           Lvmax == 0.5*rho0*Vne**2*Svt*CLvmax,
                           # Max load for structural sizing
                          ]

            # For linking to other models
            linking_constraints = [plamv == p,
                                   cvt == croot]

            CG_constraint = [TCS([xCGvt >= xCG+(dxlead+dxtrail)/2],
                                 raiseerror=False),
                             xCGvt <= lfuse]

            self.linking_constraints = linking_constraints
            self.CG_constraint = CG_constraint


        # Incorporate the structural model
        wb = WingBox()
        wb.subinplace({'L_{max}': Lmax,
                       'p': p,
                       'q': q,
                       'S': S,
                       'taper': taper,
                       '\\tau': tau})
        lc = LinkedConstraintSet([constraints, wb])

        CostedConstraintSet.__init__(self, objective, lc)

    @classmethod
    def standalone_737(cls):

        ccs = cls()

        substitutions = {
                         'C_{D_{wm}}': 0.5, # [2]
                         'C_{L_{vmax}}': 2.6, # [2]
                         'T_e': 1.29e5, # [4]
                         'V_1': 65,
                         'V_c': 234, # [7]
                         'V_{ne}': 144, # [2]
                         '\\mu': 1.4e-5, # [5]
                         '\\rho_c': 0.38, # [6]
                         '\\rho_{TO}': 1.225,
                         '\\tan(\\Lambda_{LE})': np.tan(40*np.pi/180),
                         'c_{l_{vt}}': 0.5, # [2]
                         'd_{fan}': 1.75, # [1]
                         'e': 0.8,
                         'l_{fuse}': 39,
                         'x_{CG}': 18,
                         'y_{eng}': 4.83, # [3]
                         }

        m = Model(ccs.cost, ccs, substitutions)
        return m

    @classmethod
    def aircraft_737(cls):

        ccs = cls()

        constraints = ccs + ccs.CG_constraint + ccs.linking_constraints

        substitutions = {
                         'C_{D_{wm}}': 0.5, # [2]
                         'C_{L_{vmax}}': 2.6, # [2]
                         'T_e': 1.29e5, # [4]
                         'V_1': 65,
                         'V_c': 234, # [7]
                         'V_{ne}': 144, # [2]
                         '\\mu': 1.4e-5, # [5]
                         '\\rho_c': 0.38, # [6]
                         '\\rho_{TO}': 1.225,
                         '\\tan(\\Lambda_{LE})': np.tan(40*np.pi/180),
                         'c_{l_{vt}}': 0.5, # [2]
                         'd_{fan}': 1.75, # [1]
                         'e': 0.8,
                         'y_{eng}': 4.83, # [3]
                         }

        m = Model(ccs.cost, constraints, substitutions, name='VerticalTail')
        return m

    @classmethod
    def test(cls):
        vt = cls.standalone_737()
        sol = vt.localsolve()

if __name__ == "__main__":
    VerticalTail.test()
