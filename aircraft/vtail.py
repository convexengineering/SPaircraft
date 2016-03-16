# coding=utf-8
from gpkit import Model, Variable, SignomialsEnabled, units, LinkConstraint
import numpy as np
import numpy.testing as npt
from wing.wingbox import WingBox
from gpkit.small_scripts import mag

class VerticalTail(Model):
    """
    Vertical tail sizing
    """
    def __init__(self, **kwargs):
        # Variables
        Avt    = Variable('A_{vt}', '-', 'Vertical tail aspect ratio')
        b      = Variable('b', 'm', 'Vertical tail full span')
        bvt    = Variable('b_{vt}', 'm', 'Vertical tail half span')
        CDvis  = Variable('C_{D_{vis}}', '-', 'Viscous drag coefficient')
        CLvt   = Variable('C_{L_{vt}}', '-', 'Vertical tail lift coefficient')
        cma    = Variable('\\bar{c}', 'm', 'Vertical tail mean aero chord')
        ctip   = Variable('c_{tip}', 'm', 'Vertical tail tip chord')
        croot  = Variable('c_{root}', 'm', 'Vertical tail root chord')
        Dvis   = Variable('D_{vis}', 'N', 'Vertical tail viscous drag, cruise')
        Dwm    = Variable('D_{wm}', 'N', 'Engine out windmill drag')
        dxlead = Variable('\\Delta x_{lead}', 'm',
                          'Distance from CG to vertical tail leading edge')
        dxtrail= Variable('\\Delta x_{trail}', 'm',
                          'Distance from CG to vertical tail trailing edge')
        lvt    = Variable('l_{vt}', 'm', 'Vertical tail moment arm')
        Lvt    = Variable('L_{vt}', 'N', 'Vertical tail lift in engine out')
        Lmax   = Variable('L_{max}','N', 'Maximum load for structural sizing')
        p      = Variable('p', '-', 'Substituted variable = 1 + 2*taper')
        q      = Variable('q', '-', 'Substituted variable = 1 + taper')
        Rec    = Variable('Re_c', '-', 'Vertical tail reynolds number, cruise')
        S      = Variable('S', 'm^2', 'Vertical tail reference area (full span)')
        Svt    = Variable('S_{vt}', 'm^2', 'Vertical tail ref. area (half span)')
        taper  = Variable('\\lambda', '-', 'Vertical tail taper ratio')
        tau    = Variable('\\tau', '-', 'Vertical tail thickness/chord ratio')
        Wvt    = Variable('W_{vt}', 'N', 'Vertical tail weight')
        Wstruct= Variable('W_{struct}', 'N', 'Full span weight')
        zmac   = Variable('z_{\\bar{c}}', 'm',
                          'Vertical location of mean aerodynamic chord')

        # Constants
        # References:
        # 1: LEAP engine specs (1B)
        # 2: TASOPT 737 code
        # 3: Boeing 737 airport doc
        # 4: Boeing 737 Max airport doc
        # 5: http://www.digitaldutch.com/atmoscalc
        # 6: Engineering toolbox
        # 7: Boeing.com
        Aeng  = Variable('A_{eng}', np.pi*(1.75/2)**2, 'm^2', 
                         'Engine reference area') # [1]
        CDwm  = Variable('C_{D_{wm}}', 0.5, '-', 'Windmill drag coefficient') # [2]
        CLvmax= Variable('C_{L_{vmax}}', 2.6, '-', 'Max lift coefficient')
        clvt  = Variable('c_{l_{vt}}', 0.5, '-',
                         'Sectional lift force coefficient (engine out)') # [2]
        e     = Variable('e', 0.8, '-', 'Span efficiency of vertical tail')
        le    = Variable('l_e', 4.83, 'm', 'Engine moment arm') # [3]
        Lfuse = Variable('L_{fuse}', 39, 'm', 'Length of fuselage') # [4]
        mu    = Variable('\\mu', 1.4e-5, 'N*s/m^2', 'Dynamic viscosity (35,000ft)') # [5]
        mu0   = Variable('\\mu_0', 1.8e-5, 'N*s/m^2', 'Dynamic viscosity (SL)') # [5]
        rho_c = Variable('\\rho_c', 0.38, 'kg/m^3', 'Air density (35,000ft)') # [6]
        rho0  = Variable('\\rho_{TO}', 1.225, 'kg/m^3', 'Air density (SL))')
        tanL  = Variable('\\tan(\\Lambda_{LE})', np.tan(40*np.pi/180), '-',
                         'Tangent of leading edge sweep (40 deg)')
        Te    = Variable('T_e', 1.29e5, 'N', 'Thrust per engine at takeoff') # [4]
        V1    = Variable('V_1', 65, 'm/s', 'Minimum takeoff velocity')
        Vc    = Variable('V_c', 234, 'm/s', 'Cruise velocity') # [7]
        Vne   = Variable('V_{ne}', 144, 'm/s', 'Never exceed velocity')
        xCG   = Variable('x_{CG}', 18, 'm', 'x-location of CG')


        with SignomialsEnabled():
            objective = Dvis + 0.5*Wvt
            constraints = [
                           Lvt*lvt >= Te*le + Dwm*le,
                           # Force moment balance for one engine out condition
                           # TASOPT 2.0 p45

                           dxlead + zmac*tanL + 0.25*cma >= lvt, # [SP]
                           # Tail moment arm

                           Lvt == 0.5*rho0*V1**2*Svt*CLvt,
                           # Vertical tail force (y-direction) for engine out

                           Avt == bvt**2/Svt,
                           CLvt*(1 + clvt/(np.pi*e*Avt)) <= clvt,
                           # Finite wing theory
                           # people.clarkson.edu/~pmarzocc/AE429/AE-429-4.pdf
                           # Valid because tail is untwisted and uncambered
                           # (lift curve slope passes through origin)

                           Dwm >= 0.5*rho0*V1**2*Aeng*CDwm,
                           # Drag of a windmilling engine

                           Svt <= bvt*(croot + ctip)/2, # [SP]
                           # Tail geometry relationship

                           dxtrail >= croot + dxlead,
                           # Tail geometry constraint

                           Lfuse >= dxtrail + xCG,
                           # Fuselage length constrains the tail trailing edge

                           p >= 1 + 2*taper,
                           2*q >= 1 + p,
                           zmac == (bvt/3)*q/p,
                           (2./3)*(1 + taper + taper**2)*croot/q >= cma, # [SP]
                           taper == ctip/croot,
                           # Define vertical tail geometry

                           taper >= 0.27,
                           # TODO: Constrain taper by tip Reynolds number
                           # source: b737.org.uk

                           Dvis >= 0.5*rho_c*Vc**2*Svt*CDvis,
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
                           # Relate vertical tail geometry/weight to generic
                           # wing used in structural model

                           Lmax == 0.5*rho0*Vne**2*Svt*CLvmax,
                           # Max load for structural sizing
                           ]

        # Incorporate the structural model
        wb = WingBox()
        wb.subinplace({"taper": taper})
        lc = LinkConstraint([constraints, wb])

        Model.__init__(self, objective, lc, **kwargs)

    def test(self):
        sol = self.localsolve()

        npt.assert_almost_equal(mag(sol('\\Delta x_{trail}')), mag(sol('c_{root}'))
                                + mag(sol('\\Delta x_{lead}')), decimal=4)
        npt.assert_almost_equal(mag(sol('\\Delta x_{lead}')) + mag(sol('z_{\\bar{c}}'))
                                *mag(sol('\\tan(\\Lambda_{LE})'))
                                + 0.25*mag(sol('\\bar{c}')), mag(sol('l_{vt}')),
                                decimal=3)
        npt.assert_almost_equal(mag(sol('L_{fuse}')), mag(sol('\\Delta x_{trail}'))
                                + mag(sol('x_{CG}')), decimal=4)
        npt.assert_almost_equal(mag(sol('C_{L_{vt}}'))*(1 + mag(sol('c_{l_{vt}}'))
                                /(np.pi*mag(sol('e'))*mag(sol('A_{vt}')))),
                                mag(sol('c_{l_{vt}}')), decimal=5)
        npt.assert_almost_equal((2./3)*(1 + mag(sol('\\lambda'))
                                + mag(sol('\\lambda'))**2)*mag(sol('c_{root}'))/mag(sol('q')),
                                mag(sol('\\bar{c}')))

if __name__ == "__main__":
    vt = VerticalTail()
    vt.test()
