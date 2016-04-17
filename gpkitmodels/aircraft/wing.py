from gpkit import Variable, Model, SignomialsEnabled, LinkedConstraintSet
from gpkit.constraints.costed import CostedConstraintSet
from gpkit.constraints.tight import TightConstraintSet as TCS
from numpy import pi, tan
from wingbox import WingBox
# pylint:disable=bad-whitespace

class Wing(CostedConstraintSet):
    """
    Wing sizing
    """
    def __init__(self, **kwargs):
        alpha   = Variable('\\alpha_w', '-', 'Wing angle of attack')
        AR      = Variable('AR', '-', 'Wing aspect ratio')
        amax    = Variable('\\alpha_{max,w}', '-', 'Max angle of attack')
        b       = Variable('b_w', 'm', 'Wing span')
        CDw     = Variable('C_{D_w}', '-', 'Drag coefficient')
        CD0w    = Variable('C_{D_{0_w}}', '-',
                           'Wing parasitic drag coefficient')
        CLaw    = Variable('C_{L_{aw}}', '-', 'Lift curve slope (wing)')
        CLw     = Variable('C_{L_w}', '-', 'Lift coefficient (wing)')
        CLwmax  = Variable('C_{L_{wmax}}', '-', 'Lift coefficient (wing)')
        croot   = Variable('c_{root}', 'm', 'Wing root chord')
        ctip    = Variable('c_{tip}', 'm', 'Wing tip chord')
        cwma    = Variable('\\bar{c}_{wing}', 'm', 'Mean aerodynamic chord (wing)')
        D       = Variable('D_{wing}', 'N', 'Wing drag')
        e       = Variable('e_w', '-', 'Oswald efficiency factor')
        eta     = Variable('\\eta_w', '-',
                           'Lift efficiency (diff between sectional and actual lift)')
        fl      = Variable('f(\\lambda_w)', '-', 'Empirical efficiency function of taper')
        Lmax    = Variable('L_{max_{w}}', 'N', 'Maximum load')
        mu      = Variable('\\mu', 'N*s/m^2', 'Dynamic viscosity (35,000ft)')
        p       = Variable('p_w', '-', 'Substituted variable = 1 + 2*taper')
        q       = Variable('q_w', '-', 'Substituted variable = 1 + taper')
        Rec     = Variable('Re_w', '-', 'Cruise Reynolds number (Wing)')
        rho     = Variable('\\rho', 'kg/m^3', 'Air density (35,000 ft)')
        rho0    = Variable('\\rho_0', 'kg/m^3', 'Air density (0 ft)')
        Sw      = Variable('S_w', 'm^2', 'Wing area')
        tanL    = Variable('\\tan(\\Lambda)', '-', 'tangent of wing sweep')
        taper   = Variable('\\lambda', '-', 'Wing taper ratio')
        tau     = Variable('\\tau_w', '-', 'Wing thickness/chord ratio')
        Vinf    = Variable('V_{\\infty}', 'm/s', 'Freestream velocity')
        Vne     = Variable('V_{ne}', 'm/s', 'Never exceed velocity')
        W       = Variable('W', 'N', 'Aircraft weight')
        W0      = Variable('W_0', 'N', 'Weight excluding wing')
        Ww      = Variable('W_{wing}', 'N', 'Wing weight')
        #xw      = Variable('x_w', 'm', 'Position of wing aerodynamic center')
        ymac    = Variable('y_{\\bar{c}}', 'm',
                           'Spanwise location of mean aerodynamic chord')

        objective = D + 0.5*Ww

        with SignomialsEnabled():
            constraints = [
                           W == 0.5*rho*Vinf**2*Sw*CLw,

        #                   p >= 1 + 2*taper,
        #                   2*q >= 1 + p,
        #                   ymac == (b/3)*q/p,
        #                   TCS([(2./3)*(1 + taper + taper**2)*croot/q >= cwma]), # [SP]
                           taper == ctip/croot,
                           TCS([Sw <= b*(croot + ctip)/2]), # [SP]

                           # DATCOM formula (Mach number makes it SP)
                           TCS([(AR/eta)**2 * (1 + tanL**2) + 8*pi*AR/CLaw
                                <= (2*pi*AR/CLaw)**2]),
                           CLw == CLaw*alpha,
                           alpha <= amax,

                           # Drag
                           D == 0.5*rho*Vinf**2*Sw*CDw,
                           CDw >= CD0w + CLw**2/(pi*e*AR),
                           #Rec == rho*Vinf*chma/mu,

                           # Oswald efficiency
                           # Nita, Scholz, "Estimating the Oswald factor from basic
                           # aircraft geometrical parameters"
                           TCS([fl >= 0.0524*taper**4 - 0.15*taper**3 + 0.1659*taper**2
                                    - 0.0706*taper + 0.0119], reltol=1E-2),
                           TCS([e*(1 + fl*AR) <= 1]),
                           taper >= 0.2, # TODO

                           Lmax == 0.5*rho0*Vne**2*Sw*CLwmax,
                          ]

            standalone_constraints = [W >= W0 + Ww,]
            self.standalone_constraints = standalone_constraints

        wb = WingBox()
        wb.subinplace({'A': AR,
                       'b': b,
                       'L_{max}': Lmax,
                       'p': p,
                       'q': q,
                       'S': Sw,
                       'taper': taper,
                       '\\tau': tau,
                       'W_{struct}': Ww})

        lc = LinkedConstraintSet([constraints, wb])

        CostedConstraintSet.__init__(self, objective, lc, **kwargs)

    @classmethod
    def standalone_737(cls):
        ccs = cls()

        constraints = ccs + ccs.standalone_constraints

        substitutions = {
                         '\\alpha_{max,w}': 0.1, # (6 deg)
                         '\\eta_w': 0.97,
                         'C_{D_{0_w}}': 0.05,
                         'C_{L_{wmax}}': 2.5,
                         '\\mu': 1.4E-5,
                         '\\rho': 0.38,
                         '\\rho_0': 1.225,
                         '\\tan(\\Lambda)': tan(30*pi/180),
                         'V_{\\infty}': 240,
                         'V_{ne}': 144,
                         'W_0': 1E5,
                        }

        m = Model(ccs.cost, constraints, substitutions)
        return m

    @classmethod
    def aircraft_737(cls):
        ccs = cls()

        substitutions = {
                         '\\alpha_{max,w}': 0.1, # (6 deg)
                         '\\eta_w': 0.97,
                         'C_{D_{0_w}}': 0.05,
                         'C_{L_{wmax}}': 2.5,
                         '\\mu': 1.4E-5,
                         '\\rho': 0.38,
                         '\\rho_0': 1.225,
                         '\\tan(\\Lambda)': tan(30*pi/180),
                         'V_{ne}': 144,
                        }

        m = Model(ccs.cost, ccs, substitutions, name='Wing')
        return m

    @classmethod
    def test(cls):
        w = cls.standalone_737()
        sol = w.localsolve()

if __name__ == "__main__":
    Wing.test()
