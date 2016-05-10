"Implement HorizontalTail model"
from numpy import pi, tan
from gpkit import Variable, Model, SignomialsEnabled, LinkedConstraintSet
from gpkit.constraints.costed import CostedConstraintSet
from gpkit.constraints.tight import TightConstraintSet as TCS
from wingbox import WingBox
# pylint:disable=bad-whitespace

class HorizontalTail(CostedConstraintSet):
    """
    Horizontal tail sizing

    References:
    [1] TASOPT code
    [2] http://adg.stanford.edu/aa241/stability/staticstability.html
    """
    def __init__(self, **kwargs):
        ARh     = Variable('AR_h', '-', 'Horizontal tail aspect ratio')
        ARw     = Variable('AR', '-', 'Wing aspect ratio')
        CD0h    = Variable('C_{D_{0_h}}', '-',
                           'Horizontal tail parasitic drag coefficient')
        CDh     = Variable('C_{D_h}', '-', 'Horizontal tail drag coefficient')
        CLah    = Variable('C_{L_{ah}}', '-', 'Lift curve slope (htail)')
        CLah0   = Variable('C_{L_{ah_0}}', '-',
                           'Isolated lift curve slope (htail)')
        CLaw    = Variable('C_{L_{aw}}', '-', 'Lift curve slope (wing)')
        CLh     = Variable('C_{L_h}', '-', 'Lift coefficient (htail)')
        CLhmax  = Variable('C_{L_{hmax}}', '-', 'Max lift coefficient')
        CLw     = Variable('C_{L_w}', '-', 'Lift coefficient (wing)')
        Cmac    = Variable('|C_{m_{ac}}|', '-', # Absolute value of CMwing
                           'Moment coefficient about aerodynamic centre (wing)')
        Cmfu    = Variable('C_{m_{fuse}}', '-', 'Moment coefficient (fuselage)')
        D       = Variable('D_{ht}', 'N', 'Horizontal tail drag')
        Kf      = Variable('K_f', '-',
                           'Empirical factor for fuselage-wing interference')
        Lh      = Variable('L_h', 'N', 'Horizontal tail downforce')
        Lmax    = Variable('L_{{max}_h}', 'N', 'Maximum load')
        M       = Variable('M', '-', 'Mach number')
        Rec     = Variable('Re_{c_h}', '-',
                           'Cruise Reynolds number (Horizontal tail)')
        SM      = Variable('S.M.', '-', 'Stability margin')
        SMmin   = Variable('S.M._{min}', '-', 'Minimum stability margin')
        Sh      = Variable('S_h', 'm^2', 'Horizontal tail area')
        Sw      = Variable('S_w', 'm^2', 'Wing area')
        Vinf    = Variable('V_{\\infty}', 'm/s', 'Freestream velocity')
        Vne     = Variable('V_{ne}', 'm/s', 'Never exceed velocity')
        W       = Variable('W_{ht}', 'N', 'Horizontal tail weight')
        alpha   = Variable('\\alpha', '-', 'Horizontal tail angle of attack')
        amax    = Variable('\\alpha_{max,h}', '-', 'Max angle of attack, htail')
        bht     = Variable('b_{ht}', 'm', 'Horizontal tail span')
        chma    = Variable('\\bar{c}_{ht}', 'm', 'Mean aerodynamic chord (ht)')
        croot   = Variable('c_{root_h}', 'm', 'Horizontal tail root chord')
        ctip    = Variable('c_{tip_h}', 'm', 'Horizontal tail tip chord')
        cwma    = Variable('\\bar{c}_{wing}', 'm',
                           'Mean aerodynamic chord (wing)')
        dxlead  = Variable('\\Delta x_{{lead}_h}', 'm',
                           'Distance from CG to horizontal tail leading edge')
        dxtrail = Variable('\\Delta x_{{trail}_h}', 'm',
                           'Distance from CG to horizontal tail trailing edge')
        dxw     = Variable('\\Delta x_w', 'm',
                           'Distance from aerodynamic centre to CG')
        e       = Variable('e_h', '-', 'Oswald efficiency factor')
        eta     = Variable('\\eta_h', '-',
                           ("Lift efficiency (diff between sectional and "
                            "actual lift)"))
        etaht   = Variable('\\eta_{ht}', '-', 'Tail efficiency')
        fl      = Variable(r"f(\lambda_h)", '-',
                           'Empirical efficiency function of taper')
        lfuse   = Variable('l_{fuse}', 'm', 'Fuselage length')
        lh      = Variable('l_h', 'm', 'Horizontal tail moment arm')
        mu      = Variable(r'\mu', 'N*s/m^2', 'Dynamic viscosity (35,000ft)')
        p       = Variable('p_{ht}', '-', 'Substituted variable = 1 + 2*taper')
        q       = Variable('q_{ht}', '-', 'Substituted variable = 1 + taper')
        rho     = Variable(r'\rho', 'kg/m^3', 'Air density (35,000 ft)')
        rho0    = Variable(r'\rho_0', 'kg/m^3', 'Air density (0 ft)')
        tanLh   = Variable(r'\tan(\Lambda_h)', '-',
                           'tangent of horizontal tail sweep')
        taper   = Variable(r'\lambda_h', '-', 'Horizontal tail taper ratio')
        tau     = Variable(r'\tau_h', '-',
                           'Horizontal tail thickness/chord ratio')
        wf      = Variable('w_{fuse}', 'm', 'Fuselage width')
        xcg     = Variable('x_{CG}', 'm', 'CG location')
        xcght   = Variable('x_{CG_{ht}}', 'm', 'Horizontal tail CG location')
        xw      = Variable('x_w', 'm', 'Position of wing aerodynamic center')
        ymac    = Variable('y_{\\bar{c}_{ht}}', 'm',
                           'Vertical location of mean aerodynamic chord')

        objective = D + 0.5*W

        with SignomialsEnabled():
                           # Stability from UMich AE-481 course notes
            constraints = [TCS([SM + dxw/cwma + Kf*wf**2*lfuse/(CLaw*Sw*cwma)
                                <= CLah*Sh*lh/(CLaw*Sw*cwma)]),
                           SM >= SMmin,

                           # Trim from UMich AE-481 course notes
                           TCS([CLh*Sh*lh/(Sw*cwma) + Cmac >=
                                CLw*dxw/cwma + Cmfu], reltol=0.02),
                           Lh == 0.5*rho*Vinf**2*Sh*CLh,

                           # Moment arm and geometry -- same as for vtail
                           TCS([dxlead + croot <= dxtrail]),
                           TCS([xcg + dxtrail <= lfuse], reltol=0.002),
                           TCS([dxlead + ymac*tanLh + 0.25*chma >= lh],
                               reltol=1e-2), # [SP]
                           p >= 1 + 2*taper,
                           2*q >= 1 + p,
                           ymac == (bht/3)*q/p,
                           TCS([(2./3)*(1 + taper + taper**2)*croot/q >=
                                chma]), # [SP]
                           taper == ctip/croot,
                           TCS([Sh <= bht*(croot + ctip)/2]), # [SP]

                           # DATCOM formula (Mach number makes it SP)
                           TCS([(ARh/eta)**2*(1+tanLh**2-M**2) + 8*pi*ARh/CLah0
                                <= (2*pi*ARh/CLah0)**2]),
                           # Loss of tail effectiveness due to wing downwash
                           CLah + (2*CLaw/(pi*ARw))*etaht*CLah0 <= CLah0*etaht,
                           CLh == CLah*alpha,
                           alpha <= amax,

                           # K_f as f(wing position) -- (fitted posynomial)
                           # from from UMich AE-481 course notes Table 9.1
                           Kf >= (1.5012*(xw/lfuse)**2 +
                                  0.538*(xw/lfuse) +
                                  0.0331),
                           TCS([xw >= xcg + dxw]),

                           # Drag
                           D == 0.5*rho*Vinf**2*Sh*CDh,
                           CDh >= CD0h + CLh**2/(pi*e*ARh),
                           # same drag model as vtail
                           CD0h**0.125 >= 0.19*(tau)**0.0075 *(Rec)**0.0017
                                        + 1.83e+04*(tau)**3.54*(Rec)**-0.494
                                        + 0.118*(tau)**0.0082 *(Rec)**0.00165
                                        + 0.198*(tau)**0.00774*(Rec)**0.00168,
                           Rec == rho*Vinf*chma/mu,

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

                           Lmax == 0.5*rho0*Vne**2*Sh*CLhmax,

                           # TODO: add correction due to wing downwash (kroo)
                          ]

        CG_constraint = [TCS([xcght >= xcg+(dxlead+dxtrail)/2]),
                         xcght <= lfuse
                        ]
        self.CG_constraint = CG_constraint

        wb = WingBox()
        wb.subinplace({'A': ARh,
                       'b': bht,
                       'L_{max}': Lmax,
                       'p': p,
                       'q': q,
                       'S': Sh,
                       'taper': taper,
                       r'\tau': tau,
                       'W_{struct}': W})

        lc = LinkedConstraintSet([constraints, wb])

        CostedConstraintSet.__init__(self, objective, lc, **kwargs)

    def default737subs(self):

        substitutions = {
                         'AR': 9,
                         'C_{L_w}': 0.5,
                         'C_{L_{aw}}': 2*pi,
                         'C_{L_{hmax}}': 2.6,
                         'C_{m_{fuse}}': 0.05, # [1]
                         'M': 0.8,
                         'S.M._{min}': 0.05,
                         'S_w': 125,
                         'V_{\\infty}': 240,
                         'V_{ne}': 144,
                         '\\Delta x_w': 2,
                         '\\alpha_{max,h}': 0.1, # (6 deg)
                         '\\bar{c}_{wing}': 5,
                         '\\eta_h': 0.97, # [2]
                         '\\eta_{ht}': 0.9, # [2]
                         '\\mu': 1.4E-5,
                         '\\rho': 0.38,
                         '\\rho_0': 1.225,
                         '\\tan(\\Lambda_h)': tan(30*pi/180),
                         'l_{fuse}': 40,
                         'w_{fuse}': 6,
                         'x_{CG}': 20,
                         '|C_{m_{ac}}|': 0.1,
                        }
        return substitutions

    @classmethod
    def standalone737(cls):

        ccs = cls()

        substitutions = ccs.default737subs()

        m = Model(ccs.cost, ccs, substitutions)
        return m

    @classmethod
    def coupled737(cls):

        ccs = cls()

        constraints = ccs + ccs.CG_constraint

        dsubs = ccs.default737subs()
        linkedsubs = ['AR', 'C_{L_w}', 'C_{L_{aw}}', 'S_w', 'V_{\\infty}',
                      '\\bar{c}_{wing}', 'l_{fuse}', 'w_{fuse}', 'x_{CG}']
        substitutions = {key: value for key, value in dsubs.items()
                                    if key not in linkedsubs}

        m = Model(ccs.cost, constraints, substitutions, name='HorizontalTail')
        return m

    @classmethod
    def test(cls):
        ht = cls.standalone737()
        sol = ht.localsolve()

if __name__ == "__main__":
    HorizontalTail.test()
