from gpkit import Variable, Model, SignomialsEnabled, LinkConstraint, units
from gpkit.constraints.tight import TightConstraintSet as TCS
from numpy import pi, tan
from wing.wingbox import WingBox

alpha   = Variable('\\alpha', '-', 'Horizontal tail angle of attack')
ARh     = Variable('AR_h', '-', 'Horizontal tail aspect ratio')
amax    = Variable('\\alpha_{max}', '-', 'Max angle of attack (htail)')
bht     = Variable('b_{ht}', 'm', 'Horizontal tail span')
CDh     = Variable('C_{D_h}', '-', 'Horizontal tail drag coefficient')
CD0h    = Variable('C_{D_{0_h}}', '-',
                   'Horizontal tail parasitic drag coefficient')
CLah    = Variable('C_{L_{ah}}', '-', 'Lift curve slope (htail)')
CLaw    = Variable('C_{L_{aw}}', '-', 'Lift curve slope (wing)')
CLh     = Variable('C_{L_h}', '-', 'Lift coefficient (htail)')
CLw     = Variable('C_{L_w}', '-', 'Lift coefficient (wing)')
Cmac    = Variable('|C_{m_{ac}}|', '-', # Absolute value of CMwing
                   'Moment coefficient about aerodynamic centre (wing)')
chma    = Variable('\\bar{c}_{ht}', 'm', 'Mean aerodynamic chord (ht)')
croot   = Variable('c_{root}', 'm', 'Horizontal tail root chord')
ctip    = Variable('c_{tip}', 'm', 'Horizontal tail tip chord')
cwma    = Variable('\\bar{c}_{wing}', 'm', 'Mean aerodynamic chord (wing)')
Cmfu    = Variable('C_{m_{fuse}}', '-', 'Moment coefficient (fuselage)')
D       = Variable('D_h', 'N', 'Horizontal tail drag')
dxlead  = Variable('\\Delta x_{lead}', 'm',
                   'Distance from CG to horizontal tail leading edge')
dxtrail = Variable('\\Delta x_{trail}', 'm',
                   'Distance from CG to horizontal tail trailing edge')
dxw     = Variable('\\Delta x_w', 'm', 'Distance from aerodynamic centre to CG')
e       = Variable('e_h', '-', 'Oswald efficiency factor')
eta     = Variable('\\eta', '-',
                   'Lift efficiency (diff between sectional and actual lift)')
fl      = Variable('f\(\\lambda\)', '-', 'Empirical efficiency function of taper')
Kf      = Variable('K_f', '-', 'Empirical factor for fuselage-wing interference')
lh      = Variable('l_h', 'm', 'Horizontal tail moment arm')
lfuse   = Variable('l_{fuse}', 'm', 'Fuselage length')
Lmax    = Variable('L_{max}', 'N', 'Maximum load')
mu      = Variable('\\mu', 'N*s/m^2', 'Dynamic viscosity (35,000ft)')
p       = Variable('p', '-', 'Substituted variable = 1 + 2*taper')
q       = Variable('q', '-', 'Substituted variable = 1 + taper')
Rec     = Variable('Re_c', '-', 'Cruise Reynolds number (Horizontal tail)')
rho     = Variable('\\rho', 'kg/m^3', 'Air density (35,000 ft)')
Sh      = Variable('S_h', 'm^2', 'Horizontal tail area')
Sw      = Variable('S_w', 'm^2', 'Wing area')
SM      = Variable('S.M.', '-', 'Stability margin')
SMmin   = Variable('S.M._{min}', '-', 'Minimum stability margin')
tanLh   = Variable('\\tan(\\Lambda_h)', '-', 'tangent of horizontal tail sweep')
taper   = Variable('\\lambda', '-', 'Horizontal tail taper ratio')
tau     = Variable('\\tau', '-', 'Horizontal tail thickness/chord ratio')
Vinf    = Variable('V_{\\infty}', 'm/s', 'Freestream velocity')
W       = Variable('W', 'N', 'Horizontal tail weight')
wf      = Variable('w_f', 'm', 'Fuselage width')
xCG     = Variable('x_{CG}', 'm', 'CG location')
xw      = Variable('x_w', 'm', 'Position of wing aerodynamic center')
ymac    = Variable('y_{\\bar{c}}', 'm',
                   'Vertical location of mean aerodynamic chord')

objective = D + 0.5*W

with SignomialsEnabled():
    constraints = [
                   # Stability
                   TCS([SM + dxw/cwma + Kf*wf**2*lfuse/(CLaw*Sw*cwma)
                        <= CLah*Sh*lh/(CLaw*Sw*cwma)]),
                   SM >= SMmin,

                   # Trim
                   TCS([CLh*Sh*lh/(Sw*cwma) + Cmac >= CLw*dxw/cwma + Cmfu]),

                   # Moment arm and geometry
                   TCS([dxlead + croot <= dxtrail]),
                   TCS([xCG + dxtrail <= lfuse]),
                   TCS([dxlead + ymac*tanLh + 0.25*chma >= lh], reltol=1e-5), # [SP]
                   p >= 1 + 2*taper,
                   2*q >= 1 + p,
                   ymac == (bht/3)*q/p,
                   TCS([(2./3)*(1 + taper + taper**2)*croot/q >= chma]), # [SP]
                   taper == ctip/croot,
                   TCS([Sh <= bht*(croot + ctip)/2]), # [SP]

                   # DATCOM formula (Mach number makes it SP)
                   TCS([(ARh/eta)**2 * (1 + tanLh**2) + 8*pi*ARh/CLah
                        <= (2*pi*ARh/CLah)**2]),
                   CLh == CLah*alpha,
                   alpha <= amax,

                   # K_f as function of wing position (fitted posynomial)
                   Kf >= 1.5012*(xw/lfuse)**2 + 0.538*(xw/lfuse) + 0.0331,
                   TCS([xw >= xCG + dxw]),

                   # Drag
                   D == 0.5*rho*Vinf**2*Sh*CDh,
                   CDh >= CD0h + CLh**2/(pi*e*ARh),
                   CD0h**0.125 >= 0.19*(tau)**0.0075 *(Rec)**0.0017
                                + 1.83e+04*(tau)**3.54*(Rec)**-0.494
                                + 0.118*(tau)**0.0082 *(Rec)**0.00165
                                + 0.198*(tau)**0.00774*(Rec)**0.00168,
                   Rec == rho*Vinf*chma/mu,

                   # Oswald efficiency
                   # Nita, Scholz, "Estimating the Oswald factor from basic
                   # aircraft geometrical parameters"
                   TCS([fl >= 0.0524*taper**4 - 0.15*taper**3 + 0.1659*taper**2
                            - 0.0706*taper + 0.0119]),
                   TCS([e*(1 + fl*ARh) <= 1]),
                   taper >= 0.2, # TODO
                  ]

# References
# [1] TASOPT code

substitutions = {
                 '\\alpha_{max}': 0.1, # (6 deg)
                 '\\eta': 0.97,
                 'C_{L_{aw}}': 2*pi,
                 'C_{L_w}': 0.5,
                 '|C_{m_{ac}}|': 0.1,
                 'C_{m_{fuse}}': 0.05, # [1]
                 '\\bar{c}_{wing}': 5,
                 'l_{fuse}': 40,
                 'L_{max}': 1E6, # TODO
                 '\\mu': 1.4E-5,
                 '\\rho': 0.38,
                 'S_w': 125,
                 'S.M._{min}': 0.05,
                 '\\tan(\\Lambda_h)': tan(30*pi/180),
                 'V_{\\infty}': 240,
                 'x_{CG}': 20,
                 '\\Delta x_w': 2,
                 'w_f': 6,
                }

wb = WingBox()
wb.subinplace({'A': ARh,
               'b': bht,
               'p': p,
               'q': q,
               'S': Sh,
               'taper': taper,
               '\\tau': tau,
               'W_{struct}': W})

lc = LinkConstraint([constraints, wb])

m = Model(objective, lc, substitutions)

m.localsolve()
