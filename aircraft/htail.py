from gpkit import Variable, Model, SignomialsEnabled, units
from gpkit.constraints.tight import TightConstraintSet as TCS
from numpy import pi, tan

alpha = Variable('\\alpha', '-', 'Horizontal tail angle of attack')
ARh  = Variable('AR_h', '-', 'Horizontal tail aspect ratio')
amax = Variable('\\alpha_{max}', '-', 'Max angle of attack (htail)')
CDh  = Variable('C_{D_h}', '-', 'Horizontal tail drag coefficient')
CD0h = Variable('C_{D_{0_h}}', '-',
                'Horizontal tail parasitic drag coefficient')
CLah = Variable('C_{L_{ah}}', '-', 'Lift curve slope (htail)')
CLaw = Variable('C_{L_{aw}}', '-', 'Lift curve slope (wing)')
CLh  = Variable('C_{L_h}', '-', 'Lift coefficient (htail)')
CLw  = Variable('C_{L_w}', '-', 'Lift coefficient (wing)')
cma  = Variable('\\bar{c}', 'm', 'Mean aerodynamic chord (wing)')
Cmac = Variable('C_{m_{ac}}', '-',
                'Moment coefficient about aerodynamic centre (wing)')
Cmfu = Variable('C_{m_{fuse}}', '-', 'Moment coefficient (fuselage)')
D    = Variable('D_h', 'N', 'Horizontal tail drag')
e    = Variable('e_h', '-', 'Oswald efficiency factor')
eta  = Variable('\\eta', '-',
                'Lift efficiency (diff between sectional and actual lift)')
Kf   = Variable('K_f', '-', 'Empirical factor for fuselage-wing interference')
lh   = Variable('l_h', 'm', 'Horizontal tail moment arm')
lfuse= Variable('l_{fuse}', 'm', 'Fuselage length')
rho  = Variable('\\rho', 'kg/m^3', 'Air density (35,000 ft)')
Sh   = Variable('S_h', 'm^2', 'Horizontal tail area')
Sw   = Variable('S_w', 'm^2', 'Wing area')
SM   = Variable('S.M.', '-', 'Stability margin')
SMmin= Variable('S.M._{min}', '-', 'Minimum stability margin')
tanLh= Variable('\\tan(\\Lambda_h)', '-', 'tangent of horizontal tail sweep')
Vinf = Variable('V_{\\infty}', 'm/s', 'Freestream velocity')
W    = Variable('W', 'N', 'Horizontal tail weight')
wf   = Variable('w_f', 'm', 'Fuselage width')
xCG  = Variable('x_{CG}', 'm', 'CG location')
xw   = Variable('x_w', 'm', 'Distance from aerodynamic centre to CG')


objective = D + 0.5*W

with SignomialsEnabled():
    constraints = [
                   # Stability
                   TCS([SM + xw/cma + Kf*wf**2*lfuse/(CLaw*Sw*cma)
                        <= CLah*Sh*lh/(CLaw*Sw*cma)]),
                   SM >= SMmin,

                   # Trim
                   CLh*Sh*lh/(Sw*cma) >= Cmac + CLw*xw/cma + Cmfu,

                   TCS([xCG + lh <= lfuse]),

                   # DATCOM formula (Mach number makes it SP)
                   TCS([(ARh/eta)**2 * (1 + tanLh**2) + 8*pi*ARh/CLah
                        <= (2*pi*ARh/CLah)**2]),
                   CLh == CLah*alpha,
                   alpha <= amax,

                   W >= 100*Sh*units.Pa,

                   D == 0.5*rho*Vinf**2*Sh*CDh,
                   CDh >= CD0h + CLh**2/(pi*e*ARh)
                  ]

substitutions = {
                 '\\alpha_{max}': 0.1, # (6 deg)
                 'AR_h': 4,
                 '\\eta': 0.97,
                 'C_{D_{0_h}}': 0.03,
                 'C_{L_{aw}}': 2*pi,
                 'C_{L_w}': 0.5,
                 'C_{m_{ac}}': 0.1,
                 'C_{m_{fuse}}': 0.1,
                 'e_h': 0.8,
                 '\\bar{c}': 5,
                 'K_f': 0.7,
                 'l_{fuse}': 40,
                 '\\rho': 0.38,
                 'S_w': 125,
                 'S.M._{min}': 0.05,
                 '\\tan(\\Lambda_h)': tan(30*pi/180),
                 'V_{\\infty}': 240,
                 'x_{CG}': 20,
                 'x_w': 2,
                 'w_f': 6,
                }

m = Model(objective, constraints, substitutions)

m.solve()
