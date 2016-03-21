from gpkit import Variable, Model, SignomialsEnabled, units
from gpkit.constraints.tight import TightConstraintSet as TCS
from numpy import pi, tan

alpha = Variable('\\alpha', '-')
ARh  = Variable('AR_h', '-')
amax = Variable('\\alpha_{max}', '-')
eta  = Variable('\\eta', '-')
CLah = Variable('C_{L_{ah}}', '-')
CLaw = Variable('C_{L_{aw}}', '-')
CLh  = Variable('C_{L_h}', '-')
CLw  = Variable('C_{L_w}', '-')
cma  = Variable('\\bar{c}', 'm')
Cmac = Variable('C_{m_{ac}}', '-')
Kf   = Variable('K_f', '-')
lh   = Variable('l_h', 'm')
lfuse= Variable('l_{fuse}', 'm')
Sh   = Variable('S_h', 'm^2')
Sw   = Variable('S_w', 'm^2')
SM   = Variable('S.M.', '-')
SMmin= Variable('S.M._{min}', '-')
tanLh= Variable('\\tan(\\Lambda_h)', '-')
xCG  = Variable('x_{CG}', 'm')
xw   = Variable('x_w', 'm')
W    = Variable('W', 'N')
wf   = Variable('w_f', 'm')


objective = W

with SignomialsEnabled():
    constraints = [
                   # Stability
                   TCS([SM + xw/cma + Kf*wf**2*lfuse/(CLaw*Sw*cma)
                        <= CLah*Sh*lh/(CLaw*Sw*cma)]),
                   SM >= SMmin,

                   # Trim
                   CLh*Sh*lh/(Sw*cma) >= Cmac + CLw*xw/cma,

                   TCS([xCG + lh <= lfuse]),

                   # DATCOM formula (Mach number makes it SP)
                   TCS([(ARh/eta)**2 * (1 + tanLh**2) + 8*pi*ARh/CLah
                        <= (2*pi*ARh/CLah)**2]),
                   CLh == CLah*alpha,
                   alpha <= amax,

                   W >= 100*Sh*units.Pa
                  ]

substitutions = {
                 '\\alpha_{max}': 0.1, # (6 deg)
                 'AR_h': 4,
                 '\\eta': 0.97,
                 'C_{m_{ac}}': 0.1,
                 'C_{L_{aw}}': 2*pi,
                 'C_{L_w}': 0.5,
                 '\\bar{c}': 5,
                 'K_f': 0.7,
                 'l_{fuse}': 40,
                 'S_w': 125,
                 'S.M._{min}': 0.05,
                 '\\tan(\\Lambda_h)': tan(30*pi/180),
                 'x_{CG}': 20,
                 'x_w': 2,
                 'w_f': 6,
                }

m = Model(objective, constraints, substitutions)

m.solve()
