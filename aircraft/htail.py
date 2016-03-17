from gpkit import Variable, Model, SignomialsEnabled, units
from numpy import pi


AR = Variable('AR')
CLa = Variable('C_{L_\\alpha}')
eta = Variable('\\eta')
Sh = Variable('Sh')
SM = Variable('S.M.')
dCG = Variable('\\Delta CG')
de_da = Variable('\\frac{de}{da}')
dcmf_dcl = Variable('\\frac{dCmf}{dCL}')


CLh  = Variable('C_{L_h}', '-')
CLw  = Variable('C_{L_w}', '-')
cma  = Variable('\\bar{c}', 'm')
Cmac = Variable('C_{m_{ac}}', '-')
lh   = Variable('l_h', 'm')
lfuse= Variable('l_{fuse}', 'm')
Sh   = Variable('S_h', 'm^2')
Sw   = Variable('S_w', 'm^2')
xw   = Variable('x_w', 'm')
W    = Variable('W', 'N')


objective = W

constraints = [
               # Stability


               # Trim
               CLh*Sh*lh/(Sw*cma) >= Cmac + CLw*xw/cma,

#               lh <= lfuse,

               # DATCOM formula (TODO: add sweep, Mach)
#               CLa * (2 + ((AR/eta)**2 + 4)**0.5) <= 2*pi*AR,

               W >= 100*Sh*units.Pa
              ]

substitutions = {
                 'C_{m_{ac}}': 0.1,
                 'C_{L_h}': 0.5,
                 'C_{L_w}': 0.5,
                 '\\bar{c}': 5,
#                 'l_{fuse}': 40,
                 'l_h': 20,
                 'S_w': 100,
                 'x_w': 2,
                }

m = Model(objective, constraints, substitutions)

m.solve()
