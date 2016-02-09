""" A simple aircraft performance model to be used as GPkit tutorial"""
from gpkit import Variable, Model
from gpkit.tools import te_exp_minus1
import numpy as np

# Free Variables
CD     = Variable('C_D', '-', 'Drag coefficient')
CL     = Variable('C_L', '-', 'Lift coefficient')
M      = Variable('M', '-', 'Mach number')
Mdd    = Variable('M_{dd}', '-', 'Drag divergence Mach number')
R      = Variable('R', 'nautical_miles', 'Range')
V      = Variable('V', 'knots', 'Velocity')
W_fuel = Variable('W_{fuel}', 'lbf', 'Fuel weight')
W_init = Variable('W_{init}', 'lbf', 'Initial gross weight')
z_bre  = Variable('z_{bre}', '-', 'Breguet parameter')

# Fixed Parameters
a     = Variable('a', 300, 'm/s', 'Speed of sound at 33,000 ft')
AR    = Variable('AR', 11 , '-', 'Aspect Ratio')
cosL  = Variable('\\cos(\\Lambda)', np.cos(25*np.pi/180), '-',
                 'Cosine of quarter chord sweep')
CD0   = Variable('CD0', 0.02, '-', 'Parasitic drag coefficient')
e     = Variable('e', 0.8, '-', 'Oswald efficiency factor')
k     = Variable('\\kappa', 0.95, '-', 'Technology factor')
MTOW  = Variable('MTOW', 120000, 'lbf', 'Max takeoff weight')
rho   = Variable('\\rho', 0.4, 'kg/m^3', 'Air density')
S     = Variable('S', 'm^2', 'Wing reference area')
tau   = Variable('\\tau', 0.13, '-', 'Wing thickness-to-chord ratio')
TSFC  = Variable('TSFC', 0.307, 'lb/lbf/hr',
                 'Thrust specific fuel consumption')
W_e   = Variable('W_e', 60000, 'lbf', 'Operating empty weight')
W_pay = Variable('W_{pay}', 40000, 'lbf', 'Payload weight')

# Constants
g     = Variable('g', 9.81, 'm/s^2', 'Gravitational acceleration')

# Model
objective = 1/R  # Maximize range

constraints = [# Aircraft and fuel weight
               W_init >= W_e + W_pay + W_fuel,
               
               # Performance constraints
               W_init <= MTOW,

               # Breguet range
               R <= z_bre*CL/CD*V/(TSFC*g),
               # Taylor series expansion of exp(z_bre) - 1
               W_fuel/W_e >= te_exp_minus1(z_bre, nterm=3),

               # Steady level flight
               W_init == 0.5*rho*V**2*S*CL,

               # Drag
               CD >= CD0 + CL**2/(np.pi*e*AR),

               # Speed limited by compressible drag rise
               M == V/a,
               M <= Mdd,
               # Korn equation
               Mdd + tau/cosL**2 + CL/(10*cosL**3) <= k/cosL, 

              ]

m = Model(objective, constraints)
m.solve()
