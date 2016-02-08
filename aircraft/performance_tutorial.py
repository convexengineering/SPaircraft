""" A simple aircraft performance model to be used as GPkit tutorial"""
from gpkit import Variable, Model
from gpkit.tools import te_exp_minus1

# Fixed Parameters
LD_max = Variable('\\left(\\frac{L}{D}\\right)_{max}', 15, '-',
                  'Maximum lift-to-drag ratio')
MTOW   = Variable('MTOW', 10000, 'lbf', 'Max takeoff weight')
TSFC   = Variable('TSFC', 0.307, 'lb/lbf/hr',
                  'Thrust specific fuel consumption')
V_max  = Variable('V_{max}', 420, 'knots', 'Maximum velocity')
W_e    = Variable('W_{e}', 7000, 'lbf', 'Operating empty weight')

# Constants
g      = Variable('g', 9.81, 'm/s^2', 'Gravitational acceleration')

# Free Variables
LD     = Variable('\\frac{L}{D}', '-', 'Lift-to-drag ratio')
R      = Variable('R', 'nautical_miles', 'Range')
V      = Variable('V', 'knots', 'Velocity')
W_fuel = Variable('W_{fuel}', 'lbf', 'Fuel weight')
W_init = Variable('W_{init}', 'lbf', 'Initial gross weight')
z_bre  = Variable('z_{bre}', '-', 'Breguet parameter')

# Model
objective = 1/R  # Maximize range

constraints = [# Aircraft and fuel weight
               W_init >= W_e + W_fuel,
               
               # Performance constraints
               W_init <= MTOW,
               LD <= LD_max,
               V <= V_max,

               # Breguet range
               R <= z_bre*LD*V/(TSFC*g),
               # Taylor series expansion of exp(z_bre) - 1
               W_fuel/W_e >= te_exp_minus1(z_bre, nterm=3)
              ]

m = Model(objective, constraints)
m.solve()
