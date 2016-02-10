""" A simple aircraft performance model to be used as GPkit tutorial"""
from gpkit import Variable, Model
from gpkit.tools import te_exp_minus1
import numpy as np

# Free Variables
b      = Variable('b', 'm', 'Wing span')
CD     = Variable('C_D', '-', 'Drag coefficient')
CL     = Variable('C_L', '-', 'Lift coefficient')
LD_la  = Variable('\\frac{L}{D}_{land}', '-', 'Lift-to-drag in final descent')
M      = Variable('M', '-', 'Mach number')
Mdd    = Variable('M_{dd}', '-', 'Drag divergence Mach number')
R      = Variable('R', 'nautical_miles', 'Range')
S      = Variable('S', 'm^2', 'Wing reference area')
s_l    = Variable('s_L', 'm', 'Total landing distance')
s_l_gr = Variable('s_{L_{gr}}', 'm', 'Ground roll distance (landing)')
s_l_oc = Variable('s_{L_{oc}}', 'm' , 'Obstacle clearance distance (landing)')
s_t    = Variable('s_{TO}', 'm', 'Total takeoff distance')
s_t_gr = Variable('s_{TO_{gr}}', 'm', 'Ground roll distance (takeoff)')
T      = Variable('T', 'N', 'Thrust')
T_t    = Variable('T_{TO}', 'N', 'Takeoff thrust')
V      = Variable('V', 'knots', 'Velocity')
Vst    = Variable('V_{stall}', 'm/s', 'Stall speed')
W_fuel = Variable('W_{fuel}', 'lbf', 'Fuel weight')
W      = Variable('W', 'lbf', 'Initial gross weight')
z_bre  = Variable('z_{bre}', '-', 'Breguet parameter')

# Fixed Parameters
a     = Variable('a', 300, 'm/s', 'Speed of sound at 33,000 ft')
AR    = Variable('AR', 11 , '-', 'Aspect Ratio')
bmax  = Variable('b_{max}', 40, 'm', 'Gate constraint')
cosL  = Variable('\\cos(\\Lambda)', np.cos(25*np.pi/180), '-',
                 'Cosine of quarter chord sweep')
CD0   = Variable('CD0', 0.02, '-', 'Parasitic drag coefficient')
CDl   = Variable('C_D_{land}', 0.08, '-',
                 'Drag coefficient in landing configuration')
CLmaxt= Variable('C_{L_{max}}', 2.0, '-', 'Maximum lift coefficient (takeoff)')
CLmaxl= Variable('C_{L_{max}}', 3.0, '-', 'Maximum lift coefficient (landing)')
e     = Variable('e', 0.8, '-', 'Oswald efficiency factor')
h_oc  = Variable('h_{oc}', 50, 'ft', 'Obstacle clearance height (landing)')
k     = Variable('\\kappa', 0.95, '-', 'Technology factor')
l_rw  = Variable('l_{rw}', 5000, 'ft', 'Runway Length (DCA)')
mu_r  = Variable('\\mu_r', 0.4, '-', 'Rolling friction coefficient')
rho   = Variable('\\rho', 0.4, 'kg/m^3', 'Air density')
s_t_oc= Variable('s_{TO_{oc}}', 200, 'm',
                  'Obstacle clearance distance (takeoff)')
tau   = Variable('\\tau', 0.13, '-', 'Wing thickness-to-chord ratio')
T_max = Variable('T_{max}', 1E6, 'N', 'Maximum thrust')
TSFC  = Variable('TSFC', 0.307, 'lb/lbf/hr',
                 'Thrust specific fuel consumption')
Vstmax= Variable('V_{stall_{max}}', 55, 'm/s', 'Maximum allowable stall speed')
W_e   = Variable('W_e', 60000, 'lbf', 'Operating empty weight')
W_pay = Variable('W_{pay}', 40000, 'lbf', 'Payload weight')

# Constants
g      = Variable('g', 9.81, 'm/s^2', 'Gravitational acceleration')
rho_sl = Variable('\\rho_{SL}', 1.225, 'kg/m^3', 'Air density at sea level')

# Model
objective = 1/R  # Maximize range

constraints = [# Weight buildup
               W >= W_e + W_pay + W_fuel,

               # Breguet range
               R <= z_bre*CL/CD*V/(TSFC*g),
               # Taylor series expansion of exp(z_bre) - 1
               W_fuel/W_e >= te_exp_minus1(z_bre, nterm=3),

               # Wing geometry
               AR == b**2/S,
               b <= bmax,

               # Steady level flight
               W == 0.5*rho*V**2*S*CL,
               T == 0.5*rho*V**2*S*CD,
               T <= T_max,

               # Drag buildup
               CD >= CD0 + CL**2/(np.pi*e*AR),

               # Speed limited by compressible drag rise
               M == V/a,
               M <= Mdd,
               # Korn equation
               Mdd + tau/cosL**2 + CL/(10*cosL**3) <= k/cosL, 

               # Takeoff
               s_t_gr >= 1.44*W**2/(g*rho_sl*S*CLmaxt*T_t),
               s_t >= s_t_gr +s_t_oc,
               s_t <= l_rw,
               T_t <= T_max,

               # Landing
               W <= 0.5*rho_sl*Vst**2*S*CLmaxl,
               Vst <= Vstmax,
               s_l_oc >= h_oc*LD_la, # assumes half fuel weight
               LD_la >= (W_e + 0.5*W_fuel)/(0.5*rho*(1.3*Vst)**2*S*CDl),
               s_l_gr == 1.69*W**2/(g*rho_sl*S*CLmaxl*mu_r*W), # very conservative
               s_l >= s_l_oc + s_l_gr,
               1.67*s_l <= l_rw, # FAR 25.125

              ]

m = Model(objective, constraints)
m.solve()
