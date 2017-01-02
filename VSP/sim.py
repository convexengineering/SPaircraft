from gpkit import Model, Variable
import numpy as np
import gpkit

def updateOpenVSP(inputDict):
	filename = 'design.des'
	with open(filename,'r+') as f:
		result = f.read()
		a = result.split('\n')
		outputLines = []
		for line in a:
			words = line.split(':')
			if len(words) > 1:
				key = words[0]
				value = float(words[-1])
				if key in inputDict:
					value = " " + str(inputDict[key])
				words[-1] = value
			outputLine = ":".join(words)
			outputLines += [outputLine]
		output = '\n'.join(outputLines)
		print('OpenVSP .des output:')
		print(output)
		f.seek(0)
		f.write(output)
		f.truncate()
		f.close()

z = Variable('z','-','Breguet range parameter')
theta_fuel = Variable('theta_fuel','-','Fuel fraction W_fuel/W_zfw')
g = Variable('g',3.711,'m/s/s','Surface gravity on Mars')
rho = Variable('rho',0.02,'kg/m^3','Surface air density')
mu = Variable("\\mu", 1.08e-5, "kg/m/s", "viscosity of air")

V = Variable('V','m/s','Cruise velocity')
Vmin = Variable('Vmin','m/s','Landing velocity')
R = Variable('R','m','Maximum theoretical range')
L = Variable('L','N','Total lift in cruise')
D = Variable('D','N','Total drag in cruise')
h_fuel =  Variable('h_fuel',5.55e7,'J/kg','Fuel heating value') #Nitromethane (just the methane)
n_0 = Variable('n_0',0.3,'-','Whole chain propulsion efficiency')

# Masses
m_payload = Variable('m_payload',1,'kg','payload (max requirement)')
m_fuel = Variable('m_fuel','kg','full fuel mass')
rho_fuel = Variable('rho_fuel',1137.1,'kg/m^3','density of fuel (nitromethane)')
m_limit = Variable('m_limit',5,'kg','max mass')
V_fuel = Variable('V_fuel','m^3','Volume of fuel tank')

# Weights
W_i = Variable('W_i','N','Total weight on mission takeoff')
W_f = Variable('W_f','N','Final weight, mission end')



V_wind = Variable('V_wind',20,'m/s','Worst case wind')

constraints = [
	theta_fuel == m_fuel*g/W_f,
	V_fuel == m_fuel/rho_fuel,
	z >= (g*R*D)/(h_fuel*n_0*L),
	theta_fuel >= z + (z**2)/2 + (z**3)/6 + (z**4)/24,
	L >= W_i,
	# L/D <= 17,
	# m_fuel >= theta_fuel*W_f/g,
]

# Wing
W_W_coeff1 = Variable("W_{W_{coeff1}}", 8.71e-5, "1/m",
                      "Wing Weight Coefficent 1")
W_W_coeff2 = Variable("W_{W_{coeff2}}", 45.24, "Pa",
                      "Wing Weight Coefficent 2")
A = Variable("A", "-", "aspect ratio")
S = Variable("S", "m^2", "total wing area")
W_0 = Variable("W_0","N","Aircraft weight excluding wings")
W_w = Variable("W_w", "N", "wing weight")

N_ult = Variable('N_ult',3,'-','Ultimate load factor')
tau = Variable("\\tau", 0.12, "-", "airfoil thickness to chord ratio")

k = Variable("k", 1.2, "-", "form factor")
e = Variable("e", 0.95, "-", "Oswald efficiency factor")
C_f = Variable("C_f", "-", "skin friction coefficient")
S_wetratio = Variable("(\\frac{S}{S_{wet}})", 2.05, "-", "wetted area ratio")
C_L = Variable("C_L", "-", "Lift coefficent of wing")
C_L_max = Variable("C_L_max",0.8,'-','Max lift coefficient of wing')
C_D = Variable("C_D", "-", "Drag coefficient of wing")
C_D_fuse = Variable('C_D_fuse',0.1,'-','Drag coefficient of fuselage')
pi = Variable("\\pi", np.pi, "-", "half of the circle constant")
Re = Variable("Re", "-", "Reynold's number")

# Wing weight model
W_w_strc = W_W_coeff1*(N_ult*A**1.5*(W_0*W_i*S)**0.5)/tau
W_w_surf = W_W_coeff2 * S

constraints += [
	W_0 >= g*(m_payload+m_fuel),
	W_w >= W_w_surf + W_w_strc
]

# Wing aero model
C_D_wpar = k*C_f*S_wetratio
C_D_ind = C_L**2/(pi*A*e)

b = Variable("b","m","span")
b = (S*A)**0.5
c = b/A
# Engine propulsion model
# Select to OS 120AX
P = Variable('P','W','Engine cruise power')
P_max = Variable('P_max',2311.67,'W','Engine max power')
# PSFC = Variable('PSFC','kg/s/W','Engine specific fuel consumption') Covered by n_0
SpecificP_max = Variable('SpecificP_max',2600,'W/kg','Engine maximum specific power')
W_eng = Variable('W_eng',3.711*0.647,'N','Weight of engine')

# Propeller model
k_t = Variable('k_t',0.07,'-','Propeller thrust coefficient')
C_p = Variable('C_p',0.07,'-','Propeller power coefficient')
prop_diam = Variable('prop_diam',0.2,'m','Diameter of propeller')
n = Variable('n','1/s','Revolutions per second for propeller')
T = Variable('T','N','Prop thrust')
constraints+=[
	T == k_t*rho*(n**2)*(prop_diam**4),
	P == C_p*rho*(n**3)*(prop_diam**5),
	P <= 0.9*P_max,
	T >= D,
	# W_eng >= g*P_max/SpecificP_max
]

# Fuselage
W_fuse = Variable('W_fuse','N','Weight of fuselage')
constraints+=[W_fuse>=0.2*W_i]

# Boom
l_b = Variable('l_b','m','Length of boom from center of gravity')
I_b = Variable('I_b','m^4','Moment of inertia of boom')
t_b = Variable('t_b','m','Thickness of boom')
r_b = Variable('r_b','m','Radius of boom')
s_b = Variable('s_b','N/m^2','Tensile stress on boom')
# http://www.performance-composites.com/carbonfibre/mechanicalproperties_2.asp
s_b_limit = Variable('s_b_limit',600e6,'N/m^2','Boom tensile stress limit')
rho_CFRP = Variable('rho_CFRP',1600,'kg/m^3')
A_b = Variable('A_b','m^2','Cross sectional area of boom')
M_b = Variable('M_b','N*m','Max moment on boom at root')
D_b = Variable('D_b','N','Drag on boom')
W_b = Variable('W_b','N','Weight of boom')

# Tail surfaces
rho_XPS = Variable('rho_XPS',20,'kg/m^3')

W_v = Variable('W_v','N','Weight of vertical stabilizer')
W_h = Variable('W_h','N','Weight of horizontal stabilizer')

A_v = Variable('A_v',3,'-','Aspect ratio v')
S_v = Variable('S_v','m^2','Reference area v')
thickness_v = Variable('thickness_v',0.1,'-','Thickness v')
C_D_v = Variable('C_D_v',0.02,'-','C_d for vertical stabilizer') #NACA 0012
D_v = Variable('D_v','N','Drag of vertical tail')
b_v = (S_v*A_v)**0.5
c_v = b_v/A_v

A_h = Variable('A_h',5,'-','Aspect ratio h')
S_h = Variable('S_h','m^2','Reference area h')
thickness_h = Variable('thickness_h',0.1,'-','Thickness h')
C_D_h = Variable('C_D_h',0.02,'-','C_d for horizontal stabilizer') #NACA 0012
D_h = Variable('D_h','N','Drag of horizontal tail')
b_h = (S_h*A_h)**0.5
c_h = b_h/A_h

V_h = (S_h*l_b)/(S*c)
V_v = (S_v*l_b)/(S*b)

constraints+=[
	V_v >= 0.035,
	V_h >= 0.45,
	D_v >= 0.5*rho*S_v*C_D_v*V**2,
	D_h >= 0.5*rho*S_h*C_D_h*V**2,
	W_h >= g*rho_XPS*S_h*thickness_h*((S_h*A_h)**0.5)/A_h,
	W_v >= g*rho_XPS*S_v*thickness_v*((S_v*A_v)**0.5)/A_v
]

constraints += [
	t_b >= Variable('t_b_lower',0.001,'m','Lower limit for boom thickness'),
	I_b <= 3.141*(r_b**3)*t_b,
	M_b >= N_ult*l_b*0.5*0.8*S_h*rho*V**2, #Based on max Cl for NACA 0012
	s_b >= M_b*r_b/I_b,
	s_b <= s_b_limit/1.5
]

with gpkit.SignomialsEnabled():
	constraints+=[W_b >= g*rho_CFRP*l_b * 3.141* (t_b**2 + 2*r_b*t_b)]

constraints += [C_D >=C_D_fuse + C_D_wpar + C_D_ind]
constraints+=[
				C_L <= C_L_max,
				L <= 0.5*rho*S*C_L*V**2,
				W_i <= 0.5*rho*S*C_L_max*Vmin**2,
				D >= 0.5*rho*S*C_D*V**2 + D_h + D_v,
				Re <= (rho/mu)*V*(S/A)**0.5,
				C_f >= 0.074/Re**0.2,
]

# All weights
constraints += [
	W_i >= W_b + W_h + W_v+ W_fuse + W_w + W_eng + g*(m_payload + m_fuel),
	W_f >= W_b + W_h + W_v + W_fuse + W_w + W_eng + g*(m_payload),
	# (m_struc*g)/W_i >= 0.7,
	# (m_struc*g) <= m_limit*g
]


objective = Vmin**6/R

m = Model(objective,constraints)
m.debug()
sol = m.solve()
print sol.table()
print sol(R)/1000
print sol(b)
print sol(b/A)
print sol(Vmin)
print('horiz')
print(sol(A_h))
print(sol(b_h))
print(sol(c_h))
print(sol(S_h))
print('vert')
print(sol(A_v))
print(sol(b_v))
print(sol(c_v))
print(sol(S_v))
print('fuel')
print(sol(V_fuel))
tankCrossA = 0.15*0.12;
print(sol(V_fuel)/tankCrossA)
print float(sol(S).magnitude)
resultsDict = {'VHNJACDDXEA':float(sol(S).magnitude),'FDGQUUBYWFT':float(sol(b/A).magnitude),'ZVJTAUAEWVE':float(sol(b).magnitude),'CSWCUOQMTDT':float(sol(b).magnitude)}
updateOpenVSP(resultsDict)
