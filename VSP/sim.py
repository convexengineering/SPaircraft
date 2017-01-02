from gpkit import Model, Variable
import numpy as np
from numpy import tan, cos, pi
import gpkit
from CFP_Fuselage_Performance_int_HT import Mission

def updateOpenVSP(inputDict):
	filename = 'VSP/design.des'
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

sweep = 30
substitutions = {
        # 'V_{stall}'   : 120,
        '\\delta_P_{over}': 12,
        'N_{land}': 6,
        'SPR': 8,
        'p_s': 81.,
        'ReqRng': 1000,
        # '\\theta_{db}' : 0.366,
        'CruiseAlt': 30000,
        'numeng': 2,
        'n_{pax}': 150,
        'W_{avg. pass}': 180,
        'W_{carry on}': 15,
        'W_{cargo}': 10000,
        'W_{checked}': 40,
        'w_{aisle}': 0.51,
        'w_{seat}': 0.5,
        'w_{sys}': 0.1,
        'W_{cargo}': 10000,
        'r_E': 1,  # [TAS]
        '\\lambda_{cone}': 0.4,  # [Philippe]
        '\\rho_{cone}': 2700,  # [TAS]
        '\\rho_{bend}': 2700,  # [TAS]
        '\\rho_{floor}': 2700,  # [TAS]
        '\\rho_{skin}': 2700,  # [TAS]
        'W\'\'_{floor}': 60,  # [TAS]
        'W\'\'_{insul}': 22,  # [TAS]
        'W\'_{seat}': 150,  # [TAS]
        'W\'_{window}': 145. * 3,  # [TAS]
        'f_{fadd}': 0.2,  # [TAS]
        'f_{frame}': 0.25,  # [Philippe]
        'f_{lugg,1}': 0.4,  # [Philippe]
        'f_{lugg,2}': 0.1,  # [Philippe]
        'f_{string}': 0.1,  #TODO remove
        'f_{padd}': 0.4,  # [TAS]

        # wing subs

        'C_{L_{wmax}}': 2.5,
        '\\tan(\\Lambda)': tan(sweep * pi / 180),
        '\\alpha_{max,w}': 0.1,  # (6 deg)
        '\\cos(\\Lambda)': cos(sweep * pi / 180),
        '\\eta': 0.97,
        '\\rho_0': 1.225,
        '\\rho_{fuel}': 817,  # Kerosene [TASOPT]

        #VT subs
       'C_{D_{wm}}': 0.5, # [2]
       'C_{L_{vmax}}': 2.6, # [2]
       'V_1': 70,
       '\\rho_{TO}': 1.225,
        '\\tan(\\Lambda_{vt})': tan(40*pi/180),
        'c_{l_{vtEO}}': 0.5,
        'e_v': 0.8,
        'y_{eng}': 4.83, # [3]

        'V_{land}': 72,
        'I_{z}': 12495000, #estimate for late model 737 at max takeoff weight (m l^2/12)
        '\\dot{r}_{req}': 0.174533, #10 deg/s yaw rate
        'N_{spar}': 2,

        #HT subs
        '\\alpha_{max,h}': 2.5,
        '\\tan(\\Lambda_{ht})': tan(30*pi/180),
        'C_{L_{hmax}}': 2.5,
        'SM_{min}': 0.05,
}

m = Mission()
m.substitutions.update(substitutions)
sol = m.localsolve('mosek',verbosity=2)

print sol.table()
b = sol('b')
croot = sol('c_{root}')
S = sol('S')
xwing = sol('x_{wing}')
# print(sol(S_v))
resultsDict = {'FLUGWUJQBVD':float(S.magnitude),'UOBOGEWYYZZ':float((xwing + 0.5*croot).magnitude)}
updateOpenVSP(resultsDict)
