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

# m = Mission()
# m.substitutions.update(substitutions)
# sol = m.localsolve('mosek',verbosity=2)

# print sol.table()

# Vehicle descriptors
xCGmin = sol('x_{CG_{min}}')
xAC = sol('x_{AC}')

# Wing descriptors
b = sol('b')
croot = sol('c_{root}')
ctip = sol('c_{tip}')
S = sol('S')
xwing = sol('x_{wing}')

# Fuselage descriptors
lnose = sol('l_{nose}')
lshell = sol('l_{shell}')
lcone = sol('l_{cone}')
lfuse = sol('l_{fuse}')
hfuse = sol('h_{fuse}')
wfuse = sol('w_{fuse}')
wdb = sol('w_{db}')
Rfuse = sol('R_{fuse}')

# Horizontal Tail descriptors
xCGht = sol('x_{CG_{ht}}')
crootht = sol('c_{root_h}')
ctipht = sol('c_{tip_h}')
# dxleadht = sol('\\Delta x_{{lead}_h}')
lht = sol('l_{ht}')
tanht = sol('\\tan(\\Lambda_{ht})')

# Vertical Tail descriptors
xCGvt = sol('x_{CG_{vt}}')
xtail = sol('x_{tail}')
Svt = sol('S_{vt}')
bvt = sol('b_{vt}')
lvt = sol('l_{vt}')
crootvt = sol('c_{root_{vt}}')
ctipvt = sol('c_{tip_{vt}}')
dxleadvt = sol('\\Delta x_{lead_v}')
dxtrailvt = sol('\\Delta x_{trail_v}')
tanvt = sol('\\tan(\\Lambda_{vt})')

# Engine descriptors
A_2 = sol('A_2') # Engine frontal area


# List of variables still to integrate
# Root chord of wing (c_0)
# Length of nose cone (l_{nose})
# Length of tail cone (l_{cone})
#


# Things to integrate later
# Web half-weight (h_{db})
# h_{floor}
# n_{rows}
# wingbox (x_b) (x_f)



resultsDict = {
    # Wing variables
    'FLUGWUJQBVD':float(S.magnitude), # Wing surface area
    'UOBOGEWYYZZ':float((xwing - 0.5*croot).magnitude), # Wing x-location
    'MOGKYBMVMPD':float(-1*hfuse.magnitude + 0.2), # Wing z-location
    'NNIHPEXRTCP':float(croot.magnitude), # Wing root chord
    'HGZBRNOPIRD':float(ctip.magnitude), # Wing tip chord

    # Fuselage variables
    'KBKZBHMUHEP':float((lnose/lfuse).magnitude), # Nose location as % of fuse length
    'OVEJIBRDSBJ':float(1. - (lcone/lfuse).magnitude), # Tailcone location as % of fuse length

    # VT variables
    'LLYTEYDPDID':float(xCGmin.magnitude + dxleadvt.magnitude),# VT x location (LE location)
    # VT TE location (LE location + chord)
    'BFZDOVCXTAV':float(wfuse.magnitude/2), # VT y location (as wide as fuselage)
    'FQDVQTUBLUX':float(0.5),# VT z location (0.5 m off the widest point of the fuselage)

    # HT Variables
    'USGQFZQKJWC':float((xCGmin + dxleadvt + tanvt*bvt - 1.5*tanht*0.5*wfuse).magnitude),
    'BLMHVDOLAQJ':float(0.5 + 0.8*bvt.magnitude),

    # Engine Variables
    'REBAHPKXPRR':float(xCGvt.magnitude), # Engine x location
    'GKMTRGNCEVD':float((wdb + 0.4*Rfuse).magnitude), #Engine y location
    'XFTWTTHLVRI':float(hfuse.magnitude - (2*(A_2/pi)**0.5/10).magnitude), # Engine z location
    'JTPPOOJVVPE':float((2*(A_2/pi)**0.5).magnitude),# Engine length
    'QRBDHPAPDFX':float(2) # Engine fineness ratio (set at 2 for now)


}
updateOpenVSP(resultsDict)
