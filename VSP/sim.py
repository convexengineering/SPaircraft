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
hfloor = sol('h_{floor}')
lnose = sol('l_{nose}')
lshell = sol('l_{shell}')
lcone = sol('l_{cone}')
lfloor = sol('l_{floor}')
lfuse = sol('l_{fuse}')
hfuse = sol('h_{fuse}')
wfuse = sol('w_{fuse}')
wdb = sol('w_{db}')
Rfuse = sol('R_{fuse}')

# Horizontal Tail descriptors
xCGht = sol('x_{CG_{ht}}')
crootht = sol('c_{root_h}')
ctipht = sol('c_{tip_h}')
bht = sol('b_{ht}')
xCGht = sol('x_{CG_{ht}}')
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
A_2 = sol('A_2_Mission, Aircraft, Engine') # Engine frontal area


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
    # Engine Variables
    'REBAHPKXPRR':float(xCGvt.magnitude), # Engine x location
    'GKMTRGNCEVD':float((wdb + 0.4*Rfuse).magnitude), #Engine y location
    'XFTWTTHLVRI':float(hfuse.magnitude - (2*(A_2/pi)**0.5/10).magnitude), # Engine z location
    'JTPPOOJVVPE':float((2*(A_2/pi)**0.5).magnitude),# Engine length
    'QRBDHPAPDFX':float(2), # Engine fineness ratio (set at 2 for now)

    # Floor Variables
    'MCVUEHMJBGG':float(hfloor.magnitude),  # Floor height
    'EVDJZIXRYSR':float(lfloor.magnitude), # Floor length
    'SKXVOFXEYEZ':float(wfuse.magnitude), # Floor width
    'KNNNINRHVVJ':float(lnose.magnitude-Rfuse.magnitude), # Floor x location (beginning of cyl section)
    'AFIOFOUHMWM':float(-0.5 - 0.5*hfloor.magnitude), # Floor z location (offset from thickest section)

    # Fuselage variables
    'HOVDTKFGFQC':float(lfuse.magnitude), # Fuselage length
    'KBKZBHMUHEP':float((lnose/lfuse).magnitude), # Nose location as % of fuse length
    'OVEJIBRDSBJ':float(1. - (lcone/lfuse).magnitude), # Tailcone location as % of fuse length
    'JMWPVNGZBYQ':float(2.0*hfuse.magnitude), # Fuselage height
    'KFWNCSRQOCQ':float(wfuse.magnitude), # Fuselage width
    'WKRLDITVGSF':float(2.0*hfuse.magnitude), # Fuselage height
    'TBCZTWFMJDM':float(wfuse.magnitude), # Fuselage width
    'JOBWSWPMZIB':float(2.0*hfuse.magnitude), # Fuselage height
    'HPKOTUWYSIY':float(wfuse.magnitude), # Fuselage width

    # HT Variables
    'USGQFZQKJWC':float(xCGht.magnitude - 0.5*crootht.magnitude - 1.5*tanht*0.5*wfuse.magnitude), # HT x location
    'BLMHVDOLAQJ':float(0.5 + bvt.magnitude),                                         # HT z location
    'IFZAMYYJPRP':float(30.),                                                             # HT sweep
    'CHYQUCYJMPS':float(bht.magnitude*0.5),                                               # HT half-span
    'LQXJHZEHDRX':float(crootht.magnitude),                                               # HT root chord
    'AYFSAELIRAY':float(ctipht.magnitude),                                                # HT tip chord

    # VT variables
    'LLYTEYDPDID':float(xCGvt.magnitude - 0.5*crootvt.magnitude), # VT x location (LE location)
    'BFZDOVCXTAV':float(wfuse.magnitude/2),                    # VT y location (as wide as fuselage)
    'FQDVQTUBLUX':float(0.5),                                  # VT z location (0.5 m off the widest point of the fuselage)
    'JXFRWSLYWDH':float(bvt.magnitude),                        # VT span
    'MBZGSEIYFGW':float(crootvt.magnitude),                    # VT root chord
    'CUIMIUZJQMS':float(ctipvt.magnitude),                     # VT tip chord

    # Wing variables
    'AYJHHOVUHBI':float(b.magnitude*0.5), # Wing half-span
    'UOBOGEWYYZZ':float((xwing - 0.5*croot).magnitude), # Wing x-location
    'MOGKYBMVMPD':float(-1*hfuse.magnitude + 0.2), # Wing z-location
    'NNIHPEXRTCP':float(croot.magnitude), # Wing root chord
    'HGZBRNOPIRD':float(ctip.magnitude), # Wing tip chord
}
updateOpenVSP(resultsDict)
