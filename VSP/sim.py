from gpkit import Model, Variable
import numpy as np
from numpy import tan, cos, pi, arctan
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
# xCGmin = sol('x_{CG_{min}}').to('m')
# xAC = sol('x_{AC}').to('m')

# varkeys = [k.str_without('models') for k in sol.program.gps[-1].varlocs if "Relax" not in k.models]
solNR = sol
# varkeys = [i.str_without]
# for i in varkeys:
#     if '(' not in i:
#         for j in i.iterkeys:
#             if 'Relax' not in j:
#                 solNR.update({i:sol(i).magnitude})

# for i in solNR:



# Wing descriptors
b = sol('b').to('m')
croot = sol('c_{root}').to('m')
ctip = sol('c_{tip}').to('m')
S = sol('S').to('m^2')
xwing = sol('x_{wing}').to('m')
dihedral = 7.

# Fuselage descriptors
hfloor = sol('h_{floor}_Mission, Aircraft, Fuselage').to('m')
lnose = sol('l_{nose}_Mission, Aircraft, Fuselage').to('m')
lshell = sol('l_{shell}').to('m')
lcone = sol('l_{cone}').to('m')
lfloor = sol('l_{floor}').to('m')
lfuse = sol('l_{fuse}').to('m')
hfuse = sol('h_{fuse}').to('m')
wfuse = sol('w_{fuse}').to('m')
wfloor = sol('w_{floor}').to('m')
wdb = sol('w_{db}_Mission, Aircraft, Fuselage').to('m')
Rfuse = sol('R_{fuse}_Mission, Aircraft, Fuselage').to('m')

# Horizontal Tail descriptors
xCGht = sol('x_{CG_{ht}}').to('m')
crootht = sol('c_{root_{ht}}').to('m')
ctipht = sol('c_{tip_{ht}}').to('m')
bht = sol('b_{ht}').to('m')
xCGht = sol('x_{CG_{ht}}')
lht = sol('l_{ht}').to('m')
tanht = sol('\\tan(\Lambda_{ht})_Mission, Aircraft, HorizontalTail, HorizontalTailNoStruct')

# Vertical Tail descriptors
xCGvt = sol('x_{CG_{vt}}').to('m')
xtail = sol('x_{tail}').to('m')
Svt = sol('S_{vt}').to('m^2')
bvt = sol('b_{vt}').to('m')
lvt = sol('l_{vt}').to('m')
crootvt = sol('c_{root_{vt}}').to('m')
ctipvt = sol('c_{tip_{vt}}').to('m')
dxleadvt = sol('\\Delta x_{lead_v}').to('m')
dxtrailvt = sol('\\Delta x_{trail_v}').to('m')
tanvt = sol('\\tan(\Lambda_{vt})_Mission, Aircraft, VerticalTail, VerticalTailNoStruct')

# Engine descriptors
df = sol('d_{f}') # Engine frontal area
lnace = sol('l_{nacelle}')
yeng = sol('y_{eng}')


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
    'GKMTRGNCEVD':float(yeng.magnitude), #Engine y location
    'XFTWTTHLVRI':float(hfuse.magnitude - (df/5.).magnitude), # Engine z location
    'JTPPOOJVVPE':float(lnace.magnitude),# Engine length
    'QRBDHPAPDFX':float(2.*lnace.magnitude/df.magnitude), # Engine fineness ratio (set at 2 for now)

    # Floor Variables
    'MCVUEHMJBGG':float(hfloor.magnitude),  # Floor height
    'EVDJZIXRYSR':float(lfloor.magnitude), # Floor length
    'SKXVOFXEYEZ':float(2*wfloor.magnitude), # Floor width
    'KNNNINRHVVJ':float(lnose.magnitude-Rfuse.magnitude), # Floor x location (beginning of cyl section)
    'AFIOFOUHMWM':float(-0.5 - 0.5*hfloor.magnitude), # Floor z location (offset from thickest section)

    # Fuselage variables
    'HOVDTKFGFQC':float(lfuse.magnitude), # Fuselage length
    'KBKZBHMUHEP':float((lnose/lfuse).magnitude), # Nose location as % of fuse length
    'OVEJIBRDSBJ':float(1. - (lcone/lfuse).magnitude), # Tailcone location as % of fuse length
    'JMWPVNGZBYQ':float(2.0*hfuse.magnitude), # Fuselage height
    'KFWNCSRQOCQ':float(2*wfuse.magnitude), # Fuselage width
    'WKRLDITVGSF':float(2.0*hfuse.magnitude), # Fuselage height
    'TBCZTWFMJDM':float(2*wfuse.magnitude), # Fuselage width
    'JOBWSWPMZIB':float(2.0*hfuse.magnitude), # Fuselage height
    'HPKOTUWYSIY':float(2*wfuse.magnitude), # Fuselage width

    # HT Variables
    'USGQFZQKJWC':float(float(xCGht.magnitude) - 0.5*crootht.magnitude - 1.0*tanvt*bvt.magnitude), # HT x location
    'BLMHVDOLAQJ':float(0.5 + bvt.magnitude),                                             # HT z location
    'IFZAMYYJPRP':float(arctan(tanht)*180/pi),                                                             # HT sweep
    'CHYQUCYJMPS':float(bht.magnitude*0.5),                                               # HT half-span
    'LQXJHZEHDRX':float(crootht.magnitude),                                               # HT root chord
    'AYFSAELIRAY':float(ctipht.magnitude),                                                # HT tip chord
    'IFZAMYYJPRP':float(np.arctan(tanht)*180/np.pi),                            # HT sweep angle

    # VT variables
    'LLYTEYDPDID':float(xCGvt.magnitude - 0.5*crootvt.magnitude), # VT x location (LE location)
    'BFZDOVCXTAV':float(2*wfuse.magnitude/2),                    # VT y location (as wide as fuselage)
    'FQDVQTUBLUX':float(0.5),                                  # VT z location (0.5 m off the widest point of the fuselage)
    'JXFRWSLYWDH':float(bvt.magnitude),                        # VT span
    'MBZGSEIYFGW':float(crootvt.magnitude),                    # VT root chord
    'CUIMIUZJQMS':float(ctipvt.magnitude),                     # VT tip chord
    'XLPAIOGKILI':float(arctan(tanvt)*180/pi),                 # VT sweep angle
    'GWTZZGTPXQU':-10,                                          # VT dihedral

    # Wing variables
    'AYJHHOVUHBI':float(b.magnitude*0.5), # Wing half-span
    'UOBOGEWYYZZ':float((xwing - 0.25*croot).magnitude), # Wing x-location
    'MOGKYBMVMPD':float(-1*hfuse.magnitude + 0.2), # Wing z-location
    'NNIHPEXRTCP':float(croot.magnitude), # Wing root chord
    'HGZBRNOPIRD':float(ctip.magnitude), # Wing tip chord
    'AGOKGLSLBTO':float(sweep), # Wing sweep angle
    'SMCAVCZXJSG':float(+dihedral), # Wing dihedral
}

#Differentiating between b737800 and D8
if b737800:
    resultsDict.update({
     # Engine Variables
    'REBAHPKXPRR':float((xwing-0.6*lnace).magnitude), # Engine x location
    'GKMTRGNCEVD':float(yeng.magnitude), #Engine y location
    'XFTWTTHLVRI':float(-hfuse.magnitude - 0.6*df.magnitude), # Engine z location
    # Floor Variables
    # Fuselage variables
    # HT Variables
    'USGQFZQKJWC':float(float(xCGht.magnitude) - 0.5*crootht.magnitude), # HT x location
    'BLMHVDOLAQJ':float(0.),                                             # HT z location
    # VT variables
    'LLYTEYDPDID':float(xCGvt.magnitude - 0.5*crootvt.magnitude),        # VT x location (LE location)
    'BFZDOVCXTAV':float(0.),                                             # VT y location (as wide as fuselage)
    'FQDVQTUBLUX':float(hfuse.magnitude),                                # VT z location (0.5 m off the widest point of the fuselage)
    'GWTZZGTPXQU':float(0.),                                             # VT dihedral
    # Wing variables
    })
# if D80 or D82:

updateOpenVSP(resultsDict)
