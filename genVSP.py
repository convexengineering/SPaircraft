from gpkit import Model, Variable
import numpy as np
from numpy import tan, cos, pi, arctan, arccos
import gpkit

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

def genDesFile(sol, i, b737800=True):

    sweep = arccos(sol('\cos(\Lambda)_Mission, Aircraft, Wing, WingNoStruct')[i])*180/np.pi

    # Wing descriptors
    b = sol('b')[i].to('m')
    croot = sol('c_{root}')[i].to('m')
    ctip = sol('c_{tip}')[i].to('m')
    S = sol('S')[i].to('m^2')
    xwing = sol('x_{wing}')[i].to('m')
    dihedral = 7.

    # Fuselage descriptors
    hfloor = sol('h_{floor}_Mission, Aircraft, Fuselage')[i].to('m')
    lnose = sol('l_{nose}_Mission, Aircraft, Fuselage')[i].to('m')
    lshell = sol('l_{shell}')[i].to('m')
    lcone = sol('l_{cone}')[i].to('m')
    lfloor = sol('l_{floor}')[i].to('m')
    lfuse = sol('l_{fuse}')[i].to('m')
    hfuse = sol('h_{fuse}')[i].to('m')
    wfuse = sol('w_{fuse}')[i].to('m')
    wfloor = sol('w_{floor}')[i].to('m')
    wdb = sol('w_{db}_Mission, Aircraft, Fuselage')[i].to('m')
    Rfuse = sol('R_{fuse}_Mission, Aircraft, Fuselage')[i].to('m')

    # Horizontal Tail descriptors
    xCGht = sol('x_{CG_{ht}}')[i].to('m')
    crootht = sol('c_{root_h}')[i].to('m')
    ctipht = sol('c_{tip_h}')[i].to('m')
    bht = sol('b_{ht}')[i].to('m')
    xCGht = sol('x_{CG_{ht}}')[i]
    lht = sol('l_{ht}')[i].to('m')
    tanht = sol('\\tan(\Lambda_{ht})_Mission, Aircraft, HorizontalTail, HorizontalTailNoStruct')[i]

    # Vertical Tail descriptors
    xCGvt = sol('x_{CG_{vt}}')[i].to('m')
    xtail = sol('x_{tail}')[i].to('m')
    Svt = sol('S_{vt}')[i].to('m^2')
    bvt = sol('b_{vt}')[i].to('m')
    lvt = sol('l_{vt}')[i].to('m')
    crootvt = sol('c_{root_{vt}}')[i].to('m')
    ctipvt = sol('c_{tip_{vt}}')[i].to('m')
    dxleadvt = sol('\\Delta x_{lead_v}')[i].to('m')
    dxtrailvt = sol('\\Delta x_{trail_v}')[i].to('m')
    tanvt = sol('\\tan(\Lambda_{vt})_Mission, Aircraft, VerticalTail, VerticalTailNoStruct')[i]

    # Engine descriptors
    df = sol('d_{f}')[i].to('m') # Engine frontal area
    lnace = sol('l_{nacelle}')[i].to('m')
    yeng = sol('y_{eng}')[i].to('m')

    # Things to integrate later
    # n_{rows}
    # wingbox (x_b) (x_f)

    resultsDict = {
        # Engine Variables
        'REBAHPKXPRR':float(xCGvt.magnitude), # Engine x location
        'GKMTRGNCEVD':float(yeng.magnitude), #Engine y location
        'XFTWTTHLVRI':float(hfuse.magnitude - (df/5.).magnitude), # Engine z location
        'JTPPOOJVVPE':float(lnace.magnitude),# Engine length
        'QRBDHPAPDFX':float(lnace.magnitude/df.magnitude), # Engine fineness ratio (set at 2 for now)

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
    print('File generation successful!')

def genDesFileSweep(sol, n, b737800=True):
    for i in range(0,1):
        solSingle = SolutionArray
        for j in sol['variables']:
            solSingle.update({j:sol(j)[i]})
