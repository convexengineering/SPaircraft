from gpkit import Model, Variable
import numpy as np
from numpy import tan, cos, pi, arctan, arccos
import gpkit

def updateOpenVSP(inputDict, i = 0):
    filenameOpen = 'VSP/design.des'
    filenameWrite = 'VSP/design' + str(i) + '.des'

    with open(filenameOpen,'r+') as f:
        g = open(filenameWrite,'w+')
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
        f.close()
        print('OpenVSP .des output:')
        print(output)
        g.seek(0)
        g.write(output)
        g.truncate()
        g.close()

def genDesFile(sol, aircraft = 'D82', i = 0, swpt = False):
    if swpt:
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
        crootht = sol('c_{root_{ht}}')[i].to('m')
        ctipht = sol('c_{tip_{ht}}')[i].to('m')
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
        # dxleadvt = sol('\\Delta x_{lead_v}')[i].to('m')
        # dxtrailvt = sol('\\Delta x_{trail_v}')[i].to('m')
        tanvt = sol('\\tan(\Lambda_{vt})_Mission, Aircraft, VerticalTail, VerticalTailNoStruct')[i]

        # Engine descriptors
        df = sol('d_{f}_Mission, Aircraft, Engine')[i].to('m') # Engine frontal area
        lnace = sol('l_{nacelle}')[i].to('m')
        yeng = sol('y_{eng}_Mission, Aircraft, VerticalTail, VerticalTailNoStruct')[i].to('m')

        # Things to integrate later
        # n_{rows}
        # wingbox (x_b) (x_f)
    else:
        sweep = arccos(sol('\cos(\Lambda)_Mission, Aircraft, Wing, WingNoStruct'))*180/np.pi

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
        # dxleadvt = sol('\\Delta x_{lead_v}').to('m')
        # dxtrailvt = sol('\\Delta x_{trail_v}').to('m')
        tanvt = sol('\\tan(\Lambda_{vt})_Mission, Aircraft, VerticalTail, VerticalTailNoStruct')

        # Engine descriptors
        df = sol('d_{f}_Mission, Aircraft, Engine').to('m') # Engine frontal area
        lnace = sol('l_{nacelle}').to('m')
        yeng = sol('y_{eng}_Mission, Aircraft, VerticalTail, VerticalTailNoStruct').to('m')

    # Creating the default (D82) resultsDict
    resultsDict = {
        # Engine Variables
        'OOWZWGGROQZ':float(lnace.magnitude),    # Engine length (chord)
        'TTRJCLVSWWP':float(df.magnitude),       # Engine height
        'YUWFYBTYKTL':float(0.1625),             # Engine airfoil thickness/chord
        'TVQVWMMVRYB':float(df.magnitude),       # Engine width
        'EGCVYPSLWEZ':float(xCGvt.magnitude+0.25*crootvt.magnitude),    # Engine x location
        'RJLYSBJAFOT':float(0.5*wfuse.magnitude),     #Engine y location
        'GBGVQARDEVD':float(hfuse.magnitude - (df/5.).magnitude), # Engine z location
        'HKVDGHIEXRW':float(15.),                                  # Engine up-rotation (degrees)

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
        'USGQFZQKJWC':float(float(xCGvt.magnitude) - 0.5*crootvt.magnitude + 1.0*tanvt*bvt.magnitude - wfuse.magnitude*tanht), # HT x location
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
    # if aircraft in ['D8big', 'D82_73eng', 'D8_eng_wing', 'optimalD8', 'M08D8', 'M08_D8_eng_wing']:

    # Wing mounted engines
    if aircraft in ['b737800', 'b777300ER','D8_eng_wing','M08_D8_eng_wing','optimal737','optimal777','optimalRJ']:
        resultsDict.update({
         # Engine Variables
        'EGCVYPSLWEZ':float((xwing).magnitude - 0.25*croot.magnitude + yeng.magnitude*tan(sweep*pi/180)), # Engine x location
        'RJLYSBJAFOT':float(yeng.magnitude), #Engine y location
        'GBGVQARDEVD':float(-hfuse.magnitude - 0.2*df.magnitude), # Engine z location
        'HKVDGHIEXRW':float(0.),                                  # Engine up-rotation (degrees)
        })
    # Conventional tail
    if aircraft in ['b737800', 'b777300ER','optimal737','optimal777','optimalRJ']:
        resultsDict.update({
        # HT Variables
        'USGQFZQKJWC':float(float(xCGht.magnitude) - crootht.magnitude),     # HT x location
        'BLMHVDOLAQJ':float(0.),                                             # HT z location
        'CHYQUCYJMPS':float(bht.magnitude*0.5 + wfuse.magnitude),            # HT half-span

        # VT variables
        'LLYTEYDPDID':float(xCGvt.magnitude - 0.5*crootvt.magnitude),        # VT x location (LE location)
        'BFZDOVCXTAV':float(0.),                                             # VT y location (as wide as fuselage)
        'FQDVQTUBLUX':float(hfuse.magnitude),                                # VT z location (0.5 m off the widest point of the fuselage)
        'GWTZZGTPXQU':float(0.),                                             # VT dihedral
        })
    # Rear mounted non-BLI D8 engines
    if aircraft in ['M08_D8_no_BLI', 'D8_no_BLI']:
        resultsDict.update({
            'RJLYSBJAFOT':float(yeng.magnitude), #Engine y location
            'GBGVQARDEVD':float(0.0), # Engine z location
            'EGCVYPSLWEZ':float(xCGvt.magnitude), # Engine x location
        })

    updateOpenVSP(resultsDict,i)
    print('File generation successful!')

def genDesFileSweep(sol, aircraft, n):
    for i in range(0,n):
        genDesFile(sol,True,i,aircraft)

