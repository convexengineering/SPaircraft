from __future__ import absolute_import
# only needed for plotting
import matplotlib.pyplot as plt
import numpy as np

# Percent diffs
from .percent_diff import percent_diff
from gpkit import units
from gpkit.small_scripts import mag

# Solution saving
from .saveSol import genSolOut

# Models and substitutions
from .aircraft import Mission
from .SPaircraft import optimize_aircraft
from .subs.optimalD8 import get_optimalD8_subs
from .subs.optimal737 import get_optimal737_subs
from .subs.M072_737 import get_M072_737_subs
from .subs.D8_no_BLI import get_D8_no_BLI_subs
from .subs.D8_eng_wing import get_D8_eng_wing_subs

def standard_killer_plot():
    """
    Generates the standard killer plots from the TASOPT paper
    """
    Nclimb = 3; Ncruise = 2; Nmission = 1;
    subsList = [get_optimal737_subs(), get_M072_737_subs(), get_D8_eng_wing_subs(), get_D8_no_BLI_subs(), get_optimalD8_subs(), get_optimalD8_subs()]
    configList = ['optimal737', 'M072_737', 'D8_eng_wing', 'D8_no_BLI', 'optimalD8', 'optimalD8']
    fixedBPRList = [True, True, True, True, True, False]
    pRatOptList = [False, False, False, False, False, True]
    mutategpargList = [False, False, False, False, False, False]
    sol = {}; wf = [];
    for i in range(0,6):
        m = Mission(Nclimb, Ncruise, configList[i], Nmission)
        m.cost = m['W_{f_{total}}'].sum()
        substitutions = subsList[i]
        substitutions.update({'R_{req}': 3000.*units('nmi'),
                              'n_{pass}': 180.})
        sol[i] = optimize_aircraft(m, substitutions, fixedBPRList[i], pRatOptList[i], mutategpargList[i])
        wf.append(sol[i]('W_{f_{total}}'))

    wing_sens = [sol[i]['sensitivities']['constants']['C_{wing}'] for i in range(0,6)]
    HT_sens = [sol[i]['sensitivities']['constants']['C_{ht}'] for i in range(0,6)]
    VT_sens = [sol[i]['sensitivities']['constants']['C_{VT}'] for i in range(0,6)]
    fuse_sens = [sol[i]['sensitivities']['constants']['C_{fuse}'] for i in range(0,6)]
    engine_sens = [sol[i]['sensitivities']['constants']['C_{engsys}'] for i in range(0,6)]
    lg_sens = [sol[i]['sensitivities']['constants']['C_{lg}'] for i in range(0,6)]
    Mmin_sens = [sol[i]['sensitivities']['constants']['M_{min}'] for i in range(0,6)]

    ytest = [mag(wf[i]/wf[0])[0] for i in range(0,6)]
    xtest = [0, 1, 2, 3, 4, 5]
    xlabels = ['Optimized 737-800 M = 0.8', 'Slow to M = 0.72', 'D8 fuselage, Pi tail', 'Rear podded engines', 'Integrated engines, BLI = D8', 'Optimize engine', '2020 Engines']

    plt.plot(xtest, ytest, "o--")
    plt.plot([0, 1, 2, 3, 4, 5, 6], [1, .88, .81, .82, .67, .66, .63], "o--")
    plt.plot([0, 1, 2, 3, 6], [1, .868, .871, .865, .602], "o--")
    plt.plot([0, 1, 2, 3, 4, 5, 6], [1, 41129./43843, 38402./43843, 37180./43843, 32987./43843, 32383./43843, 29753./43843], "o--")
    plt.xticks(np.linspace(0,6,7), xlabels,  rotation='vertical')
    plt.ylim([0,1.1])
    plt.xlim([-.5, 6.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('$W_{\mathrm{f}}/W_{\mathrm{f}_\mathrm{0}}$', fontsize = 20)
    plt.title('D8 Morphing Chart')
    plt.legend(['SP Model', 'TASOPT', 'NASA', 'Aurora'], loc=3)
    plt.savefig('Morphing_Chart_Figs/D8_standard_morphing_chart.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    xtest = [0, 1, 2, 3, 4, 5]

    plt.plot(xtest, wing_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.25])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Wing Weight', fontsize = 20)
    plt.title('Wing Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_standard_killer_chart_max_opt_eng_wing_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, HT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.02])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Horizontal Tail Weight', fontsize = 20)
    plt.title('Horizontal Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_standard_killer_chart_max_opt_eng_HT_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, VT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.15])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Vertical Tail Weight', fontsize = 20)
    plt.title('Vertical Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_standard_killer_chart_max_opt_eng_VT_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, fuse_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.55])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Fuselage Weight', fontsize = 20)
    plt.title('Fuselage Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_standard_killer_chart_max_opt_eng_fuse_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, engine_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.3])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Engine Weight', fontsize = 20)
    plt.title('Engine Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_standard_killer_chart_max_opt_eng_engine_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, lg_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.07])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Landing Gear Weight', fontsize = 20)
    plt.title('Landing Gear Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_standard_killer_chart_max_opt_eng_lg_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()
