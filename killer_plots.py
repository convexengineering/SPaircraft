# only needed for plotting
import matplotlib.pyplot as plt

# Percent diffs
from percent_diff import percent_diff

from gpkit.small_scripts import mag

# Solution saving
from saveSol import genSolOut

# Models and substitutions
from SPaircraft import optimize_aircraft
from subs.optimalD8 import get_optimalD8_subs
from subs.optimal737 import get_optimal737_subs
from subs.M072_737 import get_M072_737_subs
from subs.D8_no_BLI import get_D8_no_BLI_subs
from subs.D8_eng_wing import get_D8_eng_wing_subs

def standard_killer_plot_max_opt_eng():
    # Note that the killer plot is designed for the following user inputs
    # Nclimb = 3, Ncruise = 2, Nmission = 1,
    # Rreq = 3000nmi, Npass = 180
    sol0, m0, m_relax0  = optimize_aircraft('W_{f_{total}}', 'optimal737', get_optimal737_subs(), True, False, False)
    wf0 = sol0('W_{f_{total}}')

    sol1, m1, m_relax1 = optimize_aircraft('W_{f_{total}}', 'M072_737', get_M072_737_subs(), True, False, False)
    wf1 = sol1('W_{f_{total}}')

    sol2, m2, m_relax2 = optimize_aircraft('W_{f_{total}}', 'D8_eng_wing', get_D8_eng_wing_subs(), True, False, False)
    wf2 = sol2('W_{f_{total}}')

    sol3, m3, m_relax3 = optimize_aircraft('W_{f_{total}}', 'D8_no_BLI', get_D8_no_BLI_subs(), True, False, False)
    wf3 = sol3('W_{f_{total}}')

    sol4, m4, m_relax4 = optimize_aircraft('W_{f_{total}}', 'optimalD8', get_optimalD8_subs(), True, False, False)
    wf4 = sol4('W_{f_{total}}')

    sol5, m5, m_relax5 = optimize_aircraft('W_{f_{total}}', 'optimalD8', get_optimalD8_subs(), False, True, False)
    wf5 = sol5('W_{f_{total}}')

    wing_sens = [sol0['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], sol1['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], \
                 sol2['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], sol3['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], \
                 sol4['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], sol5['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing']]
    HT_sens = [sol0['sensitivities']['constants']['C_{ht}_Mission/Aircraft/HorizontalTail'], sol1['sensitivities']['constants']['C_{ht}_Mission/Aircraft/HorizontalTail'], \
                 sol2['sensitivities']['constants']['C_{ht}_Mission/Aircraft/HorizontalTail'], sol3['sensitivities']['constants']['C_{ht}_Mission/Aircraft/HorizontalTail'], \
                 sol4['sensitivities']['constants']['C_{ht}_Mission/Aircraft/HorizontalTail'], sol5['sensitivities']['constants']['C_{ht}_Mission/Aircraft/HorizontalTail']]
    VT_sens = [sol0['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail'], sol1['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail'], \
                 sol2['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail'], sol3['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail'], \
                 sol4['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail'], sol5['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail']]
    fuse_sens = [sol0['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage'], sol1['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage'], \
                 sol2['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage'], sol3['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage'], \
                 sol4['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage'], sol5['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage']]
    engine_sens = [sol0['sensitivities']['constants']['C_{engsys}_Mission/Aircraft'], sol1['sensitivities']['constants']['C_{engsys}_Mission/Aircraft'], \
                 sol2['sensitivities']['constants']['C_{engsys}_Mission/Aircraft'], sol3['sensitivities']['constants']['C_{engsys}_Mission/Aircraft'], \
                 sol4['sensitivities']['constants']['C_{engsys}_Mission/Aircraft'], sol5['sensitivities']['constants']['C_{engsys}_Mission/Aircraft']]
    lg_sens = [sol0['sensitivities']['constants']['C_{lg}_Mission/Aircraft'], sol1['sensitivities']['constants']['C_{lg}_Mission/Aircraft'], \
                sol2['sensitivities']['constants']['C_{lg}_Mission/Aircraft'], sol3['sensitivities']['constants']['C_{lg}_Mission/Aircraft'], \
                sol4['sensitivities']['constants']['C_{lg}_Mission/Aircraft'], sol5['sensitivities']['constants']['C_{lg}_Mission/Aircraft']]
    Mmin_sens = [sol0['sensitivities']['constants']['M_{min}_Mission/Aircraft'], sol1['sensitivities']['constants']['M_{min}_Mission/Aircraft'], \
             sol2['sensitivities']['constants']['M_{min}_Mission/Aircraft'], sol3['sensitivities']['constants']['M_{min}_Mission/Aircraft'], \
             sol4['sensitivities']['constants']['M_{min}_Mission/Aircraft'], sol5['sensitivities']['constants']['M_{min}_Mission/Aircraft']]

    ytest = [1, wf1/wf0, wf2/wf0, wf3/wf0, wf4/wf0, wf5/wf0, wf5/wf0]
    for i in range(len(ytest)):
        if i == 0:
            ytest[i] = ytest[i]
        if i != 0:
            ytest[i] = mag(ytest[i][0])
        

    xtest = [0, 1, 2, 3, 4, 5, 6]
    xlabels = ['Optimized 737-800 M = 0.8', 'Slow to M = 0.72', 'D8 fuselage, Pi tail', 'Rear podded engines', 'Integrated engines, BLI = D8', 'Optimize engine', '2010 Engines']

    plt.plot(xtest, ytest, "o--")
    plt.plot([0, 1, 2, 3, 4, 5, 6], [1, .88, .81, .82, .67, .66, .63], "o--")
    plt.plot([0, 1, 2, 3, 6], [1, .868, .871, .865, .602], "o--")
    plt.plot([0, 1, 2, 3, 4, 5, 6], [1, 41129./43843, 38402./43843, 37180./43843, 32987./43843, 32383./43843, 29753./43843], "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
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
