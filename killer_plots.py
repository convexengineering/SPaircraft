# only needed for plotting
import matplotlib.pyplot as plt

from SP_aircraft import run_737800
from SP_aircraft import run_D82
from SP_aircraft import run_M072_737
from SP_aircraft import run_D8_eng_wing
from SP_aircraft import run_optimal_737
from SP_aircraft import run_optimal_D8
from SP_aircraft import run_D8_no_BLI
from SP_aircraft import run_M08_D8_eng_wing
from SP_aircraft import run_M08_D8_no_BLI
from SP_aircraft import run_M08_D8
from SP_aircraft import run_optimal_777
from SP_aircraft import run_D8_big_eng_wing
from SP_aircraft import run_D8_big_no_BLI
from SP_aircraft import run_D8_big
from SP_aircraft import run_D8_big_M072
from SP_aircraft import run_D8_big_M08
from SP_aircraft import run_optimal_RJ
from SP_aircraft import run_M072_optimal_RJ
from SP_aircraft import run_small_D8_eng_wing
from SP_aircraft import run_small_D8_no_BLI
from SP_aircraft import run_small_D8
from SP_aircraft import run_small_M08_D8_eng_wing
from SP_aircraft import run_small_M08_D8_no_BLI
from SP_aircraft import run_small_M08_D8

def RJ_M08_killer_plot_standard():
    sol0 = run_optimal_RJ(True)
    wf0 = sol0('W_{f_{total}}')

    sol1 = run_small_M08_D8_eng_wing(True)
    wf1 = sol1('W_{f_{total}}')

    sol2 = run_small_M08_D8_no_BLI(True)
    wf2 = sol2('W_{f_{total}}')

    sol3 = run_small_M08_D8(True)
    wf3 = sol3('W_{f_{total}}')

    sol4 = run_small_D8(True)
    wf4 = sol4('W_{f_{total}}')

    sol5 = run_small_D8(False, True)
    wf5 = sol5('W_{f_{total}}')

    ytest = [1, wf1/wf0, wf2/wf0, wf3/wf0, wf4/wf0, wf5/wf0]
    xtest = [0, 1, 2, 3, 4, 5]
    xlabels = ['Optimized RJ M = 0.8', 'D8 fuselage, Pi tail', 'Rear podded engines', 'Integrated engines, BLI = D8', 'Slow to M = 0.72',  'Optimize Engine']

    wing_sens = [sol0['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], sol1['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], \
                 sol2['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], sol3['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], \
                 sol4['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], sol5['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing']]
    HT_sens = [sol0['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], sol1['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], \
                 sol2['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], sol3['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], \
                 sol4['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], sol5['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail']]
    VT_sens = [sol0['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], sol1['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], \
                 sol2['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], sol3['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], \
                 sol4['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], sol5['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail']]
    fuse_sens = [sol0['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], sol1['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], \
                 sol2['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], sol3['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], \
                 sol4['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], sol5['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage']]
    engine_sens = [sol0['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], sol1['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], \
                 sol2['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], sol3['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], \
                 sol4['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], sol5['sensitivities']['constants']['C_{engsys}_Mission, Aircraft']]
    lg_sens = [sol0['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], sol1['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], \
                 sol2['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], sol3['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], \
                 sol4['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], sol5['sensitivities']['constants']['C_{lg}_Mission, Aircraft']]

    plt.plot(xtest, ytest, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,1.1])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('$W_{f}/W_{f_0}$', fontsize = 20)
    plt.title('RJ Class D8 Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_M08_standard_killer_chart.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, wing_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.2])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Wing Weight', fontsize = 20)
    plt.title('Wing Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_M08_standard_killer_chart_wing_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, HT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.015])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Horizontal Tail Weight', fontsize = 20)
    plt.title('Horizontal Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_M08_standard_killer_chart_HT_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, VT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.06])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Vertical Tail Weight', fontsize = 20)
    plt.title('Vertical Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_M08_standard_killer_chart_VT_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, fuse_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.55])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Fuselage Weight', fontsize = 20)
    plt.title('Fuselage Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_M08_standard_killer_chart_fuse_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, engine_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.25])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Engine Weight', fontsize = 20)
    plt.title('Engine Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_M08_standard_killer_chart_engine_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, lg_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.07])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Landing Gear Weight', fontsize = 20)
    plt.title('Landing Gear Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_M08_standard_killer_chart_lg_sens.pdf', bbox_inches="tight")
    plt.show()

def RJ_killer_plot_standard():
    sol0 = run_optimal_RJ(True)
    wf0 = sol0('W_{f_{total}}')

    sol1 = run_M072_optimal_RJ(True)
    wf1 = sol1('W_{f_{total}}')

    sol2 = run_small_D8_eng_wing(True)
    wf2 = sol2('W_{f_{total}}')

    sol3 = run_small_D8_no_BLI(True)
    wf3 = sol3('W_{f_{total}}')

    sol4 = run_small_D8(True)
    wf4 = sol4('W_{f_{total}}')

    sol5 = run_small_D8(False, True)
    wf5 = sol5('W_{f_{total}}')

    wing_sens = [sol0['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], sol1['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], \
                 sol2['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], sol3['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], \
                 sol4['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], sol5['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing']]
    HT_sens = [sol0['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], sol1['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], \
                 sol2['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], sol3['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], \
                 sol4['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], sol5['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail']]
    VT_sens = [sol0['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], sol1['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], \
                 sol2['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], sol3['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], \
                 sol4['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], sol5['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail']]
    fuse_sens = [sol0['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], sol1['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], \
                 sol2['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], sol3['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], \
                 sol4['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], sol5['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage']]
    engine_sens = [sol0['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], sol1['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], \
                 sol2['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], sol3['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], \
                 sol4['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], sol5['sensitivities']['constants']['C_{engsys}_Mission, Aircraft']]
    lg_sens = [sol0['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], sol1['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], \
                 sol2['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], sol3['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], \
                 sol4['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], sol5['sensitivities']['constants']['C_{lg}_Mission, Aircraft']]

    ytest = [1, wf1/wf0, wf2/wf0, wf3/wf0, wf4/wf0, wf5/wf0]
    xtest = [0, 1, 2, 3, 4, 5]
    xlabels = ['Optimized RJ M = 0.8', 'Slow to M = 0.72', 'D8 fuselage, Pi tail', 'Rear podded engines', 'Integrated engines, BLI = D8', 'Optimize Engine']

    plt.plot(xtest, ytest, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,1.1])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('$W_{f}/W_{f_0}$', fontsize = 20)
    plt.title('RJ Class D8 Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_standard_killer_chart.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, wing_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.25])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Wing Weight', fontsize = 20)
    plt.title('Wing Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_killer_chart_standard_wing_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, HT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.02])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Horizontal Tail Weight', fontsize = 20)
    plt.title('Horizontal Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_killer_chart_standard_HT_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, VT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.1])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Vertical Tail Weight', fontsize = 20)
    plt.title('Vertical Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_killer_chart_standard_VT_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, fuse_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.55])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Fuselage Weight', fontsize = 20)
    plt.title('Fuselage Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_killer_chart_standard_fuse_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, engine_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.25])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Engine Weight', fontsize = 20)
    plt.title('Engine Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_killer_chart_standard_engine_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, lg_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.07])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Landing Gear Weight', fontsize = 20)
    plt.title('Landing Gear Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_killer_chart_standard_lg_sens.pdf', bbox_inches="tight")
    plt.show()

def RJ_killer_plot_max_optimal_engines():
    #TURN ON THE RIGHT CONSTRIANTS INSIDE D8.py
    sol0 = run_optimal_RJ(False, True)
    wf0 = sol0('W_{f_{total}}')

    sol1 = run_M072_optimal_RJ(False, True)
    wf1 = sol1('W_{f_{total}}')

    sol2 = run_small_D8_eng_wing(False, True)
    wf2 = sol2('W_{f_{total}}')

    sol3 = run_small_D8_no_BLI(False, True)
    wf3 = sol3('W_{f_{total}}')

    sol4 = run_small_D8(False, True)
    wf4 = sol4('W_{f_{total}}')

    wing_sens = [sol0['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], sol1['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], \
                 sol2['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], sol3['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], \
                 sol4['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing']]
    HT_sens = [sol0['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], sol1['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], \
                 sol2['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], sol3['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], \
                 sol4['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail']]
    VT_sens = [sol0['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], sol1['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], \
                 sol2['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], sol3['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], \
                 sol4['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail']]
    fuse_sens = [sol0['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], sol1['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], \
                 sol2['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], sol3['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], \
                 sol4['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage']]
    engine_sens = [sol0['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], sol1['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], \
                 sol2['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], sol3['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], \
                 sol4['sensitivities']['constants']['C_{engsys}_Mission, Aircraft']]
    lg_sens = [sol0['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], sol1['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], \
                 sol2['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], sol3['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], \
                 sol4['sensitivities']['constants']['C_{lg}_Mission, Aircraft']]

    ytest = [1, wf1/wf0, wf2/wf0, wf3/wf0, wf4/wf0]
    xtest = [0, 1, 2, 3, 4]
    xlabels = ['Optimized RJ M = 0.8', 'Slow to M = 0.72', 'D8 fuselage, Pi tail', 'Rear podded engines', 'Integrated engines, BLI = D8']

    plt.plot(xtest, ytest, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,1.1])
    plt.xlim([-.5, 4.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('$W_{f}/W_{f_0}$', fontsize = 20)
    plt.title('RJ Class D8 Morphing Chart - Maximum Engine Optimization')
    plt.savefig('Morphing_Chart_Figs/RJ_killer_chart_max_optimal_engines.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, wing_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.25])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Wing Weight', fontsize = 20)
    plt.title('Wing Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_killer_chart_max_optimal_engines_wing_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, HT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.02])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Horizontal Tail Weight', fontsize = 20)
    plt.title('Horizontal Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_killer_chart_max_optimal_engines_HT_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, VT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.15])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Vertical Tail Weight', fontsize = 20)
    plt.title('Vertical Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_killer_chart_max_optimal_engines_VT_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, fuse_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.55])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Fuselage Weight', fontsize = 20)
    plt.title('Fuselage Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_killer_chart_max_optimal_engines_fuse_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, engine_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.3])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Engine Weight', fontsize = 20)
    plt.title('Engine Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_killer_chart_max_optimal_engines_engine_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, lg_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.07])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Landing Gear Weight', fontsize = 20)
    plt.title('Landing Gear Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_killer_chart_max_optimal_engines_lg_sens.pdf', bbox_inches="tight")
    plt.show()

def 777_killer_plot_standard():
    sol0 = run_optimal_777(True)
    wf0 = sol0('W_{f_{total}}')

    sol1 = run_D8_big_eng_wing(True)
    wf1 = sol1('W_{f_{total}}')

    sol2 = run_D8_big_no_BLI(True)
    wf2 = sol2('W_{f_{total}}')

    sol3 = run_D8_big(True)
    wf3 = sol3('W_{f_{total}}')

    sol4 = run_D8_big(False)
    wf4 = sol4('W_{f_{total}}')

    sol5 = run_D8_big_M08(False)
    wf5 = sol5('W_{f_{total}}')

    sol6 = run_D8_big_M072(False)
    wf6 = sol6('W_{f_{total}}')

    wing_sens = [sol0['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], sol1['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], \
                 sol2['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], sol3['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], \
                 sol4['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], sol5['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], \
                 sol6['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail']]
    HT_sens = [sol0['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], sol1['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], \
                 sol2['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], sol3['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], \
                 sol4['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], sol5['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], \
                 sol6['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail']]
    VT_sens = [sol0['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], sol1['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], \
                 sol2['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], sol3['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], \
                 sol4['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], sol5['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], \
                 sol6['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage']]
    fuse_sens = [sol0['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], sol1['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], \
                 sol2['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], sol3['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], \
                 sol4['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], sol5['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], \
                 sol6['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage']]
    engine_sens = [sol0['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], sol1['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], \
                 sol2['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], sol3['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], \
                 sol4['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], sol5['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], \
                 sol6['sensitivities']['constants']['C_{engsys}_Mission, Aircraft']]
    lg_sens = [sol0['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], sol1['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], \
                 sol2['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], sol3['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], \
                 sol4['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], sol5['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], \
                 sol6['sensitivities']['constants']['C_{lg}_Mission, Aircraft']]

    ytest = [1, wf1/wf0, wf2/wf0, wf3/wf0, wf4/wf0, wf5/wf0, wf6/wf0]
    xtest = [0, 1, 2, 3, 4, 5, 6]
    xlabels = ['Optimized 777-300ER M = 0.84', 'D8 Fuselage, Pi Tail', 'Rear Podded Engines', 'Integrated Engines, BLI', 'Optimize BPR', 'Slow to M = 0.8','Slow to M = 0.72']

    plt.plot(xtest, ytest, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,1.1])
    plt.xlim([-.5, 6.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('$W_{f}/W_{f_0}$', fontsize = 20)
    plt.title('777 Class D8 Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_777_morphing_chart_.pdf', bbox_inches="tight")
    plt.show()


    plt.plot(xtest, wing_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,1.5])
    plt.xlim([-.5, 6.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Wing Weight', fontsize = 20)
    plt.title('Wing Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_777_morphing_chart_wing_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, HT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.08])
    plt.xlim([-.5, 6.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Horizontal Tail Weight', fontsize = 20)
    plt.title('Horizontal Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_777_morphing_chart_HT_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, VT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.45])
    plt.xlim([-.5, 6.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Vertical Tail Weight', fontsize = 20)
    plt.title('Vertical Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_777_morphing_chart_VT_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, fuse_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.8])
    plt.xlim([-.5, 6.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Fuselage Weight', fontsize = 20)
    plt.title('Fuselage Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_777_morphing_chart_fuse_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, engine_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,.6])
    plt.xlim([-.5, 6.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Engine Weight', fontsize = 20)
    plt.title('Engine Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_777_morphing_chart_engine_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, lg_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.15])
    plt.xlim([-.5, 6.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Landing Gear Weight', fontsize = 20)
    plt.title('Landing Gear Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_777_morphing_chart_lg_sens.pdf', bbox_inches="tight")
    plt.show()

def 777_killer_plot_max_optimal():
    #TURN ON THE RIGHT CONSTRIANTS INSIDE D8.py
    sol0 = run_optimal_777(False, True)
    wf0 = sol0('W_{f_{total}}')

    sol1 = run_D8_big_eng_wing(False, True)
    wf1 = sol1('W_{f_{total}}')

    sol2 = run_D8_big_no_BLI(False, True)
    wf2 = sol2('W_{f_{total}}')

    sol3 = run_D8_big(False, True)
    wf3 = sol3('W_{f_{total}}')

    sol4 = run_D8_big_M08(False, True)
    wf4 = sol4('W_{f_{total}}')

    sol5 = run_D8_big_M072(False, True)
    wf5 = sol5('W_{f_{total}}')

    wing_sens = [sol0['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], sol1['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], \
                 sol2['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], sol3['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], \
                 sol4['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], sol5['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing']]
    HT_sens = [sol0['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], sol1['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], \
                 sol2['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], sol3['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], \
                 sol4['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], sol5['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail']]
    VT_sens = [sol0['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], sol1['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], \
                 sol2['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], sol3['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], \
                 sol4['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], sol5['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail']]
    fuse_sens = [sol0['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], sol1['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], \
                 sol2['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], sol3['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], \
                 sol4['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], sol5['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage']]
    engine_sens = [sol0['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], sol1['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], \
                 sol2['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], sol3['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], \
                 sol4['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], sol5['sensitivities']['constants']['C_{engsys}_Mission, Aircraft']]
    lg_sens = [sol0['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], sol1['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], \
                 sol2['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], sol3['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], \
                 sol4['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], sol5['sensitivities']['constants']['C_{lg}_Mission, Aircraft']]

    ytest = [1, wf1/wf0, wf2/wf0, wf3/wf0, wf4/wf0, wf5/wf0]
    xtest = [0, 1, 2, 3, 4, 5]
    xlabels = ['Optimized 777-300ER M = 0.84', 'D8 Fuselage, Pi Tail', 'Rear Podded Engines', 'Integrated Engines, BLI', 'Slow to M = 0.8','Slow to M = 0.72']

    plt.plot(xtest, ytest, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,1.1])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('$W_{f}/W_{f_0}$', fontsize = 20)
    plt.title('777 Class D8 Morphing Chart - Maximum Engine Optimization')
    plt.savefig('Morphing_Chart_Figs/D8_777_morphing_chart_max_optimal_engines.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, wing_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.25])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Wing Weight', fontsize = 20)
    plt.title('Wing Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_777_morphing_chart_max_optimal_engines_wing_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, HT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.02])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Horizontal Tail Weight', fontsize = 20)
    plt.title('Horizontal Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_777_morphing_chart_max_optimal_engines_HT_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, VT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.15])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Vertical Tail Weight', fontsize = 20)
    plt.title('Vertical Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_777_morphing_chart_max_optimal_engines_VT_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, fuse_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.55])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Fuselage Weight', fontsize = 20)
    plt.title('Fuselage Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_777_morphing_chart_max_optimal_engines_fuse_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, engine_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.3])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Engine Weight', fontsize = 20)
    plt.title('Engine Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_777_morphing_chart_max_optimal_engines_engine_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, lg_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.07])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Landing Gear Weight', fontsize = 20)
    plt.title('Landing Gear Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_777_morphing_chart_max_optimal_engines_lg_sens.pdf', bbox_inches="tight")
    plt.show()

def M_08_killer_plot_max_optimal():
    #TURN ON THE RIGHT CONSTRIANTS INSIDE D8.py
    sol0 = run_optimal_737(False, True)
    wf0 = sol0('W_{f_{total}}')

    sol1 = run_M08_D8_eng_wing(False, True)
    wf1 = sol1('W_{f_{total}}')

    sol2 = run_M08_D8_no_BLI(False, True)
    wf2 = sol2('W_{f_{total}}')

    sol3 = run_M08_D8(False, True)
    wf3 = sol3('W_{f_{total}}')

    sol4 = run_optimal_D8(False, True)
    wf4 = sol4('W_{f_{total}}')

    wing_sens = [sol0['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], sol1['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], \
                 sol2['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], sol3['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], \
                 sol4['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing']]
    HT_sens = [sol0['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], sol1['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], \
                 sol2['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], sol3['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], \
                 sol4['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail']]
    VT_sens = [sol0['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], sol1['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], \
                 sol2['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], sol3['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], \
                 sol4['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail']]
    fuse_sens = [sol0['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], sol1['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], \
                 sol2['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], sol3['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], \
                 sol4['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage']]
    engine_sens = [sol0['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], sol1['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], \
                 sol2['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], sol3['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], \
                 sol4['sensitivities']['constants']['C_{engsys}_Mission, Aircraft']]
    lg_sens = [sol0['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], sol1['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], \
                 sol2['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], sol3['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], \
                 sol4['sensitivities']['constants']['C_{lg}_Mission, Aircraft']]

    ytest = [1, wf1/wf0, wf2/wf0, wf3/wf0, wf4/wf0]
    xtest = [0, 1, 2, 3, 4]
    xlabels = ['Optimized 737-800 M = 0.8', 'D8 Fuselage, Pi Tail', 'Rear Podded Engines', 'Integrated Engines, BLI', 'Slow to M = 0.72']

    plt.plot(xtest, ytest, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,1.1])
    plt.xlim([-.5, 4.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('$W_{f}/W_{f_0}$', fontsize = 20)
    plt.title('D8 Morphing Chart - Maximum Engine Optimization')
    plt.savefig('Morphing_Chart_Figs/D8_M08_morphing_chart_max_optimal_engines.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, wing_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.25])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Wing Weight', fontsize = 20)
    plt.title('Wing Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_M08_morphing_chart_max_optimal_engines_wing_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, HT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.02])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Horizontal Tail Weight', fontsize = 20)
    plt.title('Horizontal Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_M08_morphing_chart_max_optimal_engines_HT_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, VT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.15])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Vertical Tail Weight', fontsize = 20)
    plt.title('Vertical Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_M08_morphing_chart_max_optimal_engines_VT_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, fuse_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.55])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Fuselage Weight', fontsize = 20)
    plt.title('Fuselage Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/DD8_M08_morphing_chart_max_optimal_engines_fuse_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, engine_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.3])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Engine Weight', fontsize = 20)
    plt.title('Engine Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_M08_morphing_chart_max_optimal_engines_engine_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, lg_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.07])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Landing Gear Weight', fontsize = 20)
    plt.title('Landing Gear Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_M08_morphing_chart_max_optimal_engines_lg_sens.pdf', bbox_inches="tight")
    plt.show()

def M_08_killer_plot_optimal():
    #invalid due to lack of FPR optimization along with BPR
    sol0 = run_optimal_737(False)
    wf0 = sol0('W_{f_{total}}')

    sol1 = run_M08_D8_eng_wing(False)
    wf1 = sol1('W_{f_{total}}')

    sol2 = run_M08_D8_no_BLI(False)
    wf2 = sol2('W_{f_{total}}')

    sol3 = run_M08_D8(False)
    wf3 = sol3('W_{f_{total}}')

    sol4 = run_optimal_D8(False)
    wf4 = sol4('W_{f_{total}}')

    wing_sens = [sol0['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], sol1['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], \
                 sol2['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], sol3['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], \
                 sol4['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing']]
    HT_sens = [sol0['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], sol1['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], \
                 sol2['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], sol3['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], \
                 sol4['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail']]
    VT_sens = [sol0['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], sol1['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], \
                 sol2['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], sol3['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], \
                 sol4['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail']]
    fuse_sens = [sol0['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], sol1['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], \
                 sol2['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], sol3['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], \
                 sol4['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage']]
    engine_sens = [sol0['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], sol1['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], \
                 sol2['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], sol3['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], \
                 sol4['sensitivities']['constants']['C_{engsys}_Mission, Aircraft']]
    lg_sens = [sol0['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], sol1['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], \
                 sol2['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], sol3['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], \
                 sol4['sensitivities']['constants']['C_{lg}_Mission, Aircraft']]

    ytest = [1, wf1/wf0, wf2/wf0, wf3/wf0, wf4/wf0]
    xtest = [0, 1, 2, 3, 4]
    xlabels = ['Optimized 737-800 M = 0.8', 'D8 Fuselage, Pi Tail', 'Rear Podded Engines', 'Integrated Engines, BLI', 'Slow to M = 0.72']

    plt.plot(xtest, ytest, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,1.1])
    plt.xlim([-.5, 4.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('$W_{f}/W_{f_0}$', fontsize = 20)
    plt.savefig('Morphing_Chart_Figs/D8_M08_morphing_chart_optimal_engines.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, wing_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.25])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Wing Weight', fontsize = 20)
    plt.title('Wing Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_M08_morphing_chart_optimal_engines_wing_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, HT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.02])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Horizontal Tail Weight', fontsize = 20)
    plt.title('Horizontal Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_M08_morphing_chart_optimal_engines_HT_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, VT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.15])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Vertical Tail Weight', fontsize = 20)
    plt.title('Vertical Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_M08_morphing_chart_optimal_engines_VT_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, fuse_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.55])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Fuselage Weight', fontsize = 20)
    plt.title('Fuselage Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/DD8_M08_morphing_chart_optimal_engines_fuse_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, engine_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.3])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Engine Weight', fontsize = 20)
    plt.title('Engine Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_M08_morphing_chart_optimal_engines_engine_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, lg_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.07])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Landing Gear Weight', fontsize = 20)
    plt.title('Landing Gear Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_M08_morphing_chart_optimal_engines_lg_sens.pdf', bbox_inches="tight")
    plt.show()


def M_08_killer_plot_fixed_BPR():
    sol0 = run_optimal_737(True)
    wf0 = sol0('W_{f_{total}}')

    sol1 = run_M08_D8_eng_wing(True)
    wf1 = sol1('W_{f_{total}}')

    sol2 = run_M08_D8_no_BLI(True)
    wf2 = sol2('W_{f_{total}}')

    sol3 = run_M08_D8(True)
    wf3 = sol3('W_{f_{total}}')

    sol4 = run_M08_D8(False)
    wf4 = sol4('W_{f_{total}}')

    sol5 = run_optimal_D8(False)
    wf5 = sol5('W_{f_{total}}')

    wing_sens = [sol0['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], sol1['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], \
                 sol2['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], sol3['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], \
                 sol4['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], sol5['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing']]
    HT_sens = [sol0['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], sol1['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], \
                 sol2['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], sol3['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], \
                 sol4['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], sol5['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail']]
    VT_sens = [sol0['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], sol1['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], \
                 sol2['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], sol3['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], \
                 sol4['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], sol5['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail']]
    fuse_sens = [sol0['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], sol1['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], \
                 sol2['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], sol3['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], \
                 sol4['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], sol5['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage']]
    engine_sens = [sol0['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], sol1['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], \
                 sol2['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], sol3['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], \
                 sol4['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], sol5['sensitivities']['constants']['C_{engsys}_Mission, Aircraft']]
    lg_sens = [sol0['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], sol1['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], \
                 sol2['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], sol3['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], \
                 sol4['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], sol5['sensitivities']['constants']['C_{lg}_Mission, Aircraft']]

    ytest = [1, wf1/wf0, wf2/wf0, wf3/wf0, wf4/wf0, wf5/wf0]
    xtest = [0, 1, 2, 3, 4, 5]
    xlabels = ['Optimized 737-800 M = 0.8', 'D8 Fuselage, Pi Tail', 'Rear Podded Engines', 'Integrated Engines, BLI', 'Optimize BPR', 'Slow to M = 0.72']

    plt.plot(xtest, ytest, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,1.2])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('$W_{f}/W_{f_0}$', fontsize = 20)
    plt.title('D8 Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_M08_morphing_chart_fixed_BPR.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, wing_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.25])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Wing Weight', fontsize = 20)
    plt.title('Wing Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_M08_morphing_chart_fixed_BPR_wing_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, HT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.02])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Horizontal Tail Weight', fontsize = 20)
    plt.title('Horizontal Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_M08_morphing_chart_fixed_BPR_HT_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, VT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.15])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Vertical Tail Weight', fontsize = 20)
    plt.title('Vertical Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_M08_morphing_chart_fixed_BPR_VT_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, fuse_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.55])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Fuselage Weight', fontsize = 20)
    plt.title('Fuselage Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/DD8_M08_morphing_chart_fixed_BPR_fuse_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, engine_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.3])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Engine Weight', fontsize = 20)
    plt.title('Engine Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_M08_morphing_chart_fixed_BPR_engine_sens.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, lg_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.07])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Landing Gear Weight', fontsize = 20)
    plt.title('Landing Gear Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_M08_morphing_chart_fixed_BPR_lg_sens.pdf', bbox_inches="tight")
    plt.show()

def standard_killer_plot():
    sol0 = run_optimal_737(True)
    wf0 = sol0('W_{f_{total}}')

    sol1 = run_M072_737(True)
    wf1 = sol1('W_{f_{total}}')

    sol2 = run_D8_eng_wing(True)
    wf2 = sol2('W_{f_{total}}')

    sol3 = run_D8_no_BLI(True)
    wf3 = sol3('W_{f_{total}}')

    sol4 = run_optimal_D8(True)
    wf4 = sol4('W_{f_{total}}')

    sol5 = run_optimal_D8(False, True)
    wf5 = sol5('W_{f_{total}}')

    wing_sens = [sol0['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], sol1['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], \
                 sol2['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], sol3['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], \
                 sol4['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], sol5['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing']]
    HT_sens = [sol0['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], sol1['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], \
                 sol2['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], sol3['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], \
                 sol4['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], sol5['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail']]
    VT_sens = [sol0['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], sol1['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], \
                 sol2['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], sol3['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], \
                 sol4['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], sol5['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail']]
    fuse_sens = [sol0['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], sol1['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], \
                 sol2['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], sol3['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], \
                 sol4['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], sol5['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage']]
    engine_sens = [sol0['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], sol1['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], \
                 sol2['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], sol3['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], \
                 sol4['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], sol5['sensitivities']['constants']['C_{engsys}_Mission, Aircraft']]
    lg_sens = [sol0['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], sol1['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], \
                 sol2['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], sol3['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], \
                 sol4['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], sol5['sensitivities']['constants']['C_{lg}_Mission, Aircraft']]



    ytest = [1, wf1/wf0, wf2/wf0, wf3/wf0, wf4/wf0, wf5/wf0, wf5/wf0]
    xtest = [0, 1, 2, 3, 4, 5, 6]
    xlabels = ['Optimized 737-800 M = 0.8', 'Slow to M = 0.72', 'D8 fuselage, Pi tail', 'Rear podded engines', 'Integrated engines, BLI = D8', 'Optimize engine', '2010 Engines']

    plt.plot(xtest, ytest, "o--")
    plt.plot([0, 1, 2, 3, 4, 5, 6], [1, .88, .81, .82, .67, .66, .63], "o--")
    plt.plot([0, 1, 2, 3, 6], [1, .868, .871, .865, .602], "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,1.1])
    plt.xlim([-.5, 6.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('$W_{f}/W_{f_0}$', fontsize = 20)
    plt.title('D8 Morphing Chart')
    plt.legend(['SP Model', 'TASOPT', 'NASA'], loc=3)
    plt.savefig('Morphing_Chart_Figs/D8_standard_killer_chart.pdf', bbox_inches="tight")
    plt.show()

    xtest = [0, 1, 2, 3, 4, 5]

    plt.plot(xtest, wing_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.25])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Wing Weight', fontsize = 20)
    plt.title('Wing Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_standard_wing_sens_killer_chart.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, HT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.02])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Horizontal Tail Weight', fontsize = 20)
    plt.title('Horizontal Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_standard_HT_sens_killer_chart.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, VT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.15])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Vertical Tail Weight', fontsize = 20)
    plt.title('Vertical Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_standard_VT_sens_killer_chart.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, fuse_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.55])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Fuselage Weight', fontsize = 20)
    plt.title('Fuselage Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_standard_fuse_sens_killer_chart.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, engine_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.3])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Engine Weight', fontsize = 20)
    plt.title('Engine Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_standard_engine_sens_killer_chart.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, lg_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.07])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Landing Gear Weight', fontsize = 20)
    plt.title('Landing Gear Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_standard_lg_sens_killer_chart.pdf', bbox_inches="tight")
    plt.show()

def standard_killer_plot_optimal_engine():
    #invalid due to lack of FPR optimization along with BPR
    sol0 = run_optimal_737(False, False)
    wf0 = sol0('W_{f_{total}}')

    sol1 = run_M072_737(False, False)
    wf1 = sol1('W_{f_{total}}')

    sol2 = run_D8_eng_wing(False, False)
    wf2 = sol2('W_{f_{total}}')

    sol3 = run_D8_no_BLI(False, False)
    wf3 = sol3('W_{f_{total}}')

    sol4 = run_optimal_D8(False, False)
    wf4 = sol4('W_{f_{total}}')

    wing_sens = [sol0['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], sol1['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], \
                 sol2['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], sol3['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], \
                 sol4['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing']]
    HT_sens = [sol0['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], sol1['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], \
                 sol2['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], sol3['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], \
                 sol4['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail']]
    VT_sens = [sol0['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], sol1['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], \
                 sol2['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], sol3['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], \
                 sol4['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail']]
    fuse_sens = [sol0['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], sol1['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], \
                 sol2['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], sol3['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], \
                 sol4['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage']]
    engine_sens = [sol0['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], sol1['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], \
                 sol2['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], sol3['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], \
                 sol4['sensitivities']['constants']['C_{engsys}_Mission, Aircraft']]
    lg_sens = [sol0['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], sol1['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], \
                 sol2['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], sol3['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], \
                 sol4['sensitivities']['constants']['C_{lg}_Mission, Aircraft']]

    ytest = [1, wf1/wf0, wf2/wf0, wf3/wf0, wf4/wf0]
    xtest = [0, 1, 2, 3, 4]
    xlabels = ['Optimized 737-800 M = 0.8', 'Slow to M = 0.72', 'D8 fuselage, Pi tail', 'Rear podded engines', 'Integrated engines, BLI = D8']

    plt.plot(xtest, ytest, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,1.1])
    plt.xlim([-.5, 4.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('$W_{f}/W_{f_0}$', fontsize = 20)
    plt.title('D8 Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_standard_killer_chart_optimal_engines.pdf', bbox_inches="tight")
    plt.show()

def standard_killer_plot_max_optimal_engine():
    #TURN ON THE RIGHT CONSTRIANTS INSIDE D8.py
    sol0 = run_optimal_737(False, True)
    wf0 = sol0('W_{f_{total}}')

    sol1 = run_M072_737(False, True)
    wf1 = sol1('W_{f_{total}}')

    sol2 = run_D8_eng_wing(False, True)
    wf2 = sol2('W_{f_{total}}')

    sol3 = run_D8_no_BLI(False, True)
    wf3 = sol3('W_{f_{total}}')

    sol4 = run_optimal_D8(False, True)
    wf4 = sol4('W_{f_{total}}')

    wing_sens = [sol0['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], sol1['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], \
                 sol2['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], sol3['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing'], \
                 sol4['sensitivities']['constants']['C_{wing}_Mission, Aircraft, Wing']]
    HT_sens = [sol0['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], sol1['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], \
                 sol2['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], sol3['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail'], \
                 sol4['sensitivities']['constants']['C_{HT}_Mission, Aircraft, HorizontalTail']]
    VT_sens = [sol0['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], sol1['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], \
                 sol2['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], sol3['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail'], \
                 sol4['sensitivities']['constants']['C_{VT}_Mission, Aircraft, VerticalTail']]
    fuse_sens = [sol0['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], sol1['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], \
                 sol2['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], sol3['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage'], \
                 sol4['sensitivities']['constants']['C_{fuse}_Mission, Aircraft, Fuselage']]
    engine_sens = [sol0['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], sol1['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], \
                 sol2['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], sol3['sensitivities']['constants']['C_{engsys}_Mission, Aircraft'], \
                 sol4['sensitivities']['constants']['C_{engsys}_Mission, Aircraft']]
    lg_sens = [sol0['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], sol1['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], \
                 sol2['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], sol3['sensitivities']['constants']['C_{lg}_Mission, Aircraft'], \
                 sol4['sensitivities']['constants']['C_{lg}_Mission, Aircraft']]

    ytest = [1, wf1/wf0, wf2/wf0, wf3/wf0, wf4/wf0]
    xtest = [0, 1, 2, 3, 4]
    xlabels = ['Optimized 737-800 M = 0.8', 'Slow to M = 0.72', 'D8 fuselage, Pi tail', 'Rear podded engines', 'Integrated engines, BLI = D8']

    plt.plot(xtest, ytest, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,1.1])
    plt.xlim([-.5, 4.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('$W_{f}/W_{f_0}$', fontsize = 20)
    plt.title('D8 Morphing Chart - Maximum Engine Optimization')
    plt.savefig('Morphing_Chart_Figs/D8_standard_killer_chart_max_optimal_engines.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, wing_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.25])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Wing Weight', fontsize = 20)
    plt.title('Wing Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_standard_max_opt_wing_sens_killer_chart.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, HT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.02])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Horizontal Tail Weight', fontsize = 20)
    plt.title('Horizontal Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_standard_max_opt_HT_sens_killer_chart.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, VT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.15])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Vertical Tail Weight', fontsize = 20)
    plt.title('Vertical Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_standard_max_opt_VT_sens_killer_chart.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, fuse_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.55])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Fuselage Weight', fontsize = 20)
    plt.title('Fuselage Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_standard_max_opt_fuse_sens_killer_chart.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, engine_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.3])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Engine Weight', fontsize = 20)
    plt.title('Engine Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_standard_max_opt_engine_sens_killer_chart.pdf', bbox_inches="tight")
    plt.show()

    plt.plot(xtest, lg_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.07])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Landing Gear Weight', fontsize = 20)
    plt.title('Landing Gear Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_standard_max_opt_lg_sens_killer_chart.pdf', bbox_inches="tight")
    plt.show()
