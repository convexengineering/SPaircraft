# only needed for plotting
import matplotlib.pyplot as plt

# Solution saving
from saveSol import genSolOut

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
from SP_aircraft import run_M072_optimal_777
from SP_aircraft import run_M08_optimal_777
from SP_aircraft import run_M072_D8_big_eng_wing
from SP_aircraft import run_M072_D8_big_no_BLI

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

    wing_sens = [sol0['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], sol1['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], \
                 sol2['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], sol3['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], \
                 sol4['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], sol5['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing']]
    HT_sens = [sol0['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], sol1['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], \
                 sol2['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], sol3['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], \
                 sol4['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], sol5['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail']]
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

    plt.plot(xtest, ytest, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,1.1])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('$W_{f}/W_{f_0}$', fontsize = 20)
    plt.title('RJ Class D8 Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_M08_standard_killer_chart.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, wing_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.2])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Wing Weight', fontsize = 20)
    plt.title('Wing Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_M08_standard_killer_chart_wing_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, HT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.015])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Horizontal Tail Weight', fontsize = 20)
    plt.title('Horizontal Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_M08_standard_killer_chart_HT_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, VT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.15])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Vertical Tail Weight', fontsize = 20)
    plt.title('Vertical Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_M08_standard_killer_chart_VT_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, fuse_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.55])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Fuselage Weight', fontsize = 20)
    plt.title('Fuselage Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_M08_standard_killer_chart_fuse_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, engine_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.3])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Engine Weight', fontsize = 20)
    plt.title('Engine Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_M08_standard_killer_chart_engine_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, lg_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.07])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Landing Gear Weight', fontsize = 20)
    plt.title('Landing Gear Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_M08_standard_killer_chart_lg_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

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

    wing_sens = [sol0['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], sol1['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], \
                 sol2['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], sol3['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], \
                 sol4['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], sol5['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing']]
    HT_sens = [sol0['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], sol1['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], \
                 sol2['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], sol3['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], \
                 sol4['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], sol5['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail']]
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
    plt.show(), plt.close()

    plt.plot(xtest, wing_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.25])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Wing Weight', fontsize = 20)
    plt.title('Wing Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_standard_killer_chart_wing_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, HT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.02])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Horizontal Tail Weight', fontsize = 20)
    plt.title('Horizontal Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_standard_killer_chart_HT_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, VT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.1])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Vertical Tail Weight', fontsize = 20)
    plt.title('Vertical Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_standard_killer_chart_VT_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, fuse_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.55])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Fuselage Weight', fontsize = 20)
    plt.title('Fuselage Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_standard_killer_chart_fuse_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, engine_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.25])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Engine Weight', fontsize = 20)
    plt.title('Engine Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_standard_killer_chart_engine_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, lg_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.07])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Landing Gear Weight', fontsize = 20)
    plt.title('Landing Gear Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_standard_killer_chart_lg_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

def RJ_killer_plot_max_opt_eng():
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

    wing_sens = [sol0['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], sol1['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], \
                 sol2['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], sol3['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], \
                 sol4['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing']]
    HT_sens = [sol0['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], sol1['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], \
                 sol2['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], sol3['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], \
                 sol4['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail']]
    VT_sens = [sol0['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail'], sol1['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail'], \
                 sol2['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail'], sol3['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail'], \
                 sol4['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail']]
    fuse_sens = [sol0['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage'], sol1['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage'], \
                 sol2['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage'], sol3['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage'], \
                 sol4['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage']]
    engine_sens = [sol0['sensitivities']['constants']['C_{engsys}_Mission/Aircraft'], sol1['sensitivities']['constants']['C_{engsys}_Mission/Aircraft'], \
                 sol2['sensitivities']['constants']['C_{engsys}_Mission/Aircraft'], sol3['sensitivities']['constants']['C_{engsys}_Mission/Aircraft'], \
                 sol4['sensitivities']['constants']['C_{engsys}_Mission/Aircraft']]
    lg_sens = [sol0['sensitivities']['constants']['C_{lg}_Mission/Aircraft'], sol1['sensitivities']['constants']['C_{lg}_Mission/Aircraft'], \
                 sol2['sensitivities']['constants']['C_{lg}_Mission/Aircraft'], sol3['sensitivities']['constants']['C_{lg}_Mission/Aircraft'], \
                 sol4['sensitivities']['constants']['C_{lg}_Mission/Aircraft']]

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
    plt.savefig('Morphing_Chart_Figs/RJ_killer_chart_max_opt_eng.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, wing_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.25])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Wing Weight', fontsize = 20)
    plt.title('Wing Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_killer_chart_max_opt_eng_wing_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, HT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.02])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Horizontal Tail Weight', fontsize = 20)
    plt.title('Horizontal Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_killer_chart_max_opt_eng_HT_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, VT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.15])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Vertical Tail Weight', fontsize = 20)
    plt.title('Vertical Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_killer_chart_max_opt_eng_VT_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, fuse_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.55])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Fuselage Weight', fontsize = 20)
    plt.title('Fuselage Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_killer_chart_max_opt_eng_fuse_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, engine_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.3])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Engine Weight', fontsize = 20)
    plt.title('Engine Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_killer_chart_max_opt_eng_engine_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, lg_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.07])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Landing Gear Weight', fontsize = 20)
    plt.title('Landing Gear Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/RJ_killer_chart_max_opt_eng_lg_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

def b777_killer_plot_standard():
    sol0 = run_optimal_777(True)
    wf0 = sol0('W_{f_{total}}')

    sol1 = run_M08_optimal_777(True)
    wf1 = sol1('W_{f_{total}}')

    sol2 = run_M072_optimal_777(True)
    wf2 = sol2('W_{f_{total}}')

    sol3 = run_M072_D8_big_eng_wing(True)
    wf3 = sol3('W_{f_{total}}')

    sol4 = run_M072_D8_big_no_BLI(True)
    wf4 = sol4('W_{f_{total}}')

    sol5 = run_D8_big_M072(True)
    wf5 = sol5('W_{f_{total}}')

    sol6 = run_D8_big_M072(False, True)
    wf6 = sol6('W_{f_{total}}')

    wing_sens = [sol0['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], sol1['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], \
                 sol2['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], sol3['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], \
                 sol4['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], sol5['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], \
                 sol6['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail']]
    HT_sens = [sol0['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], sol1['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], \
                 sol2['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], sol3['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], \
                 sol4['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], sol5['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], \
                 sol6['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail']]
    VT_sens = [sol0['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail'], sol1['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail'], \
                 sol2['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail'], sol3['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail'], \
                 sol4['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail'], sol5['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail'], \
                 sol6['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage']]
    fuse_sens = [sol0['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage'], sol1['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage'], \
                 sol2['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage'], sol3['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage'], \
                 sol4['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage'], sol5['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage'], \
                 sol6['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage']]
    engine_sens = [sol0['sensitivities']['constants']['C_{engsys}_Mission/Aircraft'], sol1['sensitivities']['constants']['C_{engsys}_Mission/Aircraft'], \
                 sol2['sensitivities']['constants']['C_{engsys}_Mission/Aircraft'], sol3['sensitivities']['constants']['C_{engsys}_Mission/Aircraft'], \
                 sol4['sensitivities']['constants']['C_{engsys}_Mission/Aircraft'], sol5['sensitivities']['constants']['C_{engsys}_Mission/Aircraft'], \
                 sol6['sensitivities']['constants']['C_{engsys}_Mission/Aircraft']]
    lg_sens = [sol0['sensitivities']['constants']['C_{lg}_Mission/Aircraft'], sol1['sensitivities']['constants']['C_{lg}_Mission/Aircraft'], \
                 sol2['sensitivities']['constants']['C_{lg}_Mission/Aircraft'], sol3['sensitivities']['constants']['C_{lg}_Mission/Aircraft'], \
                 sol4['sensitivities']['constants']['C_{lg}_Mission/Aircraft'], sol5['sensitivities']['constants']['C_{lg}_Mission/Aircraft'], \
                 sol6['sensitivities']['constants']['C_{lg}_Mission/Aircraft']]

    ytest = [1, wf1/wf0, wf2/wf0, wf3/wf0, wf4/wf0, wf5/wf0, wf6/wf0]
    xtest = [0, 1, 2, 3, 4, 5, 6]
    xlabels = ['Optimized 777-300ER M = 0.84', 'Slow to M = 0.8','Slow to M = 0.72', 'D8 Fuselage, Pi Tail', 'Rear Podded Engines', 'Integrated Engines, BLI', 'Optimize Engine']

    plt.plot(xtest, ytest, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,1.1])
    plt.xlim([-.5, 6.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('$W_{f}/W_{f_0}$', fontsize = 20)
    plt.title('777 Class D8 Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_777_morphing_chart_.pdf', bbox_inches="tight")
    plt.show(), plt.close()


    plt.plot(xtest, wing_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,1.5])
    plt.xlim([-.5, 6.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Wing Weight', fontsize = 20)
    plt.title('Wing Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_777_morphing_chart_wing_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, HT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.08])
    plt.xlim([-.5, 6.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Horizontal Tail Weight', fontsize = 20)
    plt.title('Horizontal Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_777_morphing_chart_HT_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, VT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.45])
    plt.xlim([-.5, 6.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Vertical Tail Weight', fontsize = 20)
    plt.title('Vertical Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_777_morphing_chart_VT_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, fuse_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.8])
    plt.xlim([-.5, 6.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Fuselage Weight', fontsize = 20)
    plt.title('Fuselage Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_777_morphing_chart_fuse_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, engine_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,.6])
    plt.xlim([-.5, 6.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Engine Weight', fontsize = 20)
    plt.title('Engine Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_777_morphing_chart_engine_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, lg_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.15])
    plt.xlim([-.5, 6.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Landing Gear Weight', fontsize = 20)
    plt.title('Landing Gear Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_777_morphing_chart_lg_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

def b777_killer_plot_reordered():
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

    wing_sens = [sol0['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], sol1['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], \
                 sol2['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], sol3['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], \
                 sol4['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], sol5['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], \
                 sol6['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail']]
    HT_sens = [sol0['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], sol1['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], \
                 sol2['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], sol3['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], \
                 sol4['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], sol5['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], \
                 sol6['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail']]
    VT_sens = [sol0['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail'], sol1['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail'], \
                 sol2['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail'], sol3['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail'], \
                 sol4['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail'], sol5['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail'], \
                 sol6['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage']]
    fuse_sens = [sol0['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage'], sol1['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage'], \
                 sol2['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage'], sol3['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage'], \
                 sol4['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage'], sol5['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage'], \
                 sol6['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage']]
    engine_sens = [sol0['sensitivities']['constants']['C_{engsys}_Mission/Aircraft'], sol1['sensitivities']['constants']['C_{engsys}_Mission/Aircraft'], \
                 sol2['sensitivities']['constants']['C_{engsys}_Mission/Aircraft'], sol3['sensitivities']['constants']['C_{engsys}_Mission/Aircraft'], \
                 sol4['sensitivities']['constants']['C_{engsys}_Mission/Aircraft'], sol5['sensitivities']['constants']['C_{engsys}_Mission/Aircraft'], \
                 sol6['sensitivities']['constants']['C_{engsys}_Mission/Aircraft']]
    lg_sens = [sol0['sensitivities']['constants']['C_{lg}_Mission/Aircraft'], sol1['sensitivities']['constants']['C_{lg}_Mission/Aircraft'], \
                 sol2['sensitivities']['constants']['C_{lg}_Mission/Aircraft'], sol3['sensitivities']['constants']['C_{lg}_Mission/Aircraft'], \
                 sol4['sensitivities']['constants']['C_{lg}_Mission/Aircraft'], sol5['sensitivities']['constants']['C_{lg}_Mission/Aircraft'], \
                 sol6['sensitivities']['constants']['C_{lg}_Mission/Aircraft']]

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
    plt.savefig('Morphing_Chart_Figs/D8_777_morphing_chart_reordered.pdf', bbox_inches="tight")
    plt.show(), plt.close()


    plt.plot(xtest, wing_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,1.5])
    plt.xlim([-.5, 6.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Wing Weight', fontsize = 20)
    plt.title('Wing Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_777_morphing_chart_reordered_wing_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, HT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.08])
    plt.xlim([-.5, 6.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Horizontal Tail Weight', fontsize = 20)
    plt.title('Horizontal Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_777_morphing_chart_reordered_HT_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, VT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.7])
    plt.xlim([-.5, 6.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Vertical Tail Weight', fontsize = 20)
    plt.title('Vertical Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_777_morphing_chart_reordered_VT_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, fuse_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,1])
    plt.xlim([-.5, 6.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Fuselage Weight', fontsize = 20)
    plt.title('Fuselage Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_777_morphing_chart_reordered_fuse_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, engine_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,.6])
    plt.xlim([-.5, 6.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Engine Weight', fontsize = 20)
    plt.title('Engine Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_777_morphing_chart_reordered_engine_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, lg_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.11])
    plt.xlim([-.5, 6.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Landing Gear Weight', fontsize = 20)
    plt.title('Landing Gear Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_777_morphing_chart_reordered_lg_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

def b777_killer_plot_max_opt_eng():
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

    wing_sens = [sol0['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], sol1['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], \
                 sol2['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], sol3['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], \
                 sol4['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], sol5['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing']]
    HT_sens = [sol0['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], sol1['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], \
                 sol2['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], sol3['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], \
                 sol4['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], sol5['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail']]
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
    plt.savefig('Morphing_Chart_Figs/D8_777_morphing_chart_max_opt_eng.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, wing_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.25])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Wing Weight', fontsize = 20)
    plt.title('Wing Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_777_morphing_chart_max_opt_eng_wing_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, HT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.02])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Horizontal Tail Weight', fontsize = 20)
    plt.title('Horizontal Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_777_morphing_chart_max_opt_eng_HT_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, VT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.15])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Vertical Tail Weight', fontsize = 20)
    plt.title('Vertical Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_777_morphing_chart_max_opt_eng_VT_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, fuse_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.55])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Fuselage Weight', fontsize = 20)
    plt.title('Fuselage Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_777_morphing_chart_max_opt_eng_fuse_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, engine_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.3])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Engine Weight', fontsize = 20)
    plt.title('Engine Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_777_morphing_chart_max_opt_eng_engine_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, lg_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.07])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Landing Gear Weight', fontsize = 20)
    plt.title('Landing Gear Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_777_morphing_chart_max_opt_eng_lg_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

def M_08_killer_plot_max_opt_eng():
    #TURN ON THE RIGHT CONSTRAINTS INSIDE D8.py
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

    wing_sens = [sol0['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], sol1['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], \
                 sol2['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], sol3['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], \
                 sol4['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing']]
    HT_sens = [sol0['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], sol1['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], \
                 sol2['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], sol3['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], \
                 sol4['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail']]
    VT_sens = [sol0['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail'], sol1['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail'], \
                 sol2['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail'], sol3['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail'], \
                 sol4['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail']]
    fuse_sens = [sol0['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage'], sol1['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage'], \
                 sol2['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage'], sol3['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage'], \
                 sol4['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage']]
    engine_sens = [sol0['sensitivities']['constants']['C_{engsys}_Mission/Aircraft'], sol1['sensitivities']['constants']['C_{engsys}_Mission/Aircraft'], \
                 sol2['sensitivities']['constants']['C_{engsys}_Mission/Aircraft'], sol3['sensitivities']['constants']['C_{engsys}_Mission/Aircraft'], \
                 sol4['sensitivities']['constants']['C_{engsys}_Mission/Aircraft']]
    lg_sens = [sol0['sensitivities']['constants']['C_{lg}_Mission/Aircraft'], sol1['sensitivities']['constants']['C_{lg}_Mission/Aircraft'], \
                 sol2['sensitivities']['constants']['C_{lg}_Mission/Aircraft'], sol3['sensitivities']['constants']['C_{lg}_Mission/Aircraft'], \
                 sol4['sensitivities']['constants']['C_{lg}_Mission/Aircraft']]

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
    plt.savefig('Morphing_Chart_Figs/D8_M08_morphing_chart_max_opt_eng.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, wing_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.25])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Wing Weight', fontsize = 20)
    plt.title('Wing Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_M08_morphing_chart_max_opt_eng_wing_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, HT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.02])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Horizontal Tail Weight', fontsize = 20)
    plt.title('Horizontal Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_M08_morphing_chart_max_opt_eng_HT_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, VT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.15])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Vertical Tail Weight', fontsize = 20)
    plt.title('Vertical Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_M08_morphing_chart_max_opt_eng_VT_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, fuse_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.55])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Fuselage Weight', fontsize = 20)
    plt.title('Fuselage Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/DD8_M08_morphing_chart_max_opt_eng_fuse_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, engine_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.3])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Engine Weight', fontsize = 20)
    plt.title('Engine Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_M08_morphing_chart_max_opt_eng_engine_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, lg_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.07])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Landing Gear Weight', fontsize = 20)
    plt.title('Landing Gear Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_M08_morphing_chart_max_opt_eng_lg_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

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

    wing_sens = [sol0['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], sol1['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], \
                 sol2['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], sol3['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], \
                 sol4['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing']]
    HT_sens = [sol0['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], sol1['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], \
                 sol2['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], sol3['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], \
                 sol4['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail']]
    VT_sens = [sol0['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail'], sol1['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail'], \
                 sol2['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail'], sol3['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail'], \
                 sol4['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail']]
    fuse_sens = [sol0['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage'], sol1['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage'], \
                 sol2['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage'], sol3['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage'], \
                 sol4['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage']]
    engine_sens = [sol0['sensitivities']['constants']['C_{engsys}_Mission/Aircraft'], sol1['sensitivities']['constants']['C_{engsys}_Mission/Aircraft'], \
                 sol2['sensitivities']['constants']['C_{engsys}_Mission/Aircraft'], sol3['sensitivities']['constants']['C_{engsys}_Mission/Aircraft'], \
                 sol4['sensitivities']['constants']['C_{engsys}_Mission/Aircraft']]
    lg_sens = [sol0['sensitivities']['constants']['C_{lg}_Mission/Aircraft'], sol1['sensitivities']['constants']['C_{lg}_Mission/Aircraft'], \
                 sol2['sensitivities']['constants']['C_{lg}_Mission/Aircraft'], sol3['sensitivities']['constants']['C_{lg}_Mission/Aircraft'], \
                 sol4['sensitivities']['constants']['C_{lg}_Mission/Aircraft']]

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
    plt.show(), plt.close()

    plt.plot(xtest, wing_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.25])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Wing Weight', fontsize = 20)
    plt.title('Wing Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_M08_morphing_chart_optimal_engines_wing_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, HT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.02])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Horizontal Tail Weight', fontsize = 20)
    plt.title('Horizontal Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_M08_morphing_chart_optimal_engines_HT_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, VT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.15])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Vertical Tail Weight', fontsize = 20)
    plt.title('Vertical Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_M08_morphing_chart_optimal_engines_VT_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, fuse_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.55])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Fuselage Weight', fontsize = 20)
    plt.title('Fuselage Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/DD8_M08_morphing_chart_optimal_engines_fuse_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, engine_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.3])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Engine Weight', fontsize = 20)
    plt.title('Engine Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_M08_morphing_chart_optimal_engines_engine_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, lg_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.07])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Landing Gear Weight', fontsize = 20)
    plt.title('Landing Gear Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_M08_morphing_chart_optimal_engines_lg_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()


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

    wing_sens = [sol0['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], sol1['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], \
                 sol2['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], sol3['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], \
                 sol4['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], sol5['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing']]
    HT_sens = [sol0['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], sol1['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], \
                 sol2['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], sol3['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], \
                 sol4['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], sol5['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail']]
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
    plt.show(), plt.close()

    plt.plot(xtest, wing_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.25])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Wing Weight', fontsize = 20)
    plt.title('Wing Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_M08_morphing_chart_fixed_BPR_wing_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, HT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.02])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Horizontal Tail Weight', fontsize = 20)
    plt.title('Horizontal Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_M08_morphing_chart_fixed_BPR_HT_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, VT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.15])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Vertical Tail Weight', fontsize = 20)
    plt.title('Vertical Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_M08_morphing_chart_fixed_BPR_VT_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, fuse_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.55])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Fuselage Weight', fontsize = 20)
    plt.title('Fuselage Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/DD8_M08_morphing_chart_fixed_BPR_fuse_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, engine_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.3])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Engine Weight', fontsize = 20)
    plt.title('Engine Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_M08_morphing_chart_fixed_BPR_engine_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, lg_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.07])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Landing Gear Weight', fontsize = 20)
    plt.title('Landing Gear Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_M08_morphing_chart_fixed_BPR_lg_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

def standard_killer_plot():
    sol0 = run_optimal_737(True)
    wf0 = sol0('W_{f_{total}}')
    genSolOut(sol0.table(),0)

    sol1 = run_M072_737(True)
    wf1 = sol1('W_{f_{total}}')
    genSolOut(sol1.table(),1)

    sol2 = run_D8_eng_wing(True)
    wf2 = sol2('W_{f_{total}}')
    genSolOut(sol2.table(),2)

    sol3 = run_D8_no_BLI(True)
    wf3 = sol3('W_{f_{total}}')
    genSolOut(sol3.table(),3)

    sol4 = run_optimal_D8(True)
    wf4 = sol4('W_{f_{total}}')
    genSolOut(sol4.table(),4)

    sol5 = run_optimal_D8(False, True)
    wf5 = sol5('W_{f_{total}}')
    genSolOut(sol5.table(),5)

    wing_sens = [sol0['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], sol1['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], \
                 sol2['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], sol3['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], \
                 sol4['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], sol5['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing']]
    HT_sens = [sol0['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], sol1['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], \
                 sol2['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], sol3['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], \
                 sol4['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], sol5['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail']]
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
    CrTt41_sens = [sol0['sensitivities']['constants']['CruiseTt41max_Mission'], sol1['sensitivities']['constants']['CruiseTt41max_Mission'], \
             sol2['sensitivities']['constants']['CruiseTt41max_Mission'], sol3['sensitivities']['constants']['CruiseTt41max_Mission'], \
             sol4['sensitivities']['constants']['CruiseTt41max_Mission'], sol5['sensitivities']['constants']['CruiseTt41max_Mission']]
    CruiseAltBegin = [sol0('h_Mission/CruiseSegment/FlightState/Altitude')[0].to('ft').magnitude,sol1('h_Mission/CruiseSegment/FlightState/Altitude')[0].to('ft').magnitude,
                 sol2('h_Mission/CruiseSegment/FlightState/Altitude')[0].to('ft').magnitude,sol3('h_Mission/CruiseSegment/FlightState/Altitude')[0].to('ft').magnitude,
                 sol4('h_Mission/CruiseSegment/FlightState/Altitude')[0].to('ft').magnitude,sol5('h_Mission/CruiseSegment/FlightState/Altitude')[0].to('ft').magnitude]
    CruiseAltEnd = [sol0('h_Mission/CruiseSegment/FlightState/Altitude')[-1].to('ft').magnitude,sol1('h_Mission/CruiseSegment/FlightState/Altitude')[-1].to('ft').magnitude,
                 sol2('h_Mission/CruiseSegment/FlightState/Altitude')[-1].to('ft').magnitude,sol3('h_Mission/CruiseSegment/FlightState/Altitude')[-1].to('ft').magnitude,
                 sol4('h_Mission/CruiseSegment/FlightState/Altitude')[-1].to('ft').magnitude,sol5('h_Mission/CruiseSegment/FlightState/Altitude')[-1].to('ft').magnitude]


    ytest = [1, wf1/wf0, wf2/wf0, wf3/wf0, wf4/wf0, wf5/wf0, wf5/wf0]
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

    plt.plot(xtest, CruiseAltBegin, "o--")
    plt.plot(xtest, CruiseAltEnd, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.xlim([-.5, 6.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Cruise Altitude (ft)', fontsize = 20)
    plt.title('Cruise Altitude Morphing Chart')
    plt.legend(['Beginning', 'End'])
    plt.savefig('Morphing_Chart_Figs/D8_standard_morphing_chart_CruiseAlt.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, wing_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.25])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Wing Weight', fontsize = 20)
    plt.title('Wing Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_standard_morphing_chart_wing_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, HT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.02])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Horizontal Tail Weight', fontsize = 20)
    plt.title('Horizontal Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_standard_morphing_chart_HT_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, VT_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.15])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Vertical Tail Weight', fontsize = 20)
    plt.title('Vertical Tail Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_standard_morphing_chart_VT_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close()

    plt.plot(xtest, fuse_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.55])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Fuselage Weight', fontsize = 20)
    plt.title('Fuselage Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_standard_morphing_chart_fuse_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close(), plt.close()

    plt.plot(xtest, engine_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.3])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Engine Weight', fontsize = 20)
    plt.title('Engine Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_standard_morphing_chart_engine_sens.pdf', bbox_inches="tight")
    plt.show, plt.close()

    plt.plot(xtest, lg_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,0.12])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Landing Gear Weight', fontsize = 20)
    plt.title('Landing Gear Weight Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_standard_morphing_chart_lg_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close(),

    plt.plot(xtest, Mmin_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0, 1.1])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Minimum Cruise Mach Number')
    plt.title('Minimum Cruise Mach Number Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_standard_morphing_chart_Mmin_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close(),

    plt.plot(xtest, CrTt41_sens, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([-1, 0.1])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('Sensitivity to Max Cruise $T_{t_{4.1}}$')
    plt.title('Max Cruise $T_{t_{4.1}}$ Sensitivity Morphing Chart')
    plt.savefig('Morphing_Chart_Figs/D8_standard_morphing_chart_CrTt41max_sens.pdf', bbox_inches="tight")
    plt.show(), plt.close(),

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

    wing_sens = [sol0['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], sol1['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], \
                 sol2['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], sol3['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], \
                 sol4['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing']]
    HT_sens = [sol0['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], sol1['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], \
                 sol2['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], sol3['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], \
                 sol4['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail']]
    VT_sens = [sol0['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail'], sol1['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail'], \
                 sol2['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail'], sol3['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail'], \
                 sol4['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail']]
    fuse_sens = [sol0['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage'], sol1['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage'], \
                 sol2['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage'], sol3['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage'], \
                 sol4['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage']]
    engine_sens = [sol0['sensitivities']['constants']['C_{engsys}_Mission/Aircraft'], sol1['sensitivities']['constants']['C_{engsys}_Mission/Aircraft'], \
                 sol2['sensitivities']['constants']['C_{engsys}_Mission/Aircraft'], sol3['sensitivities']['constants']['C_{engsys}_Mission/Aircraft'], \
                 sol4['sensitivities']['constants']['C_{engsys}_Mission/Aircraft']]
    lg_sens = [sol0['sensitivities']['constants']['C_{lg}_Mission/Aircraft'], sol1['sensitivities']['constants']['C_{lg}_Mission/Aircraft'], \
                 sol2['sensitivities']['constants']['C_{lg}_Mission/Aircraft'], sol3['sensitivities']['constants']['C_{lg}_Mission/Aircraft'], \
                 sol4['sensitivities']['constants']['C_{lg}_Mission/Aircraft']]

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
    plt.show(), plt.close()

def standard_killer_plot_max_opt_eng():
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

    wing_sens = [sol0['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], sol1['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], \
                 sol2['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], sol3['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing'], \
                 sol4['sensitivities']['constants']['C_{wing}_Mission/Aircraft/Wing']]
    HT_sens = [sol0['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], sol1['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], \
                 sol2['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], sol3['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail'], \
                 sol4['sensitivities']['constants']['C_{HT}_Mission/Aircraft/HorizontalTail']]
    VT_sens = [sol0['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail'], sol1['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail'], \
                 sol2['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail'], sol3['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail'], \
                 sol4['sensitivities']['constants']['C_{VT}_Mission/Aircraft/VerticalTail']]
    fuse_sens = [sol0['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage'], sol1['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage'], \
                 sol2['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage'], sol3['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage'], \
                 sol4['sensitivities']['constants']['C_{fuse}_Mission/Aircraft/Fuselage']]
    engine_sens = [sol0['sensitivities']['constants']['C_{engsys}_Mission/Aircraft'], sol1['sensitivities']['constants']['C_{engsys}_Mission/Aircraft'], \
                 sol2['sensitivities']['constants']['C_{engsys}_Mission/Aircraft'], sol3['sensitivities']['constants']['C_{engsys}_Mission/Aircraft'], \
                 sol4['sensitivities']['constants']['C_{engsys}_Mission/Aircraft']]
    lg_sens = [sol0['sensitivities']['constants']['C_{lg}_Mission/Aircraft'], sol1['sensitivities']['constants']['C_{lg}_Mission/Aircraft'], \
                 sol2['sensitivities']['constants']['C_{lg}_Mission/Aircraft'], sol3['sensitivities']['constants']['C_{lg}_Mission/Aircraft'], \
                 sol4['sensitivities']['constants']['C_{lg}_Mission/Aircraft']]

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
    plt.savefig('Morphing_Chart_Figs/D8_standard_killer_chart_max_opt_eng.pdf', bbox_inches="tight")
    plt.show(), plt.close()

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
