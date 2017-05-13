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
    plt.savefig('D8_M08_morphing_chart_optimal_engines.pdf', bbox_inches="tight")
    plt.show()

def M_08_killer_plot_optimal():
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
    plt.savefig('D8_M08_morphing_chart_optimal_engines.pdf', bbox_inches="tight")
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

    sol4 = run_M08_D8(True)
    wf4 = sol4('W_{f_{total}}')

    sol5 = run_optimal_D8(True)
    wf5 = sol5('W_{f_{total}}')

    ytest = [1, wf1/wf0, wf2/wf0, wf3/wf0, wf4/wf0, wf5/wf0]
    xtest = [0, 1, 2, 3, 4, 5]
    xlabels = ['Optimized 737-800 M = 0.8', 'D8 Fuselage, Pi Tail', 'Rear Podded Engines', 'Integrated Engines, BLI', 'Optimize BPR', 'Slow to M = 0.72']

    plt.plot(xtest, ytest, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,1.1])
    plt.xlim([-.5, 5.5])
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('$W_{f}/W_{f_0}$', fontsize = 20)
    plt.savefig('D8_M08_morphing_chart_fixed_BPR.pdf', bbox_inches="tight")
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

    ytest = [1, wf1/wf0, wf2/wf0, wf3/wf0, wf4/wf0, wf5/wf0]
    xtest = [0, 1, 2, 3, 4, 5]
    xlabels = ['Optimized 737-800 M = 0.8', 'Slow to M = 0.72', 'D8 fuselage, Pi tail', 'Rear podded engines', 'Integrated engines, BLI = D8', 'Optimize engine']

    plt.plot(xtest, ytest, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,1.1])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('$W_{f}/W_{f_0}$', fontsize = 20)
    plt.title('D8 Morphing Chart')
    plt.savefig('D8_standard_killer_chart.pdf', bbox_inches="tight")
    plt.show()

def standard_killer_plot_optimal_engine():
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

    ytest = [1, wf1/wf0, wf2/wf0, wf3/wf0, wf4/wf0]
    xtest = [0, 1, 2, 3, 4]
    xlabels = ['Optimized 737-800 M = 0.8', 'Slow to M = 0.72', 'D8 fuselage, Pi tail', 'Rear podded engines', 'Integrated engines, BLI = D8']

    plt.plot(xtest, ytest, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,1.1])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('$W_{f}/W_{f_0}$', fontsize = 20)
    plt.title('D8 Morphing Chart')
    plt.savefig('D8_standard_killer_chart_optimal_engines.pdf', bbox_inches="tight")
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

    ytest = [1, wf1/wf0, wf2/wf0, wf3/wf0, wf4/wf0]
    xtest = [0, 1, 2, 3, 4]
    xlabels = ['Optimized 737-800 M = 0.8', 'Slow to M = 0.72', 'D8 fuselage, Pi tail', 'Rear podded engines', 'Integrated engines, BLI = D8']

    plt.plot(xtest, ytest, "o--")
    plt.xticks(xtest, xlabels,  rotation='vertical')
    plt.ylim([0,1.1])
    plt.xlim([-.5, 5.5])
    plt.grid()
    plt.xlabel('Design Step', fontsize = 20)
    plt.ylabel('$W_{f}/W_{f_0}$', fontsize = 20)
    plt.title('D8 Morphing Chart')
    plt.savefig('D8_standard_killer_chart_max_optimal_engines.pdf', bbox_inches="tight")
    plt.show()
