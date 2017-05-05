# only needed for plotting
import matplotlib.pyplot as plt

from SP_aircraft import run_737800
from SP_aircraft import run_D82
from SP_aircraft import run_M072_737
from SP_aircraft import run_optimal_737_fixedBPR
from SP_aircraft import run_D8_eng_wing_opt
from SP_aircraft import run_optimal_737
from SP_aircraft import run_optimal_D8
from SP_aircraft import run_D8_no_BLI

sol0 = run_optimal_737_fixedBPR()
wf0 = sol0('W_{f_{total}}')

sol1 = run_optimal_737()
wf1 = sol1('W_{f_{total}}')

sol2 = run_M072_737()
wf2 = sol2('W_{f_{total}}')

sol3 = run_D8_eng_wing_opt()
wf3 = sol3('W_{f_{total}}')

sol4 = run_D8_no_BLI()
wf4 = sol4('W_{f_{total}}')

sol5 = run_optimal_D8()
wf5 = sol5('W_{f_{total}}')

ytest = [1, wf1/wf0, wf2/wf0, wf3/wf0, wf4/wf0, wf5/wf0]
xtest = [0, 1, 2, 3, 4, 5]
xlabels = ['Optimized 737-800 M = 0.8', 'Optomized BPR', 'Slow to M = 0.72', 'D8 fuselage, Pi tail', 'Rear podded egnines', 'Integrated engines, BLI = D8']

plt.plot(xtest, ytest, "o--")
plt.xticks(xtest, xlabels,  rotation='vertical')
plt.ylim([0,1.1])
plt.xlim([-.5, 5.5])
plt.xlabel('Design Step', fontsize = 20)
plt.ylabel('$W_{f}/W_{f_0}$', fontsize = 20)
plt.savefig('killer_chart_test.pdf', bbox_inches="tight")
plt.show()
