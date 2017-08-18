"""
Script to run the SP aircraft model
"""

# GPkit tools
from gpkit import units, Model
from gpkit import Variable, Model, units, SignomialsEnabled, SignomialEquality, Vectorize
from gpkit.constraints.bounded import Bounded as BCS

# Constant relaxation heuristic for SP solve
from relaxed_constants import relaxed_constants, post_process

# Mission model
from aircraft import Mission

# Substitution dictionaries for different aircraft
from subs.optimal_737 import get737_optimal_subs
# Plotting tools
import matplotlib.pyplot as plt
from gpkit.small_scripts import mag
import numpy as np

# File for calculating values after solution
from post_compute import post_compute

# Solution check tool relative to TASOPT
from percent_diff_TASOPT import percent_diff

def run_optimal_737(objective = 'fuel'):
    # User definitions
    Ncruise = 2
    Nmission = 1
    aircraft = 'optimal737'

    m = Mission(Ncruise, objective, aircraft, Nmission)
    
    substitutions = get737_optimal_subs()

    if Nmission == 3:
        substitutions.update({
            'ReqRng': [3000.*units('nmi'),3000.*units('nmi'),3000.*units('nmi')],
            'n_{pass}': [180., 160., 120.],
        })
    else:
        substitutions.update({
           'n_{pass}': 180.,
           'ReqRng': 3000.*units('nmi'),
        })        

    m.substitutions.update(substitutions)

    m = Model(m.cost, BCS(m))
    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

    sol = m_relax.localsolve(verbosity=0, iteration_limit=50, reltol=0.01,
                             modifylastgp=False)
    post_process(sol)

    return sol

def run_sweeps_optimal_737(objective = 'fuel',variable = 'n_{pass}', xrange = np.linspace(150,210,20)):
    # User definitions
    Ncruise = 4
    Nmission = 1
    aircraft = 'optimal737'

    m  = Mission(Ncruise, objective, aircraft, Nmission)

    substitutions = get737_optimal_subs()

    if Nmission == 3:
        substitutions.update({
            'ReqRng': [3000.*units('nmi'),3000.*units('nmi'),3000.*units('nmi')],
            'n_{pass}': [180., 160., 120.],
        })
    else:
        substitutions.update({
           'n_{pass}': 180.,
           'ReqRng': 3000.*units('nmi'),
        })

    m = Model(m.cost, BCS(m))
    m.substitutions.update(substitutions)
    m.substitutions.update({variable : ('sweep', xrange)})
    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])
    sol = m_relax.localsolve(verbosity=4, iteration_limit=50, reltol=0.01, skipsweepfailures = True)
    return sol

def test():
    run_optimal_737()

if __name__ == "__main__":
    sol = run_optimal_737()
    percent_diff(sol, 'optimal737')
    print sol.table()

    # sol = run_sweeps_optimal_737('fuel','V_1',np.linspace(60,80,20))
    # sol = run_sweeps_optimal_737('fuel','M_{min}',np.linspace(0.7,0.83,20))


    # NPASS PLOTS
    #
    # sol = run_sweeps_optimal_737('fuel','n_{pass}',np.linspace(150,210,20))
    #
    # plt.plot(sol('n_{pass}'),sol['sensitivities']['constants']['C_{wing}'],label='Wing')
    # plt.plot(sol('n_{pass}'),sol['sensitivities']['constants']['C_{VT}'], label = 'Vertical Tail')
    # plt.plot(sol('n_{pass}'),sol['sensitivities']['constants']['C_{HT}'], label = 'Horizontal Tail')
    # plt.plot(sol('n_{pass}'),sol['sensitivities']['constants']['C_{fuse}'], label = 'Fuselage')
    # plt.plot(sol('n_{pass}'),sol['sensitivities']['constants']['C_{lg}'], label = 'Landing Gear')
    # plt.xlabel('Number of passengers')
    # plt.ylabel('Sensitivity to component weights')
    # plt.grid()
    # plt.legend(loc = "upper left")
    # plt.savefig('figs/sweep-Wsens_vs_npass.png')
    # plt.close()
    #
    # plt.plot(sol('n_{pass}'),sol('W_{f_{total}}'))
    # plt.xlabel('Number of passengers')
    # plt.ylabel('Fuel weight [lbf]')
    # plt.grid()
    # plt.savefig('figs/sweep-Wftotal_vs_npass.png')
    # plt.close()
    #
    # plt.plot(sol('n_{pass}'),sol('W_{dry}'))
    # plt.xlabel('Number of passengers')
    # plt.ylabel('Dry weight [lbf]')
    # plt.grid()
    # plt.savefig('figs/sweep-Wdry_vs_npass.png')
    # plt.close()
    #
    # plt.plot(sol('n_{pass}'),sol('W_{f_{total}}')/sol('W_{total}'))
    # plt.xlabel('Number of passengers')
    # plt.ylabel('Fuel weight fraction')
    # plt.ylabel('Fuel weight fraction')
    # plt.savefig('figs/sweep-FuelFraction_vs_npass.png')
    # plt.grid()
    # plt.close()
    #
    # plt.plot(sol('n_{pass}'),sol('W_{f_{total}}')/sol('W_{payload}'))
    # plt.xlabel('Number of passengers')
    # plt.ylabel('Fuel payload weight ratio')
    # plt.savefig('figs/sweep-FPWR_vs_npass.png')
    # plt.grid()
    # plt.close()

    #V1 PLOTS

    # plt.plot(sol('V_1'),sol['sensitivities']['constants']['C_{wing}'],label='Wing')
    # plt.plot(sol('V_1'),sol['sensitivities']['constants']['C_{VT}'], label = 'Vertical Tail')
    # plt.plot(sol('V_1'),sol['sensitivities']['constants']['C_{HT}'], label = 'Horizontal Tail')
    # plt.plot(sol('V_1'),sol['sensitivities']['constants']['C_{fuse}'], label = 'Fuselage')
    # plt.plot(sol('V_1'),sol['sensitivities']['constants']['C_{lg}'], label = 'Landing Gear')
    # plt.xlabel('Minimum takeoff velocity $V_1$')
    # plt.ylabel('Sensitivity to component weights')
    # plt.grid()
    # plt.legend(loc = "upper left")

    #Mmin PLOTS

    # plt.plot(sol('M_{min}'),sol['sensitivities']['constants']['C_{wing}'],label='Wing')
    # plt.plot(sol('M_{min}'),sol['sensitivities']['constants']['C_{VT}'], label = 'Vertical Tail')
    # plt.plot(sol('M_{min}'),sol['sensitivities']['constants']['C_{HT}'], label = 'Horizontal Tail')
    # plt.plot(sol('M_{min}'),sol['sensitivities']['constants']['C_{fuse}'], label = 'Fuselage')
    # plt.plot(sol('M_{min}'),sol['sensitivities']['constants']['C_{lg}'], label = 'Landing Gear')
    # plt.xlabel('Minimum Mach number')
    # plt.ylabel('Sensitivity to component weights')
    # plt.grid()
    # plt.legend(loc = "upper left")

