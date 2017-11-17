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
            'R_{req}': [3000.*units('nmi'),3000.*units('nmi'),3000.*units('nmi')],
            'n_{pass}': [180., 160., 120.],
        })
    else:
        substitutions.update({
           'n_{pass}': 180.,
           'R_{req}': 3000.*units('nmi'),
        })        

    m.substitutions.update(substitutions)

    m = Model(m.cost, BCS(m))
    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

    sol = m_relax.localsolve(verbosity=0, iteration_limit=50, reltol=0.01)
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
            'R_{req}': [3000.*units('nmi'),3000.*units('nmi'),3000.*units('nmi')],
            'n_{pass}': [180., 160., 120.],
        })
    else:
        substitutions.update({
           'n_{pass}': 180.,
           'R_{req}': 3000.*units('nmi'),
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

    #NPASS PLOTS
    #
    # sol = run_sweeps_optimal_737('fuel','n_{pass}',np.linspace(150,210,20))
    #
    # plt.plot(sol('n_{pass}'),sol['sensitivities']['constants']['C_{wing}'],label='Wing',ls='-',color='red',linewidth=1.5)
    # plt.plot(sol('n_{pass}'),sol['sensitivities']['constants']['C_{fuse}'], label = 'Fuselage',ls='--',color='purple',linewidth=1.5)
    # plt.plot(sol('n_{pass}'),sol['sensitivities']['constants']['C_{VT}'], label = 'Vertical Tail',ls=':',color='black',linewidth=1.5)
    # plt.plot(sol('n_{pass}'),sol['sensitivities']['constants']['C_{HT}'], label = 'Horizontal Tail',ls='-.',color='blue',linewidth=1.5)
    # plt.plot(sol('n_{pass}'),sol['sensitivities']['constants']['C_{lg}'], label = 'Landing Gear',ls='--',color='green',linewidth=1.5)
    # plt.xlabel('Number of passengers')
    # plt.ylabel('Sensitivity to component weights')
    # plt.grid()
    # plt.legend(loc = "center right")
    # plt.savefig('figs/sweep-Wsens_vs_npass.png')
    # plt.close()
    #
    # plt.plot(sol('n_{pass}'),sol('W_{f_{total}}'))
    # plt.xlabel('Number of passengers')
    # plt.ylabel('Fuel weight [lbf]')
    # plt.axis([150, 210, 37000, 43000])
    # plt.grid()
    # plt.savefig('figs/sweep-Wftotal_vs_npass.png')
    # plt.close()
    # plt.plot(sol('n_{pass}'),sol('W_{dry}'))
    # plt.xlabel('Number of passengers')
    # plt.ylabel('Dry weight [lbf]')
    # plt.axis([150, 210, 89500, 95500])
    # plt.grid()
    # plt.savefig('figs/sweep-Wdry_vs_npass.png')
    # plt.close()
    #
    # plt.plot(sol('n_{pass}'),sol('W_{f_{total}}')/sol('W_{total}'))
    # plt.xlabel('Number of passengers')
    # plt.ylabel('Fuel weight fraction')
    # plt.ylabel('Fuel weight fraction')
    # plt.grid()
    # plt.savefig('figs/sweep-FuelFraction_vs_npass.png')
    # plt.close()
    #
    # plt.plot(sol('n_{pass}'),sol('W_{f_{total}}')/sol('W_{payload}'))
    # plt.xlabel('Number of passengers')
    # plt.ylabel('Fuel payload weight ratio')
    # plt.grid()
    # plt.savefig('figs/sweep-FPWR_vs_npass.png')
    # plt.close()

    #V1 PLOTS
    #
    # sol = run_sweeps_optimal_737('fuel','V_1',np.linspace(60,80,20))
    #
    # plt.plot(sol('V_1'),sol['sensitivities']['constants']['C_{wing}'],label='Wing')
    # plt.plot(sol('V_1'),sol['sensitivities']['constants']['C_{VT}'], label = 'Vertical Tail')
    # plt.plot(sol('V_1'),sol['sensitivities']['constants']['C_{HT}'], label = 'Horizontal Tail')
    # plt.plot(sol('V_1'),sol['sensitivities']['constants']['C_{fuse}'], label = 'Fuselage')
    # plt.plot(sol('V_1'),sol['sensitivities']['constants']['C_{lg}'], label = 'Landing Gear')
    # plt.xlabel('Minimum takeoff velocity $V_1$')
    # plt.ylabel('Sensitivity to component weights')
    # plt.grid()
    # plt.legend(loc = "upper left")
    # plt.savefig('figs/sweep-Wsens_vs_V1.png')
    # plt.close()
    #
    # plt.plot(sol('V_1'),sol('W_{f_{total}}'))
    # plt.xlabel('Minimum takeoff velocity $V_1$')
    # plt.ylabel('Fuel weight [lbf]')
    # plt.grid()
    # plt.savefig('figs/sweep-Wftotal_vs_V1.png')
    # plt.close()
    #
    # plt.plot(sol('V_1'),sol('W_{dry}'))
    # plt.xlabel('Minimum takeoff velocity $V_1$')
    # plt.ylabel('Dry weight [lbf]')
    # plt.grid()
    # plt.savefig('figs/sweep-Wdry_vs_V1.png')
    # plt.close()
    #
    # plt.plot(sol('V_1'),sol('W_{f_{total}}')/sol('W_{total}'))
    # plt.xlabel('Minimum takeoff velocity $V_1$')
    # plt.ylabel('Fuel weight fraction')
    # plt.ylabel('Fuel weight fraction')
    # plt.grid()
    # plt.savefig('figs/sweep-FuelFraction_vs_V1.png')
    # plt.close()
    #
    # plt.plot(sol('V_1'),sol('W_{f_{total}}')/sol('W_{payload}'))
    # plt.xlabel('Minimum takeoff velocity $V_1$')
    # plt.ylabel('Fuel payload weight ratio')
    # plt.grid()
    # plt.savefig('figs/sweep-FPWR_vs_V1.png')
    # plt.close()
    #
    # plt.plot(sol('V_1'),sol('AR'))
    # plt.xlabel('Minimum takeoff velocity $V_1$')
    # plt.ylabel('Wing aspect ratio')
    # plt.grid()
    # plt.savefig('figs/sweep-AR_vs_V1.png')
    # plt.close()
    #
    # plt.plot(sol('V_1'),sol('S'))
    # plt.xlabel('Minimum takeoff velocity $V_1$')
    # plt.ylabel('Wing area [m$^2$]')
    # plt.grid()
    # plt.savefig('figs/sweep-S_vs_V1.png')
    # plt.close()
    #
    # plt.plot(sol('V_1'),sol('W_{wing}'))
    # plt.xlabel('Minimum takeoff velocity $V_1$')
    # plt.ylabel('Wing weight [lbf]')
    # plt.grid()
    # plt.savefig('figs/sweep-Wwing_vs_V1.png')
    # plt.close()
    #
    # plt.plot(sol('V_1'),sol('W_{total}').transpose()[0]/sol('S').to('in^2'))
    # plt.xlabel('Minimum takeoff velocity $V_1$')
    # plt.ylabel('Wing loading [lbf/m$^2$]')
    # plt.grid()
    # plt.savefig('figs/sweep-WingLoading_vs_V1.png')
    # plt.close()

    #Mmin PLOTS
    #
    # sol = run_sweeps_optimal_737('fuel','M_{min}',np.linspace(0.7,0.83,20))
    #
    # plt.plot(sol('M_{min}'),sol['sensitivities']['constants']['C_{wing}'],label='Wing')
    # plt.plot(sol('M_{min}'),sol['sensitivities']['constants']['C_{VT}'], label = 'Vertical Tail')
    # plt.plot(sol('M_{min}'),sol['sensitivities']['constants']['C_{HT}'], label = 'Horizontal Tail')
    # plt.plot(sol('M_{min}'),sol['sensitivities']['constants']['C_{fuse}'], label = 'Fuselage')
    # plt.plot(sol('M_{min}'),sol['sensitivities']['constants']['C_{lg}'], label = 'Landing Gear')
    # plt.xlabel('Minimum cruise Mach number')
    # plt.ylabel('Sensitivity to component weights')
    # plt.axis([0.7,0.83,0.0,0.3])
    # plt.grid()
    # plt.legend(loc = "center left")
    # plt.savefig('figs/sweep-Wsens_vs_Mmin.png')
    # plt.close()
    #
    # plt.plot(sol('M_{min}'),sol['sensitivities']['constants']['M_{min}'])
    # plt.xlabel('Minimum cruise Mach number')
    # plt.ylabel('Minimum cruise Mach number sensitivity')
    # plt.axis([0.7,0.83,0.0,0.6])
    # plt.grid()
    # plt.savefig('figs/sweep-Mminsens_vs_Mmin.png')
    # plt.close()
    #
    # Marray = np.zeros(20)
    # for i in range(0,20):
    #     Marray[i] = sol('M_Mission/FlightState')[i][0]
    # plt.plot(sol('M_{min}'),Marray)
    # plt.xlabel('Minimum cruise Mach number')
    # plt.ylabel('Beginning of Cruise Mach number')
    # plt.axis([0.7,0.83,0.7,0.83])
    # plt.grid()
    # plt.savefig('figs/sweep-M_vs_Mmin.png')
    # plt.close()
    #
    # plt.plot(sol('M_{min}'),sol('W_{f_{total}}'))
    # plt.xlabel('Minimum cruise Mach number')
    # plt.ylabel('Fuel weight [lbf]')
    # plt.axis([0.7,0.83,27500,30000])
    # plt.grid()
    # plt.savefig('figs/sweep-Wftotal_vs_Mmin.png')
    # plt.close()
    #
    # plt.plot(sol('M_{min}'),sol('W_{dry}'))
    # plt.xlabel('Minimum cruise Mach number')
    # plt.ylabel('Dry weight [lbf]')
    # plt.axis([0.7,0.83,89500,92000])
    # plt.grid()
    # plt.savefig('figs/sweep-Wdry_vs_Mmin.png')
    # plt.close()
    #
    # plt.plot(sol('M_{min}'),sol('W_{f_{total}}')/sol('W_{total}'))
    # plt.xlabel('Minimum cruise Mach number')
    # plt.ylabel('Fuel weight fraction')
    # plt.ylabel('Fuel weight fraction')
    # plt.grid()
    # plt.savefig('figs/sweep-FuelFraction_vs_Mmin.png')
    # plt.close()
    #
    # plt.plot(sol('M_{min}'),sol('W_{f_{total}}')/sol('W_{payload}'))
    # plt.xlabel('Minimum cruise Mach number')
    # plt.ylabel('Fuel payload weight ratio')
    # plt.grid()
    # plt.savefig('figs/sweep-FPWR_vs_Mmin.png')
    # plt.close()
    #
    # plt.plot(sol('M_{min}'),sol('AR'))
    # plt.xlabel('Minimum cruise Mach number')
    # plt.ylabel('Wing aspect ratio')
    # plt.axis([0.7,0.83,12.,13.5])
    # plt.grid()
    # plt.savefig('figs/sweep-AR_vs_Mmin.png')
    # plt.close()
    #
    # plt.plot(sol('M_{min}'),sol('S'))
    # plt.xlabel('Minimum cruise Mach number')
    # plt.ylabel('Wing area [m$^2$]')
    # plt.grid()
    # plt.savefig('figs/sweep-S_vs_Mmin.png')
    # plt.close()
    #
    # plt.plot(sol('M_{min}'),sol('W_{wing}'))
    # plt.xlabel('Minimum cruise Mach number')
    # plt.ylabel('Wing weight [lbf]')
    # plt.axis([0.7,0.83,22000,24500])
    # plt.grid()
    # plt.savefig('figs/sweep-Wwing_vs_Mmin.png')
    # plt.close()
    #
    # plt.plot(sol('M_{min}'),sol('W_{total}').transpose()[0]/sol('S').to('in^2'))
    # plt.xlabel('Minimum cruise Mach number')
    # plt.ylabel('Wing loading [lbf/m$^2$]')
    # plt.grid()
    # plt.savefig('figs/sweep-WingLoading_vs_Mmin.png')
    # plt.close()

