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
            'n_{pax}': [180., 160., 120.],
        })
    else:
        substitutions.update({
           'n_{pax}': 180.,
           'ReqRng': 3000.*units('nmi'),
        })        

    m.substitutions.update(substitutions)

    m = Model(m.cost, BCS(m))
    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

    sol = m_relax.localsolve(verbosity=4, iteration_limit=50, reltol=0.01)
    post_process(sol)

    return sol

def run_sweeps_optimal_737(objective = 'fuel',variable = 'n_{pax}', range = np.linspace(150,210,10)):
    # User definitions
    Ncruise = 4
    Nmission = 1
    aircraft = 'optimal737'

    m  = Mission(Ncruise, objective, aircraft, Nmission)

    substitutions = get737_optimal_subs()

    if Nmission == 3:
        substitutions.update({
            'ReqRng': [3000.*units('nmi'),3000.*units('nmi'),3000.*units('nmi')],
            'n_{pax}': [180., 160., 120.],
        })
    else:
        substitutions.update({
           # 'n_{pax}': 180.,
           'ReqRng': 3000.*units('nmi'),
        })

    m = Model(m.cost, BCS(m))
    m.substitutions.update(substitutions)

    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])
    m_relax.substitutions.update({variable : ('sweep', range)})


    sol = m_relax.localsolve(verbosity=0, iteration_limit=200, reltol=0.01)

    return sol

def test():
    run_optimal_737()

if __name__ == "__main__":
    sol = run_optimal_737()
    percent_diff(sol, 'optimal737')
    print sol.table()

    # sol = run_sweeps_optimal_737(objective = 'fuel',variable = 'n_{pax}', range = np.linspace(150,210,11)):
    # sol = run_sweeps_optimal_737(objective = 'fuel',variable = 'V_1', range = np.linspace(60,80,11)):
    # sol = run_sweeps_optimal_737(objective = 'fuel',variable = 'M_{min}', range = np.linspace(0.7,0.88,10))

