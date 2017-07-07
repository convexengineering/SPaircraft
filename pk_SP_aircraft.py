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
from D8 import Mission

# Substitution dictionaries for different aircraft
from subs_optimal_737 import get737_optimal_subs
# Plotting tools
import matplotlib.pyplot as plt
from gpkit.small_scripts import mag
import numpy as np

# File for calculating values after solution
from post_compute import post_compute

# Solution check tool relative to TASOPT
from D8_TASOPT_percent_diff import percent_diff

def run_optimal_737(fixedBPR, pRatOpt = False, objective = 'fuel'):
    # User definitions
    Ncruise = 4
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
#            'n_{paxx}': 180.,
            'ReqRng': 3000.*units('nmi'),
        })        

    if fixedBPR:
        substitutions.update({
            '\\alpha_{max}': 6.97,
        })
        
    if pRatOpt:
        del substitutions['\pi_{f_D}']
##        del substitutions['\pi_{lc_D}']
##        del substitutions['\pi_{hc_D}']

    m.substitutions.update(substitutions)

    m = Model(m.cost, BCS(m))
    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
    post_process(sol)

    percent_diff(sol, 'optimal737')

    return sol

