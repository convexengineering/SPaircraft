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
from subs.optimalD8 import get_optimalD8_subs
from subs.optimal777 import get_optimal777_subs
from subs.optimal737 import get_optimal737_subs
from subs.M072_737 import get_M072_737_subs
from subs.D8_no_BLI import get_D8_no_BLI_subs
from subs.D8_eng_wing import get_D8_eng_wing_subs

# Plotting tools
from gen_plots import *

# File for calculating values after solution
from post_compute import post_compute

# Solution check tool relative to TASOPT
from percent_diff import percent_diff

# VSP visualization tools
from saveSol import updateOpenVSP, genDesFile, genDesFileSweep

# Aircraft options:
# 'D8_eng_wing', 'optimal737', 'optimal777', 'optimalD8', 'D8_no_BLI', 'M072_737'

def optimize_aircraft(objective, aircraft, substitutions, fixedBPR, pRatOpt, mutategparg):
    # User definitions
    Nclimb = 3
    Ncruise = 2
    Nmission = 1

    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)

    if Nmission == 4:
        substitutions.update({
            'R_{req}': [3000.*units('nmi'),2500.*units('nmi'),2000.*units('nmi'),1000.*units('nmi')], #,2500.*units('nmi')
            'n_{pass}': [180., 180., 180., 180.], #, 140.
        })
    elif Nmission == 2:
        substitutions.update({
            'R_{req}': [3000.*units('nmi'),2000.*units('nmi'),], #,2500.*units('nmi')
            'n_{pass}': [180., 180.], #, 140.
        })
    else:
        substitutions.update({
                # Commented substitutions are for 777 class aircraft
           'R_{req}': 3000.*units('nmi'), #6000*units('nmi'),
           'n_{pass}': 180.,              #450.,
        })

    if fixedBPR:
        substitutions.update({
            '\\alpha_{max}': 6.97, #8.62,
        })

    if pRatOpt:
        del substitutions['\pi_{f_D}']
        del substitutions['\pi_{lc_D}']
        del substitutions['\pi_{hc_D}']

    m.substitutions.update(substitutions)
    m = Model(m.cost, BCS(m))
    m_relax = relaxed_constants(m)

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01, mutategp=mutategparg)
    post_process(sol)

    percent_diff(sol, aircraft, Nclimb)

    post_compute(sol, Nclimb)

    return sol

def test():
    objective = 'W_{f_{total}}'
    aircraft = 'optimal737'
    substitutions = get_optimal737_subs()
    fixedBPR = True
    pRatOpt = False
    mutategparg = True
    optimize_aircraft(objective, aircraft, substitutions, fixedBPR, pRatOpt, mutategparg)
