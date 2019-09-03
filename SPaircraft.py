"""
Script to run the SP aircraft model
"""
from __future__ import print_function
from __future__ import absolute_import

# GPkit tools
from gpkit import units, Model
from gpkit import Variable, Model, units, SignomialsEnabled, SignomialEquality, Vectorize
from gpkit.constraints.bounded import Bounded

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
#from gen_plots import *

# File for calculating values after solution
from post_compute import post_compute

# Solution check tool relative to TASOPT
from percent_diff import percent_diff

# VSP visualization tools
from saveSol import updateOpenVSP, gendes, gencsm

# Aircraft options:
# currently one of: 'D8_eng_wing', 'optimal737', 'optimal777', 'optimalD8', 'D8_no_BLI', 'M072_737'

def optimize_aircraft(m, substitutions, fixedBPR=False, pRatOpt=True, mutategparg=False, x0 = None):
    """
    Optimizes an aircraft of a given configuration
    :param m: aircraft model with objective and configuration
    :param fixedBPR: boolean specifying whether or not BPR is fixed (depends on config)
    :param pRatOpt: boolean specifying whether or not pressure ratio is optimized (depends on config)
    :param mutategparg: boolean whether to keep each GP solve intact
    :return: solution of aircraft model
    """

    if fixedBPR:
        substitutions.update({
            '\\alpha_{max}': 6.97, #8.62,
        })

    if pRatOpt:
        del substitutions['\pi_{f_D}']
        del substitutions['\pi_{lc_D}']
        del substitutions['\pi_{hc_D}']

    m.substitutions.update(substitutions)
    m_relax = Model(m.cost, Bounded(m), m.substitutions)
    m_relax = relaxed_constants(m_relax)
    sol = m.localsolve(verbosity=2, iteration_limit=200, reltol=0.01, mutategp=mutategparg)#, x0 = x0)
    post_process(sol)
    return sol

def test():
    Nclimb = 3 # number of climb segments
    Ncruise = 2 # number of cruise segments
    Nmission = 1 # number of missions
    config = 'optimalD8' # String describing configuration:
    # currently one of: 'D8_eng_wing', 'optimal737', 'optimal777', 'optimalD8', 'D8_no_BLI', 'M072_737'
    m = Mission(Nclimb, Ncruise, config, Nmission)

    # Objective
    m.cost = m['W_{f_{total}}'].sum()

    # Inputs to the model
    substitutions = get_optimalD8_subs()
    substitutions.update({'R_{req}': 3000.*units('nmi'), #6000*units('nmi'),
                         'n_{pass}': 180.})              #450.,)

    # Additional options
    fixedBPR = False
    pRatOpt = True
    mutategparg = False
    sol = optimize_aircraft(m, substitutions, fixedBPR, pRatOpt, mutategparg)
    Nclimb = m.Nclimb
    percent_diff(sol, config, Nclimb)
    post_compute(sol, Nclimb)
    sol.savetxt()
    return sol

sol = test()
