"""
Script to run the SP aircraft model
"""
# Import constant relaxation tool
from relaxed_constants import relaxed_constants, post_process

# Import tool to check solution relative to TASOPT
from D8_TASOPT_percent_diff import percent_diff

# Import VSP generation tools
from genVSP import updateOpenVSP, genDesFile, genDesFileSweep

#import substitution dict files
from subsD80 import getD80subs
from subsD82 import getD82subs
from subsD82_73eng import getD82_73engsubs
from subs_D8_eng_wing import get_D8_eng_wing_subs
from subsD8big import getD8bigsubs
from subsb737800 import getb737800subs
from subsb777300ER import getb777300ERsubs
from subs_optimal_737 import get737_optimal_subs
from subs_optimal_D8 import get_optimal_D8_subs
from subs_M08_D8 import subs_M08_D8
from subs_M08_d8_eng_wing import getM08_D8_eng_wing_subs
from subs_D8_eng_wing import get_D8_eng_wing_subs
from subsM072737 import get_M072_737_subs
from subs_D8_no_BLI import get_D8_no_BLI_subs
from subs_M08_D8_noBLI import get_subs_M08_D8_noBLI
from subs_optimal_777300ER import get_optimal_777300ER_subs
from subs_M08_D8_big import getD8big_M08_subs
from subs_M072_D8_big import getD8big_M072_subs
from subs_D8big_eng_wing import getD8big_eng_wing_subs
from subs_D8big_noBLI import getD8big_noBLI_subs
from subs_optimal_RJ import get_optimal_RJ_subs
from subs_D8_small import get_small_D8_subs
from subs_D8_small_no_BLI import get_small_D8_no_BLI_subs
from subs_D8_small_eng_wing import get_small_D8_eng_wing_subs

from gpkit import units, Model
from gpkit import Variable, Model, units, SignomialsEnabled, SignomialEquality, Vectorize
from gpkit.constraints.bounded import Bounded as BCS

from D8 import Mission

#import needed for plotting
import matplotlib.pyplot as plt
from gpkit.small_scripts import mag
import numpy as np

def gen_plots(sol):
    """
    function to generate plots of interesting values
    """

    #generate an altitude profile plot
    rng = []
    alt = []
    for i in range(len(sol('RngClimb'))):
           rng.append(mag(sol('RngClimb')[i][0]))
    for i in range(len(sol('Rng'))):
           rng.append(mag(sol('Rng')[i][0]))
    for i in range(len(sol('hft')['hft_Mission, FlightState, Altitude'])):
           alt.append(sol('hft')['hft_Mission, FlightState, Altitude'][i][0])
    rng = np.cumsum(rng)
    plt.plot(rng, alt)
    plt.ylabel('Altitude [feet]', fontsize=18)
    plt.xlabel('Down Range Distnace', fontsize=18)
    plt.title('Aircraft Altitude Profile')
#    plt.savefig('M08_D8_wing_profile_drag.pdf', bbox_inches="tight")
    plt.show()

def run_737800():
    # User definitions
    Nclimb = 3
    Ncruise = 2
    Nmission = 1
    objective = 'fuel'
    aircraft = 'b737800'

    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
    
    substitutions = getb737800subs()

    substitutions.update({
#                'n_{paxx}': 180.,
        'ReqRng': 3000.*units('nmi'),
    })

    m.substitutions.update(substitutions)

    m = Model(m.cost, BCS(m))
    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
    post_process(sol)

    percent_diff(sol, aircraft, Nclimb)

    return sol

def run_D82():
    # User definitions
    Nclimb = 3
    Ncruise = 2
    Nmission = 1
    objective = 'fuel'
    aircraft = 'D82'

    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
    
    substitutions = getD82subs()

    substitutions.update({
#                'n_{paxx}': 180.,
        'ReqRng': 3000.*units('nmi'),
    })

    m.substitutions.update(substitutions)

    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
    post_process(sol)

    percent_diff(sol, aircraft, Nclimb)

    return sol

def run_D8_no_BLI(fixedBPR, pRatOpt = False):
    # User definitions
    Nclimb = 3
    Ncruise = 2
    Nmission = 1
    objective = 'fuel'
    aircraft = 'D8_no_BLI'

    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
    
    substitutions = get_D8_no_BLI_subs()

    substitutions.update({
#                'n_{paxx}': 180.,
        'ReqRng': 3000.*units('nmi'),
    })

    if fixedBPR:
        substitutions.update({
            '\\alpha_{max}': 6.97,
        })

    if pRatOpt:
        del substitutions['\pi_{f_D}']
        del substitutions['\pi_{lc_D}']
        del substitutions['\pi_{hc_D}']

    m.substitutions.update(substitutions)
    m = Model(m.cost, BCS(m))
    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
    post_process(sol)

    percent_diff(sol, 'D82', Nclimb)

    return sol

def run_M072_737(fixedBPR, pRatOpt = False):
    # User definitions
    Nclimb = 3
    Ncruise = 2
    Nmission = 1
    objective = 'fuel'
    aircraft = 'M072_737'

    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
    
    substitutions = get_M072_737_subs()

    substitutions.update({
#                'n_{paxx}': 180.,
        'ReqRng': 3000.*units('nmi'),
    })

    if fixedBPR:
        substitutions.update({
            '\\alpha_{max}': 6.97,
        })

    if pRatOpt:
        del substitutions['\pi_{f_D}']
        del substitutions['\pi_{lc_D}']
        del substitutions['\pi_{hc_D}']

    m.substitutions.update(substitutions)

    m = Model(m.cost, BCS(m))
    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
    post_process(sol)

    percent_diff(sol, 'b737800', Nclimb)

    return sol

def run_D8_eng_wing(fixedBPR, pRatOpt = False):
    # User definitions
    Nclimb = 3
    Ncruise = 2
    Nmission = 1
    objective = 'fuel'
    aircraft = 'D8_eng_wing'

    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
    
    substitutions = get_D8_eng_wing_subs()

    substitutions.update({
#                'n_{paxx}': 180.,
        'ReqRng': 3000.*units('nmi'),
    })

    if fixedBPR:
        substitutions.update({
            '\\alpha_{max}': 6.97,
        })

    if pRatOpt:
        del substitutions['\pi_{f_D}']
        del substitutions['\pi_{lc_D}']
        del substitutions['\pi_{hc_D}']

    m.substitutions.update(substitutions)
    m = Model(m.cost, BCS(m))
    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
    post_process(sol)

    percent_diff(sol, 2, Nclimb)

    return sol

def run_optimal_D8(fixedBPR, pRatOpt = False):
    # User definitions
    Nclimb = 3
    Ncruise = 2
    Nmission = 1
    objective = 'fuel'
    aircraft = 'optimalD8'

    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
    
    substitutions = get_optimal_D8_subs()

    substitutions.update({
#                'n_{paxx}': 180.,
        'ReqRng': 3000.*units('nmi'),
    })

    if fixedBPR:
        substitutions.update({
            '\\alpha_{max}': 6.97,
        })

    if pRatOpt:
        del substitutions['\pi_{f_D}']
        del substitutions['\pi_{lc_D}']
        del substitutions['\pi_{hc_D}']

    m = Model(m.cost, BCS(m))
    m.substitutions.update(substitutions)

    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
    post_process(sol)

    percent_diff(sol, 'D82', Nclimb)

    return sol

def run_optimal_737(fixedBPR, pRatOpt = False):
    # User definitions
    Nclimb = 3
    Ncruise = 2
    Nmission = 1
    objective = 'fuel'
    aircraft = 'optimal737'

    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
    
    substitutions = get737_optimal_subs()

    substitutions.update({
#                'n_{paxx}': 180.,
        'ReqRng': 3000.*units('nmi'),
    })

    if fixedBPR:
        substitutions.update({
            '\\alpha_{max}': 6.97,
        })
        
    if pRatOpt:
        del substitutions['\pi_{f_D}']
        del substitutions['\pi_{lc_D}']
        del substitutions['\pi_{hc_D}']

    m.substitutions.update(substitutions)

    m = Model(m.cost, BCS(m))
    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
    post_process(sol)

    percent_diff(sol, 'b737800', Nclimb)

    return sol

def run_M08_D8_eng_wing(fixedBPR, pRatOpt = False):
    # User definitions
    Nclimb = 3
    Ncruise = 2
    Nmission = 1
    objective = 'fuel'
    aircraft = 'M08_D8_eng_wing'

    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
    
    substitutions = getM08_D8_eng_wing_subs()

    substitutions.update({
#                'n_{paxx}': 180.,
        'ReqRng': 3000.*units('nmi'),
    })

    if fixedBPR:
        substitutions.update({
            '\\alpha_{max}': 6.97,
        })

    if pRatOpt:
        del substitutions['\pi_{f_D}']
        del substitutions['\pi_{lc_D}']
        del substitutions['\pi_{hc_D}']

    m.substitutions.update(substitutions)

    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
    post_process(sol)

    percent_diff(sol, 'D82', Nclimb)

    return sol

def run_M08_D8(fixedBPR, pRatOpt = False):
    # User definitions
    Nclimb = 3
    Ncruise = 2
    Nmission = 1
    objective = 'fuel'
    aircraft = 'M08D8'

    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
    
    substitutions = subs_M08_D8()

    substitutions.update({
#                'n_{paxx}': 180.,
        'ReqRng': 3000.*units('nmi'),
    })

    if fixedBPR:
        substitutions.update({
            '\\alpha_{max}': 6.97,
        })

    if pRatOpt:
        del substitutions['\pi_{f_D}']
        del substitutions['\pi_{lc_D}']
        del substitutions['\pi_{hc_D}']
        
    m.substitutions.update(substitutions)
    m = Model(m.cost, BCS(m))
    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
    post_process(sol)

    percent_diff(sol, 'D82', Nclimb)

    return sol

def run_M08_D8_no_BLI(fixedBPR, pRatOpt = False):
    # User definitions
    Nclimb = 3
    Ncruise = 2
    Nmission = 1
    objective = 'fuel'
    aircraft = 'M08D8_noBLI'

    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
    
    substitutions = get_subs_M08_D8_noBLI()

    substitutions.update({
#                'n_{paxx}': 180.,
        'ReqRng': 3000.*units('nmi'),
    })

    if fixedBPR:
        substitutions.update({
            '\\alpha_{max}': 6.97,
        })

    if pRatOpt:
        del substitutions['\pi_{f_D}']
        del substitutions['\pi_{lc_D}']
        del substitutions['\pi_{hc_D}']
        
    m.substitutions.update(substitutions)
    m = Model(m.cost, BCS(m))
    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
    post_process(sol)

    percent_diff(sol, 'D82', Nclimb)

    return sol

def run_777300ER():
    # User definitions
    Nclimb = 3
    Ncruise = 2
    Nmission = 1
    objective = 'fuel'
    aircraft = 'b777300ER'

    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
    
    substitutions = getb777300ERsubs()

    substitutions.update({
##        'n_{pax}': [450.],
        'ReqRng': [6000.],
    })

    m.substitutions.update(substitutions)
    m = Model(m.cost, BCS(m))
    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
    post_process(sol)

    percent_diff(sol, aircraft, Nclimb)

    return sol

def run_optimal_777(fixedBPR, pRatOpt = False):
    # User definitions
    Nclimb = 3
    Ncruise = 2
    Nmission = 1
    objective = 'fuel'
    aircraft = 'optimal777'

    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
    
    substitutions = get_optimal_777300ER_subs()

    substitutions.update({
##        'n_{pax}': [450.],
        'ReqRng': [6000.],
    })

    if fixedBPR:
        substitutions.update({
            '\\alpha_{max}': 8.62, 
        })

    if pRatOpt:
        del substitutions['\pi_{f_D}']
        del substitutions['\pi_{lc_D}']
        del substitutions['\pi_{hc_D}']

    m.substitutions.update(substitutions)
    m = Model(m.cost, BCS(m))
    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
    post_process(sol)

    percent_diff(sol, 'b777300ER', Nclimb)

    return sol

def run_D8_big(fixedBPR, pRatOpt = False):
    # User definitions
    Nclimb = 3
    Ncruise = 2
    Nmission = 1
    objective = 'fuel'
    aircraft = 'D8big'

    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
    
    substitutions = getD8bigsubs()

    substitutions.update({
#                'n_{paxx}': 180.,
        'ReqRng': 6000.*units('nmi'),
    })

    if fixedBPR:
        substitutions.update({
            '\\alpha_{max}': 8.62, 
        })

    if pRatOpt:
        del substitutions['\pi_{f_D}']
        del substitutions['\pi_{lc_D}']
        del substitutions['\pi_{hc_D}']
        
    m.substitutions.update(substitutions)
    m = Model(m.cost, BCS(m))
    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
    post_process(sol)

    percent_diff(sol, 'b777300ER', Nclimb)

    return sol

def run_D8_big_no_BLI(fixedBPR, pRatOpt = False):
    # User definitions
    Nclimb = 3
    Ncruise = 2
    Nmission = 1
    objective = 'fuel'
    aircraft = 'D8big_no_BLI'

    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
    
    substitutions = getD8big_noBLI_subs()

    substitutions.update({
#                'n_{paxx}': 180.,
        'ReqRng': 6000.*units('nmi'),
    })

    if fixedBPR:
        substitutions.update({
            '\\alpha_{max}': 8.62, 
        })

    if pRatOpt:
        del substitutions['\pi_{f_D}']
        del substitutions['\pi_{lc_D}']
        del substitutions['\pi_{hc_D}']
        
    m.substitutions.update(substitutions)
    m = Model(m.cost, BCS(m))
    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
    post_process(sol)

    percent_diff(sol, 'b777300ER', Nclimb)

    return sol

def run_D8_big_eng_wing(fixedBPR, pRatOpt = False):
    # User definitions
    Nclimb = 3
    Ncruise = 2
    Nmission = 1
    objective = 'fuel'
    aircraft = 'D8big_eng_wing'

    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
    
    substitutions = getD8big_eng_wing_subs()

    substitutions.update({
#                'n_{paxx}': 180.,
        'ReqRng': 6000.*units('nmi'),
    })

    if fixedBPR:
        substitutions.update({
            '\\alpha_{max}': 8.62, 
        })

    if pRatOpt:
        del substitutions['\pi_{f_D}']
        del substitutions['\pi_{lc_D}']
        del substitutions['\pi_{hc_D}']
        
    m.substitutions.update(substitutions)
    m = Model(m.cost, BCS(m))
    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
    post_process(sol)

    percent_diff(sol, 'b777300ER', Nclimb)

    return sol


def run_D8_big_M072(fixedBPR, pRatOpt = False):
    # User definitions
    Nclimb = 3
    Ncruise = 2
    Nmission = 1
    objective = 'fuel'
    aircraft = 'D8big'

    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
    
    substitutions = getD8big_M072_subs()

    substitutions.update({
#                'n_{paxx}': 180.,
        'ReqRng': 6000.*units('nmi'),
    })

    if fixedBPR:
        substitutions.update({
            '\\alpha_{max}': 8.62, 
        })

    if pRatOpt:
        del substitutions['\pi_{f_D}']
        del substitutions['\pi_{lc_D}']
        del substitutions['\pi_{hc_D}']
        
    m.substitutions.update(substitutions)
    m = Model(m.cost, BCS(m))
    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
    post_process(sol)

    percent_diff(sol, 'b777300ER', Nclimb)

    return sol

def run_D8_big_M08(fixedBPR, pRatOpt = False):
    # User definitions
    Nclimb = 3
    Ncruise = 2
    Nmission = 1
    objective = 'fuel'
    aircraft = 'D8big'

    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
    
    substitutions = getD8big_M08_subs()

    substitutions.update({
#                'n_{paxx}': 180.,
        'ReqRng': 6000.*units('nmi'),
    })

    if fixedBPR:
        substitutions.update({
            '\\alpha_{max}': 8.62, 
        })

    if pRatOpt:
        del substitutions['\pi_{f_D}']
        del substitutions['\pi_{lc_D}']
        del substitutions['\pi_{hc_D}']
        
    m.substitutions.update(substitutions)
    m = Model(m.cost, BCS(m))
    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
    post_process(sol)

    percent_diff(sol, 'b777300ER', Nclimb)

    return sol

def run_optimal_RJ(fixedBPR, pRatOpt = False):
    # User definitions
    Nclimb = 3
    Ncruise = 2
    Nmission = 1
    objective = 'fuel'
    aircraft = 'optimalRJ'

    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
    
    substitutions = get_optimal_RJ_subs()

    substitutions.update({
##        'n_{pax}': [90.],
        'ReqRng': [1500.],
    })

    if fixedBPR:
        substitutions.update({
            '\\alpha_{max}': 5.6958, 
        })

    if pRatOpt:
        del substitutions['\pi_{f_D}']
        del substitutions['\pi_{lc_D}']
        del substitutions['\pi_{hc_D}']

    m.substitutions.update(substitutions)
    m = Model(m.cost, BCS(m))
    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
    post_process(sol)

    percent_diff(sol, 'b737800', Nclimb)

    return sol

def run_small_D8(fixedBPR, pRatOpt = False):
    # User definitions
    Nclimb = 3
    Ncruise = 2
    Nmission = 1
    objective = 'fuel'
    aircraft = 'smallD8'

    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
    
    substitutions = get_small_D8_subs()

    substitutions.update({
##        'n_{pax}': [90.],
        'ReqRng': [1500.],
    })

    if fixedBPR:
        substitutions.update({
            '\\alpha_{max}': 5.6958, 
        })

    if pRatOpt:
        del substitutions['\pi_{f_D}']
        del substitutions['\pi_{lc_D}']
        del substitutions['\pi_{hc_D}']

    m.substitutions.update(substitutions)
    m = Model(m.cost, BCS(m))
    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
    post_process(sol)

    percent_diff(sol, 'b737800', Nclimb)

    return sol

def run_small_D8_eng_wing(fixedBPR, pRatOpt = False):
    # User definitions
    Nclimb = 3
    Ncruise = 2
    Nmission = 1
    objective = 'fuel'
    aircraft = 'smallD8_eng_wing'

    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
    
    substitutions = get_small_D8_eng_wing_subs()

    substitutions.update({
##        'n_{pax}': [90.],
        'ReqRng': [1500.],
    })

    if fixedBPR:
        substitutions.update({
            '\\alpha_{max}': 5.6958, 
        })

    if pRatOpt:
        del substitutions['\pi_{f_D}']
        del substitutions['\pi_{lc_D}']
        del substitutions['\pi_{hc_D}']

    m.substitutions.update(substitutions)
    m = Model(m.cost, BCS(m))
    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
    post_process(sol)

    percent_diff(sol, 'b737800', Nclimb)

    return sol

def run_small_D8_no_BLI(fixedBPR, pRatOpt = False):
    # User definitions
    Nclimb = 3
    Ncruise = 2
    Nmission = 1
    objective = 'fuel'
    aircraft = 'smallD8_no_BLI'

    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
    
    substitutions = get_small_D8_no_BLI_subs()

    substitutions.update({
##        'n_{pax}': [90.],
        'ReqRng': [1500.],
    })

    if fixedBPR:
        substitutions.update({
            '\\alpha_{max}': 5.6958, 
        })

    if pRatOpt:
        del substitutions['\pi_{f_D}']
        del substitutions['\pi_{lc_D}']
        del substitutions['\pi_{hc_D}']

    m.substitutions.update(substitutions)
    m = Model(m.cost, BCS(m))
    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
    post_process(sol)

    percent_diff(sol, 'b737800', Nclimb)

    return sol

def test():
    run_737800()

if __name__ == '__main__':
    # User definitions
    Nclimb = 3
    Ncruise = 2
    Nmission = 1
    objective = 'fuel'
    aircraft = 'D8_eng_wing'

    genVSP = True
    sweeps = False
    nsweep = 5

    if Nmission == 1:
        multimission = True
    else:
        multimission = False

    # Mission definition
    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)

    # Aircraft types
    if aircraft == 'D80':
        print('D80 executing...')
        substitutions = getD80subs()
        if Nmission == 1:
                substitutions.update({
##                 'n_{pax}': 180.,
                'ReqRng': 3000.*units('nmi'),
                })

        if Nmission != 1:
                substitutions.update({
                'n_{pax}': [180.],
                'ReqRng': [3000.],
                })

    if aircraft == 'D8_eng_wing':
        print('D8_eng_wing executing...')
        substitutions = get_D8_eng_wing_subs()
        if Nmission == 1:
                substitutions.update({
##                 'n_{pax}': 180.,
                'ReqRng': 3000.*units('nmi'),
                })

        if Nmission != 1:
                substitutions.update({
                'n_{pax}': [180.],
                'ReqRng': [3000.],
                })

    if aircraft == 'D8_no_BLI':
        print('D8_no_BLI executing...')
        substitutions = get_D8_no_BLI_subs()
        if Nmission == 1:
                substitutions.update({
##                 'n_{pax}': 180.,
                'ReqRng': 3000.*units('nmi'),
                })

        if Nmission != 1:
                substitutions.update({
                'n_{pax}': [180.],
                'ReqRng': [3000.],
                })


    if aircraft == 'D82_73eng':
        print('D82_73eng executing...')
        substitutions = getD82_73engsubs()
        if Nmission == 1:
                substitutions.update({
##                'n_{pax}': 180.,
                'ReqRng': 3000.*units('nmi'),
                })
        if Nmission != 1:
                substitutions.update({
                'n_{pax}': [180.],
                'ReqRng': [3000.],
                })


    m.substitutions.update(substitutions)

    if objective == 'PRFC':
        # PRFC optimization chooses optimal mission for a given configuration
        # TO RUN: MUST REMOVE EQUATIONS SETTING N_{PAX} IN D8.
        # m.substitutions.__delitem__('n_{pax}')
        m.substitutions.__delitem__('ReqRng')

    if aircraft in ['D80','D82','D8_no_BLI']:
        # m = Model(m.cost,BCS(m))
        m_relax = relaxed_constants(m, None, ['ReqRng'])
    if aircraft in ['D8big', 'D82_73eng', 'D8_eng_wing', 'optimalD8', 'M08D8', 'M08_D8_eng_wing']:
        m = Model(m.cost,BCS(m))
        m_relax = relaxed_constants(m, None, ['ReqRng'])
    if aircraft in ['b737800', 'optimal737', 'M072_737']:
        m = Model(m.cost, BCS(m))
        m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])
    if aircraft in ['b777300ER']:
        m = Model(m.cost, BCS(m))
        m_relax = relaxed_constants(m, None)

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
    post_process(sol)
    
    if Nmission == 1:
         if aircraft in ['D82', 'D8_eng_wing', 'optimalD8']:
              percent_diff(sol, 2, Nclimb)

         if aircraft in ['b737800','optimal737']:
              percent_diff(sol, 801, Nclimb)

         if aircraft in ['b777300ER']:
              percent_diff(sol, 777, Nclimb)
    if genVSP:
        if sweeps:
            genDesFileSweep(sol,aircraft,nsweep)
        else:
            genDesFile(sol,aircraft)
