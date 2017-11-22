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
from subs.optimal_D8 import get_optimal_D8_subs
from subs.optimal_777300ER import get_optimal_777300ER_subs
from subs.optimal_737 import get737_optimal_subs
from subs.b737_M072 import get_M072_737_subs
from subs.D8_no_BLI import get_D8_no_BLI_subs
from subs.D8_eng_wing import get_D8_eng_wing_subs

#out of date sub files that don't currenlty work
##from subs.D80 import getD80subs
##from subs.D82 import getD82subs
##from subs.D82_737_engine import getD82_73engsubs
##from subs.D8_big import getD8bigsubs
##from subs.b737800 import getb737800subs
##from subs.b777300ER import getb777300ERsubs
##from subs.D8_M08 import subs_M08_D8
##from subs.D8_eng_wing_M08 import getM08_D8_eng_wing_subs
##from subs.D8_eng_wing import get_D8_eng_wing_subs
##from subs.D8_no_BLI_M08 import get_subs_M08_D8_noBLI
##from subs.D8_big_M08 import getD8big_M08_subs
##from subs.D8_big_M072 import getD8big_M072_subs
##from subs.D8_big_eng_wing import getD8big_eng_wing_subs
##from subs.D8_big_no_BLI import getD8big_noBLI_subs
##from subs.optimal_RJ import get_optimal_RJ_subs
##from subs.optimal_RJ_M072 import get_M072_optimal_RJ_subs
##from subs.D8_small import get_small_D8_subs
##from subs.D8_small_no_BLI import get_small_D8_no_BLI_subs
##from subs.D8_small_eng_wing import get_small_D8_eng_wing_subs
##from subs.D8_small_M08 import get_small_M08_D8_subs
##from subs.D8_small_no_BLI_M08 import get_small_M08_D8_no_BLI_subs
##from subs.D8_small_eng_wing_M08  import get_small_M08_D8_eng_wing_subs
##from subs.D8_big_eng_wing_M072 import getD8big_M072_eng_wing_subs
##from subs.D8_big_no_BLI_M072 import getD8big_M072_noBLI_subs
##from subs.optimal_777_M072 import get_optimal_777300ER_M072_subs
##from subs.optimal_777_M08 import get_optimal_777300ER_M08_subs
##from subs.D12  import get_D12_subs

# Plotting tools
import matplotlib.pyplot as plt
from gpkit.small_scripts import mag
import numpy as np

# File for calculating values after solution
from post_compute import post_compute

# Solution check tool relative to TASOPT
from percent_diff import percent_diff

# VSP visualization tools
from saveSol import updateOpenVSP, genDesFile, genDesFileSweep

def gen_plots(sol):
    """
    function to generate plots of interesting values
    """

    #generate an altitude profile plot
    rng = []
    alt = []
    for i in range(len(sol('R_{climb}'))):
           rng.append(mag(sol('R_{climb}')[i][0]))
    for i in range(len(sol('R_{cruise}'))):
           rng.append(mag(sol('R_{cruise}')[i][0]))
    for i in range(len(sol('hft')['hft_Mission/FlightState/Altitude'])):
           alt.append(mag(sol('hft')['hft_Mission/FlightState/Altitude'][i][0]))
    rng = np.cumsum(rng)
    plt.plot(rng, alt)
    plt.ylabel('Altitude [ft]', fontsize=18)
    plt.xlabel('Down Range Distance', fontsize=18)
    plt.title('Aircraft Altitude Profile')
#    plt.savefig('M08_D8_wing_profile_drag.pdf', bbox_inches="tight")
    plt.show()

def gen_D82_plots(sol):
    """
    function to generate plots of interesting values
    """

    #generate an altitude profile plot
    rng = [0]
    alt = [0]
    for i in range(len(sol('R_{climb}'))):
           rng.append(mag(sol('R_{climb}')[i][0]))
    for i in range(len(sol('R_{cruise}'))):
           rng.append(mag(sol('R_{cruise}')[i][0]))
    for i in range(len(sol('hft')['hft_Mission/FlightState/Altitude'])):
           alt.append(mag(sol('hft')['hft_Mission/FlightState/Altitude'][i][0]))
    rng = np.cumsum(rng)

    tasrng = [0, 11.51, 27.36, 52.64, 103.28, 103.28, 2825.41, 2825.41, 2869.08, 2912.76, 2956.43, 3000]
    tasalt = [0, 9619.5, 19239.0, 28858.5, 38478.0, 38478.0, 41681.3, 41681.3, 32129.3, 21998.5, 11288.7, 0]
    
    plt.plot(rng, alt)
    plt.plot(tasrng, tasalt)
    plt.legend(['SP Model', 'TASOPT'], loc=4)
    plt.ylabel('Altitude [ft]', fontsize=18)
    plt.xlabel('Down Range Distance [nm]', fontsize=18)
    plt.title('D8.2 Altitude Profile')
    plt.savefig('D8_altitude_profile.eps', bbox_inches="tight")
    plt.show()

def gen_D8_737_plots(solD8, sol737):
    """
    function to generate plots of interesting values
    """

    #generate an altitude profile plot
    rngD8 = []
    altD8 = []
    for i in range(len(solD8('R_{climb}'))):
           rngD8.append(mag(solD8('R_{climb}')[i][0]))
    for i in range(len(solD8('Rng'))):
           rngD8.append(mag(solD8('Rng')[i][0]))
    for i in range(len(solD8('hft')['hft_Mission/FlightState/Altitude'])):
           altD8.append(mag(solD8('hft')['hft_Mission/FlightState/Altitude'][i][0]))
    rngD8 = np.cumsum(rngD8)

    rng73 = []
    alt73 = []
    for i in range(len(sol737('R_{climb}'))):
           rng73.append(mag(sol737('R_{climb}')[i][0]))
    for i in range(len(sol737('Rng'))):
           rng73.append(mag(sol737('Rng')[i][0]))
    for i in range(len(sol737('hft')['hft_Mission/FlightState/Altitude'])):
           alt73.append(mag(sol737('hft')['hft_Mission/FlightState/Altitude'][i][0]))
    rng73 = np.cumsum(rng73)

    tasrngD8 = [0, 11.51, 27.36, 52.64, 103.28, 103.28, 2825.41, 2825.41, 2869.08, 2912.76, 2956.43, 3000]
    tasaltD8 = [0, 9619.5, 19239.0, 28858.5, 38478.0, 38478.0, 41681.3, 41681.3, 32129.3, 21998.5, 11288.7, 0]

    tasrng73 = [0, 13.68, 31.34, 59.96, 115.05, 115.05, 2875.38, 2875.38, 2906.56, 2937.74, 2968.92, 3000]
    tasalt73 = [0, 8750, 17500, 26250, 35000, 35000, 39677.3, 39677.3, 29758., 19838.6, 9919.3, 0]

    plt.plot(rngD8, altD8)
    plt.plot(tasrngD8, tasaltD8)
    plt.plot(rng73, alt73)
    plt.plot(tasrng73, tasalt73)
    plt.legend(['D8 SP Model', 'D8 TASOPT', '737 SP Model', '737 TASOPT'], loc=4)
    plt.ylabel('Altitude [ft]', fontsize=18)
    plt.xlabel('Down Range Distance [nm]', fontsize=18)
    plt.title('737 and D8 Altitude Profile', fontsize=18)
    plt.savefig('737_D8_altitude_profile.eps', bbox_inches="tight")
    plt.show()

def gen_D8_D8_no_BLI_plots(solD8, solno_BLI):
    """
    function to generate plots of interesting values
    """

    #generate an altitude profile plot
    rngD8 = []
    altD8 = []
    for i in range(len(solD8('R_{climb}'))):
           rngD8.append(mag(solD8('R_{climb}')[i][0]))
    for i in range(len(solD8('Rng'))):
           rngD8.append(mag(solD8('Rng')[i][0]))
    for i in range(len(solD8('hft')['hft_Mission/FlightState/Altitude'])):
           altD8.append(mag(solD8('hft')['hft_Mission/FlightState/Altitude'][i][0]))
    rngD8 = np.cumsum(rngD8)

    rngno_BLI = []
    altno_BLI = []
    for i in range(len(solno_BLI('R_{climb}'))):
           rngno_BLI.append(mag(solno_BLI('R_{climb}')[i][0]))
    for i in range(len(solno_BLI('Rng'))):
           rngno_BLI.append(mag(solno_BLI('Rng')[i][0]))
    for i in range(len(solno_BLI('hft')['hft_Mission/FlightState/Altitude'])):
           altno_BLI.append(mag(solno_BLI('hft')['hft_Mission/FlightState/Altitude'][i][0]))
    rngno_BLI = np.cumsum(rngno_BLI)

    plt.plot(rngD8, altD8)
    plt.plot(rngno_BLI, altno_BLI)
    plt.legend(['D8', 'D8 w/out BLI (rear podded engines)'], loc=4)
    plt.ylabel('Altitude [ft]', fontsize=18)
    plt.xlabel('Down Range Distance [nm]', fontsize=18)
    plt.title('D8 Altitude Profile with and without BLI')
    plt.savefig('D8_D8_no_BLI_altitude_profile.pdf', bbox_inches="tight")
    plt.show()

def gen_737_plots(sol):
    """
    function to generate plots of interesting values
    """

    #generate an altitude profile plot
    tasrng = [0, 13.68, 31.34, 59.96, 115.05, 115.05, 2875.38, 2875.38, 2906.56, 2937.74, 2968.92, 3000]
    tasalt = [0, 8750, 17500, 26250, 35000, 35000, 39677.3, 39677.3, 29758., 19838.6, 9919.3, 0]

    rng = [0]
    alt = [0]

    for i in range(len(sol('R_{climb}'))):
           rng.append(mag(sol('R_{climb}')[i][0]))
    for i in range(len(sol('R_{cruise}'))):
           rng.append(mag(sol('R_{cruise}')[i][0]))
    for i in range(len(sol('hft')['hft_Mission/FlightState/Altitude'])):
           alt.append(mag(sol('hft')['hft_Mission/FlightState/Altitude'][i][0]))
    rng = np.cumsum(rng)
    plt.plot(rng, alt)
    plt.plot(tasrng, tasalt)
    plt.legend(['SP Model', 'TASOPT'], loc=4)
    plt.ylabel('Altitude [ft]', fontsize=18)
    plt.xlabel('Down Range Distance [nm]', fontsize=18)
    plt.title('737 Altitude Profile', fontsize=18)
    plt.savefig('737_altitude_profile.eps', bbox_inches="tight")
    plt.show()
    
def gen_777_plots(sol):
    """
    function to generate plots of interesting values
    """

    #generate an altitude profile plot
    tasrng = [0, 15.6, 33.24, 60.40, 107.98, 107.98, 5850.37, 5850.37, 5887.82, 5925.28, 5962.74, 6000]
    tasalt = [0, 7994.2, 15988.5, 23982.8, 31977.0, 31977.0, 39723.4, 39723.4, 31282.2, 21847.9, 11420.5, 0]

    rng = [0]
    alt = [0]
    
    for i in range(len(sol('R_{climb}'))):
           rng.append(mag(sol('R_{climb}')[i][0]))
    for i in range(len(sol('R_{cruise}'))):
           rng.append(mag(sol('R_{cruise}')[i][0]))
    for i in range(len(sol('hft')['hft_Mission/FlightState/Altitude'])):
           alt.append(mag(sol('hft')['hft_Mission/FlightState/Altitude'][i][0]))
    
    rng = np.cumsum(rng)
    plt.plot(rng, alt)
    plt.plot(tasrng, tasalt)
    plt.legend(['SP Model', 'TASOPT'], loc=4)
    plt.ylabel('Altitude [ft]', fontsize=18)
    plt.xlabel('Down Range Distance [nm]', fontsize=18)
    plt.title('777 Altitude Profile')
    plt.savefig('777_altitude_profile.eps', bbox_inches="tight")
    plt.show()

def run_optimal_777(objective, fixedBPR, pRatOpt, mutategparg):
    # User definitions
    Nclimb = 3
    Ncruise = 2
    Nmission = 1
    aircraft = 'optimal777'

    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
    
    substitutions = get_optimal_777300ER_subs()

    substitutions.update({
        'R_{req}': 6000.,
        'n_{pass}': 450.,
    })

    if fixedBPR:
        substitutions.update({
            '\\alpha_{max}': 8.62, 
        })

    if pRatOpt:
        del substitutions['\pi_{f_D}']
##        del substitutions['\pi_{lc_D}']
##        del substitutions['\pi_{hc_D}']

    m.substitutions.update(substitutions)
    m = Model(m.cost, BCS(m))
    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01, mutategp=mutategparg)
    post_process(sol)

    percent_diff(sol, 'b777300ER', Nclimb)

    post_compute(sol, Nclimb)

    return sol

def run_optimal_D8(objective, fixedBPR, pRatOpt, mutategparg):
    # User definitions
    Nclimb = 3
    Ncruise = 2
    Nmission = 1
    aircraft = 'optimalD8'

    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
    
    substitutions = get_optimal_D8_subs()

    if Nmission == 2:
        substitutions.update({
            'R_{req}': [3000.*units('nmi'), 2000.*units('nmi')], #,2500.*units('nmi')
            'n_{pass}': [180., 180.], #,140.
        })
    else:
        substitutions.update({
            'R_{req}': 3000.*units('nmi'),
            'n_{pass}': 180.,
        })

    if fixedBPR:
        substitutions.update({
            '\\alpha_{max}': 6.97,
        })

    if pRatOpt:
        #del substitutions['\pi_{f_D}']
        del substitutions['\pi_{lc_D}']
        del substitutions['\pi_{hc_D}']

    m = Model(m.cost, BCS(m))
    m.substitutions.update(substitutions)

    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01, mutategp=mutategparg, skipsweepfailures=True)
    post_process(sol)

    percent_diff(sol, 'D82', Nclimb)

    post_compute(sol, Nclimb)

    return sol

def run_optimal_737(objective, fixedBPR, pRatOpt, mutategparg):
    # User definitions
    Nclimb = 3
    Ncruise = 2
    Nmission = 1
    aircraft = 'optimal737'

    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
    
    substitutions = get737_optimal_subs()

    if Nmission == 2:
        substitutions.update({
            'R_{req}': [3000.*units('nmi'),2000.*units('nmi')], #,2500.*units('nmi')
            'n_{pass}': [180., 180.], #, 140.
        })
    else:
        substitutions.update({
           'R_{req}': 3000.*units('nmi'),
           'n_{pass}': 180.,
        })

    if fixedBPR:
        substitutions.update({
            '\\alpha_{max}': 6.97,
        })
        
    if pRatOpt:
        #del substitutions['\pi_{f_D}']
        del substitutions['\pi_{lc_D}']
        del substitutions['\pi_{hc_D}']

    m.substitutions.update(substitutions)

    m = Model(m.cost, BCS(m))
    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01, mutategp=mutategparg)
    post_process(sol)

    percent_diff(sol, 'b737800', Nclimb)

    post_compute(sol, Nclimb)

    return sol

def run_M072_737(objective, fixedBPR, pRatOpt, mutategparg):
    # User definitions
    Nclimb = 3
    Ncruise = 2
    Nmission = 1
    aircraft = 'M072_737'

    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
    
    substitutions = get_M072_737_subs()

    if Nmission > 1:
        substitutions.update({
            'R_{req}': [3000.*units('nmi'),2500.*units('nmi')],
            'n_{pass}': [180., 180.],
        })
    else:
        substitutions.update({
            'R_{req}': 3000.*units('nmi'),
            'n_{pass}': 180.,
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

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01, mutategp=mutategparg)
    post_process(sol)

    percent_diff(sol, aircraft, Nclimb)

    post_compute(sol, Nclimb)

    return sol

def run_D8_eng_wing(objective, fixedBPR, pRatOpt, mutategparg):
    # User definitions
    Nclimb = 3
    Ncruise = 2
    Nmission = 1
    aircraft = 'D8_eng_wing'

    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
    
    substitutions = get_D8_eng_wing_subs()

    if Nmission > 1:
        substitutions.update({
            'R_{req}': [3000.*units('nmi'),2500.*units('nmi')],
            'n_{pass}': [180., 180.],
        })
    else:
        substitutions.update({
            'n_{pass}': 180.,
            'R_{req}': 3000.*units('nmi'),
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

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01, mutategp=mutategparg)
    post_process(sol)

    percent_diff(sol, 'b737800', Nclimb)

    post_compute(sol, Nclimb)

    return sol

def run_D8_no_BLI(objective, fixedBPR, pRatOpt, mutategparg):
    # User definitions
    Nclimb = 3
    Ncruise = 2
    Nmission = 1
    aircraft = 'D8_no_BLI'

    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
    
    substitutions = get_D8_no_BLI_subs()

    if Nmission > 1:
        substitutions.update({
            'R_{req}': [3000.*units('nmi'),2500.*units('nmi')],
            'n_{pass}': [180., 180.],
        })
    else:
        substitutions.update({
            'R_{req}': 3000.*units('nmi'),
            'n_{pass}': 180.,
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

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01, mutategp=mutategparg)
    post_process(sol)

    percent_diff(sol, 'D82', Nclimb)

    post_compute(sol, Nclimb)

    return sol

def test():
    run_optimal_737('W_{f_{total}}', True, False, True)

##    OLD CODE DOESN'T CURRENTLY WORK
    
##def run_737800(objective = 'W_{f_{total}}'):
##    # User definitions
##    Nclimb = 3
##    Ncruise = 2
##    Nmission = 1
##    aircraft = 'b737800'
##
##    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
##    
##    substitutions = getb737800subs()
##
##    substitutions.update({
###                'n_{pass}': 180.,
##        'R_{req}': 3000.*units('nmi'),
##    })
##
##    m.substitutions.update(substitutions)
##
##    m = Model(m.cost, BCS(m))
##    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])
##
##    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
##
##    percent_diff(sol, aircraft, Nclimb)
##
##    post_compute(sol, Nclimb)
##
##    return sol
##
##def run_D82(objective = 'W_{f_{total}}'):
##    # User definitions
##    Nclimb = 3
##    Ncruise = 2
##    Nmission = 1
##    aircraft = 'D82'
##
##    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
##    
##    substitutions = getD82subs()
##
##    substitutions.update({
###                'n_{pass}': 180.,
##        'R_{req}': 3000.*units('nmi'),
##    })
##
##    m.substitutions.update(substitutions)
##    m = Model(m.cost, BCS(m))
##    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])
##
##    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
##    post_process(sol)
##
##    percent_diff(sol, aircraft, Nclimb)
##
##    post_compute(sol, Nclimb)
##
##    return sol
##
##def run_D12(objective, fixedBPR, pRatOpt, mutategparg):
##    # User definitions
##    Nclimb = 3
##    Ncruise = 2
##    Nmission = 1
##    aircraft = 'D12'
##
##    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
##    
##    substitutions = get_D12_subs()
##
##    substitutions.update({
###                'n_{pass}': 500.,
##        'R_{req}': 5600.*units('nmi'),
##    })
##
##    if fixedBPR:
##        substitutions.update({
##            '\\alpha_{max}': 8.009,
##        })
##
##    if pRatOpt:
##        del substitutions['\pi_{f_D}']
####        del substitutions['\pi_{lc_D}']
####        del substitutions['\pi_{hc_D}']
##
##    m.substitutions.update(substitutions)
##    m = Model(m.cost, BCS(m))
##    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])
##
##    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
##    post_process(sol)
##
##    percent_diff(sol, aircraft, Nclimb)
##
##    post_compute(sol, Nclimb)
##
##    return sol
##

##
    

##
##def run_M08_D8_eng_wing(objective, fixedBPR, pRatOpt, mutategparg):
##    # User definitions
##    Nclimb = 3
##    Ncruise = 2
##    Nmission = 1
##    aircraft = 'M08_D8_eng_wing'
##
##    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
##    
##    substitutions = getM08_D8_eng_wing_subs()
##
##    substitutions.update({
###                'n_{pass}': 180.,
##        'R_{req}': 3000.*units('nmi'),
##    })
##
##    if fixedBPR:
##        substitutions.update({
##            '\\alpha_{max}': 6.97,
##        })
##
##    if pRatOpt:
##        del substitutions['\pi_{f_D}']
####        del substitutions['\pi_{lc_D}']
####        del substitutions['\pi_{hc_D}']
##
##    m.substitutions.update(substitutions)
##
##    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])
##
##    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
##    post_process(sol)
##
##    percent_diff(sol, 'D82', Nclimb)
##
##    post_compute(sol, Nclimb)
##
##    return sol
##
##def run_M08_D8(objective, fixedBPR, pRatOpt, mutategparg):
##    # User definitions
##    Nclimb = 3
##    Ncruise = 2
##    Nmission = 1
##    aircraft = 'M08D8'
##
##    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
##    
##    substitutions = subs_M08_D8()
##
##    substitutions.update({
###                'n_{pass}': 180.,
##        'R_{req}': 3000.*units('nmi'),
##    })
##
##    if fixedBPR:
##        substitutions.update({
##            '\\alpha_{max}': 6.97,
##        })
##
##    if pRatOpt:
##        del substitutions['\pi_{f_D}']
####        del substitutions['\pi_{lc_D}']
####        del substitutions['\pi_{hc_D}']
##        
##    m.substitutions.update(substitutions)
##    m = Model(m.cost, BCS(m))
##    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])
##
##    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
##    post_process(sol)
##
##    percent_diff(sol, 'D82', Nclimb)
##
##    post_compute(sol, Nclimb)
##
##    return sol
##
##def run_M08_D8_no_BLI(objective, fixedBPR, pRatOpt, mutategparg):
##    # User definitions
##    Nclimb = 3
##    Ncruise = 2
##    Nmission = 1
##    aircraft = 'M08D8_noBLI'
##
##    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
##    
##    substitutions = get_subs_M08_D8_noBLI()
##
##    substitutions.update({
###                'n_{pass}': 180.,
##        'R_{req}': 3000.*units('nmi'),
##    })
##
##    if fixedBPR:
##        substitutions.update({
##            '\\alpha_{max}': 6.97,
##        })
##
##    if pRatOpt:
##        del substitutions['\pi_{f_D}']
####        del substitutions['\pi_{lc_D}']
####        del substitutions['\pi_{hc_D}']
##        
##    m.substitutions.update(substitutions)
##    m = Model(m.cost, BCS(m))
##    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])
##
##    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
##    post_process(sol)
##
##    percent_diff(sol, 'D82', Nclimb)
##
##    post_compute(sol, Nclimb)
##
##    return sol
##
##def run_777300ER(objective = 'W_{f_{total}}'):
##    # User definitions
##    Nclimb = 3
##    Ncruise = 2
##    Nmission = 1
##    aircraft = 'b777300ER'
##
##    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
##    
##    substitutions = getb777300ERsubs()
##
##    substitutions.update({
####        'n_{pass}': [450.],
##        'R_{req}': 6000.*units('nmi'),
##    })
##
##    m.substitutions.update(substitutions)
##    m = Model(m.cost, BCS(m))
##    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])
##
##    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
##    post_process(sol)
##
##    percent_diff(sol, aircraft, Nclimb)
##
##    post_compute(sol, Nclimb)
##
##    return sol
##
##
##
##def run_M08_optimal_777(objective, fixedBPR, pRatOpt, mutategparg):
##    # User definitions
##    Nclimb = 3
##    Ncruise = 2
##    Nmission = 1
##    aircraft = 'optimal777_M08'
##
##    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
##    
##    substitutions = get_optimal_777300ER_M08_subs()
##
##    substitutions.update({
####        'n_{pass}': [450.],
##        'R_{req}': [6000.],
##    })
##
##    if fixedBPR:
##        substitutions.update({
##            '\\alpha_{max}': 8.62, 
##        })
##
##    if pRatOpt:
##        del substitutions['\pi_{f_D}']
####        del substitutions['\pi_{lc_D}']
####        del substitutions['\pi_{hc_D}']
##
##    m.substitutions.update(substitutions)
##    m = Model(m.cost, BCS(m))
##    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])
##
##    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
##    post_process(sol)
##
##    percent_diff(sol, 'b777300ER', Nclimb)
##
##    post_compute(sol, Nclimb)
##
##    return sol
##
##def run_M072_optimal_777(objective, fixedBPR, pRatOpt, mutategparg):
##    # User definitions
##    Nclimb = 3
##    Ncruise = 2
##    Nmission = 1
##    aircraft = 'optimal777_M072'
##
##    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
##    
##    substitutions = get_optimal_777300ER_M072_subs()
##
##    substitutions.update({
####        'n_{pass}': [450.],
##        'R_{req}': [6000.],
##    })
##
##    if fixedBPR:
##        substitutions.update({
##            '\\alpha_{max}': 8.62, 
##        })
##
##    if pRatOpt:
##        del substitutions['\pi_{f_D}']
####        del substitutions['\pi_{lc_D}']
####        del substitutions['\pi_{hc_D}']
##
##    m.substitutions.update(substitutions)
##    m = Model(m.cost, BCS(m))
##    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])
##
##    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
##    post_process(sol)
##
##    percent_diff(sol, 'b777300ER', Nclimb)
##
##    post_compute(sol, Nclimb)
##
##    return sol
##
##def run_D8_big(objective, fixedBPR, pRatOpt, mutategparg):
##    # User definitions
##    Nclimb = 3
##    Ncruise = 2
##    Nmission = 1
##    aircraft = 'D8big'
##
##    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
##    
##    substitutions = getD8bigsubs()
##
##    substitutions.update({
###                'n_{pass}': 180.,
##        'R_{req}': 6000.*units('nmi'),
##    })
##
##    if fixedBPR:
##        substitutions.update({
##            '\\alpha_{max}': 8.62, 
##        })
##
##    if pRatOpt:
##        del substitutions['\pi_{f_D}']
####        del substitutions['\pi_{lc_D}']
####        del substitutions['\pi_{hc_D}']
##        
##    m.substitutions.update(substitutions)
##    m = Model(m.cost, BCS(m))
##    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])
##
##    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
##    post_process(sol)
##
##    percent_diff(sol, 'b777300ER', Nclimb)
##
##    post_compute(sol, Nclimb)
##
##    return sol
##
##def run_D8_big_no_BLI(objective, fixedBPR, pRatOpt, mutategparg):
##    # User definitions
##    Nclimb = 3
##    Ncruise = 2
##    Nmission = 1
##    aircraft = 'D8big_no_BLI'
##
##    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
##    
##    substitutions = getD8big_noBLI_subs()
##
##    substitutions.update({
###                'n_{pass}': 180.,
##        'R_{req}': 6000.*units('nmi'),
##    })
##
##    if fixedBPR:
##        substitutions.update({
##            '\\alpha_{max}': 8.62, 
##        })
##
##    if pRatOpt:
##        del substitutions['\pi_{f_D}']
####        del substitutions['\pi_{lc_D}']
####        del substitutions['\pi_{hc_D}']
##        
##    m.substitutions.update(substitutions)
##    m = Model(m.cost, BCS(m))
##    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])
##
##    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
##    post_process(sol)
##
##    percent_diff(sol, 'b777300ER', Nclimb)
##
##    post_compute(sol, Nclimb)
##
##    return sol
##
##def run_M072_D8_big_no_BLI(objective, fixedBPR, pRatOpt, mutategparg):
##    # User definitions
##    Nclimb = 3
##    Ncruise = 2
##    Nmission = 1
##    aircraft = 'D8big_no_BLI'
##
##    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
##    
##    substitutions = getD8big_M072_noBLI_subs()
##
##    substitutions.update({
###                'n_{pass}': 180.,
##        'R_{req}': 6000.*units('nmi'),
##    })
##
##    if fixedBPR:
##        substitutions.update({
##            '\\alpha_{max}': 8.62, 
##        })
##
##    if pRatOpt:
##        del substitutions['\pi_{f_D}']
####        del substitutions['\pi_{lc_D}']
####        del substitutions['\pi_{hc_D}']
##        
##    m.substitutions.update(substitutions)
##    m = Model(m.cost, BCS(m))
##    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])
##
##    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
##    post_process(sol)
##
##    percent_diff(sol, 'b777300ER', Nclimb)
##
##    post_compute(sol, Nclimb)
##
##    return sol
##
##def run_D8_big_eng_wing(objective, fixedBPR, pRatOpt, mutategparg):
##    # User definitions
##    Nclimb = 3
##    Ncruise = 2
##    Nmission = 1
##    aircraft = 'D8big_eng_wing'
##
##    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
##    
##    substitutions = getD8big_eng_wing_subs()
##
##    substitutions.update({
###                'n_{pass}': 180.,
##        'R_{req}': 6000.*units('nmi'),
##    })
##
##    if fixedBPR:
##        substitutions.update({
##            '\\alpha_{max}': 8.62, 
##        })
##
##    if pRatOpt:
##        del substitutions['\pi_{f_D}']
####        del substitutions['\pi_{lc_D}']
####        del substitutions['\pi_{hc_D}']
##        
##    m.substitutions.update(substitutions)
##    m = Model(m.cost, BCS(m))
##    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])
##
##    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
##    post_process(sol)
##
##    percent_diff(sol, 'b777300ER', Nclimb)
##
##    post_compute(sol, Nclimb)
##
##    return sol
##
##def run_M072_D8_big_eng_wing(objective, fixedBPR, pRatOpt, mutategparg):
##    # User definitions
##    Nclimb = 3
##    Ncruise = 2
##    Nmission = 1
##    aircraft = 'D8big_eng_wing_M072'
##
##    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
##    
##    substitutions = getD8big_M072_eng_wing_subs()
##
##    substitutions.update({
###                'n_{pass}': 180.,
##        'R_{req}': 6000.*units('nmi'),
##    })
##
##    if fixedBPR:
##        substitutions.update({
##            '\\alpha_{max}': 8.62, 
##        })
##
##    if pRatOpt:
##        del substitutions['\pi_{f_D}']
####        del substitutions['\pi_{lc_D}']
####        del substitutions['\pi_{hc_D}']
##        
##    m.substitutions.update(substitutions)
##    m = Model(m.cost, BCS(m))
##    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])
##
##    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
##    post_process(sol)
##
##    percent_diff(sol, 'b777300ER', Nclimb)
##
##    post_compute(sol, Nclimb)
##
##    return sol
##
##
##def run_D8_big_M072(objective, fixedBPR, pRatOpt, mutategparg):
##    # User definitions
##    Nclimb = 3
##    Ncruise = 2
##    Nmission = 1
##    aircraft = 'D8big'
##
##    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
##    
##    substitutions = getD8big_M072_subs()
##
##    substitutions.update({
###                'n_{pass}': 180.,
##        'R_{req}': 6000.*units('nmi'),
##    })
##
##    if fixedBPR:
##        substitutions.update({
##            '\\alpha_{max}': 8.62, 
##        })
##
##    if pRatOpt:
##        del substitutions['\pi_{f_D}']
####        del substitutions['\pi_{lc_D}']
####        del substitutions['\pi_{hc_D}']
##        
##    m.substitutions.update(substitutions)
##    m = Model(m.cost, BCS(m))
##    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])
##
##    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
##    post_process(sol)
##
##    percent_diff(sol, 'b777300ER', Nclimb)
##
##    post_compute(sol, Nclimb)
##
##    return sol
##
##def run_D8_big_M08(objective, fixedBPR, pRatOpt, mutategparg):
##    # User definitions
##    Nclimb = 3
##    Ncruise = 2
##    Nmission = 1
##    aircraft = 'D8big'
##
##    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
##    
##    substitutions = getD8big_M08_subs()
##
##    substitutions.update({
###                'n_{pass}': 180.,
##        'R_{req}': 6000.*units('nmi'),
##    })
##
##    if fixedBPR:
##        substitutions.update({
##            '\\alpha_{max}': 8.62, 
##        })
##
##    if pRatOpt:
##        del substitutions['\pi_{f_D}']
####        del substitutions['\pi_{lc_D}']
####        del substitutions['\pi_{hc_D}']
##        
##    m.substitutions.update(substitutions)
##    m = Model(m.cost, BCS(m))
##    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])
##
##    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
##    post_process(sol)
##
##    percent_diff(sol, 'b777300ER', Nclimb)
##
##    post_compute(sol, Nclimb)
##
##    return sol
##
##def run_optimal_RJ(objective, fixedBPR, pRatOpt, mutategparg):
##    # User definitions
##    Nclimb = 3
##    Ncruise = 2
##    Nmission = 1
##    aircraft = 'optimalRJ'
##
##    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
##    
##    substitutions = get_optimal_RJ_subs()
##
##    substitutions.update({
####        'n_{pass}': [90.],
##        'R_{req}': [1500.],
##    })
##
##    if fixedBPR:
##        substitutions.update({
##            '\\alpha_{max}': 5.6958, 
##        })
##
##    if pRatOpt:
##        del substitutions['\pi_{f_D}']
####        del substitutions['\pi_{lc_D}']
####        del substitutions['\pi_{hc_D}']
##
##    m.substitutions.update(substitutions)
##    m = Model(m.cost, BCS(m))
##    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])
##
##    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
##    post_process(sol)
##
##    percent_diff(sol, 'b737800', Nclimb)
##
##    post_compute(sol, Nclimb)
##
##    return sol
##
##def run_M072_optimal_RJ(objective, fixedBPR, pRatOpt, mutategparg):
##    # User definitions
##    Nclimb = 3
##    Ncruise = 2
##    Nmission = 1
##    aircraft = 'optimalRJ'
##
##    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
##    
##    substitutions = get_M072_optimal_RJ_subs()
##
##    substitutions.update({
####        'n_{pass}': [90.],
##        'R_{req}': [1500.],
##    })
##
##    if fixedBPR:
##        substitutions.update({
##            '\\alpha_{max}': 5.6958, 
##        })
##
##    if pRatOpt:
##        del substitutions['\pi_{f_D}']
####        del substitutions['\pi_{lc_D}']
####        del substitutions['\pi_{hc_D}']
##
##    m.substitutions.update(substitutions)
##    m = Model(m.cost, BCS(m))
##    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])
##
##    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
##    post_process(sol)
##
##    percent_diff(sol, 'b737800', Nclimb)
##
##    post_compute(sol, Nclimb)
##
##    return sol
##
##def run_small_D8(objective, fixedBPR, pRatOpt, mutategparg):
##    # User definitions
##    Nclimb = 3
##    Ncruise = 2
##    Nmission = 1
##    aircraft = 'smallD8'
##
##    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
##    
##    substitutions = get_small_D8_subs()
##
##    substitutions.update({
####        'n_{pass}': [90.],
##        'R_{req}': [1500.],
##    })
##
##    if fixedBPR:
##        substitutions.update({
##            '\\alpha_{max}': 5.6958, 
##        })
##
##    if pRatOpt:
##        del substitutions['\pi_{f_D}']
####        del substitutions['\pi_{lc_D}']
####        del substitutions['\pi_{hc_D}']
##
##    m.substitutions.update(substitutions)
##    m = Model(m.cost, BCS(m))
##    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])
##
##    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
##    post_process(sol)
##
##    percent_diff(sol, 'b737800', Nclimb)
##
##    post_compute(sol, Nclimb)
##
##    return sol
##
##def run_small_M08_D8(objective, fixedBPR, pRatOpt, mutategparg):
##    # User definitions
##    Nclimb = 3
##    Ncruise = 2
##    Nmission = 1
##    aircraft = 'smallD8'
##
##    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
##    
##    substitutions = get_small_M08_D8_subs()
##
##    substitutions.update({
####        'n_{pass}': [90.],
##        'R_{req}': [1500.],
##    })
##
##    if fixedBPR:
##        substitutions.update({
##            '\\alpha_{max}': 5.6958, 
##        })
##
##    if pRatOpt:
##        del substitutions['\pi_{f_D}']
####        del substitutions['\pi_{lc_D}']
####        del substitutions['\pi_{hc_D}']
##
##    m.substitutions.update(substitutions)
##    m = Model(m.cost, BCS(m))
##    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])
##
##    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
##    post_process(sol)
##
##    percent_diff(sol, 'b737800', Nclimb)
##
##    return sol
##
##def run_small_D8_eng_wing(objective, fixedBPR, pRatOpt, mutategparg):
##    # User definitions
##    Nclimb = 3
##    Ncruise = 2
##    Nmission = 1
##    aircraft = 'smallD8_eng_wing'
##
##    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
##    
##    substitutions = get_small_D8_eng_wing_subs()
##
##    substitutions.update({
####        'n_{pass}': [90.],
##        'R_{req}': [1500.],
##    })
##
##    if fixedBPR:
##        substitutions.update({
##            '\\alpha_{max}': 5.6958, 
##        })
##
##    if pRatOpt:
##        del substitutions['\pi_{f_D}']
####        del substitutions['\pi_{lc_D}']
####        del substitutions['\pi_{hc_D}']
##
##    m.substitutions.update(substitutions)
##    m = Model(m.cost, BCS(m))
##    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])
##
##    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
##    post_process(sol)
##
##    percent_diff(sol, 'b737800', Nclimb)
##
##    post_compute(sol, Nclimb)
##
##    return sol
##
##def run_small_M08_D8_eng_wing(objective, fixedBPR, pRatOpt, mutategparg):
##    # User definitions
##    Nclimb = 3
##    Ncruise = 2
##    Nmission = 1
##    aircraft = 'smallD8_M08_eng_wing'
##
##    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
##    
##    substitutions = get_small_M08_D8_eng_wing_subs()
##
##    substitutions.update({
####        'n_{pass}': [90.],
##        'R_{req}': [1500.],
##    })
##
##    if fixedBPR:
##        substitutions.update({
##            '\\alpha_{max}': 5.6958, 
##        })
##
##    if pRatOpt:
##        del substitutions['\pi_{f_D}']
####        del substitutions['\pi_{lc_D}']
####        del substitutions['\pi_{hc_D}']
##
##    m.substitutions.update(substitutions)
##    m = Model(m.cost, BCS(m))
##    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])
##
##    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
##    post_process(sol)
##
##    percent_diff(sol, 'b737800', Nclimb)
##
##    post_compute(sol, Nclimb)
##
##    return sol
##
##def run_small_D8_no_BLI(objective, fixedBPR, pRatOpt, mutategparg):
##    # User definitions
##    Nclimb = 3
##    Ncruise = 2
##    Nmission = 1
##    aircraft = 'smallD8_no_BLI'
##
##    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
##    
##    substitutions = get_small_D8_no_BLI_subs()
##
##    substitutions.update({
####        'n_{pass}': [90.],
##        'R_{req}': [1500.],
##    })
##
##    if fixedBPR:
##        substitutions.update({
##            '\\alpha_{max}': 5.6958, 
##        })
##
##    if pRatOpt:
##        del substitutions['\pi_{f_D}']
####        del substitutions['\pi_{lc_D}']
####        del substitutions['\pi_{hc_D}']
##
##    m.substitutions.update(substitutions)
##    m = Model(m.cost, BCS(m))
##    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])
##
##    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
##    post_process(sol)
##
##    percent_diff(sol, 'b737800', Nclimb)
##
##    post_compute(sol, Nclimb)
##
##    return sol
##
##def run_small_M08_D8_no_BLI(objective, fixedBPR, pRatOpt, mutategparg):
##    # User definitions
##    Nclimb = 3
##    Ncruise = 2
##    Nmission = 1
##    aircraft = 'smallD8_M08_no_BLI'
##
##    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
##    
##    substitutions = get_small_M08_D8_no_BLI_subs()
##
##    substitutions.update({
####        'n_{pass}': [90.],
##        'R_{req}': [1500.],
##    })
##
##    if fixedBPR:
##        substitutions.update({
##            '\\alpha_{max}': 5.6958, 
##        })
##
##    if pRatOpt:
##        del substitutions['\pi_{f_D}']
####        del substitutions['\pi_{lc_D}']
####        del substitutions['\pi_{hc_D}']
##
##    m.substitutions.update(substitutions)
##    m = Model(m.cost, BCS(m))
##    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])
##
##    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
##    post_process(sol)
##
##    percent_diff(sol, 'b737800', Nclimb)
##
##    post_compute(sol, Nclimb)
##
##    return sol
