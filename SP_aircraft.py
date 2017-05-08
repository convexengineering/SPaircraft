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
from subs_optimal_737_fixedBPR import get737_optimal_fixedBPR_subs
from subs_optimal_737 import get737_optimal_subs
from subs_optimal_D8 import get_optimal_D8_subs
from subs_M08_D8 import subs_M08_D8
from subs_M08_d8_eng_wing import getM08_D8_eng_wing_subs
from subs_D8_eng_wing import get_D8_eng_wing_subs
from subsM072737 import get_M072_737_subs
from subs_D8_no_BLI import get_D8_no_BLI_subs
from subs_M08_D8_noBLI import get_subs_M08_D8_noBLI

from gpkit import units, Model
from gpkit import Variable, Model, units, SignomialsEnabled, SignomialEquality, Vectorize
from gpkit.constraints.bounded import Bounded as BCS

from D8 import Mission

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

def run_D8_no_BLI():
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

    m.substitutions.update(substitutions)

    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
    post_process(sol)

    percent_diff(sol, 'D82', Nclimb)

    return sol

def run_M072_737():
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

    m.substitutions.update(substitutions)

    m = Model(m.cost, BCS(m))
    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
    post_process(sol)

    percent_diff(sol, 'b737800', Nclimb)

    return sol

def run_D8_eng_wing():
    # User definitions
    Nclimb = 3
    Ncruise = 2
    Nmission = 1
    objective = 'fuel'
    aircraft = 'D8_eng_wing'

    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
    
    substitutions = getD8_eng_wing_subs()

    substitutions.update({
#                'n_{paxx}': 180.,
        'ReqRng': 3000.*units('nmi'),
    })

    m.substitutions.update(substitutions)

    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
    post_process(sol)

    percent_diff(sol, 2, Nclimb)

    return sol

def run_D8_eng_wing():
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

    m.substitutions.update(substitutions)

    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
    post_process(sol)

    percent_diff(sol, 2, Nclimb)

    return sol

def run_optimal_D8():
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

    m.substitutions.update(substitutions)

    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
    post_process(sol)

    percent_diff(sol, 'D82', Nclimb)

    return sol

def run_optimal_737():
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

    m.substitutions.update(substitutions)

    m = Model(m.cost, BCS(m))
    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
    post_process(sol)

    percent_diff(sol, 'b737800', Nclimb)

    return sol

def run_optimal_737_fixedBPR():
    # User definitions
    Nclimb = 3
    Ncruise = 2
    Nmission = 1
    objective = 'fuel'
    aircraft = 'optimal737'

    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
    
    substitutions = get737_optimal_fixedBPR_subs()

    substitutions.update({
#                'n_{paxx}': 180.,
        'ReqRng': 3000.*units('nmi'),
    })

    m.substitutions.update(substitutions)

    m = Model(m.cost, BCS(m))
    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
    post_process(sol)

    percent_diff(sol, 'b737800', Nclimb)

    return sol

def run_M08_D8_eng_wing():
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

    m.substitutions.update(substitutions)

    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
    post_process(sol)

    percent_diff(sol, 'D82', Nclimb)

    return sol

def run_M08_D8():
    # User definitions
    Nclimb = 3
    Ncruise = 2
    Nmission = 1
    objective = 'fuel'
    aircraft = 'M08D8'

    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
    
    substitutions = subs_M08_D8()

    substitutions.update({
               # 'n_{pax}': 180.,
        'ReqRng': 3000.*units('nmi'),
    })

    m.substitutions.update(substitutions)

    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
    post_process(sol)

    percent_diff(sol, 'D82', Nclimb)

    return sol

def run_M08_D8_no_BLI():
    # User definitions
    Nclimb = 3
    Ncruise = 2
    Nmission = 1
    objective = 'fuel'
    aircraft = 'M08D8_no_BLI'

    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)
    
    substitutions = get_subs_M08_D8_no_BLI()

    substitutions.update({
#                'n_{paxx}': 180.,
        'ReqRng': 3000.*units('nmi'),
    })

    m.substitutions.update(substitutions)
    m = Model(m.cost, BCS(m))
    m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
    post_process(sol)

    percent_diff(sol, 'D82', Nclimb)

    return sol

def test():
    run_737800()

if __name__ == '__main__':
    # User definitions
    Nclimb = 3
    Ncruise = 2
    Nmission = 1
    objective = 'fuel'
    aircraft = 'M08_D8_noBLI'

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

    if aircraft == 'D8big':
        print('D8big executing...')
        substitutions = getD8bigsubs()
        if Nmission == 1:
                substitutions.update({
##                 'n_{pax}': 450.,
                'ReqRng': 6000.*units('nmi'),
                })
                
        if Nmission != 1:
                substitutions.update({
                 'n_{pax}': [450.],
                'ReqRng': [6000.],
                })

    if aircraft == 'b777300ER':
           print('777-300ER executing...')
           substitutions = getb777300ERsubs()
           if Nmission == 1:
                substitutions.update({
##                 'n_{pax}': 180.,
                'ReqRng': 6000.*units('nmi'),
                })

           if Nmission != 1:
                substitutions.update({
                 'n_{pax}': [450.],
                'ReqRng': [6000.],
                })

    if aircraft == 'b737800':
        print('b737800 executing...')
        substitutions = getb737800subs()
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

    m.substitutions.update(substitutions)

    if aircraft in ['D80','D82','D8_no_BLI']:
        # m = Model(m.cost,BCS(m))
        m_relax = relaxed_constants(m, None, ['ReqRng'])
    if aircraft in ['D8big', 'D82_73eng', 'D8_eng_wing', 'optimalD8', 'M08D8', 'M08_D8_eng_wing','M08_D8_noBLI']:
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
