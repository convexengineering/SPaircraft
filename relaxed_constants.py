from gpkit.constraints.relax import ConstantsRelaxed
from gpkit import Model

"""
Methods to precondition an SP so that it solves with a relaxed constants algorithim
and postcondition an SP to ensure all relax values are 1
"""

def relaxed_constants(model, include_only=None, exclude=None):
    """
    Method to precondition an SP so it solves with a relaxed constants algorithim

    ARGUMENTS
    ---------
    model: the model to solve with relaxed constants

    RETURNS
    -------
    feas: the input model but with relaxed constants and a new objective
    """

    if model.substitutions:
        constsrelaxed = ConstantsRelaxed(model, include_only, exclude)
        feas = Model(constsrelaxed.relaxvars.prod()**20 * model.cost,
                     constsrelaxed)
        # NOTE: It hasn't yet been seen but might be possible that
        #       the model.cost component above could cause infeasibility
    else:
        feas = Model(model.cost, model)

    return feas

def post_process(sol):
    """
    Model to print relevant info for a solved model with relaxed constants
    
    ARGUMENTS
    --------
    sol: the solution to the solved model
    """

    for i in range(len(sol.program.gps)):
        varkeys = [k for k in sol.program.gps[i].varlocs if "Relax" in k.models and sol.program.gps[i].result(k) >= 1.00001]
        if varkeys:
            print "GP iteration %s has relaxed constants" % i
            print sol.program.gps[i].result.table(varkeys)
            if i == len(sol.program.gps) - 1:
                print  "WARNING: The final GP iteration had relaxation values greater than 1"

