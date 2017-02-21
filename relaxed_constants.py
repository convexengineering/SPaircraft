from gpkit.constraints.relax import ConstantsRelaxed
from gpkit import Model

"""
Methods to precondition an SP so that it solves with a relaxed constants algorithim
and postcondition an SP to ensure all relax values are 1
"""

def relaxed_constants(model):
    """
    Method to precondition an SP so it solves with a relaxed constants algorithim
    """

    if model.substitutions:
        constsrelaxed = ConstantsRelaxed(model)
        feas = Model(constsrelaxed.relaxvars.prod()**20 * model.cost,
                     constsrelaxed)
        # NOTE: It hasn't yet been seen but might be possible that
        #       the model.cost component above could cause infeasibility
    else:
        feas = Model(model.cost, model)

    return feas
