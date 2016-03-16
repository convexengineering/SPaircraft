### NPV MODEL ###

from gpkit import Model, Variable, units, VectorVariable
from gpkit.tools import te_exp_minus1
import gpkit
import numpy as np

class NPV(Model):
    def __init__(self):

        N = 3 # number of payments
        NPV = Variable('NPV', 'USD', 'net present value')
        PV = Variable('PV', 'USD', 'present value')
        C = VectorVariable(N, 'C', 'USD', 'cash flow')
        C_r = Variable('C_r', 4e6, 'USD/years', 'cash flow rate')
        r = Variable('r', 0.1, '1/years', 'interest rate')
        t = VectorVariable(N, 't', 'years', 'time')
       
        constraints = [N*PV >= NPV, 
                       C >= PV*(te_exp_minus1(r*t,4) + 1),
                       t >= t.left + C/C_r]

        cost = 1/NPV

        Model.__init__(self, cost, constraints)

    def test(self):
        seft.solve()

if __name__ == "__main__":
    M = NPV()
    sol = M.solve()
