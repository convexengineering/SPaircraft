### NPV MODEL ###

from gpkit import Model, Variable, units, VectorVariable
import gpkit
import numpy as np

class NPV(Model):
    def setup(self):

        N = 3 # number of payments
        NPV = Variable('NPV', 'USD', 'net present value')
        PV = Variable('PV', 'USD', 'present value')
        C = VectorVariable(N, 'C', 'USD', 'cash flow')
        C_r = Variable('C_r', 4e6, 'USD', 'cash flow rate')
        r = Variable('r', 0.02, '-', 'interest rate')
        t = VectorVariable(N, 't', '-', 'time')
       
        constraints = [3*PV >= NPV, 
                       C >= PV*(1 + r*t + r*t**2/2 ),
                       t >= t.left + C/C_r]

        cost = 1/NPV

        return cost, constraints

    def test(self):
        seft.solve()

if __name__ == "__main__":
    M = NPV()
    sol = M.solve()
