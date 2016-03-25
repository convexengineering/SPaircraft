"Net present value model"
from gpkit import Model, Variable, VectorVariable
from gpkit.tools import te_exp_minus1

class NPV(Model):
    "Net present value"
    def __init__(self, **kwargs):

        N = 3 # number of payments
        npv = Variable('NPV', 'USD', 'net present value')
        PV = Variable('PV', 'USD', 'present value')
        C = VectorVariable(N, 'C', 'USD', 'cash flow')
        C_r = Variable('C_r', 4e6, 'USD/years', 'cash flow rate')
        r = Variable('r', 0.1, '1/years', 'interest rate')
        t = VectorVariable(N, 't', 'years', 'time')

        constraints = [N*PV >= npv,
                       C >= PV*(te_exp_minus1(r*t, 4) + 1),
                       t >= t.left + C/C_r]

        cost = 1/npv

        Model.__init__(self, cost, constraints, **kwargs)

    def test(self):
        "test the model by solving it"
        self.solve()

if __name__ == "__main__":
    M = NPV()
    SOL = M.solve()
