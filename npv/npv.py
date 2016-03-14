### NPV MODEL ###

from gpkit import Model, Variables, units
import gpkit
import numpy as np

class NPV(Model):
    def setup(self):

        NPV = Variable('NPV', 'USD', 'net present value')
        C = Variable('C', 'USD', 'cash flow')
        C_r = Variable('C_r', 4e6, 'USD/year', 'cash flow rate')
        r = Variable('r', 0.02 '1/year', 'interest rate')


