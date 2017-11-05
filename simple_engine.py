"""
simple engine for pkjournal model
"""
from gpkit import Variable, Model, Vectorize

class SimpleEngine(Model):
    """
    most simple engine model imaginable
    """
    def setup(self, Nsegments):
        A2 = Variable('A_{2}', 'm^2', 'Fan Area')
        with Vectorize(Nsegments):
            TSFC = Variable('TSFC', '1/hr', 'Thrust Specific Fuel Consumption')
            F = Variable('F', 'N', 'Cruise Thrust')
