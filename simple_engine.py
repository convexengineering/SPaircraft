"""
simple engine for pkjournal model
"""
from gpkit import Variable, Model, units, Vectorize

class Engine(Model):
    """
    most simple engine model imaginable
    """
    def setup(self):
        TSFC = Variable('TSFC', '1/hr', 'Thrust Specific Fuel Consumption')
        A2 = Variable('A_{2}', 'm^2', 'Fan Area')
        F = Variable('F', 'N', 'Cruise Thrust')
        F_TO = Variable('F_TO', 'N', 'Takeoff Thrust for VT Sizing')

