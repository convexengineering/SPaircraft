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
        F = Variable('F', 'N', 'Cruise Thrust')
        F_TO = Variable('F_TO', 'N', 'Takeoff Thrust for VT Sizing')

