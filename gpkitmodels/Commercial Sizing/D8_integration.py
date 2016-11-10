"""temporary file containing classes in flux for full aircraft integration"""
from numpy import pi, tan, cos
import numpy as np
from gpkit import Variable, Model, units, SignomialsEnabled, VectorVariable
from gpkit.constraints.sigeq import SignomialEqualityConstraint as SignomialEquality
from gpkit.constraints.tight import TightConstraintSet as TCS

class Engine(Model):
    """
    place holder engine model
    """
    def __init__(self, **kwargs):
        #new variables
        W_engine = Variable('W_{engine}', 'N', 'Weight of a Single Turbofan Engine')
        A2 = Variable('A_2', 'm^2', 'Fan Area')
        
        constraints = []

        constraints.extend([
            W_engine == 1000 * units('N'),

            A2 == A2,
            ])

        Model.__init__(self, None, constraints)

    def dynamic(self, state):
            """
            returns an engine performance model
            """
            return EnginePerformance(self, state)

class EnginePerformance(Model):
    """
    place holder engine perofrmacne model
    """
    def __init__(self, engine, state, **kwargs):
        #new variables
        TSFC = Variable('TSFC', '1/hr', 'Thrust Specific Fuel Consumption')
        thrust = Variable('thrust', 'N', 'Thrust')
        
        #constraints
        constraints = []

        constraints.extend([
            TSFC == TSFC,

            thrust == thrust, #want thrust to enter the model
            ])

        Model.__init__(self, None, constraints)


class Fuselage(Model):
    """
    place holder fuselage model
    """
    def __init__(self, **kwargs):
        #new variables
        npax = Variable('n_{pax}', '-', 'Number of Passengers to Carry')
                           
        #weight variables
        Wpay = Variable('W_{payload}', 'N', 'Aircraft Payload Weight')
        W_e = Variable('W_{e}', 'N', 'Empty Weight of Aircraft')
        Wpax = Variable('W_{pax}', 'N', 'Estimated Average Passenger Weight, Includes Baggage')

        Afuse = Variable('A_{fuse}', 'm^2', 'Estimated Fuselage Area')
        pax_area = Variable('pax_{area}', 'm^2', 'Estimated Fuselage Area per Passenger')

        lfuse   = Variable('l_{fuse}', 'm', 'Fuselage length')
        wfuse   = Variable('w_{fuse}', 6, 'm', 'Fuselage width')

        constraints = []
        
        constraints.extend([
            #compute fuselage area for drag approximation
            Afuse == pax_area * npax,

            Afuse == lfuse * wfuse,
            
            #constraints on the various weights
            Wpay == npax * Wpax,
            
            #estimate based on TASOPT 737 model
            W_e == .75*Wpay,
            ])

        Model.__init__(self, None, constraints)

    def dynamic(self, state):
        """
        returns a fuselage performance model
        """
        return FuselagePerformance(self, state)

class FuselagePerformance(Model):
    """
    Fuselage performance model
    """
    def __init__(self, fuse, state, **kwargs):
        #new variables
        Cdfuse = Variable('C_{D_{fuse}}', '-', 'Fuselage Drag Coefficient')
        Dfuse = Variable('D_{fuse}', 'N', 'Total Fuselage Drag')
        Cmfu    = Variable('C_{m_{fuse}}', '-', 'Moment coefficient (fuselage)')

        xcg     = VectorVariable(2, 'x_{CG}', 'm', 'CG location')
        
        #constraints
        constraints = []

        constraints.extend([
            Dfuse == Cdfuse * (.5 * fuse['A_{fuse}'] * state.atm['\\rho'] * state['V']**2),

            Cdfuse == .005,

            Cmfu == .05,

            xcg == 17 * units('m'),
            ])

        Model.__init__(self, None, constraints)
