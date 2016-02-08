"""Breguet range model"""
from gpkit.shortcuts import Var, Model
from gpkit.tools import te_exp_minus1


class BreguetRange(Model):

    """Breguet Range Model

    Assumptions
    -----------
    Fuel burn varies linearily with thrust and is independent of velocity.
    """
    def setup(self):

        # Fixed Parameters
        LD_max   = Var('\\left(\\frac{L}{D}\\right)_{max}', 15, '-',
                       'Maximum lift-to-drag ratio')
        MTOW     = Var('MTOW', 10000, 'lbf', 'Max takeoff weight')
        TSFC_min = Var('TSFC_{min}', 0.307, 'lb/lbf/hr', 'Minimum TSFC')
        V_max    = Var('V_{max}', 420, 'knots', 'Maximum Velocity')
        W_e      = Var('W_{e}', 7000, 'lbf', 'Operating Empty Weight')

        # Constants
        g        = Var('g', 9.81, 'm/s^2', 'Gravitational acceleration')

        # Free Variables
        LD       = Var('\\frac{L}{D}', '-', 'Lift-to-drag ratio')
        R        = Var('R', 'nautical_miles', 'Range')
        t        = Var('t', 'hr', 'Flight time')
        TSFC     = Var('TSFC', 'lb/lbf/hr', 'Thrust specific fuel consumption')
        V        = Var('V', 'knots', 'Velocity')
        W_fuel   = Var('W_{fuel}', 'lbf', 'Fuel weight')
        W_init   = Var('W_{init}', 'lbf', 'Initial gross weight')
        z_bre    = Var('z_{bre}', '-', 'Breguet parameter')

        # Model
        objective = 1/R  # Maximize range
        constraints = [W_init >= W_e + W_fuel,
                       W_init <= MTOW,
                       LD <= LD_max,
                       TSFC >= TSFC_min,
                       V <= V_max,
                       t >= R/V,
                       z_bre >= t*TSFC*g/LD,
                       W_fuel/W_e >= te_exp_minus1(z_bre, nterm=3)
                      ]
        return objective, constraints

    def test(self):
        _ = self.solve()

if __name__ == '__main__':
    BreguetRange().test()
