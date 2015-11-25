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

        TSFC_min = Var('TSFC_{min}', 0.307, "lb/lbf/hr", "Minimum TSFC")
        MTOW = Var('MTOW', 10000, "lbf", "Max take off weight")
        W_e = Var('W_{e}', 7000, "lbf", "Operating Empty Weight")
        LoverD_max = Var('LoverD_{max}', 15, "-", "Maximum Lift to Drag Ratio")
        V_max = Var('V_{max}', 420, "knots", "Maximum Velocity")

        #Constants
        g = Var('g', 9.81, "m/s^2", "gravity")

        #Free Variables
        R = Var('R', "nautical_miles", "range")
        V = Var('V', "knots", "Velocity")
        LoverD = Var('LoverD', "-", "life to drag ratio")
        TSFC = Var('TSFC', "lb/lbf/hr", "thrust specific fuel consuption")
        W_init = Var('W_{init}', "lbf", "initial weight")
        W_fuel = Var('W_{fuel}', "lbf", "fuel weight")
        z_bre = Var('z_{bre}', "-", "Breguet parameter")
        t = Var('t', "hr", "time")

        # Set up Model Equations
        objective = 1/R  # Maximize range
        constraints = [W_init >= W_e + W_fuel,
                       W_init <= MTOW,
                       LoverD <= LoverD_max,
                       TSFC >= TSFC_min,
                       V <= V_max,
                       t >= R/V,
                       z_bre >= t*TSFC*g/LoverD,
                       W_fuel/W_e >= te_exp_minus1(z_bre, nterm=3)
                      ]
        return objective, constraints

    def test(self):
        _ = self.solve()


if __name__ == "__main__":
    BreguetRange().test()
