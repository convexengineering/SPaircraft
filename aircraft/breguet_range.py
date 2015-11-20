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
        W_oew = Var('W_{oew}', 7000, "lbf", "Operating Empty Weight")
        LoverD_max = Var('LoverD_{max}', 15, "-", "Maximum Lift to Drag Ratio")
        M_max = Var('M_max', 0.78, "-", "Maximum Mach")

        #Constants
        g = Var('g', 9.81, "m/s^2", "gravity")
        a0 = Var('a0', 340.29, "m/s", "speed of sound at sea level")

        #Free Variables
        R = Var('R', "nautical_miles", "range")
        M = Var('M', "-", "Mach number")
        LoverD = Var('LoverD', "-", "life to drag ratio")
        TSFC = Var('TSFC', "lb/lbf/hr", "thrust specific fuel consuption")
        W_init = Var('W_{init}', "lbf", "initial weight")
        W_fuel = Var('W_{fuel}', "lbf", "fuel weight")
        z_bre = Var('z_{bre}', "-", "Breguet parameter")
        t = Var('t', "hr", "time")

        # Set up Model Equations
        objective = 1/R  # Maximize range
        constraints = [W_init >= W_oew + W_fuel,
                       W_init <= MTOW,
                       LoverD <= LoverD_max,
                       TSFC >= TSFC_min,
                       M <= M_max,
                       t >= R/M/a0,
                       z_bre >= t*TSFC*g/LoverD,
                       W_fuel/W_oew >= te_exp_minus1(z_bre, nterm=3)
                       ]
        return objective, constraints
    
    def test(self):
        sol = self.solve()
        R = sol('R')
        M = sol('M')
        LoverD = sol('LoverD')
        TSFC = sol('TSFC')
        W_init = sol('W_{init}')
        W_fuel = sol('W_{fuel}')
        z_bre = sol('z_{bre}')
        t = sol('t')
        g = sol('g')
        a0 = sol('a0')

if __name__ == "__main__":
    BreguetRange().test()
