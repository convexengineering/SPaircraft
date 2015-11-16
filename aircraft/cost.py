"""Aircraft cost models"""
from gpkit.shortcuts import Var, Model
import gpkit
gpkit.disable_units()


class DAPCA4Cost(Model):
    """
    Cost Model from Raymer, Aircraft Design: A Conceptual Approach 3rd Edition,
    pp 713-714

    Note: this is also one instance of the RAND Corp cost models.

    The total cost is Research, Development, Test and Evaluation Cost

    Cost is in 2012 dollars

    Engineering hours considers all airframe, propulsion and avionics
    integration and design.  It does not include engineering of avionics,
    propulsion, tooling or production planning.

    Tooling hours covers all the preparation for production.  It also covers
    ongoing tooling cost through production.

    Manufacturing hours is the direct labor to fabricate the aircraft. It also
    includes cost of labor of subcontractors.

    Quality Control hours covers all inspections.

    Flight test cost covers all cost to prove airworthiness.

    Manufacturing Materials is the raw materials cost and assumes that airframe
    and main materials are made of aluminum.

    The Engine Cost assumes that the cost of the engine is not known.

    Uses 2012 pricing and economics

    Units are not included in variable definitions

    Labor Rates are found from Raymer
    Avionics rate taken from Raymer and based on 15% cost of unit
    Default values based on Honda Jet specs

    Arguments
    ---------
    engine_defined : bool
        If true, must specify T_max, M_max, and T_tur
        If false, must define engine_cost
        Todo: this should really be replaced with an engine cost (sub) model.
    """
    def setup(self, W_e=9200, V=420, Q=5, FTA=1, N_eng=10, T_max=2050,
              M_max=0.63, T_tur=2050, R_avn=5000, engine_cost=500000,
              eng_defined=False):

        # User Definied Variables
        W_e = Var('W_e', W_e, "lb", "Empty Weight")
        V = Var('V', V, "knots", "Maximum Velocity")
        Q = Var('Q', Q, "-", "Number produced in 5 years")
        FTA = Var('FTA', FTA, "-", "Number of Flight Test Aircraft")
        N_eng = Var('N_{eng}', N_eng, "-",
                    "Number of engines or Q*num of engines per aircraft")
        T_max = Var('T_{max}', T_max, "lb", "Maximum Thrust")
        M_max = Var('M_{max}', M_max, "-", "Maximum Mach # of engine")
        T_tur = Var('T_{tur}', T_tur, "R", "Turbine Inlet Temperature")
        R_avn = Var('C_{avn}', R_avn, "$/lb", "Avionics Cost")

        # Constants Hourly Rates
        R_E = Var('R_{E}', 115, "$/hr", "Enginering Hourly Rate")
        R_T = Var('R_{T}', 118, "$/hr", "Tooling Hourly Rate")
        R_M = Var('R_{M}', 108, "$/hr", "Manufacturing Hourly Rate")
        R_Q = Var('R_{Q}', 98, "$/hr", "Quality Check Hourly Rate")

        # Free Variables

        # Hourly costs in hrs. Multiply by hourly work rate to get cost
        H_E = Var('H_{E}', "hrs", "Engineering Time")
        H_T = Var('H_{T}', "hrs", "Tooling Time")
        H_M = Var('H_{M}', "hrs", "Manufacturing Time")
        H_Q = Var('H_{Q}', "hrs", "Quality Check Time")

        # Cost Varibles
        C_D = Var('C_D', "$", "Development Cost")
        C_F = Var('C_F', "$", "Flight Test Cost")
        C_M = Var('C_M', "$", "Manufacturing Materials Cost")
        C_fly = Var('C_{tot}', "$", "RDT&E Flyaway Cost")
        C_avn = Var('C_{avn}', "$", "Cost of avionics")
        C_max = Var('C_{max}', 1000000000000, "$", "Maximum Cost")

        eng = (3112*(0.043*T_max + 243.35*M_max + 0.969*T_tur - 2228)
               if eng_defined else engine_cost)
        C_eng = Var('C_{eng}', eng, '$', "Engine Cost")

        objective = C_fly
        constraints = [C_fly >= (H_E*R_E + H_T*R_T + H_M*R_M + H_Q*R_Q + C_D +
                                 C_F + C_M + C_eng*N_eng + C_avn),
                       H_E >= 4.86 * W_e**0.777 * V**0.894 * Q**0.163,
                       H_T >= 5.99 * W_e**0.777 * V**0.696 * Q**0.263,
                       H_M >= 7.37 * W_e**0.82 * V**0.484 * Q**0.641,
                       H_Q >= 0.133*H_E,
                       C_D >= 91.3 * W_e**0.630 * V**1.3,
                       C_F >= 2498 * W_e**0.325 * V**0.822 * FTA**1.21,
                       C_M >= 22.1 * W_e**0.921 * V**0.621 * Q**0.799,
                       C_avn >= R_avn*W_e,
                       C_max >= C_fly
                      ]
        return objective, constraints


if __name__ == "__main__":
    M = DAPCA4Cost()
    SOL = M.solve()
