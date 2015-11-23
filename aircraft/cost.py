"""Aircraft cost models"""
from gpkit.shortcuts import Var, Model
import gpkit
gpkit.disable_units()

class DAPCA4Cost(Model):
    """
    Cost Model from Raymer, Aircraft Design: A Conceptual Approach 3rd Edition,
    pp 713-714

    Note: this is also one instance of the RAND Corp cost models.

    Cost is in 2012 dollars

    Units disabled because '$' is not defined as a unit 

    Labor Rates are found from Raymer

    Default values based on Honda Jet specs

    Arguments
    ---------
    engine_defined : bool
        If true, must specify T_max, M_max, and T_tur
        If false, must define engine_cost
        Todo: this should really be replaced with an engine cost (sub) model.
    """
    def setup(self, eng_defined=False):

        # User Definied Variables
        W_e = Var('W_e', 9200, "lb", "Empty Weight")
        V = Var('V', 420, "knots", "Maximum Velocity")
        Q = Var('Q', 5, "-", "Number produced in 5 years")
        FTA = Var('FTA', 1, "-", "Number of Flight Test Aircraft")
        N_eng = Var('N_{eng}', 10, "-",
                    "Number of engines or Q*num of engines per aircraft")
        T_max = Var('T_{max}', 2050, "lb", "Maximum Thrust")
        M_max = Var('M_{max}', 0.63, "-", "Maximum Mach # of engine")
        T_tur = Var('T_{tur}', 2050, "R", "Turbine Inlet Temperature")
        R_avn = Var('C_{avn}', 5000, "$/lb", "Avionics Cost")

        # Constants Hourly Rates
        R_E = Var('R_{E}', 115, "$/hr", "Enginering Hourly Rate")
        R_T = Var('R_{T}', 118, "$/hr", "Tooling Hourly Rate")
        R_M = Var('R_{M}', 108, "$/hr", "Manufacturing Hourly Rate")
        R_Q = Var('R_{Q}', 98, "$/hr", "Quality Check Hourly Rate")

        # Free Variables

        # Hourly costs in hrs. Multiply by hourly work rate to get cost
        H_E = Var('H_{E}', "hrs",
                  "Engineering hours of airframe and component integration")
        #Eng Hours does not include cost of avionics or engine engineering
        H_T = Var('H_{T}', "hrs", "Tooling hours for production preparation")
        #H_T also covers tooling cost through ongoing production
        H_M = Var('H_{M}', "hrs",
                  "Manufacturing hours of main and subcontractors")
        H_Q = Var('H_{Q}', "hrs", "Quality control hours for inspection")

        # Cost Varibles
        C_D = Var('C_D', "$", "Development Cost")
        C_F = Var('C_F', "$", "Flight Test Cost to prove airworthiness")
        C_M = Var('C_M', "$", "Materials cost of aluminum airframe")
        C_fly = Var('C_{tot}', "$", "RDT&E Flyaway Cost")
        C_avn = Var('C_{avn}', "$", "Cost of avionics")
        #Raymer suggests C_avn = 0.15*C_fly

        eng = (3112*(0.043*T_max + 243.35*M_max + 0.969*T_tur - 2228)
               if eng_defined else 500000)
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
                      ]
        return objective, constraints

    def test(self):
        sol = self.solve()

if __name__ == "__main__":
    DAPCA4Cost().test()
