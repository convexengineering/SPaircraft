"""Aircraft cost models"""
from gpkit import units, Model, Variable


class DAPCA4Cost(Model):
    """
    Cost Model from Raymer, Aircraft Design: A Conceptual Approach 3rd Edition,
    pp 713-714

    Note: this is also one instance of the RAND Corp cost models.

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

        # User Defined Variables
        W_e = Variable("W_e", 9200, "lbf", "Empty Weight")
        V = Variable("V", 420, "knots", "Maximum Velocity")
        Q = Variable("Q", 5, "-", "Number produced in 5 years")
        FTA = Variable("FTA", 1, "-", "Number of Flight Test Aircraft")
        N_eng = Variable("N_{eng}", 10, "-",
                    "Number of engines or Q*num of engines per aircraft")
        T_max = Variable("T_{max}", 2050, "lbf", "Maximum Thrust")
        M_max = Variable("M_{max}", 0.63, "-", "Maximum Mach # of engine")
        T_tur = Variable("T_{tur}", 2050, "R", "Turbine Inlet Temperature")
        R_avn = Variable("C_{avn}", 5000, "USD2012/lbf", "Avionics Cost")

        # Constants Hourly Rates
        R_E = Variable("R_{E}", 115, "USD2012/hr", "Enginering Hourly Rate")
        R_T = Variable("R_{T}", 118, "USD2012/hr", "Tooling Hourly Rate")
        R_M = Variable("R_{M}", 108, "USD2012/hr", "Manufacturing Hourly Rate")
        R_Q = Variable("R_{Q}", 98, "USD2012/hr", "Quality Check Hourly Rate")

        # Free Variables

        # Hourly costs in hrs. Multiply by hourly work rate to get cost
        H_E = Variable("H_{E}", "hrs",
                  "Engineering hours of airframe and component integration")
        # Eng Hours does not include cost of avionics or engine engineering
        H_T = Variable("H_{T}", "hrs", "Tooling hours for production preparation")
        # H_T also covers tooling cost through ongoing production
        H_M = Variable("H_{M}", "hrs",
                  "Manufacturing hours of main and subcontractors")
        H_Q = Variable("H_{Q}", "hrs", "Quality control hours for inspection")

        # Cost Variableibles
        C_D = Variable("C_D", "USD2012", "Development Cost")
        C_F = Variable("C_F", "USD2012", "Flight Test Cost to prove airworthiness")
        C_M = Variable("C_M", "USD2012", "Materials cost of aluminum airframe")
        C_fly = Variable("C_{fly}", "USD2012", "RDT&E Flyaway Cost")
        C_avn = Variable("C_{avn}", "USD2012", "Cost of avionics")
        #Raymer suggests C_avn = 0.15*C_fly

        eng = (3112*(0.043*T_max/units.lbf +
                     243.35*M_max + 0.969*T_tur/units.R - 2228)
               if eng_defined else 500000)
        C_eng = Variable("C_{eng}", eng, "USD2012", "Engine Cost")
        H_E_const = Variable("H_{E_const}", 4.86,
                        units.hr/(units.lbf**0.777 * units.knots**0.894),
                        "H_E power law proportionality constant")
        H_T_const = Variable("H_{T_const}", 5.99,
                        units.hr/(units.lbf**0.777 * units.knots**0.696),
                        "H_T power law proportionality constant")
        H_M_const = Variable("H_{M_const}", 7.37,
                        units.hr/(units.lbf**0.82 *  units.knots**0.484),
                        "H_M power law proportionality constant")
        H_Q_const = Variable("H_{Q_const}", 0.133, "-",
                        "H_Q power law proportionality constant")
        C_D_const = Variable("C_{D_const}", 91.3,
                        units.USD2012/(units.lbf**0.630 * units.knots**1.3),
                        "C_D power law proportionality constant")
        C_F_const = Variable("C_{F_const}", 2498,
                        units.USD2012/(units.lbf**0.325 * units.knots**0.822),
                        "C_F power law proportionality constant")
        C_M_const = Variable("C_{M_const}", 22.1,
                        units.USD2012/(units.lbf**0.921 * units.knots**0.621),
                        "C_M power law proportionality constant")

        objective = C_fly
        constraints = [C_fly >= (H_E*R_E + H_T*R_T + H_M*R_M + H_Q*R_Q + C_D +
                                 C_F + C_M + C_eng*N_eng + C_avn),
                       H_E >= H_E_const * W_e**0.777 * V**0.894 * Q**0.163,
                       H_T >= H_T_const * W_e**0.777 * V**0.696 * Q**0.263,
                       H_M >= H_M_const * W_e**0.82 * V**0.484 * Q**0.641,
                       H_Q >= H_Q_const  *H_E,
                       C_D >= C_D_const * W_e**0.630 * V**1.3,
                       C_F >= C_F_const * W_e**0.325 * V**0.822 * FTA**1.21,
                       C_M >= C_M_const * W_e**0.921 * V**0.621 * Q**0.799,
                       C_avn >= R_avn*W_e,
                      ]
        return objective, constraints

    @classmethod
    def test(cls):
        cls().solve()

if __name__ == "__main__":
    DAPCA4Cost.test()
