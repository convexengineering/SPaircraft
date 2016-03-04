from gpkit import Variable, Model, SignomialsEnabled

TS = 110.4    # Sutherland temperature
C1 = 1.458E-6 # Sutherland coefficient

class Sutherland(Model):
    """
    Dynamic viscosity (mu) as a function of temperature

    Assumptions:
    Upward pressure on viscosity

    References:
    http://www-mdp.eng.cam.ac.uk/web/library/enginfo/aerothermal_dvd_only/aero/
           atmos/atmos.html
    http://www.cfd-online.com/Wiki/Sutherland's_law
    """
    def setup(self):
        # Free variables
        mu = Variable("mu", "kg/(m*s)", "Dynamic viscosity")
        T  = Variable("T", "K", "Temperature")

        # Constants
        T_S   = Variable("T_S", TS, "K", "Sutherland Temperature")
        C_1 = Variable("C_1", C1, "kg/(m*s*K^0.5)",
                       "Sutherland coefficient")

        objective = mu
        constraints = [((T+T_S).mono_approximation({T:288.15,
                                                    T_S:T_S.value.magnitude}))
                       * mu == C_1 * T**1.5]

        return objective, constraints
