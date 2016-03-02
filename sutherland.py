from gpkit import Variable, Model, SignomialsEnabled

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
        T_S   = Variable("T_S", 110.4, "K", "Sutherland Temperature")
        C_1 = Variable("C_1", 1.458E-6, "kg/(m*s*K^0.5)",
                       "Sutherland coefficient")

        objective = 1/mu
        constraints = [(1 + (T_S/T)) * mu <= C_1 * T**0.5]

        return objective, constraints

class SutherlandSP(Model):
    """
    Dynamic viscosity (mu) as a function of temperature

    Assumptions:
    Downward pressure on viscosity

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
        T_S   = Variable("T_S", 110.4, "K", "Sutherland Temperature")
        C_1 = Variable("C_1", 1.458E-6, "kg/(m*s*K^0.5)",
                       "Sutherland coefficient")

        objective = 1/mu
        with SignomialsEnabled():
            constraints = [(1 + (T_S/T)) * mu >= C_1 * T**0.5]

        return objective, constraints
