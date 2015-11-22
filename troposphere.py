from gpkit import Variable, Model, units
import numpy.testing as npt


class Troposphere(Model):
    """
    Density as a function of altitude based on standard atmosphere model

    Required pressures: upward pressure on rho and h
    Assumptions: only valid to the top of the troposphere (~20km)

    source: https://en.wikipedia.org/wiki/Density_of_air#Altitude

    Arguments
    ---------
    g : float [m/s^2]
        Acceleration due to gravity
    L : float  [K/m]
        Temperature lapse rate
    M : float [kg/mol]
        Molar mass of dry air
    R : float  [J/(mol*K)]
        Universal gas constant
    """
    def setup(self, g=9.81, L=0.0065, M=0.0289644, R=8.31447):
        th = (g*M)/(R*L)  # dimensionless
        self.g, self.L, self.M, self.R, self.th = g, L, M, R, th

        # Free variables
        h = Variable("h", "m", "Altitude")
        rho = Variable("\\rho", "kg/m^3", "Density")
        T = Variable("T", "K", "Temperature")

        # Constants
        L = Variable("L", L, "K/m", "Temperature lapse rate")
        MoR = Variable("(M/R)", M/R, "kg*K/J", "Air property")
        p_sl = Variable("p_{sl}", 101325, "Pa", "Pressure at sea level")
        T_sl = Variable("T_{sl}", 288.15, "K", "Temperature at sea level")

        objective = rho**-1 * h**-1  # maximize height * density
        constraints = [  # Model only valid up to top of the troposphere
                         h <= 20000*units.m,
                         # Temperature decreases with altitude at a rate of L
                         T_sl >= T + L*h,
                         # Pressure-altitude relation and ideal gas law
                         rho <= p_sl*T**(th-1)*MoR/(T_sl**th)]
        return objective, constraints

    def test(self):
        sol = self.solve()
        h, T, rho, p_sl, T_sl, L, MoR = map(sol, ["h", "T", "\\rho", "p_{sl}",
                                                   "T_{sl}", "L", "(M/R)"])
        npt.assert_almost_equal(T_sl, T + L*h, decimal=2)
        calculated_rho = p_sl*T**(self.th-1)*MoR/(T_sl**self.th)
        npt.assert_almost_equal(calculated_rho, rho, decimal=5)

if __name__ == "__main__":
    Troposphere().test()
