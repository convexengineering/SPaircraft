from gpkit import Variable, Model, SignomialsEnabled, units, disable_units
import numpy.testing as npt


class TroposphereSP(Model):
    """
    Density as a function of altitude based on standard atmosphere model

    Required pressures: downward pressure on rho
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
        h   = Variable('h', 'm', 'Altitude')
        rho = Variable('\\rho', 'kg/m^3', 'Density')
        p   = Variable('p', 'Pa', 'Pressure')
        T   = Variable('T', 'K', 'Temperature')

        # Constants
        L   = Variable('L', L, 'K/m', 'Temperature lapse rate')
        MoR = Variable('(M/R)', M/R, 'kg*K/J', 'Air property')
        p_sl = Variable('p_{sl}', 101325, 'Pa', 'Pressure at sea level')
        T_sl = Variable('T_{sl}', 288.15, 'K', 'Temperature at sea level')

        with SignomialsEnabled():
            objective = rho  # minimize density
            constraints = [  # Model only valid up to top of the troposphere
                             h <= 20000*units.m,
                             # Temperature decreases with height at a rate of L
                             T_sl >= T + L*h,
                             # Pressure-altitude relation and ideal gas law
                             (p/p_sl)**(1/th) + L*h/T_sl >= 1,
                             rho == p*MoR/T]
        return objective, constraints

    def test(self):
        sol = self.localsolve()
        h    = sol('h')
        p    = sol('p')
        T    = sol('T')
        rho  = sol('\\rho')
        p_sl = sol('p_{sl}')
        T_sl = sol('T_{sl}')
        L    = sol('L')
        MoR  = sol('(M/R)')

        npt.assert_almost_equal(T_sl, T + L*h, decimal=2)
        npt.assert_almost_equal((p/p_sl)**(1/self.th) + L*h/T_sl, 1, decimal=5)
        npt.assert_almost_equal(rho, p*MoR/T, decimal=5)

if __name__ == "__main__":
    AtmosphereSP().test()
