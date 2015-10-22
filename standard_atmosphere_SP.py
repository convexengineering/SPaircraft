from gpkit import Variable, Model, SignomialsEnabled, units, disable_units
import numpy.testing as npt

class Air_Density_SP(Model):
    """
    Density as a function of altitude based on standard atmosphere model

    Required pressures: downward pressure on rho
    Assumptions: only valid to the top of the troposphere (~20km)

    source: https://en.wikipedia.org/wiki/Density_of_air#Altitude
    """
    def setup(self):
        disable_units() # TODO Remove this
        g  = 9.81        # Acceleration due to gravity [m/s^2]
        L  = 0.0065      # Temperature lapse rate [K/m]
        M  = 0.0289644   # Molar mass of dry air [kg/mol]
        R  = 8.31447     # Universal gas constant [J/(mol*K)]
        th = (g*M)/(R*L) # [-]

        # Free variables
        h   = Variable('h', 'm')             # Altitude
        p   = Variable('p', 'Pa')            # Pressure
        rho = Variable('rho', 'kg/m^3')      # Air density
        T   = Variable('T', 'K')             # Temperature

        # Constants
        L   = Variable('L', L, 'K/m')        # Temperature lapse rate
        MoR = Variable('MoR', M/R, 'kg*K/J') # M/R (see above)
        p0  = Variable('p0', 101325, 'Pa')   # Sea level standard pressure
        T0  = Variable('T0', 288.15, 'K')    # Sea level standard temperature

        with SignomialsEnabled():
            objective = rho
            constraints = [
                            # Model only valid up to top of the troposphere
                            h <= 20000,#*units.m,

                            # Temperature decreases with altitude at a rate of L
                            T0 >= T + L*h,

                            # Combination of pressure altitude relationship and
                            # ideal gas law
                            (p/p0)**(1/th) + L*h/T0 >= 1,
                            rho == p*MoR/T
                          ]
        return objective, constraints

if __name__ == "__main__":
    m = Air_Density_SP()
    sol = m.localsolve()

    th = 5.2575767

    h   = sol['variables']['h']
    p   = sol['variables']['p']
    T   = sol['variables']['T']
    rho = sol['variables']['rho']

    p0  = sol['constants']['p0']
    T0  = sol['constants']['T0']
    L   = sol['constants']['L']
    MoR = sol['constants']['MoR']

    npt.assert_almost_equal(T0, T + L*h, decimal=2)
    npt.assert_almost_equal((p/p0)**(1/th) + L*h/T0, 1, decimal=5)
    npt.assert_almost_equal(rho, p*MoR/T, decimal=5)
