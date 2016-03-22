from gpkit import Variable, Model, units, SignomialsEnabled
from gpkit.tools import te_exp_minus1
from gpkit.constraints.tight import TightConstraintSet as TCS
import numpy as np
import numpy.testing as npt

class Troposphere(Model):
    """
    Density and dynamic viscosity as a function of altitude based on
    standard atmosphere model for the Troposphere

    Required pressures: downward pressure on rho
    Assumptions: only valid to the top of the Troposphere (~11km)

    Sources:
    https://en.wikipedia.org/wiki/Density_of_air#Altitude
    http://www.digitaldutch.com/atmoscalc/index.htm
    http://www.engineeringtoolbox.com/standard-atmosphere-d_604.html
    http://www-mdp.eng.cam.ac.uk/web/library/enginfo/aerothermal_dvd_only/aero/
           atmos/atmos.html
    http://www.cfd-online.com/Wiki/Sutherland's_law

    Arguments
    ---------
    g : float [m/s^2]
        Acceleration due to gravity
    L : float  [K/m]
        Temperature lapse rate
    R : float  [J/(kg*K)]
        Specific gas constant (air)
    """
    def __init__(self, g=9.81, L=0.0065, R=287, **kwargs):
        th = g/(R*L)  # dimensionless
        self.g, self.L, self.R, self.th = g, L, R, th

        # Free variables
        h   = Variable('h', 'm', 'Altitude')
        mu  = Variable("mu", "kg/(m*s)", "Dynamic viscosity")
        p   = Variable('p', 'Pa', 'Pressure')
        rho = Variable('\\rho', 'kg/m^3', 'Density')
        T   = Variable('T', 'K', 'Temperature')

        # Constants
        C_1  = Variable("C_1", 1.458E-6, "kg/(m*s*K^0.5)",
                        "Sutherland coefficient")
        L    = Variable('L', L, 'K/m', 'Temperature lapse rate')
        p_0  = Variable('p_0', 101325, 'Pa', 'Pressure at sea level')
        R    = Variable('R', R, 'J/(kg*K)', 'Specific gas constant (air)')
        T_0  = Variable('T_0', 288.15, 'K', 'Temperature at sea level')
        T_S  = Variable("T_S", 110.4, "K", "Sutherland Temperature")

        with SignomialsEnabled():
            objective = rho  # minimize density
            constraints = [  # Model only valid up to top of the troposphere
                             h <= 11000*units.m,

                             # Temperature decreases with height at a rate of L
                             TCS([T_0 <= T + L*h]),

                             # Pressure-altitude relation
                             (p/p_0)**(1/th) == T/T_0,
 
                             # Ideal gas law
                             rho == p/(R*T),

                             # Sutherland viscosity model
                             C_1*T**1.5/mu == (T+T_S).mono_approximation(
                                              {T: 288.15, T_S: T_S.value})
                          ]
        Model.__init__(self, objective, constraints, **kwargs)

    def test(self):
        sol = self.localsolve()
        h, p, T, rho, p_0, T_0, L, R = map(sol, ("h p T \\rho p_0 "
                                                   "T_0 L R").split())
        npt.assert_almost_equal(T_0/(T + L*h), 1, decimal=5)

class Tropopause(Model):
    """
    Density and dynamic viscosity as a function of altitude based on
    standard atmosphere model for the Tropopause

    Assumptions: Only valid in the Tropopause (11 km - 20 km)

    Sources:
    http://www.digitaldutch.com/atmoscalc/index.htm
    http://www.engineeringtoolbox.com/standard-atmosphere-d_604.html
    http://www-mdp.eng.cam.ac.uk/web/library/enginfo/aerothermal_dvd_only/aero/
           atmos/atmos.html
    http://www.cfd-online.com/Wiki/Sutherland's_law

    Arguments
    ---------
    g : float [m/s^2]
        Acceleration due to gravity
    R : float  [J/(kg*K)]
        Specific gas constant (air)
    T : float [K]
        Temperature in the tropopause
    """
    def __init__(self, g=9.81, T=216.65, R=287, **kwargs):
        k = g/(R*T)

        # Free variables
        h   = Variable('h', 'm', 'Altitude')
        mu  = Variable("mu", "kg/(m*s)", "Dynamic viscosity")
        p   = Variable('p', 'Pa', 'Pressure')
        rho = Variable('\\rho', 'kg/m^3', 'Density')

        # Constants
        C_1  = Variable('C_1', 1.458E-6, "kg/(m*s*K^0.5)",
                        "Sutherland coefficient")
        g    = Variable('g', g, 'm/s^2', 'Graivational acceleration')
        p11  = Variable('p_{11}', 22630, 'Pa', 'Pressure at 11 km')
        R    = Variable('R', R, 'J/(kg*K)', 'Specific gas constant (air)')
        T    = Variable('T', T, 'K', 'Temperature')
        T_S  = Variable("T_S", 110.4, "K", "Sutherland Temperature")

        objective = 1/rho  # maximize density
        constraints = [  # Model only valid up to top of the troposphere
                         h >= 11*units.km,
                         h <= 20*units.km,

                         # Pressure-altitude relation, using taylor series exp
                         TCS([np.exp(k*11000)*p11/p >=
                              1 + te_exp_minus1(g/(R*T)*h, 15)]),

                         # Ideal gas law
                         rho == p/(R*T),

                         # Sutherland viscosity model
                         C_1*T**1.5/mu == (T+T_S).mono_approximation(
                                          {T: 288.15, T_S: T_S.value})
                      ]

        substitutions = {'T': 216.65} # Temperature is constant 

        Model.__init__(self, objective, constraints, substitutions, **kwargs)

if __name__ == "__main__":
    TS = Troposphere()
    TS.test()
    TP = Tropopause()
    TP.solve()
