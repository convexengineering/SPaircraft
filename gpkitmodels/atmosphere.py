"Models for atmospheric quantities"
import numpy as np
from gpkit import Variable, Model, units, SignomialsEnabled, LinkConstraint
from gpkit.tools import te_exp_minus1
from gpkit.constraints.tight import TightConstraintSet as TCS
# pylint: disable=bad-whitespace

GAS_CONSTANT = 287  # [J/(kg*K)]
GRAVITATIONAL_ACCEL = 9.81 # [m/s^2]

g   = Variable('g', GRAVITATIONAL_ACCEL, 'm/s^2', 'Gravitational acceleration')
h   = Variable('h', 'm', 'Altitude')
mu  = Variable('\\mu', 'kg/(m*s)', 'Dynamic viscosity')
p   = Variable('p', 'Pa', 'Pressure')
R   = Variable('R', GAS_CONSTANT, 'J/(kg*K)', 'Specific gas constant (air)')
rho = Variable('\\rho', 'kg/m^3', 'Density')
T   = Variable('T', 'K', 'Temperature')


class Troposphere(Model):
    """
    Density and dynamic viscosity as a function of altitude based on
    standard atmosphere model for the Troposphere

    Assumptions: only valid to the top of the Troposphere (~11km)

    References:
    Anderson, Introduction to Flight
    https://en.wikipedia.org/wiki/Density_of_air#Altitude
    http://www.digitaldutch.com/atmoscalc/index.htm

    Arguments
    ---------
    min_rho: bool
        If true, model expects downward external pressure on rho
        If false, model expects upward external pressure on rho
    """
    def __init__(self, min_rho=True, **kwargs):
        self.min_rho = min_rho

        Lval = 0.0065   # [K/m]
        th = GRAVITATIONAL_ACCEL/(GAS_CONSTANT*Lval) # [-]

        L    = Variable('L', Lval, 'K/m', 'Temperature lapse rate')
        p_0  = Variable('p_0', 101325, 'Pa', 'Pressure at sea level')
        T_0  = Variable('T_0', 288.15, 'K', 'Temperature at sea level')

        if min_rho:
            objective = rho  # minimize density
        else:
            objective = 1/rho # maximize density

        # Temperature lapse rate constraint
        if min_rho:
            with SignomialsEnabled():
                constraints = TCS([T_0 <= T + L*h])
        else:
            constraints = TCS([T_0 >= T + L*h])

        constraints += [h <= 11000*units.m,
                        h >= 1E-6*units.m,

                        # Pressure-altitude relation
                        (p/p_0)**(1/th) == T/T_0,

                        # Ideal gas law
                        rho == p/(R*T),
                       ]

        su = Sutherland()
        lc = LinkConstraint([constraints, su])

        Model.__init__(self, objective, lc, **kwargs)

    def test(self):
        if self.min_rho:
            self.localsolve()
        else:
            self.solve()

class Tropopause(Model):
    """
    Density and dynamic viscosity as a function of altitude based on
    standard atmosphere model for the Tropopause

    Assumptions: Only valid in the Tropopause (11 km - 20 km)

    References:
    Anderson, Introduction to Flight
    http://www.digitaldutch.com/atmoscalc/index.htm
    """
    def __init__(self, **kwargs):
        T_tp = 216.65
        k = GRAVITATIONAL_ACCEL/(GAS_CONSTANT*T_tp)

        p11  = Variable('p_{11}', 22630, 'Pa', 'Pressure at 11 km')

        objective = 1/rho  # maximize density
        constraints = [h >= 11*units.km,
                       h <= 20*units.km,

                       # Temperature is constant in the tropopause
                       T == T_tp,

                       # Pressure-altitude relation, using taylor series exp
                       TCS([np.exp(k*11000)*p11/p >=
                            1 + te_exp_minus1(g/(R*T)*h, 15)], reltol=1E-4),

                       # Ideal gas law
                       rho == p/(R*T),
                      ]

        su = Sutherland()
        lc = LinkConstraint([constraints, su])

        Model.__init__(self, objective, lc, **kwargs)


class Sutherland(Model):
    """
    Dynamic viscosity (mu) as a function of temperature

    References:
    http://www-mdp.eng.cam.ac.uk/web/library/enginfo/aerothermal_dvd_only/aero/
        atmos/atmos.html
    http://www.cfd-online.com/Wiki/Sutherland's_law
    """
    def __init__(self, **kwargs):

        T_s = Variable('T_s', 110.4, "K", "Sutherland Temperature")
        C_1 = Variable('C_1', 1.458E-6, "kg/(m*s*K^0.5)",
                       'Sutherland coefficient')

        t_plus_ts_approx = (T + T_s).mono_approximation({T: 288.15,
                                                         T_s: T_s.value})

        objective = mu
        constraints = [t_plus_ts_approx * mu == C_1 * T**1.5]

        Model.__init__(self, objective, constraints, **kwargs)


if __name__ == "__main__":
    TS = Troposphere()
    TS.test()
    TP = Tropopause()
    TP.solve()
    SU = Sutherland()
    # TODO: move this to a test() method
    SU.substitutions.update({SU["T"]: 288})
    SU.solve()
