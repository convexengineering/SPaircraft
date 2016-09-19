from gpkit import Variable, VectorVariable, Model, units, SignomialsEnabled, SignomialEquality
from gpkit.constraints.set import ConstraintSet

class Atmosphere(Model):
    def __init__(self, **kwargs):
        g = Variable('g', 'm/s^2', 'Gravitational acceleration')
        p_sl = Variable("p_{sl}", "Pa", "Pressure at sea level")
        T_sl = Variable("T_{sl}", "K", "Temperature at sea level")
        L_atm = Variable("L_{atm}", "K/m", "Temperature lapse rate")
        T_atm = Variable("T_{atm}", "K", "air temperature")
        M_atm = Variable("M_{atm}", "kg/mol",
                         "Molar mass of dry air")
        R_atm = Variable("R_{atm}", "J/mol/K",
                         "air specific heating value")
        p_atm = Variable("P_{atm}", "Pa", "air pressure")
        TH = (g*M_atm/R_atm/L_atm).value

        rho = Variable('\\rho', 'kg/m^3', 'Density of air')

        h = Variable("h", "ft", "Altitude")

        """
        Dynamic viscosity (mu) as a function of temperature
        References:
        http://www-mdp.eng.cam.ac.uk/web/library/enginfo/aerothermal_dvd_only/aero/
            atmos/atmos.html
        http://www.cfd-online.com/Wiki/Sutherland's_law
        """
        mu  = Variable('\\mu', 'kg/(m*s)', 'Dynamic viscosity')

        T_s = Variable('T_s', 110.4, "K", "Sutherland Temperature")
        C_1 = Variable('C_1', 1.458E-6, "kg/(m*s*K^0.5)",
                       'Sutherland coefficient')

        t_plus_ts_approx = (T_atm + T_s).mono_approximation({T_atm: 288.15,
                                                         T_s: T_s.value})

        with SignomialsEnabled():
            constraints = [
                # Pressure-altitude relation
                (p_atm/p_sl)**(1/5.257) == T_atm/T_sl,

                # Ideal gas law
                rho == p_atm/(R_atm/M_atm*T_atm),

                #temperature equation
                SignomialEquality(T_sl, T_atm + L_atm*h),

                #constraint on mu
                t_plus_ts_approx * mu == C_1 * T_atm**1.5,
                ]


        Model.__init__(self, T_atm, constraints, **kwargs)
        
if __name__ == "__main__":
    M = Atmosphere()
    sol = M.localsolve("mosek",iteration_limit=500)
    print sol.table()
