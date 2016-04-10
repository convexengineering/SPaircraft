"""Gas-powered high-altitude-long-endurange UAV model"""
from numpy import pi
from gpkit import Variable, Model, units
from gpkit.tools import te_exp_minus1

class GasPoweredHALE(Model):
    def setup(self):
        constraints = []

        # Steady level flight relations
        CD = Variable('C_D', '-', 'Drag coefficient')
        CL = Variable('C_L', '-', 'Lift coefficient')
        P_shaft = Variable('P_{shaft}', 'W', 'Shaft power')
        S = Variable('S', 'm^2', 'Wing reference area')
        V = Variable('V', 'm/s', 'Cruise velocity')
        W = Variable('W', 'lbf', 'Aircraft weight')

        eta_prop = Variable(r'\eta_{prop}', 0.7, '-', 'Propulsive efficiency')
        rho = Variable(r'\rho', 'kg/m^3')

        constraints.extend([P_shaft == V*W*CD/CL/eta_prop,   # eta*P = D*V
                            W == 0.5*rho*V**2*CL*S])

        # Aerodynamics model
        Cd0 = Variable('C_{d0}', 0.01, '-', "non-wing drag coefficient")
        CLmax = Variable('C_{L-max}', 1.5, '-', 'maximum lift coefficient')
        e = Variable('e', 0.9, '-', "spanwise efficiency")
        A = Variable('A', 20, '-', "aspect ratio")
        b = Variable('b', 'ft', 'span')
        mu = Variable(r'\mu', 1.5e-5, 'N*s/m^2', "dynamic viscosity")
        Re = Variable("Re", '-', "Reynolds number")
        Cf = Variable("C_f", "-", "wing skin friction coefficient")
        Kwing = Variable("K_{wing}", 1.3, "-", "wing form factor")
        cl_16 = Variable("cl_{16}", 0.0001, "-", "profile stall coefficient")
        constraints.extend([CD >= Cd0 + 2*Cf*Kwing + CL**2/(pi*e*A) + cl_16*CL**16,
                            b**2 == S*A,
                            CL <= CLmax,
                            Re == rho*V/mu*(S/A)**0.5,
                            Cf >= 0.074/Re**0.2])


        # Engine Weight Model
        W_eng = Variable('W_{eng}', 'lbf', 'Engine weight')
        W_engmin = Variable('W_{min}', 1.3, 'lbf', 'min engine weight')
        W_engmax = Variable('W_{max}', 275, 'lbf', 'max engine weight')
        eta_t = Variable('\\eta_t', 0.75, '-', 'percent throttle')
        eng_cnst = Variable('eng_{cnst}', 0.0011, '-',
                            'engine constant based off of power weight model')
        constraints.extend([W_eng >= W_engmin,
                            W_eng <= W_engmax,
                            W_eng >= P_shaft*eng_cnst/eta_t * units('lbf/watt')])

        # Weight model
        W_airframe = Variable('W_{airframe}', 'lbf', 'Airframe weight')
        W_pay = Variable(r'W_{pay}', 4, 'lbf', 'Payload weight')
        W_fuel = Variable('W_{fuel}', 'lbf', 'Fuel Weight')
        W_zfw = Variable('W_{zfw}', 'lbf', 'Zero fuel weight')
        wl = Variable('wl', 'lbf/ft^2', 'wing loading')

        f_airframe = Variable('f_{airframe}', 0.3, '-',
                              'Airframe weight fraction')
        g = Variable('g', 9.81, 'm/s^2', 'Gravitational acceleration')

        constraints.extend([W_airframe >= W*f_airframe,
                            W_zfw >= W_airframe + W_eng + W_pay,
                            wl == W/S,
                            W >= W_pay + W_eng + W_airframe + W_fuel])

        # Breguet Range
        z_bre = Variable("z_bre", "-", "breguet coefficient")
        h_fuel = Variable("h_{fuel}", 42e6, "J/kg", "heat of combustion")
        eta_0 = Variable("\\eta_0", 0.2, "-", "overall efficiency")
        t = Variable('t', 7, 'days', 'time on station')

        constraints.extend([z_bre >= g*t*V*CD/(h_fuel*eta_0*CL),
                            W_fuel/W_zfw >= te_exp_minus1(z_bre, 3)])

        # Atmosphere model
        h = Variable("h", "ft", "Altitude")
        p_sl = Variable("p_{sl}", 101325, "Pa", "Pressure at sea level")
        T_sl = Variable("T_{sl}", 288.15, "K", "Temperature at sea level")
        L_atm = Variable("L_{atm}", 0.0065, "K/m", "Temperature lapse rate")
        T_atm = Variable("T_{atm}", "K", "air temperature")
        M_atm = Variable("M_{atm}", 0.0289644, "kg/mol",
                         "Molar mass of dry air")
        R_atm = Variable("R_{atm}", 8.31447, "J/mol/K")
        TH = (g*M_atm/R_atm/L_atm).value.magnitude  # dimensionless
        constraints.extend([h <= 20000*units.m,  # Model valid to top of troposphere
                            T_sl >= T_atm + L_atm*h,     # Temp decreases w/ altitude
                            rho <= p_sl*T_atm**(TH-1)*M_atm/R_atm/(T_sl**TH)])
            # http://en.wikipedia.org/wiki/Density_of_air#Altitude

        # station keeping requirement
        footprint = Variable("d_{footprint}", 100, 'km',
                             "station keeping footprint diameter")
        lu = Variable(r"\theta_{look-up}", 5, '-', "look up angle")
        R_earth = Variable("R_{earth}", 6371, "km", "Radius of earth")
        tan_lu = lu*pi/180. + (lu*pi/180.)**3/3.  # Taylor series expansion
        # approximate earth curvature penalty as distance^2/(2*Re)
        constraints.extend([
            h >= tan_lu*0.5*footprint + footprint**2/8./R_earth])

        # wind speed model
        V_wind = Variable('V_{wind}', 'm/s', 'wind speed')
        wd_cnst = Variable('wd_{cnst}', [0.002, 0.0015], 'm/s/ft',
                           'wind speed constant predited by model')
        wd_ln = Variable('wd_{ln}', [13.009, 8.845], 'm/s',
                         'linear wind speed variable')
        h_min = Variable('h_{min}', 11800, 'ft', 'minimum height')
        h_max = Variable('h_{max}', 20800, 'ft', 'maximum height')
        constraints.extend([V_wind >= wd_cnst*h + wd_ln, # model predicting worse case scenario at 45 deg latitude
                            V >= V_wind,
                            h >= h_min,
                            h <= h_max])

        objective = W

        return objective, constraints

if __name__ == "__main__":
    M = GasPoweredHALE()
    M.solve()
