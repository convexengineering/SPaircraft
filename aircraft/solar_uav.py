from gpkit import Variable, Model
Var = Variable


class SolarUAV(Model):
    """LLhale by whoburg, inspired by Drela"""

    def setup(self):
        # Performance
        P = Var("P", "W", "flight power")
        eta_p = Var(r"\eta_p", 0.75, "-", "motor+prop efficiency")
        L15D = Var(r"(C_L^{1.5}/C_D)", 25, "-", "power parameter")
        rho = Var(r"\rho", ('sweep', [0.238, 0.449]), "kg/m^3", "air density")
        wl = Var("(m/S)", "kg/m^2", "wing loading")
        g = Var("g", 9.8, "m/s^2")

        # solar cells
        Psolar = Var("P_{solar}", 200, "W/m^2", "solar cell output")
        rho_solar = Var(r"\rho_{solar}", 1.2, "kg/m^2",
                        "solar cell area density")
        eta_daylight = Var(r"\eta_{daylight}", 0.4, "-",
                           "winter daylight percentage")
        eta_sun_angle = Var(r"\eta_{sun angle}", 0.5, "-",
                            "average utilization from low sun angles")

        # battery
        eta_batt = Var(r"\eta_{batt}", 0.9, "-",
                       "round-trip battery efficiency")
        m_batt = Var("m_{batt}", "kg", "battery mass")
        h_batt = Var("h_{batt}", 0.85, "MJ/kg", "battery specific energy")

        # vehicle / airframe / weights
        m_struct = Var("m_{struct}", "kg", "airframe mass")
        f_struct = Var("f_{struct}", 0.2, "-",
                       "structural weight fraction, m_struct/m_tot")
        S = Var("S", "m^2", "wing area")
        m_tot = Var("m_{tot}", "kg", "airplane mass")
        A = Var("A", 25, "-", "aspect ratio")
        b = Var("b", "ft", "span")

        # operations / requirements
        m_pay = Var("m_{pay}", 3, "kg", "payload mass")
        P_non_flight = Var(r"P_{non flight}", 20, "W",
                           "non-flight power requirement")
        t_batt = Var("t_{batt}", 12, "hr", "time on battery")

        return Model(S,
                     [P/S >= 1/eta_p*(1/L15D)*(2/rho)**0.5*(wl*g)**1.5,
                      P + P_non_flight <=
                      Psolar*S*eta_daylight*eta_sun_angle*eta_batt,
                      m_tot >= S*rho_solar + m_pay + m_batt + m_struct,
                      m_batt*h_batt >= (P + P_non_flight)*t_batt,
                      m_struct >= f_struct*m_tot,
                      wl*S >= m_tot,
                      b**2 == S*A])


if __name__ == "__main__":
    m = SolarUAV()
    m.solve()
