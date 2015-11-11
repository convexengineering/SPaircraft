from gpkit import Variable, Model


class SolarUAV(Model):
    """LLhale by whoburg, inspired by Drela"""

    def setup(self):
        # Performance
        P = Variable("P", "W", "flight power")
        eta_p = Variable(r"\eta_p", 0.75, "-", "motor+prop efficiency")
        L15D = Variable(r"(C_L^{1.5}/C_D)", 25, "-", "power parameter")
        rho = Variable(r"\rho", ('sweep', [0.238, 0.449]), "kg/m^3",
                       "air density")
        wl = Variable("(m/S)", "kg/m^2", "wing loading")
        g = Variable("g", 9.8, "m/s^2")

        # solar cells
        Psolar = Variable("P_{solar}", 200, "W/m^2", "solar cell output")
        rho_solar = Variable(r"\rho_{solar}", 1.2, "kg/m^2",
                             "solar cell area density")
        eta_daylight = Variable(r"\eta_{daylight}", 0.4, "-",
                                "winter daylight percentage")
        eta_sun_angle = Variable(r"\eta_{sun angle}", 0.5, "-",
                                 "average utilization from low sun angles")

        # battery
        eta_batt = Variable(r"\eta_{batt}", 0.9, "-",
                            "round-trip battery efficiency")
        m_batt = Variable("m_{batt}", "kg", "battery mass")
        h_batt = Variable("h_{batt}", 0.85, "MJ/kg", "battery specific energy")

        # vehicle / airframe / weights
        m_struct = Variable("m_{struct}", "kg", "airframe mass")
        f_struct = Variable("f_{struct}", 0.2, "-",
                            "structural weight fraction, m_struct/m_tot")
        S = Variable("S", "m^2", "wing area")
        m_tot = Variable("m_{tot}", "kg", "airplane mass")
        A = Variable("A", 25, "-", "aspect ratio")
        b = Variable("b", "ft", "span")

        # operations / requirements
        m_pay = Variable("m_{pay}", 3, "kg", "payload mass")
        P_non_flight = Variable(r"P_{non flight}", 20, "W",
                                "non-flight power requirement")
        t_batt = Variable("t_{batt}", 12, "hr", "time on battery")

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
