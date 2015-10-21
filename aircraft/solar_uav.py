"""LLhale by whoburg, inspired by Drela"""
from gpkit import Variable, Model


def solar_uav():
    """Returns the Model"""
    # Performance
    P = Variable("P", units="W", label="flight power")
    eta_p = Variable(r"\eta_p", value=0.75, label="motor+prop efficiency")
    L15D = Variable(r"(C_L^{1.5}/C_D)", value=25, label="power parameter")
    rho = Variable(r"\rho", units="kg/m^3", label="air density",
                   value=('sweep', [0.238, 0.449]))
    wl = Variable("(m/S)", units="kg/m^2", label="wing loading")
    g = Variable("g", units="m/s^2", value=9.8)

    # solar cells
    Psolar = Variable("P_{solar}", units="W/m^2", value=200,
                      label="solar cell output")
    rho_solar = Variable(r"\rho_{solar}", units="kg/m^2", value=1.2,
                         label="solar cell area density")
    eta_daylight = Variable(r"\eta_{daylight}", value=0.4,
                            label="winter daylight percentage")
    eta_sun_angle = Variable(r"\eta_{sun angle}", value=0.5,
                             label="average utilization from low sun angles")

    # battery
    eta_batt = Variable(r"\eta_{batt}", value=0.9,
                        label="round-trip battery efficiency")
    m_batt = Variable("m_{batt}", units="kg", label="battery mass")
    h_batt = Variable("h_{batt}", units="MJ/kg", value=0.85,
                      label="battery specific energy")

    # vehicle / airframe / weights
    m_struct = Variable("m_{struct}", units="kg", label="airframe mass")
    f_struct = Variable("f_{struct}", value=0.2,
                        label="structural weight fraction, m_struct/m_tot")
    S = Variable("S", units="m^2", label="wing area")
    m_tot = Variable("m_{tot}", units="kg", label="airplane mass")
    A = Variable("A", value=25, label="aspect ratio")
    b = Variable("b", units="ft", label="span")

    # operations / requirements
    m_pay = Variable("m_{pay}", units="kg", value=3, label="payload mass")
    P_non_flight = Variable(r"P_{non flight}", units="W", value=20,
                            label="non-flight power requirement")
    t_batt = Variable("t_{batt}", value=12, units="hr", label="time on batt")

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
    M = solar_uav()
    M.solve()
