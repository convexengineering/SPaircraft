# coding=utf-8
from gpkit import Variable, Model, SignomialsEnabled, units
import numpy as np
import numpy.testing as npt
from gpkit.small_scripts import mag

class LandingGear(Model):
    """
    Landing gear sizing model

    Sources:
    AERO 481 notes
    www.fzt.haw-hamburg.de/pers/Scholz/dglr/hh/text_2002_04_11_Fahrwerk.pdf
    http://faculty.dwc.edu/sadraey/Chapter%209.%20Landing%20Gear%20Design.pdf
    http://www.dept.aoe.vt.edu/~mason/Mason_f/M96SC03.pdf
    http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=M434AE
    http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=M434AE
    https://en.wikipedia.org/wiki/Buckling
    http://arc.aiaa.org/doi/pdf/10.2514/3.7771
    """
    def __init__(self):
        # Variables
        B       = Variable('B', 'm', 'Landing gear base')
        d_oleo  = Variable('d_{oleo}', 'm', 'Diameter of oleo shock absorber')
        dtm     = Variable('d_{t_m}', 'in', 'Diameter of main gear tires')
        dtn     = Variable('d_{t_n}', 'in', 'Diameter of nose gear tires')
        dxm     = Variable('\\Delta x_m', 'm', 'Distance b/w main gear and CG')
        dxn     = Variable('\\Delta x_n', 'm', 'Distance b/w nose gear and CG')
        Eland   = Variable('E_{land}', 'J', 'Max KE to be absorbed in landing')
        Fwm     = Variable('F_{w_m}', '-', 'Weight factor (main)')
        Fwn     = Variable('F_{w_n}', '', 'Weight factor (nose)')
        I_m     = Variable('I_m', 'm^4', 'Area moment of inertia (main strut)')
        I_n     = Variable('I_n', 'm^4', 'Area moment of inertia (nose strut)')
        l_m     = Variable('l_m', 'm', 'Length of main gear')
        l_n     = Variable('l_n', 'm', 'Length of nose gear')
        l_oleo  = Variable('l_{oleo}', 'm', 'Length of oleo shock absorber')
        L_m     = Variable('L_m', 'N', 'Max static load through main gear')
        L_n     = Variable('L_n', 'N', 'Min static load through nose gear')
        L_n_dyn = Variable('L_{n_{dyn}}', 'N', 'Dyn. braking load, nose gear')
        Lwm     = Variable('L_{w_m}', 'N', 'Static load per wheel (main)')
        Lwn     = Variable('L_{w_n}', 'N', 'Static load per wheel (nose)')
        r_m     = Variable('r_m', 'm', 'Radius of main gear struts')
        r_n     = Variable('r_n', 'm', 'Radius of nose gear struts')
        S       = Variable('S', 'm', 'Stroke of the shock absorber')
        S_t     = Variable('S_t', 'm', 'Tire deflection')
        t_m     = Variable('t_m', 'm', 'Thickness of main gear strut wall')
        t_n     = Variable('t_n', 'm', 'Thickness of nose gear strut wall')
        T       = Variable('T', 'm', 'Main landing gear track')
        tan_phi = Variable('\\tan(\\phi)', '-', 'Angle b/w main gear and CG')
        tan_psi = Variable('\\tan(\\psi)', '-', 'Tip over angles')
        Waddm   = Variable('W_{add_m}', 'N', 'Proportional added weight, main')
        Waddn   = Variable('W_{add_n}', 'N', 'Proportional added weight, nose')
        WAWm    = Variable('W_{wa,m}', 'lbf',
                           'Wheel assembly weight for single main gear wheel')
        WAWn    = Variable('W_{wa,n}', 'lbf',
                           'Wheel assembly weight for single nose gear wheel')
        W_lg    = Variable('W_{lg}', 'N', 'Weight of landing gear')
        W_ng    = Variable('W_{ng}', 'N', 'Weight of nose gear')
        W_ns    = Variable('W_{ns}', 'N', 'Weight of nose strut')
        W_nw    = Variable('W_{nw}', 'N', 'Weight of nose wheels (total)')
        W_mg    = Variable('W_{mg}', 'N', 'Weight of main gear')
        W_ms    = Variable('W_{ms}', 'N', 'Weight of main struts')
        W_mw    = Variable('W_{mw}', 'N', 'Weight of main wheels (per strut)')
        wtm     = Variable('w_{t_m}', 'm', 'Width of main tires')
        wtn     = Variable('w_{t_n}', 'm', 'Width of nose tires')
        x_m     = Variable('x_m',  'm', 'x-location of main gear')
        x_n     = Variable('x_n', 'm', 'x-location of nose gear')
        xcg     = Variable('x_{CG}', 'm', 'x-location of CG incl. LG')
        y_m     = Variable('y_m', 'm', 'y-location of main gear (symmetric)')

        # Constants
        d_nac     = Variable('d_{nacelle}', 2, 'm', 'Nacelle diameter')
        E         = Variable('E', 205, 'GPa',
                             'Modulus of elasticity, 4340 steel') # ASM Matweb
        eta_s     = Variable('\\eta_s', 0.8, '-', 'Shock absorber efficiency')
        eta_t     = Variable('\\eta_t', 0.47, '-',
                             'Efficiency of tire in shock absorption')
        g         = Variable('g', 9.81, 'm/s^2', 'Gravitational acceleration')
        h_nac     = Variable('h_{nacelle}', 0.5, 'm', 'Min. nacelle clearance')
        hhold     = Variable('h_{hold}', 1, 'm', 'Hold height')
        K         = Variable('K', 2, '-', 'Column effective length factor')
        lam       = Variable('\\lambda', 2.5, '-',
                             'Ratio of max to static load') # Torenbeek p360
        n_mg      = Variable('n_{mg}', 2, '-', 'Number of main gear struts')
        N_s       = Variable('N_s', 2, '-', 'Factor of safety')
        nwps      = Variable('n_{wps}', 2, '-', 'Number of wheels per strut')
        p_oleo    = Variable('p_{oleo}', 1800, 'lbf/in^2', 'Oleo pressure')
#        p_t       = Variable('p_t', 170, 'lbf/in^2', 'Tyre pressure')
        rho_st    = Variable('\\rho_{st}', 7850, 'kg/m^3',
                             'Density of 4340 Steel') # ASM Matweb
        sig_y_c   = Variable('\\sigma_{y_c}', 470E6, 'Pa',
                             'Compressive yield strength 4340 steel ') #  AZOM
        tan_15    = Variable('\\tan(\\phi_{min})', np.tan(15*np.pi/180), '-',
                             'Lower bound on phi')
        tan_63    = Variable('\\tan(\\psi_{max})', np.tan(63*np.pi/180), '-',
                             'Upper bound on psi')
        tangam    = Variable('\\tan(\\gamma)', np.tan(5*np.pi/180), '-',
                             'Tangent of dihedral angle')
        tan_th_TO = Variable('\\tan(\\theta_{TO})', np.tan(15*np.pi/180), '-',
                             'Takeoff pitch angle')
        w         = Variable('w', 10, 'ft/s', 'Ultimate velocity of descent')
        W_0       = Variable('W_0', 82000*9.81, 'N', # 737 Max airport doc
                             'Weight of aircraft excluding landing gear')
        xcg0      = Variable('x_{CG_0}', 18, 'm', 'x-location of CG excl. LG')
        x_upswp   = Variable('x_{up}', 28, 'm', 'Fuselage upsweep point')
        y_eng     = Variable('y_{eng}', 4.83, 'm', 'Spanwise loc. of engines')
        z_CG_0    = Variable('z_{CG}', 2, 'm',
                             'CG height relative to bottom of fuselage')
        zwing     = Variable('z_{wing}', 0.5, 'm',
                             'Height of wing relative to base of fuselage')

        with SignomialsEnabled():

            objective = W_lg

            constraints = [
                           # Track and Base geometry definitions
                           l_n + zwing + y_m*tangam >= l_m, # [SP]
                           T == 2*y_m,
                           x_n + B <= x_m,
                           x_n >= 5*units.m, # nose gear after nose

                           # Geometric constraints relating gear placement with
                           # fore/aft CG locations
                           dxn + x_n >= xcg, # [SP] Needs to be tight
                           dxm + xcg >= x_m, # [SP] Needs to be tight
                           # TODO forward and aft CG

                           # CG location affected by landing gear position
                           xcg*(W_0 + W_lg) <= W_0*xcg0 + W_ng*x_n + W_mg*x_m,

                           # Maximum static loads through main and nose gears
                           L_n == W_0*dxm/B,
                           L_m == W_0*dxn/B,

                           # Dynamic braking load through nose gear
                           # (assumes deceleration of 10 ft/s^2)
                           L_n_dyn >= 0.31*((z_CG_0+l_m)/B)*W_0,

                           # For steering don't want too much or too little
                           # load on nose gear
                           L_n/W_0 >= 0.05,
                           L_n/W_0 <= 0.20,

                           # Longitudinal tip over (static)
                           x_m >= tan_phi*(z_CG_0+l_m) + xcg,
                           tan_phi >= tan_15,

                           # Lateral tip over in turn (dynamic)
                           # www.dept.aoe.vt.edu/~mason/Mason_f/M96SC03.pdf
                           # stricter constraint uses forward CG
                           # cos(arctan(y/x))) = x/sqrt(x^2 + y^2)
                           1 >= (z_CG_0 + l_m)**2 * (y_m**2 + B**2) /
                                (dxn * y_m * tan_psi)**2,
                           tan_psi <= tan_63,

                           # Tail strike: Longitudinal ground clearance in
                           # takeoff, landing (Raymer says 10-15 degrees)
                           # TODO?: 2 cases:(i) upsweep angle > rotation angle,
                           # (ii) upsweep angle < rotation ang
                           x_upswp - x_m <= l_m/tan_th_TO, # [SP]

                           # Engine ground clearance
                           d_nac + h_nac <= l_m + (y_eng-y_m)*tangam, # [SP]

                           # Size/Volume for retraction
                           y_m >= l_m,

                           # Constrains/is constrained by engine location
                           y_m <= y_eng,

                           # Brake sizing for stopping aircraft

                           # Hard landing
                           # http://www.boeing.com/commercial/aeromagazine/...
                           # articles/qtr_3_07/AERO_Q307_article3.pdf
                           # sink rate of 10 feet per second at the maximum
                           # design landing weight
                           # Landing condition from Torenbeek p360
                           Eland == W_0/(2*g)*w**2, # Torenbeek (10-26)
                           # S_t == 0.5*lam*Lwm/(p*(dtm*bt)**0.5), # (10-30)
                           S == (1/eta_s)*(Eland/(L_m*lam)),# - eta_t*S_t),
                           # [SP] Torenbeek (10-28)

                           l_oleo == 2.5*S, # Raymer 244
                           d_oleo == 1.3*(4*lam*L_m/n_mg/(np.pi*p_oleo))**0.5,
                           l_m >= l_oleo + dtm/2,

                           # Wheel weights
                           Fwm == Lwm*dtm/(1000*Lwm.units*dtm.units),
                           WAWm == 1.2*Fwm**0.609*units.lbf,# Currey p145
                           Fwn == Lwn*dtn/(1000*units.lbf*units.inches),
                           WAWn == 1.2*Fwn**0.609*units.lbf,# Currey p145
                           Lwm == L_m/(n_mg*nwps),
                           Lwn == L_n/nwps,

                           # Main wheel diameter/width (Raymer p233)
                           dtm == 1.63*(Lwm/(4.44*units.N))**0.315*units.inch,
                           wtm == 0.1043*(Lwm/(4.44*units.N))**0.48*units.inch,
                           dtn == 0.8*dtm,
                           wtn == 0.8*wtm,

                           # TODO: Beam sizing and max bending (limits track,
                           # downward pressure on y_m)

                           # Weight is a function of height and load through
                           # each strut as well as tyre size (and obviously
                           # number of struts)
                           # Main gear strut weight (for a single strut) is a
                           # function of length and load passing through strut
                           W_ms >= 2*np.pi*r_m*t_m*l_m * rho_st * g,
                           # Compressive yield in hard landing condition
                           N_s * lam * L_m/n_mg <= sig_y_c * (2*np.pi*r_m*t_m),
                           W_mw == nwps*WAWm,

                           # Nose gear strut weight is a function of length and
                           # load passing through strut
                           W_ns >= 2*np.pi*r_n*t_n*l_n * rho_st * g,
                           # find cross sectional area based on compressive yield
                           N_s * (L_n + L_n_dyn) <= sig_y_c*(2*np.pi*r_n*t_n),
                           W_nw >= nwps*WAWn,

                           # Buckling constraint on main gear
                           L_m <= np.pi**2*E*I_m/(K*l_m)**2,
                           I_m <= np.pi*r_m**3*t_m,

                           # Buckling constraint on nose gear
                           # source: https://en.wikipedia.org/wiki/Buckling
                           L_n <= np.pi**2*E*I_n/(K*l_n)**2,
                           I_n <= np.pi*r_n**3*t_n,

                           # Machining constraint # p89 Mason
                           # www.dept.aoe.vt.edu/~mason/Mason_f/M96SC08.pdf
                           2*r_m/t_m <= 40,
                           2*r_m/t_n <= 40,

                           # Retraction constraint on strut diameter
                           2*wtm + 2*r_m <= hhold,
                           2*wtn + 2*r_n <= 0.8*units.m, #TODO improve this

                           # Weight accounting
                           Waddm >= 1.5*W_mw, # Currey p264
                           Waddn >= 1.5*W_nw,
                           W_mg >= n_mg*(W_ms + W_mw + Waddm),
                           W_ng >= W_ns + W_nw + Waddn,
                           W_lg >= W_mg + W_ng,
                          ]

        Model.__init__(self, objective, constraints)

    @classmethod
    def test(cls):
        sol = cls().localsolve()
        npt.assert_almost_equal(mag(sol('x_n')) + mag(sol('B')), mag(sol('x_m')), decimal=5)
        npt.assert_almost_equal(mag(sol('\Delta x_n')) + mag(sol('x_n')), mag(sol('x_{CG}')),
                                decimal=1)
        npt.assert_almost_equal(mag(sol('\Delta x_m')) + mag(sol('x_{CG}')) , mag(sol('x_m')),
                                decimal=4)
        npt.assert_almost_equal(mag(sol('l_n')) + mag(sol('z_{wing}')) + mag(sol('y_m')) *
                                mag(sol('\\tan(\\gamma)')), mag(sol('l_m')), decimal=4)
        npt.assert_almost_equal(mag(sol('d_{nacelle}')) + mag(sol('h_{nacelle}')),
                                mag(sol('l_m')) + (mag(sol('y_{eng}'))-mag(sol('y_m'))) *
                                mag(sol('\\tan(\\gamma)')), decimal=4)

if __name__ == "__main__":
    LandingGear.test()
