# coding=utf-8
"Implements Landing Gear model"
import numpy as np
from gpkit import Variable, Model, SignomialsEnabled, units
from gpkit.constraints.costed import CostedConstraintSet
from gpkit.constraints.tight import TightConstraintSet as TCS

class LandingGear(CostedConstraintSet):
    """
    Landing gear sizing constraint set
    
    Sources:
    AERO 481 notes
    www.fzt.haw-hamburg.de/pers/Scholz/dglr/hh/text_2002_04_11_Fahrwerk.pdf
    http://faculty.dwc.edu/sadraey/Chapter%209.%20Landing%20Gear%20Design.pdf
    http://www.dept.aoe.vt.edu/~mason/Mason_f/M96SC03.pdf
    http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=M434AE
    https://en.wikipedia.org/wiki/Buckling
    """
    def __init__(self):
        B       = Variable('B', 'm', 'Landing gear base')
        E       = Variable('E', 'GPa', 'Modulus of elasticity, 4340 steel')
        Eland   = Variable('E_{land}', 'J', 'Max KE to be absorbed in landing')
        Fwm     = Variable('F_{w_m}', '-', 'Weight factor (main)')
        Fwn     = Variable('F_{w_n}', '-', 'Weight factor (nose)')
        I_m     = Variable('I_m', 'm^4', 'Area moment of inertia (main strut)')
        I_n     = Variable('I_n', 'm^4', 'Area moment of inertia (nose strut)')
        K       = Variable('K', '-', 'Column effective length factor')
        L_m     = Variable('L_m', 'N', 'Max static load through main gear')
        L_n     = Variable('L_n', 'N', 'Min static load through nose gear')
        L_n_dyn = Variable('L_{n_{dyn}}', 'N', 'Dyn. braking load, nose gear')
        Lwm     = Variable('L_{w_m}', 'N', 'Static load per wheel (main)')
        Lwn     = Variable('L_{w_n}', 'N', 'Static load per wheel (nose)')
        N_s     = Variable('N_s', '-', 'Factor of safety')
        S_sa    = Variable('S_sa', 'm', 'Stroke of the shock absorber')
        S_t     = Variable('S_t', 'm', 'Tire deflection')
        T       = Variable('T', 'm', 'Main landing gear track')
        W       = Variable('W', 'N', 'Total aircraft weight')
        WAWm    = Variable('W_{wa,m}', 'lbf',
                           'Wheel assembly weight for single main gear wheel')
        WAWn    = Variable('W_{wa,n}', 'lbf',
                           'Wheel assembly weight for single nose gear wheel')
        W_0     = Variable('W_{0_{lg}}', 'N',
                           'Weight of aircraft excluding landing gear')
        W_lg    = Variable('W_{lg}', 'N', 'Weight of landing gear')
        W_mg    = Variable('W_{mg}', 'N', 'Weight of main gear')
        W_ms    = Variable('W_{ms}', 'N', 'Weight of main struts')
        W_mw    = Variable('W_{mw}', 'N', 'Weight of main wheels (per strut)')
        W_ng    = Variable('W_{ng}', 'N', 'Weight of nose gear')
        W_ns    = Variable('W_{ns}', 'N', 'Weight of nose strut')
        W_nw    = Variable('W_{nw}', 'N', 'Weight of nose wheels (total)')
        d_fan   = Variable('d_{fan}', 'm', 'Fan diameter')
        d_nac   = Variable('d_{nacelle}', 'm', 'Nacelle diameter')
        d_oleo  = Variable('d_{oleo}', 'm', 'Diameter of oleo shock absorber')
        dtm     = Variable('d_{t_m}', 'in', 'Diameter of main gear tires')
        dtn     = Variable('d_{t_n}', 'in', 'Diameter of nose gear tires')
        dxm     = Variable('\\Delta x_m', 'm', 'Distance b/w main gear and CG')
        dxn     = Variable('\\Delta x_n', 'm', 'Distance b/w nose gear and CG')
        eta_s   = Variable('\\eta_s', '-', 'Shock absorber efficiency')
        eta_t   = Variable('\\eta_t', '-',
                           'Efficiency of tire in shock absorption')
        faddm   = Variable('f_{add,m}', '-', 'Proportional added weight, main')
        faddn   = Variable('f_{add,n}', '-', 'Proportional added weight, nose')
        g       = Variable('g', 9.81, 'm/s^2', 'Gravitational acceleration')
        h_nac   = Variable('h_{nacelle}', 'm', 'Min. nacelle clearance')
        hhold   = Variable('h_{hold}', 'm', 'Hold height')
        l_m     = Variable('l_m', 'm', 'Length of main gear')
        l_n     = Variable('l_n', 'm', 'Length of nose gear')
        l_oleo  = Variable('l_{oleo}', 'm', 'Length of oleo shock absorber')
        lam     = Variable('\\lambda_{LG}', '-',
                           'Ratio of max to static load') # Torenbeek p360
        n_mg    = Variable('n_{mg}', '-', 'Number of main gear struts')
        nwps    = Variable('n_{wps}', '-', 'Number of wheels per strut')
        p_oleo  = Variable('p_{oleo}', 'lbf/in^2', 'Oleo pressure')
       #p_t     = Variable('p_t', 170, 'lbf/in^2', 'Tyre pressure')
        r_m     = Variable('r_m', 'm', 'Radius of main gear struts')
        r_n     = Variable('r_n', 'm', 'Radius of nose gear struts')
        rho_st  = Variable('\\rho_{st}', 'kg/m^3', 'Density of 4340 Steel')
        sig_y_c = Variable('\\sigma_{y_c}', 'Pa',
                           'Compressive yield strength 4340 steel') #  AZOM
        t_m     = Variable('t_m', 'm', 'Thickness of main gear strut wall')
        t_n     = Variable('t_n', 'm', 'Thickness of nose gear strut wall')
        t_nac   = Variable('t_{nacelle}', 'm', 'Nacelle thickness')
        tan_15  = Variable('\\tan(\\phi_{min})', '-', 'Lower bound on phi')
        tan_63  = Variable('\\tan(\\psi_{max})', '-', 'Upper bound on psi')
        tan_gam = Variable('\\tan(\\gamma)', '-', 'Tangent, dihedral angle')
        tan_phi = Variable('\\tan(\\phi)', '-', 'Angle b/w main gear and CG')
        tan_psi = Variable('\\tan(\\psi)', '-', 'Tip over angles')
        tan_th0 = Variable('\\tan(\\theta_{TO})', '-', 'Takeoff pitch angle')
        w_ult   = Variable('w_{ult}', 'ft/s', 'Ultimate velocity of descent')
        wtm     = Variable('w_{t_m}', 'm', 'Width of main tires')
        wtn     = Variable('w_{t_n}', 'm', 'Width of nose tires')
        x_m     = Variable('x_m', 'm', 'x-location of main gear')
        x_n     = Variable('x_n', 'm', 'x-location of nose gear')
        x_upswp = Variable('x_{up}', 'm', 'Fuselage upsweep point')
        xcg     = Variable('x_{CG}', 'm', 'x-location of CG incl. LG')
        xcglg   = Variable('x_{CG_{lg}}', 'm', 'Landing gear CG')
        xcg0    = Variable('x_{CG_0}', 'm', 'x-location of CG excl. LG')
        y_eng   = Variable('y_{eng}', 'm', 'Spanwise loc. of engines')
        y_m     = Variable('y_m', 'm', 'y-location of main gear (symmetric)')
        z_CG_0  = Variable('z_{CG}', 'm',
                           'CG height relative to bottom of fuselage')
        zwing   = Variable('z_{wing}', 'm',
                           'Height of wing relative to base of fuselage')


        with SignomialsEnabled():

            objective = W_lg

            constraints = [
                           # Track and Base geometry definitions
                           TCS([l_n+zwing+y_m*tan_gam>=l_m], reltol=1E-3), #[SP]
                           T == 2*y_m,
                           TCS([x_n + B <= x_m]),
                           x_n >= 5*units.m, # nose gear after nose

                           # Geometric constraints relating gear placement with
                           # fore/aft CG locations
                           TCS([dxn + x_n >= xcg], reltol=1E-3), # [SP]
                           TCS([dxm + xcg >= x_m], reltol=1E-4), # [SP]
                           # TODO forward and aft CG

                           # Maximum static loads through main and nose gears
                           L_n == W*dxm/B,
                           L_m == W*dxn/B,

                           # Dynamic braking load through nose gear
                           # (assumes deceleration of 10 ft/s^2)
                           L_n_dyn >= 0.31*((z_CG_0+l_m)/B)*W,

                           # For steering don't want too much or too little
                           # load on nose gear
                           L_n/W >= 0.05,
                           L_n/W <= 0.20,

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
                           x_upswp - x_m <= l_m/tan_th0, # [SP]

                           # Engine ground clearance
                           d_nac >= d_fan + 2*t_nac,
                           d_nac + h_nac <= l_m + (y_eng-y_m)*tan_gam, # [SP]

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
                           Eland >= W/(2*g)*w_ult**2, # Torenbeek (10-26)
                           # S_t == 0.5*lam*Lwm/(p*(dtm*bt)**0.5), # (10-30)
                           S_sa == (1/eta_s)*(Eland/(L_m*lam)),# - eta_t*S_t),
                           # [SP] Torenbeek (10-28)

                           l_oleo == 2.5*S_sa, # Raymer 244
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
                           I_m == np.pi*r_m**3*t_m,

                           # Buckling constraint on nose gear
                           # source: https://en.wikipedia.org/wiki/Buckling
                           L_n <= np.pi**2*E*I_n/(K*l_n)**2,
                           I_n == np.pi*r_n**3*t_n,

                           # Machining constraint # p89 Mason
                           # www.dept.aoe.vt.edu/~mason/Mason_f/M96SC08.pdf
                           2*r_m/t_m <= 40,
                           2*r_m/t_n <= 40,

                           # Retraction constraint on strut diameter
                           2*wtm + 2*r_m <= hhold,
                           2*wtn + 2*r_n <= 0.8*units.m, #TODO improve this

                           # Weight accounting
                           W_mg >= n_mg*(W_ms + W_mw*(1 + faddm)),# Currey p264
                           W_ng >= W_ns + W_nw*(1 + faddn),
                           W_lg >= W_mg + W_ng,
                          ]

            standaloneCG = [
                           # CG location affected by landing gear position
                           TCS([xcg*W <= W_0*xcg0 + W_ng*x_n + W_mg*x_m]),
                           W >= W_0 + W_lg,
                           ]

            coupledCG = [
                         TCS([W_lg*xcglg >= W_ng*x_n + W_mg*x_m], reltol=1E-2,
                             raiseerror=False),
                         x_m >= xcglg,
                        ]

            self.standaloneCG = standaloneCG
            self.coupledCG = coupledCG

        CostedConstraintSet.__init__(self, objective, constraints)


    def default737subs(self):

        substitutions = {
                         'E': 205,
                         'K': 2,
                         'N_s': 2,
                         'W_{0_{lg}}': 82000*9.81,
                         '\\eta_s': 0.8,
                         '\\eta_t': 0.47,
                         '\\lambda_{LG}': 2.5,
                         '\\rho_{st}': 7850,
                         '\\tan(\\gamma)': np.tan(5*np.pi/180),
                         '\\tan(\\phi_{min})': np.tan(15*np.pi/180),
                         '\\tan(\\psi_{max})': np.tan(63*np.pi/180),
                         '\\tan(\\theta_{TO})': np.tan(15*np.pi/180),
                         '\\sigma_{y_c}': 470E6,
                         'd_{fan}': 1.75,
                         'f_{add,m}': 1.5,
                         'f_{add,n}': 1.5,
                         'h_{hold}': 1,
                         'h_{nacelle}': 0.5,
                         'n_{mg}': 2,
                         'n_{wps}': 2,
                         'p_{oleo}': 1800,
                         't_{nacelle}': 0.15,
                         'w_{ult}': 10,
                         'x_{CG_0}': 18,
                         'x_{up}': 28,
                         'y_{eng}': 4.83,
                         'z_{CG}': 2,
                         'z_{wing}': 0.5,
                        } 

        return substitutions

    @classmethod
    def standalone737(cls):
        """Returns a standalone landing gear model"""
        ccs = cls()

        constraints = ccs + ccs.standaloneCG
 
        substitutions = ccs.default737subs()

        m =  Model(ccs.cost, constraints, substitutions)
        return m

    @classmethod
    def coupled737(cls):
        """Returns a landing gear model for use in a coupled aircraft model"""
        ccs = cls()

        constraints = ccs + ccs.coupledCG

        dsubs = ccs.default737subs()
        linkedsubs = ['h_{hold}', 'x_{up}']
        substitutions = {key: value for key, value in dsubs.items()
                                    if key not in linkedsubs}

        m =  Model(ccs.cost, constraints, substitutions, name='LandingGear')
        return m
        

    @classmethod
    def test(cls):
        """Tests the standalone landing gear model"""
        m = cls.standalone737()
        sol = m.localsolve()

if __name__ == "__main__":
    LandingGear.test()
