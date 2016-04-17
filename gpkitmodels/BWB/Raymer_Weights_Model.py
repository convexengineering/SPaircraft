# Below is a Weights breakdown for the Blended wing body.  This method is GP Compatible!

from gpkit import VectorVariable, Variable, Model, units
import numpy as np
from gpkit.tools import te_exp_minus1
import gpkit
gpkit.disable_units()
import matplotlib.pyplot as plt

def secant(x, nterm):
    """Taylor expansion of secant(x).

    Arguments
    ---------
    x : gpkit.monomial
      Variable or expression argument
    nterm : int
        Number of terms in resulting Taylor expansion

    Returns
    -------
    gpkit.Posynomial
        Taylor expansion of secant(x), carried to nterm terms
    """
    if nterm < 1:
        raise ValueError("Unexpected number of terms, nterm=%s" % nterm)

    # The first 12 Euler Numbers
    E2n = np.asarray([ 1, 
                       5,
                       61,
                       1385,
                       50521,
                       2702765,
                       199360981,
                       19391512145,
                       2404879675441, 
                       370371188237525,
                       69348874393137901,
                       15514534163557086905] )
    if nterm > 12:
        n_extend = np.asarray(range(13, nterm+1))
        E2n_add = (8 * np.sqrt(n_extend/np.pi) 
                      * (4*n_extend/(np.pi * np.exp(1)))**(2*n_extend))
        E2n = np.append(E2n, E2n_add)

    res = 1
    factorial_denom = 1
    for i in range(1, nterm + 1):
        factorial_denom *= ((2*i)*(2*i-1))
        res +=  E2n[i-1]/ factorial_denom * x**(2*i)
    return res

def tangent(x, nterm):
    """Taylor expansion of tangent(x).

    Arguments
    ---------
    x : gpkit.monomial
      Variable or expression argument
    nterm : int
        Number of terms in resulting Taylor expansion

    Returns
    -------
    gpkit.Posynomial
        Taylor expansion of tangent(x), carried to nterm terms
    """
    if nterm < 1:
        raise ValueError("Unexpected number of terms, nterm=%s" % nterm)

    if nterm > 15:
        raise ValueError("Tangent expansion not implemented above 15 terms")

    # The first 15 Bernoulli Numbers
    B2n = np.asarray([  1/6.,
                        -1/30.,
                        1/42.,
                        -1/30.,
                        5/66.,
                        -691/2730.,
                        7/6.,
                        -3617/510.,
                        43867/798.,
                        -174611/330.,
                        854513/138.,
                        -236364091/2730.,
                        8553103/6.,
                        -23749461029/870.,
                        8615841276005/14322. ])

    res = 0
    factorial_denom = 1
    for i in range(1, nterm + 1):
        factorial_denom *= ((2*i)*(2*i-1))
        res +=  (-1)**(i-1) * 2**(2*i) * (2**(2*i) - 1) * B2n[i-1] / factorial_denom * x**(2*i-1)
    return res

def Weights_pie(sol):
    Wts = sol.program.result["variables"]

    plt.figure()
    plt.pie([Wts["W_{operating-empty}"],Wts["W_{c}"],Wts["W_{fuel}"]], 
            labels=['OEW','Payload','Fuel'],
            startangle=90,
            explode=(0.05, 0, 0),
            autopct='%1.1f%%', shadow=True, radius=5.0)
    plt.title("MTOW: %.0f lbs" %(Wts["W_{maximum-gross}"]), y = 1.03)
    plt.axis('equal')
    plt.savefig('MTOW_breakdown.eps')

    plt.figure()
    W_body =  Wts["W_{centerbody}"] + Wts["W_{afterbody}"]
    W_gear =  Wts["W_{main-gear}"] + Wts["W_{nose-gear}"]
    W_prop =  Wts["W_{engine-contents}"] + Wts["W_{nacelle}"] + Wts["W_{engine-controls}"] + Wts["N_{en}"]*Wts["W_{en}"]
    W_aux  =  ( Wts["W_{starter-penumatic}"] + 
                Wts["W_{fuel-system}"] + 
                Wts["W_{flight-controls}"] + 
                Wts["W_{APU_{installed}}"] + 
                Wts["W_{instruments}"] + 
                Wts["W_{hydraulics}"] + 
                Wts["W_{electrical}"] + 
                Wts["W_{avionics}"] + 
                Wts["W_{furnishings}"] + 
                Wts["W_{air-conditioning}"] + 
                Wts["W_{handling-gear}"] + 
                Wts["W_{anti-ice}"] + 
                Wts["W_{military-cargo-system}"] )
    plt.pie([W_body,
              Wts["W_{wing}"],
              Wts["W_{vertical-tail}"],
              W_gear,
              W_prop,
              W_aux],
            labels = ['Body','Wing','Vertical Tail','Landing Gear','Propulsion','Auxillary'],
            startangle=90,
            explode=(0,0,0,0,0,0),
            autopct='%1.1f%%', shadow=True, radius=5.0)
    plt.title("OEW: %.0f lbs" %(Wts["W_{operating-empty}"]), y = 1.03)
    plt.axis('equal')
    plt.savefig('OEW_breakdown.eps')



# Fixed
Bref                = Variable("B_{ref}",      6.25,       "ft",       "Reference Wing Span for Wing Weight")
Kieg                = Variable("K_{ieg}",      0.575,      "-",        "Instrumentation Constant")
Kmp                 = Variable("K_{mp}",       1.0,        "-",        "1.126 for kneeling gear; 1.0 Otherwise")
Kng                 = Variable("K_{ng}",       1.107,      "-",        "1.107 for pylon mounted nacelle; 1.0 Otherwise")
Knp                 = Variable("K_{np}",       1.0,        "-",        "1.15 for kneeling gear; 1.0 Otherwise")
Kp                  = Variable("K_{p}",        1.0,        "-",        "1.4 for Engine With Propeller; 1.0 Otherwise")
Ktr                 = Variable("K_{tr}",       1.18,       "-",        "1.18 for Jet with Thrust Reverser; 1.0 Otherwise")
Kw                  = Variable("K_{w}",        1.7e-3,     "-",        "Proportionality Factor for Wing")


# Possibly changing, but should remained fixed
npax = 450
Wpax = 210
Nc                  = Variable("N_{c}",        15.,         "-",        "Number of Crew")
Nen                 = Variable("N_{en}",       2.,          "-",        "Number of Engines")
Nf                  = Variable("N_{f}",        7.,          "-",        "Number of Functions Performed by the Controls (Typically 4-7")
Ngen                = Variable("N_{gen}",      2.,          "-",        "Number of Generators (Typically = N_{en}")
Nl                  = Variable("N_{l}",        4.5,         "-",        "Ultimate Landing Loading Factor; = N_{gear}x1.5")
Nmss                = Variable("N_{mss}",      2.,          "-",        "Number of Main Gear Shock Struts")
Nmw                 = Variable("N_{mw}",       8.,          "-",        "Number of Main Wheels")
Nnw                 = Variable("N_{nw}",       2.,          "-",        "Number of Nose Wheels")
Np                  = Variable("N_{p}",        npax,        "-",        "Number of Personnel Onboard (Crew and Passengers")
Nt                  = Variable("N_{t}",        2.,          "-",        "Number of Fuel Tanks")
Nz                  = Variable("N_{z}",        2.5,         "-",        "Ultimate Load Factor; = 1.5 x Limit Load Factor")
Rkva                = Variable("R_{kva}",      50.,         "kV*A",     "System Electrical Rating (Typically 40-60) ")


# Performance
R                   = Variable("R",           8000.,        "nmi",      "Range")
V_stall             = Variable("V_{stall}",   120.,         "m/s",      "Stall Speed")
SFC                 = Variable("SFC",         0.5,          "1/s",      "Specific Fuel Consumption")
V                   = Variable("V_{cruise}",                "m/s",      "Cruise Speed")
zbre                = Variable("Z_{bre}",                   "-",        "Breguet Range Term")

# Aero
L = Variable("L", "lbf", "Lift")
CL = Variable("C_L", "-", "Lift Coefficient")
D = Variable("D", "lbf", "Drag")
CD = Variable("C_D", "-", "Drag Coefficient")
CDp = Variable("C_{D_p}", "-", "Profile Drag Coefficient")
CDp0 = Variable("C_{D_p0}", "-", "Profile Drag Coefficient, No Mach Correction")
rho = Variable("\\rho",7.38e-4, "slg/ft^3", "Air Density")
e = Variable("e", 0.9, "-", "Oswald Efficiency")
AR = Variable("AR", 10.0, "-", "Aspect Ratio")
# k_w = Variable("k_w", 1.7e-3, "-", "Proportionality factor for wing")
# n_ult = Variable("n_{ult}", 2.5, "-", "Maximum Load Factor")
mu = Variable("\\mu", 2.995e-7, "lb*s/ft^2", "Air Viscosity")
a = Variable("a", 915, "ft/s", "Speed of Sound")
Re   = Variable("Re", "-", "Reynolds Number")
M    = Variable("M", 0.86, "-", "Mach Number")
chord  = Variable("chord", 15.0, "ft", "Wing Chord")
S = Variable("S", "ft^2", "Planform Area")
pi = Variable("\\pi", np.pi, "-", "pi") 

# Geometry Dependent
Av                  = Variable("A_{v}",                   "-",        "Vertical Tail Aspect Ratio")
Bw                  = Variable("B_{w}",                   "ft",       "Wing Span")
Kz                  = Variable("K_{z}",                   "ft",       "Aircraft Yawing Radius of Gyration (\\approx L_{t})")
La                  = Variable("L_{a}",                   "ft",       "Electrical Routing Distance, Generators to Avionics to Cockpit")
Lec                 = Variable("L_{ec}",                  "ft",       "Length From Engine Front to Cockpit--Total if Multi-Engine")
Lf                  = Variable("L_{f}",                   "ft",       "Aircraft Length Measured on Centerline")
Lm                  = Variable("L_{m}",                   "in",       "Extended Length of Main Landing Gear")
Ln                  = Variable("L_{n}",                   "in",       "Extended Nose Gear Length")
Lt                  = Variable("L_{t}",                   "ft",       "Tail Length; Wing Quarter-MAC to Tail Quarter-MAC")
NLt                 = Variable("N_{Lt}",                  "ft",       "Nacelle Length")
Nw                  = Variable("N_{w}",                   "ft",       "Nacelle Width")
Saft                = Variable("S_{aft}",                 "ft^2",     "Total Area of Afterbody")
Scab                = Variable("S_{cab}",                 "ft^2",     "Total Area of Cabin")
Sn                  = Variable("S_{n}",                   "ft^2",     "Nacelle Wetted Area")
Svt                 = Variable("S_{vt}",                  "ft^2",     "Vertical Tail Area")
Sw                  = Variable("S_{w}",                   "ft^2",     "Trapezoidal Wing Area")
tau_w               = Variable("\\tau_w",                 "-",        "Thickness to Chord Ratio of Main Wing at the Root")
tau_v               = Variable("\\tau_v",                 "-",        "Thickness to Chord Ratio of Vertical Tail at the Root")
Vpr                 = Variable("V_{pr}",                  "ft^3",     "Volume of Pressurized Section")
Vt                  = Variable("V_{t}",                   "ft^3",     "Total Fuel Volume")
Lambda_vt           = Variable("\\Lambda_{vt}",           "-",        "Vertical Tail Sweep at 0.25 MAC")
lambda_aft          = Variable("\\lambda_{aft}",          "-",        "Afterbody Taper Ratio")

# Preset Weights
W_APU_uninstalled   = Variable("W_{APU_{uninstalled}}",   2000.,     "lb",       "Uninstalled APU Weight")
Wc                  = Variable("W_{c}",                   npax*Wpax, "lb",       "Maximum Cargo Weight")
Wen                 = Variable("W_{en}",                  12000.,    "lb",       "Weight of Engine")
Wuav                = Variable("W_{uav}",                 1400.,    "-",        "Uninstalled Avionics Weight (Typically 800-1400")

# Calculated Weights
W_operating_empty                 = Variable("W_{operating-empty}",        "lbf",       "Operating Empty Weight")
W_maximum_gross                   = Variable("W_{maximum-gross}",          "lbf",       "Flight Design Gross Weight")
W_fuel                            = Variable("W_{fuel}",                   "lbf",       "Fuel Weight")
W_engine_contents                 = Variable("W_{engine-contents}",        "lbf",       "Weight of Engine and Contents (per nacelle)")
W_landing_gross                   = Variable("W_{landing_gross}",          "lbf",       "Landing Design Gross Weight")
W_centerbody                      = Variable("W_{centerbody}",             "lbf",       "Weight of Centerbody")
W_afterbody                       = Variable("W_{afterbody}",              "lbf",       "Weight of Afterbody")
W_wing                            = Variable("W_{wing}",                   "lbf",       "Weight of Main Wing")
W_vertical_tail                   = Variable("W_{vertical-tail}",          "lbf",       "Weight of Vertical Tail")
W_main_landing_gear               = Variable("W_{main-gear}",              "lbf",       "Weight of Main Landing Gear")
W_nose_landing_gear               = Variable("W_{nose-gear}",              "lbf",       "Weight of Nose Landing Gear")
W_nacelle_group                   = Variable("W_{nacelle}",                "lbf",       "Weight of Nacelle Group")
W_engine_controls                 = Variable("W_{engine-controls}",        "lbf",       "Weight of Engine Controls")
W_starter_penumatic               = Variable("W_{starter-penumatic}",      "lbf",       "Weight of Pneumatic Starter")
W_fuel_system                     = Variable("W_{fuel-system}",            "lbf",       "Weight of Fuel System")
W_flight_controls                 = Variable("W_{flight-controls}",        "lbf",       "Weight of Flight Controls")
W_APU_installed                   = Variable("W_{APU_{installed}}",        "lbf",       "Weight of Installed APU")
W_instruments                     = Variable("W_{instruments}",            "lbf",       "Weight of Instruments")
W_hydraulics                      = Variable("W_{hydraulics}",             "lbf",       "Weight of Hydraulic System")
W_electrical                      = Variable("W_{electrical}",             "lbf",       "Weight of Electrical System")
W_avionics                        = Variable("W_{avionics}",               "lbf",       "Weight of Avionics")
W_furnishings                     = Variable("W_{furnishings}",            "lbf",       "Weight of Furnishings (No seats or handling gear")
W_air_conditioning                = Variable("W_{air-conditioning}",       "lbf",       "Weight of Air Conditioning System")
W_anti_ice                        = Variable("W_{anti-ice}",               "lbf",       "Weight of Anti-Icing System")
W_handling_gear                   = Variable("W_{handling-gear}",          "lbf",       "Weight of Cargo Handling Gear")
W_military_cargo_handling_system  = Variable("W_{military-cargo-system}",  "lbf",       "Weight of Military Cargo Handling System")


constraints = [ # Bradley
                W_centerbody >= 5.698865 * 0.316422 * (W_maximum_gross ** 0.166552) * Scab ** 1.061158,
                W_afterbody >= (1.0 + 0.05*Nen) * 0.53 * Saft * (W_maximum_gross**0.2) * (lambda_aft + 0.5),
                # Raymer page 459
                W_wing >= W_maximum_gross * Kw * (Bw)**(0.75) * (1+(Bref/Bw)**(0.5)) * Nz**(0.55) * (Bw/(tau_w*chord))**(0.3) * (W_maximum_gross/(Sw*.2))**(-0.3), #Torenbeek for GP compatibility
                W_vertical_tail >= 12000,
                W_main_landing_gear >= 0.0106 * Kmp * W_landing_gross**(0.888) * Nl**(0.25) * Lm**(0.4) * Nmw**(0.321) * Nmss**(-0.5) * V_stall**(0.1),
                W_nose_landing_gear >= 0.032 * Knp * W_landing_gross**(0.646) * Nl**(0.2) * Ln**(0.5) * Nnw**(0.45),
                W_engine_contents >= 2.331 * Wen**(0.901) * Kp * Ktr,
                W_nacelle_group >= 0.6724 * Kng * NLt**(0.10) * Nw**(0.294) * Nz**(0.119) * W_engine_contents**(0.611) * Nen**(0.984) * Sn**(0.224),
                W_engine_controls >= 5.0 * Nen + 0.80 * Lec,
                W_starter_penumatic >= 49.19 * (Nen * Wen / 1000.)**(0.541),
                W_fuel_system >= 2.405 * (Vt*7.48052)**(0.606) * Nt**(0.5),
                W_flight_controls >= 0.005 * W_maximum_gross,
                W_APU_installed >= 2.2 * W_APU_uninstalled,
                W_instruments >= Kieg * W_operating_empty**(0.555) * R**(0.25), # Torenbeek for GP compatibility
                W_hydraulics >= 0.2673 * Nf * (Lf + Bw),
                W_electrical >= 7.291 * Rkva**(0.782) * La**(0.346) * Ngen**(0.10),
                W_avionics >= 1.73 * Wuav**(0.983),
                W_furnishings >= 0.0577 * Nc**(0.1) * Wc**(0.393) * Scab**(0.75),
                W_air_conditioning >= 62.36 * Np**(0.25) * (Vpr/1000.)**(0.604) * Wuav**(0.01),
                W_handling_gear >= 0.002 * W_maximum_gross,
                W_anti_ice >= 0.0003 * W_maximum_gross,
                W_military_cargo_handling_system  >= 2.4 * Scab,
              ]

# OEW
constraints += [
                W_operating_empty  >= (W_centerbody
                                    + W_afterbody 
                                    + W_wing 
                                    + W_vertical_tail 
                                    + W_main_landing_gear
                                    + W_nose_landing_gear
                                    + W_engine_contents
                                    + W_nacelle_group
                                    + W_engine_controls 
                                    + W_starter_penumatic
                                    + W_fuel_system
                                    + W_flight_controls
                                    + W_APU_installed
                                    + W_instruments
                                    + W_hydraulics 
                                    + W_electrical
                                    + W_avionics
                                    + W_furnishings
                                    + W_air_conditioning
                                    + W_handling_gear
                                    + W_anti_ice
                                    + W_military_cargo_handling_system
                                    + Nen * Wen
                                      ) 
                    ]

# Weights
constraints += [
                W_maximum_gross >= W_operating_empty + Wc + W_fuel,
                zbre >= R * (SFC/V) * (D/L),
                W_fuel >= W_operating_empty*te_exp_minus1(zbre,4),
                W_landing_gross >= .7 * W_maximum_gross
               ]

# # Performance Model
constraints += [
               L == 1./2 * rho * V**2 * CL * S,
               L == W_maximum_gross,
               D == 1./2 * rho * V**2 * CD * S,
               S >= Saft + Scab + Sw ,
               ]

# Aerodynamic Model
constraints += [
                Re == rho*V*chord/mu,
                M == V/a,
                CDp0 >= 0.005,
                CDp**2 >= CDp**2 * M**2 + CDp0**2,
                CD >= CDp + 1.*(CL**2)/(pi*e*AR),
                ]


# Geometry
constraints += [
                Av                  == 4.0,
                Bw                  >= 200.0,
                Kz                  == 150.0,
                La                  >= 200.0,
                Lec                 >= 300.0,
                Lf                  >= 170.0,
                Lm                  >= 50.0,
                Ln                  >= 45.0,
                Lt                  == 150.0,
                NLt                 >= 12.0,
                Nw                  >= 8.0,
                Saft                >= 500.0,
                Scab                >= 2000.0,
                Sn                  >= np.pi * (Nw/2.)**2 * NLt,
                Svt                 == 800.,
                Sw                  >= 8000.,
                tau_w               <= .35,
                Vpr                 >= 3000.,
                Vt                  >= 2500.,
                Lambda_vt           == 45*np.pi/180,
                lambda_aft          >= 0.80,
               ]

objective = W_maximum_gross 

m=Model(objective, constraints)

sol=m.solve()

Weights_pie(sol)



