# Below is a faithful reproduction of the weight breakdown models from Raymer.
# Note that this is NOT a GP compatible model and should not be viewed as such.

from gpkit import VectorVariable, Variable, Model, units
import numpy as np
from gpkit.tools import te_secant as secant
from gpkit.tools import te_tangent as tangent

Aw                  = Variable("A_{w}",                   "-",        "Wing Aspect Ratio")
Ah                  = Variable("A_{h}",                   "-",        "Horizontal Tail Aspect Ratio")
Av                  = Variable("A_{v}",                   "-",        "Vertical Tail Aspect Ratio")
Bh                  = Variable("B_{h}",                   "ft",       "Horizontal Tail Span")
Bw                  = Variable("B_{w}",                   "ft",       "Wing Span")
D                   = Variable("D",                       "ft",       "Fuselage Structural Depth")
Fw                  = Variable("F_{w}",                   "ft",       "Fuselage Width at Horizontal Tail Intersection")
Ht_Hv               = Variable("\\frac{H_t}{h_v}",        "-",        "0.0 for Conventional Tail; 1.0 for T Tail")
Iy                  = Variable("I_{y}",                   "lb*ft^2",  "Yawing Moment of Inertia")
K_door              = Variable("K_{door}",                "-",        "1.0 if no cargo door; 1.06 if one side cargo door; " +
                                                                      "1.12 if two side cargo doors; " +
                                                                      "1.2 if aft clamshell door; " +
                                                                      "1.25 if two side cargo doors and aft clamshell door")
KLg                 = Variable("K_{Lg}",                  "-",        "1.12 if fuselage-mounted main landing gear; 1.0 Otherwise")
Kmp                 = Variable("K_{mp}",                  "-",        "1.126 for kneeling gear; 1.0 Otherwise")
Kng                 = Variable("K_{ng}",                  "-",        "1.107 for pylon mounted nacelle; 1.0 Otherwise")
Knp                 = Variable("K_{np}",                  "-",        "1.15 for kneeling gear; 1.0 Otherwise")
Kp                  = Variable("K_{p}",                   "-",        "1.4 for Engine With Propeller; 1.0 Otherwise")
Kr                  = Variable("K_{r}",                   "-",        "1.133 if Reciprocating Engine; 1.0 Otherwise")
Ktpg                = Variable("K_{tpg}",                 "-",        "0.826 for Tripod (A-7) Gear; 1.0 Otherwise")
Ktr                 = Variable("K_{tr}",                  "-",        "1.18 for Jet with Thrust Reverser; 1.0 Otherwise")
Kuht                = Variable("K_{uht}",                 "-",        "1.143 for Unit (All-Moving) Horizontal Tail; 1.0 Otherwise")
Kws                 = Variable("K_{ws}",                  "-",        "0.75[(1+2lambda)/(1+lambda)](Bw tanLambda / L")
Ky                  = Variable("K_{y}",                   "ft",       "Aircraft Pitching Radius of Gyration (\\approx 0.3 L_{t")
Kz                  = Variable("K_{z}",                   "ft",       "Aircraft Yawing Radius of Gyration (\\approx L_{t")
L                   = Variable("L",                       "ft",       "Fuselage Structural Length (excludes radome cowling, tail cap")
La                  = Variable("L_{a}",                   "ft",       "Electrical Routing Distance, Generators to Avionics to Cockpit")
Lf                  = Variable("L_{f}",                   "ft",       "Total Fuselage Length")
Lm                  = Variable("L_{m}",                   "in",       "Extended Length of Main Landing Gear")
Ln                  = Variable("L_{n}",                   "in",       "Extended Nose Gear Length")
Lt                  = Variable("L_{t}",                   "ft",       "Tail Length; Wing Quarter-MAC to Tail Quarter-MAC")
Nc                  = Variable("N_{c}",                   "-",        "Number of Crew")
Nen                 = Variable("N_{en}",                  "-",        "Number of Engines")
Nf                  = Variable("N_{f}",                   "-",        "Number of Functions Performed by the Controls (Typically 4-7")
Ngen                = Variable("N_{gen}",                 "-",        "Number of Generators (Typically = N_{en")
NLt                 = Variable("N_{Lt}",                  "ft",       "Nacelle Length")
Nl                  = Variable("N_{l}",                   "-",        "Ultimate Landing Loading Factor; = N_{gear}x1.5")
Nm                  = Variable("N_{m}",                   "-",        "Number of Mechanical Functions (Typically 0-2")
Nmss                = Variable("N_{mss}",                 "-",        "Number of Main Gear Shock Struts")
Nmw                 = Variable("N_{mw}",                  "-",        "Number of Main Wheels")
Np                  = Variable("N_{p}",                   "-",        "Number of Personnel Onboard (Crew and Passengers")
Nt                  = Variable("N_{t}",                   "-",        "Number of Fuel Tanks")
Nw                  = Variable("N_{w}",                   "ft",       "Nacelle Width")
Nz                  = Variable("N_{z}",                   "-",        "Ultimate Load Factor; = 1.5 x Limit Load Factor")
Rkva                = Variable("R_{kva}",                 "kV*A",     "System Electrical Rating (Typically 40-60) ")
Scs                 = Variable("S_{cs}",                  "ft^2",     "Total Area of Control Surfaces")
Scsw                = Variable("S_{csw}",                 "ft^2",     "Control Surface Area (Wing Mounted")
Se                  = Variable("S_{e}",                   "ft^2",     "Elevator Area")
Sf                  = Variable("S_{f}",                   "ft^2",     "Fuselage Wetted Area")
S_floor             = Variable("S_{floor}",               "ft^2",     "Cargo Floor Area")
Sht                 = Variable("S_{ht}",                  "ft^2",     "Horizontal Tail Area")
Sn                  = Variable("S_{n}",                   "ft^2",     "Nacelle Wetted Area")
Svt                 = Variable("S_{vt}",                  "ft^2",     "Vertical Tail Area")
Sw                  = Variable("S_{w}",                   "ft^2",     "Trapezoidal Wing Area")
t_c_root_w          = Variable("\\frac{t}{c}_{root_w}",   "-",        "Thickness to Chord Ratio of Main Wing at the Root")
t_c_root_v          = Variable("\\frac{t}{c}_{root_v}",   "-",        "Thickness to Chord Ratio of Vertical Tail at the Root")
Vi                  = Variable("V_{i}",                   "gal",      "Volume of Integral Tanks")
Vp                  = Variable("V_{p}",                   "gal",      "Volume of Self Sealing (Protected) Tanks")
Vpr                 = Variable("V_{pr}",                  "ft^3",     "Volume of Pressurized Section")
V_stall             = Variable("V_{stall}",               "m/s",      "Stall Speed")
Vt                  = Variable("V_{t}",                   "gal",      "Total Fuel Volume")
W_APU_uninstalled   = Variable("W_{APU_{uninstalled}}",   "lb",       "Uninstalled APU Weight")
Wc                  = Variable("W_{c}",                   "lb",       "Maximum Cargo Weight")
Wdg                 = Variable("W_{dg}",                  "lb",       "Flight Design Gross Weight")
Wec                 = Variable("W_{ec}",                  "lb",       "Weight of Engine and Contents (per nacelle); \\approx 2.331 W_en^0.901 K_p K_tr")
Wl                  = Variable("W_{l}",                   "lb",       "Landing Design Gross Weight")
Wuav                = Variable("W_{uav}",                 "-",        "Uninstalled Avionics Weight (Typically 800-1400")
Lambda_w            = Variable("\\Lambda_{w}",            "-",        "Main Wing Sweep at 0.25 MAC")
Lambda_ht           = Variable("\\Lambda_{ht}",           "-",        "Horizontal Tail Sweep at 0.25 MAC")
Lambda_vt           = Variable("\\Lambda_{vt}",           "-",        "Vertical Tail Sweep at 0.25 MAC")
lambda_w            = Variable("\\lambda_{w}",            "-",        "Main Wing Taper Ratio")

W_wing                            = Variable("W_{wing}",                   "lbf", "Weight of Main Wing")
W_horizontal_tail                 = Variable("W_{horizontal-tail}",        "lbf", "Weight of Horizontal Tail")
W_vertical_tail                   = Variable("W_{vertical-tail}",          "lbf", "Weight of Vertical Tail")
W_fuselage                        = Variable("W_{fuselage}",               "lbf", "Weight of Fuselage")
W_main_landing_gear               = Variable("W_{main-gear}",              "lbf", "Weight of Main Landing Gear")
W_nose_landing_gear               = Variable("W_{nose-gear}",              "lbf", "Weight of Nose Landing Gear")
W_nacelle_group                   = Variable("W_{nacelle}",                "lbf", "Weight of Nacelle Group")
W_engine_controls                 = Variable("W_{engine-controls}",        "lbf", "Weight of Engine Controls")
W_starter_penumatic               = Variable("W_{starter_{penumatic}}",    "lbf", "Weight of Pneumatic Starter")
W_fuel_system                     = Variable("W_{fuel-system}",            "lbf", "Weight of Fuel System")
W_flight_controls                 = Variable("W_{flight-controls}",        "lbf", "Weight of Flight Controls")
W_APU_installed                   = Variable("W_{APU_{installed}}",        "lbf", "Weight of Installed APU")
W_instruments                     = Variable("W_{instruments}",            "lbf", "Weight of Instruments")
W_hydraulics                      = Variable("W_{hydraulics}",             "lbf", "Weight of Hydraulic System")
W_electrical                      = Variable("W_{electrical}",             "lbf", "Weight of Electrical System")
W_avionics                        = Variable("W_{avionics}",               "lbf", "Weight of Avionics")
W_furnishings                     = Variable("W_{furnishings}",            "lbf", "Weight of Furnishings (No seats or handling gear")
W_air_conditioning                = Variable("W_{air-conditioning}",       "lbf", "Weight of Air Conditioning System")
W_anti_ice                        = Variable("W_{anti-ice}",               "lbf", "Weight of Anti-Icing System")
W_handling_gear                   = Variable("W_{handling-gear}",          "lbf", "Weight of Cargo Handling Gear")
W_military_cargo_handling_system  = Variable("W_{military-cargo-system}",  "lbf", "Weight of Military Cargo Handling System")



constraints = [ # Raymer page 459
                W_wing >= 0.0051 * (Wdg*Nz)**(0.557) * Sw**(0.649) * Aw**(0.5) * t_c_root_w**(-0.4) * (1+lambda_w)**(0.1) * secant(Lambda_w, 6) * Scsw**(0.1),
                W_horizontal_tail >= 0.0379 * Kuht * (1 + Fw/Bh)**(-0.25) * Wdg**(0.639) * Nz**(0.1) * Sht**(0.75) * Lt**(-1.0) * Ky**(0.704) * secant(Lambda_ht, 6) * Ah**(0.166) * (1 + Se/Sht)**(0.1),
                W_vertical_tail >= 0.0026 * (1 + Ht_Hv)**(0.225) * Wdg**(0.556) * Nz**(0.536) * Lt**(-0.5) * Svt**(0.5) * Kz**(0.875) * secant(Lambda_vt, 6) * Av**(0.35) * t_c_root_vt**(-0.5),
                Kws >= 0.75 * (1+2*lambda_w)/(1+lambda_w) * Bw * tangent(Lambda_w) / L,
                W_fuselage >= 0.3280 * K_door * KLg * (Wdg*Nz)**(0.5) * L**(0.25) * Sf**(0.302) * (1 + Kws)**(0.04) * (L/D)**(0.10),
                W_main_landing_gear >= 0.0106 * Kmp * Wl**(0.888) * Nl**(0.25) * Lm**(0.4) * Nmw**(0.321) * Nmss**(-0.5) * V_stall**(0.1),
                W_nose_landing_gear >= 0.032 * Knp * Wl**(0.646) * Nl**(0.2) * Ln**(0.5) * Nnw**(0.45),
                Wec == 2.331 * W_en**(0.901) * K_p * K_tr,
                W_nacelle_group >= 0.6724 * Kng * NLt**(0.10) * Nw**(0.294) * Nz**(0.119) * Wec**(0.611) * Nen**(0.984) * Sn**(0.224),
                W_engine_controls >= 5.0 * Nen + 0.80 * Lec,
                W_starter_penumatic >= 49.19 * (Nen * Wen / 1000.)**(0.541),
                W_fuel_system >= 2.405 * Vt**(0.606) * (1 + Vi/Vt)**(-1.0) * (1 + Vp/Vt) * Nt**(0.5),
                W_flight_controls >= 145.9 * Nf**(0.554) * (1 + Nm/Nf)**(-1.0) * Scs**(0.20) * (Iy*10**(-6))**(0.07),
                W_APU_installed >= 2.2 * W_APU_uninstalled,
                W_instruments >= 4.509 * Kr * Ktp * Nc**(0.541) * Nen * (Lf + Bw)**(0.5),
                W_hydraulics >= 0.2673 * Nf * (Lf + Bw)**(0.937),
                W_electrical >= 7.291 * Rkva**(0.782) * La**(0.346) * Ngen**(0.10),
                W_avionics >= 1.73 * Wuav**(0.983),
                W_furnishings >= 0.0577 * Nc**(0.1) * Wc**(0.393) * Sf**(0.75),
                W_air_conditioning >= 62.36 * Np**(0.25) * (Vpr/1000.)**(0.604) * Wuav**(0.01),
                W_handling_gear >= 0.002 * Wdg,
                W_anti_ice >= 0.0003 * Wdg,
                W_military_cargo_handling_system  >= 2.4 * S_floor,
              ]

objective = Wdg

m=Model(objective, constraints)





