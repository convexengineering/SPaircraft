from gpkit import units
from numpy import cos, tan, pi

def get_M072_737_subs():
    """
    returns substitution dic for the a M 0.72 Boeing 737-800 model
    """
    sweep = 13.327 # [deg]
    VTsweep = 25. #[deg]
    HTsweep = 8. #[deg]
    M4a = .2
    fan = 1.685
    lpc  = 8./1.685
    hpc = 30./8.

    Mcruisemin = 0.72

    substitutions = {
            'N_{land}': 6.,
            'p_s': 81.*units('cm'),
            'numeng': 2.,
            'W_{avg. pass}': 180.*units('lbf'),
            'W_{carry on}': 15.*units('lbf'),
            'W_{checked}':40.*units('lbf'),
            'W_{fix}': 3000.*units('lbf'),
            'w_{aisle}': 0.51*units('m'),
            'w_{seat}': 0.5*units('m'),
            'w_{sys}': 0.1*units('m'),
            'r_E': 1.,  # [TAS]
            'p_{\\lambda_v}':1.6,
            '\\lambda_{cone}': 0.3,  # [TAS]
            '\\rho_{cone}': 2700.,#*units('kg/m^3'),  # [TAS]
            '\\rho_{bend}': 2700.,#*units('kg/m^3'),  # [TAS]
            '\\rho_{floor}': 2700.,#*units('kg/m^3'),  # [TAS]
            '\\rho_{skin}': 2700.,#*units('kg/m^3'),  # [TAS]
            '\\sigma_{floor}': 30000. / 0.000145, # [TAS] [Al]
            '\\sigma_{skin}': 15000. / 0.000145,  # [TAS] [Al]
            '\\sigma_{bend}': 30000. / 0.000145, # [TAS] [Al]
            '\\tau_{floor}': 30000. / 0.000145, # [TAS] [Al]
            'W\'\'_{floor}': 60.,  # [TAS]
            'W\'\'_{insul}': 22.,  # [TAS]
            'W\'_{window}': 145.*3.*units('N/m'),  # [TAS]
            'V_{mn}': 133.76*units('m/s'),

            # TASOPT Fuselage substitutions
            'L_{total/wing}': 1.127,
            'l_{nose}':20.*units('ft'),
            'numaisle': 1.,
            'SPR': 6.,
            'f_{seat}': 0.1,
            'W\'_{seat}': 1. * units('N'),  # Seat weight determined by weight fraction instead
            'W_{cargo}': 0.1*units('N'), # Cargo weight determined by W_{avg. pass_{total}}
            'W_{avg. pass_{total}}':215.*units('lbf'),
            'f_{string}': 0.35,
            'h_{floor}': 5. * units('in'),
            # 'R_{fuse}' : 1.715*units('m'),
            '\\delta_P_{over}': 8.382 * units('psi'),
            # fuselage subs that make fuse circular
            '\\delta R_{fuse}': 0.0001 * units('m'),
            '\\theta_{db}': 0.0001,

            # Power system and landing gear and engine weight fraction subs
            'f_{hpesys}': 0.01, # [TAS]
            'f_{lgmain}':0.044, # [TAS]
            'f_{lgnose}':0.011, # [TAS]
            'f_{pylon}': 0.10,

            # Fractional weights
            'f_{fadd}': 0.2,  # [TAS]
            'f_{frame}': 0.25,  # [Philippe]
            'f_{lugg,1}': 0.4,  # [Philippe]
            'f_{lugg,2}': 0.1,  # [Philippe]
            'f_{padd}': 0.35,  # [TAS]
            'f_{hpesys}': 0.01, # [TAS]
            'f_{lgmain}':0.03, # [TAS]
            'f_{lgnose}':0.0075, # [TAS]

            # Wing substitutions
            'C_{L_{wmax}}': 2.25/(cos(sweep)**2), # [TAS]
            '\\tan(\\Lambda)': tan(sweep * pi / 180.),
            '\\cos(\\Lambda)': cos(sweep * pi / 180.),
            '\\eta': 0.97,
            '\\rho_0': 1.225*units('kg/m^3'),
##            'AR': 10.1,
            'b_{max}': 117.5 * units('ft'),
            '\\rho_{fuel}': 817.*units('kg/m^3'),  # Kerosene [TASOPT]
            'f_{wingfuel}': .35,
            '\\tau_{max_w}': 0.1267,
            'TipReduct': 0.96,

            # Wing fractional weights
            'FuelFrac': 0.9,
            'f_{flap}': 0.2,
            'f_{slat}': 0.1,
            'f_{aileron}': 0.04,
            'f_{lete}': 0.1,
            'f_{ribs}': 0.15,
            'f_{spoiler}': 0.02,
            'f_{watt}': 0.03,

            # VT substitutions
            'C_{D_{wm}}': 0.5, # [2]
            'C_{L_{vmax}}': 2.6, # [TAS]
            'V_1': 70.*units('m/s'),
            '\\rho_{TO}': 1.225*units('kg/m^3'),
            'c_{l_{vtEO}}': 0.5, # [TAS]
            'e_v': 0.8,
            'V_{land}': 72.*units('m/s'),
            'N_{spar}': 1.,
            'f_{VT}': 0.4,
            'y_{eng}': 4.8768*units('m'),
            'numVT': 1.,
            'A_{vt}': 2.,
            '\\lambda_{vt}': 0.3,
            '\\tan(\\Lambda_{vt})': tan(VTsweep * pi / 180.),  # tangent of VT sweep
            'N_{spar}': 1.,
            '\\dot{r}_{req}': 0.001, # 10 deg/s/s yaw rate acceleration #NOTE: Constraint inactive
            '\\cos(\\Lambda_{vt})^3': cos(VTsweep * pi / 180.)**3,
            'c_{d_{fv}}': 0.0060,
            'c_{d_{pv}}': 0.0030,
            'V_{vt_{min}}': 0.1,
            
            # HT substitutions
            '\\alpha_{max,h}': 2.5,
            'f_{HT}': 0.3,
##            'AR_{ht}': 6.,
            '\\lambda_{ht}': 0.25,
            '\\tan(\\Lambda_{ht})': tan(HTsweep * pi / 180.),  # tangent of HT sweep
            'C_{L_{hmax}}': 2.0,  # [TAS]
            'C_{L_{hfcG}}': 0.7,
            '\\Delta x_{CG}': 7.68 * units('ft'),
            'x_{CG_{min}}': 56.75 * units('ft'),
            'SM_{min}': .31,
            '\\cos(\\Lambda_{ht})^3': cos(HTsweep * pi / 180.)**3,
            'c_{d_{fh}}': 0.0060,
            'c_{d_{ph}}': 0.0030,

            #engine system subs
            'f_{pylon}': 0.1,
            'f_{eadd}': 0.1,

            # Cabin air substitutions in AircraftP

            #set the fuel reserve fraction
            'ReserveFraction': .20,

            # Engine substitutions
            '\\pi_{tn}': .989,
            '\pi_{b}': .94,
            '\pi_{d}': .998,
            '\pi_{fn}': .98,
            'T_{ref}': 288.15,
            'P_{ref}': 101.325,
            '\eta_{HPshaft}': .99,
            '\eta_{LPshaft}': .978,
            'eta_{B}': .985,

            '\pi_{f_D}': fan,
            '\pi_{hc_D}': hpc,
            '\pi_{lc_D}': lpc,

##            '\\alpha_{OD}': 5.1,
##            '\\alpha_{max}': 5.6958,

            'hold_{4a}': 1. + .5 * (1.313 - 1.) * M4a ** 2.,
            'r_{uc}': .05,
            '\\alpha_c': .19036,
            'T_{t_f}': 435.,

            'M_{takeoff}': .9709,

            'G_f': 1.,

            'h_f': 43.003,

            'Cp_t1': 1253.,
            'Cp_t2': 1223.,
            'Cp_c': 1283.,

            'HTR_{f_SUB}': 1. - .3 ** 2.,
            'HTR_{lpc_SUB}': 1. - 0.6 ** 2.,

            # Minimum Cruise Mach Number
            'M_{min}': Mcruisemin,

            # engine system subs
            'rSnace': 16.,
            # nacelle drag calc parameter
            'r_{vnace}': 1.02,
            'T_{t_{4.1_{max}}}': 1833.*units('K'),
    }

    return substitutions
