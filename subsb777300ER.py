from gpkit import units
from numpy import cos, tan, pi

def getb777300ERsubs():
    """
    returns substitution dic for the TASOPT Boeing 777-300ER model
    """

    sweep = 32.583 # [deg]
    M4a = .1025
    fan = 1.5958
    lpc  = 1.26
    hpc = 20.033

    substitutions = {
            'N_{land}': 6.,
            'SPR': 10.,
            'p_s': 81.*units('cm'),
            'numeng': 2.,
            'numaisle':2.,
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

            # Fuselage subs
            'l_{nose}': 40*units('ft'),
            'f_{seat}': 0.1,
            'W\'_{seat}': 1. * units('N'),  # Seat weight determined by weight fraction instead
            'W_{cargo}': 0.1*units('N'), # Cargo weight determined by W_{avg. pass_{total}}
            'W_{avg. pass_{total}}':230.*units('lbf'),
            'f_{string}': 0.34,
            'h_{floor}': 5. * units('in'),
##            'R_{fuse}' : 3.0988*units('m'),
            '\\delta_P_{over}': 8.382 * units('psi'),

            # fuselage subs that make fuse circular
            '\\delta R_{fuse}': 0.0001 * units('m'),
            '\\theta_{db}': 0.0001,

            # TASOPT Fuselage substitutions
            'L_{total/wing}': 1.127,

            # Power system and landing gear and engine weight fraction subs
            'f_{hpesys}': 0.01, # [TAS]
            'f_{lgmain}':0.04, # [TAS]
            'f_{lgnose}':0.01, # [TAS]
            'f_{pylon}': 0.05,

            # Fractional weights
            'f_{fadd}': 0.2,  # [TAS]
            'f_{frame}': 0.24,  # [Philippe]
            'f_{lugg,1}': 0.4,  # [Philippe]
            'f_{lugg,2}': 0.4,  # [Philippe]
            'f_{padd}': 0.35,  # [TAS]

            # Wing substitutions
            'C_{L_{wmax}}': 2.25/(cos(sweep)**2), # [TAS]
            '\\tan(\\Lambda)': tan(sweep * pi / 180.),
            '\\cos(\\Lambda)': cos(sweep * pi / 180.),
            '\\eta': 0.97,
            '\\rho_0': 1.225*units('kg/m^3'),
           'b_{max}': 200 * units('ft'),

            '\\rho_{fuel}': 817.*units('kg/m^3'),  # Kerosene [TASOPT]
##            'AR': 8.455,

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
            'C_{D_{wm}}': 0.01, # [2] #NOTE: TAIL SIZED BY VOLUME COEFF!!
            'C_{L_{vmax}}': 2.6, # [TAS]
            'V_1': 70.*units('m/s'),
            '\\rho_{TO}': 1.225*units('kg/m^3'),
            '\\tan(\\Lambda_{vt})': tan(28*pi/180),
            'c_{l_{vtEO}}': 0.5, # [TAS]
            'e_v': 0.8,
            'V_{land}': 72.*units('m/s'),
            '\\dot{r}_{req}': 0.0001, # 10 deg/s/s yaw rate acceleration
            'y_{eng}': 32*units('ft'),
            'numVT': 1.,
           'A_{vt}': 2.35,
            '\\lambda_{vt}': 0.25,
            '\\tan(\\Lambda_{vt})': tan(28. * pi / 180.),  # tangent of VT sweep
            'N_{spar}': 1.,
            'f_{VT}': 0.4,
            'V_{vt}':0.06,

            # HT substitutions
            '\\alpha_{max,h}': 2.5,
           'AR_{ht}': 4.8,
            '\\lambda_{ht}': 0.32,
            '\\tan(\\Lambda_{ht})': tan(33. * pi / 180.),  # tangent of HT sweep
            'C_{L_{hmax}}': 2.0,  # [TAS]
            'C_{L_{hfcG}}': 0.7,
            '\\Delta x_{CG}': 11.97 * units('ft'),
            'x_{CG_{min}}': 117.31*units('ft'),
            'SM_{min}': .05,
            'f_{HT}': 0.3,

            #engine system subs
            'f_{pylon}': 0.05,
            'f_{eadd}': 0.1,

            #new engine params
            '\pi_{tn}': .995,
            '\pi_{b}': .94,
            '\pi_{d}': .995,
            '\pi_{fn}': .985,
            'T_{ref}': 288.15,
            'P_{ref}': 101.325,
            '\eta_{HPshaft}': .98,
            '\eta_{LPshaft}': .98,
            'eta_{B}': .9970,

            '\pi_{f_D}': fan,
            '\pi_{hc_D}': hpc,
            '\pi_{lc_D}': lpc,

            '\\alpha_{OD}': 8.62,
            '\\alpha_{max}': 8.62,

            'hold_{4a}': 1+.5*(1.313-1)*M4a**2,
            'r_{uc}': .01,
            '\\alpha_c': .156,
            'T_{t_f}': 435,

            'M_{takeoff}': .9539,

            'G_f': 1,

            'h_f': 43.003,

            'Cp_t1': 1280,
            'Cp_t2': 1184,
            'Cp_c': 1216,

            'HTR_{f_SUB}': 1-.3**2,
            'HTR_{lpc_SUB}': 1 - 0.6**2,

            'T_{t_{4.1_{max}}}':1828.5*units('K'),

            # Minimum Cruise Mach Number
            'M_{min}': 0.8,

            # engine system subs
            'rSnace': 16.,
            # nacelle drag calc parameter
            'r_{vnace}': 1.02,
           'T_{t_{4.1_{max}}}': 1860.*units('K'),
            'ReserveFraction': .05,
    }

    return substitutions
