from gpkit import units
from numpy import cos, tan, pi

def getD8bigsubs():
        """
        returns substitution dic for a 777esque D8
        """

        sweep = 13.237
        M4a = .1025
        fan = 1.58
        lpc  = 1.26
        hpc = 20.033

        substitutions = {
                'N_{land}': 6.,
                'p_s': 81.*units('cm'),
                '\\theta_{db}' : 0.366,
                'numeng': 2.,
                'numVT': 2.,
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

                # TASOPT Fuselage substitutions
                'l_{nose}': 29.*units('ft'),
                'L_{total/wing}': 1.179,

                # Fuselage subs
                'f_{seat}': 0.1,
                'W\'_{seat}': 1.,  # Seat weight determined by weight fraction instead
                'W_{cargo}': 0.1*units('N'), # Cargo weight determined by W_{avg. pass_{total}}
                'W_{avg. pass_{total}}':215.*units('lbf'),
                'f_{string}': 0.35,
                ##            'AR':15.749,
                'h_{floor}': 5.12*units('in'),
                ##            'R_{fuse}': 1.715*units('m'),
                ##            '\\delta R_{fuse}': 0.43*units('m'),
                'w_{db}': 0.93*units('m'),
                ##            'b_{max}': 140.0 * 0.3048*units('m'),
                # 'c_0': 17.4*0.3048,#units('ft'),
                '\\delta_P_{over}': 8.382 * units('psi'),
                'SPR': 10.,


                # Power system and landing gear subs
                'f_{hpesys}': 0.01, # [TAS]
                'f_{lgmain}':0.03, # [TAS]
                'f_{lgnose}':0.0075, # [TAS]
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
                ##        '\\alpha_{max,w}': 0.1,  # (6 deg)
                '\\cos(\\Lambda)': cos(sweep * pi / 180.),
                '\\eta': 0.97,
                '\\rho_0': 1.225*units('kg/m^3'),

                '\\rho_{fuel}': 817.*units('kg/m^3'),  # Kerosene [TASOPT]

                # Wing fractional weights
                'FuelFrac': 0.9,
                'f_{flap}': 0.2,
                'f_{slat}': 0.001,
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
                '\\tan(\\Lambda_{vt})': tan(40*pi/180),
                'c_{l_{vtEO}}': 0.5, # [TAS]
                'e_v': 0.8,
                # 'y_{eng}': 4.83*units('m'), # [3]
                'V_{land}': 72.*units('m/s'),
                '\\dot{r}_{req}': 0.1475, # 10 deg/s/s yaw rate acceleration
                'N_{spar}': 2.,
                'f_{VT}': 0.4,

                # HT substitutions
                '\\alpha_{max,h}': 2.5,
                'C_{L_{hmax}}': 1.225,#2.0, # [TAS]
                'SM_{min}': 0.05,
                '\\Delta x_{CG}': 2.0*units('m'),
                'x_{CG_{min}}' : 10.0*units('m'),
                'C_{L_{hfcG}}': 0.85,
                'f_{HT}': 0.3,

                # HT subs
                ##            'AR_{ht}': 12.,
                '\\lambda_{ht}': 0.3,
                '\\tan(\\Lambda_{ht})': tan(8. * pi / 180.),  # tangent of HT sweep
                '\\lambda_{vt}': 0.3,
                '\\tan(\\Lambda_{vt})': tan(25. * pi / 180.),  # tangent of VT sweep

                #engine system subs
                'rSnace': 6.,
                'f_{pylon}': 0.05,
                'f_{eadd}': 0.1,

                #nacelle drag calc parameter
                'r_{vnace}': 0.925,

                # Cabin air substitutions in AircraftP

                #set the fuel reserve fraction
                'ReserveFraction': .05,

                # Minimum Cruise Mach Number
                'M_{min}': 0.72,

                #new engine params
                '\pi_{tn}': .98,
                '\pi_{b}': .94,
                '\pi_{d}': .98,
                '\pi_{fn}': .98,
                'T_{ref}': 288.15,
                'P_{ref}': 101.325,
                '\eta_{HPshaft}': .98,
                '\eta_{LPshaft}': .98,
                'eta_{B}': .9970,

                '\pi_{f_D}': 1.58,
                '\pi_{hc_D}': 20.033,
                '\pi_{lc_D}': 1.26,

                '\\alpha_{OD}': 8.7877,
                '\\alpha_{max}': 8.7877,

                'hold_{4a}': 1+.5*(1.313-1)*M4a**2,#sol('hold_{4a}'),
                'r_{uc}': .01,
                '\\alpha_c': .14,
                'T_{t_f}': 435,

                'M_{takeoff}': .955,

                'G_f': 1,

                'h_f': 43.003,

                'Cp_t1': 1280,
                'Cp_t2': 1184,
                'Cp_c': 1216,

                'HTR_{f_SUB}': 1-.3**2,
                'HTR_{lpc_SUB}': 1 - 0.6**2,
        }

##                if not multimission:
##                    m.substitutions.update({
##        ##                 'n_{pax}': 180.,
##                         'ReqRng': 7370.*units('nmi'),
##                         })
        return substitutions
