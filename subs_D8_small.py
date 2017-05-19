from gpkit import units
from numpy import cos, tan, pi

def get_small_D8_subs():
        """
        returns substitution dic for the D8.2
        """
        sweep = 13.237  # [deg]
        VTsweep = 25.0 #[deg]
        HTsweep = 8.0 #[deg]
        M4a = .2
        fan = 1.60474
        lpc  = 4.98
        hpc = 35./8.

        #Percent of velocity loss from BL smeared across entire fan
        BLIVloss = 0.088197
        #Min cruise mach number
        Mcruisemin = 0.72

        fBLI = 0.4
        wake_benefit = 0.1
        wake_percent = 0.33

        substitutions = {
                'N_{land}': 6.,
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
                'V_{mn}': 133.76*units('m/s'),

                #BLI drag reduction factor
                'D_{reduct}': 1-((1-fBLI)*wake_benefit*wake_percent),

                # Fuselage subs
                'f_{seat}': 0.1,
                'W\'_{seat}': 1.,  # Seat weight determined by weight fraction instead
                'W_{cargo}': 10000*units('lbf'), # Cargo weight determined by W_{avg. pass_{total}}
                'W_{avg. pass_{total}}':215.*units('lbf'),
                'f_{string}': 0.35,

                'h_{floor}': 5.12*units('in'),
##                'R_{fuse}': 1.715*units('m'),
##                '\\delta R_{fuse}': 0.43*units('m'),
                'w_{db}': 0.93*units('m'),
                '\\delta_P_{over}': 8.382 * units('psi'),
                'SPR': 6.,

                # TASOPT Fuselage substitutions
                'l_{nose}': 18.*units('ft'),
                'L_{total/wing}': 1.179,

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
                'C_{L_{wmax}}': 2.15/(cos(sweep)**2), # [TAS]
                '\\tan(\\Lambda)': tan(sweep * pi / 180.),
                '\\cos(\\Lambda)': cos(sweep * pi / 180.),
                '\\eta': 0.97,
                '\\rho_0': 1.225*units('kg/m^3'),
                '\\rho_{fuel}': 817.*units('kg/m^3'),  # Kerosene [TASOPT]
##                'AR':15.749,
                'b_{max}': 140.0 * 0.3048*units('m'),
                '\\tau_{max_w}': 0.14733,
                'TipReduct': 1.0,
                'f_{wingfuel}': .5,

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
                'c_{l_{vtEO}}': 0.5, # [TAS]
                'e_v': 0.8,
                'V_{land}': 72.*units('m/s'),
                '\\dot{r}_{req}': 0.0001, #0.1475 10 deg/s/s yaw rate acceleration
                'N_{spar}': 1,
                'f_{VT}': 0.4,
                'numVT': 2.,
                'A_{vt}' : 2.2,
                '\\lambda_{vt}': 0.3,
                '\\tan(\\Lambda_{vt})': tan(VTsweep * pi / 180.),  # tangent of VT sweep
                '\\cos(\\Lambda_{vt})^3': cos(VTsweep * pi / 180.)**3,
                'c_{d_{fv}}': 0.0060,
                'c_{d_{pv}}': 0.0035,
                'V_{vt_{min}}': 0.03,

                # HT substitutions
                '\\alpha_{max,h}': 2.5,
                'C_{L_{hmax}}': 2.0, # [TAS]
                'SM_{min}': 0.05,
                '\\Delta x_{CG}': 6.8*units('ft'),
                'x_{CG_{min}}' : 56.02*units('ft'),
                'C_{L_{hfcG}}': 0.85,
                'f_{HT}': 0.3,
##                'AR_{ht}': 12.,
                '\\lambda_{ht}': 0.3,
                '\\tan(\\Lambda_{ht})': tan(HTsweep * pi / 180.),  # tangent of HT sweep
                '\\cos(\\Lambda_{ht})^3': cos(HTsweep * pi / 180.)**3,
                'c_{d_{fh}}': 0.0060,
                'c_{d_{ph}}': 0.0035,
                
                #engine system subs
                'rSnace': 6.,
                'f_{pylon}': 0.05,
                'f_{eadd}': 0.1,

                #nacelle drag calc parameter
                'r_{vnace}': 0.925,

                # Cabin air substitutions in AircraftP

                #set the fuel reserve fraction
                'ReserveFraction': .20,

                # Minimum Cruise Mach Number
                'M_{min}': Mcruisemin,

                # Engine substitutions
                '\\pi_{tn}': .995,
                '\pi_{b}': .94,
                '\pi_{d}': .995,
                '\pi_{fn}': .985,
                'T_{ref}': 288.15,
                'P_{ref}': 101.325,
                '\eta_{HPshaft}': .978,
                '\eta_{LPshaft}': .99,
                'eta_{B}': .985,

                '\pi_{f_D}': fan,
                '\pi_{hc_D}': hpc,
                '\pi_{lc_D}': lpc,

##                '\\alpha_{OD}': 6.97,
##                '\\alpha_{max}': 6.97,

                'hold_{4a}': 1.+.5*(1.313-1.)*M4a**2.,
                'r_{uc}': .01,
                '\\alpha_c': .16,
                'T_{t_f}': 435.,

                'M_{takeoff}': .9556,

                'G_f': 1.,

                'h_f': 43.003,

                'Cp_t1': 1236.5,
                'Cp_t2': 1200.4,
                'Cp_c': 1257.9,

                'HTR_{f_SUB}': 1.-.3**2.,
                'HTR_{lpc_SUB}': 1. - 0.6**2.,

                'T_{t_{4.1_{max}}}': 1567.*units('K'),

                #BLI factors
                #compute the cruise stagnation pressure loss factor given a min
                #mach number and the BLI velocity loss factor
                'f_{BLI_P}': (10400. +.2546*1.225*((1-BLIVloss)*Mcruisemin*295.)**2)/ \
                              (10400. +.2546*1.225*(Mcruisemin*295.)**2),
                'f_{BLI_V}': 1 - BLIVloss,
        }

        return substitutions
