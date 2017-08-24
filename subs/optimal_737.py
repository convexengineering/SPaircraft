from gpkit import units
from numpy import cos, tan, pi
import numpy as np

def get737_optimal_subs():
    """
    returns substitution dic for the TASOPT Boeing 737-800 model
    """
    sweep = 26.0 # [deg]
    VTsweep = 25. #[deg]
    HTsweep = 25. #[deg]
    M4a = .2
    fan = 1.60474
    lpc  = 4.98
    hpc = 35./8.

    Mcruisemin = 0.8

    substitutions = {
            'N_{land}': 6.,
            'p_s': 81.*units('cm'),
            'n_{eng}': 2.,
            'W_{avg. pass}': 180.*units('lbf'),
            'W_{carry on}': 15.*units('lbf'),
            'W_{checked}':40.*units('lbf'),
            'W_{fix}': 3000.*units('lbf'),
            'w_{aisle}': 0.51*units('m'),
            'w_{seat}': 0.5*units('m'),
            'w_{sys}': 0.1*units('m'),
            'r_E': 1.,  # [TAS]
            'p_{\\lambda_{vt}}':1.6,
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
            'V_{mn}': 133.76*units('m/s'), 'V_{ne}':143.92*units('m/s'),

            'D_{reduct}': 1,

            # TASOPT Fuselage substitutions
            'f_{L_{total/wing}}': 1.127,
            'l_{nose}':20.*units('ft'),
            'n_{aisle}': 1.,
            'SPR': 6.,
            'f_{seat}': 0.1,
            'W\'_{seat}': 1. * units('N'),  # Seat weight determined by weight fraction instead
            'W_{cargo}': 0.1*units('N'), # Cargo weight determined by W_{avg. pass_{total}}
            'W_{avg. pass_{total}}':215.*units('lbf'),
            'f_{string}': 0.35,
            'h_{floor}': 5. * units('in'),
            # 'R_{fuse}' : 1.715*units('m'),
            '\\Delta P_{over}': 8.382 * units('psi'),
            # fuselage subs that make fuse circular
            '\\Delta R_{fuse}': 0.0001 * units('m'),
            '\\theta_{db}': 0.0001,

            # Power system and landing gear and engine weight fraction subs
            'f_{hpesys}': 0.01, # [TAS]

            # Fractional weights
            'f_{fadd}': 0.2,  # [TAS]
            'f_{frame}': 0.25,  # [Philippe]
            'f_{lugg,1}': 0.4,  # [Philippe]
            'f_{lugg,2}': 0.1,  # [Philippe]
            'f_{padd}': 0.35,  # [TAS]
            'f_{hpesys}': 0.01, # [TAS]

            # Wing substitutions
            'C_{L_{w,max}}': 2.25/(cos(sweep* pi / 180.)**2), # [TAS]
            '\\tan(\\Lambda)': tan(sweep * pi / 180.),
            '\\cos(\\Lambda)': cos(sweep * pi / 180.),
            '\\eta': 0.97,
            '\\rho_0': 1.225*units('kg/m^3'),
##            'AR': 10.1,
            'b_{max}': 117.5 * units('ft'),
            '\\rho_{fuel}': 817.*units('kg/m^3'),  # Kerosene [TASOPT]
            'f_{wingfuel}': .5,  #.35 may be a better answer
            '\\tau_{max_w}': 0.1267,
            'TipReduct': 1,

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
            'C_{L_{vt,max}}': 2.6, # [TAS]
            'V_1': 70.*units('m/s'),
            '\\rho_{TO}': 1.225*units('kg/m^3'),
            'c_{l_{vt,EO}}': 0.5, # [TAS]
            'e_{vt}': 0.8,
            'V_{land}': 72.*units('m/s'),
            'N_{spar}': 1.,
            'f_{VT}': 0.4,
            'y_{eng}': 4.8768*units('m'),
            'n_{vt}': 1.,
            'A_{vt}': 2.,
            '\\lambda_{vt}': 0.3,
            '\\tan(\\Lambda_{vt})': tan(VTsweep * pi / 180.),  # tangent of VT sweep
            'N_{spar}': 1.,
            '\\dot{r}_{req}': 0.0001, # 10 deg/s/s yaw rate acceleration #NOTE: Constraint inactive
            '\\cos(\\Lambda_{vt})^3': cos(VTsweep * pi / 180.)**3,
            'c_{d_{fv}}': 0.0060,
            'c_{d_{pv}}': 0.0030,
            'V_{vt_{min}}': 0.001, #0.1,
            'Fsafetyfac': 1.8,
            
            # HT substitutions
            '\\alpha_{ht,max}': 2.5,
            'f_{ht}': 0.3,
##            'AR_{ht}': 6.,
            '\\lambda_{ht}': 0.25,
            '\\tan(\\Lambda_{ht})': tan(HTsweep * pi / 180.),  # tangent of HT sweep
            'C_{L_{ht,max}}': 2.0,  # [TAS]
            'C_{L_{ht,fCG}}': 0.7,
            '\\Delta x_{CG}': 7.68 * units('ft'),
            'x_{CG_{min}}': 56.75 * units('ft'),
            'SM_{min}': .15,
            '\\cos(\\Lambda_{ht})^3': cos(HTsweep * pi / 180.)**3,
            'c_{d_{fh}}': 0.0060,
            'c_{d_{ph}}': 0.0030,
            '\\eta_{ht}': 0.9, 

            #engine system subs
            'f_{pylon}': 0.1,
            'f_{eadd}': 0.1,

            # Cabin air substitutions in AircraftP

            #set the fuel reserve fraction
            'ReserveFraction': .20,

            #min altitude for start of cruise
            'MinCruiseAlt': 35000*units('ft'),

            # Engine substitutions
            '\\pi_{tn}': .995,
            '\\pi_{b}': .94,
            '\\pi_{d}': .995,
            '\\pi_{fn}': .985,
            'T_{ref}': 288.15,
            'P_{ref}': 101.325,
            '\\eta_{HPshaft}': .978,
            '\\eta_{LPshaft}': .99,
            '\\eta_{B}': .985,

            '\\pi_{f_D}': fan,
            '\\pi_{hc_D}': hpc,
            '\\pi_{lc_D}': lpc,

##            '\\alpha_{OD}': 6.97,
##            '\\alpha_{max}': 6.97,

            'hold_{4a}': 1.+.5*(1.313-1.)*M4a**2.,
            'r_{uc}': .01,
            '\\alpha_c': .16,
            'T_{t_f}': 435.,

            'M_{takeoff}': .9556,

            'G_{f}': 1.,

            'h_{f}': 42.5,

            'C_{p_{t1}}': 1236.5,
            'C_{p_{t2}}': 1200.4,
            'C_{p_{c}}': 1253.9,

            'HTR_{f_{SUB}}': 1. - .3 ** 2.,
            'HTR_{lpc_{SUB}}': 1. - 0.6 ** 2.,

            # Minimum Cruise Mach Number
            'M_{min}': Mcruisemin,

            # engine system subs
            'r_{S_{nacelle}}': 16.,
            # nacelle drag calc parameter
            'r_{v_{nacelle}}': 1.02,
            'T_{t_{4.1_{max}}}': 1567.*units('K'),#('sweep', np.linspace(1250, 1700, 11)),#

            'T_{t_{4.1_{max-Cruise}}}': 1125.*units('K'),

          'MaxClimbTime': 13*units('min'),
##          'MaxClimbDistance': 180*units('nautical_miles')

            #LG subs
            'E': 205,
            'K': 2,
            'N_s': 2,
            '\\eta_s': 0.8,
            '\\lambda_{LG}': 2.5,
            '\\rho_{st}': 7850,
            '\\tan(\\gamma)': np.tan(5*np.pi/180),
            '\\tan(\\phi_{min})': np.tan(15*np.pi/180),
            '\\tan(\\psi_{max})': np.tan(63*np.pi/180),
            '\\tan(\\theta_{max})': np.tan(15*np.pi/180),
            '\\sigma_{y_c}': 470E6,
            'f_{add,m}': 1.5,
            'f_{add,n}': 1.5,
            'h_{hold}': 1,
            'h_{nacelle}': 0.5,
            'n_{mg}': 2,
            'n_{wps}': 2,
            'p_{oleo}': 1800,
            't_{nacelle}': 0.15,
            'w_{ult}': 10,
            'z_{CG}': 2,
            'z_{wing}': 0.5,
    }

    return substitutions
