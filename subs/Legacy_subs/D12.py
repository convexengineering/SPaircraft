from gpkit import units
from numpy import cos, tan, pi

def get_D12_subs():
        """
        returns substitution dic for a D12
        """

        sweep = 32.0 # [deg]
        VTsweep = 25. #[deg]
        HTsweep = 20. #[deg]
        M4a = .1025
        fan = 1.65
        lpc  = 8./fan
        hpc = 48./(fan*lpc)

        #Percent of velocity loss from BL smeared across entire fan
        BLIVloss = 0.072712
        #Min cruise mach number
        Mcruisemin = 0.83

        fBLI = 0.4
        wake_benefit = 0.08
        wake_percent = 0.33

        substitutions = {
                'N_{land}': 6.,
                'p_s': 81.*units('cm'),
                '\\theta_{db}' : 0.366,
                'n_{eng}': 2.,
                'n_{vt}': 2.,
                'n_{aisle}':2.,
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

                #BLI drag reduction factor
                'D_{reduct}': 1-((1-fBLI)*wake_benefit*wake_percent),

                # TASOPT Fuselage substitutions
                'l_{nose}': 29.*units('ft')*1.75,   #1.75 is estimated length increase factor from D8.2
                'f_{L_{total/wing}}': 1.179,

                # Fuselage subs
                'f_{seat}': 0.1,
                'W\'_{seat}': 1.,  # Seat weight determined by weight fraction instead
                'W_{cargo}': 0.1*units('N'), # Cargo weight determined by W_{avg. pass_{total}}
                'W_{avg. pass_{total}}':230.*units('lbf'),
                'f_{string}': 0.35,
                'h_{floor}': 5 *units('in'),
                'R_{fuse}': 2.1*units('m'),
##                '\\Delta R_{fuse}': 0.43*units('m'),
                'w_{db}': 1.5*units('m'),
                '\\Delta P_{over}': 8.382 * units('psi'),
                'SPR': 12.,

                # Power system and landing gear subs
                'f_{hpesys}': 0.01, # [TAS]
                'f_{lgmain}':0.03, # [TAS]
                'f_{lgnose}':0.0075, # [TAS]
                'f_{pylon}': 0.05,

                # Fractional weights
                'f_{fadd}': 0.2,  # [TAS]
                'f_{frame}': 0.25,  # [Philippe]
                'f_{lugg,1}': 0.4,  # [Philippe]
                'f_{lugg,2}': 0.4,  # [Philippe]
                'f_{padd}': 0.35,  # [TAS]
                'f_{hpesys}': 0.01, # [TAS]
                'f_{lgmain}':0.03, # [TAS]
                'f_{lgnose}':0.0075, # [TAS]

                # Wing substitutions
                'AR':10,
                'b_{max}': 220.0 * units('ft'),
                'C_{L_{w,max}}': 2.25/(cos(sweep)**2), # [TAS]
                '\\tan(\\Lambda)': tan(sweep * pi / 180.),
                '\\cos(\\Lambda)': cos(sweep * pi / 180.),
                '\\eta': 0.97,
                '\\rho_0': 1.225*units('kg/m^3'),

                '\\rho_{fuel}': 817.*units('kg/m^3'),  # Kerosene [TASOPT]
                '\\tau_{max_w}': 0.145,
                'f_{wingfuel}': 1.0,
                'TipReduct': 0.945,

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
                'A_{vt}' : 2.0,
                'C_{D_{wm}}': 0.5, # [2]
                'C_{L_{vt,max}}': 2.6, # [TAS]
                'V_1': 70.*units('m/s'),
                '\\rho_{TO}': 1.225*units('kg/m^3'),
                '\\tan(\\Lambda_{vt})': tan(VTsweep*pi/180),
                'c_{l_{vt,EO}}': 0.5, # [TAS]
                'e_{vt}': 0.8,
                # 'y_{eng}': 4.83*units('m'), # [3]
                'V_{land}': 72.*units('m/s'),
                '\\dot{r}_{req}': 0.00001, # 10 deg/s/s yaw rate acceleration
                'N_{spar}': 1.,
                'f_{VT}': 0.4,
                '\\cos(\\Lambda_{vt})^3': cos(VTsweep * pi / 180.)**3,
                'c_{d_{fv}}': 0.0060,
                'c_{d_{pv}}': 0.0035,
                'V_{vt_{min}}': 0.04,
                '\\lambda_{vt}': 0.3,

                # HT substitutions
                '\\alpha_{ht,max}': 2.5,
                'C_{L_{ht,max}}': 2.0, # [TAS]
                'SM_{min}': 0.05,
                '\\Delta x_{CG}': 11.97 * units('ft'),
                'x_{CG_{min}}' : 117.31*units('ft'),
                'C_{L_{ht,fCG}}': 0.85,
                'f_{ht}': 0.3,
                '\\cos(\\Lambda_{ht})^3': cos(HTsweep * pi / 180.)**3,
                'c_{d_{fh}}': 0.0060,
                'c_{d_{ph}}': 0.0035,
                'AR_{ht}': 10.,
                '\\lambda_{ht}': 0.3,
                '\\tan(\\Lambda_{ht})': tan(HTsweep * pi / 180.),  # tangent of HT sweep
                

                #engine system subs
                'f_{pylon}': 0.05,
                'f_{eadd}': 0.1,

                # Cabin air substitutions in AircraftP

                #set the fuel reserve fraction
                'f_{fuel_{res}}': .05,

                # Minimum Cruise Mach Number
                'M_{min}': Mcruisemin,

                #new engine params
                '\pi_{tn}': .995,
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

##                '\\alpha_{OD}': 8.62,
##                '\\alpha_{max}': 8.62,     #place holder value not active during analysis

                'hold_{4a}': 1+.5*(1.313-1)*M4a**2,
                'r_{uc}': .01,
                '\\alpha_c': 0.14757,
                'T_{t_f}': 435,

                'M_{takeoff}': .93997,

                'G_{f}': 1,

                'h_{f}': 43.003,

                'C_{p_{t1}}': 1257.3,
                'C_{p_{t2}}': 1158.35,
                'C_{p_{c}}': 1278.5,

                'HTR_{f_{SUB}}': 1-.3**2,
                'HTR_{lpc_{SUB}}': 1 - 0.6**2,

                # engine system subs
                'r_{S_{nacelle}}': 8.,
                # nacelle drag calc parameter
                'r_{v_{nacelle}}': .925,
                'T_{t_{4.1_{max}}}': 1860.*units('K'),
                
                #BLI factors
                #compute the cruise stagnation pressure loss factor given a min
                #mach number and the BLI velocity loss factor
                'f_{BLI_{P}}': (18753.9 +.5*0.301559*((1-BLIVloss)*Mcruisemin*295.)**2)/ \
                              (18753.9 +.5*0.301559*(Mcruisemin*295.)**2),
                'f_{BLI_{V}}': 1 - BLIVloss,
        }

        return substitutions
