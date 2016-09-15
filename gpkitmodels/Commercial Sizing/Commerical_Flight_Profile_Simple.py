"""Simple commercial aircraft flight profile and aircraft model"""
from numpy import pi
import gpkit
import numpy as np
from gpkit import VectorVariable, Variable, Model, units, ConstraintSet, LinkedConstraintSet, SignomialEquality
from gpkit.tools import te_exp_minus1
from gpkit.constraints.tight import TightConstraintSet as TCS
import matplotlib.pyplot as plt
from gpkit.small_scripts import mag
from atm_test import Atmosphere

#packages just needed for plotting since this is for sweeps
import matplotlib.pyplot as plt

#only needed for the local bounded debugging tool
from collections import defaultdict

"""
minimizes the aircraft total weight, must specify the number of passtengers,
the fusealge area per passenger (recommended to use 1 m^2 based on research), and the
engine weight
Rate of climb equation taken from John Anderson's Aircraft Performance and Design (eqn 5.85)
"""

class CommericalMissionConstraints(Model):
    """
    class that is general constraints that apply across the mission
    """
    def __init__(self, Nclimb, Ncruise, **kwargs):

        #altitude variables
        dhTakeoff = Variable('dh_{takeoff}', 1500, 'feet', 'Total Altitude Change During Takeoff Profile')
        dhftClimb = VectorVariable(Nclimb, 'dhftClimb', 'feet', 'Change in Altitude Per Climb Segment [feet]')
        hClimb = VectorVariable(Nclimb, 'hClimb', 'm', 'Altitude [meters]')
        hftClimb = VectorVariable(Nclimb, 'hftClimb', 'feet', 'Altitude [feet]')
        hCruise = VectorVariable(Ncruise, 'hCruise', 'm', 'Altitude [meters]')
        hftCruise = VectorVariable(Ncruise, 'hftCruise', 'feet', 'Altitude [feet]')

        #alttidue at the top of climb
        htoc = Variable('h_{toc}', 'feet', 'Altitude at Top of Climb')

        #weight variables
        W_payload = Variable('W_{payload}', 'N', 'Aircraft Payload Weight')
        W_e = Variable('W_{e}', 'N', 'Empty Weight of Aircraft')
        W_engine = Variable('W_{engine}', 'N', 'Weight of a Single Turbofan Engine')
        W_ftotal = Variable('W_{f_{total}}', 'N', 'Total Fuel Weight')
        W_wing = Variable('W_{wing}', 'N', 'Wing Weight')
        W_total = Variable('W_{total}', 'N', 'Total Aircraft Weight')
        W_startClimb = VectorVariable(Nclimb, 'W_{startClimb}', 'N', 'Segment Start Weight')
        W_fuelClimb = VectorVariable(Nclimb, 'W_{fuelClimb}', 'N', 'Segment Fuel Weight')
        W_endClimb = VectorVariable(Nclimb, 'W_{enDClimb}', 'N', 'Segment End Weight')
        W_avgClimb = VectorVariable(Nclimb, 'W_{avgClimb}', 'N', 'Geometric Average of Segment Start and End Weight')
        W_startCruise = VectorVariable(Ncruise, 'W_{startCruise}', 'N', 'Segment Start Weight')
        W_fuelCruise = VectorVariable(Ncruise, 'W_{fuelCruise}', 'N', 'Segment Fuel Weight')
        W_endCruise = VectorVariable(Ncruise, 'W_{endCruise}', 'N', 'Segment End Weight')
        W_avgCruise = VectorVariable(Ncruise, 'W_{avgCruise}', 'N', 'Geometric Average of Segment Start and End Weight')
        W_pax = Variable('W_{pax}', 'N', 'Estimated Average Passenger Weight, Includes Baggage')

        #number of passengers
        n_pax = Variable('n_{pax}', '-', 'Number of Passengers to Carry')

        #aircraft geometry
        S = Variable('S', 'm^2', 'Wing Planform Area')
        AR = Variable('AR', '-', 'Aspect Ratio')
        span = Variable('span', 'm', 'Wing Span')
        span_max = Variable('span_{max}', 'm', 'Max Wing Span')

        #number of engines
        numeng = Variable('numeng', '-', 'Number of Engines')
        
        #physical constants
        g = Variable('g', 9.81, 'm/s^2', 'Gravitational Acceleration')
        gamma = Variable('\gamma', 1.4, '-', 'Air Specific Heat Ratio')
        R = Variable('R', 287, 'J/kg/K', 'Gas Constant for Air')

        #Fuselage area
        A_fuse = Variable('A_{fuse}', 'm^2', 'Estimated Fuselage Area')
        pax_area = Variable('pax_{area}', 'm^2', 'Estimated Fuselage Area per Passenger')
        
        constraints = []
        constraints.extend([
            #convert m to ft
            hftClimb  == hClimb,
            hftCruise  == hCruise,
            
            #constraints on the various weights
            W_payload == n_pax * W_pax,
            
            TCS([W_e + W_payload + W_ftotal + numeng * W_engine + W_wing <= W_total]),
 
            W_startClimb[0] == W_total,
            W_startCruise[0] == W_endClimb[Nclimb-1],
   
            TCS([W_e + W_payload + numeng * W_engine + W_wing <= W_endCruise[Ncruise-1]]),
            TCS([W_ftotal >= sum(W_fuelClimb) + sum(W_fuelCruise)]),

            #wing weight constraint
            #based off of a raymer weight and 737 data from TASOPT output file
            (S/(124.58*units('m^2')))**.65 == W_wing/(105384.1524*units('N')),

            #compute wing span and aspect ratio, subject to a span constraint
            AR == (span**2)/S,
            span <= span_max,

            #estimate based on TASOPT 737 model
            W_e == .75*W_payload,

            #compute fuselage area for drag approximation
            A_fuse == pax_area * n_pax,
            
            #altitude buildup constraints for climb segment 2
            TCS([hftClimb[0] >= dhftClimb[0]]),
            
            hftClimb <= htoc,
            ])

        for i in range(1, Nclimb):
            constraints.extend([
                TCS([W_startClimb[i] == W_endClimb[i-1]]),
                TCS([hftClimb[i] >= hftClimb[i - 1] + dhftClimb[i]]),
                ])

        for i in range(1, Ncruise):
            constraints.extend([
                TCS([W_startCruise[i] == W_endCruise[i-1]]),
                ])

        for i in range(0, Nclimb):
            constraints.extend([
                W_avgClimb[i] == (W_startClimb[i]*W_endClimb[i])**.5,
                TCS([W_startClimb[i] >= W_endClimb[i] + W_fuelClimb[i]]),
                ])

        for i in range(0, Ncruise):
            constraints.extend([
                W_avgCruise[i] == (W_startCruise[i]*W_endCruise[i])**.5,
                TCS([W_startCruise[i] >= W_endCruise[i] + W_fuelCruise[i]]),
                ])

        Model.__init__(self, W_ftotal, constraints, **kwargs)
            
class Climb(Model):
    """
    class to model the climb portion of the flight (from zero alt to cruise alt)
    """
    def __init__(self, Nclimb, Ncruise, **kwargs):
        #aero
        CLClimb = VectorVariable(Nclimb, 'C_{L_{Climb}}', '-', 'Lift Coefficient')
        WLoadClimb = VectorVariable(Nclimb, 'W_{Load_{Climb}}', 'N/m^2', 'Wing Loading')
        WLoadmax = Variable('W_{Load_max}', 'N/m^2', 'Max Wing Loading')
        DClimb = VectorVariable(Nclimb, 'DragClimb', 'N', 'Drag')
        Vstall = Variable('V_{stall}', 'knots', 'Aircraft Stall Speed')
        Cdwc = VectorVariable(Nclimb, 'C_{d_wc}', '-', 'Cd for a NC130 Airfoil at Re=2e7')
        cdp_rc = VectorVariable(Nclimb, 'cdp_rc', '-', 'Hold Variable for Drag Fit')
        Cdfuse = Variable('C_{d_fuse}', '-', 'Fuselage Drag Coefficient')
        K = Variable('K', '-', 'K for Parametric Drag Model')
        e = Variable('e', '-', 'Oswald Span Efficiency Factor')
        AR = Variable('AR', '-', 'Aspect Ratio')
        
        #atmosphere
        aClimb = VectorVariable(Nclimb, 'aClimb', 'm/s', 'Speed of Sound')
        rhoClimb = VectorVariable(Nclimb, '\\rhoClimb', 'kg/m^3', 'Air Density')
        pClimb = VectorVariable(Nclimb, 'pClimb', 'kPa', 'Pressure')
        TClimb = VectorVariable(Nclimb, 'TClimb', 'K', 'Air Temperature')
 
        #number of engines
        numeng = Variable('numeng', '-', 'Number of Engines')

        #aircraft geometry
        S = Variable('S', 'm^2', 'Wing Planform Area')

        #physical constants
        g = Variable('g', 9.81, 'm/s^2', 'Gravitational Acceleration')
        gamma = Variable('\gamma', 1.4, '-', 'Air Specific Heat Ratio')
        R = Variable('R', 287, 'J/kg/K', 'Gas Constant for Air')

        #thrust
        thrustClimb = Variable('thrustClimb', 'N', 'Engine Thrust')

        #excess power during climb
        excesspClimb = VectorVariable(Nclimb, 'Excess Power Climb', 'W', 'Excess Power During Climb')

        #climb rate and angle
        RCClimb = VectorVariable(Nclimb, 'RCClimb', 'feet/min', 'Rate of Climb/Decent')
        thetaClimb = VectorVariable(Nclimb, '\\thetaClimb', '-', 'Aircraft Climb Angle')

        #time
        tmiNclimb = VectorVariable(Nclimb, 'tmiNclimb', 'min', 'Flight Time in Minutes')
        thoursClimb = VectorVariable(Nclimb, 'thrClimb', 'hour', 'Flight Time in Hours')

        #range
        RngClimb = VectorVariable(Nclimb, 'RngClimb', 'miles', 'Segment Range During Climb2')

        #velocitites and mach numbers
        VClimb = VectorVariable(Nclimb, 'VClimb', 'knots', 'Aircraft Flight Speed')
        MClimb = VectorVariable(Nclimb, 'MClimb', '-', 'Aircraft Mach Number')

        #altitude
        dhftClimb = VectorVariable(Nclimb, 'dhftClimb', 'feet', 'Change in Altitude Per Climb Segment [feet]') 
        htoc = Variable('h_{toc}', 'feet', 'Altitude at Top of Climb')

        W_fuelClimb = VectorVariable(Nclimb, 'W_{fuelClimb}', 'N', 'Segment Fuel Weight')
        W_avgClimb = VectorVariable(Nclimb, 'W_{avgClimb}', 'N', 'Geometric Average of Segment Start and End Weight')
        
        TSFCcl = VectorVariable(Nclimb, 'TSFC_{cl}', '1/hr', 'Thrust Specific Fuel Consumption During Climb2')
        thrustcl = VectorVariable(Nclimb, 'thrust_{cl}', 'N', 'Thrust During Climb Segment #2')

        thrustcr = VectorVariable(Ncruise, 'thrust_{cr}', 'N', 'Thrust During Cruise Segment #2')

        #Fuselage area
        A_fuse = Variable('A_{fuse}', 'm^2', 'Estimated Fuselage Area')

        #non-GPkit variables
        #climb 2 lsit
        icl = map(int, np.linspace(0, Nclimb - 1, Nclimb))

        constraints = []
        
        constraints.extend([            
            #set the velocity limits
            #needs to be replaced by an actual Vne and a mach number
            VClimb[icl] >= Vstall,

            VClimb == MClimb * aClimb,

            #constraint on drag and thrust
            numeng*thrustcl[icl] >= DClimb[icl] + W_avgClimb[icl] * thetaClimb[icl],
            
            #climb rate constraints
            TCS([excesspClimb[icl]+VClimb[icl]*DClimb[icl] <= VClimb[icl]*numeng*thrustcl[icl]]),
            
            TCS([DClimb[icl] >= (.5*S*rhoClimb[icl]*VClimb[icl]**2)*(Cdwc[icl] + K*CLClimb[icl]**2) + Cdfuse * (.5 * A_fuse * rhoClimb[icl] * VClimb[icl]**2)]),
            
            K == (pi * e * AR)**-1,
            
            cdp_rc[icl] >= (1.02458748e10 * CLClimb[icl]**15.587947404823325 * MClimb[icl]**156.86410659495155 +
                2.85612227e-13 * CLClimb[icl]**1.2774976672501526 * MClimb[icl]**6.2534328002723703 +
                2.08095341e-14 * CLClimb[icl]**0.8825277088649582 * MClimb[icl]**0.0273667615730107 +
                1.94411925e+06 * CLClimb[icl]**5.6547413360261691 * MClimb[icl]**146.51920742858428),

            Cdwc[icl]**6.5 >= cdp_rc[icl],
            
            RCClimb[icl] == excesspClimb[icl]/W_avgClimb[icl],
            RCClimb[icl] >= 500*units('ft/min'),
            
            #make the small angle approximation and compute theta
            thetaClimb[icl]*VClimb[icl]  == RCClimb[icl],
           
            dhftClimb[icl]  == tmiNclimb[icl] * RCClimb[icl],

            tmiNclimb == thoursClimb,
            #compute the distance traveled for each segment
            #takes into account two terms of a cosine expansion
            RngClimb[icl] == thoursClimb[icl]*VClimb[icl],

            W_avgClimb[icl] == .5*CLClimb[icl]*S*rhoClimb[icl]*VClimb[icl]**2,      
            WLoadClimb[icl] == .5*CLClimb[icl]*S*rhoClimb[icl]*VClimb[icl]**2/S,
            
            #compute fuel burn from TSFC
            W_fuelClimb[icl]  == numeng*TSFCcl[icl] * thoursClimb[icl] * thrustcl[icl],

            #compute the dh
            dhftClimb[icl] == htoc/Nclimb,

            #constrain the max wing loading
            WLoadClimb <= WLoadmax,

            thrustcl <= 2 * thrustcr[0],

            TSFCcl[0] == .7*units('1/hr'),
            TSFCcl[1] == .6*units('1/hr'),
            ])

        for i in range(0, Nclimb):
            constraints.extend([
                #speed of sound
                aClimb[i]  == (gamma * R * TClimb[i])**.5,
                ])
        
        Model.__init__(self, None, constraints, **kwargs)
        
class Cruise(Model):
    """
    class to model the second cruise portion of a flight (if the flight is
    long enough to mandate two cruise portions)
    Model is based off of a discretized Breguet Range equation
    """
    def __init__(self, Nclimb, Ncruise, **kwargs):
        """
    class to model the second cruise portion of a flight (if the flight is
    long enough to mandate two cruise portions)
    Model is based off of a discretized Breguet Range equation
    """
        #aero
        CLCruise = VectorVariable(Ncruise, 'C_{L_{Cruise}}', '-', 'Lift Coefficient')
        WLoadCruise = VectorVariable(Ncruise, 'W_{Load_{Cruise}}', 'N/m^2', 'Wing Loading')
        WLoadmax = Variable('W_{Load_max}', 'N/m^2', 'Max Wing Loading')
        DCruise = VectorVariable(Ncruise, 'DragCruise', 'N', 'Drag')
        Vstall = Variable('V_{stall}', 'knots', 'Aircraft Stall Speed')
        Cdwcr = VectorVariable(Ncruise, 'C_{d_wcr}', '-', 'Cd for a NC130 Airfoil at Re=2e7')
        cdp_rcr = VectorVariable(Ncruise, 'cdp_rcr', '-', 'Hold Variable for Drag Fit')
        Cdfuse = Variable('C_{d_fuse}', '-', 'Fuselage Drag Coefficient')
        K = Variable('K', '-', 'K for Parametric Drag Model')
        e = Variable('e', '-', 'Oswald Span Efficiency Factor')
        AR = Variable('AR', '-', 'Aspect Ratio')

        #atmosphere
        aCruise = VectorVariable(Ncruise, 'aCruise', 'm/s', 'Speed of Sound')
        rhoCruise = VectorVariable(Ncruise, '\\rhoCruise', 'kg/m^3', 'Air Density')
        pCruise = VectorVariable(Ncruise, 'pCruise', 'kPa', 'Pressure')
        TCruise = VectorVariable(Ncruise, 'TCruise', 'K', 'Air Temperature')

        #aircraft geometry
        S = Variable('S', 'm^2', 'Wing Planform Area')

        #number of engines
        numeng = Variable('numeng', '-', 'Number of Engines')

        #time
        tmiNcruise = VectorVariable(Ncruise, 'tmiNcruise', 'min', 'Flight Time in Minutes')
        thoursCruise = VectorVariable(Ncruise, 'thrCruise', 'hour', 'Flight Time in Hours')

        #velocitites and mach numbers
        VCruise = VectorVariable(Ncruise, 'VCruise', 'knots', 'Aircraft Flight Speed')
        MCruise = VectorVariable(Ncruise, 'MCruise', '-', 'Aircraft Mach Number')

        #physical constants
        g = Variable('g', 9.81, 'm/s^2', 'Gravitational Acceleration')
        gamma = Variable('\gamma', 1.4, '-', 'Air Specific Heat Ratio')
        R = Variable('R', 287, 'J/kg/K', 'Gas Constant for Air')
        
        TSFCcr = VectorVariable(Ncruise, 'TSFC_{cr}', '1/hr', 'Thrust Specific Fuel Consumption During Cruise2')
        D = Variable('D', 'N', 'Drag for Cruise')

        thrustcr = VectorVariable(Ncruise, 'thrust_{cr}', 'N', 'Thrust During Cruise')

        constraints = []
        
        #defined here for linking purposes
        T0 = Variable('T_0', 'K', 'Free Stream Stagnation Temperature')
        P0 = Variable('P_0', 'kPa', 'Free Stream Static Pressure')

        #parameter to make breguent range eqn gp compatible
        z_bre = VectorVariable(Ncruise, 'z_{bre}', '-', 'Breguet Parameter')

        #range variables
        ReqRng = Variable('ReqRng', 'miles', 'Required Mission Range')
        RngCruise = VectorVariable(Ncruise, 'RngCruise', 'miles', 'Segment Range During Cruise2')

        #altitude
        hCruise = VectorVariable(Ncruise, 'hCruise', 'm', 'Altitude [meters]')
        hftCruise = VectorVariable(Ncruise, 'hftCruise', 'feet', 'Altitude [feet]')

        #alttidue at the top of climb
        htoc = Variable('h_{toc}', 'feet', 'Altitude at Top of Climb')
 
        W_fuelCruise = VectorVariable(Ncruise, 'W_{fuelCruise}', 'N', 'Segment Fuel Weight')
        W_avgCruise = VectorVariable(Ncruise, 'W_{avgCruise}', 'N', 'Geometric Average of Segment Start and End Weight')
        W_avgClimb = VectorVariable(Nclimb, 'W_{avgClimb}', 'N', 'Geometric Average of Segment Start and End Weight')
        W_endCruise = VectorVariable(Ncruise, 'W_{endCruise}', 'N', 'Segment End Weight')
            
        #Fuselage area
        A_fuse = Variable('A_{fuse}', 'm^2', 'Estimated Fuselage Area')

        #non-GPkit variables
        #cruise 2 lsit
        icr = map(int, np.linspace(0, Ncruise - 1, Ncruise))
            
        constraints.extend([
##            MCruise[icr] == 0.8,

            MCruise * aCruise == VCruise,
            
            TCS([DCruise[icr] >= (.5 * S * rhoCruise[icr] * VCruise[icr]**2) *
                 (Cdwcr[icr] + K * (W_avgCruise[icr] / (.5 * S * rhoCruise[icr]* VCruise[icr]**2))**2)
                 + Cdfuse * (.5 * A_fuse * rhoCruise[icr] * VCruise[icr]**2)]),

            K == (pi * e * AR)**-1,

            cdp_rcr[icr] >= (1.02458748e10 * CLCruise[icr]**15.587947404823325 * MCruise[icr]**156.86410659495155 +
                2.85612227e-13 * CLCruise[icr]**1.2774976672501526 * MCruise[icr]**6.2534328002723703 +
                2.08095341e-14 * CLCruise[icr]**0.8825277088649582 * MCruise[icr]**0.0273667615730107 +
                1.94411925e+06 * CLCruise[icr]**5.6547413360261691 * MCruise[icr]**146.51920742858428),

            Cdwcr[icr]**6.5 >= cdp_rcr[icr],
            
            DCruise[icr] == numeng * thrustcr[icr],

            W_avgCruise[icr] == .5 * CLCruise[icr] * S * rhoCruise[icr] * VCruise[icr]**2,
            WLoadCruise[icr] == .5 * CLCruise[icr] * S * rhoCruise[icr] * VCruise[icr]**2 / S,
            
            #taylor series expansion to get the weight term
            TCS([W_fuelCruise[icr]/W_endCruise[icr] >= te_exp_minus1(z_bre[icr], nterm=3)]),

            #breguet range eqn
            TCS([z_bre[icr] >= (numeng * TSFCcr[icr] * thoursCruise[icr] * DCruise[icr]) / W_avgCruise[icr]]),

            #time
            thoursCruise[icr]*VCruise[icr]  == RngCruise[icr],
            tmiNcruise == thoursCruise,

            #constrain the max wing loading
            WLoadCruise <= WLoadmax,

            #altitude constraint
            hCruise == 40000*units('ft'),

            TSFCcr[0] == .5*units('1/hr'),
            TSFCcr[1] == .5*units('1/hr'),
            ])
        
        #constraint on the aircraft meeting the required range
        for i in range(min(icr), max(icr)+1):
            constraints.extend([
                TCS([RngCruise[i] == ReqRng/(Ncruise)])
                ])

        for i in range(0, Ncruise):
            constraints.extend([
                #speed of sound
                aCruise[i]  == (gamma * R * TCruise[i])**.5,
                ])
            
        Model.__init__(self, None, constraints, **kwargs)

#-------------------------------------
#build the linked model

class CommercialAircraft(Model):
    """
    class to link all models needed to simulate a commercial flight
    """
    def __init__(self, **kwargs):
        #defining the number of segments
        Nclimb = 2
        Ncruise = 2

        #make the segment range
        Nseg = Nclimb + Ncruise

        #define all the submodels
        cmc = CommericalMissionConstraints(Nclimb, Ncruise)
        climb = Climb(Nclimb, Ncruise)
        cruise = Cruise(Nclimb, Ncruise)

        atmvec = []
        
        for i in range(Nseg):
            atmvec.append(Atmosphere())
            
        substitutions = {      
            'V_{stall}': 120,
            'ReqRng': 2000,
            'K': 0.05,
##            'h_{toc}': 40000, #('sweep', np.linspace(20000,40000,4)),
            'numeng': 2,
            'W_{Load_max}': 6664,
            'W_{engine}': 1000,
            'W_{pax}': 91 * 9.81,
            'n_{pax}': 150,
            'pax_{area}': 1,
            'C_{d_fuse}': .005, #assumes turbulent flow, from wikipedia
            'e': .9,
            'span_{max}': 35,

            #atm subs
            "p_{sl}": 101325,
            "T_{sl}": 288.15,
            "L_{atm}": .0065,
            "M_{atm}":.0289644,
            "R_{atm}": 8.31447
            }
        
        submodels = [cmc, climb, cruise]

        for i in range(len(atmvec)):
            submodels.extend(atmvec[i])

        constraints = ConstraintSet([submodels])

        subs= {}

        for i in range(Nclimb):
            subs.update({
                climb["\\rhoClimb"][i]: atmvec[i]["\\rho"], climb["TClimb"][i]: atmvec[i]["T_{atm}"], cmc['hClimb'][i]: atmvec[i]["h"]
                })

        for i in range(Ncruise):
            subs.update({
                cruise["\\rhoCruise"][i]: atmvec[i + Nclimb]["\\rho"], cruise["TCruise"][i]:atmvec[i + Nclimb]["T_{atm}"],
                cruise['hCruise'][i]: atmvec[i + Nclimb]["h"]
                })

        constraints.subinplace(subs)
        
        lc = LinkedConstraintSet(constraints, exclude={"T_{atm}", "P_{atm}", '\\rho', "h"})

        Model.__init__(self, cmc.cost, lc, substitutions, **kwargs)

    def bound_all_variables(self, model, eps=1e-30, lower=None, upper=None):
        "Returns model with additional constraints bounding all free variables"
        lb = lower if lower else eps
        ub = upper if upper else 1/eps
        constraints = []
        freevks = tuple(vk for vk in model.varkeys if "value" not in vk.descr)
        for varkey in freevks:
            units = varkey.descr.get("units", 1)
            constraints.append([ub*units >= Variable(**varkey.descr),
                                Variable(**varkey.descr) >= lb*units])
        m = Model(model.cost, [constraints, model], model.substitutions)
        m.bound_all = {"lb": lb, "ub": ub, "varkeys": freevks}
        return m


    # pylint: disable=too-many-locals
    def determine_unbounded_variables(self, model, solver=None, verbosity=0,
                                      eps=1e-30, lower=None, upper=None, **kwargs):
        "Returns labeled dictionary of unbounded variables."
        m = self.bound_all_variables(model, eps, lower, upper)
        sol = m.localsolve(solver, verbosity, **kwargs)
        solhold = sol
        lam = sol["sensitivities"]["la"][1:]
        out = defaultdict(list)
        for i, varkey in enumerate(m.bound_all["varkeys"]):
            lam_gt, lam_lt = lam[2*i], lam[2*i+1]
            if abs(lam_gt) >= 1e-7:  # arbitrary threshold
                out["sensitive to upper bound"].append(varkey)
            if abs(lam_lt) >= 1e-7:  # arbitrary threshold
                out["sensitive to lower bound"].append(varkey)
            value = mag(sol["variables"][varkey])
            distance_below = np.log(value/m.bound_all["lb"])
            distance_above = np.log(m.bound_all["ub"]/value)
            if distance_below <= 3:  # arbitrary threshold
                out["value near lower bound"].append(varkey)
            elif distance_above <= 3:  # arbitrary threshold
                out["value near upper bound"].append(varkey)
        return out, solhold

    
if __name__ == '__main__':
    m = CommercialAircraft()
##    sol = m.localsolve(solver="mosek", verbosity = 4, iteration_limit=100, skipsweepfailures=True)
    
    sol, solhold = m.determine_unbounded_variables(m, solver="mosek",verbosity=4, iteration_limit=100)
    
