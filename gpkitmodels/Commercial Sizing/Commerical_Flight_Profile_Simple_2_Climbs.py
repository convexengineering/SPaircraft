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
Minimizes the aircraft total fuel weight. Rate of climb equation taken from John
Anderson's Aircraft Performance and Design (eqn 5.85).

Inputs
-----

- Number of passtengers,
- Fusealge area per passenger (recommended to use 1 m^2 based on research)
- Engine weight
"""

class CommericalMissionConstraints(Model):
    """
    General constraints that apply across the mission
    """
    def __init__(self, Nclimb1, Nclimb2, Ncruise2, signomial=0, **kwargs):
        #variable local to this model
        #10,000 foot altitude
        alt10k = Variable('alt10k', 10000, 'feet', '10,000 feet')

        #altitude variables
        dhTakeoff = Variable('dh_{takeoff}', 1500, 'feet', 'Total Altitude Change During Takeoff Profile')
        dhClimb1 = Variable('dh_{climb1}', 8500, 'feet', 'Total Altitude Change Required in Climb 1')
        dhClimb2 = Variable('dh_{climb2}', 'feet', 'Total Altitude Change Required in Climb 2')
        dhftClimb1 = VectorVariable(Nclimb1, 'dhftClimb1', 'feet', 'Change in Altitude Per Climb Segment [feet]')
        dhftClimb2 = VectorVariable(Nclimb2, 'dhftClimb2', 'feet', 'Change in Altitude Per Climb Segment [feet]')
        hClimb1 = VectorVariable(Nclimb1, 'hClimb1', 'm', 'Altitude [meters]')
        hftClimb1 = VectorVariable(Nclimb1, 'hftClimb1', 'feet', 'Altitude [feet]')
        hClimb2 = VectorVariable(Nclimb2, 'hClimb2', 'm', 'Altitude [meters]')
        hftClimb2 = VectorVariable(Nclimb2, 'hftClimb2', 'feet', 'Altitude [feet]')
        hCruise2 = Variable('hCruise2', 'm', 'Altitude [meters]')
        hftCruise2 = Variable('hftCruise2', 'feet', 'Altitude [feet]')

        #weight variables
        W_payload = Variable('W_{payload}', 'N', 'Aircraft Payload Weight')
        W_e = Variable('W_{e}', 'N', 'Empty Weight of Aircraft')
        W_engine = Variable('W_{engine}', 'N', 'Weight of a Single Turbofan Engine')
        W_ftotal = Variable('W_{f_{total}}', 'N', 'Total Fuel Weight')
        W_wing = Variable('W_{wing}', 'N', 'Wing Weight')
        W_total = Variable('W_{total}', 'N', 'Total Aircraft Weight')
        W_startClimb1 = VectorVariable(Nclimb1, 'W_{startClimb1}', 'N', 'Segment Start Weight')
        W_fuelClimb1 = VectorVariable(Nclimb1, 'W_{fuelClimb1}', 'N', 'Segment Fuel Weight')
        W_endClimb1 = VectorVariable(Nclimb1, 'W_{endClimb1}', 'N', 'Segment End Weight')
        W_avgClimb1 = VectorVariable(Nclimb1, 'W_{avgClimb1}', 'N', 'Geometric Average of Segment Start and End Weight')
        W_startClimb2 = VectorVariable(Nclimb2, 'W_{startClimb2}', 'N', 'Segment Start Weight')
        W_fuelClimb2 = VectorVariable(Nclimb2, 'W_{fuelClimb2}', 'N', 'Segment Fuel Weight')
        W_endClimb2 = VectorVariable(Nclimb2, 'W_{endClimb2}', 'N', 'Segment End Weight')
        W_avgClimb2 = VectorVariable(Nclimb2, 'W_{avgClimb2}', 'N', 'Geometric Average of Segment Start and End Weight')
        W_startCruise2 = VectorVariable(Ncruise2, 'W_{startCruise2}', 'N', 'Segment Start Weight')
        W_fuelCruise2 = VectorVariable(Ncruise2, 'W_{fuelCruise2}', 'N', 'Segment Fuel Weight')
        W_endCruise2 = VectorVariable(Ncruise2, 'W_{endCruise2}', 'N', 'Segment End Weight')
        W_avgCruise2 = VectorVariable(Ncruise2, 'W_{avgCruise2}', 'N', 'Geometric Average of Segment Start and End Weight')
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

        #range variables
        ReqRngCruise = Variable('ReqRngCruise', 'miles', 'Required Cruise Range')
        ReqRng = Variable('ReqRng', 'miles', 'Required Mission Range')
        RngClimb1 = VectorVariable(Nclimb1, 'RngClimb1', 'miles', 'Segment Range During Climb1')
        RngClimb2 = VectorVariable(Nclimb2, 'RngClimb2', 'miles', 'Segment Range During Climb2')

        #physical constants
        g = Variable('g', 9.81, 'm/s^2', 'Gravitational Acceleration')
        gamma = Variable('\gamma', 1.4, '-', 'Air Specific Heat Ratio')
        R = Variable('R', 287, 'J/kg/K', 'Gas Constant for Air')

        #Fuselage area
        A_fuse = Variable('A_{fuse}', 'm^2', 'Estimated Fuselage Area')
        pax_area = Variable('pax_{area}', 'm^2', 'Estimated Fuselage Area per Passenger')

        test= VectorVariable(2, 'test', 'm', 'test')
        
        constraints = []
        constraints.extend([
            #convert m to ft
            hftClimb1  == hClimb1,
            hftClimb2  == hClimb2,
            hftCruise2  == hCruise2,
            
            
            #constraints on the various weights
            W_payload == n_pax * W_pax,
            
            TCS([W_e + W_payload + W_ftotal + numeng * W_engine + W_wing <= W_total]),
 
            W_startClimb1[0]  == W_total,
   
            TCS([W_e + W_payload + numeng * W_engine + W_wing <= W_endCruise2[Ncruise2-1]]),
            TCS([W_ftotal >= sum(W_fuelClimb1) + sum(W_fuelClimb2) + sum(W_fuelCruise2)]),

            W_startClimb2[0] == W_endClimb1[Nclimb1-1],
            W_startCruise2[0] == W_endClimb2[Nclimb2-1],

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

            #altitude buildup constraints for climb segment 1
            TCS([hftClimb1[0] >= 1500 * units('ft') + dhftClimb1[0]]),
            hftClimb1 <= alt10k,
            
            #altitude buildup constraints for climb segment 2
            TCS([hftClimb2[0] >= hftClimb1[Nclimb1-1] + dhftClimb2[0]]),
            hftClimb2 <= hftCruise2,
            ])
                        
        for i in range(0, Nclimb1):
            constraints.extend([
                #constrain the geometric weight average
                W_avgClimb1[i] == (W_startClimb1[i]*W_endClimb1[i])**.5,
                TCS([W_startClimb1[i] >= W_endClimb1[i] + W_fuelClimb1[i]]),                
                ])

        for i in range(1, Nclimb1):
            constraints.extend([
                TCS([W_startClimb1[i] == W_endClimb1[i-1]]),
                TCS([hftClimb1[i] >= hftClimb1[i-1] + dhftClimb1[i]]),
                ])

        for i in range(1, Nclimb2):
            constraints.extend([
                TCS([W_startClimb2[i] == W_endClimb2[i-1]]),
                TCS([hftClimb2[i] >= hftClimb2[i-1] + dhftClimb2[i]]),
                ])

        for i in range(1, Ncruise2):
            constraints.extend([
                TCS([W_startCruise2[i] == W_endCruise2[i-1]]),
                ])

        for i in range(0, Nclimb2):
            constraints.extend([
                W_avgClimb2[i] == (W_startClimb2[i]*W_endClimb2[i])**.5,
                TCS([W_startClimb2[i] >= W_endClimb2[i] + W_fuelClimb2[i]]),
                ])

        for i in range(0, Ncruise2):
            constraints.extend([
                W_avgCruise2[i] == (W_startCruise2[i]*W_endCruise2[i])**.5,
                TCS([W_startCruise2[i] >= W_endCruise2[i] + W_fuelCruise2[i]]),
                ])

        with gpkit.SignomialsEnabled():
            if signomial == True:
                constraints.extend([
                    #range constraints
                    TCS([sum(RngClimb1) + sum(RngClimb2) + ReqRngCruise >= ReqRng]),
                    ])
                
            if signomial == False:
                 constraints.extend([
                     TCS([ReqRngCruise >= ReqRng]),
                    ])
                 
            constraints.extend([
                #compute hdclimb2
                hftCruise2 <= dhClimb2 + alt10k,
                ])

        Model.__init__(self, W_ftotal, constraints, **kwargs)
        
#---------------------------------------
#takeoff

class Climb1(Model):
    """
    class to model the climb portion of a flight, applies to all climbs below
    10,000'
    """
    def __init__(self, Nclimb1, Ncruise2, **kwargs):
        #set the speed limit under 10,000'
        speedlimit = Variable('speedlimit', 'kts', 'Speed Limit Under 10,000 ft')

        #atmosphere
        aClimb1 = VectorVariable(Nclimb1, 'aClimb1', 'm/s', 'Speed of Sound')
        rhoClimb1 = VectorVariable(Nclimb1, '\\rhoClimb1', 'kg/m^3', 'Air Density')
        pClimb1 = VectorVariable(Nclimb1, 'pClimb1', 'kPa', 'Pressure')
        muClimb1 = VectorVariable(Nclimb1, '\muClimb1', 'kg/m/s', 'Air Kinematic Viscosity')
        TClimb1 = VectorVariable(Nclimb1, 'TClimb1', 'K', 'Air Temperature')
        
        #aero
        CLClimb1 = VectorVariable(Nclimb1, 'C_{L_{Climb1}}', '-', 'Lift Coefficient')
        WLoadClimb1 = VectorVariable(Nclimb1, 'W_{Load_{Climb1}}', 'N/m^2', 'Wing Loading')
        WLoadmax = Variable('W_{Load_max}', 'N/m^2', 'Max Wing Loading')
        DClimb1 = VectorVariable(Nclimb1, 'DragClimb1', 'N', 'Drag')
        Vstall = Variable('V_{stall}', 'knots', 'Aircraft Stall Speed')
        Cdwc1 = VectorVariable(Nclimb1, 'C_{d_wc1}', '-', 'Cd for a NC130 Airfoil at Re=2e7')
        Cdfuse = Variable('C_{d_fuse}', '-', 'Fuselage Drag Coefficient')
        K = Variable('K', '-', 'K for Parametric Drag Model')
        e = Variable('e', '-', 'Oswald Span Efficiency Factor')
        AR = Variable('AR', '-', 'Aspect Ratio')

        #aircraft geometry
        S = Variable('S', 'm^2', 'Wing Planform Area')

        #number of engines
        numeng = Variable('numeng', '-', 'Number of Engines')

        #physical constants
        g = Variable('g', 9.81, 'm/s^2', 'Gravitational Acceleration')
        gamma = Variable('\gamma', 1.4, '-', 'Air Specific Heat Ratio')
        R = Variable('R', 287, 'J/kg/K', 'Gas Constant for Air')

        #excess power during climb
        excessPclimb1 = VectorVariable(Nclimb1, 'Excess Power Climb1', 'W', 'Excess Power During Climb')

        #climb rate and angle
        RCClimb1 = VectorVariable(Nclimb1, 'RCClimb1', 'feet/min', 'Rate of Climb/Decent')
        thetaClimb1 = VectorVariable(Nclimb1, '\\thetaClimb1', '-', 'Aircraft Climb Angle')

        #time
        tminClimb1 = VectorVariable(Nclimb1, 'tminClimb1', 'min', 'Flight Time in Minutes')
        thoursClimb1 = VectorVariable(Nclimb1, 'thrClimb1', 'hour', 'Flight Time in Hours')

        #range
        RngClimb1 = VectorVariable(Nclimb1, 'RngClimb1', 'miles', 'Segment Range During Climb1')

        #velocitites and mach numbers
        VClimb1 = VectorVariable(Nclimb1, 'VClimb1', 'knots', 'Aircraft Flight Speed')
        MClimb1 = VectorVariable(Nclimb1, 'MClimb1', '-', 'Aircraft Mach Number')

        #altitude
        dhftClimb1 = VectorVariable(Nclimb1, 'dhftClimb1', 'feet', 'Change in Altitude Per Climb Segment [feet]') 
        dhClimb1 = Variable('dh_{climb1}', 8500, 'feet', 'Total Altitude Change Required in Climb 1')

        #Weights
        W_avgClimb1 = VectorVariable(Nclimb1, 'W_{avgClimb1}', 'N', 'Geometric Average of Segment Start and End Weight')
        W_fuelClimb1 = VectorVariable(Nclimb1, 'W_{fuelClimb1}', 'N', 'Segment Fuel Weight')

        #Engine
        TSFCc1 = VectorVariable(Nclimb1, 'TSFC_{c1}', '1/hr', 'Thrust Specific Fuel Consumption During Climb1')
        thrustc1 = VectorVariable(Nclimb1, 'thrust_{c1}', 'N', 'Thrust During Climb Segment #1')

        #Cruise engine thrust
        thrustcr2 = VectorVariable(Ncruise2, 'thrust_{cr2}', 'N', 'Thrust During Cruise Segment #2')

        #Fuselage area
        A_fuse = Variable('A_{fuse}', 'm^2', 'Estimated Fuselage Area')
        
        #non-GPkit variables
        #cruise 2 lsit
        icl1 = map(int, np.linspace(0, Nclimb1 - 1, Nclimb1))
        
        constraints = []
            
        constraints.extend([            
            #set the velocity limits
            VClimb1[icl1] <= speedlimit,
            VClimb1[icl1] >= Vstall,

            MClimb1 * aClimb1 == VClimb1,

            #constraint on drag and thrust
            numeng*thrustc1[icl1] >= DClimb1[icl1] + W_avgClimb1[icl1]*thetaClimb1[icl1],

            #climb rate constraints
            TCS([excessPclimb1[icl1]+VClimb1[icl1]*DClimb1[icl1] <= VClimb1[icl1]*numeng*thrustc1[icl1]]),
            
            TCS([DClimb1[icl1] >= (.5*S*rhoClimb1[icl1]*VClimb1[icl1]**2)*(Cdwc1[icl1] + K*CLClimb1[icl1]**2) + Cdfuse * (.5 * A_fuse * rhoClimb1[icl1] * VClimb1[icl1]**2)]),
            
            K == (pi * e * AR)**-1,
            
            Cdwc1[icl1]**6.5 >= (1.02458748e10 * CLClimb1[icl1]**15.587947404823325 * MClimb1[icl1]**156.86410659495155 +
                2.85612227e-13 * CLClimb1[icl1]**1.2774976672501526 * MClimb1[icl1]**6.2534328002723703 +
                2.08095341e-14 * CLClimb1[icl1]**0.8825277088649582 * MClimb1[icl1]**0.0273667615730107 +
                1.94411925e+06 * CLClimb1[icl1]**5.6547413360261691 * MClimb1[icl1]**146.51920742858428),
            
            RCClimb1[icl1] == excessPclimb1[icl1]/W_avgClimb1[icl1],
            RCClimb1[icl1] >= 500*units('ft/min'),
            
            #make the small angle approximation and compute theta
            thetaClimb1[icl1]*VClimb1[icl1]  == RCClimb1[icl1],
           
            dhftClimb1[icl1]  == tminClimb1[icl1] * RCClimb1[icl1],

            tminClimb1 == thoursClimb1,
            #compute the distance traveled for each segment
            #takes into account two terms of a cosine expansion
            RngClimb1[icl1] == thoursClimb1[icl1]*VClimb1[icl1],

            W_avgClimb1[icl1] == .5*CLClimb1[icl1]*S*rhoClimb1[icl1]*VClimb1[icl1]**2,
            WLoadClimb1[icl1] == .5*CLClimb1[icl1]*S*rhoClimb1[icl1]*VClimb1[icl1]**2/S,
            
            #compute fuel burn from TSFC
            W_fuelClimb1[icl1]  == numeng*TSFCc1[icl1] * thoursClimb1[icl1] * thrustc1[icl1],

            #compute the dh required for each climb 1 segment
            dhftClimb1[icl1] == dhClimb1/Nclimb1,

            #constrain the max wing loading
            WLoadClimb1 <= WLoadmax,

            #make a realistic constraint on thrust
            thrustc1 <= 2 * thrustcr2[0],
            ])

        for i in range(0, Nclimb1):
            constraints.extend([
                #speed of sound
                aClimb1[i]  == (gamma * R * TClimb1[i])**.5,
                TSFCc1[i] == .7*units('1/hr'),
                ])
        
        Model.__init__(self, None, constraints, **kwargs)
            
class Climb2(Model):
    """
    class to model the climb portion above 10,000'
    """
    def __init__(self, Nclimb2, Ncruise2, **kwargs):
        #aero
        CLClimb2 = VectorVariable(Nclimb2, 'C_{L_{Climb2}}', '-', 'Lift Coefficient')
        WLoadClimb2 = VectorVariable(Nclimb2, 'W_{Load_{Climb2}}', 'N/m^2', 'Wing Loading')
        WLoadmax = Variable('W_{Load_max}', 'N/m^2', 'Max Wing Loading')
        DClimb2 = VectorVariable(Nclimb2, 'DragClimb2', 'N', 'Drag')
        Vstall = Variable('V_{stall}', 'knots', 'Aircraft Stall Speed')
        Cdwc2 = VectorVariable(Nclimb2, 'C_{d_wc2}', '-', 'Cd for a NC130 Airfoil at Re=2e7')
        Cdfuse = Variable('C_{d_fuse}', '-', 'Fuselage Drag Coefficient')
        K = Variable('K', '-', 'K for Parametric Drag Model')
        e = Variable('e', '-', 'Oswald Span Efficiency Factor')
        AR = Variable('AR', '-', 'Aspect Ratio')
        
        #atmosphere
        aClimb2 = VectorVariable(Nclimb2, 'aClimb2', 'm/s', 'Speed of Sound')
        rhoClimb2 = VectorVariable(Nclimb2, '\\rhoClimb2', 'kg/m^3', 'Air Density')
        pClimb2 = VectorVariable(Nclimb2, 'pClimb2', 'kPa', 'Pressure')
        muClimb2 = VectorVariable(Nclimb2, '\muClimb2', 'kg/m/s', 'Air Kinematic Viscosity')
        TClimb2 = VectorVariable(Nclimb2, 'TClimb2', 'K', 'Air Temperature')
 
        #number of engines
        numeng = Variable('numeng', '-', 'Number of Engines')

        #aircraft geometry
        S = Variable('S', 'm^2', 'Wing Planform Area')

        #physical constants
        g = Variable('g', 9.81, 'm/s^2', 'Gravitational Acceleration')
        gamma = Variable('\gamma', 1.4, '-', 'Air Specific Heat Ratio')
        R = Variable('R', 287, 'J/kg/K', 'Gas Constant for Air')

        #thrust
        thrustClimb2 = Variable('thrustClimb2', 'N', 'Engine Thrust')

        #excess power during climb
        excessPclimb2 = VectorVariable(Nclimb2, 'Excess Power Climb2', 'W', 'Excess Power During Climb')

        #climb rate and angle
        RCClimb2 = VectorVariable(Nclimb2, 'RCClimb2', 'feet/min', 'Rate of Climb/Decent')
        thetaClimb2 = VectorVariable(Nclimb2, '\\thetaClimb2', '-', 'Aircraft Climb Angle')

        #time
        tminClimb2 = VectorVariable(Nclimb2, 'tminClimb2', 'min', 'Flight Time in Minutes')
        thoursClimb2 = VectorVariable(Nclimb2, 'thrClimb2', 'hour', 'Flight Time in Hours')

        #range
        RngClimb2 = VectorVariable(Nclimb2, 'RngClimb2', 'miles', 'Segment Range During Climb2')

        #velocitites and mach numbers
        VClimb2 = VectorVariable(Nclimb2, 'VClimb2', 'knots', 'Aircraft Flight Speed')
        MClimb2 = VectorVariable(Nclimb2, 'MClimb2', '-', 'Aircraft Mach Number')

        #altitude
        dhftClimb2 = VectorVariable(Nclimb2, 'dhftClimb2', 'feet', 'Change in Altitude Per Climb Segment [feet]') 
        dhClimb2 = Variable('dh_{climb2}', 'feet', 'Total Altitude Change Required in Climb 2')

        #Weights
        W_fuelClimb2 = VectorVariable(Nclimb2, 'W_{fuelClimb2}', 'N', 'Segment Fuel Weight')
        W_avgClimb2 = VectorVariable(Nclimb2, 'W_{avgClimb2}', 'N', 'Geometric Average of Segment Start and End Weight')

        #Engine
        TSFCc2 = VectorVariable(Nclimb2, 'TSFC_{c2}', '1/hr', 'Thrust Specific Fuel Consumption During Climb2')
        thrustc2 = VectorVariable(Nclimb2, 'thrust_{c2}', 'N', 'Thrust During Climb Segment #2')

        #Cruise engine thrust
        thrustcr2 = VectorVariable(Ncruise2, 'thrust_{cr2}', 'N', 'Thrust During Cruise Segment #2')

        #Fuselage area
        A_fuse = Variable('A_{fuse}', 'm^2', 'Estimated Fuselage Area')

        #non-GPkit variables
        #climb 2 lsit
        icl2 = map(int, np.linspace(0, Nclimb2 - 1, Nclimb2))

        constraints = []
        
        constraints.extend([            
            #set the velocity limits
            #needs to be replaced by an actual Vne and a mach number
            VClimb2[icl2] >= Vstall,

            VClimb2 == MClimb2 * aClimb2,

            #constraint on drag and thrust
            numeng*thrustc2[icl2] >= DClimb2[icl2] + W_avgClimb2[icl2] * thetaClimb2[icl2],
            
            #climb rate constraints
            TCS([excessPclimb2[icl2]+VClimb2[icl2]*DClimb2[icl2] <= VClimb2[icl2]*numeng*thrustc2[icl2]]),
            
            TCS([DClimb2[icl2] >= (.5*S*rhoClimb2[icl2]*VClimb2[icl2]**2)*(Cdwc2[icl2] + K*CLClimb2[icl2]**2) + Cdfuse * (.5 * A_fuse * rhoClimb2[icl2] * VClimb2[icl2]**2)]),
            
            K == (pi * e * AR)**-1,
            
            Cdwc2[icl2]**6.5 >= (1.02458748e10 * CLClimb2[icl2]**15.587947404823325 * MClimb2[icl2]**156.86410659495155 +
                2.85612227e-13 * CLClimb2[icl2]**1.2774976672501526 * MClimb2[icl2]**6.2534328002723703 +
                2.08095341e-14 * CLClimb2[icl2]**0.8825277088649582 * MClimb2[icl2]**0.0273667615730107 +
                1.94411925e+06 * CLClimb2[icl2]**5.6547413360261691 * MClimb2[icl2]**146.51920742858428),
            
            RCClimb2[icl2] == excessPclimb2[icl2]/W_avgClimb2[icl2],
            RCClimb2[icl2] >= 500*units('ft/min'),
            
            #make the small angle approximation and compute theta
            thetaClimb2[icl2]*VClimb2[icl2]  == RCClimb2[icl2],
           
            dhftClimb2[icl2]  == tminClimb2[icl2] * RCClimb2[icl2],

            tminClimb2 == thoursClimb2,
            #compute the distance traveled for each segment
            #takes into account two terms of a cosine expansion
            RngClimb2[icl2] == thoursClimb2[icl2]*VClimb2[icl2],

            W_avgClimb2[icl2] == .5*CLClimb2[icl2]*S*rhoClimb2[icl2]*VClimb2[icl2]**2,      
            WLoadClimb2[icl2] == .5*CLClimb2[icl2]*S*rhoClimb2[icl2]*VClimb2[icl2]**2/S,
            
            #compute fuel burn from TSFC
            W_fuelClimb2[icl2]  == numeng*TSFCc2[icl2] * thoursClimb2[icl2] * thrustc2[icl2],

            #compute the dh required for each climb 1 segment
            dhftClimb2[icl2] == dhClimb2/Nclimb2,

            #constrain the max wing loading
            WLoadClimb2 <= WLoadmax,

            thrustc2 <= 1.5 * thrustcr2[0],
            ])

        for i in range(0, Nclimb2):
            constraints.extend([
                #speed of sound
                aClimb2[i]  == (gamma * R * TClimb2[i])**.5,
                TSFCc2[i] == .6*units('1/hr'),
                ])
        
        Model.__init__(self, None, constraints, **kwargs)
        
#--------------------------------------
#cruise #1
#Breguet Range discretized to model the cruise


#---------------------------------------
#cruise climb/decent...might switch altitude in cruise

class Cruise2(Model):
    """
    class to model the second cruise portion of a flight (if the flight is
    long enough to mandate two cruise portions)
    Model is based off of a discretized Breguet Range equation
    """
    def __init__(self, Nclimb2, Ncruise2, **kwargs):
 
        #aero
        CLCruise2 = VectorVariable(Ncruise2, 'C_{L_{Cruise2}}', '-', 'Lift Coefficient')
        WLoadCruise2 = VectorVariable(Ncruise2, 'W_{Load_{Cruise2}}', 'N/m^2', 'Wing Loading')
        WLoadmax = Variable('W_{Load_max}', 'N/m^2', 'Max Wing Loading')
        DCruise2 = VectorVariable(Ncruise2, 'DragCruise2', 'N', 'Drag')
        Vstall = Variable('V_{stall}', 'knots', 'Aircraft Stall Speed')
        Cdwcr2 = VectorVariable(Ncruise2, 'C_{d_wcr2}', '-', 'Cd for a NC130 Airfoil at Re=2e7')
        cdp_rcr2 = VectorVariable(Ncruise2, 'cdp_rcr2', '-', 'Hold Variable for Drag Fit')
        Cdfuse = Variable('C_{d_fuse}', '-', 'Fuselage Drag Coefficient')
        K = Variable('K', '-', 'K for Parametric Drag Model')
        e = Variable('e', '-', 'Oswald Span Efficiency Factor')
        AR = Variable('AR', '-', 'Aspect Ratio')

        #atmosphere
        aCruise2 = VectorVariable(Ncruise2, 'aCruise2', 'm/s', 'Speed of Sound')
        rhoCruise2 = VectorVariable(Ncruise2, '\\rhoCruise2', 'kg/m^3', 'Air Density')
        pCruise2 = VectorVariable(Ncruise2, 'pCruise2', 'kPa', 'Pressure')
        muCruise2 = VectorVariable(Ncruise2, '\muCruise2', 'kg/m/s', 'Air Kinematic Viscosity')
        TCruise2 = VectorVariable(Ncruise2, 'TCruise2', 'K', 'Air Temperature')

        #aircraft geometry
        S = Variable('S', 'm^2', 'Wing Planform Area')

        #number of engines
        numeng = Variable('numeng', '-', 'Number of Engines')

        #time
        tminCruise2 = VectorVariable(Ncruise2, 'tminCruise2', 'min', 'Flight Time in Minutes')
        thoursCruise2 = VectorVariable(Ncruise2, 'thrCruise2', 'hour', 'Flight Time in Hours')

        #velocitites and mach numbers
        VCruise2 = VectorVariable(Ncruise2, 'VCruise2', 'knots', 'Aircraft Flight Speed')
        MCruise2 = VectorVariable(Ncruise2, 'MCruise2', '-', 'Aircraft Mach Number')

        #physical constants
        g = Variable('g', 9.81, 'm/s^2', 'Gravitational Acceleration')
        gamma = Variable('\gamma', 1.4, '-', 'Air Specific Heat Ratio')
        R = Variable('R', 287, 'J/kg/K', 'Gas Constant for Air')

        #Engine
        TSFCcr2 = VectorVariable(Ncruise2, 'TSFC_{cr2}', '1/hr', 'Thrust Specific Fuel Consumption During Cruise2')
        thrustcr2 = VectorVariable(Ncruise2, 'thrust_{cr2}', 'N', 'Thrust During Cruise Segment #2')
        
        #defined here for linking purposes
        T0 = Variable('T_0', 'K', 'Free Stream Stagnation Temperature')
        P0 = Variable('P_0', 'kPa', 'Free Stream Static Pressure')

        #parameter to make breguent range eqn gp compatible
        z_brec2 = VectorVariable(Ncruise2, 'z_{brec2}', '-', 'Breguet Parameter')

        #range variables
        ReqRngCruise = Variable('ReqRngCruise', 'miles', 'Required Cruise Range')
        RngCruise2 = VectorVariable(Ncruise2, 'RngCruise2', 'miles', 'Segment Range During Cruise2')
        
        #Weights
        W_fuelCruise2 = VectorVariable(Ncruise2, 'W_{fuelCruise2}', 'N', 'Segment Fuel Weight')
        W_avgCruise2 = VectorVariable(Ncruise2, 'W_{avgCruise2}', 'N', 'Geometric Average of Segment Start and End Weight')
        W_avgClimb2 = VectorVariable(Nclimb2, 'W_{avgClimb2}', 'N', 'Geometric Average of Segment Start and End Weight')
        W_endCruise2 = VectorVariable(Ncruise2, 'W_{endCruise2}', 'N', 'Segment End Weight')

        #Fuselage area
        A_fuse = Variable('A_{fuse}', 'm^2', 'Estimated Fuselage Area')

        #non-GPkit variables
        #cruise 2 lsit
        izbre = map(int, np.linspace(0, Ncruise2 - 1, Ncruise2))
        
        constraints = []
            
        constraints.extend([
##            MCruise2[izbre] == 0.8,

            MCruise2 * aCruise2 == VCruise2,
            
            TCS([DCruise2[izbre] >= (.5 * S * rhoCruise2[izbre] * VCruise2[izbre]**2) *
                 (Cdwcr2[izbre] + K * (W_avgCruise2[izbre] / (.5 * S * rhoCruise2[izbre]* VCruise2[izbre]**2))**2)
                 + Cdfuse * (.5 * A_fuse * rhoCruise2[izbre] * VCruise2[izbre]**2)]),

            K == (pi * e * AR)**-1,

            cdp_rcr2[izbre] >= (1.02458748e10 * CLCruise2[izbre]**15.587947404823325 * MCruise2[izbre]**156.86410659495155 +
                2.85612227e-13 * CLCruise2[izbre]**1.2774976672501526 * MCruise2[izbre]**6.2534328002723703 +
                2.08095341e-14 * CLCruise2[izbre]**0.8825277088649582 * MCruise2[izbre]**0.0273667615730107 +
                1.94411925e+06 * CLCruise2[izbre]**5.6547413360261691 * MCruise2[izbre]**146.51920742858428),

            Cdwcr2[izbre]**6.5 >= cdp_rcr2[izbre],
            
            DCruise2[izbre] == numeng * thrustcr2[izbre],

            W_avgCruise2[izbre] == .5 * CLCruise2[izbre] * S * rhoCruise2[izbre] * VCruise2[izbre]**2,
            WLoadCruise2[izbre] == .5 * CLCruise2[izbre] * S * rhoCruise2[izbre] * VCruise2[izbre]**2 / S,
            
            #taylor series expansion to get the weight term
            TCS([W_fuelCruise2[izbre]/W_endCruise2[izbre] >= te_exp_minus1(z_brec2[izbre], nterm=3)]),

            #breguet range eqn
            TCS([z_brec2[izbre] >= (numeng * TSFCcr2[izbre] * thoursCruise2[izbre] * DCruise2[izbre]) / W_avgCruise2[izbre]]),

            #time
            thoursCruise2[izbre]*VCruise2[izbre]  == RngCruise2[izbre],
            tminCruise2 == thoursCruise2,

            #constrain the max wing loading
            WLoadCruise2 <= WLoadmax,
            ])
        
        #constraint on the aircraft meeting the required range
        for i in range(min(izbre), max(izbre)+1):
            constraints.extend([
                TCS([RngCruise2[i] == ReqRngCruise/(Ncruise2)])
                ])

        for i in range(0, Ncruise2):
            constraints.extend([
                #speed of sound
                aCruise2[i]  == (gamma * R * TCruise2[i])**.5,
                TSFCcr2[i] == .5*units('1/hr'),
                ])
            
        Model.__init__(self, None, constraints, **kwargs)

#---------------------------------------
#decent

#----------------------------------------
#landing


#-------------------------------------
#build the linked model
class CommercialAircraft(Model):
    """
    class to link all models needed to simulate a commercial flight
    """
    def __init__(self, **kwargs):
        #defining the number of segments
        Nclimb1 = 2
        Nclimb2 = 2
        Ncruise2 = 3

        #total number of flight segments in each category, possibly useful
        Nclimb = Nclimb1 + Nclimb2
        Ncruise = Ncruise2

        Nseg = Nclimb + Ncruise

        #define all the submodels
        cmc = CommericalMissionConstraints(Nclimb1, Nclimb2, Ncruise2, False)
        climb1 = Climb1(Nclimb1, Ncruise2)
        climb2 = Climb2(Nclimb2, Ncruise2)
        cruise2 = Cruise2(Nclimb2, Ncruise2)

        atmvec = []
        
        for i in range(Nseg):
            atmvec.append(Atmosphere())
            
        substitutions = {      
            'V_{stall}': 120,
            'ReqRng': 2000,
            'hftCruise2': 40000,#('sweep', np.linspace(20000,40000,4)),
            'speedlimit': 250,
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
            "p_{sl}":101325,
            "T_{sl}":288.15,
            "L_{atm}":.0065,
            "M_{atm}":.0289644,
            "R_{atm}":8.31447
            }
        #for engine on design must link T0, P0, F_D,TSFC w/TSFC from icruise 2
        
        self.submodels = [cmc, climb1, climb2, cruise2]

        for i in range(len(atmvec)):
            self.submodels.extend(atmvec[i])

        constraints = ConstraintSet([self.submodels])

        subs= {}

        for i in range(Nclimb1):
            subs.update({
                climb1["\\rhoClimb1"][i]: atmvec[i]["\\rho"], climb1["TClimb1"][i]: atmvec[i]["T_{atm}"], cmc['hClimb1'][i]: atmvec[i]["h"]
                })

        for i in range(Nclimb2):
            subs.update({
                climb2["\\rhoClimb2"][i]: atmvec[i + Nclimb1]["\\rho"], climb2["TClimb2"][i]: atmvec[i + Nclimb1]["T_{atm}"], cmc['hClimb2'][i]: atmvec[i + Nclimb1]["h"]
                })

        for i in range(Ncruise2):
            subs.update({
                cruise2["\\rhoCruise2"][i]: atmvec[i + Nclimb1 + Nclimb2]["\\rho"], cruise2["TCruise2"][i]: atmvec[i + Nclimb1 + Nclimb2]["T_{atm}"],
                atmvec[i + Nclimb1 + Nclimb2]["h"]: cmc['hCruise2']
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
    sol = m.localsolve(solver="mosek", verbosity = 4, iteration_limit=100, skipsweepfailures=True)
    
##    sol, solhold = m.determine_unbounded_variables(m, solver="mosek",verbosity=4, iteration_limit=100)
    
