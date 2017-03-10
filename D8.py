"""Simple commercial aircraft flight profile and D8 aircraft model"""
""" Combines Wing, VerticalTail, and Fuselage models for D8"""

import numpy as np
from gpkit import Variable, Model, units, SignomialsEnabled, SignomialEquality, Vectorize
from gpkit.constraints.tight import Tight as TCS
from gpkit.tools import te_exp_minus1
from numpy import pi, cos, tan

TCS.reltol = 1e-3

# only needed for plotting
import matplotlib.pyplot as plt

# importing from D8_integration
from stand_alone_simple_profile import FlightState
from TASOPT_VT_yaw_rate_and_EO_simple_profile import VerticalTail
from D8_HT_simple_profile import HorizontalTail
from D8_Wing_simple_profile import Wing
from turbofan.engine_validation import Engine
from D8_Fuselage import Fuselage

#import constant relaxtion tool
from relaxed_constants import relaxed_constants, post_process

#import tool to check solution relative to TASOPT
from D8_TASOPT_percent_diff import percent_diff

"""
Models required to minimize the aircraft total fuel weight. Rate of climb equation taken from John
Anderson's Aircraft Performance and Design (eqn 5.85).
Inputs
-----
- Number of passengers
- Passenger weight [N]
- Fuselage area per passenger (recommended to use 1 m^2 based on research) [m^2]
- Engine weight [N]
- Number of engines
- Required mission range [nm]
- Oswald efficiency factor
- Max allowed wing span [m]
- Cruise altitude [ft]
Sources for substitutions:
-[b757 freight doc]
-[Boeing]
-[Philippe]
-[stdAtm]
-[TAS]
Other markers:
-[SP]
-[SPEquality]
"""

# Script for doing sweeps
n = 5
sweeps = False
sweepSMmin = False
sweepdxCG = True
sweepReqRng = True
sweepthetadb = True
sweepxCG = True
sweepCruiseAlt = True
sweepMmin = True
sweepnpax = True
sweepResFuel = True

plot = True

# Only one active at a time
D80 = False
D82 = False
b737800 = True

sweep = 27.566#30 [deg]

if D82:
     sweep = 13.237  # [deg]

if b737800:
     sweep = 26.0 # [deg]

g = 9.81 * units('m*s**-2')

class Aircraft(Model):
    "Aircraft class"

    def setup(self, Nclimb, Ncruise, enginestate, eng, Nfleet=0, **kwargs):
        # create submodels
        self.fuse = Fuselage()
        self.wing = Wing()
        if Nfleet != 0:
            self.engine = Engine(0, True, Nclimb+Ncruise, enginestate, eng, Nfleet)
        else:
           self.engine = Engine(0, True, Nclimb+Ncruise, enginestate, eng)
        self.VT = VerticalTail()
        self.HT = HorizontalTail()

        # variable definitions
        numaisle = Variable('numaisle','-','Number of Aisles')
        numeng = Variable('numeng', '-', 'Number of Engines')
        numVT = Variable('numVT','-','Number of Vertical Tails')
        Vne = Variable('V_{ne}',144, 'm/s', 'Never-exceed speed')  # [Philippe]
        rhoTO = Variable('\\rho_{T/O}',1.225,'kg*m^-3','Air density at takeoff')
        ReserveFraction = Variable('ReserveFraction', '-', 'Fuel Reserve Fraction')
        
        SMmin = Variable('SM_{min}','-', 'Minimum Static Margin')
        dxCG = Variable('\\Delta x_{CG}', 'm', 'Max CG Travel Range')
        xCGmin = Variable('x_{CG_{min}}','m','Maximum Forward CG')

        Izwing = Variable('I_{z_{wing}}','kg*m**2','Wing moment of inertia')
        Iztail = Variable('I_{z_{tail}}','kg*m**2','Tail moment of inertia')
        Izfuse = Variable('I_{z_{fuse}}','kg*m**2','Fuselage moment of inertia')

        Mmin = Variable('M_{min}','-','Minimum Cruise Mach Number')

        # Weights
        W_total = Variable('W_{total}', 'lbf', 'Total Aircraft Weight')
        W_dry = Variable('W_{dry}', 'lbf', 'Zero Fuel Aircraft Weight')
        W_ftotal = Variable('W_{f_{total}}', 'lbf', 'Total Fuel Weight')
        W_fclimb = Variable('W_{f_{climb}}', 'lbf',
                            'Fuel Weight Burned in Climb')
        W_fcruise = Variable('W_{f_{cruise}}', 'lbf',
                             'Fuel Weight Burned in Cruise')
        W_fprimary = Variable('W_{f_{primary}}', 'lbf', 'Total Fuel Weight Less Fuel Reserves')

        Wwing = Variable('W_{wing}','lbf','Wing Weight')
        WHT = Variable('W_{HT}','lbf','Horizontal Tail Weight')
        WVT = Variable('W_{VT}','lbf','Vertical Tail Weight')

        #engine system weight variables
        rSnace = Variable('rSnace', '-', 'Nacelle and Pylon Wetted Area')
        lnace = Variable('l_{nacelle}', 'm', 'Nacelle Length')
        fSnace = Variable('f_{S_nacelle}', '-', 'Non-dimensional Nacelle Area')
        Snace = Variable('S_{nace}', 'm^2', 'Nacelle Surface Area')
        Ainlet = Variable('A_{inlet}','m^2', 'Inlet Area')
        Afancowl = Variable('A_{fancowl}', 'm^2', 'Fan Cowling Area')
        Aexh = Variable('A_{exh}', 'm^2', 'Exhaust Area')
        Acorecowl = Variable('A_{corecowl}', 'm^2', 'Core Cowling Area')
        Wnace = Variable('W_{nace}', 'lbf', 'Nacelle Weight')
        Wpylon = Variable('W_{pylon}', 'lbf','Engine Pylon Weight')
        fpylon = Variable('f_{pylon}', '-', 'Pylong Weight Fraction')
        feadd = Variable('f_{eadd}', '-', 'Additional Engine Weight Fraction')
        Weadd = Variable('W_{eadd}', 'lbf', 'Additional Engine System Weight')
        Wengsys = Variable('W_{engsys}', 'lbf', 'Total Engine System Weight')
        rvnace = Variable('r_{vnace}', '-', 'Incoming Nacelle Velocity Ratio')
     
        constraints = []
        with SignomialsEnabled():
            constraints.extend([
                            self.wing['c_{root}'] == self.fuse['c_0'],
                            self.wing.wb['wwb'] == self.fuse['wtc'],
                            self.wing['x_w'] == self.fuse['x_{wing}'],
                            self.wing['V_{ne}'] == 144*units('m/s'),
                            self.VT['V_{ne}'] == 144*units('m/s'),

                            #compute the aircraft's zero fuel weight
                            TCS([self.fuse['W_{fuse}'] + numeng \
                                * Wengsys + self.fuse['W_{tail}'] + Wwing <= W_dry]),

                            # Total takeoff weight constraint
                            TCS([W_ftotal + W_dry + self.fuse['W_{payload}'] <= W_total]),
                            TCS([W_ftotal >= W_fprimary + ReserveFraction * W_fprimary]),
                            TCS([W_fprimary >= W_fclimb + W_fcruise]),

                            # Load factor matching
                            self.fuse['N_{lift}'] == self.wing['N_{lift}'],
                            #set the wing lift
                            self.wing['L_{max}'] >= self.wing['N_{lift}'] * W_total + self.HT['L_{{max}_h}'],

                            # Wing fuel constraints
                            self.wing['W_{fuel_{wing}}'] == W_ftotal/self.wing['FuelFrac'],

                            # Lifting surface weights
                            Wwing == self.wing['W_{wing_system}'],
                            WHT == self.HT['W_{HT_system}'],
                            WVT == self.VT['W_{VT_system}'],

                            # Tail cone sizing
                            3 * self.VT['M_r'] * self.VT['c_{root_{vt}}'] * \
                                (self.fuse['p_{\\lambda_v}'] - 1) >= numVT*self.VT[
                                    'L_{v_{max}}'] * self.VT['b_{vt}'] * (self.fuse['p_{\\lambda_v}']),
                            TCS([self.fuse['V_{cone}'] * (1 + self.fuse['\\lambda_{cone}']) * \
                             (pi + 4 * self.fuse['\\theta_{db}']) >= numVT*self.VT[
                                'M_r'] * self.VT['c_{root_{vt}}'] / self.fuse['\\tau_{cone}'] * \
                                 (pi + 2 * self.fuse['\\theta_{db}']) * \
                                  (self.fuse['l_{cone}'] / self.fuse['R_{fuse}'])]), #[SP]

                            # Lift curve slope ratio for HT and Wing
                            SignomialEquality(self.HT['m_{ratio}']*(1+2/self.wing['AR']), 1 + 2/self.HT['AR_h']),

                            # HT Location and Volume Coefficient
                            self.HT['x_{CG_{ht}}'] <= self.fuse['l_{fuse}'],
                            self.fuse['x_{tail}'] == self.VT['x_{CG_{vt}}'],
                            TCS([self.HT['V_{h}'] == self.HT['S_h']*self.HT['l_{ht}']/(self.wing['S']*self.wing['mac'])]),
                            # self.HT['V_{h}'] >= 0.4,

                            # HT Max Loading
                            TCS([self.HT['L_{{max}_h}'] >= 0.5*rhoTO*Vne**2*self.HT['S_h']*self.HT['C_{L_{hmax}}']]),

                            # Tail weight
                            self.fuse['W_{tail}'] >= numVT*WVT + WHT + self.fuse['W_{cone}'],

                            # VT root chord constraint #TODO find better constraint
                            self.VT['c_{root_{vt}}'] <= 1.5*self.fuse['l_{cone}'],

                            # VT volume coefficient
                            self.VT['V_{vt}'] == numVT*self.VT['S_{vt}'] * self.VT['x_{CG_{vt}}']/(self.wing['S']*self.wing['b']),

                            # VT sizing constraints
                            # Yaw rate constraint at flare
                            numVT*.5*self.VT['\\rho_{TO}']*self.VT['V_{land}']**2*self.VT['S_{vt}']*self.VT['l_{vt}']* \
                                            self.VT['C_{L_{vyaw}}'] >= self.VT['\\dot{r}_{req}']*self.VT['I_{z}'],

                            # Force moment balance for one engine out condition
                            numVT*self.VT['L_{vtEO}']*self.VT['l_{vt}'] >= self.VT['T_e']*self.VT['y_{eng}'] + \
                                        self.VT['D_{wm}']*self.VT['y_{eng}'],
                            # TASOPT 2.0 p45

                            # Vertical bending material coefficient (VT aero loads)
                            self.fuse['B1v'] == self.fuse['r_{M_v}']*numVT*self.VT['L_{v_{max}}']/(self.fuse['w_{fuse}']*self.fuse['\\sigma_{M_v}']),

                            # Moment of inertia around z-axis
                            SignomialEquality(self.VT['I_{z}'], Izwing + Iztail + Izfuse),

                            # Fuselage width (numaisle comes in)
                            TCS([2.*self.fuse['w_{fuse}'] >= self.fuse['SPR'] * self.fuse['w_{seat}'] + \
                                 numaisle*self.fuse['w_{aisle}'] + 2 * self.fuse['w_{sys}'] + self.fuse['t_{db}']]),

                            #engine system weight constraints
                            Snace == rSnace * np.pi * 0.25 * self.engine['d_{f}']**2,
                            lnace == 0.15 * self.engine['d_{f}'] * rSnace,
                            fSnace == Snace * self.wing['S']**-1,
                            Ainlet == 0.4 * Snace,
                            Afancowl == 0.2 * Snace,
                            Aexh == 0.4 * Snace,
                            Acorecowl == 3 * np.pi * self.engine['d_{LPC}']**2,
                            TCS([Wnace >= ((2.5+ 0.238*self.engine['d_{f}']/units('in')) * Ainlet + 1.9*Afancowl \
                                + (2.5+ 0.0363*self.engine['d_{f}']/units('in'))*Aexh + 1.9*Acorecowl)*units('lbf/ft^2'),
                            Weadd == feadd * self.engine['W_{engine}']]),
                            TCS([Wpylon >= (Wnace + Weadd + self.engine['W_{engine}']) * fpylon]),
                            TCS([Wengsys >= Wpylon + Wnace + Weadd + self.engine['W_{engine}']]),
                            ])

        #d8 only constraints
        if D80 or D82:
            with SignomialsEnabled():
                constraints.extend([

                    # VT height constraint (4*engine radius)
                    self.VT['b_{vt}'] >= 2 * self.engine['d_{f}'],

                    # Engine out moment arm,
                    self.VT['y_{eng}'] == 0.5*self.fuse['w_{fuse}'],

                    # HT root moment
                    self.HT['M_r']*self.HT['c_{root_h}'] >= self.HT['N_{lift}']*self.HT['L_{h_{rect}}']*(self.HT['b_{ht}']/4) \
                        + self.HT['N_{lift}']*self.HT['L_{h_{tri}}']*(self.HT['b_{ht}']/6) - self.HT['N_{lift}']*self.fuse['w_{fuse}']*self.HT['L_{{max}_h}']/2., #[SP]


                    # Pin VT joint moment constraint #TODO may be problematic, should check
                    SignomialEquality(self.HT['L_{h_{rect}}']*(self.HT['b_{ht}']/4 - self.fuse['w_{fuse}']),
                                    self.HT['L_{h_{tri}}']*(self.fuse['w_{fuse}'] - self.HT['b_{ht}']/6)), #[SP] #[SPEquality]

                    # HT/VT joint constraint
                    self.HT['b_{ht}']/(2.*self.fuse['w_{fuse}'])*self.HT['\lambda_h']*self.HT['c_{root_h}'] == self.HT['c_{attach}'],


                # Moment of inertia
                Izwing >= (self.wing['W_{fuel_{wing}}'] + Wwing)/(self.wing['S']*g)* \
                                    self.wing['c_{root}']*self.wing['b']**3*(1./12.-(1-self.wing['\\lambda'])/16), #[SP]
                Iztail >= (self.fuse['W_{apu}'] + numeng*Wengsys + self.fuse['W_{tail}'])*self.VT['l_{vt}']**2/g,
                            #NOTE: Using xwing as a CG surrogate. Reason: xCG moves during flight; want scalar Izfuse
                Izfuse >= (self.fuse['W_{fuse}'] + self.fuse['W_{payload}'])/self.fuse['l_{fuse}'] * \
                                    (self.fuse['x_{wing}']**3 + self.VT['l_{vt}']**3)/(3.*g),

                # Floor loading
                self.fuse['S_{floor}'] == (5. / 16.) * self.fuse['P_{floor}'],
                self.fuse['M_{floor}'] == 9. / 256. * self.fuse['P_{floor}']*self.fuse['w_{floor}'],


                # Horizontal tail aero+landing loads constants A1h
                self.fuse['A1h_{Land}'] >= (self.fuse['N_{land}'] * \
                                (self.fuse['W_{tail}'] + numeng*Wengsys + self.fuse['W_{apu}'])) / \
                                 (self.fuse['h_{fuse}'] * self.fuse['\\sigma_{M_h}']),

                self.fuse['A1h_{MLF}'] >= (self.fuse['N_{lift}'] * \
                                (self.fuse['W_{tail}'] + numeng*Wengsys + self.fuse['W_{apu}']) \
                                + self.fuse['r_{M_h}'] * self.HT['L_{{max}_h}']) / \
                                 (self.fuse['h_{fuse}'] * self.fuse['\\sigma_{M_h}']),
                ])

          #737 only constraints
        if b737800:
            with SignomialsEnabled():
               constraints.extend([
                   # Engine out moment arm,
                    self.VT['y_{eng}'] == 4.83*units('m'),

                    # HT root moment
                    # TCS([self.HT['M_r'] >= self.HT['L_{{max}_h}']*self.HT['AR_h']*self.HT['p_{ht}']/24]),
                    TCS([self.HT['M_r']*self.HT['c_{root_h}'] >= 1./6.*self.HT['L_{h_{tri}}']*self.HT['b_{ht}'] + \
                         1./4.*self.HT['L_{h_{rect}}']*self.HT['b_{ht}']]),

                    # HT joint constraint
                   self.HT['c_{attach}'] == self.HT['c_{root_h}'],

                   # Moment of inertia
                    Izwing >= numeng*Wengsys*self.VT['y_{eng}']**2/g + \
                                    (self.wing['W_{fuel_{wing}}'] + Wwing)/(self.wing['S']*g)* \
                                    self.wing['c_{root}']*self.wing['b']**3*(1./12.-(1-self.wing['\\lambda'])/16), #[SP]
                    Iztail >= (self.fuse['W_{apu}'] + self.fuse['W_{tail}'])*self.VT['l_{vt}']**2/g,
                            #NOTE: Using xwing as a CG surrogate. Reason: xCG moves during flight; want scalar Izfuse
                    Izfuse >= (self.fuse['W_{fuse}'] + self.fuse['W_{payload}'])/self.fuse['l_{fuse}'] * \
                                    (self.fuse['x_{wing}']**3 + self.VT['l_{vt}']**3)/(3.*g),

                   # Floor loading
                    self.fuse['S_{floor}'] == 1./2. * self.fuse['P_{floor}'],
                    self.fuse['M_{floor}'] == 1./4. * self.fuse['P_{floor}']*self.fuse['w_{floor}'],

                    # Horizontal tail aero+landing loads constants A1h
                    self.fuse['A1h_{Land}'] >= (self.fuse['N_{land}'] * \
                                (self.fuse['W_{tail}'] + self.fuse['W_{apu}'])) / \
                                 (self.fuse['h_{fuse}'] * self.fuse['\\sigma_{M_h}']),

                    self.fuse['A1h_{MLF}'] >= (self.fuse['N_{lift}'] * \
                                (self.fuse['W_{tail}'] + self.fuse['W_{apu}']) \
                                + self.fuse['r_{M_h}'] * self.HT['L_{{max}_h}']) / \
                                 (self.fuse['h_{fuse}'] * self.fuse['\\sigma_{M_h}']),
                    ])

        self.components = [self.fuse, self.wing, self.engine, self.VT, self.HT]

        return self.components, constraints

    def climb_dynamic(self, state, Nclimb):  # creates an aircraft climb performance model, given a state
        return ClimbP(self, state, Nclimb)

    def cruise_dynamic(self, state, Nclimb): # creates an aircraft cruise performance model, given a state
        return CruiseP(self, state, Nclimb)


class AircraftP(Model):
    """
    Aircraft performance models superclass, contains constraints true for
    all flight segments
    """

    def setup(self, aircraft, state):
        # make submodels
        self.aircraft = aircraft
        self.wingP = aircraft.wing.dynamic(state)
        self.fuseP = aircraft.fuse.dynamic(state)
        self.VTP = aircraft.VT.dynamic(aircraft.fuse,state)
        self.HTP = aircraft.HT.dynamic(aircraft.fuse,aircraft.wing,state)
        self.Pmodels = [self.wingP, self.fuseP, self.VTP, self.HTP]

        # variable definitions
        Vstall = Variable('V_{stall}',120, 'knots', 'Aircraft Stall Speed')
        D = Variable('D', 'N', 'Total Aircraft Drag')
        C_D = Variable('C_D', '-', 'Total Aircraft Drag Coefficient')
        LoD = Variable('L/D','-','Lift-to-Drag Ratio')
        W_avg = Variable(
            'W_{avg}', 'lbf', 'Geometric Average of Segment Start and End Weight')
        W_start = Variable('W_{start}', 'lbf', 'Segment Start Weight')
        W_end = Variable('W_{end}', 'lbf', 'Segment End Weight')
        W_burn = Variable('W_{burn}', 'lbf', 'Segment Fuel Burn Weight')
        WLoadmax = Variable('W_{Load_max}',6664, 'N/m^2', 'Max Wing Loading')
        WLoad = Variable('W_{Load}', 'N/m^2', 'Wing Loading')
        t = Variable('tmin', 'min', 'Segment Flight Time in Minutes')
        thours = Variable('thr', 'hour', 'Segment Flight Time in Hours')

        xAC = Variable('x_{AC}','m','Aerodynamic Center of Aircraft')
        xCG = Variable('x_{CG}','m','Center of Gravity of Aircraft')
        xNP = Variable('x_{NP}','m','Neutral Point of Aircraft')
        # SM = Variable('SM','-','Stability Margin of Aircraft')
        PCFuel = Variable('PCFuel','-','Percent Fuel Remaining (end of segment)')

        Pcabin = Variable('P_{cabin}','Pa','Cabin Air Pressure')
        W_buoy = Variable('W_{buoy}','lbf','Buoyancy Weight')
        Tcabin = Variable('T_{cabin}','K','Cabin Air Temperature')
        rhocabin = Variable('\\rho_{cabin}','kg/m^3','Cabin Air Density')

        #variables for nacelle drag calcualation
        Vnace = Variable('V_{nacelle}', 'm/s', 'Incoming Nacelle Flow Velocity')
        V2 = Variable('V_2', 'm/s', 'Interior Nacelle Flow Velcoity')
        Vnacrat = Variable('V_{nacelle_ratio}', '-', 'Vnle/Vinf')
        rvnsurf = Variable('rvnsurf', '-', 'Intermediate Nacelle Drag Parameter')
        Cfnace = Variable('C_{f_nacelle}', '-', 'Nacelle Drag Coefficient')
        Renace = Variable('R_{e_nacelle}', '-', 'Nacelle Reynolds Number')
        Cfturb = Variable('C_{f_nacelle}', '-', 'Turbulent Nacelle Skin Friction Coefficient')
        Cdnace = Variable('C_{d_nacelle}', '-', 'Nacelle Drag Coeffecient')
        Dnace = Variable('D_{nacelle}', 'N', 'Drag On One Nacelle')

        constraints = []

        with SignomialsEnabled():
            constraints.extend([
            W_burn == W_burn,
            PCFuel == PCFuel,

            # Cabin Air properties
            rhocabin == Pcabin/(state['R']*Tcabin),
            Pcabin == 75000*units('Pa'),
            Tcabin == 297*units('K'),

            # Buoyancy weight #TODO relax the equality
            SignomialEquality(W_buoy,(rhocabin - state['\\rho'])*g*aircraft['V_{cabin}']),  #[SP] #[SPEquality]
            # W_buoy >= (rhocabin - state['\\rho'])*g*aircraft['V_{cabin}'], # [SP]

            # speed must be greater than stall speed
            state['V'] >= Vstall,

            # compute the drag
            D >= self.wingP['D_{wing}'] + self.fuseP['D_{fuse}'] + self.aircraft['numVT']*self.VTP['D_{vt}'] + self.HTP['D_{ht}'] + aircraft['numeng'] * Dnace,
            C_D == D/(.5*state['\\rho']*state['V']**2 * self.aircraft.wing['S']),
            LoD == W_avg/D,

            # Wing looading
            WLoad == .5 * self.wingP['C_{L}'] * self.aircraft['S'] * state.atm['\\rho'] * state['V']**2 / self.aircraft.wing['S'],

            # Geometric average of start and end weights of flight segment
            W_avg >= (W_start * W_end)**.5 + W_buoy, # Buoyancy weight included in Breguet Range

            # Maximum wing loading constraint
            WLoad <= WLoadmax,

            # Flight time unit conversion
            t == thours,

            #VTP constraints
            TCS([aircraft.fuse['l_{fuse}'] >= aircraft.VT['\\Delta x_{trail_v}'] + xCG]),
            TCS([aircraft.VT['x_{CG_{vt}}'] <= xCG + (aircraft.VT['\\Delta x_{lead_v}']+aircraft.VT['\\Delta x_{trail_v}'])/2]),
            aircraft.VT['x_{CG_{vt}}'] <= aircraft.fuse['l_{fuse}'],

            # Drag of a windmilling engine (VT sizing)
            TCS([aircraft.VT['D_{wm}'] >= 0.5*aircraft.VT['\\rho_{TO}']*aircraft.VT['V_1']**2*aircraft.engine['A_2']*aircraft.VT['C_{D_{wm}}']]),


            # Aircraft trim conditions
            TCS([xAC/aircraft.wing['mac'] <= xCG/aircraft.wing['mac'] + \
                 self.wingP['c_{m_{w}}']/self.wingP['C_{L}']  +\
                              aircraft.HT['V_{h}']*(self.HTP['C_{L_h}']/self.wingP['C_{L}'])]),

            # Tail aspect ratio and lift constraints
            aircraft.HT['AR_h'] >= 6, #TODO change to tip Re constraint
            self.HTP['C_{L_h}'] >= 0.01, #TODO remove

            # HT/VT moment arm constraints
            aircraft.HT['l_{ht}'] <= aircraft.HT['x_{CG_{ht}}'] - xCG,
            aircraft.VT['l_{vt}'] <= aircraft.VT['x_{CG_{vt}}'] - xCG,

           # Tail downforce penalty to wing lift
            self.wingP['L_w'] >= W_avg + self.HTP['L_h'],

            # Wing location and AC constraints
            TCS([xCG + self.HTP['\\Delta x_{{trail}_h}'] <= aircraft.fuse['l_{fuse}']]), #TODO tighten
            # xAC == aircraft['x_{wing}'], #TODO improve, only works because cmw == 0.1
            TCS([xAC >= aircraft['x_{wing}'] + xNP]),
            SignomialEquality(xAC,xCG + self.HTP['\\Delta x_w']),

            # Neutral point approximation (taken from Basic Aircraft Design Rules, Unified)
            # TODO improve
            SignomialEquality(xNP/aircraft['mac']/aircraft['V_{h}']*(aircraft['AR']+2.)*(1.+2./aircraft['AR_h']),
                              (1.+2./aircraft['AR'])*(aircraft['AR']-2)),

            TCS([aircraft.HT['x_{CG_{ht}}'] <= xCG + 0.5*(self.HTP['\\Delta x_{{trail}_h}'] + self.HTP['\\Delta x_{{lead}_h}'])]), #TODO tighten

            # Static margin constraint with and without dxCG
            self.wingP['c_{m_{w}}'] == .65,
            # SM >= aircraft['SM_{min}'],
            TCS([aircraft['SM_{min}'] + aircraft['\\Delta x_{CG}']/aircraft.wing['mac'] \
                 + self.wingP['c_{m_{w}}']/aircraft.wing['C_{L_{wmax}}'] <= \
                                            aircraft.HT['V_{h}']*aircraft.HT['m_{ratio}'] +\
                                            aircraft.HT['V_{h}']*aircraft.HT['C_{L_{hmax}}']/aircraft.wing['C_{L_{wmax}}']]), # [SP]

          #nacelle drag
          Renace == state['\\rho']*state['V'] * aircraft['l_{nacelle}']/state['\\mu'],
          Cfnace == 4*0.0743/(Renace**(0.2)), #from http://www.calpoly.edu/~kshollen/ME347/Handouts/Friction_Drag_Coef.pdf
          Vnace == aircraft['r_{vnace}'] * state['V'],
          Vnacrat >= 2*Vnace/state['V'] - V2/state['V'],
          rvnsurf**3 >= 0.25*(Vnacrat + aircraft['r_{vnace}'])*(Vnacrat**2 + aircraft['r_{vnace}']**2),
          Cdnace == aircraft['f_{S_nacelle}'] * Cfnace[0] * rvnsurf **3,
          Dnace == Cdnace * 0.5 * state['\\rho'] * state['V']**2 * aircraft['S'],
           ])


            # CG CONSTRAINT #TODO improve; how to account for decreasing fuel volume?
            if D80 or D82:
                constraints.extend([
                    # TCS([xCG*W_start >= 0.5*(aircraft.fuse['W_{fuse}']+aircraft.fuse['W_{payload}'])*aircraft.fuse['l_{fuse}'] \
                    # + (aircraft['W_{tail}']+aircraft['numeng']*aircraft['W_{engsys}'])*aircraft['x_{tail}'] \
                    # + (aircraft['W_{wing_system}']+(PCFuel + aircraft['ReserveFraction'])*aircraft['W_{f_{primary}}'])*aircraft.fuse['x_{wing}']]),
                    # # Center of gravity constraints #TODO Refine
                    # xCG >= aircraft['x_{CG_{min}}'],

                ])

            if b737800:
                with SignomialsEnabled():
                    constraints.extend([
                         # SignomialEquality(xCG*aircraft['W_{dry}'], 0.5*(aircraft.fuse['W_{fuse}']+aircraft.fuse['W_{payload}'])*aircraft.fuse['l_{fuse}'] \
                         #    + (aircraft['W_{tail}'])*aircraft['x_{tail}'] \
                         #    + aircraft['W_{wing_system}']*aircraft.fuse['x_{wing}'] \
                         #    + aircraft['numeng']*aircraft['W_{engsys}'] \
                         #    *(aircraft.fuse['x_{wing}'] + 0.3*aircraft.VT['y_{eng}']*aircraft.wing['\\tan(\\Lambda)'])),
                         #    # Center of gravity constraints #TODO Refine
                         #    xCG >= aircraft['x_{CG_{min}}'],
                #cap max rear wing position
                # aircraft.fuse['x_{wing}'] <= 18.94*units('m'),
                    ])
        return self.Pmodels, constraints

class ClimbP(Model): # Climb performance constraints

    def setup(self, aircraft, state, Nclimb, **kwargs):
        # submodels
        self.aircraft = aircraft
        self.aircraftP = AircraftP(aircraft, state)
        self.wingP = self.aircraftP.wingP
        self.fuseP = self.aircraftP.fuseP
        self.engine = aircraft.engine

        # variable definitions
        theta = Variable('\\theta', '-', 'Aircraft Climb Angle')
        excessP = Variable('excessP', 'W', 'Excess Power During Climb')
        RC = Variable('RC', 'feet/min', 'Rate of Climb/Decent')
        dhft = Variable(
            'dhft', 'feet', 'Change in Altitude Per Climb Segment [feet]')
        RngClimb = Variable('RngClimb', 'nautical_miles',
                            'Down Range Covered in Each Climb Segment')

        # constraints
        constraints = []

        constraints.extend([
            # Thrust >= Drag + Vertical Potential Energy
            self.aircraft['numeng'] * self.engine['F'][:Nclimb] >= self.aircraftP[
                'D'] + self.aircraftP['W_{avg}'] * theta,

            # Excess power for climb
           TCS([excessP + state['V'] * self.aircraftP['D'] <= state['V']
                 * aircraft['numeng'] * self.engine['F'][:Nclimb]]),

            RC == excessP / self.aircraftP['W_{avg}'],
            RC >= 500 * units('ft/min'),

            # Climb angle and rate constraint
            theta * state['V'] == RC,

            dhft == self.aircraftP['tmin'] * RC,

            # Small angle assumption during climb
            RngClimb == self.aircraftP['thr'] * state['V'],
        ])

        return constraints + self.aircraftP


class CruiseP(Model): # Cruise performance constraints

    def setup(self, aircraft, state, Nclimb, **kwargs):
        self.aircraft = aircraft
        self.aircraftP = AircraftP(aircraft, state)
        self.wingP = self.aircraftP.wingP
        self.fuseP = self.aircraftP.fuseP
        self.engine = aircraft.engine
        
        # variable definitions
        z_bre = Variable('z_{bre}', '-', 'Breguet Parameter')
        Rng = Variable('Rng', 'nautical_miles', 'Cruise Segment Range')

        constraints = []

        constraints.extend([
            # Steady level flight constraint on D
            self.aircraftP['D'] == aircraft[
                'numeng'] * self.engine['F'][Nclimb:],

            # Taylor series expansion to get the weight term
            TCS([self.aircraftP['W_{burn}'] / self.aircraftP['W_{end}'] >=
                 te_exp_minus1(z_bre, nterm=3)]),

            # Breguet range eqn
            TCS([z_bre >= (self.engine['TSFC'][Nclimb:] * self.aircraftP['thr'] *
                           self.aircraftP['D']) / self.aircraftP['W_{avg}']]),

            # Time
            self.aircraftP['thr'] * state['V'] == Rng,
        ])

        return constraints + self.aircraftP


class CruiseSegment(Model): # Combines FlightState and Aircraft to form a cruise flight segment
    def setup(self, aircraft, Nclimb, **kwargs):
        self.state = FlightState()
        self.cruiseP = aircraft.cruise_dynamic(self.state, Nclimb)
        return self.state, self.cruiseP


class ClimbSegment(Model): # Combines FlightState and Aircraft to form a climb flight segment
    def setup(self, aircraft, Nclimb, **kwargs):
        self.state = FlightState()
        self.climbP = aircraft.climb_dynamic(self.state, Nclimb)
        return self.state, self.climbP

class StateLinking(Model):
    """
    link all the state model variables
    """
    def setup(self, climbstate, cruisestate, enginestate, Nclimb, Ncruise):
        statevarkeys = ['p_{sl}', 'T_{sl}', 'L_{atm}', 'M_{atm}', 'P_{atm}', 'R_{atm}',
                        '\\rho', 'T_{atm}', '\\mu', 'T_s', 'C_1', 'h', 'hft', 'V', 'a', 'R', '\\gamma', 'M']
        constraints = []
        for i in range(len(statevarkeys)):
            varkey = statevarkeys[i]
            for i in range(Nclimb):
                constraints.extend([
                    climbstate[varkey][i] == enginestate[varkey][i]
                    ])
            for i in range(Ncruise):
                constraints.extend([
                    cruisestate[varkey][i] == enginestate[varkey][i+Nclimb]
                    ])           
        
        return constraints

class Mission(Model):
    """
    Mission superclass, links together all subclasses into an optimization problem
    """

    def setup(self, Nclimb, Ncruise, **kwargs):
        # define the number of each flight segment

        if D80 or D82:
             eng = 3

        if b737800:
             eng = 1
        
        # vectorize
        with Vectorize(Nclimb + Ncruise):
            enginestate = FlightState()

        # build required submodels
        aircraft = Aircraft(Nclimb, Ncruise, enginestate, eng)

        # vectorize
        with Vectorize(Nclimb):
            climb = ClimbSegment(aircraft, Nclimb)

        with Vectorize(Ncruise):
            cruise = CruiseSegment(aircraft, Nclimb)

        statelinking = StateLinking(climb.state, cruise.state, enginestate, Nclimb, Ncruise)

        # declare new variables
        CruiseAlt = Variable('CruiseAlt', 'ft', 'Cruise Altitude [feet]')
        ReqRng = Variable('ReqRng', 'nautical_miles', 'Required Cruise Range')
        RngCruise = Variable('RngCruise', 'nautical_miles', 'Total Cruise Range')

        #hold variables
        hftCruiseHold = Variable('hftCruiseHold', 'ft', 'Temp Variable to Avoid Dimension Mismatch')

        h = climb.state['h']
        hftClimb = climb.state['hft']
        dhft = climb.climbP['dhft']
        hftCruise = cruise.state['hft']

        # make overall constraints
        constraints = []

        with SignomialsEnabled():
            #CG constraints
            if D80 or D82:
                constraints.extend([
                TCS([climb['x_{CG}']*climb['W_{end}'] >=
                    0.5*(aircraft.fuse['W_{fuse}']+aircraft.fuse['W_{payload}'])*aircraft.fuse['l_{fuse}'] \
                    + (aircraft['W_{tail}']+aircraft['numeng']*aircraft['W_{engsys}'])*aircraft['x_{tail}'] \
                    + (aircraft['W_{wing_system}']+(climb['PCFuel']+aircraft['ReserveFraction'])*aircraft['W_{f_{primary}}'])*aircraft.fuse['x_{wing}'] \
                    ]), # TODO improve; using x_b as a surrogate for xeng
##               cruise['x_{CG}'][0] <= climb['x_{CG}'],
                TCS([cruise['x_{CG}']*cruise['W_{end}'] >=
                    0.5*(aircraft.fuse['W_{fuse}']+aircraft.fuse['W_{payload}'])*aircraft.fuse['l_{fuse}'] \
                    + (aircraft['W_{tail}']+aircraft['numeng']*aircraft['W_{engsys}'])*aircraft['x_{tail}'] \
                    + (aircraft['W_{wing_system}']+(cruise['PCFuel']+aircraft['ReserveFraction'])*aircraft['W_{f_{primary}}'])*aircraft.fuse['x_{wing}']
                     ]), # TODO improve; using x_b as a surrogate for xeng
              ])
            if b737800:
                constraints.extend([
                TCS([climb['x_{CG}']*climb['W_{end}'] >=
                    0.5*(aircraft.fuse['W_{fuse}']+aircraft.fuse['W_{payload}'])*aircraft.fuse['l_{fuse}'] \
                    + (aircraft['W_{tail}'])*aircraft['x_{tail}'] \
                    + (aircraft['W_{wing_system}']+(climb['PCFuel']+aircraft['ReserveFraction'])*aircraft['W_{f_{primary}}'])*aircraft.fuse['x_{wing}'] \
                    + aircraft['numeng']*aircraft['W_{engsys}']*aircraft['x_b']]), # TODO improve; using x_b as a surrogate for xeng
##               cruise['x_{CG}'][0] <= climb['x_{CG}'],
                TCS([cruise['x_{CG}']*cruise['W_{end}'] >=
                    0.5*(aircraft.fuse['W_{fuse}']+aircraft.fuse['W_{payload}'])*aircraft.fuse['l_{fuse}'] \
                    + (aircraft['W_{tail}'])*aircraft['x_{tail}'] \
                    + (aircraft['W_{wing_system}']+(cruise['PCFuel']+aircraft['ReserveFraction'])*aircraft['W_{f_{primary}}'])*aircraft.fuse['x_{wing}'] \
                    + aircraft['numeng']*aircraft['W_{engsys}']*aircraft['x_b']]), # TODO improve; using x_b as a surrogate for xeng
              ])

        constraints.extend([
            climb.climbP.aircraftP['W_{start}'][0] == aircraft['W_{total}'],
            climb.climbP.aircraftP[
                'W_{end}'][-1] == cruise.cruiseP.aircraftP['W_{start}'][0],

            # similar constraint 1
            TCS([climb.climbP.aircraftP['W_{start}'] >= climb.climbP.aircraftP[
                'W_{end}'] + climb.climbP.aircraftP['W_{burn}']]),
            # similar constraint 2
            TCS([cruise.cruiseP.aircraftP['W_{start}'] >= cruise.cruiseP.aircraftP[
                'W_{end}'] + cruise.cruiseP.aircraftP['W_{burn}']]),

            climb.climbP.aircraftP['W_{start}'][
                1:] == climb.climbP.aircraftP['W_{end}'][:-1],
            cruise.cruiseP.aircraftP['W_{start}'][
                1:] == cruise.cruiseP.aircraftP['W_{end}'][:-1],

            TCS([aircraft['W_{dry}'] + aircraft['W_{payload}'] + \
                 aircraft['ReserveFraction']*aircraft['W_{f_{primary}}'] <= cruise.cruiseP.aircraftP['W_{end}'][-1]]),
            TCS([aircraft['W_{f_{climb}}'] >= sum(climb.climbP.aircraftP['W_{burn}'])]),
            TCS([aircraft['W_{f_{cruise}}'] >= sum(cruise.cruiseP.aircraftP['W_{burn}'])]),

            # Altitude constraints
            hftCruise >= CruiseAlt,
            TCS([hftClimb[1:Nclimb] >= hftClimb[:Nclimb-1] + dhft[1:Nclimb]]),
            TCS([hftClimb[0] == dhft[0]]),
            hftClimb[-1] <= hftCruise,

##            hftCruise <= 39692*units('ft'),

            hftCruise[0] == hftCruiseHold,

            # Compute dh
            dhft == hftCruiseHold / Nclimb,

            #compute fuel burn from TSFC
            cruise.cruiseP.aircraftP['W_{burn}'] == aircraft['numeng']*aircraft.engine['TSFC'][Nclimb:] * cruise['thr'] * aircraft.engine['F'][Nclimb:],              
            climb.climbP.aircraftP['W_{burn}'] == aircraft['numeng']*aircraft.engine['TSFC'][:Nclimb] * climb['thr'] * aircraft.engine['F'][:Nclimb],

            # Thrust constraint
            aircraft.VT['T_e'] == climb.climbP.engine['F'][0],

            # Set the range for each cruise segment, doesn't take credit for
            # down range distance covered during climb
            cruise.cruiseP['Rng'] == RngCruise / (Ncruise),

            # Cruise Mach Number constraint
            cruise['M'] >= aircraft['M_{min}'],

          #nacelle drag constraint
          #elevated this constarint to mission for dimensionality
          cruise.cruiseP['V_2'] == aircraft.engine['M_2'][Nclimb:]*cruise.state['a'],
          climb.climbP['V_2'] == aircraft.engine['M_2'][:Nclimb]*climb.state['a'],

          climb['\\alpha_{max,w}'] == .18,
          cruise['\\alpha_{max,w}'] == .1,


##          climb['\\alpha_{max,w}'] == .1,
##          cruise['\\alpha_{max,w}'] == .1,
##            aircraft['e'] <= 0.91,
        ])

        # Calculating percent fuel remaining
        with SignomialsEnabled():
            for i in range(0,Nclimb):
                constraints.extend([
                    TCS([climb['PCFuel'][i] >= (sum(climb['W_{burn}'][i+1:]) + \
                                                             aircraft['W_{f_{cruise}}'])/aircraft['W_{f_{primary}}']]) ,
                    climb['PCFuel'] <= 1.0,
                ])
            for i in range(0,Ncruise):
                constraints.extend([
                    TCS([cruise['PCFuel'][i] >= (sum(cruise['W_{burn}'][i+1:]) + \
                                0.0000001*aircraft['W_{f_{primary}}'])/aircraft['W_{f_{primary}}']]),
                    cruise['PCFuel'] <= 1.0,
                    ])

        if D80 or D82:
             constraints.extend([
                    # Set the BLI Benefit
                    climb.climbP.fuseP['f_{BLI}'] == 0.91, #TODO area for improvement
                    cruise.cruiseP.fuseP['f_{BLI}'] == 0.91,
                  ])
        if b737800:
            constraints.extend([
                climb.climbP.fuseP['f_{BLI}'] == 1.0,
                cruise.cruiseP.fuseP['f_{BLI}'] == 1.0,
            ])

        with SignomialsEnabled():
            constraints.extend([
                #set the range constraints
                TCS([sum(climb['RngClimb']) + RngCruise >= ReqRng]),
                ])
        
        M2 = .6
        M25 = .6
        M4a = .1025
        M0 = .5

        engineclimb = [
            aircraft.engine.engineP['M_2'][:Nclimb] == climb['M'],
            aircraft.engine.engineP['M_{2.5}'][:Nclimb] == M25,
            aircraft.engine.engineP['hold_{2}'][:Nclimb] == 1+.5*(1.398-1)*M2**2,
            aircraft.engine.engineP['hold_{2.5}'][:Nclimb] == 1+.5*(1.354-1)*M25**2,
            aircraft.engine.engineP['c1'][:Nclimb] == 1+.5*(.401)*M0**2,

            #constraint on drag and thrust
            aircraft['numeng']*aircraft.engine['F_{spec}'][:Nclimb] >= climb['D'] + climb['W_{avg}'] * climb['\\theta'],

            #climb rate constraints
            TCS([climb['excessP'] + climb.state['V'] * climb['D'] <=  climb.state['V'] * aircraft['numeng'] * aircraft.engine['F_{spec}'][:Nclimb]]),
            ]

        if D80 or D82:
             M2 = .6
             M25 = .6
             M4a = .1025
             M0 = .72
        if b737800:
             M2 = .6
             M25 = .6
             M4a = .1025
             M0 = .8


        enginecruise = [
            aircraft.engine.engineP['M_2'][Nclimb:] == cruise['M'],
            aircraft.engine.engineP['M_{2.5}'][Nclimb:] == M25,
            aircraft.engine.engineP['hold_{2}'][Nclimb:] == 1+.5*(1.398-1)*M2**2,
            aircraft.engine.engineP['hold_{2.5}'][Nclimb:] == 1+.5*(1.354-1)*M25**2,
            aircraft.engine.engineP['c1'][Nclimb:] == 1+.5*(.401)*M0**2,
            
            #steady level flight constraint on D 
            cruise['D'] == aircraft['numeng'] * aircraft.engine['F_{spec}'][Nclimb:],

            #breguet range eqn
            TCS([cruise['z_{bre}'] >= (aircraft.engine['TSFC'][Nclimb:] * cruise['thr']*
            cruise['D']) / cruise['W_{avg}']]),
            ]

        self.cost = aircraft['W_{f_{total}}']

        return constraints, aircraft, climb, cruise, enginestate, statelinking, engineclimb, enginecruise

M4a = .1025
fan = 1.60474
lpc  = 4.98
hpc = 35/8

substitutions = {
        # 'V_{stall}'   : 120,
        '\\delta_P_{over}': 12*units('psi'),
        'N_{land}': 6,
        'SPR': 8,
        'p_s': 81.*units('cm'),
        'ReqRng': 3000*units('nmi'),
        '\\theta_{db}' : 0.366,
##        'CruiseAlt': 36632*units('ft'),
        'numeng': 2,
        'numVT': 2,
        'numaisle':2,
        'n_{pax}': 180,
        'W_{avg. pass}': 180*units('lbf'),
        'W_{carry on}': 15*units('lbf'),
        'W_{cargo}': 10000*units('N'),
        'W_{checked}':40*units('lbf'),
        'W_{fix}': 3000*units('lbf'),
        'w_{aisle}': 0.51*units('m'),
        'w_{seat}': 0.5*units('m'),
        'w_{sys}': 0.1*units('m'),
        'r_E': 1,  # [TAS]
        '\\lambda_{cone}': 0.3,  # [TAS]
        '\\rho_{cone}': 2700,#*units('kg/m^3'),  # [TAS]
        '\\rho_{bend}': 2700,#*units('kg/m^3'),  # [TAS]
        '\\rho_{floor}': 2700,#*units('kg/m^3'),  # [TAS]
        '\\rho_{skin}': 2700,#*units('kg/m^3'),  # [TAS]
        '\\sigma_{floor}': 30000 / 0.000145, # [TAS] [Al]
        '\\sigma_{skin}': 15000 / 0.000145,  # [TAS] [Al]
        '\\tau_{floor}': 30000 / 0.000145, # [TAS] [Al]
        'W\'\'_{floor}': 60,  # [TAS]
        'W\'\'_{insul}': 22,  # [TAS]
        'W\'_{seat}': 150*units('N'),  # [TAS]
        'W\'_{window}': 145.*3.*units('N/m'),  # [TAS]

        # TASOPT Fuselage substitutions
        'l_{nose}': 29*0.3048*units('m'),

        # Fractional weights
        'f_{fadd}': 0.2,  # [TAS]
        'f_{frame}': 0.25,  # [Philippe]
        'f_{lugg,1}': 0.4,  # [Philippe]
        'f_{lugg,2}': 0.1,  # [Philippe]
        'f_{padd}': 0.4,  # [TAS]

        # Wing substitutions
        'C_{L_{wmax}}': 2.25, # [TAS]
        '\\tan(\\Lambda)': tan(sweep * pi / 180),
##        '\\alpha_{max,w}': 0.1,  # (6 deg)
        '\\cos(\\Lambda)': cos(sweep * pi / 180),
        '\\eta': 0.97,
        '\\rho_0': 1.225*units('kg/m^3'),
        '\\rho_{fuel}': 817*units('kg/m^3'),  # Kerosene [TASOPT]
        'FuelFrac': 0.9,
        'f_{flap}': 0.2,
        'f_{slat}': 0.0001,
        'f_{aileron}': 0.04,
        'f_{lete}': 0.1,
        'f_{ribs}': 0.15,
        'f_{spoiler}': 0.02,
        'f_{watt}': 0.03,

        # VT substitutions
       'C_{D_{wm}}': 0.5, # [2]
       'C_{L_{vmax}}': 2.6, # [TAS]
       'V_1': 70*units('m/s'),
       '\\rho_{TO}': 1.225*units('kg/m^3'),
        '\\tan(\\Lambda_{vt})': np.tan(40*np.pi/180),
        'c_{l_{vtEO}}': 0.5, # [TAS]
        'e_v': 0.8,
        # 'y_{eng}': 4.83*units('m'), # [3]
        'V_{land}': 72*units('m/s'),
        # 'I_{z}': 12495000, # estimate for late model 737 at max takeoff weight (m l^2/12)
        '\\dot{r}_{req}': 0.1475, # 10 deg/s/s yaw rate acceleration #NOTE: Constraint inactive
        'N_{spar}': 2,
        'f_{VT}': 0.4,

        # HT substitutions
        '\\alpha_{max,h}': 2.5,
        '\\tan(\\Lambda_{ht})': tan(30*pi/180),
        'C_{L_{hmax}}': 1.225,#2.0, # [TAS]
        'SM_{min}': 0.05,
        '\\Delta x_{CG}': 2.0*units('m'),
        'x_{CG_{min}}' : 10.0*units('m'),
        'C_{L_{hfcG}}': 0.85,
        'f_{HT}': 0.3,

        #engine system subs
        'rSnace': 6,
        'f_{pylon}': 0.05,
        'f_{eadd}': 0.1,

        #nacelle drag calc parameter
       'r_{vnace}': 0.925,

        # Cabin air substitutions in AircraftP

        #set the fuel reserve fraction
        'ReserveFraction': .20,

         # Engine substitutions
        '\\pi_{tn}': .995,
        '\pi_{b}': .94,
        '\pi_{d}': .995,
        '\pi_{fn}': .985,
        'T_{ref}': 288.15,
        'P_{ref}': 101.325,
        '\eta_{HPshaft}': .97,
        '\eta_{LPshaft}': .97,
        'eta_{B}': .9827,

        '\pi_{f_D}': fan,
        '\pi_{hc_D}': hpc,
        '\pi_{lc_D}': lpc,

        '\\alpha_{OD}': 6.97,
        '\\alpha_{max}': 6.97,

        'hold_{4a}': 1+.5*(1.313-1)*M4a**2,
        'r_{uc}': .01,
        '\\alpha_c': .19036,
        'T_{t_f}': 435,

        'M_{takeoff}': .9556,

        'G_f': 1,

        'h_f': 43.003,

        'Cp_t1': 1280,
        'Cp_t2': 1184,
        'Cp_c': 1216,

        'HTR_{f_SUB}': 1-.3**2,
        'HTR_{lpc_SUB}': 1 - 0.6**2,
}

def test():
     """
     solves a D82 and b737-800
     """
     global D82, b737800, sweep

     #run the D82 case
     D82 = True

     sweep = 13.237  # [deg]

     Nclimb = 5
     Ncruise = 2

     m = Mission(Nclimb, Ncruise)
     m.substitutions.update(substitutions)

     sweep = 13.237
     m.substitutions.update({
          #Fuselage subs
          'f_{seat}':0.1,
          'W\'_{seat}':1, # Seat weight determined by weight fraction instead
          'f_{string}':0.35,
          'AR':15.749,
          'h_{floor}': 0.13,
          'R_{fuse}' : 1.715,
          '\\delta R_{fuse}': 0.43,
          'w_{db}': 0.93,
          'b_{max}':140.0*0.3048,
          # 'c_0': 17.4*0.3048,#units('ft'),
          '\\delta_P_{over}': 8.382*units('psi'),

          #HT subs
          'AR_h': 12.,
          '\\lambda_h' : 0.3,
          '\\tan(\\Lambda_{ht})': np.tan(8*np.pi/180), # tangent of HT sweep
          # 'V_{h}': 0.895,

          #VT subs
          'A_{vt}' : 2.2,
          '\\lambda_{vt}': 0.3,
          '\\tan(\\Lambda_{vt})': np.tan(25*np.pi/180), # tangent of VT sweep
          ##                'V_{vt}': .03,

          #Wing subs
          'C_{L_{wmax}}': 2.15,

          # Minimum Cruise Mach Number
          'M_{min}': 0.72,
     })
     m.substitutions.__delitem__('\\theta_{db}')

     m = relaxed_constants(m)

     sol = m.localsolve( verbosity = 4, iteration_limit=50)

if __name__ == '__main__':
    Nclimb = 3
    Ncruise = 2
    
    if sweeps == False: 
        m = Mission(Nclimb, Ncruise)
        m.substitutions.update(substitutions)
        # m = Model(m.cost,BCS(m))
        if D80:
            print('D80 executing...')
            # for constraint in m.flat(constraintsets=False):
            #         if 'l_{nose}' in constraint.varkeys:
            #             print constraint
            sweep = 27.566
            m.substitutions.update({
                #Fuselage subs
                'f_{seat}':0.1,
                'W\'_{seat}':1, # Seat weight determined by weight fraction instead
                'f_{string}':0.35,
                'AR':10.8730,
                'h_{floor}': 0.13,
                'R_{fuse}' : 1.715,
                'w_{db}': 0.93,
                # 'b':116.548*0.3048,#units('ft'),
                # 'c_0': 17.4*0.3048,#units('ft'),
                #HT subs
                'AR_h': 8.25,
                '\\lambda_h' : 0.25,
                '\\tan(\\Lambda_{ht})':np.tan(20*np.pi/180), #tangent of HT sweep

                #VT subs
                'numVT': 2,
                'A_{vt}' : 2.0,
                '\\lambda_{vt}': 0.3,
                '\\tan(\\Lambda_{vt})':np.tan(25*np.pi/180),

                # Minimum Cruise Mach Number
                'M_{min}': 0.8,
            })
            m.substitutions.__delitem__('\\theta_{db}')

        if D82:
            print('D82 executing...')
            sweep = 13.237
            m.substitutions.update({
                #Fuselage subs
                'f_{seat}':0.1,
                'W\'_{seat}':1, # Seat weight determined by weight fraction instead
                'f_{string}':0.35,
                # 'AR':15.749,
                'h_{floor}': 0.13,
                'R_{fuse}' : 1.715,
                '\\delta R_{fuse}': 0.43,
                'w_{db}': 0.93,
                'b_{max}':140.0*0.3048,
                # 'c_0': 17.4*0.3048,#units('ft'),
                '\\delta_P_{over}': 8.382*units('psi'),

                #HT subs
                'AR_h': 12.,
                '\\lambda_h' : 0.3,
                '\\tan(\\Lambda_{ht})': np.tan(8*np.pi/180), # tangent of HT sweep
                # 'V_{h}': 0.895,

                #VT subs
                'numVT': 2,
                # 'A_{vt}' : 2.2,
                '\\lambda_{vt}': 0.3,
                '\\tan(\\Lambda_{vt})': np.tan(25*np.pi/180), # tangent of VT sweep
##                'V_{vt}': .03,

                #Wing subs
                'C_{L_{wmax}}': 2.15,

                # Minimum Cruise Mach Number
                'M_{min}': 0.72,
            })
            m.substitutions.__delitem__('\\theta_{db}')

        if b737800:
               print('737-800 executing...')
               
               M4a = .1025
               fan = 1.685
               lpc  = 8/1.685
               hpc = 30/8
               
               m.substitutions.update({
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

                  '\\alpha_{OD}': 5.1,
                  '\\alpha_{max}': 5.1,

                  'hold_{4a}': 1+.5*(1.313-1)*M4a**2,
                  'r_{uc}': .01,
                  '\\alpha_c': .19036,
                  'T_{t_f}': 435,

                  'M_{takeoff}': .9556,

                  'G_f': 1,

                  'h_f': 43.003,

                  'Cp_t1': 1280,
                  'Cp_t2': 1184,
                  'Cp_c': 1216,

                  'HTR_{f_SUB}': 1-.3**2,
                  'HTR_{lpc_SUB}': 1 - 0.6**2,

                  #fuselage subs that make fuse circular
                  '\\delta R_{fuse}': 0.0001*units('m'),
                  '\\theta_{db}': 0.0001,

                    #Fuselage subs
                   'numaisle':1.,
                   'SPR':6,
                    'f_{seat}':0.1,
                    'W\'_{seat}':1*units('N'), # Seat weight determined by weight fraction instead
                    'f_{string}':0.35,
                    'h_{floor}': 0.13*units('m'),
                    # 'R_{fuse}' : 1.715,
                    'b_{max}':117.5*units('ft'),
                    # 'c_0': 17.4*0.3048,#units('ft'),
                    '\\delta_P_{over}': 8.382*units('psi'),

                    #HT subs
                    'AR_h': 6,
                    '\\lambda_h' : 0.25,
                    '\\tan(\\Lambda_{ht})': np.tan(25*np.pi/180), # tangent of HT sweep
##                    'V_{h}': .6,
                    'C_{L_{hmax}}': 2.0, # [TAS]
                    'C_{L_{hfcG}}': 0.7,
                    '\\Delta x_{CG}': 7.68*units('ft'),
                    'x_{CG_{min}}' : 56.75*units('ft'),
                    'SM_{min}': .05,

                    #VT subs
                    'numVT': 1,
                    'A_{vt}' : 2,
                    '\\lambda_{vt}': 0.3,
                    '\\tan(\\Lambda_{vt})': np.tan(25*np.pi/180), # tangent of VT sweep
##                    'V_{vt}': .07,
                  'N_{spar}': 1,

                    #Wing subs
                    'C_{L_{wmax}}': 2.15,
                   'f_{slat}': 0.1,
                   'AR':10.1,

##                   'e': .91,

                    # Minimum Cruise Mach Number
                    'M_{min}': 0.8,
                   # 'CruiseAlt':38000*units('ft'),

                    #engine system subs
                    'rSnace': 16,

                    #nacelle drag calc parameter
                    'r_{vnace}': 1.02,
                })

        if D82:
             m_relax = relaxed_constants(m)
        if b737800:
             m_relax = relaxed_constants(m, None, ['M_{takeoff}', '\\theta_{db}'])

        sol = m_relax.localsolve( verbosity = 4, iteration_limit=200)

        post_process(sol)

##        m.cost = m_relax.cost
##
##        sol = m.localsolve( verbosity = 4, iteration_limit=50, x0=sol['variables'])

        if D82:
             percent_diff(sol, 2, Nclimb)

        if b737800:
             percent_diff(sol, 801, Nclimb)

    if sweeps:
        if sweepSMmin:
            m = Mission(Nclimb, Ncruise)
            m.substitutions.update(substitutions)
            SMminArray = np.linspace(0.05,0.5,n)
            m.substitutions.update({'SM_{min}': ('sweep',SMminArray)})
            m = relaxed_constants(m)
            solSMsweep = m.localsolve(verbosity = 4, skipsweepfailures=True)

            if plot:
                plt.plot(solSMsweep('SM_{min}'), solSMsweep('S_h'), '-r')
                plt.xlabel('Minimum Allowed Static Margin')
                plt.ylabel('Horizontal Tail Area [m$^2$]')
                plt.title('Horizontal Tail Area vs Min Static Margin')
                plt.savefig('CFP_Sweeps/S_{h}-vs-SM_{min}.pdf')
                plt.show(), plt.close()

                plt.plot(solSMsweep('SM_{min}'), solSMsweep('V_{h}'), '-r')
                plt.xlabel('Minimum Allowed Static Margin')
                plt.ylabel('Horizontal Tail Volume Coefficient')
                plt.title('Horizontal Tail Volume Coefficient vs Min Static Margin')
                plt.savefig('CFP_Sweeps/V_{h}-vs-SM_{min}.pdf')
                plt.show(), plt.close()

##                plt.plot(solSMsweep('SM_{min}'), np.mean(solSMsweep('x_{CG}_Mission, CruiseSegment, CruiseP, AircraftP'),axis = 1), '-r')
##                plt.xlabel('Minimum Allowed Static Margin')
##                plt.ylabel('CG Location [m]')
##                plt.title('CG Location vs Min Static Margin')
##                plt.savefig('CFP_Sweeps/x_{CG}-vs-SM_{min}.pdf')
##                plt.show(), plt.close()

        if sweepReqRng:
            m = Mission(Nclimb, Ncruise)
            m.substitutions.update(substitutions)
            ReqRngArray = np.linspace(500,3000,n)
            m.substitutions.update({'ReqRng': ('sweep',ReqRngArray)})
            m = relaxed_constants(m)
            solReqRngsweep = m.localsolve(verbosity=2, iteration_limit=25, skipsweepfailures=True)

            if plot:
                plt.plot(solReqRngsweep('ReqRng'),solReqRngsweep('W_{f_{total}}'))
                plt.xlabel('Range [nmi]')
                plt.ylabel('Total Fuel Weight [lbf]')
                plt.title('Total Fuel Weight vs Range')
                plt.savefig('CFP_Sweeps/W_ftotal-vs-ReqRng.pdf')
                plt.show(), plt.close()

                plt.plot(solReqRngsweep('ReqRng'),solReqRngsweep('W_{dry}'))
                plt.xlabel('Range [nmi]')
                plt.ylabel('Dry Weight [lbf]')
                plt.title('Dry Weight vs Range')
                plt.savefig('CFP_Sweeps/W_dry-vs-ReqRng.pdf')
                plt.show(), plt.close()

                plt.plot(solReqRngsweep('ReqRng'),solReqRngsweep('S'))
                plt.xlabel('Range [nmi]')
                plt.ylabel('Wing Area [m^2]')
                plt.title('Wing Area vs Range')
                plt.savefig('CFP_Sweeps/S-vs-ReqRng.pdf')
                plt.show(), plt.close()

                plt.plot(solReqRngsweep('ReqRng'),solReqRngsweep('b'))
                plt.xlabel('Range [nmi]')
                plt.ylabel('Wing Span [m]')
                plt.title('Wing Span vs Range')
                plt.savefig('CFP_Sweeps/b-vs-ReqRng.pdf')
                plt.show(), plt.close()

        if sweepthetadb:
            m = Mission(Nclimb, Ncruise)
            m.substitutions.update(substitutions)
            thetadbArray = np.linspace(0,0.5,n)
            m.substitutions.update({'\\theta_{db}': ('sweep', thetadbArray)})
            m = relaxed_constants(m)
            solthetadbsweep = m.localsolve(verbosity=2, iteration_limit=25, skipsweepfailures=True)

            if plot:
                plt.plot(solthetadbsweep('\\theta_{db}'),solthetadbsweep('W_{f_{total}}'))
                plt.xlabel('Fuselage Joint Angle [radians]')
                plt.ylabel('Total Fuel Weight [lbf]')
                plt.title('Total Fuel Weight vs Fuselage Joint Angle')
                plt.savefig('CFP_Sweeps/W_ftotal-vs-thetadb.pdf')
                plt.show(), plt.close()

                plt.plot(solthetadbsweep('\\theta_{db}'),solthetadbsweep('W_{dry}'))
                plt.xlabel('Fuselage Joint Angle [radians]')
                plt.ylabel('Dry Weight [lbf]')
                plt.title('Dry Weight vs Fuselage Joint Angle')
                plt.savefig('CFP_Sweeps/W_dry-vs-thetadb.pdf')
                plt.show(), plt.close()

                plt.plot(solthetadbsweep('\\theta_{db}'),solthetadbsweep('W_{fuse}'))
                plt.xlabel('Fuselage Joint Angle [radians]')
                plt.ylabel('Fuselage Weight [lbf]')

                plt.title('Fuselage Weight vs Fuselage Joint Angle')
                plt.savefig('CFP_Sweeps/W_fuse-vs-thetadb.pdf')
                plt.show(), plt.close()

        if sweepdxCG:
            m = Mission(Nclimb, Ncruise)
            m.substitutions.update(substitutions)
            dxCGArray = np.linspace(0.5,3.5,n)
            m.substitutions.update({'\\Delta x_{CG}': ('sweep',dxCGArray)})
            m = relaxed_constants(m)
            soldxCGsweep = m.localsolve(verbosity=2,skipsweepfailures=True)

            if plot:
                plt.plot(soldxCGsweep('\\Delta x_{CG}'),soldxCGsweep('V_{h}'),'-r')
                plt.xlabel('Allowed CG shift [m]')
                plt.ylabel('Horizontal Tail Volume Coefficient')
                plt.title('Horizontal Tail Volume Coefficient vs Allowed CG Shift')
                plt.savefig('CFP_Sweeps/Vht-vs-dxCG.pdf')
                plt.show(), plt.close()

                plt.plot(soldxCGsweep('\\Delta x_{CG}'),np.mean(soldxCGsweep('x_{CG}_Mission, CruiseSegment, CruiseP, AircraftP'),axis = 1), '-r')
                plt.xlabel('Allowed CG shift [m]')
                plt.ylabel('Nominal CG location [m]')
                plt.title('CG Location vs Allowed CG Shift')
                plt.savefig('CFP_Sweeps/xCG-vs-dxCG.pdf')

                plt.plot(soldxCGsweep('\\Delta x_{CG}'),np.mean(soldxCGsweep('x_{AC}_Mission, CruiseSegment, CruiseP, AircraftP'),axis = 1), '-r')
                plt.xlabel('Allowed CG shift [m]')
                plt.ylabel('Average AC location [m]')
                plt.title('AC Location vs Allowed CG Shift')
                plt.savefig('CFP_Sweeps/xAC-vs-dxCG.pdf')

        if sweepxCG:
            m = Mission(Nclimb, Ncruise)
            m.substitutions.update(substitutions)
            xCGArray = np.linspace(8,14,n)
            m.substitutions.update({'x_{CG_{min}}': ('sweep',xCGArray)})
            m = relaxed_constants(m)
            solxCGsweep = m.localsolve(verbosity=2,skipsweepfailures=True,iteration_limit=30)

            if plot:
                plt.plot(solxCGsweep('x_{CG_{min}}'),solxCGsweep('V_{h}'),'-r')
                plt.xlabel('Max Forward CG [m]')
                plt.ylabel('Horizontal Tail Volume Coefficient')
                plt.title('Horizontal Tail Volume Coefficient vs Max Forward CG')
                plt.savefig('CFP_Sweeps/Vht-vs-xCG.pdf')
                plt.show(), plt.close()

                plt.plot(solxCGsweep('x_{CG_{min}}'),np.mean(solxCGsweep('x_{CG}_Mission, CruiseSegment, CruiseP, AircraftP'),axis = 1), '-r')
                plt.xlabel('Max Forward CG [m]')
                plt.ylabel('Nominal CG location [m]')
                plt.title('CG Location vs Max Forward CG')
                plt.savefig('CFP_Sweeps/xCG-vs-xCG.pdf')
                plt.show(), plt.close()

                plt.plot(solxCGsweep('x_{CG_{min}}'),np.mean(solxCGsweep('x_{AC}_Mission, CruiseSegment, CruiseP, AircraftP'),axis = 1), '-r')
                plt.xlabel('Max Forward CG [m]')
                plt.ylabel('Average AC location [m]')
                plt.title('AC Location vs Max Forward CG')
                plt.savefig('CFP_Sweeps/xAC-vs-xCG.pdf')
                plt.show(), plt.close()

        if sweepCruiseAlt:
            m = Mission(Nclimb, Ncruise)
            m.substitutions.update(substitutions)
            CruiseAltArray = np.linspace(25000,40000,n)
            m.substitutions.update({'CruiseAlt': ('sweep',CruiseAltArray)})
            m = relaxed_constants(m)
            solCruiseAltsweep = m.localsolve(verbosity=2,skipsweepfailures=True,iteration_limit=30)

            if plot:
                plt.plot(solCruiseAltsweep('CruiseAlt'),solCruiseAltsweep('W_{f_{total}}'))
                plt.xlabel('Cruise Altitude [ft]')
                plt.ylabel('Mission Fuel Burn [lbs]')
                plt.title('Fuel Burn vs Cruise Altitude')
                plt.savefig('CFP_Sweeps/Wftotal-vs-CruiseAlt.pdf')
                plt.show(), plt.close()

                plt.plot(solCruiseAltsweep('CruiseAlt'),solCruiseAltsweep('AR'))
                plt.xlabel('Cruise Altitude [ft]')
                plt.ylabel('Aspect Ratio')
                plt.title('Aspect Ratio vs Cruise Altitude')
                plt.savefig('CFP_Sweeps/AR-vs-CruiseAlt.pdf')
                plt.show(), plt.close()

                plt.plot(solCruiseAltsweep('CruiseAlt'),solCruiseAltsweep('S'))
                plt.xlabel('Cruise Altitude [ft]')
                plt.ylabel('Wing Area [ft^2]')
                plt.title('Wing Area vs Cruise Altitude')
                plt.savefig('CFP_Sweeps/S-vs-CruiseAlt.pdf')
                plt.show(), plt.close()

                plt.plot(solCruiseAltsweep('CruiseAlt'),np.mean(solCruiseAltsweep('M_Mission, CruiseSegment, FlightState'),axis = 1))
                plt.xlabel('Cruise Altitude [m]')
                plt.ylabel('Mach Number')
                plt.title('Average Mach Number vs Cruise Altitude')
                plt.savefig('CFP_Sweeps/M-vs-CruiseAlt.pdf')
                plt.show(), plt.close()

        if sweepMmin:
            m = Mission(Nclimb, Ncruise)
            m.substitutions.update(substitutions)
            MminArray = np.linspace(0.4,0.9,n) # No lower limit to high Mach
            m.substitutions.update({'M_{min}':('sweep',MminArray)})
            m = relaxed_constants(m)
            solMminsweep = m.localsolve(verbosity=2,skipsweepfailures=True,iteration_limit=30)

            if plot:
                plt.plot(solMminsweep('M_{min}'),solMminsweep('W_{f_{total}}'))
                plt.xlabel('Minimum Cruise Mach Number')
                plt.ylabel('Total Fuel Weight [lbf]')
                plt.title('Total Fuel Weight vs Minimum Cruise Mach Number')
                plt.savefig('CFP_Sweeps/Wftotal-vs-Mmin.pdf')
                plt.show(), plt.close()

                plt.plot(solMminsweep('M_{min}'),solMminsweep('M_Mission, CruiseSegment, FlightState')[:,0])
                plt.xlabel('Minimum Cruise Mach Number')
                plt.ylabel('TOC Mach Number')
                plt.title('TOC vs Minimum Cruise Mach Number')
                plt.savefig('CFP_Sweeps/M-vs-Mmin.pdf')
                plt.show(), plt.close()

                plt.plot(solMminsweep('M_{min}'),solMminsweep('W_{dry}'))
                plt.xlabel('Minimum Cruise Mach Number')
                plt.ylabel('Aircraft Dry Weight [lbf]')
                plt.title('Aircraft Dry Weight vs Minimum Cruise Mach Number')
                plt.savefig('CFP_Sweeps/W_dry-vs-Mmin.pdf')
                plt.show(), plt.close()

                plt.plot(solMminsweep('M_{min}'),solMminsweep('b'))
                plt.xlabel('Minimum Cruise Mach Number')
                plt.ylabel('Wingspan [m]')
                plt.title('Wing Span vs Minimum Cruise Mach Number')
                plt.savefig('CFP_Sweeps/b-vs-Mmin.pdf')
                plt.show(), plt.close()

                plt.plot(solMminsweep('M_{min}'),solMminsweep('AR'))
                plt.xlabel('Minimum Cruise Mach Number')
                plt.ylabel('Aspect Ratio')
                plt.title('Aspect Ratio vs Minimum Cruise Mach Number')
                plt.savefig('CFP_Sweeps/AR-vs-Mmin.pdf')
                plt.show(),plt.close()

                plt.plot(solMminsweep('M_{min}'),solMminsweep('S'))
                plt.xlabel('Minimum Cruise Mach Number')
                plt.ylabel('Wing Area [m^2]')
                plt.title('Wing Are vs. Minimum Cruise Mach Number')
                plt.savefig('CFP_Sweeps/S-vs-Mmin.pdf')
                plt.show(),plt.close()

##                plt.plot(solMminsweep('M_{min}'),solMminsweep('CruiseAlt'))
##                plt.xlabel('Minimum Cruise Mach Number')
##                plt.ylabel('Cruise Altitude')
##                plt.title('Cruise Altitude vs. Minimum Cruise Mach Number')
##                plt.savefig('CFP_Sweeps/CruiseAlt-vs-Mmin.pdf')
##                plt.show(),plt.close()

        if sweepnpax:
            m = Mission(Nclimb, Ncruise)
            m.substitutions.update(substitutions)
            npaxArray = np.linspace(150,400,n)
            m.substitutions.update({'n_{pax}':('sweep',npaxArray)})
            m = relaxed_constants(m)
            solnpaxsweep = m.localsolve(verbosity=2,skipsweepfailures=True,iteration_limit=30)

            if plot:
                plt.plot(solnpaxsweep('n_{pax}'),solnpaxsweep('W_{total}')),
                plt.xlabel('Number of Passengers')
                plt.ylabel('Takeoff Weight [lbf]')
                plt.title('Takeoff Weight vs Number of Passengers')
                plt.savefig('CFP_Sweeps/Wtotal-vs-npax.pdf')
                plt.show(),plt.close()

                plt.plot(solnpaxsweep('n_{pax}'),solnpaxsweep('W_{f_{total}}'))
                plt.xlabel('Number of Passengers')
                plt.ylabel('Total Fuel Weight [lbf]')
                plt.title('Total Fuel Weight vs Number of Passengers')
                plt.savefig('CFP_Sweeps/Wftotal-vs-npax.pdf')
                plt.show(),plt.close()

                plt.plot(solnpaxsweep('n_{pax}'),solnpaxsweep('AR'))
                plt.xlabel('Number of Passengers')
                plt.ylabel('Aspect Ratio')
                plt.title('Aspect Ratio vs Number of Passengers')
                plt.savefig('CFP_Sweeps/AR-vs-npax.pdf')
                plt.show(),plt.close()

                plt.plot(solnpaxsweep('n_{pax}'),solnpaxsweep('CruiseAlt'))
                plt.xlabel('Number of Passengers')
                plt.ylabel('Cruise Altitude [ft]')
                plt.title('Cruise Altitude vs Number of Passengers')
                plt.savefig('CFP_Sweeps/CruiseAlt-vs-npax.pdf')
                plt.show(),plt.close()

                plt.plot(solnpaxsweep('n_{pax}'),solnpaxsweep('S'))
                plt.xlabel('Number of Passengers')
                plt.ylabel('Wing Area [m^2]')
                plt.title('Wing Area vs Number of Passengers')
                plt.savefig('CFP_Sweeps/S-vs-npax.pdf')
                plt.show(),plt.close()

                plt.plot(solnpaxsweep('n_{pax}'),solnpaxsweep('W_{fuse}'))
                plt.xlabel('Number of Passengers')
                plt.ylabel('Fuselage Weight [lbf]')
                plt.title('Fuselage Weight vs Number of Passengers')
                plt.savefig('CFP_Sweeps/Wfuse-vs-npax.pdf')
                plt.show(),plt.close()

                plt.plot(solnpaxsweep('n_{pax}'),solnpaxsweep('W_{f_{total}}')/solnpaxsweep('W_{total}'))
                plt.xlabel('Number of Passengers')
                plt.ylabel('Fuel Mass Fraction')
                plt.title('Fuel Mass Fraction vs Number of Passengers')
                plt.savefig('CFP_Sweeps/ffuel-vs-npax.pdf')
                plt.show(),plt.close()

                plt.plot(solnpaxsweep('n_{pax}'),solnpaxsweep('W_{payload}')/solnpaxsweep('W_{total}'))
                plt.xlabel('Number of Passengers')
                plt.ylabel('Payload Mass Fraction')
                plt.title('Payload Mass Fraction vs Number of Passengers')
                plt.savefig('CFP_Sweeps/fpayload-vs-npax.pdf')
                plt.show(),plt.close()

                plt.plot(solnpaxsweep('n_{pax}'),solnpaxsweep('W_{wing}')/solnpaxsweep('W_{total}'))
                plt.xlabel('Number of Passengers')
                plt.ylabel('Wing Mass Fraction')
                plt.title('Wing Mass Fraction vs Number of Passengers')
                plt.savefig('CFP_Sweeps/fwing-vs-npax.pdf')
                plt.show(),plt.close()

        if sweepResFuel:
            m = Mission(Nclimb, Ncruise)
            m.substitutions.update(substitutions)
            ResFuelArray = np.linspace(.05,.25,n)
            m.substitutions.update({'ReserveFraction':('sweep',ResFuelArray)})
            m = relaxed_constants(m)
            solResFuelsweep = m.localsolve(verbosity=2,skipsweepfailures=True,iteration_limit=30)
            
            if plot:
                plt.plot(solResFuelSweep('n_{pax}'),solResFuelSweep('W_{total}')),
                plt.xlabel('Reserve Fuel Fraction')
                plt.ylabel('Takeoff Weight [lbf]')
                plt.title('Takeoff Weight vs Reserve Fuel Fraction')
                plt.savefig('CFP_Sweeps/Wtotal-vs-resfuel.pdf')
                plt.show(),plt.close()

                plt.plot(solResFuelSweep('n_{pax}'),solResFuelSweep('W_{f_{total}}'))
                plt.xlabel('Reserve Fuel Fraction')
                plt.ylabel('Total Fuel Weight [lbf]')
                plt.title('Total Fuel Weight vs Reserve Fuel Fraction')
                plt.savefig('CFP_Sweeps/Wftotal-vs-resfuel.pdf')
                plt.show(),plt.close()

                plt.plot(solResFuelSweep('n_{pax}'),solResFuelSweep('AR'))
                plt.xlabel('Reserve Fuel Fraction')
                plt.ylabel('Aspect Ratio')
                plt.title('Aspect Ratio vs Reserve Fuel Fraction')
                plt.savefig('CFP_Sweeps/AR-vs-resfuel.pdf')
                plt.show(),plt.close()

                plt.plot(solResFuelSweep('n_{pax}'),solResFuelSweep('CruiseAlt'))
                plt.xlabel('Reserve Fuel Fraction')
                plt.ylabel('Cruise Altitude [ft]')
                plt.title('Cruise Altitude vs Reserve Fuel Fraction')
                plt.savefig('CFP_Sweeps/CruiseAlt-vs-resfuel.pdf')
                plt.show(),plt.close()

                plt.plot(solResFuelSweep('n_{pax}'),solResFuelSweep('S'))
                plt.xlabel('Reserve Fuel Fraction')
                plt.ylabel('Wing Area [m^2]')
                plt.title('Wing Area vs Reserve Fuel Fraction')
                plt.savefig('CFP_Sweeps/S-vs-resfuel.pdf')
                plt.show(),plt.close()

                plt.plot(solResFuelSweep('n_{pax}'),solResFuelSweep('W_{fuse}'))
                plt.xlabel('Reserve Fuel Fraction')
                plt.ylabel('Fuselage Weight [lbf]')
                plt.title('Fuselage Weight vs Reserve Fuel Fraction')
                plt.savefig('CFP_Sweeps/Wfuse-vs-resfuel.pdf')
                plt.show(),plt.close()

                plt.plot(solResFuelSweep('n_{pax}'),solResFuelSweep('W_{f_{total}}')/solResFuelSweep('W_{total}'))
                plt.xlabel('Reserve Fuel Fraction')
                plt.ylabel('Fuel Mass Fraction')
                plt.title('Fuel Mass Fraction vs Reserve Fuel Fraction')
                plt.savefig('CFP_Sweeps/ffuel-vs-resfuel.pdf')
                plt.show(),plt.close()

                plt.plot(solResFuelSweep('n_{pax}'),solResFuelSweep('W_{payload}')/solResFuelSweep('W_{total}'))
                plt.xlabel('Reserve Fuel Fraction')
                plt.ylabel('Payload Mass Fraction')
                plt.title('Payload Mass Fraction vs Reserve Fuel Fraction')
                plt.savefig('CFP_Sweeps/fpayload-vs-resfuel.pdf')
                plt.show(),plt.close()

                plt.plot(solResFuelSweep('n_{pax}'),solResFuelSweep('W_{wing}')/solResFuelSweep('W_{total}'))
                plt.xlabel('Reserve Fuel Fraction')
                plt.ylabel('Wing Mass Fraction')
                plt.title('Wing Mass Fraction vs Reserve Fuel Fraction')
                plt.savefig('CFP_Sweeps/fwing-vs-resfuel.pdf')
                plt.show(),plt.close()
