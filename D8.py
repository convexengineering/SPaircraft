"""Simple commercial aircraft flight profile and D8 aircraft model"""
""" Integrates Wing, VerticalTail, HorizontalTail and Fuselage models """

import numpy as np
from gpkit import Variable, Model, units, SignomialsEnabled, SignomialEquality, Vectorize
from gpkit.constraints.tight import Tight as TCS
from gpkit.constraints.bounded import Bounded as BCS
from gpkit.tools import te_exp_minus1
from numpy import pi

TCS.reltol = 1e-3

# importing from D8_integration
from stand_alone_simple_profile import FlightState
from D8_VT_yaw_rate_and_EO_simple_profile import VerticalTail
from D8_HT_simple_profile import HorizontalTail
from D8_Wing_simple_profile import Wing
from turbofan.engine_validation import Engine
from D8_Fuselage import Fuselage


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

g = 9.81 * units('m*s**-2')

class Aircraft(Model):
    """
    Aircraft class

    ARGUMENTS
    ---------
    BLI: True = have engine stagnation pressure drop drom BLI, False = no engine stagnation pressure drop
    fitDrag: True = use Martin's tail drag fits, False = use the TASOPT tail drag model
    """

    def setup(self, Nclimb, Ncruise, enginestate, eng, fitDrag, BLI = False, Nmissions=0,  **kwargs):
        # create submodels
        self.fuse = Fuselage(Nmissions)
        self.wing = Wing()
        if Nmissions != 0:
            self.engine = Engine(0, True, Nclimb+Ncruise, enginestate, eng, Nmissions, BLI)
        else:
           self.engine = Engine(0, True, Nclimb+Ncruise, enginestate, eng, BLI)
        self.VT = VerticalTail()
        self.HT = HorizontalTail()

        #set the tail drag flag
        self.fitDrag = fitDrag

        # variable definitions
        numaisle = Variable('n_{aisle}','-','Number of Aisles')
        numeng = Variable('n_{eng}', '-', 'Number of Engines')
        numVT = Variable('n_{VT}','-','Number of Vertical Tails')
        Vne = Variable('V_{ne}', 'm/s', 'Never-exceed speed')  # [Philippe]
        Vmn = Variable('V_{mn}', 'm/s','Maneuvering speed')
        rhoTO = Variable('\\rho_{T/O}',1.225,'kg*m^-3','Air density at takeoff')
        ReserveFraction = Variable('ReserveFraction', '-', 'Fuel Reserve Fraction')

        #fraction of fuel in wings
        f_wingfuel = Variable('f_{wingfuel}', '-', 'Fraction of fuel stored in wing tanks')

        SMmin = Variable('SM_{min}','-', 'Minimum Static Margin')
        dxCG = Variable('\\Delta x_{CG}', 'm', 'Max CG Travel Range')
        xCGmin = Variable('x_{CG_{min}}','m','Maximum Forward CG')

        with Vectorize(Nmissions):
             Izwing = Variable('I_{z_{wing}}','kg*m**2','Wing moment of inertia')
             Iztail = Variable('I_{z_{tail}}','kg*m**2','Tail moment of inertia')
             Izfuse = Variable('I_{z_{fuse}}','kg*m**2','Fuselage moment of inertia')

        Mmin = Variable('M_{min}','-','Minimum Cruise Mach Number')

        # Weights
        with Vectorize(Nmissions):
             PRFC = Variable('PRFC','','Payload-Range Fuel Consumption')
             W_total = Variable('W_{total}', 'lbf', 'Total Aircraft Weight')
             W_dry = Variable('W_{dry}', 'lbf', 'Zero Fuel Aircraft Weight')
             W_ftotal = Variable('W_{f_{total}}', 'lbf', 'Total Fuel Weight')
             W_fclimb = Variable('W_{f_{climb}}', 'lbf','Fuel Weight Burned in Climb')
             W_fcruise = Variable('W_{f_{cruise}}', 'lbf','Fuel Weight Burned in Cruise')
             W_fprimary = Variable('W_{f_{primary}}', 'lbf', 'Total Fuel Weight Less Fuel Reserves')

        Wwing = Variable('W_{wing}','lbf','Wing Weight')
        WHT = Variable('W_{HT}','lbf','Horizontal Tail Weight')
        WVT = Variable('W_{VT}','lbf','Vertical Tail Weight')

        # Fuselage lift fraction variables
        Ltow = Variable('L_{total/wing}','-','Total lift as a percentage of wing lift')

        # Misc system variables
        Wmisc   = Variable('W_{misc}','lbf','Sum of Miscellaneous Weights')
        Wlgnose = Variable('W_{lgnose}','lbf','Nose Landing Gear Weight')
        Wlgmain = Variable('W_{lgmain}','lbf','Main Landing Gear Weight')
        Wlg = Variable('W_{lg}', 'lbf', 'Total Landing Gear Weight')
        Clg = Variable('C_{lg}', 1, '-', 'Landing Gear Weight Margin/Sens Factor')
        Whpesys = Variable('W_{hpesys}','lbf','Power Systems Weight')
        #
        flgnose = Variable('f_{lgnose}','-','Nose Landing Gear Weight Fraction')
        flgmain = Variable('f_{lgmain}','-','Main Landing Gear Weight Fraction')
        fhpesys = Variable('f_{hpesys}','-','Power Systems Weight Fraction')
        #
        xmisc   = Variable('x_{misc}','m','Misc Weight Centroid')
        xlgnose = Variable('x_{lgnose}','m','Nose Landing Gear Weight x-Location')
        xlgmain = Variable('x_{lgmain}','m','Main Landing Gear Weight x-Location')
        xhpesys = Variable('x_{hpesys}','m','Power Systems Weight x-Location')

        #engine system weight variables
        rSnace = Variable('r_{S_{nacelle}}', '-', 'Nacelle and Pylon Wetted Area')
        lnace = Variable('l_{nacelle}', 'm', 'Nacelle Length')
        fSnace = Variable('f_{S_{nacelle}}', '-', 'Non-dimensional Nacelle Area')
        Snace = Variable('S_{nacelle}', 'm^2', 'Nacelle Surface Area')
        Ainlet = Variable('A_{inlet}','m^2', 'Inlet Area')
        Afancowl = Variable('A_{fancowl}', 'm^2', 'Fan Cowling Area')
        Aexh = Variable('A_{exh}', 'm^2', 'Exhaust Area')
        Acorecowl = Variable('A_{corecowl}', 'm^2', 'Core Cowling Area')
        Wnace = Variable('W_{nacelle}', 'lbf', 'Nacelle Weight')
        Wpylon = Variable('W_{pylon}', 'lbf','Engine Pylon Weight')
        fpylon = Variable('f_{pylon}', '-', 'Pylong Weight Fraction')
        feadd = Variable('f_{eadd}', '-', 'Additional Engine Weight Fraction')
        Weadd = Variable('W_{eadd}', 'lbf', 'Additional Engine System Weight')
        Wengsys = Variable('W_{engsys}', 'lbf', 'Total Engine System Weight')
        rvnace = Variable('r_{v_{nacelle}}', '-', 'Incoming Nacelle Velocity Ratio')

        Ceng = Variable('C_{engsys}', 1, '-', 'Engine System Weight Margin/Sens Factor')

        #BLI total drag reduction factor
        Dreduct = Variable('D_{reduct}', '-', 'BLI Drag Reduction Factor')
        Dwakefrac = Variable('D_{wakefraction}', 0.33, '-', 'Percent of Total Drag From Wake Dissipation')
        BLI_wake_benefit = Variable('BLI_{wakebenefit}', 0.02, '-', 'Wake Drag Reduction from BLI Wake Ingestion')
     
        constraints = []
        with SignomialsEnabled():
            constraints.extend([
                            self.wing['c_{root}'] == self.fuse['c_0'],
                            self.wing.wb['r_{w/c}'] == self.fuse['r_{wb}'],
                            self.wing['x_w'] == self.fuse['x_{wing}'],
                            self.wing['V_{ne}'] == Vmn,
                            self.VT['V_{ne}'] == Vmn,

                            #compute the aircraft's zero fuel weight
                            TCS([self.fuse['W_{fuse}'] + numeng \
                                * Wengsys + self.fuse['W_{tail}'] + Wwing + Wmisc <= W_dry]),

                            # Total takeoff weight constraint
                            TCS([W_ftotal + W_dry + self.fuse['W_{payload}'] <= W_total]),
                            TCS([W_ftotal >= W_fprimary + ReserveFraction * W_fprimary]),
                            TCS([W_fprimary >= W_fclimb + W_fcruise]),

                            # Load factor matching
                            self.fuse['N_{lift}'] == self.wing['N_{lift}'], # To make sure that the loads factors match.
                            Ltow*self.wing['L_{max}'] >= self.wing['N_{lift}'] * W_total + self.HT['L_{h_{max}}'],

                            # Wing fuel constraints
                            self.wing['W_{fuel_{wing}}'] >= W_ftotal/self.wing['FuelFrac'],

                            # Lifting surface weights
                            Wwing == self.wing['W_{wing_system}'],
                            WHT == self.HT['W_{HT_system}'],
                            WVT == self.VT['W_{VT_system}'],

                            # LG and Power Systems weights
                            Wmisc >= Wlg + Whpesys,
                            Wlgnose == flgnose*W_total,
                            Wlgmain == flgmain*W_total,
                            Wlg >= Clg*(Wlgnose + Wlgmain),
                            Whpesys == fhpesys*W_total,

                            # LG and Power System locations
                            xlgnose <= self.fuse['l_{nose}'],
                            TCS([xlgnose >= 0.6*self.fuse['l_{nose}']]),
                            TCS([xlgmain >= self.fuse['x_{wing}']]),
                            xlgmain <= self.wing['\\Delta x_{AC_{wing}}'] + self.fuse['x_{wing}'],
                            xhpesys == 1.1*self.fuse['l_{nose}'],
                            xmisc*Wmisc >= xlgnose*Wlgnose + xlgmain*Wlgmain + xhpesys*Whpesys,

                            # Tail cone sizing
                            3. * (numVT*self.VT['M_r']) * self.VT['c_{root_{vt}}'] * \
                                (self.fuse['p_{\\lambda_v}'] - 1.) >= numVT*self.VT[
                                    'L_{v_{max}}'] * self.VT['b_{vt}'] * (self.fuse['p_{\\lambda_v}']),
                            TCS([self.fuse['V_{cone}'] * (1. + self.fuse['\\lambda_{cone}']) * \
                             (pi + 4. * self.fuse['\\theta_{db}']) >= numVT*self.VT[
                                'M_r'] * self.VT['c_{root_{vt}}'] / self.fuse['\\tau_{cone}'] * \
                                 (pi + 2. * self.fuse['\\theta_{db}']) * \
                                  (self.fuse['l_{cone}'] / self.fuse['R_{fuse}'])]), #[SP]

                            # Lift curve slope ratio for HT and Wing
                            SignomialEquality(self.HT['m_{ratio}']*(1+2/self.wing['AR']), 1 + 2/self.HT['AR_{ht}']),

                            # HT Location and Volume Coefficient
                            self.HT['x_{CG_{ht}}'] <= self.fuse['l_{fuse}'],
                            self.fuse['x_{tail}'] == self.VT['x_{CG_{vt}}'],
                            TCS([self.HT['V_{ht}'] == self.HT['S_{ht}']*self.HT['l_{ht}']/(self.wing['S']*self.wing['mac'])]),

                            # HT Max Loading
                            TCS([self.HT['L_{h_{max}}'] >= 0.5*rhoTO*Vmn**2*self.HT['S_{ht}']*self.HT['C_{L_{hmax}}']]),

                            # Tail weight
                            self.fuse['W_{tail}'] >= numVT*WVT + WHT + self.fuse['W_{cone}'],

                            # VT volume coefficient
                            self.VT['V_{vt}'] == numVT*self.VT['S_{vt}'] * self.VT['l_{vt}']/(self.wing['S']*self.wing['b']),

                            # VT sizing constraints
                            # Yaw rate constraint at flare
                            numVT*.5*self.VT['\\rho_{TO}']*self.VT['V_{land}']**2*self.VT['S_{vt}']*self.VT['l_{vt}']* \
                                            self.VT['C_{L_{vyaw}}'] >= self.VT['\\dot{r}_{req}']*self.VT['I_{z}'],

                            # Force moment balance for one engine out condition
                            numVT*self.VT['L_{vtEO}']*self.VT['l_{vt}'] >= self.VT['T_e']*self.VT['y_{eng}'] + \
                                        self.VT['D_{wm}']*self.VT['y_{eng}'],
                            # TASOPT 2.0 p45

                            # Vertical bending material coefficient (VT aero loads)
                            self.fuse['B_{1v}'] == self.fuse['r_{M_v}']*numVT*self.VT['L_{v_{max}}']/(self.fuse['w_{fuse}']*self.fuse['\\sigma_{M_v}']),

                            # Moment of inertia around z-axis
                            # SignomialEquality(self.VT['I_{z}'], Izwing + Iztail + Izfuse),
                            self.VT['I_{z}'] >= Izwing + Iztail + Izfuse,

                            # Fuselage width (numaisle comes in)
                            TCS([2.*self.fuse['w_{fuse}'] >= self.fuse['SPR'] * self.fuse['w_{seat}'] + \
                                 numaisle*self.fuse['w_{aisle}'] + 2. * self.fuse['w_{sys}'] + self.fuse['t_{db}']]),

                            #engine system weight constraints
                            Snace == rSnace * np.pi * 0.25 * self.engine['d_{f}']**2,
                            lnace == 0.15 * self.engine['d_{f}'] * rSnace,
                            fSnace == Snace * self.wing['S']**-1,
                            Ainlet == 0.4 * Snace,
                            Afancowl == 0.2 * Snace,
                            Aexh == 0.4 * Snace,
                            Acorecowl == 3. * np.pi * self.engine['d_{LPC}']**2,
                            TCS([Wnace >= ((2.5+ 0.238*self.engine['d_{f}']/units('in')) * Ainlet + 1.9*Afancowl \
                                + (2.5+ 0.0363*self.engine['d_{f}']/units('in'))*Aexh + 1.9*Acorecowl)*units('lbf/ft^2'),
                            Weadd == feadd * self.engine['W_{engine}']]),
                            TCS([Wpylon >= (Wnace + Weadd + self.engine['W_{engine}']) * fpylon]),
                            TCS([Wengsys >= Ceng*(Wpylon + Wnace + Weadd + self.engine['W_{engine}'])]),
                            ])

        if rearengine and BLI:
            constraints.extend({self.VT['y_{eng}'] == 0.5 * self.fuse['w_{fuse}']})# Engine out moment arm
        if rearengine and not BLI:
            constraints.extend({self.VT['y_{eng}'] >= self.fuse['w_{fuse}'] + 0.5*self.engine['d_{f}'] + 1.*units('ft')})

        ### ENGINE LOCATION RELATED CONSTRAINTS

        # Wing-engined aircraft constraints
        if wingengine:
            with SignomialsEnabled():
                constraints.extend([
                    # Wing root moment constraint, with wing and engine weight load relief
                    TCS([self.wing['M_r']*self.wing['c_{root}'] >= (self.wing['L_{max}'] - self.wing['N_{lift}']*(Wwing+f_wingfuel*W_ftotal)) * \
                        (1./6.*self.wing['A_{tri}']/self.wing['S']*self.wing['b'] + \
                                1./4.*self.wing['A_{rect}']/self.wing['S']*self.wing['b']) - \
                                        self.wing['N_{lift}']*Wengsys*self.VT['y_{eng}']]), #[SP]

                    # Horizontal tail aero+landing loads constants A1h
                    self.fuse['A_{1h_{Land}}'] >= (self.fuse['N_{land}'] * \
                                (self.fuse['W_{tail}'] + self.fuse['W_{apu}'])) / \
                                 (self.fuse['h_{fuse}'] * self.fuse['\\sigma_{bend}']),

                    self.fuse['A_{1h_{MLF}}'] >= (self.fuse['N_{lift}'] * \
                                (self.fuse['W_{tail}'] + self.fuse['W_{apu}']) \
                                + self.fuse['r_{M_h}'] * self.HT['L_{h_{max}}']) / \
                                 (self.fuse['h_{fuse}'] * self.fuse['\\sigma_{M_h}']),

                    # Moment of inertia constraints
                    Izwing >= numeng*Wengsys*self.VT['y_{eng}']**2./g + \
                                    (self.wing['W_{fuel_{wing}}'] + Wwing)/(self.wing['S']*g)* \
                                    self.wing['c_{root}']*self.wing['b']**3.*(1./12.-(1.-self.wing['\\lambda'])/16.), #[SP]
                    Iztail >= (self.fuse['W_{apu}'] + self.fuse['W_{tail}'])*self.VT['l_{vt}']**2/g,
                            #NOTE: Using xwing as a CG surrogate. Reason: xCG moves during flight; want scalar Izfuse
                    Izfuse >= (self.fuse['W_{fuse}'] + self.fuse['W_{payload}'])/self.fuse['l_{fuse}'] * \
                                    (self.fuse['x_{wing}']**3 + self.VT['l_{vt}']**3.)/(3.*g),
                ])

        # Rear-engined aircraft constraints
        if rearengine:
            with SignomialsEnabled():
                constraints.extend([
                    # Wing root moment constraint, with wing weight + fuel load relief
                    TCS([self.wing['M_r']*self.wing['c_{root}'] >= (self.wing['L_{max}'] - self.wing['N_{lift}'] * (Wwing+f_wingfuel*W_ftotal)) * \
                        (1./6.*self.wing['A_{tri}']/self.wing['S']*self.wing['b'] + \
                                1./4.*self.wing['A_{rect}']/self.wing['S']*self.wing['b'])]), #[SP]

                    # Horizontal tail aero+landing loads constants A1h
                    self.fuse['A_{1h_{Land}}'] >= (self.fuse['N_{land}'] * \
                                                (self.fuse['W_{tail}'] + numeng * Wengsys + self.fuse['W_{apu}'])) / \
                                                (self.fuse['h_{fuse}'] * self.fuse['\\sigma_{bend}']),

                    self.fuse['A_{1h_{MLF}}'] >= (self.fuse['N_{lift}'] * \
                                               (self.fuse['W_{tail}'] + numeng * Wengsys + self.fuse['W_{apu}']) \
                                               + self.fuse['r_{M_h}'] * self.HT['L_{h_{max}}']) / \
                                                (self.fuse['h_{fuse}'] * self.fuse['\\sigma_{M_h}']),

                    # Moment of inertia constraints
                    Izwing >= (self.wing['W_{fuel_{wing}}'] + Wwing) / (self.wing['S'] * g) * \
                    self.wing['c_{root}'] * self.wing['b'] ** 3. * (1. / 12. - (1. - self.wing['\\lambda']) / 16.),
                    # [SP]
                    Iztail >= (self.fuse['W_{apu}'] + numeng * Wengsys + self.fuse['W_{tail}']) * self.VT[
                        'l_{vt}'] ** 2. / g,
                    # NOTE: Using xwing as a CG surrogate. Reason: xCG moves during flight; want scalar Izfuse
                    Izfuse >= (self.fuse['W_{fuse}'] + self.fuse['W_{payload}']) / self.fuse['l_{fuse}'] * \
                    (self.fuse['x_{wing}'] ** 3. + self.VT['l_{vt}'] ** 3.) / (3. * g),
                ])

        ### FUSELAGE CONSTRAINTS

        # Double-bubble
        if doublebubble:
            with SignomialsEnabled():
                constraints.extend([
                    # Floor loading
                    self.fuse['S_{floor}'] == (5. / 16.) * self.fuse['P_{floor}'],
                    self.fuse['M_{floor}'] == 9. / 256. * self.fuse['P_{floor}'] * self.fuse['w_{floor}'],
                    self.fuse['\\Delta R_{fuse}'] == self.fuse['R_{fuse}'] * 0.43/1.75,
                ])
        # Tube
        if tube:
            with SignomialsEnabled():
                constraints.extend([
                   # Floor loading
                    self.fuse['S_{floor}'] == 1./2. * self.fuse['P_{floor}'],
                    self.fuse['M_{floor}'] == 1./4. * self.fuse['P_{floor}']*self.fuse['w_{floor}'],
                ])

        ### VERTICAL TAIL CONSTRAINTS

        ### HORIZONTAL TAIL CONSTRAINTS
        # Pi HT constraints:
        if piHT:
            with SignomialsEnabled():
                constraints.extend([
                    # Pin VT joint moment constraint #TODO may be problematic, should check
                    SignomialEquality(self.HT['L_{h_{rect}}'] * (self.HT['b_{ht}'] / 4. - self.fuse['w_{fuse}']),
                                      self.HT['L_{h_{tri}}'] * (self.fuse['w_{fuse}'] - self.HT['b_{ht}'] / 6.)), # [SP] #[SPEquality]

                    # HT outboard half-span
                    self.HT['b_{ht_{out}}'] >= 0.5*self.HT['b_{ht}'] - self.fuse['w_{fuse}'], # [SP]

                    # Triangular HT lift outboard
                    self.HT['L_{h_{tri_{out}}}'] >= self.HT['L_{h_{tri}}'] * \
                                self.HT['b_{ht_{out}}']**2 / self.fuse['w_{fuse}']**2,

                    # Rectangular HT lift outboard
                    self.HT['L_{h_{rect_{out}}}'] >= self.HT['L_{h_{rect}}'] * \
                            self.HT['b_{ht_{out}}'] / self.fuse['w_{fuse}'],

                    # HT center moment
                    self.HT['M_r'] * self.HT['c_{root_{ht}}'] >= self.HT['L_{h_{rect}}'] * (
                    self.HT['b_{ht}'] / 4.) + self.HT['L_{h_{tri}}'] * (self.HT['b_{ht}'] / 6.) - \
                    self.fuse['w_{fuse}'] * self.HT['L_{h_{max}}'] / 2.,

                    # HT joint moment
                    self.HT['M_{r_{out}}']*self.HT['c_{attach}'] >= self.HT['L_{h_{rect_{out}}}'] * (0.5*self.HT['b_{ht_{out}}']) + \
                                                                    self.HT['L_{h_{tri_{out}}}'] * (1./3.*self.HT['b_{ht_{out}}']),

                    # HT joint shear (max shear)
                    self.HT['L_{shear}'] >= self.HT['L_{h_{rect_{out}}}'] + self.HT['L_{h_{tri_{out}}}'],

                    # constraint to lower bound M_r w.r.t. a conventional tail config (0.25*M_r of conventional)
                    self.HT['M_r']*self.HT['c_{root_{ht}}'] >= 0.25*(1./6.*self.HT['L_{h_{tri}}']*self.HT['b_{ht}'] + \
                         1./4.*self.HT['L_{h_{rect}}']*self.HT['b_{ht}']), # [SP] TODO: Remove when good HT structural model is ready.

                    # HT/VT joint constraint
                    self.HT['b_{ht}'] / (2. * self.fuse['w_{fuse}']) * self.HT['\lambda_{ht}'] * self.HT['c_{root_{ht}}'] ==
                    self.HT['c_{attach}'],
                ])
        # Conventional HT constraints
        if conventional:
            with SignomialsEnabled():
                constraints.extend([
                    # HT root moment
                    TCS([self.HT['M_r']*self.HT['c_{root_{ht}}'] >= 1./6.*self.HT['L_{h_{tri}}']*self.HT['b_{ht}'] + \
                         1./4.*self.HT['L_{h_{rect}}']*self.HT['b_{ht}']]),
                    # HT joint constraint
                    self.HT['c_{attach}'] == self.HT['c_{root_{ht}}'],

                    # HT auxiliary variables
                    self.HT['b_{ht_{out}}'] == 0.5*self.HT['b_{ht}'],
                    self.HT['L_{h_{tri_{out}}}'] == self.HT['L_{h_{tri}}'],
                    self.HT['L_{h_{rect_{out}}}'] == self.HT['L_{h_{rect}}'],
                    self.HT['M_{r_{out}}'] == self.HT['M_r'],
                    self.HT['L_{shear}'] >= self.HT['L_{h_{rect_{out}}}'] + self.HT['L_{h_{tri_{out}}}'],
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
        self.VTP = aircraft.VT.dynamic(state, aircraft.fitDrag)
        self.HTP = aircraft.HT.dynamic(state, aircraft.fitDrag)
        self.Pmodels = [self.wingP, self.fuseP, self.VTP, self.HTP]

        # Variable Definitions
        Vstall = Variable('V_{stall}',120., 'knots', 'Aircraft Stall Speed')
        D = Variable('D', 'N', 'Total Aircraft Drag')
        C_D = Variable('C_D', '-', 'Total Aircraft Drag Coefficient')
        LoD = Variable('L/D','-','Lift-to-Drag Ratio')
        W_avg = Variable(
            'W_{avg}', 'lbf', 'Geometric Average of Segment Start and End Weight')
        W_start = Variable('W_{start}', 'lbf', 'Segment Start Weight')
        W_end = Variable('W_{end}', 'lbf', 'Segment End Weight')
        W_burn = Variable('W_{burn}', 'lbf', 'Segment Fuel Burn Weight')
        WLoadmax = Variable('W_{Load_max}',6664., 'N/m^2', 'Max Wing Loading')
        WLoad = Variable('W_{Load}', 'N/m^2', 'Wing Loading')
        t = Variable('tmin', 'min', 'Segment Flight Time in Minutes')
        thours = Variable('thr', 'hour', 'Segment Flight Time in Hours')

        # Longitudinal stability variables
        xAC = Variable('x_{AC}','m','Aerodynamic Center of Aircraft')
        xCG = Variable('x_{CG}','m','Center of Gravity of Aircraft')
        xNP = Variable('x_{NP}','m','Neutral Point of Aircraft')
        SM = Variable('SM','-','Stability Margin of Aircraft')
        PCFuel = Variable('F_{fuel}','-','Percent Fuel Remaining (end of segment)')

        # Buoyancy weight variables
        Pcabin = Variable('P_{cabin}','Pa','Cabin Air Pressure')
        W_buoy = Variable('W_{buoy}','lbf','Buoyancy Weight')
        Tcabin = Variable('T_{cabin}','K','Cabin Air Temperature')
        rhocabin = Variable('\\rho_{cabin}','kg/m^3','Cabin Air Density')

        # Lift fraction variables
        Ltotal = Variable('L_{total}','N','Total lift')

        #variables for nacelle drag calcualation
        Vnace = Variable('V_{nacelle}', 'm/s', 'Incoming Nacelle Flow Velocity')
        V2 = Variable('V_2', 'm/s', 'Interior Nacelle Flow Velcoity')
        Vnacrat = Variable('V_{nacelle_ratio}', '-', 'Vnle/Vinf')
        rvnsurf = Variable('r_{v_{nsurf}}', '-', 'Intermediate Nacelle Drag Parameter')
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

            #Cabin Air properties
            rhocabin == Pcabin/(state['R']*Tcabin),
            Pcabin == 75000*units('Pa'),
            Tcabin == 297*units('K'),

            # speed must be greater than stall speed
            state['V'] >= Vstall,

            # Drag calculations
            self.fuseP['D_{fuse}'] == 0.5 * state['\\rho'] * state['V']**2 * \
                                        self.fuseP['C_{D_{fuse}}'] * aircraft['l_{fuse}'] * aircraft['R_{fuse}'] * (state['M']**2/aircraft.fuse['M_{fuseD}']**2),
            D >= aircraft['D_{reduct}'] * (self.wingP['D_{wing}'] + self.fuseP['D_{fuse}'] + self.aircraft['n_{VT}']*self.VTP['D_{vt}'] + self.HTP['D_{ht}'] + aircraft['n_{eng}'] * Dnace),
            C_D == D/(.5*state['\\rho']*state['V']**2 * self.aircraft.wing['S']),
            LoD == W_avg/D,

            # Wing loading
            WLoad == .5 * self.wingP['C_{L}'] * self.aircraft['S'] * state.atm['\\rho'] * state['V']**2 / self.aircraft.wing['S'],

            # Center wing lift loss
            # self.wingP['p_{o}'] >= self.wingP['L_w']*self.wing['c_{root}']*(.5 + 0.5*self.wingP['\\eta_{o}'](/(self.wing['S']),
            self.wingP['p_{o}'] >= self.wingP['L_w']*aircraft.wing['c_{root}']/(aircraft.wing['S']), #TODO improve approx without making SP
            self.wingP['\\eta_{o}'] == aircraft['w_{fuse}']/(aircraft['b']/2),

            # Fuselage lift (just calculating)
            SignomialEquality(self.fuseP['L_{fuse}'], (self.aircraft['L_{total/wing}']-1.)*self.wingP['L_w']),

            # Geometric average of start and end weights of flight segment
            W_avg >= (W_start * W_end)**.5 + W_buoy, # Buoyancy weight included in Breguet Range

            # Maximum wing loading constraint
            WLoad <= WLoadmax,

            # Flight time unit conversion
            t == thours,

            #VTP constraints
            aircraft.VT['x_{CG_{vt}}'] + 0.5*aircraft.VT['c_{root_{vt}}'] <= aircraft.fuse['l_{fuse}'],

            # Drag of a windmilling engine (VT sizing)
            TCS([aircraft.VT['D_{wm}'] >= 0.5*aircraft.VT['\\rho_{TO}']*aircraft.VT['V_1']**2.*aircraft.engine['A_{2}']*aircraft.VT['C_{D_{wm}}']]),

            # Aircraft trim conditions
            TCS([xAC/aircraft.wing['mac'] <= xCG/aircraft.wing['mac'] + \
                 self.wingP['c_{m_{w}}']/self.wingP['C_{L}']  +\
                              aircraft.HT['V_{ht}']*(self.HTP['C_{L_h}']/self.wingP['C_{L}'])]),

            # Tail aspect ratio and lift constraints
            aircraft.HT['AR_{ht}'] >= 3., #TODO change to tip Re constraint
            self.HTP['C_{L_h}'] >= 0.01, #TODO remove

            # HT/VT moment arm constraints
            aircraft.HT['l_{ht}'] <= aircraft.HT['x_{CG_{ht}}'] - xCG,
            aircraft.VT['l_{vt}'] <= aircraft.VT['x_{CG_{vt}}'] - xCG,

           # Tail downforce penalty to total lift
            TCS([Ltotal == self.aircraft['L_{total/wing}']*self.wingP['L_w']]),
            TCS([Ltotal >= W_avg + self.HTP['L_h']]),

            # Wing location and AC constraints

            TCS([xAC <= aircraft['x_{wing}'] + 0.25*aircraft['\\Delta x_{AC_{wing}}'] + xNP]), #[SP] #TODO relax and improve
            TCS([SM <= (xAC-xCG)/aircraft['mac']]),
            SM >= aircraft['SM_{min}'],

            # Neutral point approximation (taken from Basic Aircraft Design Rules, Unified)
            # TODO improve
            SignomialEquality(xNP/aircraft['mac']/aircraft['V_{ht}']*(aircraft['AR']+2.)*(1.+2./aircraft['AR_{ht}']),
                              (1.+2./aircraft['AR'])*(aircraft['AR']-2.)),

            # HT Location constraints
            aircraft.HT['x_{CG_{ht}}'] <= aircraft.fuse['l_{fuse}'],

            # Static margin constraints
            self.wingP['c_{m_{w}}'] == 1.9,
              
            TCS([aircraft['SM_{min}'] + aircraft['\\Delta x_{CG}']/aircraft.wing['mac'] \
                 + self.wingP['c_{m_{w}}']/aircraft.wing['C_{L_{wmax}}'] <= \
                                            aircraft.HT['V_{ht}']*aircraft.HT['m_{ratio}'] +\
                                            aircraft.HT['V_{ht}']*aircraft.HT['C_{L_{hmax}}']/aircraft.wing['C_{L_{wmax}}']]), # [SP]

            #nacelle drag
            Renace == state['\\rho']*state['V'] * aircraft['l_{nacelle}']/state['\\mu'],
            Cfnace == 0.94*4.*0.0743/(Renace**(0.2)), #from http://www.calpoly.edu/~kshollen/ME347/Handouts/Friction_Drag_Coef.pdf
            Vnace == aircraft['r_{v_{nacelle}}'] * state['V'],
            Vnacrat >= 2.*Vnace/state['V'] - V2/state['V'],
            rvnsurf**3. >= 0.25*(Vnacrat + aircraft['r_{v_{nacelle}}'])*(Vnacrat**2. + aircraft['r_{v_{nacelle}}']**2.),
            Cdnace == aircraft['f_{S_{nacelle}}'] * Cfnace[0] * rvnsurf **3.,
            Dnace == Cdnace * 0.5 * state['\\rho'] * state['V']**2. * aircraft['S'],
            ])

        if not aircraft.fitDrag:
            constraints.extend([
                #set the VT drag coefficient
                self.VTP['C_{D_{vis}}'] >= (self.aircraft.VT['c_{d_{fv}}'] + self.aircraft.VT['c_{d_{pv}}']*self.aircraft.VT['\\cos(\\Lambda_{vt})^3']),

                #set the HT drag coefficient
                self.HTP['C_{D_{0_h}}'] >= (self.aircraft.HT['c_{d_{fh}}'] + self.aircraft.HT['c_{d_{ph}}']*self.aircraft.HT['\\cos(\\Lambda_{ht})^3']),
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
        RC = Variable('RC', 'feet/min', 'Rate of Climb/Descent')
        dhft = Variable(
            'dhft', 'feet', 'Change in Altitude Per Climb Segment [feet]')
        RngClimb = Variable('RngClimb', 'nautical_miles',
                            'Down Range Covered in Each Climb Segment')

        constraints = []
        constraints.extend([
            # Excess power for climb
            TCS([excessP + state['V'] * self.aircraftP['D'] <= state['V']
                 * aircraft['n_{eng}'] * self.engine['F'][:Nclimb]]),

            RC == excessP / self.aircraftP['W_{avg}'],
            RC >= 500. * units('ft/min'),

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
        # z_bre = Variable('z_{bre}', '-', 'Breguet Parameter')
        Rng = Variable('Rng', 'nautical_miles', 'Cruise Segment Range')
        RC = Variable('RC', 'feet/min', 'Rate of Climb/Descent')
        theta = Variable('\\theta','-','Climb Angle')
        dhft = Variable('dhft', 'feet', 'Change in Altitude Per Climb Segment [feet]')

        constraints = []

        constraints.extend([
            RC == theta*state['V'],
            RC >= 0.01 * units('ft/min'),

            # Time
            self.aircraftP['thr'] * state['V'] == Rng,
            dhft == self.aircraftP['tmin'] * RC,
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
        if conventional:
             statevarkeys = ['L_{atm}', 'M_{atm}', 'P_{atm}', 'R_{atm}',
                             '\\rho', 'T_{atm}', '\\mu', 'T_s', 'C_1', 'h', 'hft', 'V', 'a', 'R', '\\gamma', 'M']
        else:
             statevarkeys = ['P_{atm}', 'R_{atm}',
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
    Inputs:
    Nclimb: number of climb segments (for Brequet Range)
    Ncruise: number of cruise segments (for Brequet Range)
    objective: defines the objective function
    airplane: string representing the aircraft model
    Nmission: specifies whether single-point or multi-point optimization
              Nmission >/= 1 requires specification of range and number of passengers for each mission
    """

    def setup(self, Nclimb, Ncruise, objective, airplane, Nmission = 1):
        # define the number of each flight segment

        global D80, D82, D82, D82_73eng, D8_eng_wing, D8big, b737800, b777300ER, optimal737, \
               optimalD8, Mo8D8, M08_D8_eng_wing, M072_737, D8fam, D8_no_BLI, \
               M08D8_noBLI, optimal777, D8big_eng_wing, multimission, \
               D8bigfam, optimalRJ, RJfam, smallD8, smallD8_no_BLI, smallD8_eng_wing, D12
        global manufacturer, operator, fuel
        global wingengine, rearengine, doublebubble, tube, piHT, conventional

        wingengine = False; rearengine = False; doublebubble = False; tube = False;
        piHT = False; conventional = False

        # Objective type, only one active at once
        manufacturer = False; operator = False; fuel = False; PRFC = False
        if objective == 'manufacturer':
            manufacturer = True
        if objective == 'operator':
            operator = True
        if objective == 'fuel':
            fuel = True
        if objective == 'PRFC':
            PRFC = True

        # Aircraft type, only one active at once
        D80 = False
        D82 = False
        D82_73eng = False
        D8_eng_wing = False
        D8big = False
        D8big_eng_wing = False
        b737800 = False
        b777300ER = False
        optimal737 = False
        optimalD8 = False
        optimal777 = False
        M08D8 = False
        M08D8_noBLI = False
        M08_D8_eng_wing = False
        M072_737 = False
        D8_no_BLI = False
        D8big_no_BLI = False
        D8big_M072 = False
        D8big_M08 = False
        optimalRJ = False
        smallD8 = False
        smallD8_no_BLI = False
        smallD8_eng_wing = False
        smallD8_M08_eng_wing = False
        smallD8_M08 = False
        smallD8_M08_no_BLI = False
        D12 = False
        optimal777_M08 = False
        optimal777_M072 = False
        D8big_M072 = False
        D8big_eng_wing_M072 = False
        D8big_no_BLI_M072 = False

        if airplane == 'D80':
            D80 = True; rearengine = True; BLI = True; piHT = True; doublebubble = True;
        if airplane == 'D82':
            D82 = True; rearengine = True; BLI = True; piHT = True; doublebubble = True;
        if airplane == 'D82_73eng':
            D82_73eng = True; rearengine = True; BLI = True; piHT = True; doublebubble = True;
        if airplane == 'D8_eng_wing':
            D8_eng_wing = True; wingengine = True; piHT = True; doublebubble = True;
        if airplane == 'D8big':
            D8big = True; rearengine = True; BLI = True; piHT = True; doublebubble = True;
        if airplane == 'D8big_no_BLI':
            D8big_no_BLI = True; rearengine = True; piHT = True; doublebubble = True;
        if airplane == 'D8big_eng_wing':
            D8big_eng_wing = True; wingengine = True; piHT = True; doublebubble = True;
        if airplane == 'D8big_M072':
            D8big = True
            D8big_M072 = True; rearengine = True; BLI = True; piHT = True; doublebubble = True;
        if airplane == 'D8big_no_BLI_M072':
            D8big_no_BLI = True
            D8big_no_BLI_M072 = True; rearengine = True; piHT = True; doublebubble = True;
        if airplane == 'D8big_eng_wing_M072':
            D8big_eng_wing = True
            D8big_eng_wing_M072 = True; wingengine = True; piHT = True; doublebubble = True;
        if airplane == 'b737800':
            b737800 = True; conventional = True
        if airplane == 'b777300ER':
            b777300ER = True; conventional = True
        if airplane == 'optimal737':
            optimal737 = True; conventional = True
        if airplane == 'optimalD8':
            optimalD8 = True; rearengine = True; BLI = True; piHT = True; doublebubble = True;
        if airplane == 'optimal777':
            optimal777 = True; conventional = True
        if airplane == 'optimal777_M08':
            optimal777 = True
            optimal777_M08 = True; conventional = True
        if airplane == 'optimal777_M072':
            optimal777 = True
            optimal777_M072 = True; conventional = True
        if airplane == 'M08D8':
            M08D8 = True; rearengine = True; BLI = True; piHT = True; doublebubble = True;
        if airplane == 'M08D8_noBLI':
            M08D8_noBLI = True; rearengine = True; piHT = True; doublebubble = True;
        if airplane == 'M08_D8_eng_wing':
            M08_D8_eng_wing = True; wingengine = True; piHT = True; doublebubble = True;
        if airplane == 'M072_737':
            M072_737 = True; conventional = True
        if airplane == 'D8_no_BLI':
            D8_no_BLI = True; rearengine = True; piHT = True; doublebubble = True;
        if airplane == 'D8big_M072':
            D8big = True
            D8big_M072 = True; rearengine = True; BLI = True; piHT = True; doublebubble = True;
        if airplane == 'D8big_M08':
            D8big = True
            D8big_M08 = True; rearengine = True; BLI = True; piHT = True; doublebubble = True;
        if airplane == 'optimalRJ':
            optimalRJ = True; conventional = True
        if airplane == 'smallD8':
            smallD8 = True; rearengine = True; BLI = True; piHT = True; doublebubble = True;
        if airplane == 'smallD8_eng_wing':
            smallD8_eng_wing = True; wingengine = True; piHT = True; doublebubble = True;
        if airplane == 'smallD8_no_BLI':
            smallD8_no_BLI = True; rearengine = True; piHT = True; doublebubble = True;
        if airplane == 'smallD8_M08_no_BLI':
            smallD8_no_BLI = True
            smallD8_M08_no_BLI = True; rearengine = True; piHT = True; doublebubble = True;
        if airplane == 'smallD8_M08':
            smallD8 = True
            smallD8_M08 = True; rearengine = True; BLI = True; piHT = True; doublebubble = True;
        if airplane == 'smallD8_M08_eng_wing':
            smallD8_eng_wing = True
            smallD8_M08_eng_wing = True; wingengine = True; piHT = True; doublebubble = True;
        if airplane == 'D12':
            D12 = True; rearengine = True; BLI = True; piHT = True; doublebubble = True;

        if conventional:
            wingengine = True; tube = True;

        # Multimission?
        if Nmission == 1:
            multimission = False
        else:
            multimission = True

        # Defining fitDrag, boolean describing whether or not to use tail drag fits
        fitDrag = None

        #eng and BLI specify enginestate and BLI parameters for input into the turbofan model.

        if D80 or D82 or optimalD8 or M08D8 or smallD8:
             eng = 3
             BLI = True

        if D8_eng_wing or D8_no_BLI or M08_D8_eng_wing or optimal737 or M08D8_noBLI or M072_737 \
           or optimalRJ or smallD8_eng_wing or smallD8_no_BLI:
            eng = 3
            BLI = False

        if b737800:
             eng = 1
             BLI = False

        if D82_73eng:
             eng = 1
             BLI = True

        if D8big or D12:
             eng = 4
             BLI = True

        if b777300ER or optimal777 or D8big_eng_wing or D8big_no_BLI:
             eng = 4
             BLI = False

        if optimalD8 or D80 or D82 or D82_73eng or D8big or M08D8 or D8_no_BLI or M08D8_noBLI or D8big_no_BLI or smallD8 or smallD8_no_BLI or D12:
            D8fam = True
        else:
            D8fam = False

        if D8big_eng_wing or D8big_no_BLI or D8big:
            D8bigfam = True
        else:
            D8bigfam = False

        if optimalRJ or smallD8 or smallD8_eng_wing or smallD8_no_BLI:
            RJfam = True
        else:
            RJfam = False

        # vectorize
        with Vectorize(Nmission):
             with Vectorize(Nclimb + Ncruise):
                 enginestate = FlightState()

        # True is use xfoil fit tail drag model, False is TASOPT tail drag model
        if optimal737 or optimalD8 or M08_D8_eng_wing or M08D8_noBLI or M08D8 or M072_737 or \
           D8_eng_wing or D8_no_BLI or D8big or optimal777 or D8big_eng_wing or D8big_no_BLI or RJfam or D12:
            fitDrag = True
        else:
            fitDrag = False

        # Build required submodels
        aircraft = Aircraft(Nclimb, Ncruise, enginestate, eng, fitDrag, BLI, Nmission)

        # Vectorize dynamic variables
        with Vectorize(Nmission):
             with Vectorize(Nclimb):
                 climb = ClimbSegment(aircraft, Nclimb)

        with Vectorize(Nmission):
             with Vectorize(Ncruise):
                 cruise = CruiseSegment(aircraft, Nclimb)

        # StateLinking links the climb and cruise state variables to the engine state,
        # so that atmospheric variables match.
        statelinking = StateLinking(climb.state, cruise.state, enginestate, Nclimb, Ncruise)

        # Declare Mission variables
        if multimission:
             with Vectorize(Nmission):
                  CruiseAlt = Variable('CruiseAlt', 'ft', 'Cruise Altitude [feet]')
                  ReqRng = Variable('ReqRng', 'nautical_miles', 'Required Cruise Range')
##                  Total_Time = Variable('TotalTime', 'hr', 'Total Mission Time')
                  climb_time = Variable('ClimbTime', 'min', 'Total Time in Climb')
                  climb_distance = Variable('ClimbDistance', 'nautical_miles', 'Climb Distance')
        else:
          CruiseAlt = Variable('CruiseAlt', 'ft', 'Cruise Altitude [feet]')
          ReqRng = Variable('ReqRng', 'nautical_miles', 'Required Cruise Range')
          Total_Time = Variable('TotalTime', 'hr', 'Total Mission Time')
          climb_time = Variable('ClimbTime', 'min', 'Total Time in Climb')
          climb_distance = Variable('ClimbDistance', 'nautical_miles', 'Climb Distance')

        max_climb_time = Variable('MaxClimbTime', 'min', 'Total Time in Climb')
        max_climb_distance = Variable('MaxClimbDistance', 'nautical_miles', 'Climb Distance')
        CruiseTt41max = Variable('T_{t_{4.1_{max-Cruise}}}', 'K', 'Max Cruise Turbine Inlet Temp')

        # make overall constraints
        constraints = []

        # Setting maximum OPR
        if RJfam:
            OPRmax = 30.
        elif D8bigfam or optimal777 or b777300ER or D12:
            OPRmax = 42.
        else:
            OPRmax = 35.

        constraints.extend([
            aircraft['OPR'][Nclimb] <= OPRmax,
            ])

        with SignomialsEnabled():
            # Buoyancy weight #TODO relax the equality
            # SignomialEquality(W_buoy,(rhocabin - state['\\rho'])*g*aircraft['V_{cabin}']),  #[SP] #[SPEquality]
            # Note: Buoyancy model has been simplified, since it causes significant increases in runtime.
            constraints.extend([
                cruise['W_{buoy}'] >= (cruise['\\rho_{cabin}'])*g*aircraft['V_{cabin}'], # [SP] # - cruise['\\rho']
                climb['W_{buoy}'] >= (climb['\\rho_{cabin}'])*g*aircraft['V_{cabin}'],
                aircraft['PRFC'] == aircraft['W_{f_{primary}}']/g*aircraft.engine['h_{f}']/(ReqRng*aircraft['W_{payload}'])
            ])

            ### CG CONSTRAINTS
            if rearengine:
                constraints.extend([
                TCS([climb['x_{CG}']*climb['W_{end}'] >=
                    aircraft['x_{misc}']*aircraft['W_{misc}'] \
                    + 0.5*(aircraft.fuse['W_{fuse}']+aircraft.fuse['W_{payload}'])*aircraft.fuse['l_{fuse}'] \
                    + (aircraft['W_{tail}']+aircraft['n_{eng}']*aircraft['W_{engsys}'])*aircraft['x_{tail}'] \
                    + (aircraft['W_{wing_system}']*(aircraft.fuse['x_{wing}']+aircraft.wing['\\Delta x_{AC_{wing}}'])) \
                    + (climb['F_{fuel}']+aircraft['ReserveFraction'])*aircraft['W_{f_{primary}}'] \
                    * (aircraft.fuse['x_{wing}']+aircraft.wing['\\Delta x_{AC_{wing}}']*climb['F_{fuel}']) \
                    ]),
                TCS([cruise['x_{CG}']*cruise['W_{end}'] >=
                    aircraft['x_{misc}']*aircraft['W_{misc}'] \
                    + 0.5*(aircraft.fuse['W_{fuse}']+aircraft.fuse['W_{payload}'])*aircraft.fuse['l_{fuse}'] \
                    + (aircraft['W_{tail}']+aircraft['n_{eng}']*aircraft['W_{engsys}'])*aircraft['x_{tail}'] \
                    + (aircraft['W_{wing_system}']*(aircraft.fuse['x_{wing}']+aircraft.wing['\\Delta x_{AC_{wing}}'])) \
                    + (cruise['F_{fuel}']+aircraft['ReserveFraction'])*aircraft['W_{f_{primary}}'] \
                    * (aircraft.fuse['x_{wing}']+aircraft.wing['\\Delta x_{AC_{wing}}']*cruise['F_{fuel}'])
                     ]),
              ])
            if wingengine:
                constraints.extend([
                TCS([climb['x_{CG}']*climb['W_{end}'] >=
                    aircraft['x_{misc}']*aircraft['W_{misc}'] \
                    + 0.5*(aircraft.fuse['W_{fuse}']+aircraft.fuse['W_{payload}'])*aircraft.fuse['l_{fuse}'] \
                    + (aircraft['W_{tail}'])*aircraft['x_{tail}'] \
                    + (aircraft['W_{wing_system}']*(aircraft.fuse['x_{wing}']+aircraft.wing['\\Delta x_{AC_{wing}}'])) \
                    + (climb['F_{fuel}']+aircraft['ReserveFraction'])*aircraft['W_{f_{primary}}'] \
                    * (aircraft.fuse['x_{wing}']+aircraft.wing['\\Delta x_{AC_{wing}}']*climb['F_{fuel}']) \
                    + aircraft['n_{eng}']*aircraft['W_{engsys}']*aircraft['x_b']]), # TODO improve; using x_b as a surrogate for xeng
                TCS([cruise['x_{CG}']*cruise['W_{end}'] >=
                    aircraft['x_{misc}']*aircraft['W_{misc}'] \
                    + 0.5*(aircraft.fuse['W_{fuse}']+aircraft.fuse['W_{payload}'])*aircraft.fuse['l_{fuse}'] \
                    + (aircraft['W_{tail}'])*aircraft['x_{tail}'] \
                    + (aircraft['W_{wing_system}']*(aircraft.fuse['x_{wing}']+aircraft.wing['\\Delta x_{AC_{wing}}'])) \
                    + (cruise['F_{fuel}']+aircraft['ReserveFraction'])*aircraft['W_{f_{primary}}'] \
                    * (aircraft.fuse['x_{wing}']+aircraft.wing['\\Delta x_{AC_{wing}}']*cruise['F_{fuel}'])
                    + aircraft['n_{eng}']*aircraft['W_{engsys}']*aircraft['x_b']]), # TODO improve; using x_b as a surrogate for xeng
              ])

            #Setting fuselage drag and lift, and BLI correction
            if optimalD8 or D80 or D82 or D82_73eng or M08D8 or D8_no_BLI or M08D8_noBLI or smallD8 or smallD8_no_BLI or D8_eng_wing or smallD8_eng_wing \
               or M08_D8_eng_wing:
                constraints.extend([
                    climb.climbP.fuseP['C_{D_{fuse}}'] == 0.018081,
                    cruise.cruiseP.fuseP['C_{D_{fuse}}'] == 0.018081,
                    aircraft.fuse['M_{fuseD}'] == 0.72,
                  ])

            if D12 or D8big_eng_wing or D8big_no_BLI or D8big:
                constraints.extend([
                        climb.climbP.fuseP['C_{D_{fuse}}'] == 0.0167620,
                        cruise.cruiseP.fuseP['C_{D_{fuse}}'] == 0.0167620,
                        aircraft.fuse['M_{fuseD}'] == 0.83,
                      ])
            if conventional and not (b777300ER or optimal777 or M072_737):
                constraints.extend([
                    #Setting fuselage drag coefficient
                    climb.climbP.fuseP['C_{D_{fuse}}'] == 0.01107365,
                    cruise.cruiseP.fuseP['C_{D_{fuse}}'] == 0.01107365,
                    aircraft.fuse['M_{fuseD}'] == 0.80,
                ])
            if M072_737:
                constraints.extend([
                    #Setting fuselage drag coefficient
                    climb.climbP.fuseP['C_{D_{fuse}}'] == 0.01107365,#0.0129077,
                    cruise.cruiseP.fuseP['C_{D_{fuse}}'] == 0.01107365,#0.0129077,
                    aircraft.fuse['M_{fuseD}'] == 0.80, #0.72,
                ])
            if b777300ER or optimal777:
                constraints.extend([
                    #Setting fuselage drag coefficient
                    #additioanl 1.1 factor accounts for mach drag rise model
                    climb.climbP.fuseP['C_{D_{fuse}}'] == 0.00987663,
                    cruise.cruiseP.fuseP['C_{D_{fuse}}'] == 0.00987663,
                    aircraft.fuse['M_{fuseD}'] == 0.84,
                ])

        constraints.extend([
            climb.climbP.aircraftP['W_{start}'][0] == aircraft['W_{total}'],
            climb.climbP.aircraftP[
                'W_{end}'][-1] == cruise.cruiseP.aircraftP['W_{start}'][0],

            # Climb segment weight decreases by the fuel burn...
            TCS([climb.climbP.aircraftP['W_{start}'] >= climb.climbP.aircraftP[
                'W_{end}'] + climb.climbP.aircraftP['W_{burn}']]),
            # Cruise segment weight decreasesby the fuel burn...
            TCS([cruise.cruiseP.aircraftP['W_{start}'] >= cruise.cruiseP.aircraftP[
                'W_{end}'] + cruise.cruiseP.aircraftP['W_{burn}']]),

            climb.climbP.aircraftP['W_{start}'][
            1:] == climb.climbP.aircraftP['W_{end}'][:-1],
            cruise.cruiseP.aircraftP['W_{start}'][
            1:] == cruise.cruiseP.aircraftP['W_{end}'][:-1],

            TCS([aircraft['W_{dry}'] + aircraft['W_{payload}'] + \
                 aircraft['ReserveFraction'] * aircraft['W_{f_{primary}}'] <= cruise.cruiseP.aircraftP['W_{end}'][-1]]),
            TCS([aircraft['W_{f_{climb}}'] >= sum(climb.climbP.aircraftP['W_{burn}'])]),
            TCS([aircraft['W_{f_{cruise}}'] >= sum(cruise.cruiseP.aircraftP['W_{burn}'])]),
            ])

        with SignomialsEnabled():
            constraints.extend([
                # WARNING: Arbitrary cruise altitude constraint
                CruiseAlt >= 25000. * units('ft'),
                
                # Altitude constraints
                climb['hft'][-1] >= CruiseAlt,
                SignomialEquality(climb['hft'][1:Nclimb], climb['hft'][:Nclimb - 1] + climb['dhft'][1:Nclimb]), #[SP]
                TCS([climb['hft'][0] == climb['dhft'][0]]),

                # All climb segments have the same total altitude change
                climb['dhft'][1:Nclimb] == climb['dhft'][:Nclimb - 1],

                # compute fuel burn from TSFC
                cruise.cruiseP.aircraftP['W_{burn}'] == aircraft['n_{eng}'] * aircraft.engine['TSFC'][Nclimb:] * \
                    cruise['thr'] * aircraft.engine['F'][Nclimb:],
                climb.climbP.aircraftP['W_{burn}'] == aircraft['n_{eng}'] * aircraft.engine['TSFC'][:Nclimb] * \
                    climb['thr'] * aircraft.engine['F'][:Nclimb],

                # Thrust >= Drag + Vertical Potential Energy
                aircraft['n_{eng}'] * aircraft.engine['F'][:Nclimb] >= climb['D'] + climb['W_{avg}'] * climb['\\theta'],
                aircraft['n_{eng}'] * aircraft.engine['F'][Nclimb:] >= cruise['D'] + cruise['W_{avg}'] * cruise['\\theta'],

                # Takeoff thrust T_e calculated for engine out + vertical tail sizing.
                # Note: coeff can be varied as desired.
                aircraft.VT['T_e'] == 1.8 * climb.climbP.engine['F'][0],

                # Set the range for each cruise segment.
                # All cruise segments cover the same range.
                cruise['Rng'][:Ncruise-1] == cruise['Rng'][1:Ncruise],

                # Cruise Mach Number constraint
                cruise['M'] >= aircraft['M_{min}'],

                # Nacelle drag constraint
                # Elevated this constraint to Mission for dimensionality
                cruise.cruiseP['V_2'] == aircraft.engine['M_2'][Nclimb:] * cruise.state['a'],
                climb.climbP['V_2'] == aircraft.engine['M_2'][:Nclimb] * climb.state['a'],

                climb['\\alpha_{max,w}'] == .18,
                cruise['\\alpha_{max,w}'] == .1,

                # T/O minimum climb rate constraint
                climb['RC'][0] >= 2500. * units('ft/min'),

                # TASOPT TOC climb rate constraint
                climb['\\theta'][-1] >= 0.015, #higher than 0.015 radian climb gradient at top-of-climb

                #compute the total time
                Total_Time >= sum(cruise['thr']) + sum(climb['thr']),

                #compute the climb in time
                climb_time >= sum(climb['thr']),
                climb_distance >= sum(cruise['Rng']),
                climb_time <= max_climb_time,
                climb_distance <= max_climb_distance,

                #set the max allowed cruise Tt4.1
                aircraft['T_{t_{4.1}}'][Nclimb:] <= CruiseTt41max,
                ])

        # Calculating percent fuel remaining
        with SignomialsEnabled():
            for i in range(0,Nclimb):
                constraints.extend([
                    TCS([climb['F_{fuel}'][i] >= (sum(climb['W_{burn}'][i+1:]) + \
                                                             aircraft['W_{f_{cruise}}'])/aircraft['W_{f_{primary}}']]) ,
                    climb['F_{fuel}'] <= 1.0, #just in case, TODO remove later
                ])
            for i in range(0,Ncruise):
                constraints.extend([
                    TCS([cruise['F_{fuel}'][i] >= (sum(cruise['W_{burn}'][i+1:]) + \
                                0.0000001*aircraft['W_{f_{primary}}'])/aircraft['W_{f_{primary}}']]),
                    cruise['F_{fuel}'] <= 1.0, #just in case, TODO remove later
                    ])

        with SignomialsEnabled():
            constraints.extend([
                #set the range constraints
                TCS([sum(climb['RngClimb']) + sum(cruise['Rng']) >= ReqRng]), #[SP]

                # Cruise climb constraint
                cruise['hft'][0] <= climb['hft'][-1] + cruise['dhft'][0], #[SP]
                cruise['hft'][1:Ncruise] <=  cruise['hft'][:Ncruise-1] + cruise['dhft'][1:Ncruise], #[SP]
                ])

        if multimission and not D8bigfam and not b777300ER and not optimal777 and not RJfam:
             W_fmissions = Variable('W_{f_{missions}', 'N', 'Fuel burn across all missions')
             constraints.extend([
                  W_fmissions >= sum(aircraft['W_{f_{total}}']),
                  aircraft['n_{seat}'] == aircraft['n_{pax}'][0], # TODO find a more robust way of doing this!
                  ])

        if not multimission and not D8bigfam and not b777300ER and not optimal777 and not RJfam and not D12:
             constraints.extend([
                  aircraft['n_{pax}'] == 180.,
                  aircraft['n_{seat}'] == aircraft['n_{pax}']
                  ])

        if multimission and (D8bigfam or b777300ER or optimal777):
             W_fmissions = Variable('W_{f_{missions}', 'N', 'Fuel burn across all missions')

             constraints.extend([
                  W_fmissions >= sum(aircraft['W_{f_{total}}']),

                  aircraft['n_{seat}'] == aircraft['n_{pax}'][0], # TODO find a more robust way of doing this!
                  ])
        if not multimission and (D8bigfam or b777300ER or optimal777):
             constraints.extend([
                  aircraft['n_{pax}'] == 450.,
                  aircraft['n_{seat}'] == aircraft['n_{pax}']
                  ])

        if multimission and RJfam:
             W_fmissions = Variable('W_{f_{missions}', 'N', 'Fuel burn across all missions')

             constraints.extend([
                  W_fmissions >= sum(aircraft['W_{f_{total}}']),

                  aircraft['n_{seat}'] == aircraft['n_{pax}'][0], # TODO find a more robust way of doing this!
                  ])
        if not multimission and RJfam:
             constraints.extend([
                  aircraft['n_{pax}'] == 90.,
                  aircraft['n_{seat}'] == aircraft['n_{pax}']
                  ])

        if multimission and D12:
             W_fmissions = Variable('W_{f_{missions}', 'N', 'Fuel burn across all missions')

             constraints.extend([
                  W_fmissions >= sum(aircraft['W_{f_{total}}']),
                  aircraft['n_{seat}'] == aircraft['n_{pax}'][0], # TODO find a more robust way of doing this!
                  ])
        if not multimission and D12:
             constraints.extend([
                  aircraft['n_{pax}'] == 500.,
                  aircraft['n_{seat}'] == aircraft['n_{pax}']
                  ])

        M2 = .6
        M25 = .6
        M4a = .2
        M0 = .5

        engineclimb = [
            aircraft.engine.engineP['M_2'][:Nclimb] == climb['M'],
            aircraft.engine.engineP['M_{2.5}'][:Nclimb] == M25,
            aircraft.engine.engineP['hold_{2}'][:Nclimb] == 1.+.5*(1.398-1.)*M2**2.,
            aircraft.engine.engineP['hold_{2.5}'][:Nclimb] == 1.+.5*(1.354-1.)*M25**2.,
            
            #climb rate constraints
            TCS([climb['excessP'] + climb.state['V'] * climb['D'] <= climb.state['V'] * aircraft['n_{eng}'] * aircraft.engine['F_{spec}'][:Nclimb]]),
            ]
    
        if D8fam or M072_737 or D8big_M072 or D8_eng_wing or smallD8_eng_wing or optimal777_M072 or D8big_M072 or D8big_eng_wing_M072 or D8big_no_BLI_M072:
             M2 = .6
             M25 = .6
             M4a = .2
             M0 = .72

        if D8_no_BLI or D8big_M08 or optimalRJ or M08D8 or smallD8_M08_eng_wing or smallD8_M08_no_BLI or smallD8_M08 or optimal777_M08:
             M2 = .6
             M25 = .6
             M4a = .2
             M0 = .8

        if (b777300ER or optimal777 or D8big_eng_wing or D8big or D8big_no_BLI) and not (D8big_M08 or D8big_M072 or optimal777_M08 or optimal777_M072 or \
                                                                                         D8big_M072 or D8big_eng_wing_M072 or D8big_no_BLI_M072):
             M2 = .65
             M25 = .6
             M4a = .2
             M0 = .84

        if D12:
             M2 = .65
             M25 = .6
             M4a = .2
             M0 = .83

        enginecruise = [
            aircraft.engine.engineP['M_2'][Nclimb:] == cruise['M'],
            aircraft.engine.engineP['M_{2.5}'][Nclimb:] == M25,
            aircraft.engine.engineP['hold_{2}'][Nclimb:] == 1.+.5*(1.398-1.)*M2**2.,
            aircraft.engine.engineP['hold_{2.5}'][Nclimb:] == 1.+.5*(1.354-1.)*M25**2.,
            
            # Thrust >= Drag + Vertical Component of Weight
            cruise['D'] + cruise['W_{avg}'] * cruise['\\theta'] <= aircraft['n_{eng}'] * aircraft.engine['F_{spec}'][Nclimb:],
            ]

        if BLI or M072_737 or b737800 or b777300ER or optimal737 or D8_eng_wing or D8_no_BLI:
             with SignomialsEnabled():
                  engineclimb.extend([
                       SignomialEquality(aircraft.engine.engineP['c1'][:Nclimb], (1. + 0.5*(.401)*climb['M']**2.)),
                       ])
                  enginecruise.extend([
                       SignomialEquality(aircraft.engine.engineP['c1'][Nclimb:], (1. + 0.5*(.401)*cruise['M']**2.)),                
                       ])
        else:
             engineclimb.extend([
                  aircraft.engine.engineP['c1'][:Nclimb] <= 1. + 0.5*(.401)*0.54592**2.,
                  ])
             enginecruise.extend([
                  aircraft.engine.engineP['c1'][Nclimb:] <= 1. + 0.5*(.401)*0.8**2.,
                  ])

        if fuel:
             # Fuel burn cost model
             if not multimission:
                  self.cost = aircraft['W_{f_{total}}']
                  self.cost = self.cost.sum()
             else:
                  self.cost = W_fmissions

             return constraints, aircraft, climb, cruise, enginestate, statelinking, engineclimb, enginecruise
             
        if operator:
             # Operator cost model
             if not multimission:
                  self.cost = aircraft['W_{dry}'] + aircraft['W_{f_{total}}']
                  self.cost = self.cost.sum()
             else:
                  self.cost = aircraft['W_{dry}'] + W_fmissions

             return constraints, aircraft, climb, cruise, enginestate, statelinking, engineclimb, enginecruise

        if manufacturer:
             # Manufacturer cost model
             if not multimission:
                  self.cost = aircraft['W_{dry}'] + aircraft['W_{f_{total}}']
                  self.cost = self.cost.sum()
             else:
                  self.cost = aircraft['W_{dry}'] + W_fmissions

             return constraints, aircraft, climb, cruise, enginestate, statelinking, engineclimb, enginecruise

        if PRFC:
             # Payload-range fuel consumption optimization - CHOOSES THE OPTIMAL MISSION, DO NOT NEED TO SUB ReqRng OR n_{pax}.
             if not multimission:
                self.cost = sum(aircraft['PRFC'])
             else:
                self.cost = sum(aircraft['PRFC'])
             return constraints, aircraft, climb, cruise, enginestate, statelinking, engineclimb, enginecruise


