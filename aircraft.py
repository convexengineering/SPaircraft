"""Simple commercial aircraft flight profile and aircraft model
Integrates Wing, VerticalTail, HorizontalTail , Fuselage, and Landing Gear models """
from __future__ import absolute_import

from builtins import range
import numpy as np
from gpkit import Variable, Model, units, SignomialsEnabled, SignomialEquality, Vectorize
from gpkit.constraints.tight import Tight as TCS
from gpkit.constraints.bounded import Bounded as BCS
from numpy import pi

TCS.reltol = 1e-3

# importing from D8_integration
from stand_alone_simple_profile import FlightState
from vertical_tail import VerticalTail
from horizontal_tail import HorizontalTail
from wing import Wing
from turbofan.engine_validation import Engine
from fuselage import Fuselage
from landing_gear import LandingGear

# for ESP
from collections import OrderedDict

"""
Models required to minimize the aircraft total fuel weight.

Sources for substitutions and equations:
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

    SKIP VERIFICATION

    ARGUMENTS
    ---------
    BLI: True = have engine stagnation pressure drop drom BLI, False = no engine stagnation pressure drop
    fitDrag: True = use Martin's tail drag fits, False = use the TASOPT tail drag model
    """

    def setup(self, Nclimb, Ncruise, flightstate, eng, fitDrag, BLI = False, Nmissions=0,  **kwargs):
        # create submodels
        self.fuse = Fuselage(Nmissions)
        self.wing = Wing()
        if Nmissions != 0:
            self.engine = Engine(True, Nclimb+Ncruise, flightstate, eng, Nmissions, BLI)
        else:
           self.engine = Engine(True, Nclimb+Ncruise, flightstate, eng, BLI)
        self.VT = VerticalTail()
        self.HT = HorizontalTail()
        self.LG = LandingGear()

        #set the tail drag flag
        self.fitDrag = fitDrag

        # variable definitions
        numaisle = Variable('n_{aisle}','-','Number of Aisles')
        numeng = Variable('n_{eng}', '-', 'Number of Engines')
        numVT = Variable('n_{vt}','-','Number of Vertical Tails')
        Vne = Variable('V_{ne}',143.92,'m/s', 'Never-exceed speed')  # [Philippe]
        Vmn = Variable('V_{mn}', 'm/s','Maneuvering speed')
        rhoTO = Variable('\\rho_{T/O}',1.225,'kg*m^-3','Air density at takeoff')
        ReserveFraction = Variable('f_{fuel_{res}}', '-', 'Fuel Reserve Fraction')

        #fraction of fuel in wings
        f_wingfuel = Variable('f_{wingfuel}', '-', 'Fraction of fuel stored in wing tanks')

        SMmin = Variable('SM_{min}','-', 'Minimum Static Margin')
        dxCG = Variable('\\Delta x_{CG}', 'm', 'Max CG Travel Range')
        xCGmin = Variable('x_{CG_{min}}','m','Maximum Forward CG')

        with Vectorize(Nmissions):
             Izwing = Variable('I_{z_{wing}}','kg*m**2','Wing moment of inertia')
             Iztail = Variable('I_{z_{tail}}','kg*m**2','Tail moment of inertia')
             Izfuse = Variable('I_{z_{fuse}}','kg*m**2','Fuselage moment of inertia')
             Iz = Variable('I_{z}', 'kg*m^2', 'Aircraft Z-axis Moment of Inertia')

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

        W_totalmax = Variable('W_{total, max}', 'lbf', 'Total Aircraft Weight')

        # Fuselage lift fraction variables
        Ltow = Variable('f_{L_{total/wing}}','-','Total lift as a percentage of wing lift')

        # Misc system variables
        Wmisc   = Variable('W_{misc}','lbf','Sum of Miscellaneous Weights')
        Whpesys = Variable('W_{hpesys}','lbf','Power Systems Weight')
        fhpesys = Variable('f_{hpesys}','-','Power Systems Weight Fraction')
        xmisc   = Variable('x_{misc}','m','Misc Weight Centroid')
        xhpesys = Variable('x_{hpesys}','m','Power Systems Weight x-Location')

        # Engine system variables
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
        xeng = Variable('x_{eng}', 'm', 'Engine x-location')

        Ceng = Variable('C_{engsys}', 1, '-', 'Engine System Weight Margin/Sens Factor')

        #BLI total drag reduction factor
        Dreduct = Variable('D_{reduct}', '-', 'BLI Drag Reduction Factor')
        Dwakefrac = Variable('D_{wakefraction}', 0.33, '-', 'Percent of Total Drag From Wake Dissipation')
        BLI_wake_benefit = Variable('BLI_{wakebenefit}', 0.02, '-', 'Wake Drag Reduction from BLI Wake Ingestion')

        # Design parameters for .csm generation
        self.design_parameters = OrderedDict([
            # Wing descriptors
            ('wing:b',             self.wing['b']),
            ('wing:croot',         self.wing['c_{root}']),
            ('wing:ctip',          self.wing['c_{tip}']),
            ('wing:S',             self.wing['S']),
            ('wing:xwing',         self.fuse['x_{wing}']),
            ('wing:tanw',          self.wing['\\tan(\\Lambda)']),
            # Fuselage descriptors
            ('fuselage:hfloor',    self.fuse['h_{floor}']),
            ('fuselage:lnose',     self.fuse['l_{nose}']),
            ('fuselage:lshell',    self.fuse['l_{shell}']),
            ('fuselage:lcone',     self.fuse['l_{cone}']),
            ('fuselage:lfloor',    self.fuse['l_{floor}']),
            ('fuselage:lfuse',     self.fuse['l_{fuse}']),
            ('fuselage:hfuse',     self.fuse['h_{fuse}']),
            ('fuselage:wfuse',     self.fuse['w_{fuse}']),
            ('fuselage:wfloor',    self.fuse['w_{floor}']),
            ('fuselage:wdb',       self.fuse['w_{db}']),
            ('fuselage:Rfuse',     self.fuse['R_{fuse}']),
            ('fuselage:dRfuse',    self.fuse['\\Delta R_{fuse}']),
            # Horizontal tail descriptors
            ('ht:xCGht',           self.HT['x_{CG_{ht}}']),
            ('ht:crootht',         self.HT['c_{root_{ht}}']),
            ('ht:ctipht',          self.HT['c_{tip_{ht}}']),
            ('ht:bht',             self.HT['b_{ht}']),
            ('ht:lht',             self.HT['l_{ht}']),
            ('ht:tanht',           self.HT['\\tan(\Lambda_{ht})']),
            # Vertical tail descriptors
            ('vt:xCGvt',           self.VT['x_{CG_{vt}}']),
            ('vt:Svt',             self.VT['S_{vt}']),
            ('vt:bvt',             self.VT['b_{vt}']),
            ('vt',                 self.VT['l_{vt}']),
            ('vt:crootvt',         self.VT['c_{root_{vt}}']),
            ('vt:ctipvt',          self.VT['c_{tip_{vt}}']),
            ('vt:tanvt',           self.VT['\\tan(\Lambda_{vt})']),
            # Engine descriptors
            ('engine:df',          self.engine['d_{f}']),
            ('engine:lnace',       lnace),
            ('engine:yeng',        self.VT['y_{eng}']),
            ('engine:xeng',        xeng),
        ])

        constraints = []
        with SignomialsEnabled():
            constraints.extend([
                            #varaible linking
                            self.wing['c_{root}'] == self.fuse['c_0'],
                            self.wing.wb['r_{w/c}'] == self.fuse['r_{w/c}'],
                            self.wing['x_w'] == self.fuse['x_{wing}'],
                            # Load factor matching
                            self.fuse['N_{lift}'] == self.wing['N_{lift}'], # To make sure that the loads factors match.
                            Ltow*self.wing['L_{max}'] >= self.wing['N_{lift}'] * W_totalmax + self.HT['L_{ht_{max}}'],

                            ## ----------------- WEIGHT BUILD UPS --------------
                            #compute the aircraft's zero fuel weight
                            TCS([self.fuse['W_{fuse}'] + numeng \
                                * Wengsys + self.fuse['W_{tail}'] + self.wing['W_{wing}'] + Wmisc <= W_dry]),

                            # Total takeoff weight constraint
                            TCS([W_ftotal + W_dry + self.fuse['W_{payload}'] <= W_total]),
                            TCS([W_ftotal + W_dry + self.fuse['W_{payload}'] <= W_total]),
                            TCS([W_ftotal >= W_fprimary + ReserveFraction * W_fprimary]),
                            TCS([W_fprimary >= W_fclimb + W_fcruise]),
                            W_totalmax >= W_total,
                            ## NEEDED FOR MULTIMISSION TODO: add if statement
##                            W_totalmax >= W_total[0],
##                            W_totalmax >= W_total[1],


                            ## ---------------- WING CONSTRAINTS -------------
                            # Wing fuel constraints
                            self.wing['W_{fuel_{wing}}'] >= f_wingfuel*W_ftotal/self.wing['FuelFrac'],


                            ## --------------------LANDING GEAR and POWER SYSTEM ----------------
                            # LG and Power Systems weights
                            Wmisc >= self.LG['W_{lg}'] + Whpesys,
                            Whpesys == fhpesys*W_totalmax,

                            # LG and Power System locations
                            self.LG['x_n'] <= self.fuse['l_{nose}'],
                            TCS([self.LG['x_m'] >= self.fuse['x_{wing}']]),
                            self.LG['x_m'] <= self.wing['\\Delta x_{AC_{wing}}'] + self.fuse['x_{wing}'],
                            xhpesys == 1.1*self.fuse['l_{nose}'],
                            xmisc*Wmisc >= xhpesys*Whpesys,

                            #compute nacelle diameter
                            self.LG['d_{nacelle}'] >= self.engine['d_{f}'] + 2*self.LG['t_{nacelle}'],
                            self.LG['d_{nacelle}'] <= self.wing['b'], # upper bounding

                           # Hard landing
                           # http://www.boeing.com/commercial/aeromagazine/...
                           # articles/qtr_3_07/AERO_Q307_article3.pdf
                           # sink rate of 10 feet per second at the maximum
                           # design landing weight
                           # Landing condition from Torenbeek p360
                           self.LG['E_{land}'] >= W_totalmax/(2*self.LG['g'])*self.LG['w_{ult}']**2, # Torenbeek (10-26)

                            #setting fuselage upsweep location
                            self.LG['x_{up}'] == self.fuse['x_{shell2}'],

                           # Maximum static loads through main and nose gears
                           self.LG['L_n'] == W_totalmax*self.LG['\\Delta x_m']/self.LG['B'],
                           self.LG['L_m'] == W_totalmax*self.LG['\\Delta x_n']/self.LG['B'],

                            # (assumes deceleration of 10 ft/s^2)
                            self.LG['L_{n_{dyn}}'] >= 0.31*((self.LG['z_{CG}']+self.LG['l_m'])/self.LG['B'])*W_totalmax,
                                                         self.VT['y_{eng}'] >= self.LG['y_m'],

                            ## ------------------- FUSELAGE CONSTRAINTS ---------------
                            # Tail cone sizing
                            3. * (numVT*self.VT['M_r']) * self.VT['c_{root_{vt}}'] * \
                                (self.fuse['p_{\\lambda_{vt}}'] - 1.) >= numVT*self.VT[
                                    'L_{vt_{max}}'] * self.VT['b_{vt}'] * (self.fuse['p_{\\lambda_{vt}}']),
                            TCS([self.fuse['V_{cone}'] * (1. + self.fuse['\\lambda_{cone}']) * \
                             (pi + 4. * self.fuse['\\theta_{db}']) >= numVT*self.VT[
                                'M_r'] * self.VT['c_{root_{vt}}'] / self.fuse['\\tau_{cone}'] * \
                                 (pi + 2. * self.fuse['\\theta_{db}']) * \
                                  (self.fuse['l_{cone}'] / self.fuse['R_{fuse}'])]), #[SP]

                            # Tail weight
                            self.fuse['W_{tail}'] >= numVT*self.VT['W_{vt}'] + self.HT['W_{ht}'],

                            # Fuselage width (numaisle comes in)
                            TCS([2.*self.fuse['w_{fuse}'] >= self.fuse['SPR'] * self.fuse['w_{seat}'] + \
                                 numaisle*self.fuse['w_{aisle}'] + 2. * self.fuse['w_{sys}'] + self.fuse['t_{db}']]),

                            # Vertical bending material coefficient (VT aero loads)
                            self.fuse['B_{1v}'] == self.fuse['r_{M_v}']*numVT*self.VT['L_{vt_{max}}']/(self.fuse['w_{fuse}']*self.fuse['\\sigma_{M_v}']),


                            ## ------------------ HORIZONTAL TAIL -----------------
                            # Lift curve slope ratio for HT and Wing
                            SignomialEquality(self.HT['m_{ratio}']*(1+2/self.wing['AR']), 1 + 2/self.HT['AR_{ht}']),

                            # HT Location and Volume Coefficient
                            self.HT['x_{CG_{ht}}'] <= self.fuse['l_{fuse}'],
                            TCS([self.HT['V_{ht}'] == self.HT['S_{ht}']*self.HT['l_{ht}']/(self.wing['S']*self.wing['mac'])]),

                            # HT Max Loading
                            TCS([self.HT['L_{ht_{max}}'] >= 0.5*rhoTO*Vne**2*self.HT['S_{ht}']*self.HT['C_{L_{ht,max}}']]),


                            ## ------------------- VERTICAL TAIL -----------------
                            # VT Max Loading
                            TCS([self.VT['L_{vt_{max}}'] >= 0.5*rhoTO*Vne**2*self.VT['S_{vt}']*self.VT['C_{L_{vt,max}}']]),

                            #VT CG location
                            self.VT['x_{CG_{vt}}'] <= self.fuse['l_{fuse}'],

                            # VT volume coefficient
                            self.VT['V_{vt}'] == numVT*self.VT['S_{vt}'] * self.VT['l_{vt}']/(self.wing['S']*self.wing['b']),

                            # VT sizing constraints
                            # Yaw rate constraint at flare
                            numVT*.5*self.VT['\\rho_{TO}']*self.VT['V_{land}']**2*self.VT['S_{vt}']*self.VT['l_{vt}']* \
                                            self.VT['C_{L_{vt,yaw}}'] >= self.VT['\\dot{r}_{req}']*self.VT['I_{z, max}'],

                            # Force moment balance for one engine out condition
                            # TASOPT 2.0 p45
                            numVT*self.VT['L_{vt,EO}']*self.VT['l_{vt}'] >= self.VT['T_e']*self.VT['y_{eng}'] + \
                                        self.VT['D_{wm}']*self.VT['y_{eng}'],

                            # Drag of a windmilling engine (VT sizing)
                            TCS([self.VT['D_{wm}'] >= 0.5*self.VT['\\rho_{TO}']*self.VT['V_1']**2.*self.engine['A_{2}']*self.VT['C_{D_{wm}}']]),


                            ## ------------- MOMENT OF INERTIA ------------
                            # Moment of inertia around z-axis
                            Iz >= Izwing + Iztail + Izfuse,
                            self.VT['I_{z, max}'] >= Iz,


                            ## --------------ENGINE SYSTEM---------------F-
                            #engine system weight constraints, nacelle dimensions
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

        ### -------------- ENGINE LOCATION RELATED CONSTRAINTS -----------------
        # Wing-engined aircraft constraints
        if wingengine:
            with SignomialsEnabled():
                constraints.extend([
                    # Wing root moment constraint, with wing and engine weight load relief
                    TCS([self.wing['M_r']*self.wing['c_{root}'] >= (self.wing['L_{max}'] - self.wing['N_{lift}']*(self.wing['W_{wing}']+f_wingfuel*W_ftotal)) * \
                       (self.wing['b']**2/(12*self.wing['S'])*(self.wing['c_{root}'] + 2*self.wing['c_{tip}'])) - \
                                        self.wing['N_{lift}']*Wengsys*self.VT['y_{eng}']]), #[SP]

                    # Horizontal tail aero+landing loads constants A1h
                    self.fuse['A_{1h_{Land}}'] >= (self.fuse['N_{land}'] * \
                                (self.fuse['W_{tail}'] + self.fuse['W_{apu}'])) / \
                                 (self.fuse['h_{fuse}'] * self.fuse['\\sigma_{bend}']),

                    self.fuse['A_{1h_{MLF}}'] >= (self.fuse['N_{lift}'] * \
                                (self.fuse['W_{tail}'] + self.fuse['W_{apu}']) \
                                + self.fuse['r_{M_h}'] * self.HT['L_{ht_{max}}']) / \
                                 (self.fuse['h_{fuse}'] * self.fuse['\\sigma_{M_h}']),

                    # Moment of inertia constraints
                    Izwing >= numeng*Wengsys*self.VT['y_{eng}']**2./g + \
                                    (self.wing['W_{fuel_{wing}}'] + self.wing['W_{wing}'])/(self.wing['S']*g)* \
                                    self.wing['c_{root}']*self.wing['b']**3.*(1./12.-(1.-self.wing['\\lambda'])/16.), #[SP]
                    Iztail >= (self.fuse['W_{apu}'] + self.VT['W_{vt}']) * self.VT['l_{vt}'] ** 2. / g + \
                        self.HT['W_{ht}'] * self.HT['l_{ht}'] ** 2. / g,
                    Izfuse >= (self.fuse['W_{fuse}'] + self.fuse['W_{payload_{max}}'])/self.fuse['l_{fuse}'] * \
                                    (self.fuse['x_{wing}']**3 + self.VT['l_{vt}']**3.)/(3.*g),

                    # Engine ground clearance
                    self.LG['d_{nacelle}']  + self.LG['h_{nacelle}'] <= self.LG['l_m'] + (self.VT['y_{eng}']-self.LG['y_m'])*self.LG['\\tan(\\gamma)'], # [SP]

                    # Engine x-location (weight centroid, roughly)
                    xeng >= self.fuse['x_{wing}'] + self.wing['\\tan(\\Lambda)']*self.VT['y_{eng}'] - 0.5*lnace,
                ])

        # Rear-engined aircraft constraints
        if rearengine:
            with SignomialsEnabled():
                constraints.extend([
                    # Wing root moment constraint, with wing weight + fuel load relief
                    TCS([self.wing['M_r']*self.wing['c_{root}'] >= (self.wing['L_{max}'] - self.wing['N_{lift}'] * (self.wing['W_{wing}']+f_wingfuel*W_ftotal)) * \
                        (self.wing['b']**2/(12*self.wing['S'])*(self.wing['c_{root}'] + 2*self.wing['c_{tip}']))]), #[SP]

                    # Horizontal tail aero+landing loads constants A1h
                    self.fuse['A_{1h_{Land}}'] >= (self.fuse['N_{land}'] * \
                                                (self.fuse['W_{tail}'] + numeng * Wengsys + self.fuse['W_{apu}'])) / \
                                                (self.fuse['h_{fuse}'] * self.fuse['\\sigma_{bend}']),

                    self.fuse['A_{1h_{MLF}}'] >= (self.fuse['N_{lift}'] * \
                                               (self.fuse['W_{tail}'] + numeng * Wengsys + self.fuse['W_{apu}']) \
                                               + self.fuse['r_{M_h}'] * self.HT['L_{ht_{max}}']) / \
                                                (self.fuse['h_{fuse}'] * self.fuse['\\sigma_{M_h}']),

                    # Moment of inertia constraints
                    Izwing >= (self.wing['W_{fuel_{wing}}'] + self.wing['W_{wing}']) / (self.wing['S'] * g) * \
                    self.wing['c_{root}'] * self.wing['b'] ** 3. * (1. / 12. - (1. - self.wing['\\lambda']) / 16.),
                    # [SP]
                    Iztail >= (self.fuse['W_{apu}'] + self.VT['W_{vt}'] + numeng * Wengsys) * \
                        self.VT['l_{vt}'] ** 2. / g + \
                        self.HT['W_{ht}'] * self.HT['l_{ht}'] ** 2. / g,
                    # NOTE: Using l_{vt} as an xeng - xCG surrogate. Reason: xCG moves during flight; want scalar Izfuse
                    Izfuse >= (self.fuse['W_{fuse}'] + self.fuse['W_{payload_{max}}']) / self.fuse['l_{fuse}'] * \
                    (self.fuse['x_{wing}'] ** 3. + self.VT['l_{vt}'] ** 3.) / (3. * g),
                    # Note: Using x_{wing} as a CG surrogate. Want scalar Izfuse.

                    # Engine x-location (weight centroid, roughly)
                    xeng <= self.fuse['x_{shell2}'] + 1.00*self.fuse['l_{cone}'],
                    xeng >= self.fuse['x_{shell2}'] + 0.75*self.fuse['l_{cone}'],
                ])

        ### -------------- FUSELAGE CONSTRAINTS ----------------
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


        ### ---------------- HORIZONTAL TAIL CONSTRAINTS ------------------
        # Pi HT constraints:
        if piHT:
            with SignomialsEnabled():
                constraints.extend([
                    # Pin VT joint moment constraint #PROBLEMATIC, instead using wingtip moment
                    # SignomialEquality(self.HT['L_{ht_{rect}}'] * (self.HT['b_{ht}'] / 2. - self.fuse['w_{fuse}']),
                    #                   self.HT['L_{ht_{tri}}'] * (self.fuse['w_{fuse}'] - self.HT['b_{ht}'] / 3.)), # [SP] #[SPEquality]
                    # Pin VT constraint (wingtip moment = 0Nm) #TODO: may be problematic as well, relax if doesn't solve
                    SignomialEquality(self.HT['b_{ht}']/4.*self.HT['L_{ht_{rect}}'] + self.HT['b_{ht}']/3.*self.HT['L_{ht_{tri}}'],
                                      self.HT['b_{ht_{out}}'] * self.HT['L_{ht_{max}}']/2.), #[SP] #[SPEquality]

                    # HT outboard half-span
                    SignomialEquality(self.HT['b_{ht_{out}}'] , 0.5*self.HT['b_{ht}'] - self.fuse['w_{fuse}']), #[SP] #[SPEquality]

                    # HT center moment
                    self.HT['M_r'] * self.HT['c_{root_{ht}}'] >= self.HT['L_{ht_{rect}}'] * (
                    self.HT['b_{ht}'] / 4.) + self.HT['L_{ht_{tri}}'] * (self.HT['b_{ht}'] / 6.) - \
                    self.fuse['w_{fuse}'] * self.HT['L_{ht_{max}}'] / 2., # [SP]

                    # HT joint moment
                    self.HT['M_{r_{out}}']*self.HT['c_{attach}'] >= self.HT['L_{ht_{rect_{out}}}'] * (0.5*self.HT['b_{ht_{out}}']) + \
                                                                    self.HT['L_{ht_{tri_{out}}}'] * (1./3.*self.HT['b_{ht_{out}}']),
                    self.HT['M_{r_{out}}'] <= 1e20*units('N'), # upper bounding

                    # HT joint shear (max shear)
                    self.HT['L_{shear}'] >= self.HT['L_{ht_{rect_{out}}}'] + self.HT['L_{ht_{tri_{out}}}'],

                    # HT/VT joint constraint
                    SignomialEquality(self.HT['c_{tip_{ht}}'] + (1. - self.HT['\\lambda_{ht}']) * 2. * self.HT['b_{ht_{out}}'] / self.HT['b_{ht}'] *
                        self.HT['c_{root_{ht}}'],
                    self.HT['c_{attach}']),

                    # HT structural factor calculation
                    self.HT['\\pi_{M-fac}'] >= (0.5*(self.HT['M_{r_{out}}']*self.HT['c_{attach}'] + \
                                                    self.HT['M_r']* self.HT['c_{root_{ht}}']) * self.fuse['w_{fuse}'] / \
                                                    (0.5*self.HT['M_{r_{out}}']*self.HT['c_{attach}']*self.HT['b_{ht_{out}}']) + 1.0) * \
                                                        self.HT['b_{ht_{out}}'] / (0.5*self.HT['b_{ht}']),
                ])
        # Conventional HT constraints
        if conventional:
            with SignomialsEnabled():
                constraints.extend([
                    # HT root moment
                    TCS([self.HT['M_r']*self.HT['c_{attach}'] >= 1./3.*self.HT['L_{ht_{tri_{out}}}']*self.HT['b_{ht_{out}}'] + \
                         1./2.*self.HT['L_{ht_{rect_{out}}}']*self.HT['b_{ht_{out}}']]),
                    # HT joint constraint
                    self.HT['c_{attach}'] == self.HT['c_{root_{ht}}'],

                    # HT auxiliary variables
                    self.HT['b_{ht_{out}}'] == 0.5*self.HT['b_{ht}'],
                    self.HT['M_{r_{out}}'] == self.HT['M_r'],
                    self.HT['L_{shear}'] >= self.HT['L_{ht_{rect_{out}}}'] + self.HT['L_{ht_{tri_{out}}}'],

                    # HT structural factor calculation
                    self.HT['\\pi_{M-fac}'] == 1.0,
                ])

        self.components = [self.fuse, self.wing, self.engine, self.VT, self.HT, self.LG]

        return self.components, constraints

    def flight_dynamic(self, state, Nclimb, Ncruise): # creates an aircraft flight performance model, given a state
        return FlightP(self, state, Nclimb, Ncruise)


class AircraftP(Model):
    """
    Aircraft performance models superclass, contains constraints true for
    all flight segments
    SKIP VERIFICATION
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
        WLoadmax = Variable('W_{Load_{max}}',6664., 'N/m^2', 'Max Wing Loading')
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
        Cfnace = Variable('C_{f_{nacelle}}', '-', 'Nacelle Drag Coefficient')
        Renace = Variable('R_{e_{nacelle}}', '-', 'Nacelle Reynolds Number')
        Cfturb = Variable('C_{f_{nacelle}}', '-', 'Turbulent Nacelle Skin Friction Coefficient')
        Cdnace = Variable('C_{d_{nacelle}}', '-', 'Nacelle Drag Coeffecient')
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

            # Geometric average of start and end weights of flight segment
            W_avg >= (W_start * W_end)**.5 + W_buoy, # Buoyancy weight included in Breguet Range

            # Flight time unit conversion
            t == thours,

            # Fuselage lift (just calculating)
            SignomialEquality(self.fuseP['L_{fuse}'], (self.aircraft['f_{L_{total/wing}}']-1.)*self.wingP['L_w']),

            # Tail downforce penalty to total lift
            TCS([Ltotal == self.aircraft['f_{L_{total/wing}}']*self.wingP['L_w']]),
            TCS([Ltotal >= W_avg + self.HTP['L_{ht}']]),

            ## ----------------- DRAG CONSTRAINTS --------------
            self.fuseP['D_{fuse}'] == 0.5 * state['\\rho'] * state['V']**2 * \
                                        self.fuseP['C_{D_{fuse}}'] * aircraft['l_{fuse}'] * aircraft['R_{fuse}'] * (state['M']**2/aircraft.fuse['M_{fuseD}']**2),
            D >= aircraft['D_{reduct}'] * (self.wingP['D_{wing}'] + self.fuseP['D_{fuse}'] + self.aircraft['n_{vt}']*self.VTP['D_{vt}'] + self.HTP['D_{ht}'] + aircraft['n_{eng}'] * Dnace),
            C_D == D/(.5*state['\\rho']*state['V']**2 * self.aircraft.wing['S']),
            LoD == W_avg/D,

            ## ----------------- WING CONSTARINTS -------------
            # Maximum wing loading constraint
            WLoad <= WLoadmax,

            # Wing loading
            WLoad == .5 * self.wingP['C_{L}'] * self.aircraft['S'] * state.atm['\\rho'] * state['V']**2 / self.aircraft.wing['S'],

            # Center wing lift loss
            # self.wingP['p_{o}'] >= self.wingP['L_w']*self.wing['c_{root}']*(.5 + 0.5*self.wingP['\\eta_{o}'](/(self.wing['S']),
            self.wingP['p_{o}'] >= self.wingP['L_w']*aircraft.wing['c_{root}']/(aircraft.wing['S']), #TODO improve approx without making SP
            self.wingP['\\eta_{o}'] == aircraft['w_{fuse}']/(aircraft['b']/2),

            # Wing location and AC constraints
            TCS([xAC <= aircraft['x_{wing}'] + 0.25*aircraft['\\Delta x_{AC_{wing}}'] + xNP]), #[SP] #TODO relax and improve

            # Static margin constraints
            self.wingP['c_{m_{w}}'] == 1.9,

            # Neutral point approximation (taken from Basic Aircraft Design Rules, Unified)
            # TODO improve
            SignomialEquality(xNP/aircraft['mac']/aircraft['V_{ht}']*(aircraft['AR']+2.)*(1.+2./aircraft['AR_{ht}']),
                              (1.+2./aircraft['AR'])*(aircraft['AR']-2.)),


            ## ------------- VERTICAL TAIL CONSTRAINTS -----------
            # VT TE constraint, and CG calculation
            xCG + aircraft.VT['\\Delta x_{trail_{vt}}'] <= aircraft.fuse['l_{fuse}'],
            aircraft.VT['x_{CG_{vt}}'] >= xCG +0.5*(aircraft.VT['\\Delta x_{lead_{vt}}']+aircraft.VT['\\Delta x_{trail_{vt}}']),


            ## --------------- HORIZONTAL LOCATION GEOMETRY AND PERFORMANCE -----------
            # HT CG calculation
            aircraft.HT['x_{CG_{ht}}'] >= xCG +0.5*(aircraft.HT['\\Delta x_{lead_{ht}}']+aircraft.HT['\\Delta x_{trail_{ht}}']),

            # HT lift coefficient calc
            self.HTP['C_{L_{\\alpha,ht}}'] + (2*self.wingP['C_{L_{\\alpha,w}}']/(pi*aircraft.wing['AR']))*aircraft.HT['\\eta_{ht}']*self.HTP['C_{L_{\\alpha,ht_0}}'] <= self.HTP['C_{L_{\\alpha,ht_0}}']*aircraft.HT['\\eta_{ht}'],

            # Tail aspect ratio and lift constraints
            aircraft.HT['AR_{ht}'] >= 4., #TODO change to tip Re constraint
            self.HTP['C_{L_{ht}}'] >= 0.01, #TODO remove

            ## -------------- HORIZONTAL TAIL SIZING AND STATIC MARGIN----------------
            TCS([SM <= (xAC-xCG)/aircraft['mac']]),
            SM >= aircraft['SM_{min}'],

            #min static margin at forward and aft CG locations
            TCS([aircraft['SM_{min}'] + aircraft['\\Delta x_{CG}']/aircraft.wing['mac'] \
                 + self.wingP['c_{m_{w}}']/aircraft.wing['C_{L_{w,max}}'] <= \
                                            aircraft.HT['V_{ht}']*aircraft.HT['m_{ratio}'] +\
                                            aircraft.HT['V_{ht}']*aircraft.HT['C_{L_{ht,max}}']/aircraft.wing['C_{L_{w,max}}']]), # [SP]

            # Aircraft trim conditions
            TCS([xAC/aircraft.wing['mac'] <= xCG/aircraft.wing['mac'] + \
                 self.wingP['c_{m_{w}}']/self.wingP['C_{L}']  +\
                              aircraft.HT['V_{ht}']*(self.HTP['C_{L_{ht}}']/self.wingP['C_{L}'])]),

            ## ------------- NACELLE DRAG CONSTRAINTS ----------------
            #nacelle drag
            Renace == state['\\rho']*state['V'] * aircraft['l_{nacelle}']/state['\\mu'],
            Cfnace == 0.94*4.*0.0743/(Renace**(0.2)), #from http://www.calpoly.edu/~kshollen/ME347/Handouts/Friction_Drag_Coef.pdf
            Vnace == aircraft['r_{v_{nacelle}}'] * state['V'],
            Vnacrat >= 2.*Vnace/state['V'] - V2/state['V'],
            rvnsurf**3. >= 0.25*(Vnacrat + aircraft['r_{v_{nacelle}}'])*(Vnacrat**2. + aircraft['r_{v_{nacelle}}']**2.),
            Cdnace == aircraft['f_{S_{nacelle}}'] * Cfnace[0] * rvnsurf **3.,
            Dnace == Cdnace * 0.5 * state['\\rho'] * state['V']**2. * aircraft['S'],
            ])

        ## ------------------ HORIZONTAL TAIL TRAILING EDGE -------------
        # HT TE constraint
        if conventional:
            constraints.extend([
            aircraft['l_{fuse}'] >= xCG + aircraft.HT['\\Delta x_{trail_{ht}}']])
        if piHT:
            with SignomialsEnabled():
                constraints.extend([
                    aircraft.HT['\\Delta x_{trail_{ht}}'] <= aircraft.VT['\\Delta x_{lead_{vt}}'] + \
                        aircraft['b_{vt}']/aircraft['\\tan(\\Lambda_{vt})'] + \
                        aircraft['w_{fuse}']/aircraft['\\tan(\\Lambda_{ht})'] + aircraft['c_{root_{ht}}']])

        # use the TASOPT tail drag model if fitDrag == False
        if not aircraft.fitDrag:
            constraints.extend([
                #set the VT drag coefficient
                self.VTP['C_{D_{vis}}'] >= (self.aircraft.VT['c_{d_{fv}}'] + self.aircraft.VT['c_{d_{pv}}']*self.aircraft.VT['\\cos(\\Lambda_{vt})^3']),

                #set the HT drag coefficient
                self.HTP['C_{D_{0,ht}}'] >= (self.aircraft.HT['c_{d_{fh}}'] + self.aircraft.HT['c_{d_{ph}}']*self.aircraft.HT['\\cos(\\Lambda_{ht})^3']),
                ])

        return self.Pmodels, constraints

class FlightP(Model): # Flight segment performance constraints
    "SKIP VERIFICATION"

    def setup(self, aircraft, state, Nclimb, Ncruise, **kwargs):
        # submodels
        self.aircraft = aircraft
        self.aircraftP = AircraftP(aircraft, state)
        self.wingP = self.aircraftP.wingP
        self.fuseP = self.aircraftP.fuseP
        self.engine = aircraft.engine

        # variable definitions
        theta = Variable('\\theta', '-', 'Aircraft Climb Angle')
        excessP = Variable('P_{excess}', 'W', 'Excess Power During Climb')
        RC = Variable('RC', 'feet/min', 'Rate of Climb')
        minRC = Variable('RC_{min}',500, 'feet/min', 'Minimum Rate of Climb')
        dhft = Variable(
            'dhft', 'feet', 'Change in Altitude Per Flight Segment [feet]')
        R = Variable('R_{segment}', 'nautical_miles',
                            'Down Range Covered in Each Segment')

        constraints = []
        constraints.extend([
            # Excess power for climb
            TCS([excessP + state['V'] * self.aircraftP['D'] <= state['V']
                 * aircraft['n_{eng}'] * self.engine['F']]),

            #compute climb rate
            RC == excessP / self.aircraftP['W_{avg}'],

            # Climb angle and rate constraint
            theta * state['V'] == RC,

            # compute total altitude change
            dhft == self.aircraftP['tmin'] * RC,

            # Small angle assumption during climb
            R == self.aircraftP['thr'] * state['V'],
        ])

        return constraints, self.aircraftP

class FlightSegment(Model):
    "SKIP VERIFICATION"
    def setup(self, aircraft, flightstate, Nclimb, Ncruise, **kwargs):
        self.state = flightstate
        self.flightP = aircraft.flight_dynamic(self.state, Nclimb, Ncruise)
        return self.state, self.flightP

class Mission(Model):
    """
    Mission superclass, links together all subclasses into an optimization problem
    SKIP VERIFICATION
    Inputs:
    Nclimb: number of climb segments (for Breguet Range)
    Ncruise: number of cruise segments (for Breguet Range)
    objective: defines the objective function
    aircraft: string representing the aircraft model
    Nmission: specifies whether single-point or multi-point optimization
              Nmission >/= 1 requires specification of range and number of passengers for each mission
    """

    def setup(self, Nclimb, Ncruise, config, Nmission = 1):
        # define global variables
        global wingengine, rearengine, doublebubble, tube, piHT, conventional
        global largeAC, multimission

        self.Nclimb = Nclimb
        self.Ncruise = Ncruise
        self.Nmission = Nmission

        # aircraft geometry flags
        wingengine = False; rearengine = False; doublebubble = False; tube = False;
        piHT = False; conventional = False; BLI = False; eng = 0;

        # Aircraft type, only one active at once
        largeAC = False

        # set geometry flags based on config type
        if config == 'D8_eng_wing':
            wingengine = True; piHT = True; doublebubble = True; eng = 3; BLI = False;
        if config == 'optimal737':
            conventional = True; eng = 3; BLI = False;
        if config == 'optimalD8':
            optimalD8 = True; rearengine = True; piHT = True; doublebubble = True; eng = 3; BLI = True;
        if config == 'optimal777':
            largeAC = True; conventional = True; eng = 4; BLI = False;
        if config == 'M072_737':
            conventional = True; eng = 3; BLI = False;
        if config == 'D8_no_BLI':
            rearengine = True; piHT = True; doublebubble = True; eng = 3; BLI = False;

        # if conventional choose wing engine and tube fuselage
        if conventional:
            wingengine = True; tube = True;

        # Multimission?
        if Nmission == 1:
            multimission = False
        else:
            multimission = True

        # Defining fitDrag, boolean describing whether or not to use XFOIL tail drag fits
        # False uses TASOPT tail drag model. Currently on.
        fitDrag = True

        # vectorize
        with Vectorize(Nmission):
             with Vectorize(Nclimb + Ncruise):
                 self.flightstate = flightstate = FlightState()

        # Build required submodels
        self.aircraft = aircraft = Aircraft(Nclimb, Ncruise, flightstate, eng, fitDrag, BLI, Nmission)
        self.aircraft.config = config

        # Vectorize dynamic variables
        with Vectorize(Nmission):
             with Vectorize(Nclimb+Ncruise):
                 self.flight = flight = FlightSegment(aircraft, flightstate, Nclimb, Ncruise)

        # Declare Mission variables
        if multimission:
             with Vectorize(Nmission):
                  CruiseAlt = Variable('CruiseAlt', 'ft', 'Cruise Altitude [feet]')
                  ReqRng = Variable('R_{req}', 'nautical_miles', 'Required Cruise Range')
                  Total_Time = Variable('TotalTime', 'hr', 'Total Mission Time')
                  climb_time = Variable('ClimbTime', 'min', 'Total Time in Climb')
                  climb_distance = Variable('ClimbDistance', 'nautical_miles', 'Climb Distance')
        else:
          CruiseAlt = Variable('CruiseAlt', 'ft', 'Cruise Altitude [feet]')
          ReqRng = Variable('R_{req}', 'nautical_miles', 'Required Cruise Range')
          Total_Time = Variable('TotalTime', 'hr', 'Total Mission Time')
          climb_time = Variable('ClimbTime', 'min', 'Total Time in Climb')
          climb_distance = Variable('ClimbDistance', 'nautical_miles', 'Climb Distance')

        max_climb_time = Variable('MaxClimbTime', 'min', 'Total Time in Climb')
        max_climb_distance = Variable('MaxClimbDistance', 'nautical_miles', 'Climb Distance')
        CruiseTt41max = Variable('T_{t_{4.1_{max-Cruise}}}', 3000., 'K', 'Max Cruise Turbine Inlet Temp')
        MinCruiseAlt = Variable('MinCruiseAlt', 'ft', 'Minimum Cruise Altitude')
        Fsafetyfac = Variable('Fsafetyfac', '-', 'Safety factor on inital climb thrust')

        # make overall constraints
        constraints = []

        with SignomialsEnabled():
            ## -------------------- BUOYANCY CONSTRAINTS ----------------
            # Buoyancy weight #TODO relax the equality
            # SignomialEquality(W_buoy,(rhocabin - state['\\rho'])*g*aircraft['V_{cabin}']),  #[SP] #[SPEquality]
            # Note: Buoyancy model has been simplified, since it causes significant increases in runtime.
            constraints.extend([
                flight['W_{buoy}'] >= flight['\\rho_{cabin}']*g*aircraft['V_{cabin}'], # [SP] # - cruise['\\rho']
                aircraft['PRFC'] == aircraft['W_{f_{primary}}']/g*aircraft.engine['h_{f}']/(ReqRng*aircraft['W_{payload}'])
            ])

            ## ---------------------- CG CONSTRAINTS ----------------------
            #depends on engine location
            if rearengine:
                constraints.extend([
                TCS([flight['x_{CG}']*flight['W_{avg}'] >=
                    aircraft['x_{misc}']*aircraft['W_{misc}'] + aircraft['x_{CG_{lg}}']*aircraft['W_{lg}'] \
                    + 0.5*(aircraft.fuse['W_{fuse}']+aircraft.fuse['W_{payload}'])*aircraft.fuse['l_{fuse}'] \
                    + (aircraft['W_{ht}']*aircraft['x_{CG_{ht}}']) + (aircraft['W_{vt}'])*aircraft['x_{CG_{vt}}'] \
                    + aircraft['n_{eng}']*aircraft['W_{engsys}'] * aircraft['x_{eng}'] \
                    + (aircraft['W_{wing}']*(aircraft.fuse['x_{wing}']+aircraft.wing['\\Delta x_{AC_{wing}}'])) \
                    + (flight['F_{fuel}']+aircraft['f_{fuel_{res}}'])*aircraft['W_{f_{primary}}'] \
                    * (aircraft.fuse['x_{wing}']+aircraft.wing['\\Delta x_{AC_{wing}}']*flight['F_{fuel}']) \
                    ])
              ])
            if wingengine:
                constraints.extend([
                TCS([flight['x_{CG}']*flight['W_{avg}'] >=
                    aircraft['x_{misc}']*aircraft['W_{misc}']  + aircraft['x_{CG_{lg}}']*aircraft['W_{lg}'] \
                    + 0.5*(aircraft.fuse['W_{fuse}']+aircraft.fuse['W_{payload}'])*aircraft.fuse['l_{fuse}'] \
                    + (aircraft['W_{ht}']*aircraft['x_{CG_{ht}}'] + (aircraft['W_{vt}'])*aircraft['x_{CG_{vt}}'])  \
                    + (aircraft['W_{wing}']*(aircraft.fuse['x_{wing}']+aircraft.wing['\\Delta x_{AC_{wing}}'])) \
                    + (flight['F_{fuel}']+aircraft['f_{fuel_{res}}'])*aircraft['W_{f_{primary}}'] \
                    * (aircraft.fuse['x_{wing}']+aircraft.wing['\\Delta x_{AC_{wing}}']*flight['F_{fuel}']) \
                    + aircraft['n_{eng}']*aircraft['W_{engsys}']*aircraft['x_{eng}']])
              ])

            # ------------------ LG CG DISTANCE AND TIP OVER COMPUTATIONS ----------------------
            constraints.extend([
                TCS([aircraft['\\Delta x_n'] + aircraft['x_n'] >= flight['x_{CG}'][Nclimb]]),
                TCS([aircraft['\\Delta x_m'] + flight['x_{CG}'][Nclimb] >= aircraft['x_m']]),
                # Longitudinal tip over (static)
                aircraft['x_m'] >= aircraft['\\tan(\\phi)']*(aircraft['z_{CG}']+aircraft['l_m']) + flight['x_{CG}'][Nclimb],
                ])

            # ---------------------- FUSELAGE LIFT, BLI CORRECTION, AND DRAG ----------------
            if doublebubble:
                constraints.extend([
                    flight.flightP.fuseP['C_{D_{fuse}}'] == 0.018081,
                    aircraft.fuse['M_{fuseD}'] == 0.72,
                  ])

            # Option for a high speed double bubble fuselage
            # elif D12 or D8bigfam:
            #     constraints.extend([
            #             flight.flightP.fuseP['C_{D_{fuse}}'] == 0.0167620,
            #             aircraft.fuse['M_{fuseD}'] == 0.83,
            #     ])

            if conventional and not largeAC:
                constraints.extend([
                    #Setting fuselage drag coefficient
                    flight.flightP.fuseP['C_{D_{fuse}}'] == 0.01107365,
                    aircraft.fuse['M_{fuseD}'] == 0.80,
                ])
            elif largeAC:
                constraints.extend([
                    #Setting fuselage drag coefficient
                    #additioanl 1.1 factor accounts for mach drag rise model
                    flight.flightP.fuseP['C_{D_{fuse}}'] == 0.00987663,
                    aircraft.fuse['M_{fuseD}'] == 0.84,
                ])

        ## ------------------------ WEIGHT BUILD UP AND LINKING -------------------
        constraints.extend([
            flight.flightP.aircraftP['W_{start}'][0] == aircraft['W_{total}'],

            # Climb segment weight decreases by the fuel burn...
            TCS([flight.flightP.aircraftP['W_{start}'] >= flight.flightP.aircraftP[
                'W_{end}'] + flight.flightP.aircraftP['W_{burn}']]),

            flight.flightP.aircraftP['W_{start}'][1:] == flight.flightP.aircraftP['W_{end}'][:-1],

            TCS([aircraft['W_{dry}'] + aircraft['W_{payload}'] + \
                 aircraft['f_{fuel_{res}}'] * aircraft['W_{f_{primary}}'] <= flight.flightP.aircraftP['W_{end}'][-1]]),
            TCS([aircraft['W_{f_{climb}}'] >= sum(flight.flightP.aircraftP['W_{burn}'][:Nclimb])]),
            TCS([aircraft['W_{f_{cruise}}'] >= sum(flight.flightP.aircraftP['W_{burn}'][Nclimb:])]),
            ])

        with SignomialsEnabled():
            constraints.extend([
                ## ------------------------ FLIGHT SEGMENT ALTITUDE AND PERFORMANCE CONSTRAINTS ---------------
                # Altitude constraints
                flight['hft'][Nclimb-1] >= CruiseAlt,
                CruiseAlt >= 1e-5*units('m'), # lower bounding
                SignomialEquality(flight['hft'][1:Nclimb+Ncruise], flight['hft'][:Nclimb+Ncruise - 1] + flight['dhft'][1:Nclimb+Ncruise]), #[SP]
                TCS([flight['hft'][0] == flight['dhft'][0]]),

                # All climb segments have the same total altitude change
                flight['dhft'][1:Nclimb] == flight['dhft'][:Nclimb-1],

                # T/O minimum climb rate constraint
                flight['RC'][0] >= 2500. * units('ft/min'),

                # Overall minimum climb rate constraint
                flight['RC'][:Nclimb] >= flight['RC_{min}'][:Nclimb],

                # TASOPT TOC climb rate constraint
                flight['\\theta'][Nclimb-1] >= 0.015, #higher than 0.015 radian climb gradient at top-of-climb

                # cruise ends at or above min cruise altitude
                flight['hft'][Nclimb-1] >= MinCruiseAlt,

                # Thrust >= Drag + Vertical Potential Energy #TODO: remove?
                aircraft['n_{eng}'] * aircraft.engine['F'] >= flight['D'] + flight['W_{avg}'] * flight['\\theta'],

                # All cruise segments cover the same range.
                flight['R_{segment}'][Nclimb:Nclimb+Ncruise-1] == flight['R_{segment}'][Nclimb+1:Nclimb+Ncruise],

                # Cruise Mach Number constraint
                flight['M'][Nclimb:Nclimb+Ncruise] >= aircraft['M_{min}'],

                ## ----------------------- CONSTRAIN TOTAL RANGE -----------
                TCS([sum(flight['R_{segment}']) >= ReqRng]), #[SP]

                ## ---------------------- VERTICAL TAIL SIZING ----------
                # Takeoff thrust T_e calculated for engine out + vertical tail sizing.
                # Note: coeff can be varied as desired...this exists because we
                # don't include a takeoff simulation and climb thrust < max takeoff
                # thrust. Future work is to add a balanced field length constraint
                # and remove this.
                aircraft.VT['T_e'] == Fsafetyfac * flight.flightP.engine['F'][0],

                ## -------------------- FUEL BURN --------------------
                # compute fuel burn from TSFC
                flight.flightP.aircraftP['W_{burn}'] == aircraft['n_{eng}'] * aircraft.engine['TSFC'] * \
                    flight['thr'] * aircraft.engine['F'],

                ## --------------------- HT AND VT GEOMETRY ----------------


                # ----------------- NACELLE DRAG ---------------
                # Elevated this constraint to Mission for dimensionality
                flight.flightP['V_2'] == aircraft.engine['M_2'] * flight.state['a'],

                ## ---------------------- SET WING MAX AOA -----------------
                flight['\\alpha_{max,w}'][:Nclimb] == .18,
                flight['\\alpha_{max,w}'][Nclimb:] == .1,


                ## -------------------- VARIOUS FLIGHT TIME COMPUTATIONS -------------------
                #compute the total time
                Total_Time >= sum(flight['thr']),
                Total_Time <= 1e10*units('s'), # upper bounding
                #compute the climb in time
                climb_time >= sum(flight['thr'][:Nclimb]),
                climb_time <= max_climb_time,

                ## --------------------- ENGINE CONSTRAINTS --------------------
                #set the max allowed cruise Tt4.1
                aircraft['T_{t_{4.1}}'][Nclimb:] <= CruiseTt41max,
                ])

        ## ------------------------ PERCENT FUEL REMAINING -------------------
        with SignomialsEnabled():
            for i in range(0, Nclimb+Ncruise):
                constraints.extend([
                    TCS([flight['F_{fuel}'][i] >= (sum(flight['W_{burn}'][i+1:]) + \
                                0.0000001*aircraft['W_{f_{primary}}'])/aircraft['W_{f_{primary}}']]),
                    flight['F_{fuel}'] <= 1.0000001, #just in case, TODO remove later
                    ])

        ## ---------------------- MULTIMISSION SETUP --------------------------
        if multimission:
            W_fmissions = Variable('W_{f_{missions}}', 'lbf', 'Fuel burn across all missions')
            constraints.extend([
                  W_fmissions >= sum(aircraft['W_{f_{total}}']),
                  ])

        ## -------------------- SETTING ENGINE PARAMETERS ----------------------
        constraints.extend([
            # constrain OPR less than max OPR
            aircraft['OPR'] <= aircraft.engine['OPR_{max}'],
            # bounding
            aircraft.engine['\\alpha_{max}'] <= 100., # bypass ratio
        ])

        M2 = .6
        M25 = .6

        engineclimb = [
            aircraft.engine.engineP['M_2'][:Nclimb] == flight['M'][:Nclimb],
            aircraft.engine.engineP['M_{2.5}'][:Nclimb] == M25,
            aircraft.engine.engineP['hold_{2}'][:Nclimb] == 1.+.5*(1.398-1.)*M2**2.,
            aircraft.engine.engineP['hold_{2.5}'][:Nclimb] == 1.+.5*(1.354-1.)*M25**2.,

            #climb rate constraints (TODO: remove?)
            TCS([flight['P_{excess}'] + flight.state['V'] * flight['D'] <= flight.state['V'] * aircraft['n_{eng}'] * aircraft.engine['F']]),
            ]


        if largeAC: #(or for larger aircraft)
             M2 = .65

        enginecruise = [
            aircraft.engine.engineP['M_2'][Nclimb:] == flight['M'][Nclimb:],
            aircraft.engine.engineP['M_{2.5}'][Nclimb:] == M25,
            aircraft.engine.engineP['hold_{2}'][Nclimb:] == 1.+.5*(1.398-1.)*M2**2.,
            aircraft.engine.engineP['hold_{2.5}'][Nclimb:] == 1.+.5*(1.354-1.)*M25**2.,

            ]

        with SignomialsEnabled():
            engineclimb.extend([
                       SignomialEquality(aircraft.engine.engineP['c1'], (1. + 0.5*(.401)*flight['M']**2.)),
                       ])

        return constraints, aircraft, flight, engineclimb, enginecruise
