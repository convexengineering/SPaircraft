"""Simple commercial aircraft flight profile and aircraft model"""
from numpy import pi
import numpy as np
from gpkit import Variable, Model, units, SignomialsEnabled, Vectorize
from gpkit.constraints.sigeq import SignomialEquality
from gpkit.constraints.tight import Tight as TCS
from wingbox import WingBox

class VerticalTail(Model):
    """
    Vertical tail sizing
    SKIP VERIFICATION

    References:
    1: LEAP engine specs (1B)
    2: TASOPT 737 code
    3: Boeing 737 airport doc
    4: Boeing 737 Max airport doc
    5: http://www.digitaldutch.com/atmoscalc
    6: Engineering toolbox
    7: Boeing.com
    """
    def setup(self, **kwargs):
        self.vtns = VerticalTailNoStruct()
        self.wb = WingBox(self.vtns, "vertical_tail")

        #total weight variables
        Wvt = Variable('W_{vt}', 'N', 'Total VT System Weight')
        fVT = Variable('f_{VT}', '-', 'VT Fractional Weight')

        #Margin and Sensitivity
        CVT = Variable('C_{VT}', 1, '-', 'VT Weight Margin and Sensitivity')

        #variables only used for the TASOPT tail drag formulation
        cdfv = Variable('c_{d_{fv}}', '-', 'VT friction drag coefficient')
        cdpv = Variable('c_{d_{pv}}', '-', 'VT pressure drag coefficient')
        coslamcube = Variable('\\cos(\\Lambda_{vt})^3', '-', 'Cosine of tail sweep cubed')

        numspar = Variable('N_{spar}', '-', 'Number of Spars in Each VT Carrying Stress in 1 in 20 Case')

        constraints = [
            self.vtns['\\lambda_{vt}'] == self.wb['taper'],
            Wvt >= numspar*CVT*(self.wb['W_{struct}'] + self.wb['W_{struct}'] * fVT),
            ]

        return self.vtns, self.wb, constraints

    def dynamic(self, state, fitDrag):
        """"
        creates a horizontal tail performance model
        """
        return VerticalTailPerformance(self, state, fitDrag)

class VerticalTailNoStruct(Model):
    """
    Vertical tail sizing w/no structural model
    SKIP VERIFICATION

    References:
    1: LEAP engine specs (1B)
    2: TASOPT 737 code
    3: Boeing 737 airport doc
    4: Boeing 737 Max airport doc
    5: http://www.digitaldutch.com/atmoscalc
    6: Engineering toolbox
    7: Boeing.com
    """
    def setup(self, **kwargs):
        #define new variables
        Avt    = Variable('A_{vt}', '-', 'Vertical tail aspect ratio')
        CDwm   = Variable('C_{D_{wm}}', '-', 'Windmill drag coefficient')
        Dwm    = Variable('D_{wm}', 'N', 'Engine out windmill drag')
        Lvmax  = Variable('L_{vt_{max}}', 'N',
                          'Maximum load for structural sizing')
        CLvmax = Variable('C_{L_{vt,max}}', '-', 'Max lift coefficient')
        CLvtEO   = Variable('C_{L_{vt,EO}}', '-', 'Vertical tail lift coefficient (engine out)')
        clvtEO   = Variable('c_{l_{vt,EO}}', '-',
                            'Sectional lift force coefficient (engine out)')
        LvtEO    = Variable('L_{vt,EO}', 'N', 'Vertical tail lift in engine out')
        Svt    = Variable('S_{vt}', 'm^2', 'Vertical tail reference area')
        V1     = Variable('V_1', 'm/s', 'Minimum takeoff velocity')
        bvt    = Variable('b_{vt}', 'm', 'Vertical tail span')
        cma    = Variable('\\bar{c}_{vt}', 'm', 'Vertical tail mean aero chord')
        croot  = Variable('c_{root_{vt}}', 'm', 'Vertical tail root chord')
        ctip   = Variable('c_{tip_{vt}}', 'm', 'Vertical tail tip chord')
        e      = Variable('e_{vt}', '-', 'Span efficiency of vertical tail')
        dxlead  = Variable('\\Delta x_{lead_{vt}}', 'm',
                           'Distance from CG to vertical tail leading edge')
        dxtrail = Variable('\\Delta x_{trail_{vt}}', 'm',
                           'Distance from CG to vertical tail trailing edge')
        lvt    = Variable('l_{vt}', 'm', 'Vertical tail moment arm')
        mu0    = Variable('\\mu_0', 1.8E-5, 'N*s/m^2', 'Dynamic viscosity (SL)')
        p      = Variable('p_{vt}', '-', 'Substituted variable = 1 + 2*taper')
        q      = Variable('q_{vt}', '-', 'Substituted variable = 1 + taper')
        rho0   = Variable('\\rho_{TO}', 'kg/m^3', 'Air density (SL))')
        tanL   = Variable('\\tan(\\Lambda_{vt})', '-',
                          'Tangent of leading edge sweep')
        taper  = Variable('\\lambda_{vt}', '-', 'Vertical tail taper ratio')
        tau    = Variable('\\tau_{vt}', '-', 'Vertical tail thickness/chord ratio')
        xCGvt  = Variable('x_{CG_{vt}}', 'm', 'x-location of tail CG')
        y_eng  = Variable('y_{eng}', 'm', 'Engine moment arm')
        ymac    = Variable('y_{\\bar{c}_{vt}}', 'm',
                           'Spanwise location of mean aerodynamic chord')
        zmac   = Variable('z_{\\bar{c}_{vt}}', 'm',
                          'Vertical location of mean aerodynamic chord')
        Vvt = Variable('V_{vt}', '-', 'Vertical Tail Volume Coefficient')
        Vvtmin = Variable('V_{vt_{min}}', '-', 'Minimum Vertical Tail Volume Coefficient')

        #engine values
        Te     = Variable('T_e', 'N', 'Thrust per engine at takeoff')

        #variables specific to yaw rate sizing
        Vland = Variable('V_{land}', 'm/s', 'Aircraft Landing Speed')
        CLvyaw = Variable('C_{L_{vt,yaw}}', '-', 'VT CL at rotation')
        Iz = Variable('I_{z, max}', 'kg*m^2', 'Aircraft Z-axis Moment of Inertia')
        rreq = Variable('\\dot{r}_{req}', 's^-2', 'Required Yaw Rate at Landing')

        #constraints
        constraints = []

        with SignomialsEnabled():
            #non vectorized constraints
            constraints.extend([
               #constraint tail Cl at flare
                CLvyaw == .85*CLvmax,

                LvtEO == 0.5*rho0*V1**2*Svt*CLvtEO,
                # Vertical tail force (y-direction) for engine out

                TCS([CLvtEO*(1 + clvtEO/(np.pi*e*Avt)) <= clvtEO]),
                #engine out CL computation

                Avt == bvt**2/Svt,

                # Tail geometry relationship
                Svt <= bvt*(croot + ctip)/2, # [SP]

                # Fuselage length constrains the tail trailing edge
                TCS([p >= 1 + 2*taper]),
                TCS([2*q >= 1 + p]),
                # Mean aerodynamic chord calculations
                ymac == (bvt/3)*q/p,
                zmac == (bvt/3)*q/p,

##                TCS([(2./3)*(1 + taper + taper**2)*croot/q >= cma]), # [SP]
                SignomialEquality((2./3)*(1 + taper + taper**2)*croot/q, cma),

                # Define vertical tail geometry
                taper == ctip/croot,
                # Moment arm and geometry -- same as for htail
                dxlead + croot <= dxtrail,
                TCS([dxlead + ymac * tanL + 0.25 * cma >= lvt], reltol=1e-2),  # [SP]


                # TODO: Constrain taper by tip Reynolds number
                taper >= 0.25,

                #Enforce a minimum vertical tail volume
                Vvt >= Vvtmin,
                ])

        return constraints

class VerticalTailPerformance(Model):
    """
    Vertical tail perofrmance model
    SKIP VERIFICATION

    ARGUMENTS
    ---------
    fitDrag: True = use Martin's tail drag fits, False = use the TASOPT tail drag model
    """
    def setup(self, vt, state, fitDrag):
        self.vt = vt

        #define new variables
        Dvt    = Variable('D_{vt}', 'N', 'Vertical tail viscous drag, cruise')
        Rec    = Variable('Re_{vt}', '-', 'Vertical tail reynolds number, cruise')
        CDvis  = Variable('C_{D_{vis}}', '-', 'Viscous drag coefficient')

        #constraints
        constraints = []

         #vectorized constraints
        with SignomialsEnabled():
            constraints.extend([
                # Finite wing theory
                # people.clarkson.edu/~pmarzocc/AE429/AE-429-4.pdf
                # Valid because tail is untwisted and uncambered
                # (lift curve slope passes through origin)

                Dvt >= 0.5*state['\\rho']*state['V']**2*self.vt['S_{vt}']*CDvis,

                # Cruise Reynolds number
                Rec == state.atm['\\rho']*state['V']*self.vt['\\bar{c}_{vt}']/state.atm['\\mu'],

              ])

        if fitDrag:
            constraints.extend([
                #Martin's TASOPT tail fit
                    CDvis**1.18909 >= 2.43701e-77 * (Rec)**-0.52841 * (self.vt['\\tau_{vt}'])**133.796 * (state['M'])**1022.7
                                + 0.00304307 * (Rec)**-0.409988 * (self.vt['\\tau_{vt}'])**1.22062 * (state['M'])**1.55119
                                + 0.000196709 * (Rec)**0.214479 * (self.vt['\\tau_{vt}'])**-0.0383195 * (state['M'])**-0.137561
                                + 6.59349e-50 * (Rec)**-0.498092 * (self.vt['\\tau_{vt}'])**1.55922 * (state['M'])**-114.577

            #Philippe thesis fit
            # Vertical tail viscous drag in cruise
            # Data fit from Xfoil
##            CDvis**0.125 >= 0.19*(self.vt['\\tau_{vt}'])**0.0075 *(Rec)**0.0017
##                        + 1.83e+04*(self.vt['\\tau_{vt}'])**3.54*(Rec)**-0.494
##                        + 0.118*(self.vt['\\tau_{vt}'])**0.0082 *(Rec)**0.00165
##                        + 0.198*(self.vt['\\tau_{vt}'])**0.00774*(Rec)**0.00168,
            ])
        else:
            None
            #drag constraints found in aircraftP

        return constraints

## Old code that doens't work. Could be used for subystem testing.

##if __name__ == '__main__':
##    plot = True
##
##    substitutions = {
##            'R_{req}': 500,
##            'CruiseAlt': 30000,
##            'numeng': 2,
####            'W_{Load_{max}}': 6664,
##            'W_{pass}': 91 * 9.81,
##            'n_{pass}': 150,
##            'pax_{area}': 1,
####            'C_{D_{fuse}}': .005, #assumes flat plate turbulent flow, from wikipedia
##            'e': .9,
##            'b_{max}': 60,
##            'W_{engine}': 1000,
##
##            #VT subs
##           'C_{D_{wm}}': 0.5, # [2]
##           'C_{L_{vt,max}}': 2.6, # [2]
##           'V_1': 70,
##           'V_{ne}': 144, # [2]
##           '\\rho_{TO}': 1.225,
##           '\\tan(\\Lambda_{vt})': np.tan(40*np.pi/180),
####           'c_{l_{vt}}': 0.5, # [2]
##           'c_{l_{vt,EO}}': 0.5,
##           'A_2': np.pi*(.5*1.75)**2, # [1]
##           'e_{vt}': 0.8,
##           'l_{fuse}': 39,
####           'x_{CG}': 18,
##           'y_{eng}': 4.83, # [3]
##
##           'V_{land}': 72,
##           'I_{z}': 12495000, #estimate for late model 737 at max takeoff weight (m l^2/12)
##           '\\dot{r}_{req}': 0.174533, #10 deg/s yaw rate
##
##            'N_{spar}': 2,
##            }
##
##    if plot == True:
##        mission = Mission()
##        m = Model(mission['W_{f_{total}}'], mission)
##        m.substitutions.update(substitutions)
##        sol = m.localsolve(solver='mosek', verbosity = 4)

##import matplotlib.pyplot as plt
##from gpkit.small_scripts import mag
##from simple_ac_imports_no_engine import Wing, Fuselage, Engine, CruiseP, ClimbP, FlightState, CruiseSegment, ClimbSegment
##    class Aircraft(Model):
##    "Aircraft class"
##    def setup(self, **kwargs):
##        #create submodels
##        self.fuse = Fuselage()
##        self.wing = Wing()
##        self.engine = Engine()
##        self.VT = VerticalTail()
##
##        #variable definitions
##        numeng = Variable('numeng', '-', 'Number of Engines')
##
##        constraints = []
##
##        constraints.extend([
##            numeng == numeng, #need numeng in the model
##            ])
##
##        self.components = [self.fuse, self.wing, self.engine, self.VT]
##
##        return self.components, constraints
##
##    def dynamic(self, state):
##        """
##        creates an aircraft climb performance model, given a state
##        """
##        return AircraftP(self, state)
##
##    def climb_dynamic(self, state):
##        """
##        creates an aircraft climb performance model, given a state
##        """
##        return ClimbP(self, state)
##
##    def cruise_dynamic(self, state):
##        """
##        creates an aircraft cruise performance model, given a state
##        """
##        return CruiseP(self, state)
##
##class AircraftP(Model):
##    """
##    aircraft performance models superclass, contains constraints true for
##    all flight segments
##    """
##    def  setup(self, aircraft, state, **kwargs):
##        #make submodels
##        self.aircraft = aircraft
##        self.wingP = aircraft.wing.dynamic(state)
##        self.fuseP = aircraft.fuse.dynamic(state)
##        self.engineP = aircraft.engine.dynamic(state)
##        self.VTP = aircraft.VT.dynamic(state)
##
##        self.Pmodels = [self.wingP, self.fuseP, self.engineP, self.VTP]
##
##        #variable definitions
##        Vstall = Variable('V_{stall}', 'knots', 'Aircraft Stall Speed')
##        D = Variable('D', 'N', 'Total Aircraft Drag')
##        W_avg = Variable('W_{avg}', 'N', 'Geometric Average of Segment Start and End Weight')
##        W_start = Variable('W_{start}', 'N', 'Segment Start Weight')
##        W_end = Variable('W_{end}', 'N', 'Segment End Weight')
##        W_burn = Variable('W_{burn}', 'N', 'Segment Fuel Burn Weight')
##        WLoadmax = Variable('W_{Load_{max}}', 'N/m^2', 'Max Wing Loading')
##        WLoad = Variable('W_{Load}', 'N/m^2', 'Wing Loading')
##        t = Variable('tmin', 'min', 'Segment Flight Time in Minutes')
##        thours = Variable('thr', 'hour', 'Segment Flight Time in Hours')
##
##        xCG     = Variable('x_{CG}', 'm', 'CG location')
##
##        constraints = []
##
##        constraints.extend([
##            #speed must be greater than stall speed
##            state['V'] >= Vstall,
##
##
##            #Figure out how to delete
##            Vstall == 120*units('kts'),
##            WLoadmax == 6664 * units('N/m^2'),
##
##            #compute the drag
##            TCS([D >= self.wingP['D_{wing}'] + self.fuseP['D_{fuse}'] + self.VTP['D_{vt}']]),
##
##            #constraint CL and compute the wing loading
##            W_avg == .5*self.wingP['C_{L}']*self.aircraft['S']*state.atm['\\rho']*state['V']**2,
##            WLoad == .5*self.wingP['C_{L}']*self.aircraft['S']*state.atm['\\rho']*state['V']**2/self.aircraft.wing['S'],
##
##            #set average weight equal to the geometric avg of start and end weight
##            W_avg == (W_start * W_end)**.5,
##
##            #constrain the max wing loading
##            WLoad <= WLoadmax,
##
##            #compute fuel burn from TSFC
##            W_burn == aircraft['numeng']*self.engineP['TSFC'] * thours * self.engineP['F'],
##
##            #time unit conversion
##            t == thours,
##
##            #make lift equal weight --> small angle approx in climb
##            self.wingP['L_{wing}'] == W_avg,
##
##            xCG == 18*units('m'),
##            ])
##
##        return self.Pmodels, constraints
##
##class Mission(Model):
##    """
##    mission class, links together all subclasses
##    """
##    def setup(self):
##        #define the number of each flight segment
##        Nclimb = 2
##        Ncruise = 2
##
##        #build required submodels
##        ac = Aircraft()
##
##        #vectorize
##        with Vectorize(Nclimb):
##            climb = ClimbSegment(ac)
##
##        with Vectorize(Ncruise):
##            cruise = CruiseSegment(ac)
##
##        #declare new variables
##        W_ftotal = Variable('W_{f_{total}}', 'N', 'Total Fuel Weight')
##        W_fclimb = Variable('W_{f_{climb}}', 'N', 'Fuel Weight Burned in Climb')
##        W_fcruise = Variable('W_{f_{cruise}}', 'N', 'Fuel Weight Burned in Cruise')
##        W_total = Variable('W_{total}', 'N', 'Total Aircraft Weight')
##        CruiseAlt = Variable('CruiseAlt', 'ft', 'Cruise Altitude [feet]')
##        ReqRng = Variable('R_{req}', 'nautical_miles', 'Required Cruise Range')
##
##        h = climb.state['h']
##        hftClimb = climb.state['hft']
##        dhft = climb.climbP['dhft']
##        hftCruise = cruise.state['hft']
##
##        #make overall constraints
##        constraints = []
##
##        constraints.extend([
##            #weight constraints
##            TCS([ac['W_{e}'] + ac['W_{payload}'] + W_ftotal + ac['numeng'] * ac['W_{engine}'] + ac['W_{wing}'] + 2*ac.VT['W_{struct}'] <= W_total]),
##
##            climb.climbP.aircraftP['W_{start}'][0] == W_total,
##            climb.climbP.aircraftP['W_{end}'][-1] == cruise.cruiseP.aircraftP['W_{start}'][0],
##
##            # similar constraint 1
##            TCS([climb.climbP.aircraftP['W_{start}'] >= climb.climbP.aircraftP['W_{end}'] + climb.climbP.aircraftP['W_{burn}']]),
##            # similar constraint 2
##            TCS([cruise.cruiseP.aircraftP['W_{start}'] >= cruise.cruiseP.aircraftP['W_{end}'] + cruise.cruiseP.aircraftP['W_{burn}']]),
##
##            climb.climbP.aircraftP['W_{start}'][1:] == climb.climbP.aircraftP['W_{end}'][:-1],
##            cruise.cruiseP.aircraftP['W_{start}'][1:] == cruise.cruiseP.aircraftP['W_{end}'][:-1],
##
##            TCS([ac['W_{e}'] + ac['W_{payload}'] + ac['numeng'] * ac['W_{engine}'] + ac['W_{wing}'] + 2*ac.VT['W_{struct}'] <= cruise.cruiseP.aircraftP['W_{end}'][-1]]),
##
##            TCS([W_ftotal >=  W_fclimb + W_fcruise]),
##            TCS([W_fclimb >= sum(climb.climbP['W_{burn}'])]),
##            TCS([W_fcruise >= sum(cruise.cruiseP['W_{burn}'])]),
##
##            #altitude constraints
##            hftCruise == CruiseAlt,
##            TCS([hftClimb[1:Ncruise] >= hftClimb[:Ncruise-1] + dhft]),
##            TCS([hftClimb[0] >= dhft[0]]),
##            hftClimb[-1] <= hftCruise,
##
##            #compute the dh
##            dhft == hftCruise/Nclimb,
##
##            #constrain the thrust
##            climb.climbP.engineP['F'] <= 2 * max(cruise.cruiseP.engineP['F']),
##
##            #set the range for each cruise segment, doesn't take credit for climb
##            #down range disatnce covered
##            cruise.cruiseP['Rng'] == ReqRng/(Ncruise),
##
##            #set the TSFC
##            climb.climbP.engineP['TSFC'] == .7*units('1/hr'),
##            cruise.cruiseP.engineP['TSFC'] == .5*units('1/hr'),
##            ])
##
##        #VT constriants
##        constraints.extend([
##            ac.VT['T_e'] == climb.climbP.engineP['F'][0],
##
##            # Drag of a windmilling engine
##            ac.VT['D_{wm}'] >= 0.5*ac.VT['\\rho_{TO}']*ac.VT['V_1']**2*ac.engine['A_2']*ac.VT['C_{D_{wm}}'],
##
##            ac.VT['x_{CG_{vt}}'] <= ac.fuse['l_{fuse}'],
##
##            #VTP constraints
##            ac.fuse['l_{fuse}'] >= ac.VT['\\Delta x_{lead_v}'] + climb.climbP.aircraftP['x_{CG}'],
##            ac.VT['x_{CG_{vt}}'] >= climb.climbP.aircraftP['x_{CG}']+(ac.VT['\\Delta x_{lead_v}']+ac.VT['\\Delta x_{trail_v}'])/2,
##
##            ac.fuse['l_{fuse}'] >= ac.VT['\\Delta x_{lead_v}'] + cruise.cruiseP.aircraftP['x_{CG}'],
##            ac.VT['x_{CG_{vt}}'] >= cruise.cruiseP.aircraftP['x_{CG}']+(ac.VT['\\Delta x_{lead_v}']+ac.VT['\\Delta x_{trail_v}'])/2,
##            ])
##
##        # Model.__init__(self, W_ftotal + s*units('N'), constraints + ac + climb + cruise, subs)
##        self.cost = W_ftotal
##        return constraints, ac, cruise, climb
