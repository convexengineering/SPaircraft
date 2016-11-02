"""Simple commercial aircraft flight profile and aircraft model"""
from numpy import pi
import numpy as np
from gpkit import Variable, Model, units, SignomialsEnabled, SignomialEquality, vectorize
from gpkit.tools import te_exp_minus1
from gpkit.constraints.tight import TightConstraintSet as TCS

#only needed for the local bounded debugging tool
from collections import defaultdict
from gpkit.small_scripts import mag

"""
Models requird to minimze the aircraft total fuel weight. Rate of climb equation taken from John
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

g = 9.81*units('m*s**-2')

class Aircraft(Model):
    "Aircraft class"
    def __init__(self, **kwargs):
        #create submodels
        self.fuse = Fuselage()
        self.wing = Wing()
        self.engine = Engine()

        #variable definitions
        numeng = Variable('numeng', '-', 'Number of Engines')

        constraints = []

        constraints.extend([
            numeng == numeng, #need numeng in the model
            ])

        self.components = [self.fuse, self.wing, self.engine]

        Model.__init__(self, None, [self.components + constraints], **kwargs)
        
    def climb_dynamic(self, state):
        """
        creates an aircraft climb performance model, given a state
        """
        return ClimbP(self, state)

    def cruise_dynamic(self, state):
        """
        creates an aircraft cruise performance model, given a state
        """
        return CruiseP(self, state)

class AircraftP(Model):
    """
    aircraft performance models superclass, contains constraints true for
    all flight segments
    """
    def  __init__(self, aircraft, state, **kwargs):
        #make submodels
        self.aircraft = aircraft
        self.wingP = aircraft.wing.dynamic(state)
        self.fuseP = aircraft.fuse.dynamic(state)
        self.engineP = aircraft.engine.dynamic(state)

        self.Pmodels = [self.wingP, self.fuseP, self.engineP]

        #variable definitions
        Vstall = Variable('V_{stall}', 'knots', 'Aircraft Stall Speed')
        D = Variable('D', 'N', 'Total Aircraft Drag')
        W_avg = Variable('W_{avg}', 'N', 'Geometric Average of Segment Start and End Weight')
        W_start = Variable('W_{start}', 'N', 'Segment Start Weight')
        W_end = Variable('W_{end}', 'N', 'Segment End Weight')
        W_burn = Variable('W_{burn}', 'N', 'Segment Fuel Burn Weight')
        WLoadmax = Variable('W_{Load_max}', 'N/m^2', 'Max Wing Loading')
        WLoad = Variable('W_{Load}', 'N/m^2', 'Wing Loading')
        t = Variable('tmin', 'min', 'Segment Flight Time in Minutes')
        thours = Variable('thr', 'hour', 'Segment Flight Time in Hours')

        constraints = []

        constraints.extend([
            #speed must be greater than stall speed
            state['V'] >= Vstall,


            #Figure out how to delete
            Vstall == 120*units('kts'),
            WLoadmax == 6664 * units('N/m^2'),

            #compute the drag
            TCS([D >= self.wingP['D_{wing}'] + self.fuseP['D_{fuse}']]),

            #constraint CL and compute the wing loading
            W_avg == .5*self.wingP['C_{L}']*self.aircraft['S']*state.atm['\\rho']*state['V']**2,      
            WLoad == .5*self.wingP['C_{L}']*self.aircraft['S']*state.atm['\\rho']*state['V']**2/self.aircraft.wing['S'],

            #set average weight equal to the geometric avg of start and end weight
            W_avg == (W_start * W_end)**.5,

            #constrain the max wing loading
            WLoad <= WLoadmax,

            #compute fuel burn from TSFC
            W_burn == aircraft['numeng']*self.engineP['TSFC'] * thours * self.engineP['thrust'],
               
            #time unit conversion
            t == thours,
            ])

        Model.__init__(self, None, [self.Pmodels + constraints], **kwargs)

class ClimbP(Model):
    """
    Climb constraints
    """
    def __init__(self, aircraft, state, **kwargs):
        #submodels
        self.aircraft = aircraft
        self.aircraftP = AircraftP(aircraft, state)
        self.wingP = self.aircraftP.wingP
        self.fuseP = self.aircraftP.fuseP
        self.engineP = self.aircraftP.engineP
                                  
        #variable definitions
        theta = Variable('\\theta', '-', 'Aircraft Climb Angle')
        excessP = Variable('excessP', 'W', 'Excess Power During Climb')
        RC = Variable('RC', 'feet/min', 'Rate of Climb/Decent')
        dhft = Variable('dhft', 'feet', 'Change in Altitude Per Climb Segment [feet]')
        RngClimb = Variable('RngClimb', 'nautical_miles', 'Down Range Covered in Each Climb Segment')

        #constraints
        constraints = []
        
        constraints.extend([
           #constraint on drag and thrust
            self.aircraft['numeng']*self.engineP['thrust'] >= self.aircraftP['D'] + self.aircraftP['W_{avg}'] * theta,
            
            #climb rate constraints
            TCS([excessP + state['V'] * self.aircraftP['D'] <=  state['V'] * aircraft['numeng'] * self.engineP['thrust']]),
            
            RC == excessP/self.aircraftP['W_{avg}'],
            RC >= 500*units('ft/min'),
            
            #make the small angle approximation and compute theta
            theta * state['V']  == RC,
           
            dhft == self.aircraftP['tmin'] * RC,
        
            #makes a small angle assumption during climb
            RngClimb == self.aircraftP['thr']*state['V'],
            ])

        Model.__init__(self, None, constraints + self.aircraftP)

class CruiseP(Model):
    """
    Cruise constraints
    """
    def __init__(self, aircraft, state, **kwargs):
        self.aircraft = aircraft
        self.aircraftP = AircraftP(aircraft, state)
        self.wingP = self.aircraftP.wingP
        self.fuseP = self.aircraftP.fuseP
        self.engineP = self.aircraftP.engineP
                        
        #variable definitions
        z_bre = Variable('z_{bre}', '-', 'Breguet Parameter')
        Rng = Variable('Rng', 'nautical_miles', 'Cruise Segment Range')

        constraints = []

        constraints.extend([
             #steady level flight constraint on D 
             self.aircraftP['D'] == aircraft['numeng'] * self.engineP['thrust'],

             #taylor series expansion to get the weight term
             TCS([self.aircraftP['W_{burn}']/self.aircraftP['W_{end}'] >=
                  te_exp_minus1(z_bre, nterm=3)]),

             #breguet range eqn
             # old version -- possibly unneeded numeng
 #            TCS([z_bre >= (self.aircraft['numeng'] * self.engineP['TSFC'] * self.aircraftP['thr']*
 #                           self.aircraftP['D']) / self.aircraftP['W_{avg}']]),

            # new version -- needs to be thought through carefully
             # seems correct to me - I switched T to D below (steady level flight) but fogot
             #about the Negn term
             TCS([z_bre >= (self.engineP['TSFC'] * self.aircraftP['thr']*
                            self.aircraftP['D']) / self.aircraftP['W_{avg}']]),

             #time
             self.aircraftP['thr'] * state['V'] == Rng,
             ])

        Model.__init__(self, None, constraints + self.aircraftP)

class CruiseSegment(Model):
    """
    Combines a flight state and aircrat to form a cruise flight segment
    """
    def __init__(self, aircraft, **kwargs):
        self.state = FlightState()
        self.cruiseP = aircraft.cruise_dynamic(self.state)

        Model.__init__(self, None, [self.state, self.cruiseP], **kwargs)

class ClimbSegment(Model):
    """
    Combines a flight state and aircrat to form a cruise flight segment
    """
    def __init__(self, aircraft, **kwargs):
        self.state = FlightState()
        self.climbP = aircraft.climb_dynamic(self.state)

        Model.__init__(self, None, [self.state, self.climbP], **kwargs)

class FlightState(Model):
    """
    creates atm model for each flight segment, has variables
    such as veloicty and altitude
    """
    def __init__(self,**kwargs):
        #make an atmosphere model
        self.atm = Atmosphere()
        
        #declare variables
        V = Variable('V', 'kts', 'Aircraft Flight Speed')
        a = Variable('a', 'm/s', 'Speed of Sound')
        h = Variable('h', 'm', 'Segment Altitude [meters]')
        hft = Variable('hft', 'feet', 'Segment Altitude [feet]')
        R = Variable('R', 287, 'J/kg/K', 'Air Specific Heat')
        gamma = Variable('\\gamma', 1.4, '-', 'Air Specific Heat Ratio')
        M = Variable('M', '-', 'Mach Number')

        #make new constraints
        constraints = []

        constraints.extend([
            V == V, #required so velocity variable enters the model

            h == hft, #convert the units on altitude

            #compute the speed of sound with the state
            a  == (gamma * R * self.atm['T_{atm}'])**.5,

            #compute the mach number
            V == M * a,
            ])

        #build the model
        Model.__init__(self, None, constraints + self.atm, **kwargs)

class Atmosphere(Model):
    def __init__(self, **kwargs):
        g = Variable('g', 'm/s^2', 'Gravitational acceleration')
        p_sl = Variable("p_{sl}", "Pa", "Pressure at sea level")
        T_sl = Variable("T_{sl}", "K", "Temperature at sea level")
        L_atm = Variable("L_{atm}", "K/m", "Temperature lapse rate")
        M_atm = Variable("M_{atm}", "kg/mol",
                         "Molar mass of dry air")
        p_atm = Variable("P_{atm}", "Pa", "air pressure")
        R_atm = Variable("R_{atm}", "J/mol/K", "air specific heating value")
        TH = 5.257386998354459 #(g*M_atm/R_atm/L_atm).value
        rho = Variable('\\rho', 'kg/m^3', 'Density of air')
        T_atm = Variable("T_{atm}", "K", "air temperature")
        h = Variable("h", "m", "Altitude")

        """
        Dynamic viscosity (mu) as a function of temperature
        References:
        http://www-mdp.eng.cam.ac.uk/web/library/enginfo/aerothermal_dvd_only/aero/
            atmos/atmos.html
        http://www.cfd-online.com/Wiki/Sutherland's_law
        """
        mu  = Variable('\\mu', 'kg/(m*s)', 'Dynamic viscosity')

        T_s = Variable('T_s', 110.4, "K", "Sutherland Temperature")
        C_1 = Variable('C_1', 1.458E-6, "kg/(m*s*K^0.5)",
                       'Sutherland coefficient')
        
##        t_plus_ts_approx = (T_atm + T_s).mono_approximation({T_atm: 288.15,
##                                                         T_s: T_s.value})

        with SignomialsEnabled():
            constraints = [
                # Pressure-altitude relation
                (p_atm/p_sl)**(1/5.257) == T_atm/T_sl,

                # Ideal gas law
                rho == p_atm/(R_atm/M_atm*T_atm),

                #temperature equation
##                SignomialEquality(T_sl, T_atm + L_atm*h),
                T_atm == 218*units('K'),

                #constraint on mu
##                SignomialEquality((T_atm + T_s) * mu, C_1 * T_atm**1.5),
##                TCS([(T_atm + T_s) * mu >= C_1 * T_atm**1.5])
                T_sl == 288.15*units('K'),
                p_sl == 101325*units('Pa'),
                R_atm == 8.31447*units('J/mol/K'),
                M_atm == .0289644*units('kg/mol')
                ]

        #like to use a local subs here in the future
        subs = None

        Model.__init__(self, None, constraints, subs)

class Engine(Model):
    """
    place holder engine model
    """
    def __init__(self, **kwargs):
        #new variables
        W_engine = Variable('W_{engine}', 'N', 'Weight of a Single Turbofan Engine')
        
        constraints = []

        constraints.extend([
            W_engine == 1000 * units('N')
            ])

        Model.__init__(self, None, constraints)

    def dynamic(self, state):
            """
            returns an engine performance model
            """
            return EnginePerformance(self, state)

class EnginePerformance(Model):
    """
    place holder engine perofrmacne model
    """
    def __init__(self, engine, state, **kwargs):
        #new variables
        TSFC = Variable('TSFC', '1/hr', 'Thrust Specific Fuel Consumption')
        thrust = Variable('thrust', 'N', 'Thrust')
        
        #constraints
        constraints = []

        constraints.extend([
            TSFC == TSFC,

            thrust == thrust, #want thrust to enter the model
            ])

        Model.__init__(self, None, constraints)


class Wing(Model):
    """
    place holder wing model
    """
    def __init__(self, ** kwargs):
        #new variables
        W_wing = Variable('W_{wing}', 'N', 'Wing Weight')
                           
        #aircraft geometry
        S = Variable('S', 'm^2', 'Wing Planform Area')
        AR = Variable('AR', '-', 'Aspect Ratio')
        span = Variable('b', 'm', 'Wing Span')
        span_max = Variable('b_{max}', 'm', 'Max Wing Span')

        K = Variable('K', '-', 'K for Parametric Drag Model')
        e = Variable('e', '-', 'Oswald Span Efficiency Factor')
        
        constraints = []

        constraints.extend([
            #wing weight constraint
            #based off of a raymer weight and 737 data from TASOPT output file
            (S/(124.58*units('m^2')))**.65 == W_wing/(105384.1524*units('N')),

            #compute wing span and aspect ratio, subject to a span constraint
            AR == (span**2)/S,
            #AR == 9,

            span <= span_max,

            #compute K for the aircraft
            K == (pi * e * AR)**-1,
            ])

        Model.__init__(self, None, constraints)

    def dynamic(self, state):
        """
        creates an instance of the wing's performance model
        """
        return WingPerformance(self, state)
        

class WingPerformance(Model):
    """
    wing aero modeling
    """
    def __init__(self, wing, state, **kwargs):
        #new variables
        CL= Variable('C_{L}', '-', 'Lift Coefficient')
        Cdw = Variable('C_{d_w}', '-', 'Cd for a NC130 Airfoil at Re=2e7')
        Dwing = Variable('D_{wing}', 'N', 'Total Wing Drag')

        #constraints
        constraints = []

        constraints.extend([
            #airfoil drag constraint
            TCS([Cdw**6.5 >= (1.02458748e10 * CL**15.587947404823325 * state['M']**156.86410659495155 +
                         2.85612227e-13 * CL**1.2774976672501526 * state['M']**6.2534328002723703 +
                         2.08095341e-14 * CL**0.8825277088649582 * state['M']**0.0273667615730107 +
                         1.94411925e+06 * CL**5.6547413360261691 * state['M']**146.51920742858428)]),
            TCS([Dwing >= (.5*wing['S']*state.atm['\\rho']*state['V']**2)*(Cdw + wing['K']*CL**2)]),
            ])

        Model.__init__(self, None, constraints)

class Fuselage(Model):
    """
    place holder fuselage model
    """
    def __init__(self, **kwargs):
        #new variables
        dPover       = Variable('\\delta_P_{over}','psi','Cabin overpressure')
        npass        = Variable('n_{pass}', '-', 'Number of Passengers to Carry')
        Nland        = Variable('N_{land}',6.0,'-', 'Emergency landing load factor') #[TAS]
        VNE          = Variable('V_{NE}','m/s','Never-exceed speed') #[Philippe]
        SPR          = Variable('SPR', '-', 'Number of seats per row')
        nrows        = Variable('n_{rows}', '-', 'Number of rows')
        nseat        = Variable('n_{seat}','-','Number of seats')
        pitch        = Variable('p_s', 'cm', 'Seat pitch')


        # Cross-sectional variables
        Adb          = Variable('A_{db}', 'm^2', 'Web cross sectional area')
        Afloor       = Variable('A_{floor}', 'm^2', 'Floor beam x-sectional area')
        Afuse        = Variable('A_{fuse}', 'm^2', 'Fuselage x-sectional area')
        Askin        = Variable('A_{skin}', 'm^2', 'Skin cross sectional area')
        hdb          = Variable('h_{db}','m', 'Web half-height')
        hfloor       = Variable('h_{floor}', 'm', 'Floor beam height')        
        Rfuse        = Variable('R_{fuse}', 'm', 'Fuselage radius') # will assume for now there: no under-fuselage extension deltaR
        tdb          = Variable('t_{db}', 'm', 'Web thickness')
        thetadb      = Variable('\\theta_{db}','-','DB fuselage joining angle')
        tshell       = Variable('t_{shell}', 'm', 'Shell thickness')
        tskin        = Variable('t_{skin}', 'm', 'Skin thickness')
        waisle       = Variable('w_{aisle}', 'm', 'Aisle width') #[Boeing]
        wdb          = Variable('w_{db}','m','DB added half-width')
        wfloor       = Variable('w_{floor}', 'm', 'Floor half-width')
        wfuse        = Variable('w_{fuse}', 'm', 'Fuselage width')
        wseat        = Variable('w_{seat}','m', 'Seat width') #[Philippe]
        wsys         = Variable('w_{sys}','m', 'Width between cabin and skin for systems') #[Philippe]

        #Tail cone variables
        lamcone      = Variable('\\lambda_{cone}', '-','Tailcone radius taper ratio (xshell2->xtail)')
        lcone        = Variable('l_{cone}', 'm', 'Cone length')
        
        # Lengths (free)
        lfuse        = Variable('l_{fuse}', 'm', 'Fuselage length')
        lnose        = Variable('l_{nose}', 'm', 'Nose length')
        lshell       = Variable('l_{shell}', 'm', 'Shell length')
        lfloor       = Variable('l_{floor}', 'm', 'Floor length')       

        # Surface areas (free)
        Sbulk        = Variable('S_{bulk}', 'm^2', 'Bulkhead surface area')
        Snose        = Variable('S_{nose}', 'm^2', 'Nose surface area')
        
        # Volumes (free)
        Vbulk        = Variable('V_{bulk}', 'm^3', 'Bulkhead skin volume')
        Vcabin       = Variable('V_{cabin}', 'm^3', 'Cabin volume')
        Vcyl         = Variable('V_{cyl}', 'm^3', 'Cylinder skin volume')   
        Vdb          = Variable('V_{db}', 'm^3', 'Web volume')
        Vfloor       = Variable('V_{floor}', 'm^3', 'Floor volume')
        Vnose        = Variable('V_{nose}', 'm^3', 'Nose skin volume')


        # Loads 
        sigfloor     = Variable('\\sigma_{floor}',30000/0.000145, 'N/m^2', 'Max allowable floor stress') #[TAS]
        sigskin      = Variable('\\sigma_{skin}', 15000/0.000145,'N/m^2', 'Max allowable skin stress') #[TAS] 
        sigth        = Variable('\\sigma_{\\theta}', 'N/m^2', 'Skin hoop stress')
        sigx         = Variable('\\sigma_x', 'N/m^2', 'Axial stress in skin')

        Mfloor       = Variable('M_{floor}', 'N*m', 'Max bending moment in floor beams')
        Pfloor       = Variable('P_{floor}','N', 'Distributed floor load')
        Sfloor       = Variable('S_{floor}', 'N', 'Maximum shear in floor beams')
        sigfloor     = Variable('\\sigma_{floor}',30000/0.000145, 'N/m^2', 'Max allowable floor stress') #[TAS]
        taufloor     = Variable('\\tau_{floor}',30000/0.000145, 'N/m^2', 'Max allowable shear web stress') #[TAS]

        # Material properties
        rE           = Variable('r_E', 1,'-', 'Ratio of stringer/skin moduli') #[TAS]
        rhocone      = Variable('\\rho_{cone}','kg/m^3','Cone material density') #[TAS]
        rhobend      = Variable('\\rho_{bend}','kg/m^3', 'Stringer density') #[TAS]
        rhofloor     = Variable('\\rho_{floor}','kg/m^3', 'Floor material density') #[TAS]
        rhoskin      = Variable('\\rho_{skin}','kg/m^3', 'Skin density') #[TAS]
        Wppfloor     = Variable('W\'\'_{floor}','N/m^2', 'Floor weight/area density') #[TAS]
        Wppinsul     = Variable('W\'\'_{insul}','N/m^2', 'Weight/area density of insulation material') #[TAS]
        Wpseat       = Variable('W\'_{seat}','N', 'Weight per seat') #[TAS]
        Wpwindow     = Variable('W\'_{window}','N/m', 'Weight/length density of windows') #[TAS]
        
        # Weight fractions   
        ffadd        = Variable('f_{fadd}','-','Fractional added weight of local reinforcements') #[TAS]
        fframe       = Variable('f_{frame}','-', 'Fractional frame weight') #[Philippe]        
        flugg1       = Variable('f_{lugg,1}','-','Proportion of passengers with one suitcase') #[Philippe]
        flugg2       = Variable('f_{lugg,2}','-','Proportion of passengers with two suitcases') #[Philippe]
        fstring      = Variable('f_{string}','-','Fractional stringer weight') #[Philippe]
                             
        # Weights
        Wavgpass     = Variable('W_{avg. pass}', 'lbf', 'Average passenger weight') #[Philippe]
        Wcargo       = Variable('W_{cargo}', 'N', 'Cargo weight') #[Philippe]        
        Wcarryon     = Variable('W_{carry on}', 'lbf', 'Ave. carry-on weight') #[Philippe]
        Wchecked     = Variable('W_{checked}', 'lbf', 'Ave. checked bag weight') #[Philippe]
        Wdb          = Variable('W_{db}' , 'N', 'Web weight')
        Wfloor       = Variable('W_{floor}', 'lbf', 'Floor weight')
        Wlugg        = Variable('W_{lugg}', 'N', 'Passenger luggage weight')
        Wpass        = Variable('W_{pass}', 'N', 'Passenger weight')
        Wpay         = Variable('W_{pay}', 'N', 'Payload weight')
        Wseat        = Variable('W_{seat}', 'lbf', 'Seating weight')
        Wshell       = Variable('W_{shell}','N','Shell weight')
        Wskin        = Variable('W_{skin}', 'N', 'Skin weight')

        #weight variables
        We = Variable('W_{e}', 'N', 'Empty Weight of Aircraft')
        Wpass = Variable('W_{pass}', 'N', 'Estimated Average Passenger Weight, Includes Baggage')
        pax_area = Variable('pax_{area}', 'm^2', 'Estimated Fuselage Area per Passenger')

        # x-location variables
        xshell1      = Variable('x_{shell1}', 'm', 'Start of cylinder section')
        xshell2      = Variable('x_{shell2}', 'm', 'End of cylinder section')


        constraints = []
        with SignomialsEnabled():
            constraints.extend([
                Nland == Nland,
                VNE == VNE,

                # Passenger constraints
                Wlugg    >= flugg2*npass*2*Wchecked + flugg1*npass*Wchecked + Wcarryon,
                Wpass    == npass*Wavgpass,
                Wpay     >= Wpass + Wlugg + Wcargo,
                Wseat    == Wpseat*nseat,
                nseat    == npass,
                nrows    == nseat/SPR,
                lshell   == nrows*pitch,

                # Fuselage joint angle relations
                thetadb     == wdb/Rfuse, # first order Taylor works...
                thetadb     >= 0.05, thetadb <= 0.5, #Temporarily
                hdb         >= Rfuse*(1.0-.5*thetadb**2), #[SP]

                #compute fuselage area for drag approximation
                #Afuse >= pax_area * npass,
                

                # Cross-sectional constraints
                Adb         == (2*hdb)*tdb,
                Afuse       >= (pi + 2*thetadb + 2*thetadb*(1-thetadb**2/2))*Rfuse**2, #[SP]
                #Afuse       >= (pi + 4*thetadb)*Rfuse**2, #Bad approx, should improve
                Askin       >= (2*pi + 4*thetadb)*Rfuse*tskin + Adb, #no delta R for now
                wfloor      == .5*wfuse,
                wfuse       >= SPR*wseat + 2*waisle + 2*wsys + tdb,
                wfuse       <= 2*(Rfuse + wdb),
                SignomialEquality(tshell,tskin*(1+rE*fstring*rhoskin/rhobend)),

                            # Fuselage surface area relations
                Snose    >= (2*pi + 4*thetadb)*Rfuse**2 *(1/3 + 2/3*(lnose/Rfuse)**(8/5))**(5/8),
                Sbulk    >= (2*pi + 4*thetadb)*Rfuse**2,
            
                # Fuselage length relations
                lfuse    >= lnose+lshell+lcone, 
                lnose    == 0.3*lshell, # Temporarily
                lcone    == Rfuse/lamcone,  

                 ## Stress relations
                #Pressure shell loading
                tskin    == dPover*Rfuse/sigskin,
                tdb      == 2*dPover*wdb/sigskin,
                sigx     == dPover*Rfuse/(2*tshell),
                sigth    == dPover*Rfuse/tskin,
            
                # Volume relations
                Vcyl   == Askin*lshell,
                Vnose  == Snose*tskin,
                Vbulk  == Sbulk*tskin,
                Vdb    == Adb*lshell,
                #TODO Revert to posynomial after debugging
                SignomialEquality(Vcabin, Afuse*(lshell + 0.67*lnose + 0.67*Rfuse)), #[SP] #[SPEquality]

                # Floor loading

                lfloor   >= lshell + 2*Rfuse,            
                Pfloor   >= Nland*(Wpay + Wseat),
                Mfloor   == 9./256.*Pfloor*wfloor,
                Afloor   >= 2.*Mfloor/(sigfloor*hfloor) + 1.5*Sfloor/taufloor,
                Vfloor   == 2*wfloor*Afloor,
                Wfloor   >= rhofloor*g*Vfloor + 2*wfloor*lfloor*Wppfloor,
                Sfloor   == (5./16.)*Pfloor,
                hfloor   <= 0.1*Rfuse,

                # Weight estimate (extra items added for convergence)

                Wdb      == rhoskin*g*Vdb,
                Wskin    >= rhoskin*g*(Vcyl + Vnose + Vbulk),
                Wshell   >= Wskin*(1 + fstring + ffadd + fframe) + Wdb, #+ Whbend, #+ Wvbend,
                We       >= .75*Wpay + Wshell + Wfloor + (lfuse + wfuse+lshell+Rfuse+hdb+(tskin+tshell+tdb)*10)*10*units('N/m'),
                ])

        Model.__init__(self, None, constraints)

    def dynamic(self, state):
        """
        returns a fuselage performance model
        """
        return FuselagePerformance(self, state)

class FuselagePerformance(Model):
    """
    Fuselage performance model
    """
    def __init__(self, fuse, state, **kwargs):
        #new variables
        Cdfuse = Variable('C_{D_{fuse}}', '-', 'Fuselage Drag Coefficient')
        Dfuse = Variable('D_{fuse}', 'N', 'Total Fuselage Drag')
        
        #constraints
        constraints = []

        constraints.extend([
            Dfuse == Cdfuse * (.5 * fuse['A_{fuse}'] * state.atm['\\rho'] * state['V']**2),
            Cdfuse == .005,
            ])

        Model.__init__(self, None, constraints)
    

class Mission(Model):
    """
    mission class, links together all subclasses
    """
    def __init__(self, subs = None, **kwargs):
        #define the number of each flight segment
        Nclimb = 2
        Ncruise = 2

        #build required submodels
        ac = Aircraft()

        #vectorize
        with vectorize(Nclimb):
            cls = ClimbSegment(ac)

        with vectorize(Ncruise):
            crs = CruiseSegment(ac)

        #declare new variables
        W_ftotal = Variable('W_{f_{total}}', 'N', 'Total Fuel Weight')
        W_fclimb = Variable('W_{f_{climb}}', 'N', 'Fuel Weight Burned in Climb')
        W_fcruise = Variable('W_{f_{cruise}}', 'N', 'Fuel Weight Burned in Cruise')
        W_total = Variable('W_{total}', 'N', 'Total Aircraft Weight')
        CruiseAlt = Variable('CruiseAlt', 'ft', 'Cruise Altitude [feet]')
        ReqRng = Variable('ReqRng', 'nautical_miles', 'Required Cruise Range')

        h = cls.state['h']
        hftClimb = cls.state['hft']
        dhft = cls.climbP['dhft']
        hftCruise = crs.state['hft']

        #make overall constraints
        constraints = []

        constraints.extend([
            #weight constraints
            TCS([ac['W_{e}'] + ac['W_{pay}'] + W_ftotal + ac['numeng'] * ac['W_{engine}'] + ac['W_{wing}'] <= W_total]),

            cls.climbP.aircraftP['W_{start}'][0] == W_total,
            cls.climbP.aircraftP['W_{end}'][-1] == crs.cruiseP.aircraftP['W_{start}'][0],

            # similar constraint 1
            TCS([cls.climbP.aircraftP['W_{start}'] >= cls.climbP.aircraftP['W_{end}'] + cls.climbP.aircraftP['W_{burn}']]),
            # similar constraint 2
            TCS([crs.cruiseP.aircraftP['W_{start}'] >= crs.cruiseP.aircraftP['W_{end}'] + crs.cruiseP.aircraftP['W_{burn}']]),

            cls.climbP.aircraftP['W_{start}'][1:] == cls.climbP.aircraftP['W_{end}'][:-1],
            crs.cruiseP.aircraftP['W_{start}'][1:] == crs.cruiseP.aircraftP['W_{end}'][:-1],

            TCS([ac['W_{e}'] + ac['W_{pay}'] + ac['numeng'] * ac['W_{engine}'] + ac['W_{wing}'] <= crs.cruiseP.aircraftP['W_{end}'][-1]]),

            TCS([W_ftotal >=  W_fclimb + W_fcruise]),
            TCS([W_fclimb >= sum(cls.climbP['W_{burn}'])]),
            TCS([W_fcruise >= sum(crs.cruiseP['W_{burn}'])]),

            #altitude constraints
            hftCruise == CruiseAlt,
            TCS([hftClimb[1:Ncruise] >= hftClimb[:Ncruise-1] + dhft]),
            TCS([hftClimb[0] >= dhft[0]]),
            hftClimb[-1] <= hftCruise,

            #compute the dh
            dhft == hftCruise/Nclimb,

            #constrain the thrust
            cls.climbP.engineP['thrust'] <= 2 * max(crs.cruiseP.engineP['thrust']),

            #set the range for each cruise segment, doesn't take credit for climb
            #down range disatnce covered
            crs.cruiseP['Rng'] == ReqRng/(Ncruise),

            #set the TSFC
            cls.climbP.engineP['TSFC'] == .7*units('1/hr'),
            crs.cruiseP.engineP['TSFC'] == .5*units('1/hr'),
            ])
        
        # Model.__init__(self, W_ftotal + s*units('N'), constraints + ac + cls + crs, subs)
        Model.__init__(self, W_ftotal, constraints + ac + cls + crs, subs)

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
    substitutions = {      
##            'V_{stall}': 120,
            '\\delta_P_{over}':12,
            'N_{land}': 6,
            'V_{NE}': 144,
            'SPR':8,
            'p_s':81.,
            'ReqRng': 500, #('sweep', np.linspace(500,2000,4)),
            'CruiseAlt': 30000, #('sweep', np.linspace(20000,40000,4)),
            'numeng': 2,
##            'W_{Load_max}': 6664,
            'n_{pass}': 150,
            'W_{avg. pass}': 180,
            'W_{carry on}': 15,
            'W_{checked}': 40,
##            'pax_{area}': 1,
            'w_{aisle}':0.51,
            'w_{seat}':0.5,
            'w_{sys}':0.1,
##            'C_{D_{fuse}}': .005, #assumes flat plate turbulent flow, from wikipedia
            'e'            : .9,
            'b_{max}'      : 35,
            'W_{cargo}'    : 10000,
            'r_E'          : 1, #[TAS]
            '\\lambda_{cone}':0.4,
#            '\\rho_{cone}' : 2700, #[TAS]
            '\\rho_{bend}' : 2700, #[TAS]
            '\\rho_{floor}':2700, #[TAS]
            '\\rho_{skin}' :2700, #[TAS]
            'W\'\'_{floor}': 60, #[TAS]
#            'W\'\'_{insul}':22, #[TAS]
            'W\'_{seat}'   :150, #[TAS]
#            'W\'_{window}' : 145.*3, #[TAS]
            'f_{fadd}'     :0.2, #[TAS]
            'f_{frame}'    :0.25, #[Philippe]        
            'f_{lugg,1}'   :0.4, #[Philippe]
            'f_{lugg,2}'   :0.1, #[Philippe]
            # 'f_{string}'   :0.25 #[Philippe]
            'f_{string}'   :0.1
                             
            }
           
    m = Mission(substitutions)
    #sol = m.solve(solver='mosek', verbosity = 4)
    bounds, sol = m.determine_unbounded_variables(m, solver="mosek",verbosity=3, iteration_limit=100)
