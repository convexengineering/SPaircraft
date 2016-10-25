""" d8fuselage.py """
""" 
Sources for substitutions:
#[b757 freight doc]
#[Boeing]
#[Philippe]
#[stdAtm]
#[TAS]

Other markers:
#[SP]
#[SPEquality]

"""
from numpy import pi
import numpy as np
import matplotlib.pyplot as plt
from gpkit import VectorVariable, Variable, Model, units, SignomialsEnabled
from gpkit import LinkedConstraintSet as LSC
#from gpkit.constraints.bounded import BoundedConstraintSet as BCS
from gpkit import SignomialEquality
from gpkit.tools import te_exp_minus1
from collections import defaultdict
from gpkit.small_scripts import mag

# Note that sweep has to be True for any sweep to take place. 
sweep         = False
nsweep        = 10

sweep_thetadb = False
thetadb_bounds=[0.0,0.5]

sweep_npass   = False
npass_bounds  = [160, 232]

sweep_Shtail   = False
Shtail_bounds  = [10,50]

sweep_fstring = True
fstring_bounds= [0.0,0.3]

class Ops(Model):
    def __init__(self,**kwargs):
        # Operational, arbitrary constraints listed here. 
        Nland = Variable('N_{land}',6.,'-', 'Emergency landing load factor') #[TAS]
        VNE          = Variable('V_{NE}',144,'m/s','Never-exceed speed') #[Philippe]
        rho0         = Variable(r'\rho_0', 1.225,'kg/m^3', 'Air density (0 ft)') #[stdAtm]
        rhoinf       = Variable('\\rho_{\\infty}',0.38,'kg/m^3','Air density (35,000ft)') #[stdAtm]
        Vinf         = Variable('V_{\\infty}',234, 'm/s','Cruise velocity')
        muinf        = Variable('\\mu_{\\infty}',1.46*10**-5, 'N*s/m^2', 'Dynamic viscosity (35,000 ft)')

        constraints = [Nland == Nland,
                        VNE == VNE,
                        rho0 == rho0,
                        rhoinf == rhoinf,
                        Vinf  == Vinf,
                        muinf == muinf]

        Model.__init__(self, None, constraints, **kwargs)


class Fuselage(Model):
    def __init__(self,ops,wingbox,htail,vtail,**kwargs):
        self.ops = ops
        self.wingbox = wingbox
        self.htail = htail
        self.vtail = vtail
        constraints = []
        g = 9.81*units('m*s**-2')

        # Will try to stick to Philippe's naming methods as closely as possible
        # for cross-compatibility (to be able to switch models quickly)

        # Fixed variables
        SPR          = Variable('SPR', 8, '-', 'Number of seats per row')
        nseat        = Variable('n_{seat}','-','Number of seats')
        dPover       = Variable('\\delta_P_{over-pressure}',12,'psi','Cabin overpressure')
        

        npass        = Variable('n_{pass}',192,'-', 'Number of passengers')
        nrows        = Variable('n_{rows}', '-', 'Number of rows')
        pitch        = Variable('p_s',81., 'cm', 'Seat pitch')
        
        # Cross sectional parameters
        Adb          = Variable('A_{db}', 'm^2', 'Web cross sectional area')
        Afloor       = Variable('A_{floor}', 'm^2', 'Floor beam x-sectional area')
        Afuse        = Variable('A_{fuse}', 'm^2', 'Fuselage x-sectional area')
        Ahold        = Variable('A_{hold}', 'm^2', 'Cargo hold x-sectional area')
        Askin        = Variable('A_{skin}', 'm^2', 'Skin cross sectional area')
        dh           = Variable('\\Delta_h',1, 'm','Distance from floor to widest part of fuselage') #[Philippe]
        hdb          = Variable('h_{db}','m', 'Web half-height')
        hfloor       = Variable('h_{floor}', 'm', 'Floor beam height')
        hfuse        = Variable('h_{fuse}','m','Fuselage height')
        #hhold        = Variable('h_{hold}', 'm', 'Height of the cargo hold')        
        Rfuse        = Variable('R_{fuse}', 'm', 'Fuselage radius') # will assume for now there: no under-fuselage extension deltaR
        tdb          = Variable('t_{db}', 'm', 'Web thickness')
        thetadb      = Variable('\\theta_{db}','-','DB fuselage joining angle')
        tshell       = Variable('t_{shell}', 'm', 'Shell thickness')
        tskin        = Variable('t_{skin}', 'm', 'Skin thickness')
        waisle       = Variable('w_{aisle}',0.51, 'm', 'Aisle width') #[Boeing]
        wdb          = Variable('w_{db}','m','DB added half-width')
        wfloor       = Variable('w_{floor}', 'm', 'Floor half-width')
        #wcargofloor  = Variable('w_{cargo}','m','Cargo floor half-width')
        wfuse        = Variable('w_{fuse}', 'm', 'Fuselage width')
        wseat        = Variable('w_{seat}',0.5,'m', 'Seat width') #[Philippe]
        wsys         = Variable('w_{sys}', 0.1,'m', 'Width between cabin and skin for systems') #[Philippe]
        
        # Tail cone variables
        lamcone      = Variable('\\lambda_{cone}',0.4, '-','Tailcone radius taper ratio (xshell2->xtail)')
        lcone        = Variable('l_{cone}', 'm', 'Cone length')
        plamv        = Variable('p_{\\lambda_v}',1.4,'-', '1 + 2*Tail taper ratio')
        tcone        = Variable('t_{cone}', 'm', 'Cone thickness')
        
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
        Vcargo       = Variable('V_{cargo}', 'm^3', 'Cargo volume')
        Vcone        = Variable('V_{cone}', 'm^3', 'Cone skin volume')
        Vcyl         = Variable('V_{cyl}', 'm^3', 'Cylinder skin volume')   
        Vdb          = Variable('V_{db}', 'm^3', 'Web volume')
        Vfloor       = Variable('V_{floor}', 'm^3', 'Floor volume')
        #Vhold        = Variable('V_{hold}', 'm^3', 'Hold volume')
        Vlugg        = Variable('V_{lugg}', 'm^3', 'Luggage volume')
        Vnose        = Variable('V_{nose}', 'm^3', 'Nose skin volume')
        
        # Weights
        #Wbuoy       = Variable('W_{buoy}', 'N', 'Buoyancy weight')
        Wapu         = Variable('W_{apu}', 'N', 'APU weight')
        Wavgpass     = Variable('W_{avg. pass}', 180, 'lbf', 'Average passenger weight') #[Philippe]
        Wcargo       = Variable('W_{cargo}',10000, 'N', 'Cargo weight') #[Philippe]
        Wcarryon     = Variable('W_{carry on}', 15, 'lbf', 'Ave. carry-on weight') #[Philippe]
        Wchecked     = Variable('W_{checked}', 40, 'lbf', 'Ave. checked bag weight') #[Philippe]
        Wcone        = Variable('W_{cone}', 'N', 'Cone weight')
        Wdb          = Variable('W_{db}' , 'N', 'Web weight')
        Wfix         = Variable('W_{fix}', 3000, 'lbf', 'Fixed weights (pilots, cockpit seats, navcom)') #[Philippe]
        Wfloor       = Variable('W_{floor}', 'N', 'Floor weight')
        Wfuse        = Variable('W_{fuse}', 'N', 'Fuselage weight')
        Winsul       = Variable('W_{insul}', 'N', 'Insulation material weight')
        Wlugg        = Variable('W_{lugg}', 'N', 'Passenger luggage weight')
        Wpay         = Variable('W_{pay}', 'N', 'Payload weight')
        Wpass        = Variable('W_{pass}', 'N', 'Passenger weight')
        Wpadd        = Variable('W_{padd}', 'N', 'Misc weights (galley, toilets, doors etc.)')
        Wseat        = Variable('W_{seat}', 'N', 'Seating weight')
        Wshell       = Variable('W_{shell}','N','Shell weight')
        Wskin        = Variable('W_{skin}', 'N', 'Skin weight')
        Wwindow      = Variable('W_{window}', 'N', 'Window weight')
        Wtail        = Variable('W_{tail}','N','Total tail weight')

        # Weight fractions and crago densities
        ffadd        = Variable('f_{fadd}',0.2, '-','Fractional added weight of local reinforcements') #[TAS]
        fapu         = Variable('f_{apu}',0.035,'-','APU weight as fraction of payload weight') #[TAS]
        #ffairing     = Variable('f_{fairing}',0.151,'-','  Fractional fairing weight')
        fframe       = Variable('f_{frame}',0.25,'-', 'Fractional frame weight') #[Philippe]
        flugg1       = Variable('f_{lugg,1}',0.4,'-','Proportion of passengers with one suitcase') #[Philippe]
        flugg2       = Variable('f_{lugg,2}',0.1, '-','Proportion of passengers with two suitcases') #[Philippe]
        fpadd        = Variable('f_{padd}',0.4, '-', 'Other misc weight as fraction of payload weight')
        fstring      = Variable('f_{string}','-','Fractional stringer weight')
        rhobend      = Variable('\\rho_{bend}', 2700, 'kg/m^3', 'Stringer density') #[TAS]
        rhocargo     = Variable('\\rho_{cargo}', 150, 'kg/m^3', 'Cargo density')  #[b757 freight doc]
        rholugg      = Variable('\\rho_{lugg}',100,'kg/m^3', 'Luggage density') #[Philippe]

        # Misc free variables
        
        # Material properties
        rE           = Variable('r_E', 1,'-', 'Ratio of stringer/skin moduli') #[TAS]
        rhocone      = Variable('\\rho_{cone}',2700,'kg/m^3','Cone material density') #[TAS]
        rhofloor     = Variable('\\rho_{floor}',2700, 'kg/m^3', 'Floor material density') #[TAS]
        rhoskin      = Variable('\\rho_{skin}',2700,'kg/m^3', 'Skin density') #[TAS]
        Wppfloor     = Variable('W\'\'_{floor}', 60,'N/m^2', 'Floor weight/area density') #[TAS]
        Wppinsul     = Variable('W\'\'_{insul}',22,'N/m^2', 'Weight/area density of insulation material') #[TAS]
        Wpseat       = Variable('W\'_{seat}',150,'N', 'Weight per seat') #[TAS]
        Wpwindow     = Variable('W\'_{window}', 145.*3,'N/m', 'Weight/length density of windows') #[TAS]
        
        # Loads
        #Pcargofloor = Variable ('P_{cargo floor}','N','Distributed cargo floor load')
        Mfloor       = Variable('M_{floor}', 'N*m', 'Max bending moment in floor beams')
        Pfloor       = Variable('P_{floor}','N', 'Distributed floor load')
        Sfloor       = Variable('S_{floor}', 'N', 'Maximum shear in floor beams')
        sigfloor     = Variable('\\sigma_{floor}',30000/0.000145, 'N/m^2', 'Max allowable floor stress') #[TAS]
        sigskin      = Variable('\\sigma_{skin}', 15000/0.000145,'N/m^2', 'Max allowable skin stress') #[TAS] 
        sigth        = Variable('\\sigma_{\\theta}', 'N/m^2', 'Skin hoop stress')
        sigx         = Variable('\\sigma_x', 'N/m^2', 'Axial stress in skin')
        taucone      = Variable('\\tau_{cone}', 'N/m^2', 'Shear stress in cone')
        taufloor     = Variable('\\tau_{floor}',30000/0.000145, 'N/m^2', 'Max allowable shear web stress') #[TAS]
        
        Mhmax        = Variable('M_{h_max}', 'N*m','Maximum horizontal axis bending moment')
        Mvmax        = Variable('M_{v_max}', 'N*m','Maximum vertical axis bending moment')
        Mhaero       = Variable('M_{h_aero}','N*m','Maximum horizontal tail aero bending load')
        Mvaero       = Variable('M_{v_aero}','N*m','Maximum vertical tail aero bending load')
        
        # Bending inertias (ported from TASOPT)
        A0           = Variable('A0','m^2','Horizontal bending area constant A0') #(shell inertia contribution)
        A1           = Variable('A1','m','Horizontal bending area constant A1') #(tail impact + aero loading)
        A2           = Variable('A2','-','Horizontal bending area constant A2') #(fuselage impact)
        Ahbendb      = Variable('A_{hbendb}','m^2','Horizontal bending area at rear wingbox')
        Ahbendf      = Variable('A_{hbendf}','m^2','Horizontal bending area at front wingbox')
        Avbendb      = Variable('A_{vbendb}','m^2','Vertical bending material area at rear wingbox')
        #B0           = Variable('B0','m^2','Vertical bending area constant B0') #(shell inertia contribution)
        #B1           = Variable('B1','m','Vertical bending area constant B1') #(vertical tail bending load)
        Ihshell      = Variable('I_{hshell}','m^4','Shell horizontal bending inertia')
        #Ivshell      = Variable('I_{vshell}','m^4','Shell vertical bending inertia')
        rMh          = Variable('r_{M_h}',.4,'-','Horizontal inertial relief factor') #[TAS]
        rMv          = Variable('r_{M_v}',.7,'-','Vertical inertial relief factor') #[TAS]
        sigbend      = Variable('\\sigma_{bend}','N/m^2','Bending material stress')
        sigMh        = Variable('\\sigma_{M_h}','N/m^2','Horizontal bending material stress')
        #sigMv        = Variable('\\sigma_{M_v}','N/m^2','Vertical bending material stress')
        Vhbend       = Variable('V_{hbend}','m^3','Horizontal bending material volume')
        Vhbendb      = Variable('V_{hbendb}','m^3','Horizontal bending material volume b') #back fuselage
        Vhbendc      = Variable('V_{hbendc}','m^3','Horizontal bending material volume c') #center fuselage
        Vhbendf      = Variable('V_{hbendf}','m^3','Horizontal bending material volume f') #front fuselage
        #Vvbend       = Variable('V_{vbend}','m^3','Vertical bending material volume')
        #Vvbendb      = Variable('V_{vbendb}','m^3','Vertical bending material volume b') #back fuselage
        #Vvbendc      = Variable('V_{vbendc}','m^3','Vertical bending material volume c') #center fuselage
        Whbend       = Variable('W_{hbend}','N','Horizontal bending material weight')
        #Wvbend       = Variable('W_{vbend}','N','Vertical bending material weight')
        xhbend       = Variable('x_{hbend}','m','Horizontal zero bending location')
        #xvbend       = Variable('x_{vbend}','m','Vertical zero bending location')
        # x-location variables
        xshell1      = Variable('x_{shell1}', 'm', 'Start of cylinder section')
        xshell2      = Variable('x_{shell2}', 'm', 'End of cylinder section')
        xtail        = Variable('x_{tail}','m', 'x-location of tail')
        
        # Drag variables
        # Dfuse    = Variable('D_{fuse}', 'N', 'Total drag in cruise')
        # Dfrict   =riable('D_{friction}', 'N', 'Friction drag')
        # Dupswp   = Variable('D_{upsweep}', 'N', 'Drag due to fuse upsweep')
        # f        = Variable('f', '-', 'Fineness ratio')
        # FF       = Variable('FF', '-','Fuselage form factor')

        with SignomialsEnabled():
            constraints = [
            # Passenger constraints (assuming 737-sixed aircraft)
            #Temporarily
            Wpass    == npass*Wavgpass,
            Wpay     >= Wpass + Wlugg + Wcargo,
            nseat    == npass,
            nrows    == nseat/SPR,
            lshell   == nrows*pitch,

            
            # Fuselage joint angle relations
            thetadb     == wdb/Rfuse, # first order Taylor works...
            thetadb     >= 0.05, thetadb <= 0.5, #Temporarily
            hdb         >= Rfuse*(1.0-.5*thetadb**2), #[SP]
            
            # Fuselage cross-sectional relations
            Askin       >= (2*pi + 4*thetadb)*Rfuse*tskin + Adb, #no delta R for now
            Adb         == (2*hdb)*tdb,
            Afuse       >= (pi + 2*thetadb + thetadb)*Rfuse**2,
            SignomialEquality(tshell,tskin*(1+rE*fstring*rhoskin/rhobend)),
            #tshell     >= tskin*(1+rE*fstring*rhoskin/rhobend),
            #Rfuse       >= dh + hhold + hfloor,
            #fstring     == .1, #Temporarily, so I can see changes in horizontal bending material
            wfuse       >= SPR*wseat + 2*waisle + tdb + 2*tshell + 2*wsys,
            wfuse       <= 2*(Rfuse + wdb),
            wfloor      >= .5*wfuse, # half of the total floor width in fuselage
            #wcargofloor == 2*Rfuse*thetadb, #Temporarily #Approx
            #Ahold       >= (wfloor + wcargofloor)*hhold, #[SP]
            # Added synthetic constraint on hfloor to keep it from growing too large #Temporarily
            hfloor      <= .1*Rfuse,

            # Fuselage surface area relations
            Snose    >= (2*pi + 4*thetadb)*Rfuse**2 *(1/3 + 2/3*(lnose/Rfuse)**(8/5))**(5/8),
            Sbulk    >= (2*pi + 4*thetadb)*Rfuse**2,
            
            # Fuselage length relations
            lfuse    >= lnose+lshell+lcone, 
            lnose    == 0.3*lshell, # Temporarily
            xshell1  == lnose,
            xshell2  >= lnose + lshell,
            lcone    == Rfuse/lamcone,            
            
            # Fuselage volume relations
            Vcyl   == Askin*lshell,
            Vnose  == Snose*tskin,
            Vbulk  == Sbulk*tskin,
            Vdb    == Adb*lshell,
            Vcabin >= Afuse*(lshell + 0.67*lnose + 0.67*Rfuse),
            #Vhold  >= Vcargo + Vlugg,
            #Vhold == Ahold*lshell,
            
            # Fuselage weight relations
            Wapu     == Wpay*fapu,
            Wcargo   == Vcargo*g*rhocargo,
            Wdb      == rhoskin*g*Vdb,
            Winsul   >= Wppinsul*((1.1*pi+2*thetadb)*Rfuse*lshell + 0.55*(Snose+Sbulk)),
            Wlugg    >= flugg2*npass*2*Wchecked + flugg1*npass*Wchecked + Wcarryon,
            Wlugg    == Vlugg*g*rholugg,
            Wwindow  >= Wpwindow*lshell,
            Wpadd    == Wpay*fpadd,
            Wseat    == Wpseat*nseat,
            Wskin    >= rhoskin*g*(Vcyl + Vnose + Vbulk),
            Wshell   >= Wskin*(1 + fstring + ffadd + fframe) + Wdb + Whbend, #+ Wvbend,
            #Wfuse    >= Wapu + Wcargo + Wlugg + Wfix + Winsul + Wshell + Wfloor + Wwindow + Whbend + Wvbend + Wtail,
            Wfuse    >= Wfix + Wapu + Wpadd + Wseat + Wshell + Wwindow + Winsul + Wfloor + Wtail,

            ## Stress relations
            #Pressure shell loading
            tskin    == dPover*Rfuse/sigskin,
            tdb      == 2*dPover*wdb/sigskin,
            sigx     == dPover*Rfuse/(2*tshell),
            sigth    == dPover*Rfuse/tskin,
            
            # Floor loading 
            lfloor   >= lshell + 2*Rfuse,            
            Pfloor   >= self.ops['N_{land}']*(Wpay + Wseat),
            Mfloor   == 9./256.*Pfloor*wfloor,
            Afloor   >= 2.*Mfloor/(sigfloor*hfloor) + 1.5*Sfloor/taufloor,
            Vfloor   == 2*wfloor*Afloor,
            Wfloor   >= rhofloor*g*Vfloor + 2*wfloor*lfloor*Wppfloor,
            Sfloor   == (5./16.)*Pfloor,

            # Tail cone sizing
            taucone                         == sigskin,
            3*self.vtail['Q_v']*(plamv-1)                  >= self.vtail['L_{v_{max}}']*self.vtail['b_{vt}']*(plamv),
            Vcone*(1+lamcone)*(pi+4*thetadb)>= self.vtail['Q_v']/taucone*(pi+2*thetadb)*(lcone/Rfuse)*2,
            Wcone                           >= rhocone*g*Vcone*(1+fstring+fframe),
            Wtail                           >= self.vtail['W_{vtail}'] + self.htail['W_{htail}'] + Wcone,
            
            # Tail aero loads
            xtail    >= lnose + lshell + .5*lcone, #Temporarily
            #Mhaero   >= rMh*Lhmax*(xtail-xwing), #[SP]
            #Mvaero   >= rMv*Lvmax*(xtail-xwing), #[SP]

            hfuse    == Rfuse, # may want to consider adding deltaRfuse later...
            xhbend   >= self.wingbox['x_{wing}'],

            
            # Horizontal bending model
            # Maximum axial stress is the sum of bending and pressurization stresses
            Ihshell <= ((pi+4*thetadb)*Rfuse**2)*Rfuse*tshell + 2/3*hdb**3*tdb, # [SP]
            #Ivshell <= (pi*Rfuse**2 + 8*wdb*Rfuse + (2*pi+4*thetadb)*wdb**2)*Rfuse*tshell, #[SP] #Ivshell approximation needs to be improved
            sigbend == rE*sigskin,
        
            # Horizontal bending material model
            # Calculating xbend, the location where additional bending material is required
            SignomialEquality(A0,A2*(xshell2-xhbend)**2 + A1*(xtail-xhbend)), #[SP] #[SPEquality] 
            A2      >=  self.ops['N_{land}']*(Wpay+Wshell+Wwindow+Winsul+Wfloor+Wseat)/(2*lshell*hfuse*sigMh), # Landing loads constant A2
            A1      >= (self.ops['N_{land}']*Wtail + rMh*self.htail['L_{h_{max}}'])/(hfuse*sigMh),                                # Aero loads constant A1
            A0      == (Ihshell/(rE*hfuse**2)),                                                # Shell inertia constant A0
            Ahbendf >= A2*(xshell2-self.wingbox['x_f'])**2 + A1*(xtail-self.wingbox['x_f']) - A0, #[SP]                           # Bending area forward of wingbox
            Ahbendb >= A2*(xshell2-self.wingbox['x_b'])**2 + A1*(xtail-self.wingbox['x_b']) - A0, #[SP]                           # Bending area behind wingbox

            Vhbendf >= A2/3*((xshell2-self.wingbox['x_f'])**3 - (xshell2-xhbend)**3) \
                            + A1/2*((xtail-self.wingbox['x_f'])**2 - (xtail - xhbend)**2) \
                            + A0*(xhbend-self.wingbox['x_f']), #[SP]

            Vhbendb >= A2/3*((xshell2-self.wingbox['x_b'])**3 - (xshell2-xhbend)**3) \
                            + A1/2*((xtail-self.wingbox['x_b'])**2 - (xtail - xhbend)**2) \
                            + A0*(xhbend-self.wingbox['x_b']), #[SP]
            Vhbendc >= .5*(Ahbendf + Ahbendb)*self.wingbox['c_0']*self.wingbox['\\bar_w'],
            Vhbend  >= Vhbendc + Vhbendf + Vhbendb,
            Whbend  >= g*rhobend*Vhbend,
            
            # Vertical bending material model
            # Calculating xvbend, the location where additional bending material is required
            #xvbend  >= self.wingbox['x_{wing}'],
            #SignomialEquality(B0,B1*(xtail-xvbend)), #[SP] #[SPEquality]
            #B1      == rMv*self.vtail['L_{v_{max}}']/(wfloor*sigMv),                                               # Aero loads constant B1
            #B0      == Ivshell/(rE*wfloor**2),                                                 # Shell inertia constant B0
            #Avbendb >= B1*(xtail-self.wingbox['x_b']) - B0,                                                      # Bending area behind wingbox
            
            #Vvbendb >= .5*B1*((xtail-self.wingbox['x_b'])**2 - (xtail-xvbend)**2) - B0*(xvbend - self.wingbox['x_b']), #[SP]
            #Vvbendc >= .5*Avbendb*self.wingbox['c_0']*self.wingbox['\\bar_w'],
            #Vvbend  >= Vvbendb + Vvbendc,
            #Wvbend  >= g*rhobend*Vvbend,

            # Temporary wing variable substitutions
            self.wingbox['c_0']       == 0.1*lshell, #Temporarily
            self.wingbox['dx_{wing}']   == 0.25*self.wingbox['c_0'], #Temporarily
                        
            sigMh   <= sigbend - rE*dPover/2*Rfuse/tshell, # The stress available to the bending material reduced because of pressurization
            #sigMv   == sigMh,

            # Drag model
            # f == lfuse/((4/np.pi*Afuse)**0.5), # fineness ratio
            # FF >= 1 + 60/f**3 + f/400, # form factor
            # Dfrict >= FF * (pi*Rfuse+2*wdb) * muinf*Vinf* 0.074*(rhoinf*Vinf*lfuse/muinf)**0.8,
            # Dfuse >= Dfrict

            # Wingbox/fuselage constraints
            SignomialEquality(self.wingbox['x_{wing}'], lnose + 0.6*lshell),
            ]

        Model.__init__(self, None, constraints, **kwargs)

class HTail(Model):
    def dynamic(self,state):
        return HTailP(self,state)

    def __init__(self,ops,**kwargs):
        self.ops = ops
        Whtail       = Variable('W_{htail}',10000, 'N', 'Horizontal tail weight') #Temporarily
        Lhmax        = Variable('L_{h_{max}}',35000,'N', 'Max horizontal tail load')
        Shtail       = Variable('S_{htail}',32*0.8,'m^2','Horizontal tail area') #Temporarily
        CLhmax       = Variable('C_{L_{h_{max}}}', 2.5, '-', 'Max lift coefficient') #Temporarily
        constraints = [#Lhmax    == 0.5*self.ops['\\rho_{\\infty}']*self.ops['V_{NE}']**2*Shtail*CLhmax,
                       Lhmax    == Lhmax,
                       Whtail   == Whtail,
                       Shtail   == Shtail,
                       CLhmax   == CLhmax]
        Model.__init__(self, None, constraints, **kwargs)

class VTail(Model):
    def dynamic(self,state):
        return VTailP(self,state)

    def __init__(self,ops,**kwargs):
        bvt          = Variable('b_{vt}',7, 'm', 'Vertical tail span')
        Lvmax        = Variable('L_{v_{max}}',35000,'N', 'Max vertical tail load')
        Wvtail       = Variable('W_{vtail}',10000, 'N', 'Vertical tail weight') #Temporarily
        Qv           = Variable('Q_v', 'N*m', 'Torsion moment imparted by tail')

        constraints = [bvt == bvt, 
                       Lvmax == Lvmax,
                       Wvtail == Wvtail,
                       Qv == Qv]
        Model.__init__(self, None, constraints, **kwargs)


# class Tail(Model):
#     def dynamic(self,state):
#         return TailP(self,state)

#     def __init__(self,**kwargs):

class Wing(Model):
    def dynamic(self,state):
        return WingP(self,state)

    def __init__(self,ops,**kwargs):
        self.ops = ops
        constraints = [];
        Model.__init__(self,None,constraints,**kwargs)

class WingBox(Model):
    def dynamic(self,state):
        return WingBoxP(self,state)

    def __init__(self,ops,**kwargs):
        self.ops = ops
        xf           = Variable('x_f','m','x-location of front of wingbox')
        xb           = Variable('x_b','m','x-location of back of wingbox')
        c0           = Variable('c_0','m','Root chord of the wing')
        wbar         = Variable('\\bar_w',0.5,'-','Wingbox to chord ratio') #Temporarily
        xwing        = Variable('x_{wing}','m', 'x-location of wing')
        dxwing       = Variable('dx_{wing}','m','wing box offset')
        # Setting bending area integration bounds (defining wing box locations)
        with SignomialsEnabled():
            constraints  = [SignomialEquality(xf,xwing + dxwing + .5*c0*wbar), #[SP] [SPEquality]
                        #xf >= xwing + dxwing + .5*c0*wbar,
                        SignomialEquality(xb, xwing - dxwing + .5*c0*wbar), #[SP] [SPEquality]
                        #xb <= xwing - dxwing + .5*c0*wbar, #[SP]        
                        ];
        Model.__init__(self,None,constraints,**kwargs)

class Aircraft(Model):
    "The D8 Double Bubble"

    def dynamic(self,state):
        """Creates an instance of this component's performance model,
        given a state"""
        return AircraftP(self,state)

    def __init__(self,**kwargs):
        self.ops = Ops()
        self.wing    = Wing(self.ops)
        self.wingbox = WingBox(self.ops)
        self.htail   = HTail(self.ops)
        self.vtail   = VTail(self.ops)
        self.fuse    = Fuselage(self.ops,self.wingbox,self.htail,self.vtail)

        self.components = [self.fuse, self.wing, self.wingbox, self.htail, self.vtail] 

        with SignomialsEnabled():
            constraints = []
        objective = self.fuse["W_{fuse}"] + \
                     self.fuse["V_{cabin}"]*units('N/m^3') + \
                     self.fuse["t_{shell}"]*units('N/m') + \
                     self.fuse["l_{fuse}"]*units('N/m')
                     
        Model.__init__(self, objective, self.components + constraints, **kwargs)

        # # Free variables
        # Wfuse   = Variable('W_{fuse}', 'N', 'Fuselage weight')
        # Wlg     = Variable('W_{lg}', 'N', 'Landing gear weight')
        # Wvt     = Variable('W_{vt}', 'N', 'Vertical tail weight')
        # xCG     = Variable('x_{CG}', 'm', 'x-location of CG')
        # xCGfu   = Variable('x_{CG_{fu}}', 'm', 'x-location of fuselage CG')
        # xCGlg   = Variable('x_{CG_{lg}}', 'm', 'x-location of landing gear CG')
        # xCGvt   = Variable('x_{CG_{vt}}', 'm', 'x-location of tail CG') 

        # # Fixed variables (pulled from Philippe's aircraft model)
        # Weng    = Variable('W_{eng}', 10000, 'N', 'Engine weight')
        # Wht     = Variable('W_{ht}', 5000, 'N', 'Horizontal tail weight')
        # Wwing   = Variable('W_{wing}', 30000, 'N', 'Wing weight')
        # xCGeng  = Variable('x_{CG_{eng}}', 15, 'm', 'x-location of engine CG')
        # xCGht   = Variable('x_{CG_{ht}}', 38, 'm', 'x-location of horizontal tail CG')
        # xCGwing = Variable('x_{CG_{wing}}', 15, 'm', 'x-location of wing CG')

if __name__ == "__main__":
    M = Aircraft()
    #M = Model(M.cost, BCS(M))
    
    #M.substitutions.update({'f_{string}':0.1})
    #bounds, sol = M.determine_unbounded_variables(M, solver="mosek",verbosity=2, iteration_limit=100)
    sol = M.localsolve("mosek",tolerance = 0.01, verbosity = 1, iteration_limit=50)
    varVals = sol['variables']
