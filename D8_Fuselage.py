from gpkit import Model, Variable, units, SignomialsEnabled, SignomialEquality
from gpkit.constraints.tight import Tight as TCS
from numpy import pi
import numpy as np

"""
D8 double bubble fuselage files
"""

class Fuselage(Model):
    '''
    A double-bubble fuselage model
    '''

    def setup(self, **kwargs):
        g = Variable('g',9.81,'m*s^-2','Acceleration due to gravity')
        dPover = Variable('\\delta_P_{over}', 'psi', 'Cabin overpressure')
        npax = Variable('n_{pax}', '-', 'Number of Passengers to Carry')
        Nland = Variable('N_{land}', 6.0, '-',
                         'Emergency landing load factor')  # [TAS]
        Nlift = Variable('N_{lift}','-','Wing maximum load factor')
        SPR = Variable('SPR', '-', 'Number of seats per row')
        nrows = Variable('n_{rows}', '-', 'Number of rows')
        nseat = Variable('n_{seat}', '-', 'Number of seats')
        pitch = Variable('p_s', 'cm', 'Seat pitch')

        # Cross-sectional variables
        Adb = Variable('A_{db}', 'm^2', 'Web cross sectional area')
        Afloor = Variable('A_{floor}', 'm^2', 'Floor beam x-sectional area')
        Afuse = Variable('A_{fuse}', 'm^2', 'Fuselage x-sectional area')
        Askin = Variable('A_{skin}', 'm^2', 'Skin cross sectional area')
        hdb = Variable('h_{db}', 'm', 'Web half-height')
        hfloor = Variable('h_{floor}', 'm', 'Floor beam height')
        hfuse = Variable('h_{fuse}', 'm', 'Fuselage height')
        dRfuse = Variable('\\delta R_{fuse}','m','Fuselage extension height')
        Rfuse = Variable('R_{fuse}', 'm', 'Fuselage radius')
        tdb = Variable('t_{db}', 'm', 'Web thickness')
        thetadb = Variable('\\theta_{db}', '-', 'DB fuselage joining angle')
        tshell = Variable('t_{shell}', 'm', 'Shell thickness')
        tskin = Variable('t_{skin}', 'm', 'Skin thickness')
        waisle = Variable('w_{aisle}', 'm', 'Aisle width')
        wdb = Variable('w_{db}', 'm', 'DB added half-width')
        wfloor = Variable('w_{floor}', 'm', 'Floor half-width')
        wfuse = Variable('w_{fuse}', 'm', 'Fuselage half-width')
        wseat = Variable('w_{seat}', 'm', 'Seat width')
        wsys = Variable('w_{sys}', 'm', 'Width between cabin and skin for systems')

        # Tail cone variables
        lamcone = Variable('\\lambda_{cone}', '-', 'Tailcone radius taper ratio')
        lcone = Variable('l_{cone}', 'm', 'Cone length')
        plamv = Variable('p_{\\lambda_v}', 1.6, '-', '1 + 2*Tail taper ratio')
        # tcone = Variable('t_{cone}', 'm', 'Cone thickness') # perhaps to be added later

        # Lengths
        c0 = Variable('c_0', 'm', 'Root chord of the wing')
        lfuse = Variable('l_{fuse}', 'm', 'Fuselage length')
        lnose = Variable('l_{nose}', 'm', 'Nose length')
        lshell = Variable('l_{shell}', 'm', 'Shell length')
        lfloor = Variable('l_{floor}', 'm', 'Floor length')

        # Surface areas
        Sbulk = Variable('S_{bulk}', 'm^2', 'Bulkhead surface area')
        Snose = Variable('S_{nose}', 'm^2', 'Nose surface area')

        # Volumes
        Vbulk = Variable('V_{bulk}', 'm^3', 'Bulkhead skin volume')
        Vcabin = Variable('V_{cabin}', 'm^3', 'Cabin volume')
        Vcone = Variable('V_{cone}', 'm^3', 'Cone skin volume')
        Vcyl = Variable('V_{cyl}', 'm^3', 'Cylinder skin volume')
        Vdb = Variable('V_{db}', 'm^3', 'Web volume')
        Vfloor = Variable('V_{floor}', 'm^3', 'Floor volume')
        Vnose = Variable('V_{nose}', 'm^3', 'Nose skin volume')

        # Loads
        sigskin = Variable('\\sigma_{skin}', 'N/m^2', 'Max allowable skin stress')
        sigth = Variable('\\sigma_{\\theta}', 'N/m^2', 'Skin hoop stress')
        sigx = Variable('\\sigma_x', 'N/m^2', 'Axial stress in skin')

        # Floor loads
        Mfloor = Variable('M_{floor}', 'N*m',
                          'Max bending moment in floor beams')
        Pfloor = Variable('P_{floor}', 'N', 'Distributed floor load')
        Sfloor = Variable('S_{floor}', 'N', 'Maximum shear in floor beams')
        sigfloor = Variable('\\sigma_{floor}', 'N/m^2', 'Max allowable floor stress')
        taucone = Variable('\\tau_{cone}', 'N/m^2', 'Shear stress in cone')
        taufloor = Variable('\\tau_{floor}', 'N/m^2', 'Max allowable shear web stress')

        # Bending inertias (ported from TASOPT)
        # (shell inertia contribution)
        A0h = Variable('A0h', 'm^2', 'Horizontal bending area constant A0h')

        # (tail impact + aero loading)
        A1hLand = Variable('A1h_{Land}', 'm', 'Horizontal bending area constant A1h (landing case)')
        A1hMLF = Variable('A1h_{MLF}', 'm', 'Horizontal bending area constant A1h (max aero load case)')

        # (fuselage impact)
        A2hLand = Variable('A2h_{Land}', '-', 'Horizontal bending area constant A2h (landing case)')
        A2hMLF = Variable('A2h_{MLF}', '-', 'Horizontal bending area constant A2h (max aero load case)')

        AhbendbLand = Variable('A_{hbendb_{Land}}', 'm^2','Horizontal bending area at rear wingbox (landing case)')
        AhbendbMLF = Variable('A_{hbendb_{MLF}}', 'm^2','Horizontal bending area at rear wingbox (max aero load case)')

        AhbendfLand = Variable('A_{hbendf_{Land}}', 'm^2', 'Horizontal bending area at front wingbox (landing case)')
        AhbendfMLF = Variable('A_{hbendf_{MLF}}', 'm^2', 'Horizontal bending area at front wingbox (max aero load case)')

        Avbendb = Variable('A_{vbendb}', 'm^2', 'Vertical bending material area at rear wingbox')

        B0v           = Variable('B0v','m^2','Vertical bending area constant B0') #(shell inertia contribution)
        B1v         = Variable('B1v','m','Vertical bending area constant B1')
        # #(vertical tail bending load)
        Ihshell = Variable('I_{hshell}', 'm^4',
                           'Shell horizontal bending inertia')
        Ivshell      = Variable('I_{vshell}','m^4','Shell vertical bending inertia')
        rMh = Variable('r_{M_h}', .4, '-','Horizontal inertial relief factor')  # [TAS]
        rMv = Variable('r_{M_v}', .7, '-','Vertical inertial relief factor')  # [TAS]
        sigbend = Variable('\\sigma_{bend}', 'N/m^2',
                           'Bending material stress')
        sigMh = Variable('\\sigma_{M_h}', 'N/m^2',
                         'Horizontal bending material stress')
        sigMv        = Variable('\\sigma_{M_v}','N/m^2','Vertical bending material stress')
        Vhbend = Variable('V_{hbend}', 'm^3',
                          'Horizontal bending material volume')

        Vhbendb = Variable('V_{hbendb}', 'm^3', 'Horizontal bending material volume b')  # back fuselage
        Vhbendc = Variable('V_{hbendc}', 'm^3', 'Horizontal bending material volume c') # center fuselage
        Vhbendf = Variable('V_{hbendf}', 'm^3','Horizontal bending material volume f') # front fuselage

        Vvbend       = Variable('V_{vbend}','m^3','Vertical bending material volume')
        Vvbendb      = Variable('V_{vbendb}','m^3','Vertical bending material volume b') #back fuselage
        Vvbendc      = Variable('V_{vbendc}','m^3','Vertical bending material volume c') #center fuselage

        Whbend = Variable('W_{hbend}', 'lbf','Horizontal bending material weight')
        Wvbend = Variable('W_{vbend}','lbf','Vertical bending material weight')

        xhbendLand = Variable('x_{hbend_{Land}}', 'ft', 'Horizontal zero bending location (landing case)')
        xhbendMLF = Variable('x_{hbend_{MLF}}', 'ft', 'Horizontal zero bending location (maximum aero load case)')
        xvbend       = Variable('x_{vbend}','ft','Vertical zero bending location')

        # Material properties
        rE = Variable('r_E', 1., '-', 'Ratio of stringer/skin moduli')  # [TAS] # [b757 freight doc]
        rhocargo = Variable('\\rho_{cargo}', 150, 'kg/m^3', 'Cargo density')
        rhocone = Variable('\\rho_{cone}', 'kg/m^3',
                           'Cone material density')  # [TAS]
        rhobend = Variable('\\rho_{bend}', 'kg/m^3',
                           'Stringer density')  # [TAS]
        rhofloor = Variable('\\rho_{floor}', 'kg/m^3',
                            'Floor material density')  # [TAS]
        rholugg = Variable('\\rho_{lugg}', 100,
                           'kg/m^3', 'Luggage density')  # [Philippe]
        rhoskin = Variable('\\rho_{skin}', 'kg/m^3', 'Skin density')  # [TAS]
        Wppfloor = Variable('W\'\'_{floor}', 'N/m^2',
                            'Floor weight/area density')  # [TAS]
        Wppinsul = Variable(
            'W\'\'_{insul}', 'N/m^2', 'Weight/area density of insulation material')  # [TAS]
        Wpseat = Variable('W\'_{seat}', 'N', 'Weight per seat')  # [TAS]
        Wpwindow = Variable('W\'_{window}', 'N/m',
                            'Weight/length density of windows')  # [TAS]

        # Weight fractions
        fapu = Variable('f_{apu}', 0.035, '-',
                        'APU weight as fraction of payload weight')  # [TAS]
        ffadd = Variable(
            'f_{fadd}', '-', 'Fractional added weight of local reinforcements')  # [TAS]
        fframe = Variable('f_{frame}', '-',
                          'Fractional frame weight')  # [Philippe]
        flugg1 = Variable(
            'f_{lugg,1}', '-', 'Proportion of passengers with one suitcase')  # [Philippe]
        flugg2 = Variable(
            'f_{lugg,2}', '-', 'Proportion of passengers with two suitcases')  # [Philippe]
        fpadd = Variable('f_{padd}', 0.4, '-',
                         'Other misc weight as fraction of payload weight')
        fseat = Variable('f_{seat}','-','Fractional seat weight')
        fstring = Variable('f_{string}', '-',
                           'Fractional stringer weight')  # [Philippe]

        # Weights
        Wapu = Variable('W_{apu}', 'lbf', 'APU weight')
        Wavgpass = Variable('W_{avg. pass}', 'lbf',
                            'Average passenger weight')  # [Philippe]
        Wcargo = Variable('W_{cargo}', 'lbf', 'Cargo weight')  # [Philippe]
        Wcarryon = Variable('W_{carry on}', 'lbf',
                            'Ave. carry-on weight')  # [Philippe]
        Wchecked = Variable('W_{checked}', 'lbf',
                            'Ave. checked bag weight')  # [Philippe]
        Wcone = Variable('W_{cone}', 'lbf', 'Cone weight')
        Wdb = Variable('W_{db}', 'lbf', 'Web weight')
        Wfix = Variable(
            'W_{fix}', 'lbf', 'Fixed weights (pilots, cockpit seats, navcom)')
        Wfloor = Variable('W_{floor}', 'lbf', 'Floor weight')
        Wfuse = Variable('W_{fuse}', 'lbf', 'Fuselage weight')
        Winsul = Variable('W_{insul}', 'lbf', 'Insulation material weight')
        Wlugg = Variable('W_{lugg}', 'lbf', 'Passenger luggage weight')
        Wpadd = Variable('W_{padd}', 'lbf',
                         'Misc weights (galley, toilets, doors etc.)')
        Wpax = Variable('W_{pax}', 'lbf', 'Passenger weight')
        Wpay = Variable('W_{payload}', 'lbf', 'Payload weight')
        Wseat = Variable('W_{seat}', 'lbf', 'Seating weight')
        Wshell = Variable('W_{shell}', 'lbf', 'Shell weight')
        Wskin = Variable('W_{skin}', 'lbf', 'Skin weight')
        Wtail = Variable('W_{tail}', 'lbf', 'Total tail weight')
        Wwindow = Variable('W_{window}', 'lbf', 'Window weight')

        # x-location variables
        xshell1 = Variable('x_{shell1}', 'm', 'Start of cylinder section')
        xshell2 = Variable('x_{shell2}', 'm', 'End of cylinder section')
        xtail = Variable('x_{tail}', 'm', 'x-location of tail')
        xwing = Variable('x_{wing}', 'm', 'x-location of wing')

        # Wingbox variables
        xf = Variable('x_f', 'm', 'x-location of front of wingbox')
        xb = Variable('x_b', 'm', 'x-location of back of wingbox')
        w = Variable('wtc', 0.5, '-', 'Wingbox-width-to-chord ratio')

        constraints = []
        with SignomialsEnabled():
            constraints.extend([

                # Passenger constraints
                Wlugg >= flugg2 * npax * 2 * Wchecked + flugg1 * npax * Wchecked + Wcarryon,
                Wpax == npax * Wavgpass,
                Wpay >= Wpax + Wlugg + Wcargo,
                nseat == npax,
                nrows == nseat / SPR,
                lshell == nrows * pitch,

                # Fuselage joint angle relations
                thetadb == wdb / Rfuse,  # first order Taylor works...
                hdb >= Rfuse * (1.0 - .5 * thetadb**2),  # [SP]

                # Cross-sectional constraints
                Adb >= (2 * hdb + dRfuse) * tdb,
                Afuse >= (pi + 2 * thetadb + 2 * thetadb * \
                          (1 - thetadb**2 / 2)) * Rfuse**2 + 2*dRfuse*Rfuse,  # [SP]
                Askin >= (2 * pi + 4 * thetadb) * Rfuse * tskin + 2*dRfuse*tskin,
                wfloor == wfuse,
                wfuse <= (Rfuse + wdb),
                SignomialEquality(hfuse, Rfuse + 0.5*dRfuse), #[SP] #[SPEquality]
                TCS([tshell <= tskin * (1. + rE * fstring * rhoskin / rhobend)]), #[SP]

                # Fuselage surface area relations
                Snose >= (2 * pi + 4 * thetadb) * Rfuse**2 * \
                (1 / 3 + 2 / 3 * (lnose / Rfuse)**(8 / 5))**(5 / 8),
                Sbulk >= (2 * pi + 4 * thetadb) * Rfuse**2,

                # Fuselage length relations
                SignomialEquality(lfuse, lnose + lshell + lcone),  
                lcone == Rfuse / lamcone,
                xshell1 == lnose,
                TCS([xshell2 >= lnose + lshell]), 
                # STRESS RELATIONS
                # Pressure shell loading
                tskin == dPover * Rfuse / sigskin,
                tdb == 2 * dPover * wdb / sigskin,
                sigx == dPover * Rfuse / (2 * tshell),
                sigth == dPover * Rfuse / tskin,

                # Floor loading
                lfloor >= lshell + 2 * Rfuse,
                Pfloor >= Nland * (Wpay + Wseat),
                Afloor >= 2. * Mfloor / (sigfloor * hfloor) + 1.5 * Sfloor / taufloor,
                Vfloor == 2 * wfloor * Afloor,
                Wfloor >= rhofloor * g * Vfloor + 2 * wfloor * lfloor * Wppfloor,
                # hfloor <= 0.1 * Rfuse,

                # Tail cone sizing
                taucone == sigskin,
                Wcone >= rhocone * g * Vcone * (1 + fstring + fframe),
                xtail >= lnose + lshell + .5 * lcone,
                
                # BENDING MODEL
                # Maximum axial stress is the sum of bending and pressurization
                # stresses
                Ihshell <= ((pi + 4 * thetadb) * Rfuse**2 + 8.*(1-thetadb**2/2) * (dRfuse/2.)*Rfuse + \
                            (2*pi + 4 * thetadb)*(dRfuse/2)**2) * Rfuse * tshell + \
                    2 / 3 * (hdb + dRfuse/2.)**3 * tdb,  # [SP]
                Ivshell <= (pi*Rfuse**2 + 8*wdb*Rfuse + (2*pi+4*thetadb)*wdb**2)*Rfuse*tshell, #[SP] #Ivshell
                # approximation needs to be improved

                # Horizontal bending material model
                # Calculating xhbend, the location where additional bending
                # material is required
                xhbendLand >= xwing, xhbendLand <= lfuse,
                xhbendMLF >= xwing, xhbendMLF <= lfuse,

                SignomialEquality(A0h, A2hLand * (xshell2 - xhbendLand) ** 2 + A1hLand * (xtail - xhbendLand)), # [SP] #[SPEquality]
                SignomialEquality(A0h, A2hMLF * (xshell2 - xhbendMLF) ** 2 + A1hMLF * (xtail - xhbendMLF)), # [SP] #[SPEquality]

                A2hLand >= Nland * (Wpay + Wpadd + Wshell + Wwindow + Winsul + Wfloor + Wseat) / \
                (2 * lshell * hfuse * sigbend),  # Landing loads constant A2hLand
                A2hMLF >= Nlift * (Wpay + Wpadd + Wshell + Wwindow + Winsul + Wfloor + Wseat) / \
                (2 * lshell * hfuse * sigMh),  # Max wing aero loads constant A2hMLF

                # Shell inertia constant A0h
                A0h == (Ihshell / (rE * hfuse**2)), # [SP]

                # Bending area behind wingbox
                AhbendfLand >= A2hLand * (xshell2 - xf)**2 + A1hLand * (xtail - xf) - A0h, # [SP]
                AhbendfMLF >= A2hMLF * (xshell2 - xf)**2 + A1hMLF * (xtail - xf) - A0h, # [SP]

                # Bending area in front of wingbox
                AhbendbLand >= A2hLand * (xshell2 - xb)**2 + A1hLand * (xtail - xb) - A0h, # [SP]
                AhbendbMLF >= A2hMLF * (xshell2 - xb)**2 + A1hMLF * (xtail - xb) - A0h, # [SP]

                # Bending volume forward of wingbox
                Vhbendf >= A2hLand / 3 * ((xshell2 - xf)**3 - (xshell2 - xhbendLand)**3) \
                + A1hLand / 2 * ((xtail - xf)**2 - (xtail - xhbendLand)**2) \
                - A0h * (xhbendLand - xf),  # [SP]
                Vhbendf >= A2hMLF / 3 * ((xshell2 - xf)**3 - (xshell2 - xhbendMLF)**3) \
                + A1hMLF / 2 * ((xtail - xf)**2 - (xtail - xhbendMLF)**2) \
                - A0h * (xhbendMLF - xf),  # [SP]

                # Bending volume behind wingbox
                Vhbendb >= A2hLand / 3 * ((xshell2 - xb)**3 - (xshell2 - xhbendLand)**3) \
                + A1hLand / 2 * ((xtail - xb)**2 - (xtail - xhbendLand)**2) \
                - A0h * (xhbendLand - xb),  # [SP]
                Vhbendb >= A2hMLF / 3 * ((xshell2 - xb)**3 - (xshell2 - xhbendMLF)**3) \
                + A1hMLF / 2 * ((xtail - xb)**2 - (xtail - xhbendMLF)**2) \
                - A0h * (xhbendMLF - xb),  # [SP]

                # Bending volume over wingbox
                Vhbendc >= .5 * (AhbendfLand + AhbendbLand) * c0 * w,
                Vhbendc >= .5 * (AhbendfMLF + AhbendbMLF) * c0 * w,

                # Determining more constraining load case (landing vs. max aero horizontal bending)
                Vhbend >= Vhbendc + Vhbendf + Vhbendb,
                Whbend >= g * rhobend * Vhbend,

                # Vertical bending material model
                # Calculating xvbend, the location where additional bending material is required
                xvbend >= xwing, xvbend <= lfuse,
                SignomialEquality(B0v, B1v * (xtail - xvbend)), # [SP] #[SPEquality]
                #B1v definition in Aircraft()
                B0v == Ivshell/(rE*wfuse**2),
                Avbendb >= B1v * (xtail - xb) - B0v,
                Vvbendb >= 0.5*B1v * ((xtail-xb)**2 - (xtail - xvbend)**2) - B0v * (xvbend - xb),
                Vvbendc >= 0.5*Avbendb*c0*w,
                Vvbend >= Vvbendb + Vvbendc,
                Wvbend >= rhobend*g*Vvbend,

                # Wing variable substitutions
                SignomialEquality(xf, xwing + .5 * c0 * w),  # [SP] [SPEquality]
                SignomialEquality(xb, xwing - .5 * c0 * w),  # [SP] [SPEquality]

                sigMh <= sigbend - rE * dPover / 2 * Rfuse / tshell,
                sigMv <= sigbend - rE * dPover / 2 * Rfuse / tshell,

                # Volume relations
                Vcyl == Askin * lshell,
                Vnose == Snose * tskin,
                Vbulk == Sbulk * tskin,
                Vdb == Adb * lshell,
                # [SP] #[SPEquality]
                Vcabin >= Afuse * (lshell + 0.67 * lnose + 0.67 * Rfuse),

                # Weight relations
                Wapu == Wpay * fapu,
                Wdb == rhoskin * g * Vdb,
                Winsul >= Wppinsul * ((1.1 * pi + 2 * thetadb) * Rfuse * lshell + 0.55 * (Snose + Sbulk)),
                Wlugg >= flugg2 * npax * 2 * Wchecked + flugg1 * npax * Wchecked + Wcarryon,
                Wwindow >= Wpwindow * lshell,
                Wpadd == Wpay * fpadd,
                Wseat >= Wpseat * nseat,
                Wseat >= fseat * Wpay,

                Wskin >= rhoskin * g * (Vcyl + Vnose + Vbulk),
                Wshell >= Wskin * (1 + fstring + ffadd + fframe) + Wdb,
                Wfuse >= Wshell + Wfloor + Winsul + \
                    Wapu + Wfix + Wwindow + Wpadd + Wseat + Whbend + Wvbend,
            ])

        return constraints

    def dynamic(self, state):
        """
        returns a fuselage performance model
        """
        return FuselagePerformance(self, state)


class FuselagePerformance(Model):
    """
    Fuselage performance model
    """

    def setup(self, fuse, state, **kwargs):
        # new variables
        CDfuse = Variable('C_{D_{fuse}}', '-', 'Fuselage Drag Coefficient')
        Dfuse = Variable('D_{fuse}', 'N', 'Fuselage Drag')
        # CLfuse = Variable('C_{L_{fuse}}','-', 'Fuselage Lift Coefficient')
        Lfuse = Variable('L_{fuse}','N','Fuselage Lift')
        # Dfrict = Variable('D_{friction}', 'N', 'Friction drag')
        # Dupswp = Variable('D_{upsweep}', 'N', 'Drag due to fuse upsweep')
        # f = Variable('f', '-', 'Fineness ratio')
        # FF = Variable('FF', '-', 'Fuselage form factor')
        # phi = Variable('\\phi', '-', 'Upsweep angle')


        # BLI surrogate
        fBLI = Variable('f_{BLI}','-','1-fBLI surrogate')

        constraints = []
        constraints.extend([
            CDfuse == CDfuse,
            Dfuse == Dfuse,
            Lfuse == Lfuse,
            #Dfuse == Cdfuse * (.5 * fuse['A_{fuse}'] * state.atm['\\rho'] * state['V']**2),
            # fineness ratio
            # f == fuse['l_{fuse}'] / ((4 / np.pi * fuse['A_{fuse}'])**0.5),
            # FF >= 1 + 60 / f**3 + f / 400,  # form factor
            # Dfrict >= FF * np.pi * fuse['R_{fuse}'] * state.atm['\\mu'] * state['V'] * 0.074 * (state.atm['\\rho'] * state['V']
            #                                                                                     * fuse['l_{fuse}'] / state.atm['\\mu'])**0.8,
            # Monomial fit of tan(phi)
            # 1.13226 * phi**1.03759 == fuse['R_{fuse}'] / fuse['l_{cone}'],
            # Dupswp >= 3.83 * phi**2.5 * fuse['A_{fuse}'] * 0.5 * state.atm['\\rho'] * state['V']**2,
            # Dfuse >= Dfrict + Dupswp,
        ])

        return constraints
