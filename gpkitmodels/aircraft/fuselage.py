"Implements Fuselage model"
import numpy as np
from gpkit import Variable, Model, SignomialsEnabled, units
from gpkit.constraints.costed import CostedConstraintSet
from gpkit.constraints.tight import TightConstraintSet as TCS

class Fuselage(CostedConstraintSet):
    """
    Fuselage model
    """
    def __init__(self):

        Afloor   = Variable('A_{floor}', 'm^2', 'Floor beam x-sectional area')
        Afuse    = Variable('A_{fuse}', 'm^2', 'Fuselage x-sectional area')
        Ahold    = Variable('A_{hold}', 'm^2', 'Cargo hold x-sectional area')
        Askin    = Variable('A_{skin}', 'm^2', 'Skin cross sectional area')
        Dfuse    = Variable('D_{fuse}', 'N', 'Total drag in cruise')
        Dfrict   = Variable('D_{friction}', 'N', 'Friction drag')
        Dupswp   = Variable('D_{upsweep}', 'N', 'Drag due to fuse upsweep')
        FF       = Variable('FF', '-','Fuselage form factor')
        LF       = Variable('LF', '-', 'Load factor')
        Lvmax    = Variable('L_{v_{max}}', 'N', 'Max vertical tail load')
        Mfloor   = Variable('M_{floor}', 'N*m',
                            'Max bending moment in floor beams')
        Nland    = Variable('N_{land}', '-', 'Emergency landing load factor')
        Pfloor   = Variable('P_{floor}', 'N', 'Distributed floor load')
        Qv       = Variable('Q_v', 'N*m', 'Torsion moment imparted by tail')
        R        = Variable('R', 287, 'J/(kg*K)', 'Universal gas constant')
        Rfuse    = Variable('R_{fuse}', 'm', 'Fuselage radius')
        SPR      = Variable('SPR', 6, '-', 'Number of seats per row')
        Sbulk    = Variable('S_{bulk}', 'm^2', 'Bulkhead surface area')
        Sfloor   = Variable('S_{floor}', 'N', 'Maximum shear in floor beams')
        Snose    = Variable('S_{nose}', 'm^2', 'Nose surface area')
        Tcabin   = Variable('T_{cabin}', 300, 'K', 'Cabin temperature')
        Vbulk    = Variable('V_{bulk}', 'm^3', 'Bulkhead skin volume')
        Vcabin   = Variable('V_{cabin}', 'm^3', 'Cabin volume')
        Vcargo   = Variable('V_{cargo}', 'm^3', 'Cargo volume')
        Vcone    = Variable('V_{cone}', 'm^3', 'Cone skin volume')
        Vcyl     = Variable('V_{cyl}', 'm^3', 'Cylinder skin volume')
        Vfloor   = Variable('V_{floor}', 'm^3', 'Floor volume')
        Vhold    = Variable('V_{hold}', 'm^3', 'Hold volume')
        Vinf     = Variable('V_{\\infty}', 'm/s', 'Cruise velocity')
        Vlugg    = Variable('V_{lugg}', 'm^3', 'Luggage volume')
        Vnose    = Variable('V_{nose}', 'm^3', 'Nose skin volume')
        Wapu     = Variable('W_{apu}', 'N', 'APU weight')
        Wavgpass = Variable('W_{avg. pass}', 180, 'lbf',
                            'Average passenger weight')
        Wbuoy    = Variable('W_{buoy}', 'N', 'Buoyancy weight')
        Wcargo   = Variable('W_{cargo}', 10000, 'N', 'Cargo weight')
        Wcarryon = Variable('W_{carry on}', 15, 'lbf', 'Ave. carry-on weight')
        Wchecked = Variable('W_{checked}', 40, 'lbf', 'Ave. checked bag weight')
        Wcone    = Variable('W_{cone}', 'N', 'Cone weight')
        Wfix     = Variable('W_{fix}', 3000, 'lbf',
                            'Fixed weights (pilots, cockpit seats, navcom)')
        Wfloor   = Variable('W_{floor}', 'N', 'Floor weight')
        Wfuse    = Variable('W_{fuse}', 'N', 'Fuselage weight')
        Winsul   = Variable('W_{insul}', 'N', 'Insulation material weight')
        Wlugg    = Variable('W_{lugg}', 'N', 'Passenger luggage weight')
        Wpadd    = Variable('W_{padd}', 'N',
                            'Misc weights (galley, toilets, doors etc.)')
        Wpass    = Variable('W_{pass}', 'N', 'Passenger weight')
        Wpay     = Variable('W_{pay}', 'N', 'Payload weight')
        Wppfloor = Variable('W\'\'_{floor}', 60, 'N/m^2',
                            'Floor weight/area density')
        Wppinsul = Variable('W\'\'_{insul}', 22, 'N/m^2',
                            'Weight/area density of insulation material')
        Wpseat   = Variable('W\'_{seat}', 150, 'N', 'Weight per seat')
        Wpwindow = Variable('W\'_{window}', 145.*3, 'N/m',
                            'Weight/length density of windows')
        Wseat    = Variable('W_{seat}', 'N', 'Seating weight')
        Wshell   = Variable('W_{shell}', 'N', 'Shell weight')
        Wskin    = Variable('W_{skin}', 'N', 'Skin weight')
        Wwindow  = Variable('W_{window}', 'N', 'Window weight')
        bvt      = Variable('b_{vt}', 'm', 'Vertical tail span')
        cvt      = Variable('c_{vt}', 'm', 'Vertical tail root chord')
        dh       = Variable('\\Delta h', 'm',
                            'Distance from floor to widest part of fuselage')
        dp       = Variable('\\Delta p', 'Pa',
                            'Pressure difference across fuselage skin')
        f        = Variable('f', '-', 'Fineness ratio')
        fapu     = Variable('f_{apu}', '-',
                            'APU weight as fraction of payload weight')
        ffadd    = Variable('f_{fadd}', '-',
                            'Fractional added weight of local reinforcements')
        fframe   = Variable('f_{frame}', '-', 'Fractional frame weight')
        flugg1   = Variable('f_{lugg,1}', '-',
                            'Proportion of passengers with one suitcase')
        flugg2   = Variable('f_{lugg,2}', '-',
                            'Proportion of passengers with two suitcases')
        fpadd    = Variable('f_{padd}', '-',
                            'Other misc weight as fraction of payload weight')
        fseat    = Variable('f_{seat}', '-',
                            'Seat weight as fraction of payload weight')
        fstring  = Variable('f_{string}', '-','Fractional weight of stringers')
        g        = Variable('g', 9.81, 'm/s^2', 'Gravitational acceleration')
        hfloor   = Variable('h_{floor}', 'm', 'Floor I-beam height')
        hhold    = Variable('h_{hold}', 'm', 'Height of the cargo hold')
        lamcone  = Variable('\\lambda_{cone}', '-',
                            'Tailcone radius taper ratio (xshell2->xtail)')
        lcone    = Variable('l_{cone}', 'm', 'Cone length')
        lfloor   = Variable('l_{floor}', 'm', 'Floor length')
        lfuse    = Variable('l_{fuse}', 'm', 'Fuselage length')
        lnose    = Variable('l_{nose}', 'm', 'Nose length')
        lshell   = Variable('l_{shell}', 'm', 'Shell length')
        mu       = Variable('\\mu', 1.4E-5, 'N*s/m^2',
                            'Dynamic viscosity (35,000 ft)')
        npass    = Variable('n_{pass}', '-', 'Number of passengers')
        nrows    = Variable('n_{rows}', '-', 'Number of rows')
        nseat    = Variable('n_{seat}', '-',' Number of seats')
        pcabin   = Variable('p_{cabin}', 'Pa','Cabin air pressure (8,000ft)')
        phi      = Variable('\\phi', '-', 'Upsweep angle')
        pitch    = Variable('p_s', 'in', 'Seat pitch')
        plamv    = Variable('p_{\\lambda_v}', '-', '1 + 2*Tail taper ratio')
        rE       = Variable('r_E', '-', 'Ratio of stringer/skin moduli')
        rhobend  = Variable('\\rho_{bend}', 2700, 'kg/m^3', 'Stringer density')
        rhocabin = Variable('\\rho_{cabin}', 'kg/m^3', 'Air density in cabin')
        rhocargo = Variable('\\rho_{cargo}', 150, 'kg/m^3', 'Cargo density')
        rhocone  = Variable('\\rho_{cone}', 2700, 'kg/m^3',
                            'Cone material density')
        rhofloor = Variable('\\rho_{floor}', 2700, 'kg/m^3',
                            'Floor material density')
        rhoinf   = Variable('\\rho_{\\infty}', 0.38, 'kg/m^3',
                            'Air density (35,000ft)')
        rholugg  = Variable('\\rho_{lugg}', 100, 'kg/m^3', 'Luggage density')
        rhoskin  = Variable('\\rho_{skin}', 2700, 'kg/m^3', 'Skin density')
        sigfloor = Variable('\\sigma_{floor}', 30000/0.000145, 'N/m^2',
                            'Max allowable cap stress')
        sigskin  = Variable('\\sigma_{skin}', 15000/0.000145, 'N/m^2',
                            'Max allowable skin stress')
        sigth    = Variable('\\sigma_{\\theta}', 'N/m^2', 'Skin hoop stress')
        sigx     = Variable('\\sigma_x', 'N/m^2', 'Axial stress in skin')
        taucone  = Variable('\\tau_{cone}', 'N/m^2', 'Shear stress in cone')
        taufloor = Variable('\\tau_{floor}', 30000/0.000145, 'N/m^2',
                            'Max allowable shear web stress')
        tcone    = Variable('t_{cone}', 'm', 'Cone thickness')
        tshell   = Variable('t_{shell}', 'm', 'Shell thickness')
        tskin    = Variable('t_{skin}', 'm', 'Skin thickness')
        waisle   = Variable('w_{aisle}', 0.51, 'm', 'Aisle width')
        wfloor   = Variable('w_{floor}', 'm', 'Floor width')
        wfuse    = Variable('w_{fuse}', 'm', 'Fuselage width')
        wseat    = Variable('w_{seat}', 0.5, 'm', 'Seat width')
        wsys     = Variable('w_{sys}', 0.10, 'm',
                            'Width between cabin and skin for systems')
        xCGfu    = Variable('x_{CG_{fu}}', 'm', 'x-location of fuselage CG')
        xVbulk   = Variable('xVbulk', 'm^4', 'Volume moment of bulkhead')
        xVcyl    = Variable('xVcyl', 'm^4', 'Volume moment of cylinder')
        xVnose   = Variable('xVnose', 'm^4', 'Volume moment of nose')
        xWapu    = Variable('xWapu', 'N*m', 'Moment of APU')
        xWcone   = Variable('xWcone', 'N*m', 'Moment of cone')
        xWfix    = Variable('xWfix', 'N*m', 'Moment of fixed weights')
        xWfloor  = Variable('xWfloor', 'N*m', 'Moment of floor weight')
        xWfuse   = Variable('xWfuse', 'N*m', 'Fuselage moment')
        xWinsul  = Variable('xWinsul', 'N*m', 'Moment of insulation material')
        xWpadd   = Variable('xWpadd', 'N*m', 'Moment of misc weights')
        xWseat   = Variable('xWseat', 'N*m', 'Moment of seats')
        xWshell  = Variable('xWshell', 'N*m', 'Mass moment of shell')
        xWskin   = Variable('xWskin', 'N*m', 'Mass moment of skin')
        xWwindow = Variable('xWwindow', 'N*m', 'Mass moment of windows')
        x_upswp  = Variable('x_{up}', 'm', 'Fuselage upsweep point')
        xapu     = Variable('xapu', 120, 'ft', 'x-location of APU')
        xconend  = Variable('xconend', 'm', 'x-location of cone end')
        xfix     = Variable('xfix', 2.1, 'm', 'x-location of fixed weight')
        xshell1  = Variable('x_{shell1}', 'm', 'Start of cylinder section')
        xshell2  = Variable('x_{shell2}', 'm', 'End of cylinder section')

        with SignomialsEnabled():
            objective = Dfuse + 0.5*Wfuse

            constraints = [
                            # Geometry relations
                            lnose == xshell1,

                            # Cross section relations
                            wfuse == 2*Rfuse,
                            wfuse >= SPR*wseat + waisle + 2*wsys,
                            Askin >= 2*np.pi*Rfuse*tskin, # simplified
                            Afuse >= np.pi*Rfuse**2, # simplified
                            tshell >= tskin*(1 + rE*fstring*rhoskin/rhobend),

                            # TODO: Bending loads

                            # Pressure shell loads
                            sigx == dp/2*Rfuse/tshell,
                            sigth == dp*Rfuse/tskin,
                            sigskin >= sigth,
                            sigskin >= sigx,
                            Snose**(8./5) >= (2*np.pi*Rfuse**2)**(8./5) *
                                             (1./3 + (2./3)*(lnose/Rfuse)
                                             **(8./5)),
                            Sbulk == 2*np.pi*Rfuse**2,
                            Vcyl == Askin*lshell,
                            Vnose == Snose*tskin,
                            Vbulk == Sbulk*tskin,
                            Wskin >= rhoskin*g*(Vcyl + Vnose + Vbulk),
                            Wshell >= Wskin*(1 + fstring + fframe + ffadd),

                            # Cabin volume and buoyancy weight
                            rhocabin == (1/(R*Tcabin))*pcabin,
                            Vcabin >= Afuse*(lshell + 0.67*lnose + 0.67*Rfuse),
                            TCS([Wbuoy >= (rhocabin - rhoinf)*g*Vcabin], reltol=1E-3), # [SP]

                            # Windows and insulation
                            Wwindow == Wpwindow * lshell,
                            Winsul >= Wppinsul*(1.1*np.pi*Rfuse*lshell
                                      + 0.55*(Snose + Sbulk)),

                            # Payload-proportional weights
                            Wapu == Wpay*fapu,
                            Wseat == Wpseat*nseat,
                            Wpadd == Wpay*fpadd,

                            # Floor
                            Pfloor >= Nland*(Wpay + Wseat),
                            Sfloor == 0.5*Pfloor, # without support
                            Mfloor == 0.25*Pfloor*wfloor/2,
                            lnose >= 5.2*units.m, # TODO less arbitrary
                            Afloor >= 2*Mfloor/(sigfloor*hfloor)
                                      + 1.5*Sfloor/taufloor,
                            Vfloor >= wfloor*Afloor,
                            lfloor >= lshell + 2*Rfuse,
                            lshell >= nrows*pitch,
                            (wfloor/2)**2 + dh**2 >= Rfuse**2, # [SP]
                            Wfloor >= rhofloor*g*Vfloor
                                      + wfloor*lfloor*Wppfloor,

                            # Tail cone
                            Rfuse*taucone*(1+plamv)*Vcone*(1+lamcone)/(4*lcone)
                                >= Lvmax*bvt*plamv/3, # [SP]
                            plamv >= 1.6,
                            taucone == sigskin,
                            lamcone == 0.4, # TODO remove
                            lamcone == cvt/lcone,
                            Wcone >= rhocone*g*Vcone*(1 + fstring + fframe
                                                      + ffadd),

                            # Payload weight breakdown
                            npass == nseat*LF,
                            nseat == nrows*SPR,
                            Wpass == npass*Wavgpass,
                            Wlugg >= flugg2*npass*2*Wchecked
                                     + flugg1*npass*Wchecked + Wcarryon,
                            Wlugg == Vlugg*g*rholugg,
                            Wcargo == Vcargo*g*rhocargo,
                            Vhold >= Vcargo + Vlugg,
                            Vhold == Ahold*lshell,
                            TCS([Ahold <= (2./3)*wfloor*hhold + hhold**3/(2*wfloor)], reltol=1E-5),
                            # [SP] Harris stocker 1998 (wolfram)
                            TCS([dh + hhold + hfloor <= Rfuse]),
                            Wpay >= Wpass + Wlugg + Wcargo,

                            # Total fuselage weight
                            Wfuse >= Wfix + Wapu + Wpadd + Wseat + Wshell
                                   + Wwindow + Winsul + Wcone + Wfloor + Wbuoy,

                            # Drag
                            # sources: Raymer (p285), kfid 325 notes (p180)
                            lfuse >= lnose+lshell+lcone,
                            f == lfuse/((4/np.pi*Afuse)**0.5), # fineness ratio
                            FF >= 1 + 60/f**3 + f/400, # form factor
                            Dfrict >= FF * np.pi*Rfuse * mu*Vinf
                                      * 0.074*(rhoinf*Vinf*lfuse/mu)**0.8,

                            # Drag due to fuselage upsweep (Raymer p286)
                            xshell2 >= xshell1 + lshell,
                            xshell2 == x_upswp,
                            x_upswp + lcone <= lfuse,
                            1.13226*phi**1.03759 == Rfuse/lcone, # monomial fit
                                                                 # of tan(phi)
                            Dupswp >= 3.83*phi**2.5*Afuse * 0.5*rhoinf*Vinf**2,
                            Dfuse >= Dfrict + Dupswp
                          ]

        CG_constraints = [
                          xVcyl >= 0.5*(xshell1+xshell2)*Vcyl,
                          xVnose >= 0.5*(xshell1)*Vnose,
                          xVbulk >= xshell2*Vbulk, # simplified
                          xWskin >= rhoskin*g*(xVcyl + xVnose + xVbulk),
                          xWshell >= xWskin*(1 + fstring + fframe + ffadd),
                          xWwindow >= 0.5* (xshell1+xshell2)*Wwindow,
                          xWinsul >= 0.5*(xshell1 + xshell2)*Winsul,
                          xWapu == xapu*Wapu,
                          xWseat >= 0.5*(xshell1 + xshell2)*Wseat,
                          xWpadd >= 0.5*(xshell1 + xshell2)*Wpadd,
                          xWfix == xfix*Wfix,
                          xWfloor >= 0.5*(xshell1 + xshell2)*Wfloor,
                          xWcone >= 0.5*(xshell2 + lfuse) * Wcone,
                          xWfuse >= xWfix + xWapu + xWpadd + xWseat + xWshell
                                  + xWcone + xWwindow + xWinsul + xWfloor,
                          xCGfu == xWfuse/Wfuse,
                         ]
        self.CG_constraints = CG_constraints

        CostedConstraintSet.__init__(self, objective, constraints)


    def default737subs(self):

        substitutions = {
                         'LF': 0.898,
                         'L_{v_{max}}': 35000,
                         'N_{land}': 6.0, # [TAS]
                         'SPR': 6,
                         'T_{cabin}': 300,
                         'V_{\\infty}': 234,
                         'W\'\'_{floor}': 60, # [TAS]
                         'W\'\'_{insul}': 22, # [TAS]
                         'W\'_{seat}': 150, # Boeing
                         'W\'_{window}': 145.*3, # [TAS]
                         'W_{avg. pass}': 180,
                         'W_{cargo}': 10000,
                         'W_{carry on}': 15,
                         'W_{checked}': 40,
                         'W_{fix}': 3000,
                         '\\Delta h': 1,
                         '\\Delta p': 52000,
                         '\\rho_{\\infty}': 0.38,
                         '\\rho_{bend}': 2700, # [TAS]
                         '\\rho_{cargo}': 150, # b757 freight doc
                         '\\rho_{cone}': 2700, # [TAS]
                         '\\rho_{floor}': 2700, # [TAS]
                         '\\rho_{lugg}': 100,
                         '\\rho_{skin}': 2700, # [TAS]
                         '\\sigma_{floor}': 30000/0.000145, # [TAS]
                         '\\sigma_{skin}': 15000/0.000145, # [TAS]
                         '\\tau_{floor}': 30000/0.000145, # [TAS]
                         'b_{vt}': 7,
                         'c_{vt}': 4,
                         'f_{apu}': 0.035, # [TAS]
                         'f_{fadd}': 0.20, # [TAS]
                         'f_{frame}': 0.25,
                         'f_{lugg,1}': 0.4,
                         'f_{lugg,2}': 0.1,
                         'f_{padd}': 0.4, # [TAS]
                         'f_{seat}': 0.10,
                         'f_{string}': 0.35, # [TAS]
                         'n_{seat}': 186,
                         'p_s': 31,
                         'p_{cabin}': 75000,
                         'r_E': 1.0, # [TAS]
                         'w_{aisle}': 0.51, # Boeing
                         'w_{seat}': 0.5,
                         'w_{sys}': 0.10,
                         'xapu': 120,
                         'xfix': 2.1
                        }

        return substitutions

    @classmethod
    def standalone737(cls):
        """Create standalone instance of fuselage model"""

        ccs = cls()

        substitutions = ccs.default737subs()

        m = Model(ccs.cost, ccs, substitutions)
        return m

    @classmethod
    def coupled737(cls):
        """Creates instance of fuselage model for use in full aircraft model"""

        ccs = cls()

        constraints = ccs + ccs.CG_constraints
 
        dsubs = ccs.default737subs()
        linkedsubs = ['L_{v_{max}}', 'V_{\\infty}', 'b_{vt}', 'c_{vt}']
        substitutions = {key: value for key, value in dsubs.items()
                                    if key not in linkedsubs}

        m = Model(ccs.cost, constraints, substitutions, name='Fuselage')
        return m

    @classmethod
    def test(cls):
        fu = cls.standalone737()
        sol = fu.localsolve()


if __name__ == "__main__":
    Fuselage.test()
