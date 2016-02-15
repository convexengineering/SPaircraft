from gpkit import Variable, Model, SignomialsEnabled, units
import numpy as np
import numpy.testing as npt

class Fuselage(Model):
    """
    Fuselage model
    """
    def setup(self):

        # Variables
        Afloor    = Variable('A_{floor}', 'm^2', 'Floor beam x-sectional area')
        Afuse     = Variable('A_{fuse}', 'm^2', 'Fuselage x-sectional area')
        Ahold     = Variable('A_{hold}', 'm^2', 'Cargo hold x-sectional area')
        Askin     = Variable('A_{skin}', 'm^2', 'Skin cross sectional area')
        D         = Variable('D', 'N', 'Total drag in cruise')
        Dfrict    = Variable('D_{friction}', 'N', 'Friction drag')
        Dupswp    = Variable('D_{upsweep}', 'N', 'Drag due to fuse upsweep')
#       dRfuse    = Variable('\\Delta R_{fuse}', 'm', 'Lower cylinder offset')
        f         = Variable('f', '-', 'Fineness ratio')
        FF        = Variable('FF', '-','Fuselage form factor')
        hfloor    = Variable('h_{floor}', 'm', 'Floor I-beam height')
        hhold     = Variable('h_{hold}', 'm', 'Height of the cargo hold')
#       Ihshell   = Variable('Ihshell', 'm^4', 'Bending inertia')
#       Ivshell   = Variable('Ivshell', 'm^4', 'Bending inertia')
        lfloor    = Variable('l_{floor}', 'm', 'Floor length')
        lfuse     = Variable('l_{fuse}', 'm', 'Fuselage length')
        lcone     = Variable('l_{cone}', 'm', 'Cone length')
        lnose     = Variable('l_{nose}', 'm', 'Nose length')
        lshell    = Variable('l_{shell}', 'm', 'Shell length')
        Mfloor    = Variable('M_{floor}', 'N*m',
                             'Max bending moment in floor beams')
        nrows     = Variable('n_{rows}', '-', 'Number of rows')
        npass     = Variable('n_{pass}', '-', 'Number of passengers')
        Pfloor    = Variable('P_{floor}', 'N', 'Distributed floor load')
        phi       = Variable('\\phi', '-', 'Upsweep angle')
        lamcone   = Variable('\\lambda_{cone}', '-',
                             'Tailcone radius taper ratio (xshell2->xtail)')
        plamv     = Variable('p_{\\lambda_v}', '-', '1 + 2*Tail taper ratio')
        Qv        = Variable('Q_v', 'N*m', 'Torsion moment imparted by tail')
        Rfuse     = Variable('R_{fuse}', 'm', 'Fuselage radius')
        rhocabin  = Variable('\\rho_{cabin}', 'kg/m^3', 'Air density in cabin')
        Sbulk     = Variable('S_{bulk}', 'm^2', 'Bulkhead surface area')
        Sfloor    = Variable('S_{floor}', 'N', 'Maximum shear in floor beams')
        sigth     = Variable('\\sigma_{\\theta}', 'N/m^2', 'Skin hoop stress')
        sigx      = Variable('\\sigma_x', 'N/m^2', 'Axial stress in skin')
        Snose     = Variable('S_{nose}', 'm^2', 'Nose surface area')
        taucone   = Variable('\\tau_{cone}', 'N/m^2', 'Shear stress in cone')
        tcone     = Variable('t_{cone}', 'm', 'Cone thickness')
        tshell    = Variable('t_{shell}', 'm', 'Shell thickness')
        tskin     = Variable('t_{skin}', 'm', 'Skin thickness')
        Vbulk     = Variable('V_{bulk}', 'm^3', 'Bulkhead skin volume')
        Vcabin    = Variable('V_{cabin}', 'm^3', 'Cabin volume')
        Vcargo    = Variable('V_{cargo}', 'm^3', 'Cargo volume')
        Vcone     = Variable('V_{cone}', 'm^3', 'Cone skin volume')
        Vcyl      = Variable('V_{cyl}', 'm^3', 'Cylinder skin volume')
        Vfloor    = Variable('V_{floor}', 'm^3', 'Floor volume')
        Vhold     = Variable('V_{hold}', 'm^3', 'Hold volume')
        Vlugg     = Variable('V_{lugg}', 'm^3', 'Luggage volume')
        Vnose     = Variable('V_{nose}', 'm^3', 'Nose skin volume')
        Wapu      = Variable('W_{apu}', 'N', 'APU weight')
        Wbuoy     = Variable('W_{buoy}', 'N', 'Buoyancy weight')
        Wcone     = Variable('W_{cone}', 'N', 'Cone weight')
        wfloor    = Variable('w_{floor}', 'm', 'Floor width')
        Wfloor    = Variable('W_{floor}', 'N', 'Floor weight')
        Wfuse     = Variable('W_{fuse}', 'N', 'Fuselage weight')
        Winsul    = Variable('W_{insul}', 'N', 'Insulation material weight')
        Wlugg     = Variable('W_{lugg}', 'N', 'Passenger luggage weight')
        Wpadd     = Variable('W_{padd}', 'N',
                             'Misc weights (galley, toilets, doors etc.)')
        Wpass     = Variable('W_{pass}', 'N', 'Passenger weight')
        Wpay      = Variable('W_{pay}', 'N', 'Payload weight')
        Wseat     = Variable('W_{seat}', 'N', 'Seating weight')
        Wshell    = Variable('W_{shell}', 'N', 'Shell weight')
        Wskin     = Variable('W_{skin}', 'N', 'Skin weight')
        Wwindow   = Variable('W_{window}', 'N', 'Window weight')
#       xapu      = Variable('xapu', 120*0.3048, 'm', 'x-location of APU')
#       xconend   = Variable('xconend', 'm', 'x-location of cone end')
        xshell1   = Variable('x_{shell1}', 'm', 'Start of cylinder section')
#       xshell2   = Variable('xshell2', 102*0.3048, 'm', 'End of cylinder sect')
#       xVbulk    = Variable('xVbulk', 'm^4', 'Volume moment of bulkhead')
#       xVcyl     = Variable('xVcyl', 'm^4', 'Volume moment of cylinder')
#       xVnose    = Variable('xVnose', 'm^4', 'Volume moment of nose')
#       xWapu     = Variable('xWapu', 'N*m', 'Moment of APU')
#       xWcone    = Variable('xWcone', 'N*m', 'Moment of cone')
#       xWfix     = Variable('xWfix', 'N*m', 'Moment of fixed weights')
#       xWfloor   = Variable('xWfloor', 'N*m', 'Moment of floor weight')
#       xWfuse    = Variable('xWfuse', 'N*m', 'Fuselage moment')
#       xWinsul   = Variable('xWinsul', 'N*m', 'Moment of insulation material')
#       xWpadd    = Variable('xWpadd', 'N*m', 'Moment of misc weights')
#       xWseat    = Variable('xWseat', 'N*m', 'Moment of seats')
#       xWshell   = Variable('xWshell', 'N*m', 'Mass moment of shell')
#       xWskin    = Variable('xWskin', 'N*m', 'Mass moment of skin')
#       xWwindow  = Variable('xWwindow', 'N*m', 'Mass moment of windows')

        # Constants ([TAS] means sourced from TASOPT)
        bv       = Variable('b_v', 7, 'm', 'Vertical tail span')
        cvt      = Variable('c_{vt}', 4, 'm', 'Vertical tail chord')
        dh       = Variable('\\Delta h', 1, 'm',
                            'Distance from floor to widest part of fuselage')
        dp       = Variable('\\Delta p', 52000, 'Pa',
                            'Pressure difference across fuselage skin')
        fapu     = Variable('f_{apu}', 0.035, '-',
                            'APU weight as fraction of payload weight') # [TAS]
        ffadd    = Variable('f_{fadd}', 0.20, '-',
                            'Fractional added weight of local reinforcements') # [TAS]
        fframe   = Variable('f_{frame}', 0.25, '-', 'Fractional frame weight')
        flugg1   = Variable('f_{lugg,1}', 0.4, '-',
                            'Proportion of passengers with one suitcase')
        flugg2   = Variable('f_{lugg,2}', 0.1, '-',
                            'Proportion of passengers with two suitcases')
        fpadd    = Variable('f_{padd}', 0.4, '-',
                            'Other misc weight as fraction of payload weight') # [TAS]
        fseat    = Variable('f_{seat}', 0.10, '-',
                            'Seat weight as fraction of payload weight')
        fstring  = Variable('f_{string}', 0.35, '-',
                            'Fractional weight of stringers') # [TAS]
        g        = Variable('g', 9.81, 'm/s^2', 'Gravitational acceleration')
        LF       = Variable('LF', 0.898, '-', 'Load factor')
        Lvmax    = Variable('L_{v_{max}}', 35000, 'N', 'Max vertical tail load')
        mu       = Variable('\\mu', 1.4E-5, 'N*s/m^2',
                            'Dynamic viscosity (35,000ft)')
        Nland    = Variable('N_{land}', 6.0,  '-',
                            'Emergency landing load factor') # [TAS]
        nseat    = Variable('n_{seat}', 186, '-',' Number of seats')
        pcabin   = Variable('p_{cabin}', 75000, 'Pa',
                            'Cabin air pressure (8,000ft)')
        pitch    = Variable('p_s', 31 , 'in', 'Seat pitch') # Boeing doc
        R        = Variable('R', 287, 'J/(kg*K)', 'Universal gas constant')
        rE       = Variable('r_E', 1.0, '-', 'Ratio of stringer/skin moduli') # [TAS]
        rhobend  = Variable('\\rho_{bend}', 2700, 'kg/m^3', 'Stringer density') # [TAS]
        rhocargo = Variable('\\rho_{cargo}', 150, 'kg/m^3', 'Cargo density') # b757 freight doc
        rhocone  = Variable('\\rho_{cone}', 2700, 'kg/m^3',
                            'Cone material density') # [TAS]
        rhofloor = Variable('\\rho_{floor}', 2700, 'kg/m^3',
                            'Floor material density') # [TAS]
        rhoinf   = Variable('\\rho_{\\infty}', 0.38, 'kg/m^3',
                            'Air density (35,000ft)')
        rholugg  = Variable('\\rho_{lugg}', 100, 'kg/m^3', 'Luggage density')
        rhoskin  = Variable('\\rho_{skin}', 2700, 'kg/m^3', 'Skin density') # [TAS]
        sigfloor = Variable('\\sigma_{floor}', 30000/0.000145, 'N/m^2',
                            'Max allowable cap stress') # [TAS]
        sigskin  = Variable('\\sigma_{skin}', 15000/0.000145, 'N/m^2',
                            'Max allowable skin stress') # [TAS]
        SPR      = Variable('SPR', 6, '-', 'Number of seats per row')
        taufloor = Variable('\\tau_{floor}', 30000/0.000145, 'N/m^2',
                            'Max allowable shear web stress') # [TAS]
        Tcabin   = Variable('T_{cabin}', 300, 'K', 'Cabin temperature')
        Vinf     = Variable('V_{\\infty}', 234, 'm/s', 'Cruise velocity')
        waisle   = Variable('w_{aisle}', 0.51, 'm', 'Aisle width') # Boeing doc
        Wavgpass = Variable('W_{avg. pass}', 180, 'lbf',
                            'Average passenger weight')
        Wcargo   = Variable('W_{cargo}', 10000, 'N', 'Cargo weight')
        Wcarryon = Variable('W_{carry on}', 15, 'lbf', 'Ave. carry-on weight')
        Wchecked = Variable('W_{checked}', 40, 'lbf', 'Ave. checked bag weight')
        Wfix     = Variable('W_{fix}', 3000, 'lbf',
                            'Fixed weights (pilots, cockpit seats, navcom)')
        Wpseat   = Variable('W\'_{seat}', 150, 'N', 'Weight per seat')
        Wppfloor = Variable('W\'\'_{floor}', 60, 'N/m^2',
                            'Floor weight/area density') # [TAS]
        Wppinsul = Variable('W\'\'_{insul}', 22, 'N/m^2',
                            'Weight/area density of insulation material') # [TAS]
        Wpwindow = Variable('W\'_{window}', 145.*3, 'N/m',
                            'Weight/length density of windows') # [TAS]
        wseat    = Variable('w_{seat}', 0.5, 'm', 'Seat width') # Boeing doc
        wsys     = Variable('w_{sys}', 0.10, 'm',
                            'Width between cabin and skin for systems')
#       xfix     = Variable('xfix', 2.1, 'm', 'x-location of fixed weight')

        with SignomialsEnabled():
            objective = Wfuse*D

            constraints = [
                            # Geometry relations
                            lnose == xshell1,

                            # Cross section relations
                            2*Rfuse >= SPR*wseat + waisle + 2*wsys,
                            Askin >= 2*np.pi*Rfuse*tskin,# + 2*dRfuse*tskin,
#                            Afuse <= np.pi*Rfuse**2 + 2*Rfuse*dRfuse,
                            Afuse >= np.pi*Rfuse**2,# + (2./3)*wfloor*dRfuse
                            #        + (3./2)*hhold**2*dRfuse/wfloor
                            #        - (3./2)*hhold*dRfuse**2/wfloor,
                            tshell >= tskin*(1 + rE*fstring*rhoskin/rhobend),
#                           Ihshell >= (np.pi*Rfuse**2 + 2*np.pi*(dRfuse/2)**2)
                            #          *Rfuse*tshell,
#                           Ivshell >= np.pi*Rfuse**2 * Rfuse*tshell,

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
#                           xVcyl >= 0.5*(xshell1+xshell2)*Vcyl,
#                           xVnose >= 0.5*(xshell1)*Vnose,
#                           xVbulk >= (xshell2 + 0.5*dRfuse)*Vbulk,
                            Wskin >= rhoskin*g*(Vcyl + Vnose + Vbulk),
#                           xWskin >= rhoskin*g*(xVcyl + xVnose + xVbulk),
                            Wshell >= Wskin*(1 + fstring + fframe + ffadd),
#                           xWshell >= xWskin*(1 + fstring + fframe + ffadd),

                            # Cabin volume and buoyancy weight
                            rhocabin == (1/(R*Tcabin))*pcabin,
                            Vcabin >= Afuse*(lshell + 0.67*lnose + 0.67*Rfuse),
                            Wbuoy >= (rhocabin - rhoinf)*g*Vcabin, # [SP]

                           # Windows and insulation
                           Wwindow == Wpwindow * lshell,
#                          xWwindow >= 0.5* (xshell1+xshell2)*Wwindow,
                           Winsul >= Wppinsul*(1.1*np.pi*Rfuse*lshell
                                     + 0.55*(Snose + Sbulk)),
#                          xWinsul >= 0.5*(xshell1 + xshell2)*Winsul,

                           # Payload-proportional weights
                           Wapu == Wpay*fapu,
#                          xWapu == xapu*Wapu,
                           Wseat == Wpseat*nseat,
#                          xWseat >= 0.5*(xshell1 + xshell2)*Wseat,
                           Wpadd == Wpay*fpadd,
#                          xWpadd >= 0.5*(xshell1 + xshell2)*Wpadd,

                            # Fixed weight
#                          xWfix == xfix*Wfix,

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
#                           xWfloor >= 0.5*(xshell1 + xshell2)*Wfloor,

                            # Tail cone
                            Qv*(1 + plamv/2) >= Lvmax*bv*plamv/3,
                            plamv >= 1.6,
                            taucone == sigskin,
                            Qv == 2*Afuse*taucone*tcone, # matches TASOPT code
                            Vcone*Rfuse*(1+lamcone) >= 2*Qv/taucone*lcone, # [SP]
                            lamcone == 0.4, # TODO remove
                            lamcone == cvt/lcone,
                            Wcone >= rhocone*g*Vcone*(1 + fstring + fframe
                                                      + ffadd),
#                           xWcone >= 0.5*(xshell2 + xconend) * Wcone,

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
                            Ahold <= (2./3)*wfloor*hhold + hhold**3/(2*wfloor),
                            # [SP] Harris stocker 1998 (wolfram)
                            dh + hhold + hfloor <= Rfuse,
                            Wpay >= Wpass + Wlugg + Wcargo,

                            # Total fuselage weight
                            Wfuse >= Wfix + Wapu + Wpadd + Wseat + Wshell
                                   + Wwindow + Winsul + Wcone + Wfloor + Wbuoy,
#                           xWfuse >= xWfix + xWapu + xWpadd + xWseat + xWshell
#                                   + xWcone + xWwindow + xWinsul + xWfloor

                            # Drag
                            # sources: Raymer (p285), kfid 325 notes (p180)
                            lfuse >= lnose+lshell+lcone,
                            f == lfuse/((4/np.pi*Afuse)**0.5), # fineness ratio
                            FF >= 1 + 60/f**3 + f/400, # form factor
                            Dfrict >= FF * np.pi*Rfuse * mu*Vinf
                                      * 0.074*(rhoinf*Vinf*lfuse/mu)**0.8,

                            # Drag due to fuselage upsweep (Raymer p286)
                            1.13226*phi**1.03759 == Rfuse/lcone, # monomial fit
                                                                 # of tan(phi)
                            Dupswp >= 3.83*phi**2.5*Afuse * 0.5*rhoinf*Vinf**2, 
                            D >= Dfrict + Dupswp
                          ]
            return objective, constraints

    def test(self):
        sol = self.localsolve()

        npt.assert_almost_equal(sol('A_{hold}'), (2./3)*sol('w_{floor}')
                                *sol('h_{hold}') + sol('h_{hold}')**3
                                /(2*sol('w_{floor}')), decimal=4)
#       npt.assert_almost_equal(sol('\\Delta h') + sol('h_{hold}')
#                               + sol('h_{floor}'), sol('R_{fuse}'))
#       npt.assert_almost_equal(sol('A_{fuse}'),  np.pi*sol('R_{fuse}')**2
#                               + 2*sol('R_{fuse}')*sol('\\Delta R_{fuse}'))
#       npt.assert_almost_equal(sol('A_{fuse}'), np.pi*sol('R_{fuse}')**2
#                               + (2./3)*sol('w_{floor}')*sol('\Delta R_{fuse}')
#                               + (3./2)*sol('h_{hold}')**2
#                               *sol('\Delta R_{fuse}')/sol('w_{floor}')
#                               - (3./2)*sol('h_{hold}')*sol('\Delta R_{fuse}')
#                               **2/sol('w_{floor}'))
        npt.assert_almost_equal(sol('W_{buoy}'), (sol('\\rho_{cabin}')
                                - sol('\\rho_{\\infty}'))*sol('g')
                                * sol('V_{cabin}'), decimal=1)
        npt.assert_almost_equal(sol('V_{cone}')*sol('R_{fuse}')
                                *(1 + sol('\\lambda_{cone}')), 2*sol('Q_v')
                                /sol('\\tau_{cone}') * sol('l_{cone}'))

#        with open('fuselage.tex', 'w') as outfile:
#            outfile.write(self.latex())

#        with open('vartable.tex', 'w') as outfile:
#            outfile.write(sol.table(latex=2))

#        with open('soltable.tex', 'w') as outfile:
#            outfile.write(sol.table(latex=3))

if __name__ == "__main__":
    Fuselage().test()
