"""Simple commercial aircraft flight profile and aircraft model"""
from numpy import pi, cos, tan
import numpy as np
from gpkit import Variable, Model, units, SignomialsEnabled, Vectorize
from gpkit.constraints.sigeq import SignomialEquality
from gpkit.tools import te_exp_minus1
from gpkit.tools.tools import te_tangent, te_secant
from gpkit.constraints.tight import Tight as TCS
from wingbox import WingBox

class Wing(Model):
    """
    Philippe's thesis wing model
    SKIP VERIFICATION
    """
    def setup(self, **kwargs):
        self.wns = WingNoStruct()
        self.wb = WingBox(self.wns, "wing")

        Wwing = Variable('W_{wing}', 'N', 'Wing System Weight')

        Cwing = Variable('C_{wing}', 1, '-', 'Wing Weight Margin and Sensitivity Factor')

        # w.r.t. the quarter chord of the root of the wing.
        dxACwing = Variable('\\Delta x_{AC_{wing}}','m','Wing Aerodynamic Center Shift')

        #wing induced drag reduction due to wing tip devices
        TipReduct = Variable('TipReduct', '-', 'Induced Drag Reduction Factor from Wing Tip Devices')

        constraints = []
        with SignomialsEnabled():
            constraints.extend([
            self.wns['\\lambda'] == self.wb['taper'],

            TCS([Wwing >= Cwing * self.wb['W_{struct}'] + self.wb['W_{struct}']*(self.wns['f_{flap}'] + \
                    self.wns['f_{slat}'] + self.wns['f_{aileron}'] + self.wns['f_{lete}'] + self.wns['f_{ribs}'] + \
                    self.wns['f_{spoiler}'] + self.wns['f_{watt}'])]),
            TCS([dxACwing <= 1./24.*(self.wns['c_{root}'] + 5.*self.wns['c_{tip}'])/self.wns['S'] \
                            *self.wns['b']**2*self.wns['\\tan(\\Lambda)']]),
            ])

        return self.wns, self.wb, constraints

    def dynamic(self, state):
        """
        returns an instance of the wing perofrmance model
        """
        return WingPerformance(self, state)

class WingNoStruct(Model):
    """
    Philippe's wing model minus structure
    SKIP VERIFICATION
    """
    def setup(self, **kwargs):
        #declare variables
               #Variables
        Afuel   = Variable('\\bar{A}_{fuel, max}', '-', 'Non-dim. fuel area')
        CLwmax  = Variable('C_{L_{w,max}}', '-', 'Max lift coefficient, wing')
        Vfuel   = Variable('V_{fuel, max}', 'm^3', 'Available fuel volume')
        L       = Variable('\\Lambda', '-', 'Wing quarter-chord sweep angle')
        cosL    = Variable('\\cos(\\Lambda)', '-',
                           'Cosine of quarter-chord sweep angle')
        croot   = Variable('c_{root}', 'm', 'Wing root chord')
        ctip    = Variable('c_{tip}', 'm', 'Wing tip chord')
        eta     = Variable('\\eta', '-', 'Lift efficiency (diff b/w sectional, actual lift)')
        fl      = Variable('f(\\lambda_w)', '-', 'Empirical efficiency function of taper')
        g       = Variable('g', 9.81, 'm/s^2', 'Gravitational acceleration')
        p       = Variable('p', '-', 'Substituted variable = 1 + 2*taper')
        q       = Variable('q', '-', 'Substituted variable = 1 + taper')
        rho0    = Variable('\\rho_0', 1.225, 'kg/m^3', 'Air density (0 ft)')
        rhofuel = Variable('\\rho_{fuel}', 'kg/m^3', 'Density of fuel')
        tanL    = Variable('\\tan(\\Lambda)', '-',
                           'Tangent of quarter-chord sweep angle')
        taper   = Variable('\\lambda', '-', 'Wing taper ratio')
        tau     = Variable('\\tau', '-', 'Wing thickness/chord ratio')
        tau_max = Variable('\\tau_{max_w}', '-', 'Max allowed wing thickness')
        wwn       = Variable('wwn', 0.5, '-', 'Wingbox-width-to-chord ratio')
        xw     = Variable('x_w', 'm', 'Position of wing aerodynamic center')
        ymac    = Variable('y_{mac}', 'm',
                           'Spanwise location of mean aerodynamic chord')
        bmax = Variable('b_{max}', 'm', 'Max Wing Span')

        AR      = Variable('AR', '-', 'Wing aspect ratio')
        Lmax    = Variable('L_{max}', 'N', 'Maximum load')
        Sw      = Variable('S', 'm^2', 'Wing area')
        WfuelWing   = Variable('W_{fuel_{wing}}', 'N', 'Fuel weight')
        b       = Variable('b', 'm', 'Wing span')
        mac    = Variable('mac', 'm',
                          'Mean aerodynamic chord (wing)')
        e       = Variable('e', '-', 'Oswald efficiency factor')
        FuelFrac = Variable('FuelFrac', '-', 'Usability Factor of Max Fuel Volume')

        #fractional component weights
        fflap = Variable('f_{flap}', '-', 'Flap Fractional Weight')
        fslat = Variable('f_{slat}', '-', 'Slat Fractional Weight')
        faile = Variable('f_{aileron}', '-', 'Aileron Fractional Weight')
        flete = Variable('f_{lete}', '-', 'Lete Fractional Weight')
        fribs = Variable('f_{ribs}', '-', 'Rib Fractional Weight')
        fspoi = Variable('f_{spoiler}', '-', 'Spoiler Fractional Weight')
        fwatt = Variable('f_{watt}', '-', 'Watt Fractional Weight')

        # Area fractions
        Atri = Variable('A_{tri}','m^2','Triangular Wing Area')
        Arect = Variable('A_{rect}','m^2','Rectangular Wing Area')

        #make constraints
        constraints = []

        with SignomialsEnabled():

            constraints.extend([
                 TCS([cosL*te_secant(L,5) <= 1]),
                 TCS([tanL >= te_tangent(L,5)]),
                 CLwmax <= 2.15/(cosL**2), # [TAS]
                 Arect == ctip*b,
                 Atri >= 1./2.*(1-taper)*croot*b, #[SP]
                 p >= 1 + 2*taper,
                 2*q >= 1 + p,
                 ymac == (b/3)*q/p,
                 TCS([(2./3)*(1+taper+taper**2)*croot/q <= mac],
                                   reltol=1E-2),
                 taper == ctip/croot,

                 SignomialEquality(Sw, b*(croot + ctip)/2),

                 # Oswald efficiency
                 # Nita, Scholz, "Estimating the Oswald factor from
                 # basic aircraft geometrical parameters"
                 TCS([fl >= 0.0524*taper**4 - 0.15*taper**3
                         + 0.1659*taper**2 - 0.0706*taper + 0.0119],
                    reltol=1E-2),
                TCS([e*(1 + fl*AR) <= 1]),

                taper >= 0.15, # TODO

                # Fuel volume [TASOPT doc]
                # GP approx of the signomial constraint:
                Vfuel <= .3026*mac**2*b*tau,

                WfuelWing <= rhofuel*Vfuel*g,

                b <= bmax,
                ])

        return constraints

class WingPerformance(Model):
    """
    Wing performance model
    SKIP VERIFICATION
    """
    def setup(self, wing, state, **kwargs):
        self.wing = wing

        #declare variables
        #Vector Variables
        alpha   = Variable('\\alpha_w', '-', 'Wing angle of attack')
        CLaw    = Variable('C_{L_{\\alpha,w}}', '-', 'Lift curve slope, wing')

        Re      = Variable('Re_w', '-', 'Reynolds number (wing)')
        CDp     = Variable('C_{D_{p_w}}', '-',
                           'Wing parasitic drag coefficient')
        CDi     = Variable('C_{D_{i_w}}', '-',
                           'Wing induced drag coefficient')
        CDw     = Variable('C_{d_w}', '-', 'Drag coefficient, wing')
        CLw     = Variable('C_{L}', '-', 'Lift coefficient, wing')

        D       = Variable('D_{wing}', 'N', 'Wing drag')
        Lw      = Variable('L_w', 'N', 'Wing lift')

        # Center wing section lift reduction variables
        dLo     = Variable('\\Delta L_{o}','N','Center wing lift loss')
        etao    = Variable('\\eta_{o}','-','Center wing span coefficient')
        po      = Variable('p_{o}','N/m','Center section theoretical wing loading')
        fLo     = Variable('f_{L_{o}}',0.5,'-','Center wing lift reduction coefficient')

        # Wing tip lift reduction variables
        dLt = Variable('\\Delta L_{t}','N','Wing tip lift loss')
        fLt = Variable('f_{L_{t}}',0.05,'-','Wing tip lift reduction coefficient')

        #wing moment variables -- need a good way to model this, currently using TAT
        cmw = Variable('c_{m_{w}}', '-', 'Wing Pitching Moment Coefficient')

        amax    = Variable('\\alpha_{max,w}', '-', 'Max angle of attack')

        #make constraints
        constraints = []

        with SignomialsEnabled():
            constraints.extend([
                0.5*state['\\rho']*state['V']**2*self.wing['S']*CLw >= Lw + dLo + 2.*dLt,
                dLo == etao*fLo*self.wing['b']/2*po,
                dLt == fLt*po*self.wing['c_{root}']*self.wing['taper']**2, # TODO improve approximations croot~co and taper~gammat

                # DATCOM formula (Mach number makes it SP)
                # Swept wing lift curve slope constraint
                SignomialEquality((self.wing['AR']/self.wing['\\eta'])**2*(1 + self.wing['\\tan(\\Lambda)']**2 - state['M']**2) + 8*pi*self.wing['AR']/CLaw
                      , (2*pi*self.wing['AR']/CLaw)**2),

                CLw == CLaw*alpha,
                alpha <= amax,

                # Drag
                D == 0.5*state['\\rho']*state['V']**2*self.wing['S']*CDw,
                TCS([CDw >= CDp + CDi]),
                TCS([CDi >= self.wing['TipReduct']*CLw**2/(pi*(self.wing['e'])*self.wing['AR'])]),
                Re == state['\\rho']*state['V']*self.wing['mac']/state['\\mu'],

                #original Philippe thesis fit
##                TCS([CDp**6.5 >= (1.02458748e10 * CLw**15.587947404823325 * (self.wing['\\cos(\\Lambda)']*state['M'])**156.86410659495155 +
##                                 2.85612227e-13 * CLw**1.2774976672501526 * (self.wing['\\cos(\\Lambda)']*state['M'])**6.2534328002723703 +
##                                 2.08095341e-14 * CLw**0.8825277088649582 * (self.wing['\\cos(\\Lambda)']*state['M'])**0.0273667615730107 +
##                                 1.94411925e+06 * CLw**5.6547413360261691 * (self.wing['\\cos(\\Lambda)']*state['M'])**146.51920742858428)]),
##                ])

            #Martin's TASOPT c series airfoil fit
            TCS([
                CDp**1.6515 >= 1.61418 * (Re/1000)**-0.550434 * (self.wing['\\tau'])**1.29151 * (self.wing['\\cos(\\Lambda)']*state['M'])**3.03609 * CLw**1.77743
                    + 0.0466407 * (Re/1000)**-0.389048 * (self.wing['\\tau'])**0.784123 * (self.wing['\\cos(\\Lambda)']*state['M'])**-0.340157 * CLw**0.950763
                    + 190.811 * (Re/1000)**-0.218621 * (self.wing['\\tau'])**3.94654 * (self.wing['\\cos(\\Lambda)']*state['M'])**19.2524 * CLw**1.15233
                    + 2.82283e-12 * (Re/1000)**1.18147 * (self.wing['\\tau'])**-1.75664 * (self.wing['\\cos(\\Lambda)']*state['M'])**0.10563 * CLw**-1.44114
                ]),
            ])

        return constraints
