from gpkit import Variable, Model, units, SignomialsEnabled, Vectorize
from gpkit.constraints.sigeq import SignomialEquality
from numpy import pi
from gpkit.constraints.tight import Tight as TCS
import matplotlib.pyplot as plt
import numpy as np

class TestMission(Model):
    def setup(self):
        #create submodels
        state = TestState()
        ht = BasicHT()
        wing = TestWing()
        fuse = TestFuse()

        with Vectorize(4):
            wingP = TestWingPerformance(state, wing)
            htP = BasicHTPerformance(state, ht)
            W_fuel = Variable('W_{fuel}', 'N', 'Fuel Weight')
            W_start = Variable('W_{start}', 'N', 'Segment Start Weight')
            W_end = Variable('W_{end}', 'N', 'Segment End Weight')
            W_avg = Variable('W_{avg}', 'N', 'Average Segment Weight')
            alpha = Variable('\\alpha', '-', 'Angle of Attack')
            D = Variable('D', 'N', 'Drag')
            xcg     = Variable('x_{CG}', 'm', 'CG location')

 

        submodels = [state, ht, htP, wing, wingP, fuse]
        
        #define mission level variables
        #weights
        W_payload = Variable('W_{payload}', 'N', 'Payload Weight')
        W_dry = Variable('W_{dry}', 'N', 'Zero Fuel Aircraft Weight')
        Cd0 = Variable('C_{D_{0}}', '-', 'Profile Drag')

 

        with SignomialsEnabled():
        
            constraints = [
                    #mission level
                    TCS([D >= .5*state['\\rho']*wing['S']*state['V']**2* \
                    (Cd0 + wing['K']*wingP['C_{L}']**2) + .5*state['\\rho']*ht['Sh']*state['V']**2*(Cd0 + ht['Kh']*htP['C_{L_{h}}']**2)]),

                    #assumes a 2 hour flight w/TSFC = 0.5 hr^-1
                    W_fuel == 0.5 * D * 2,

                    #compute the lift
                    TCS([W_dry >= wing['W_{wing}']+ ht['W_{HT}'] + W_payload]),
                    TCS([W_start >= W_dry + W_fuel]),
                    TCS([W_end >= W_dry]),
                    W_avg == (W_start*W_end)**.5,
                    TCS([wingP['L'] >= W_avg + htP['L_{h}']]),

                    #lift coefficient constraints
                    wingP['C_{L}'] == 2*pi*alpha,
##                    htP['C_{L_{h}}'] <= 2.2*pi*alpha,
##                    htP['C_{L_{h}}'] >= 1.8*pi*alpha,

                    #arbitrary, sturctural model will remove the need for this constraint
                    ht['b_{h}'] <= .33*wing['b_{max}'],

                    #HT sizing constraints
                    #compute mrat, is a signomial equality
                    SignomialEquality(ht['m_{ratio}']*(1+2/wing['AR']), 1 + 2/ht['ARh']),

                    #tail volume coefficient
                    ht['V_{h}'] == ht['Sh']*ht['l_{h}']/(wing['S']*wing['MAC']),

                    #enforce max tail location is the end of the fuselage
                    ht['l_{h}'] <= fuse['l_{fuse}'],

                    #Stability constraint, is a signomial
##                    TCS([ht['SM_{min}'] + ht['\\Delta x_{CG}']/wing['MAC'] <= ht['V_{h}']*ht['m_{ratio}'] \
##                         + wingP['c_{m_{w}}']/wing['C_{L_{max}}'] + ht['V_{h}']*ht['CL_{h_{max}}']/wing['C_{L_{max}}']]),
                    SignomialEquality(ht['SM_{min}'] + ht['\\Delta x_{CG}']/wing['MAC'], ht['V_{h}']*ht['m_{ratio}'] \
                                      + wingP['c_{m_{w}}']/wing['C_{L_{max}}'] + ht['V_{h}']*ht['CL_{h_{max}}']/wing['C_{L_{max}}']),

                    # Trim condidtion for each flight segment
##                    TCS([wingP['x_{ac}']/wing['MAC'] <= wingP['c_{m_{w}}']/wingP['C_{L}'] + xcg/wing['MAC'] + ht['V_{h}']*(htP['C_{L_{h}}']/wingP['C_{L}'])]),
                    SignomialEquality(wingP['x_{ac}']/wing['MAC'], wingP['c_{m_{w}}']/wingP['C_{L}'] + xcg/wing['MAC'] + ht['V_{h}']*(htP['C_{L_{h}}']/wingP['C_{L}'])),

                    #modify this
##                    xcg == 15*units('m'),
                    wing['x_w'] == 21*units('m'),
                    fuse['x_{fuse}'] == fuse['l_{fuse}']/2,
                    fuse['W_{fuse}'] == W_payload,

                    TCS([xcg <= (wing['W_{wing}'] * wing['x_w'] + fuse['x_{fuse}'] * fuse['W_{fuse}'] + ht['l_{h}'] * ht['W_{HT}'])/W_dry])

                    ]

        return submodels, constraints

class TestState(Model):
    def setup(self):
        #atm
        rho = Variable('\\rho', 'kg/m^3', 'Air Density')
        #airspeed
        V = Variable('V', 'm/s', 'Airspeed')

class TestWingPerformance(Model):
    def setup(self, state, wing):
        #aero
        L = Variable('L', 'N', 'Wing Lift')
        CL = Variable('C_{L}', '-', 'Wing Lift Coefficient')
        
        #Moments
        cmw = Variable('c_{m_{w}}', '-', 'Wing Pitching Moment Coefficient')   #approximtaed as a constant via TAT

        xac = Variable('x_{ac}', 'm', 'Position of wing aerodynamic center')

        constraints = [
            L == .5 * state['\\rho'] * state['V']**2 * wing['S'] * CL,

            CL <= wing['C_{L_{max}}'],

            cmw == 1,

            xac == 20*units('m'),
            ]

        return constraints

class TestWing(Model):
    def setup(self):
        #weight variables
        W_wing = Variable('W_{wing}', 'N', 'Wing Weight')

        CLmax = Variable('C_{L_{max}}', '-', 'Max Wing Lift Coefficient')

        #wing geometry
        S = Variable('S', 'm^2', 'Wing Planform Area')
        AR = Variable('AR', '-', 'Aspect Ratio')
        b = Variable('b', 'm', 'Wing Span')
        b_max = Variable('b_{max}', 'm', 'Max Wing Span')
        MAC = Variable('MAC', 'm', 'Mean Aerodynamic Chord')

        K = Variable('K', '-', 'Induced Drag Parameter')
        e = Variable('e', '-', 'Oswald Efficiency')

        xw = Variable('x_w', 'm', 'Position of wing CG')

        constraints = [
            #Wing geometry
            S == b*MAC,
            AR == b/MAC,    
            b <= b_max,

            #drag constraints
            K == 1/(pi*e*AR),

            #wing weight constraint
            #based off of a raymer weight and 737 data from TASOPT output file
            (S/(124.58*units('m^2')))**.65 == W_wing/(105384.1524*units('N')),
            ]

        return constraints

class TestFuse(Model):
    def setup(self):
        #aircraft geometry
        lfuse = Variable('l_{fuse}', 'm', 'Fuselage Length')
        x_fuse = Variable('x_{fuse}', 'm', 'Fuselage CG Location')
        W_fuse = Variable('W_{fuse}', 'N', 'Fuselage Weight')

class BasicHT(Model):
    """
    Basic horizontal tail sizing model. Given a CG travel range the
    model enforces a minium static margin for aft CG and pitch trim
    at forward CG and max lift.
    """
    def setup(self):
        #define variables
        #weight variables
        W_HT = Variable('W_{HT}', 'N', 'Horizontal Tail Weight')

        #HT geometry
        Sh = Variable('Sh', 'm^2', 'HT Planform Area')
        ARh = Variable('ARh', '-', 'HT Aspect Ratio')
        Vh = Variable('V_{h}', '-', 'Horizontal Tail Volume Coefficient')
        bh = Variable('b_{h}', 'm', 'HT Span')
        MACh = Variable('MAC_{h}', 'm', 'HT Mean Aerodynamic Chord')

        #aircraft geometry
        lh = Variable('l_{h}', 'm', 'Horizontal Tail Location')

        #aero
        CLhmax = Variable('CL_{h_{max}}', '-', 'Max Tail Downforce Coefficient')
 
        Kh = Variable('Kh', '-', 'HT Induced Drag Parameter')
        eh = Variable('eh', '-', 'HT Oswald Efficiency')
        mrat = Variable('m_{ratio}', '-', 'Wing to Tail Lift Slope Ratio')

        #min static margin
        SMmin = Variable('SM_{min}', '-', 'Minimum Static Margin')
        dxcg = Variable('\\Delta x_{CG}', 'm', 'Max CG Travel Range')

        #make the constraints
        constraints = []

        with SignomialsEnabled():

            constraints.extend([
                #HT weight constraint
                #based off of a raymer weight and 737 data from TASOPT output file
                (Sh/(46.1*units('m^2')))**.65 == W_HT/(16064.7523*units('N')),

                #HT geometry
                Sh == bh*MACh,
                ARh == bh/MACh,
 
                #HT Drag
                Kh == 1/(pi*eh*ARh),
                ])

            return constraints

class BasicHTPerformance(Model):
    def setup(self, state, ht):
        L_h = Variable('L_{h}', 'N', 'Horizontal Tail Downforce')
        CLh = Variable('C_{L_{h}}', '-', 'Tail Downforce Coefficient')
        
        constraints = [
            CLh <= ht['CL_{h_{max}}'],

            L_h == .5 * state['\\rho'] * state['V']**2 * ht['Sh'] * CLh,
            ]

        return constraints

if __name__ == '__main__':
    PLOT = False
    
    substitutions = {
        'W_{payload}': 85100*9.81,
        'b_{max}': 30,
        'l_{fuse}': 30,
##        'c_{m_{w}}': 1,
        'C_{L_{max}}': 2,
        'CL_{h_{max}}': 2.5,
        '\\rho': .8,
        'V': 230,
        'SM_{min}': 0.5,
        'C_{D_{0}}': 0.05,
        'e': 0.9,
        'eh': 0.9,
        '\\Delta x_{CG}': 4,
    }

    mission = TestMission()
    m = Model(sum(mission['W_{fuel}']), mission, substitutions)

    sol = m.localsolve(solver="mosek", verbosity=4)

    if PLOT == True:
        #sweeps of cg range
        substitutions = {
            'W_{payload}': 85100*9.81,
            'b_{max}': 30,
            'l_{fuse}': 30,
            'c_{m_{w}}': 1,
            'C_{L_{max}}': 2,
            'CL_{h_{max}}': 2.5,
            '\\rho': .8,
            'V': 230,
            'SM_{min}': 0.5,
            'C_{D_{0}}': 0.05,
            'e': 0.9,
            'eh': 0.9,
            '\\Delta x_{CG}': ('sweep', np.linspace(.5,6,10)),
        }

        mission = TestMission()
        m = Model(mission['W_{fuel}'], mission, substitutions)

        solCGsweep = m.localsolve(solver="mosek", verbosity=1)

        plt.plot(solCGsweep('\\Delta x_{CG}'), solCGsweep('Sh'), '-r')
        plt.xlabel('CG Travel Range [m]')
        plt.ylabel('Horizontal Tail Area [m$^2$]')
        plt.title('Horizontal Tail Area vs CG Travel Range')
        plt.show()

        plt.plot(solCGsweep('\\Delta x_{CG}'), solCGsweep('V_{h}'), '-r')
        plt.xlabel('CG Travel Range [m]')
        plt.ylabel('Horizontal Tail Volume Coefficient')
        plt.title('Horizontal Tail Area vs CG Travel Range')
        plt.show()

        #sweeps of SMmin
        substitutions = {
            'W_{payload}': 85100*9.81,
            'b_{max}': 30,
            'l_{fuse}': 30,
            'c_{m_{w}}': 1,
            'C_{L_{max}}': 2,
            'CL_{h_{max}}': 2.5,
            '\\rho': .8,
            'V': 230,
            'SM_{min}': ('sweep', np.linspace(.05,1,10)),
            'C_{D_{0}}': 0.05,
            'e': 0.9,
            'eh': 0.9,
            '\\Delta x_{CG}': 2,
        }

        mission = TestMission()
        m = Model(mission['W_{fuel}'], mission, substitutions)

        solSMsweep = m.localsolve(solver="mosek", verbosity=1)

        plt.plot(solSMsweep('SM_{min}'), solSMsweep('Sh'), '-r')
        plt.xlabel('Minimum Allowed Static Margin')
        plt.ylabel('Horizontal Tail Area [m$^2$]')
        plt.title('Horizontal Tail Area vs Min Static Margin')
        plt.show()

        plt.plot(solSMsweep('SM_{min}'), solSMsweep('V_{h}'), '-r')
        plt.xlabel('Minimum Allowed Static Margin')
        plt.ylabel('Horizontal Tail Volume Coefficient')
        plt.title('Horizontal Tail Volume Coefficient vs Min Static Margin')
        plt.show()

        #sweeps of payload
        substitutions = {
            'W_{payload}': ('sweep', np.linspace(50000*9.81,85100*9.81,10)),
            'b_{max}': 30,
            'l_{fuse}': 30,
            'c_{m_{w}}': 1,
            'C_{L_{max}}': 2,
            'CL_{h_{max}}': 2.5,
            '\\rho': .8,
            'V': 230,
            'SM_{min}': .1,
            'C_{D_{0}}': 0.05,
            'e': 0.9,
            'eh': 0.9,
            '\\Delta x_{CG}': 2,
        }

        mission = TestMission()
        m = Model(sum(mission['W_{fuel}']), mission, substitutions)

        solPayloadsweep = m.localsolve(solver="mosek", verbosity=1)

        plt.plot(solPayloadsweep('W_{payload}'), solPayloadsweep('Sh'), '-r')
        plt.xlabel('Payload Weight [N]')
        plt.ylabel('Horizontal Tail Area [m$^2$]')
        plt.title('Horizontal Tail Area vs Payload Weight')
        plt.show()

        plt.plot(solPayloadsweep('W_{payload}'), solPayloadsweep('V_{h}'), '-r')
        plt.xlabel('Payload Weight [N]')
        plt.ylabel('Horizontal Tail Volume Coefficient')
        plt.title('Horizontal Tail Volume Coefficient vs Payload Weight')
        plt.show()
