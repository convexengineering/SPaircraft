# coding=utf-8
"Implements an aircraft model composed of multiple sub-models"
from gpkit import Model, Variable, SignomialsEnabled, LinkedConstraintSet
from gpkit.constraints.tight import TightConstraintSet as TCS
from vtail import VerticalTail
from fuselage import Fuselage
from landing_gear import LandingGear
from htail import HorizontalTail
from wing import Wing
from wingbox import WingBox
from gpkit.tools import te_exp_minus1

class Aircraft(Model):
    """
    Combined fuselage, tail, and landing gear model
    """

    def __init__(self):

        CL     = Variable('C_L', '-', 'Lift coefficient')
        CD     = Variable('C_D', '-', 'Drag coefficient')
        D      = Variable('D', 'N', 'Total aircraft drag (cruise)')
        Dfuse  = Variable('D_{fuse}', 'N', 'Fuselage drag')
        Dht    = Variable('D_{ht}', 'N', 'Horizontal tail drag')
        Dvt    = Variable('D_{vt}', 'N', 'Vertical tail drag')
        Dwing  = Variable('D_{wing}', 'N', 'Wing drag')
        LD     = Variable('\\frac{L}{D}', '-', 'Lift/drag ratio')
        Lh     = Variable('L_h', 'N', 'Horizontal tail downforce')
        Lw     = Variable('L_w', 'N', 'Wing lift')
        R      = Variable('Range', 'nautical_miles', 'Range')
        Sw     = Variable('S_w', 'm**2', 'Wing reference area')
        TSFC   = Variable('TSFC', 'lb/lbf/hr',
                          'Thrust specific fuel consumption')
        V      = Variable('V_{\\infty}', 'm/s', 'Cruise velocity')
        W      = Variable('W', 'N', 'Total aircraft weight')
        Weng   = Variable('W_{eng}', 'N', 'Engine weight')
        Wfuel  = Variable('W_{fuel}', 'N', 'Fuel weight')
        Wfuse  = Variable('W_{fuse}', 'N', 'Fuselage weight')
        Wht    = Variable('W_{ht}', 'N', 'Horizontal tail weight')
        Wlg    = Variable('W_{lg}', 'N', 'Landing gear weight')
        Wpay   = Variable('W_{pay}', 'N', 'Payload weight')
        Wvt    = Variable('W_{vt}', 'N', 'Vertical tail weight')
        Wwing  = Variable('W_{wing}', 'N', 'Wing weight')
        Wzf     = Variable('W_{zf}', 'N', 'Zero fuel weight')
        g      = Variable('g', 9.81, 'm/s^2', 'Gravitational acceleration')
        rho    = Variable('\\rho', 'kg/m^3', 'Air density')
        xCG    = Variable('x_{CG}', 'm', 'x-location of CG')
        xCGeng = Variable('x_{CG_{eng}}', 'm', 'x-location of engine CG')
        xCGfu  = Variable('x_{CG_{fu}}', 'm', 'x-location of fuselage CG')
        xCGht  = Variable('x_{CG_{ht}}', 'm', 'x-location of htail CG')
        xCGlg  = Variable('x_{CG_{lg}}', 'm', 'x-location of landing gear CG')
        xCGvt  = Variable('x_{CG_{vt}}', 'm', 'x-location of vtail CG') 
        xCGwing = Variable('x_{CG_{wing}}', 'm', 'x-location of wing CG')
        xw     = Variable('x_w', 'm', 'x-location of wing aerodynamic center')
        z_bre  = Variable('z_{bre}', '-', 'Breguet parameter')

        with SignomialsEnabled():

            objective = Wfuel
            hlc = [# High level constraints

                   # Drag and weight buildup
                   D >= Dvt + Dfuse       + Dwing + Dht,
                   Wzf >= Wvt + Wfuse + Wlg + Wwing + Wht + Weng + Wpay,
                   W >= Wzf + Wfuel,

                   # Range equation for a jet
                   D == 0.5*rho*V**2*Sw*CD,
                   W == 0.5*rho*V**2*Sw*CL,
                   LD == CL/CD,
                   Lw >= W, # TODO: add Lh
                   R <= z_bre*LD*V/(TSFC*g),
                   Wfuel/Wzf >= te_exp_minus1(z_bre, nterm=3),

                   # CG relationships
                   TCS([xCG*W >= Wvt*xCGvt + Wfuse*xCGfu + Wlg*xCGlg
                               + Wwing*xCGwing + Wht*xCGht + Weng*xCGeng],
                       reltol=1E-2, raiseerror=False),
                   TCS([0.99*xCG*W <= Wvt*xCGvt + Wfuse*xCGfu + Wlg*xCGlg
                               + Wwing*xCGwing + Wht*xCGht + Weng*xCGeng],
                       reltol=1E-2, raiseerror=False),
                   xw == xCGwing,
                   xCGeng == xCGwing,
                  ]

            # Subsystem models
            vt = VerticalTail.aircraft_737()
            fu = Fuselage.aircraft_737()
            lg = LandingGear.aircraft_737()
            ht = HorizontalTail.aircraft_737()
            wi = Wing.aircraft_737()
            wb = WingBox()

        substitutions = {
                         'Range': 3000,
                         'TSFC': 0.3,
                         'V_{\\infty}': 234,
                         'W_{eng}': 10000,
                        }

        lc = LinkedConstraintSet([hlc, vt, fu, lg, ht, wi],
                                 exclude=[vk.name for vk in wb.varkeys])
        Model.__init__(self, objective,
                             lc,
                             substitutions)

    def test(self):
        sol = self.localsolve()
        return sol

if __name__ == "__main__":
    a = Aircraft()
    sol = a.test()

