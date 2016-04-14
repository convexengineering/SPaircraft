# coding=utf-8
"Implements an aircraft model composed of multiple sub-models"
from gpkit import Model, Variable, SignomialsEnabled, LinkedConstraintSet
from gpkit.constraints.tight import TightConstraintSet as TCS
from vtail import VerticalTail
from fuselage import Fuselage
from landing_gear import LandingGear
from htail import HorizontalTail

class Aircraft(Model):
    """
    Combined fuselage, tail, and landing gear model
    """

    def __init__(self):

        D      = Variable('D', 'N', 'Total aircraft drag (cruise)')
        Dfuse  = Variable('D_{fuse}', 'N', 'Fuselage drag')
        Dht    = Variable('D_{ht}', 'N', 'Horizontal tail drag')
        Dvt    = Variable('D_{vt}', 'N', 'Vertical tail drag')
        Dwing  = Variable('D_{wing}', 'N', 'Wing drag')
        W      = Variable('W', 'N', 'Total aircraft weight')
        Weng   = Variable('W_{eng}', 'N', 'Engine weight')
        Wfuse  = Variable('W_{fuse}', 'N', 'Fuselage weight')
        Wht    = Variable('W_{ht}', 'N', 'Horizontal tail weight')
        Wlg    = Variable('W_{lg}', 'N', 'Landing gear weight')
        Wvt    = Variable('W_{vt}', 'N', 'Vertical tail weight')
        Wwing  = Variable('W_{wing}', 'N', 'Wing weight')
        xCG    = Variable('x_{CG}', 'm', 'x-location of CG')
        xCGeng = Variable('x_{CG_{eng}}', 'm', 'x-location of engine CG')
        xCGfu  = Variable('x_{CG_{fu}}', 'm', 'x-location of fuselage CG')
        xCGht  = Variable('x_{CG_{ht}}', 'm', 'x-location of htail CG')
        xCGlg  = Variable('x_{CG_{lg}}', 'm', 'x-location of landing gear CG')
        xCGvt  = Variable('x_{CG_{vt}}', 'm', 'x-location of vtail CG') 
        xCGwing = Variable('x_{CG_{wing}}', 'm', 'x-location of wing CG')

        with SignomialsEnabled():

            objective = D + 0.5*W
            hlc = [# High level constraints
                   D >= Dvt + Dfuse       + Dwing + Dht,
                   W >= Wvt + Wfuse + Wlg + Wwing + Wht + Weng,
                   TCS([xCG*W >= Wvt*xCGvt + Wfuse*xCGfu + Wlg*xCGlg
                               + Wwing*xCGwing + Wht*xCGht + Weng*xCGeng], reltol=1E-2),
                  ]

            # Subsystem models
            vt = VerticalTail.aircraft_737()
            vts = VerticalTail.standalone_737()
            fu = Fuselage.aircraft_737()
            fus = Fuselage.standalone_737()
            lg = LandingGear.aircraft_737()
            lgs = LandingGear.standalone_737()
            ht = HorizontalTail.aircraft_737()
            hts = HorizontalTail.standalone_737()

            # Need to initialize solve with solution of uncoupled models
            vt_sol = vts.localsolve(verbosity=0)
            fu_sol = fus.localsolve(verbosity=0)
            lg_sol = lgs.localsolve(verbosity=0)
            ht_sol = hts.localsolve(verbosity=0)

            init = vt_sol['variables'].copy()
            init.update(fu_sol['variables'])
            init.update(lg_sol['variables'])
            init.update(ht_sol['variables'])
            init.update({
                         'x_{CG}': 15,
                         'x_{CG_{fu}}': 15,
                         'x_{CG_{ht}}': 38,
                         'x_{CG_{lg}}': 16,
                         'x_{CG_{vt}}': 35,
                        })

            self.init = init

        substitutions = {
                         'D_{wing}': 10000,
                         'W_{eng}': 10000,
                         'W_{wing}': 30000,
                         'x_{CG_{eng}}': 15,
                         'x_{CG_{wing}}': 15,
                        }

        Model.__init__(self, objective,
                             LinkedConstraintSet([hlc, vt, fu, lg, ht]),
                             substitutions)

    def test(self):
        sol = self.localsolve(x0=self.init)

if __name__ == "__main__":
    a = Aircraft()
    a.test()

