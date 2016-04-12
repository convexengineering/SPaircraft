# coding=utf-8
"Implements an aircraft model composed of multiple sub-models"
from gpkit import Model, Variable, SignomialsEnabled, LinkedConstraintSet
from vtail import VerticalTail
from fuselage import Fuselage
from landing_gear import LandingGear

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
        xCGeng  = Variable('x_{CG_{eng}}', 'm', 'x-location of engine CG')
        xCGfu  = Variable('x_{CG_{fu}}', 'm', 'x-location of fuselage CG')
        xCGht   = Variable('x_{CG_{ht}}', 'm', 'x-location of htail CG')
        xCGlg  = Variable('x_{CG_{lg}}', 'm', 'x-location of landing gear CG')
        xCGvt  = Variable('x_{CG_{vt}}', 'm', 'x-location of vtail CG') 
        xCGwing = Variable('x_{CG_{wing}}', 'm', 'x-location of wing CG')

        with SignomialsEnabled():

            objective = D + 0.5*W
            constraints = [# High level constraints
                           D >= Dvt + Dfuse       + Dwing + Dht,
                           W >= Wvt + Wfuse + Wlg + Wwing + Wht + Weng,
                          ]
            #               xCG*W >= Wvt*xCGvt + Wfuse*xCGfu + Wlg*xCGlg
            #                      + Wwing*xCGwing + Wht*xCGht + Weng*xCGeng,

            # Subsystem models
            vt = VerticalTail.aircraft_737()
            fu = Fuselage.aircraft_737()
            lg = LandingGear.aircraft_737()
            lgs = LandingGear.standalone_737()

            # Need to initialize solve with solution of uncoupled models
            vtfu = Model(vt.cost + fu.cost, LinkedConstraintSet([vt, fu]))
            vtfu.substitutions.update({
                                       'x_{CG_{fu}}': 15,
                                       'x_{CG_{vt}}': 35,
                                       'x_{CG}': 15
                                      })
            vtfu_sol = vtfu.localsolve(verbosity=0)
            vtfu.substitutions.clear()

            #lgs.substitutions.update({'x_{CG_{lg}}': 15, 'x_{CG}': 15})
            lg_sol = lgs.localsolve(verbosity=0)

            init = vtfu_sol['variables'].copy()
            init.update(lg_sol['variables'])

            self.init = init

        substitutions = {
                         'D_{ht}': 1000,
                         'D_{wing}': 10000,
                         'W_{eng}': 10000,
                         'W_{ht}': 5000,
                         'W_{wing}': 30000,
                         'x_{CG_{eng}}': 15,
                         'x_{CG_{ht}}': 38,
                         'x_{CG_{wing}}': 15,
                        }

        Model.__init__(self, objective,
                             LinkedConstraintSet([constraints, vt, fu, lg]),
                             substitutions)

    def test(self):
        sol = self.localsolve(x0=self.init)

if __name__ == "__main__":
    a = Aircraft()
    a.test()

