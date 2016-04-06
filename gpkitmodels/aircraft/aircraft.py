# coding=utf-8
"Implements an aircraft model composed of multiple sub-models"
from gpkit import Model, Variable, SignomialsEnabled, LinkConstraint
from gpjet.aircraft.vtail import VerticalTail
from gpjet.aircraft.fuselage import Fuselage
from gpjet.aircraft.landing_gear import LandingGear

class Aircraft(Model):
    """
    Combined fuselage, tail, and landing gear model
    """

    def __init__(self):

        # Free variables
        W      = Variable('W', 'N', 'Total aircraft weight')
        Wfuse  = Variable('W_{fuse}', 'N', 'Fuselage weight')
        Wlg    = Variable('W_{lg}', 'N', 'Landing gear weight')
        Wvt    = Variable('W_{vt}', 'N', 'Vertical tail weight')
        xCG    = Variable('x_{CG}', 'm', 'x-location of CG')
        xCGfu  = Variable('x_{CG_{fu}}', 'm', 'x-location of fuselage CG')
        xCGlg  = Variable('x_{CG_{lg}}', 'm', 'x-location of landing gear CG')
        xCGvt  = Variable('x_{CG_{vt}}', 'm', 'x-location of tail CG') 

        # Fixed variables
        Weng    = Variable('W_{eng}', 10000, 'N', 'Engine weight')
        Wht     = Variable('W_{ht}', 5000, 'N', 'Horizontal tail weight')
        Wwing   = Variable('W_{wing}', 30000, 'N', 'Wing weight')
        xCGeng  = Variable('x_{CG_{eng}}', 15, 'm', 'x-location of engine CG')
        xCGht   = Variable('x_{CG_{ht}}', 38, 'm', 'x-location of horizontal tail CG')
        xCGwing = Variable('x_{CG_{wing}}', 15, 'm', 'x-location of wing CG')

        with SignomialsEnabled():

            # High level model
            #objective = W
            #constraints = [W <= Wvt + Wfuse + Wlg + Wwing + Wht + Weng,
            #               xCG*W >= Wvt*xCGvt + Wfuse*xCGfu + Wlg*xCGlg
            #                      + Wwing*xCGwing + Wht*xCGht + Weng*xCGeng,
            #              ]
            #m = Model(objective, constraints)

            # Subsystem models
            vt = VerticalTail.aircraft_737()
            fu = Fuselage.aircraft_737()
            lg = LandingGear.aircraft_737()
            lgs = LandingGear.standalone_737()

            # Need to initialize solve with solution of uncoupled models
            vtfu = Model(vt.cost + fu.cost, LinkConstraint([vt, fu]))
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

        Model.__init__(self, vt.cost + fu.cost + 0.5*lg.cost , LinkConstraint([vt, fu, lg]))

    def test(self):
        sol = self.localsolve(x0=self.init)

if __name__ == "__main__":
    a = Aircraft()
    a.test()

