from gpkit import Variable, Model, units
from gpkit.constraints.tight import Tight as TCS

class WingBox(Model):
    """
    Structural model for a wing
    source: Hoburg, "Geometric Programming for Aircraft Design Optimization"

    Note - does not have a performance model
    """

    def setup(self, surface, surfacetype):
        # Variables
        Icap    = Variable('I_{cap}', '-',
                           'Non-dim spar cap area moment of inertia')
        Mr      = Variable('M_r', 'N', 'Root moment per root chord')
        nu      = Variable('\\nu', '-',
                           'Dummy variable = $(t^2 + t + 1)/(t+1)$')
        Wcap    = Variable('W_{cap}', 'N', 'Weight of spar caps')
        Wweb    = Variable('W_{web}', 'N', 'Weight of shear web')
        Wstruct = Variable('W_{struct}', 'N', 'Structural weight')

        # Constants
        taper = Variable('taper', '-', 'Taper ratio')
        g      = Variable('g', 9.81, 'm/s^2', 'Gravitational acceleration')
        Nlift  = Variable('N_{lift}', 1.0, '-', 'Wing loading multiplier')
        rh     = Variable('r_h', 0.75, '-',
                          'Fractional wing thickness at spar web')
        rhocap = Variable('\\rho_{cap}', 2700, 'kg/m^3',
                          'Density of spar cap material')
        rhoweb = Variable('\\rho_{web}', 2700, 'kg/m^3',
                          'Density of shear web material')
        sigmax = Variable('\\sigma_{max}', 250e6, 'Pa',
                          'Allowable tensile stress')
        sigmaxshear = Variable('\\sigma_{max,shear}', 167e6, 'Pa',
                               'Allowable shear stress')
        wwb      = Variable('r_{w/c}', 0.5, '-', 'Wingbox width-to-chord ratio')
        tcap    = Variable('t_{cap}' ,'-', 'Non-dim. spar cap thickness')
        tweb    = Variable('t_{web}', '-', 'Non-dim. shear web thickness')
        tau_max = Variable('\\tau_{max_w}', '-', 'Max allowed wing thickness')

        objective = Wstruct

        if surfacetype == "wing":
            AR = surface['AR']
            b = surface['b']
            S = surface['S']
            p = surface['p']
            q = surface['q']
            tau = surface['\\tau']
            Lmax = surface['L_{max}']

        constraints = [
                       # Aspect ratio definition
                       AR == b**2/S,

                       # Upper bound on maximum thickness
                       tau <= tau_max,

                       # Root stiffness (see Hoburg 2014)
                       # Assumes rh = 0.75, so that rms box height = ~0.92*tmax
                       TCS([0.92*wwb*tau*tcap**2 + Icap <= 0.92**2/2*wwb*tau**2*tcap]),

                       # Stress limit
                       # Assumes bending stress carried by caps (Icap >> Iweb)
                       TCS([8 >= Nlift*Mr*AR*q**2*tau/(S*Icap*sigmax)]),

                       # Shear web sizing
                       # Assumes all shear loads are carried by web and rh=0.75
                       TCS([12 >= AR*Lmax*Nlift*q**2/(tau*S*tweb*sigmaxshear)]),

                       # Posynomial approximation of nu=(1+lam+lam^2)/(1+lam^2)
                       # nu**3.94 >= 0.86*p**(-2.38)+ 0.14*p**0.56, # Woody's fit
                       (nu/1.09074074)**.166 >= 0.205*(p/1.7)**0.772
                                              + 0.795*(p/1.7)**-0.125, # Berk's fit

                       # Weight of spar caps and shear webs
                       Wcap >= 8*rhocap*g*wwb*tcap*S**1.5*nu/(3*AR**0.5),
                       Wweb >= 8*rhoweb*g*rh*tau*tweb*S**1.5*nu/(3*AR**0.5),
                      ]


        if surfacetype in ("wing", "horizontal_tail"):
            constraints += [
                    # Total wing weight using an additional weight fraction
                    Wstruct >= (Wweb + Wcap),
                    ]
        elif surfacetype == "vertical_tail":
            constraints += [
                    Wstruct >= 0.5*(Wweb + Wcap),
                    # Root moment calculation (see Hoburg 2014)
                    # Assumes lift per unit span proportional to local chord
                    Mr >= Lmax*AR*p/24
                    ]
        
        return constraints
