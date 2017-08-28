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
                           'Dummy variable = $(t^2 + t + 1)/(t+1)^2$')
        Wcap    = Variable('W_{cap}', 'N', 'Weight of spar caps')
        Wweb    = Variable('W_{web}', 'N', 'Weight of shear web')
        Wstruct = Variable('W_{struct}', 'N', 'Structural weight')

        # Constants
        g      = Variable('g', 9.81, 'm/s^2', 'Gravitational acceleration')
        Nlift  = Variable('N_{lift}', '-', 'Wing loading multiplier')
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

        objective = Wstruct

        if surfacetype == "wing":
            taper = Variable('taper', '-', 'Taper ratio')
            AR = surface['AR']
            b = surface['b']
            S = surface['S']
            p = surface['p']
            q = surface['q']
            tau = surface['\\tau']
            tau_max = surface['\\tau_{max_w}']
            Lmax = surface['L_{max}']
        elif surfacetype == "vertical_tail":
            taper = Variable('taper', '-', 'Taper ratio')
            #factors of 2 required since the VT is only half span and the wing model
            #is for a full span wing
            AR = Variable('AR_{vt}','-','Vertical tail aspect ratio (double span)')
            b = 2.*surface['b_{vt}']
            S = 2.*surface['S_{vt}']
            p = surface['p_{vt}']
            q = surface['q_{vt}']
            tau = surface['\\tau_{vt}']
            Lmax = 2.*surface['L_{vt_{max}}']
            taper = surface['\\lambda_{vt}']
        elif surfacetype == "horizontal_tail":
            taper = Variable('taper', 0.3, '-', 'Taper ratio')
            AR = surface['AR_{ht}']
            b = surface['b_{ht}']
            S = surface['S_{ht}']
            p = surface['p_{ht}']
            q = surface['q_{ht}']
            tau = surface['\\tau_{ht}']
            Lmax = surface['L_{ht_{max}}']

            # Pi tail sizing variables
            # Splits the max lift into triangular and rectangular components
            # for root bending sizing.
            bhtout = Variable('b_{ht_{out}}','m','Horizontal tail outboard half-span')
            Lhtri = Variable('L_{ht_{tri}}','N','Triangular HT load')
            Lhrect = Variable('L_{ht_{rect}}','N','Rectangular HT load')
            Lhtriout = Variable('L_{ht_{tri_{out}}}','N','Triangular HT load outboard')
            Lhrectout = Variable('L_{ht_{rect_{out}}}','N','Rectangular HT load outboard')
            Mrout = Variable('M_{r_{out}}','N','Wing moment at pin joint ') #TODO
            Lshear = Variable('L_{shear}','N','Maximum shear load (at pin joint)')
            piMfac = Variable('\\pi_{M-fac}','-','Pi-tail bending structural factor')

        constraints = [
                       # Aspect ratio definition
                       AR == b**2/S,

                       # Root stiffness (see Hoburg 2014)
                       # Assumes rh = 0.75, so that rms box height = ~0.92*tmax
                       TCS([0.92*wwb*tau*tcap**2 + Icap <= 0.92**2/2*wwb*tau**2*tcap]),

                       # Posynomial approximation of nu=(1+lam+lam^2)/(1+lam)^2
                       nu**3.94 >= 0.86*p**(-2.38)+ 0.14*p**0.56, # Woody's fit

                       # Weight of spar caps and shear webs
                       Wweb >= 8*rhoweb*g*rh*tau*tweb*S**1.5*nu/(3*AR**0.5),
                      ]


        if surfacetype == "wing":
            constraints += [
                    # Upper bound on maximum thickness
                    tau <= tau_max,
                    Wstruct >= (Wweb + Wcap),
                    # Shear web sizing
                    # Assumes all shear loads are carried by web and rh=0.75
                    TCS([12 >= AR*Lmax*q**2/(tau*S*tweb*sigmaxshear)]),
                    # Weight of spar caps and shear webs
                    Wcap >= 8*rhocap*g*wwb*tcap*S**1.5*nu/(3*AR**0.5),
                    Nlift == 3,
                    # Stress limit
                    # Assumes bending stress carried by caps (Icap >> Iweb)
                    TCS([8 >= Mr*AR*q**2*tau/(S*Icap*sigmax)]),
                    ]
        elif surfacetype == "vertical_tail":
            constraints += [
                    # Upper bound on maximum thickness
                    tau <= 0.14,
                    Wstruct >= 0.5*(Wweb + Wcap),
                    # Root moment calculation (see Hoburg 2014)
                    # Assumes lift per unit span proportional to local chord
                    Mr >= Lmax*AR*p/24,
                    # Shear web sizing
                    # Assumes all shear loads are carried by web and rh=0.75
                    TCS([12 >= AR*Lmax*Nlift*q**2/(tau*S*tweb*sigmaxshear)]),
                    # Weight of spar caps and shear webs
                    Wcap >= 8*rhocap*g*wwb*tcap*S**1.5*nu/(3*AR**0.5),
                    Nlift == 1,
                    # Stress limit
                    # Assumes bending stress carried by caps (Icap >> Iweb)
                    TCS([8 >= Nlift*Mr*AR*q**2*tau/(S*Icap*sigmax)]),
                    ]
        elif surfacetype == "horizontal_tail":
            constraints += [
                    # Upper bound on maximum thickness
                    tau <= 0.14,
                    Wstruct >= (Wweb + Wcap),
                    Lhtriout >= Lhtri * bhtout**2 / (0.5*b)**2,
                    Lhrectout >= Lhrect * bhtout /(0.5*b),
                    12 >= 2*AR*Lshear*Nlift*q**2/(tau*S*tweb*sigmaxshear), #TODO
                    Wcap >= piMfac*8*rhocap*g*wwb*tcap*S**1.5*nu/(3*AR**0.5), #TODO
                    Nlift == 1,
                    # Stress limit
                    # Assumes bending stress carried by caps (Icap >> Iweb)
                    TCS([8 >= Nlift*Mr*AR*q**2*tau/(S*Icap*sigmax)]),
                    ]
        
        return constraints
