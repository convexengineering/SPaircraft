from gpkit import Variable, Model

class Structure(Model):
    """
    Structural model for a wing
    source: Hoburg, Abbeel, "Geometric Programming for Aircraft Design Optimization", 2014
    """

    def setup(self, Lmax, fwadd=0.4, Nlift=2.0, rh=0.75, rhocap=2700, rhoweb=2700,
            sigmax=250E6, sigmaxshear=167E6, taper=0.45, w=0.5):

        if type(taper) != Variable:
            taper = Variable('taper', taper, '-', 'Taper ratio')
        if type(Lmax) != Variable:
            Lmax  = Variable('L_{max}', Lmax, 'N', 'Maximum load')

        # Variables
        A       = Variable('A', '-', 'Aspect ratio')
        b       = Variable('b', 'm', 'Span')
        Icap    = Variable('I_{cap}', '-', 'Non-dim spar cap area moment of inertia')
        Mr      = Variable('M_r', 'N', 'Root moment per root chord')
        nu      = Variable('\\nu', '-', 'Substituted variable = $(t^2 + t + 1)/(t+1)$')
        p       = Variable('p', '-', 'Substituted variable = 1 + 2*taper')
        q       = Variable('q', '-', 'Substituted variable = 1 + taper')
        S       = Variable('S', 'm^2', 'Reference area')
        tau     = Variable('\\tau', '-', 'Thickness to chord ratio')
        tcap    = Variable('t_{cap}' ,'-', 'Non-dim. spar cap thickness')
        tweb    = Variable('t_{web}', '-', 'Non-dim. shear web thickness')
        Wcap    = Variable('W_{cap}', 'N', 'Weight of spar caps')
        Wweb    = Variable('W_{web}', 'N', 'Weight of shear web')
        Wstruct = Variable('W_{struct}', 'N', 'Structural weight')

        # Constants
        fwadd  = Variable('f_{w,add}', fwadd, '-', 'Wing added weight fraction') # TASOPT code (737.tas)
        g      = Variable('g', 9.81, 'm/s^2', 'Gravitational acceleration')
        Nlift  = Variable('N_{lift}', Nlift, '-', 'Wing loading multiplier')
        rh     = Variable('r_h', rh, '-', 'Fractional wing thickness at spar web')
        rhocap = Variable('\\rho_{cap}', rhocap, 'kg/m^3', 'Density of spar cap material' )
        rhoweb = Variable('\\rho_{web}', rhoweb, 'kg/m^3', 'Density of shear web material')
        sigmax = Variable('\\sigma_{max}', sigmax, 'Pa', 'Allowable tensile stress')
        sigmaxshear = Variable('\\sigma_{max,shear}', sigmaxshear, 'Pa', 'Allowable shear stress')
        w      = Variable('w', w, '-', 'Wingbox-width-to-chord ratio')

        objective = Wstruct

        constraints = [
                       # Aspect ratio definition
                       A == b**2/S,

                       # Defining taper dummy variables
                       p >= 1 + 2*taper,
                       2*q >= 1 + p,

                       # Upper bound on maximum thickness
                       tau <= 0.15,

                       # Root moment calculation (see Hoburg 2014)
                       # Depends on a given load the wing must support, Lmax
                       # Assumes net lift per unit span proportional to the local chord
                       Mr >= Lmax*A*p/24,

                       # Root stiffness (see Hoburg 2014)
                       # Assumes rh = 0.75, so that rms box height = ~0.92*max thickness
                       0.92*w*tau*tcap**2 + Icap <= 0.92**2/2*w*tau**2*tcap,

                       # Stress limit
                       # Assumes bending stress carried by caps only (Icap >> Iweb)
                       8 >= Nlift*Mr*A*q**2*tau/(S*Icap*sigmax),

                       # Shear web sizing
                       # Assumes all shear loads are carried by web and rh = 0.75
                       12 >= A*Lmax*Nlift*q**2/(tau*S*tweb*sigmaxshear),

                       # Posynomial approximation of nu = (1+lam+lam**2)/(1+lam**2)
                       nu**3.94 >= 0.86*p**(-2.38)+ 0.14*p**0.56,

                       # Weight of spar caps and shear webs
                       Wcap >= 8*rhocap*g*w*tcap*S**1.5*nu/(3*A**0.5),
                       Wweb >= 8*rhoweb*g*rh*tau*tweb*S**1.5*nu/(3*A**0.5),

                       # Total wing weight using an additional weight fraction
                       Wstruct >= fwadd*(Wweb + Wcap),
                       ]
        return objective, constraints
