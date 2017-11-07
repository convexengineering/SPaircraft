Wing Model
==========

The aircraft wing is designed generate sufficient lift such that the aircraft can
takeoff, climb, cruise, descend, and land safely. Typically the wings also carry fuel tanks and
support the engines. Unfortunately, wings are heavy and produced drag. The purpose of this
model is to relate all of these considerations.

Model Assumptions
-----------------

The wing model assumes a continuous-taper, low-wing configuration with a
modern transonic airfoil. It does not currently consider wing twist or
wing dihedral. It also does not consider roll or yaw stability.

Model Description
-----------------

Variable tables are available for download below:

* :download:`Free variables <tables/wi_freevars.pdf>`

* :download:`Fixed variables <tables/wi_fixedvars.pdf>`

Wing Geometry
~~~~~~~~~~~~~

Before considering a wing’s performance, the variables that prescribe
its geometry must be appropriately constrained.

The relationship between reference area, span and mean geometric chord
is enforced using a constraint that assumes a trapezoidal planform. This
constraint is implemented as a signomial equality constraint because
there is both upward and downward (optimization) pressure on the
reference area, and it is not possible to know a priori which will
dominate.

.. math:: {S_{w}} = {b_{w}} \frac{c_{root_{w}} + c_{tip_{w}}}{2} \label{eq:planformarea}

The mean aerodynamic chord relationship for a trapezoidal wing can be
written as a signomial constraint, and its spanwise location can be
written as a monomial equality constraint. These constraints make use of
dummy variables, :math:`p_w` and :math:`q_w`, introduced by the
structural model below.

.. math::

   \begin{aligned}
   \bar{c}_{w} &\leq \frac23 \left(\frac{1 + \lambda_{w} 
   + \lambda_{w}^2}{q_{w}}\right) c_{root_{w}} \label{eq:meanaerochord} \\
   y_{\bar{c}_w} &= \frac{b_w q_w}{3 p_w} \label{eq:spanwisemac}\end{aligned}

The wing taper ratio is defined by a monomial equality constraint. It
is necessary to lower bound taper to avoid an unacceptably small
Reynolds number at the wing tip [Kroo, 2001].
For the purpose of this work, the taper is lower bounded using the taper
ratio of the reference aircraft’s wing.

.. math::

   \begin{aligned}
   \lambda_{w} &= \frac{c_{tip_{w}}}{c_{root_{w}}} \label{eq:taperratio}\\
   {\lambda_{w}} &\geq \lambda_{w_{min}} \label{eq:mintaperratio}\end{aligned}

Finally, a maximum span constraint can be imposed to reflect, for
example, a gate size constraint.

.. math:: b_w \leq b_{w,max}

Wing Lift
~~~~~~~~~

Total lift is constrained to be greater than the weight of the aircraft
plus the downforce from the horizontal tail. The constant
:math:`f_{L_{total/wing}}` is greater than one and used to account for
fuselage lift.

.. math::

   \begin{aligned}
   L_{total} &\geq W + L_{ht}\\
   L_{total} &= f_{L_{total/wing}} L_{w}\end{aligned}

The standard equation for the lift of a wing is a natural monomial
equality constraint.

.. math::

   \begin{aligned}
   L_w = \frac12 \rho_{\infty} V_{\infty}^2 S_w C_{L_w}\end{aligned}

However, this assumes a continuous unobstructed wing planform.
Correcting for lift loss at the fuselage and at the wing tips, gives the
adjusted Equation , which can be rearranged into the posynomial
Constraint .

.. math::

   \begin{aligned}
   L_w &= \frac12 \rho_{\infty} V_{\infty}^2 S_w C_{L_w} - \Delta L_o - 2 \Delta L_t 
   \label{eq:liftadjeq} \\
   \frac12 \rho_{\infty} V_{\infty}^2 S_w C_{L_w} &\geq L_w + \Delta L_o + 2 \Delta L_t
   \label{eq:liftadjcon}\end{aligned}

The lift corrections are given as monomial equality
constraints [Drela, 2011].

.. math::

   \begin{aligned}
   \Delta L_o &= \eta_o f_{L_o} \frac{b_w}{2} p_o \\
   \Delta L_t &= f_{L_t} p_o c_{root_{w}} \lambda_w^2\end{aligned}

The lift coefficient of the wing goes linearly with the angle of
attack, which is limited by a maximum angle of attack due to stall.

.. math::

   \begin{aligned}
       C_{L_w} &= C_{L_{\alpha,w}}\alpha_w \\
       \alpha_{w} &\leq \alpha_{w,max}\end{aligned}

The DATCOM formula is an analytic function for estimating the lift
curve slope of a wing or tail, based on empirical
results [Kroo, 2001].

.. math:: C_{L_{\alpha,w}} = \frac{2 \pi AR_{w}}{2+\sqrt{(AR_{w}/\eta_w)^2(1+\tan^2\Lambda - M^2)+4}}

This relationship can be used as a signomial inequality to constrain
the lift curve slope, although some algebraic manipulation is needed.

.. math::

   \begin{aligned}
   C_{L_{\alpha,w}} &\leq \frac{2 \pi AR_{w}}{2+\sqrt{(AR_{w}/\eta_w)^2(1+\tan^2\Lambda-M^2)+4}}
    \\
   (AR_{w}/\eta_w)^2(1+\tan^2\Lambda - M^2)+4 &\leq \left( \frac{2\pi AR_{w}}{C_{L_{\alpha,w}}}
    - 2 \right)^2 \\
   (AR_{w}/\eta_w)^2(1+\tan^2\Lambda - M^2) &\leq \frac{4 \pi^2 AR_{w}^2}{C_{L_{\alpha,w}}^2}
    - \frac{8 \pi AR_{w}}{C_{L_{\alpha,w}}} \\
   \frac{C_{L_{\alpha,w}}^2}{\eta_w^2}\left(1 + \tan^2\Lambda - M^2\right) +
   \frac{8\pi C_{L_{\alpha,w}}}{AR_{w}} &\leq 4\pi^2 \end{aligned}

Maximum wing lift is constrained using an assumed load factor,
:math:`N_{lift}`.

.. math::

   \label{e:Lmax}
   f_{L_{total/wing}} L_{w_{max}} \geq N_{lift} W_{total} + L_{ht_{max}}

Finally, wing loading is constrained to be less than a user specified
maximum.

.. math::

   \begin{aligned}
   W_{S} &= \frac{1}{2} \rho_{\infty} C_{L_w} {V_{\infty}}^2 \\
   W_{S} &\leq W_{S_{max}}\end{aligned}

Wing Weight
~~~~~~~~~~~

Wing weight is constrained to be greater than the wing structural weight
plus a series of fractional weights to account for wing ribs and control
surfaces.

.. math::

   W_{wing} \geq W_{struct_{w}}(1 + f_{flap} + f_{slat} + f_{aileron}
                           + f_{lete} + f_{ribs} + f_{spoiler} + f_{watt})

Wing structural weight is constrained using an adaptation of the
structural model from [Hoburg, 2013],
which comprises 12 monomial and posynomial constraints.

.. math::

   \begin{aligned}
   {W_{struct_{w}}} &\geq ({W_{cap}} + {W_{web}}) \\
   {W_{cap}} &\geq  \frac{8{\rho_{cap}} {g} {w} {t_{cap}} {S_{w}}^{1.5} {\nu}}{3{AR_{w}}^{0.5}} \\
   {W_{web}} &\geq \frac{8{\rho_{web}}{g}{r_h}{\tau_{w}}{t_{web}}{S_{w}}^{1.5}{\nu}}{3{AR_{w}}^{0.5}} \\
   {\nu}^{3.94} &\geq 0.14{p_{w}}^{0.56} + \frac{0.86}{{p_{w}}^{2.4}} \\
   {p_{w}} &\geq 1 + 2{\lambda_{w}} \\
   2{q_{w}} &\geq 1 + {p_{w}} \\
   \frac{0.92^2}{2}{\tau_{w}}^{2}{t_{cap}}{w} &\geq 0.92{\tau_{w}}{t_{cap}}^{2}{w} + {I_{cap}} \\
   \frac{{AR_{w}} {M_r} {N_{lift}} {\tau_{w}} {q_{w}}^{2}}{{I_{cap}} {S_{w}} {\sigma_{max}}} &\leq 8 \\
   \frac{{AR_{w}}{L_{w_{max}}}{N_{lift}}{q_{w}}^{2}}{{S_{w}}{\sigma_{max,shear}}{\tau_{w}}{t_{web}}}&\leq 12 \\
   {AR_{w}} &= \frac{{b_{w}}^{2}}{{S_{w}}} \\
   {\tau_{w}} &\leq 0.14 \end{aligned}

The original root bending moment constraint,

.. math:: {M_r} \geq \frac{{AR_{w}} {L_{w_{max}}} {p_{w}}}{24},

is replaced with a more sophisticated signomial constraint that
considers the load relief effect due to the weight of the engine and the
fuel tanks. To derive the constraint, the lift per unit span of wing is
assumed to be proportional to the local chord, and the wing planform
area is partitioned into an untapered (rectangular) area
:math:`A_{rect}` and a fully tapered (triangular) area :math:`A_{tri}`.

.. math::

   \begin{aligned}
   A_{tri} &= \frac{1}{2} (1-\lambda_w) c_{root_{w}} b_w \\
   A_{rect} &= c_{tip_{w}} b_w\end{aligned}

The wing area component loads are treated as point loads to determine
the equivalent wing root moment.

.. math::

   \begin{aligned}
   \label{eq:M_rinit}
   M_r c_{root_{w}} \geq &\left(L_{w_{max}} - N_{lift}\left(W_{wing} + f_{fuel,wing}
              W_{fuel}\right)\right) \left(\frac16 A_{tri} + \frac14
              A_{rect}\right)\frac{b_{w}}{S_{w}} \\
              &- N_{lift} W_{engine} y_{eng} \nonumber\end{aligned}

This constraint can be further simplified to remove the need for
intermediary variables :math:`A_{tri}` and :math:`A_{rect}`, since

.. math::

   \begin{aligned}
   \frac{1}{6} A_{tri} + \frac{1}{4} A_{rect} &= \frac{1}{12} (c_{root_{w}} - c_{tip_{w}}) b_{w} 
   + \frac{1}{4} c_{tip_{w}} b_{w} \\
   &= \frac{b_{w}}{12} (c_{root_{w}} + 2 c_{tip_{w}}).
   \label{eq:Asub}\end{aligned}

Substituting Equation  into Constraint  yields the following wing root
moment constraint.

.. math::

   \begin{aligned}
   M_r c_{root_{w}} \geq &\left(L_{w_{max}} - N_{lift}\left(W_{wing} + f_{fuel,wing}
              W_{fuel}\right)\right) \left(\frac{b_{w}^2}{12 S_{w}} \left(c_{root_{w}} 
              + 2 c_{tip_{w}}\right)\right) \\
              & - N_{lift} W_{engine} y_{eng} \nonumber\end{aligned}

Note that this provides a conservative estimate for the root moment,
since it assumes that the lift per unit area is constant throughout the
wing, whereas in reality the lift per unit area diminishes towards the
wingtips.

Wing Drag
~~~~~~~~~

Wing drag is captured by five monomial and posynomial constraints. The
parasitic drag coefficient is constrained using a softmax affine fit of
XFOIL[Drela, 1989]simulation data for the
TASOPT[Drela, 2011] C-series airfoils, which are
representative of modern transonic
airfoils[Drela, 2011]. The fit, which considers
wing thickness, lift coefficient, Reynolds number, and Mach number, was
developed with `GPfit <https://github.com/convexengineering/gpfit>`_ and has an RMS
error of approximately 5%. Constraint  is an adaption of the standard
definition of the induced drag
coefficient [Anderson, 2001], with an
adjustment factor for wingtip devices.

.. math::

   \begin{aligned}
   D_w &= \frac12 \rho_{\infty} V_{\infty}^2 S_w C_{D_w} \label{eq:wingdrag}\\
   C_{D_w} &\geq C_{D_{p_w}} + C_{D_{i_w}} \label{eq:wingdragcoeff}\\
   \label{eq:wingpdragcoeff}
   C_{D_{p_w}}^{1.65} &\geq 1.61  \left(\frac{Re_w}{1000}\right)^{-0.550}
           (\tau_{w})^{1.29}
           (M \cos(\Lambda))^{3.04}
           C_{L_w}^{1.78} \\
           &+ 0.0466  \left(\frac{Re_w}{1000}\right)^{-0.389}
           (\tau_{w})^{0.784}
           (M \cos(\Lambda))^{-0.340}
           C_{L_w}^{0.951} \nonumber \\
             &+ 191  \left(\frac{Re_w}{1000}\right)^{-0.219}
           (\tau_{w})^{3.95}
           (M\cos(\Lambda))^{19.3}
           C_{L_w}^{1.15} \nonumber \\
           &+ 2.82e-12  \left(\frac{Re_w}{1000}\right)^{1.18}
           (\tau_{w})^{-1.76}
           (M \cos(\Lambda))^{0.105}
           C_{L_w}^{-1.44} \nonumber \\
   \label{eq:wingRe}
   Re_w &= \frac{\rho_{\infty} V_{\infty} \bar{c}_w}{\mu} \\
   C_{D_{i_w}} &\geq f_{tip} \frac{C_{L_w}^2}{\pi e AR_{w}} \label{eq:induceddrag}\end{aligned}

The Oswald efficiency is constrained by a relationship
from [Nita, 2012], in which the authors fit a
polynomial function to empirical data. Given that all polynomials are
signomials, this can easily be used in the framework.

.. math:: e\leq \frac{1}{1 + f(\lambda_w) AR_{w} }

.. math::

   \label{eq:flambda}
   f(\lambda_w) \geq 0.0524 \lambda_w^4 - 0.15 \lambda_w^3 + 0.1659 \lambda_w^2 -
   0.0706 \lambda_w + 0.0119

The Oswald efficiency is plotted as a function of taper ratio, as
imposed by this pair of constraints, in .

.. figure:: figs/e_fit.eps
   :alt: Empirical relationship for Oswald efficiency as a function of
   taper for a wing with :math:`AR_{w}`\ =10

   Empirical relationship for Oswald efficiency as a function of taper
   for a wing with :math:`AR_{w}`\ =10

Wing Aerodynamic Center
~~~~~~~~~~~~~~~~~~~~~~~

The true aerodynamic center and the of the wing are shifted in the
aircraft’s x-axis with respect to the wing root quarter chord due to the
swept geometry of the wing. This effect is captured with the variable
:math:`\Delta x_{ac_w}`. Assuming that the wing lift per unit area is
constant, and by integrating the product of the local quarter chord
offset :math:`\delta x(y)` and local chord area :math:`c(y)dy` over the
wing-half span, it can be calculated by

.. math::

   \label{eq:dXACwingDerivation}
   \Delta x_{ac_w} = \frac{2}{S} \int_{0}^{b/2} c(y) \delta x(y) dy,

where the local root chord :math:`c(y)` and the local quarter chord
offset :math:`\delta x(y)` are given by:

.. math::

   \begin{aligned}
   \label{eq:cy}
   c(y) &= \left(1 - (1-\lambda_w) \frac{2y}{b_w} \right) c_{root_{w}} \\
   \label{eq:dxy}
   \delta x(y) &= y \tan(\Lambda)\end{aligned}

By substituting Equations and into Equation , expanding out the
integral and relaxing the equality, :math:`\Delta x_{ac_w}` can be
constrained as follows.

.. math:: \Delta x_{ac_w} \geq \frac{1}{4} \tan(\Lambda) AR_{w} c_{root_{w}} \left(\frac{1}{3} + \frac{2}{3} \lambda_w \right)

Fuel Volume
~~~~~~~~~~~

Fuel tanks are typically located inside the wingbox. Using the geometry
of a TASOPT-optimized 737-800[Drela, 2011], a
constraint on the maximum fuel volume in the wing was developed. For a
wing of the same mean aerodynamic chord, thickness, and span as a TASOPT
737-800, the maximum available fuel volumes in the wing will match
exactly. To allow for the possibility of auxiliary tanks in the
horizontal tail or fuselage the user-specified value
:math:`f_{fuel, usable}` is introduced.

.. math::

   \begin{aligned}
   \label{e:V_fuel}
   V_{fuel, max} &\leq 0.303 {\bar{c}_w}^2 b_{w} \tau_{w} \\
   W_{fuel_{wing}} &\leq \rho_{fuel} V_{fuel, max} g  \\
   W_{fuel_{wing}} &\geq \frac{f_{fuel, wing} W_{fuel_{total}}}{ f_{fuel, usable}}\end{aligned}
