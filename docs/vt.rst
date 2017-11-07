Vertical Tail Model
===================

At a conceptual design level, the purpose of an aircraft’s vertical tail
is two-fold. Firstly, it must provide stability in yaw. Secondly, it
must provide adequate yaw control authority in critical flight
conditions. For a multi-engine aircraft, the critical flight condition
is typically an engine failure at low speeds. The vertical tail must be
capable of providing sufficient sideforce in this
case [Raymer, 1992]. The vertical tail must also
provide adequate yaw rate acceleration during landing flare in crosswind
conditions. The design of the vertical tail is therefore coupled to the
size of the fuselage, the position of the engines, and the aircraft’s
moment of inertia.

Model Assumptions
-----------------

The high-level assumptions for this model are that the the horizontal
tail is mounted on the fuselage, so as to not require a reinforced
vertical tail structure, and that the aircraft has two engines.

Model Description
-----------------

Variable tables are available for download below:

* :download:`Free variables <tables/vt_freevars.pdf>`

* :download:`Fixed variables <tables/vt_fixedvars.pdf>`

Vertical Tail Geometry and Structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The variables that define geometry are illustrated in . The moment arm
of the vertical tail is the distance from the aircraft to the
aerodynamic center of the vertical tail, which is assumed to be at the
quarter chord. The moment arm is therefore upper bounded by the distance
from the to the leading edge of the tail at the root, the height of the
mean aerodynamic chord above the fuselage, the sweep angle, and the mean
aerodynamic chord.

.. math::

   l_{vt}\leq\Delta x_{lead_{vt}}+z_{\bar{c}_{vt}}{\tan(\Lambda_{LE})}+0.25\bar{c}_{vt}
   \label{eq:vtmomentarm}

The x-coordinates of the leading and trailing edge at the root are
related by the root chord. The tail trailing edge is upper bounded by
imposing a constraint that the tail root cannot extend beyond the end of
the fuselage. Together these constraints put an upper bound on the
moment arm of the tail based on the length of the fuselage.

.. math::

   \begin{aligned}
   {\Delta x_{trail_{vt}}} &\geq {\Delta x_{lead_{vt}}} + {c_{root_{vt}}} \label{eq:vtleading}\\
   {l_{fuse}} &\geq {x_{CG}} + {\Delta x_{trail_{vt}}} \label{eq:vttrailing}\end{aligned}

The vertical tail structure is sized by its maximum lift coefficient and
the never-exceed speed.

.. math::

   \begin{aligned}
   L_{vt_{max}} &= \frac12 \rho_{TO} V_{ne}^2 S_{vt}C_{L_{v,max}}\end{aligned}

The remaining geometry and structural constraints were already
introduced in the wing model.
Constraints [eq:planformarea,eq:meanaerochord,eq:spanwisemac,eq:taperratio,eq:mintaperratio]
are adapted to the vertical tail model to constrain its geometry, with
two minor modifications. Constraint  can be relaxed from a signomial
equality to a signomial inequality constraint, meanwhile Constraint 
needs to be implemented as a signomial equality constraint. The wing
structure model from [Hoburg, 2013] is also
reused, however, given that the vertical tail only has a half-span, the
definitions of :math:`b_{vt}`, :math:`S_{vt}`, and :math:`W_{vt}` differ
from those of their wing counterparts.

Engine-out Condition
~~~~~~~~~~~~~~~~~~~~

The first performance constraint specifies that the maximum moment
exerted by the tail must be greater than or equal to the moment exerted
by the engines in an engine-out condition, exacerbated by the windmill
drag of the engine that is
inoperative [Drela, 2011].

.. math:: {L_{vt,EO}}{l_{vt}} \geq {D_{wm}} {y_{eng}} + {T_e} {y_{eng}}

The worst case engine out condition is likely to occur during takeoff,
when the velocity is lowest but the engine force required to safely
complete takeoff is highest. The force exerted by the vertical tail in
this critical low speed case is constrained by its maximum lift
coefficient, its reference area and the minimum dynamic pressure. As a
conservative estimate, the :math:`V_1` speed is used because it is the
minimum speed after which a takeoff can be completed, following a
critical engine failure.

.. math:: {L_{vt,EO}} = \frac12{\rho_{TO}}{V_1}^{2} {S_{vt}} {C_{L_{vt,EO}}}

The 3D lift coefficient is constrained by the airfoil sectional lift
coefficient using finite wing
theory [Anderson, 2001].

.. math:: C_{L_{vt,EO}}\left(1 + \frac{c_{l_{vt,EO}}}{\pi e_{vt} AR_{vt}}\right) \leq c_{l_{vt,EO}}

The windmill drag can, to a first approximation, be lower bounded using
a drag coefficient and a reference
area [Drela, 2011], in this case the area of the
engine fan.

.. math::

   \begin{aligned}
   {D_{wm}} &\geq \frac12{\rho_{TO}}{V_1}^{2}  {A_{fan}} {C_{D_{wm}}}\end{aligned}

Crosswind Landing Condition
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The second performance constraint ensures the vertical tail can provide
adequate yaw rate acceleration in a crosswind landing, where the moment
of inertia was constrained at the system level (Section
[chap:full\_aircraft]). To provide a safety margin during cross-wind
landing, :math:`C_{L_{vt,landing}}` is taken to be 85% of takeoff
:math:`{C_{L_{vt}}}`.

.. math::

   \frac12{\rho_{TO}{V_{land}}^{2}} S_{vt} l_{vt} C_{L_{vt, landing}} \geq 
   \frac{\dot{r}_{req}}{I_{z}}

Vertical Tail Drag
~~~~~~~~~~~~~~~~~~

The vertical tail produces drag, regardless of the flight condition.
Neglecting any induced drag, the parasitic drag coefficient,
:math:`C_{D_{p_{vt}}}`, is set by a softmax affine fit of
XFOIL[Drela, 1989]data for the symmetric NACA 0008
through 0020 airfoils. The fit considers airfoil thickness, Mach number,
and Reynolds number. It was developed with
`GPfit <https://github.com/convexengineering/gpfit>`_ and has an RMS error of 1.31%.

.. math::

   \begin{aligned}
   {D_{vt}} &\geq \frac12  {\rho_{\infty}} {V_\infty}^{2}{S_{vt}}{C_{D_{p_{vt}}}} \\
   {C_{D_{p_{vt}}}}^{1.189} &\geq 2.44\times10^{-77} (Re_{vt})^{-0.528} (\tau_{vt})^{133.8} (M)^{1022.7} \\
   &+ 0.003 (Re_{vt})^{-0.410}  (\tau_{vt})^{1.22} (M)^{1.55} \nonumber\\
   &+ 1.967\times10^{-4} (Re_{vt})^{0.214} (\tau_{vt})^{-0.04} (M)^{-0.14}
       \nonumber\\
   &+ 6.590\times10^{-50} (Re_{vt})^{-0.498} (\tau_{vt})^{1.56} (M)^{-114.6}
       \nonumber\\
   {Re_{vt}} &= \frac{ {\rho_\infty} {V_\infty} {\bar{c}_{vt}}}{{\mu}}\end{aligned}
