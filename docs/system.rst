System-level Model
==================

The objective of the optimization problem presented in this work is to
minimize fuel consumption, or equivalently fuel weight,
:math:`W_{fuel}`, using an adaptation of the Breguet range formulation
introduced in [Hoburg, 2013]. The purpose of
the system-level model is threefold: it enforces system-level
performance constraints such as required range and minimum cruise speed,
it encodes weight and drag buildups, and it constrains system-level
properties such as the aircraft’s and moment of inertia. In doing these
things, it also couples the subsystem models.

Model Assumptions
-----------------

The model presented in this work is a set of constraints that describe
the performance and design of a conventional-configuration narrowbody
aircraft, with a simple cruise-only mission profile. A more
sophisticated mission profile is left for future work.

Model Description
-----------------

Variable tables are available for download below:

* :download:`Free variables <tables/ac_freevars.pdf>`

* :download:`Fixed variables <tables/ac_fixedvars.pdf>`

Flight Performance
~~~~~~~~~~~~~~~~~~

The Breguet range formulation is discretized over multiple cruise
segments to improve accuracy, meaning the constraints
from [Hoburg, 2013] apply during each of the
:math:`N` flight segments. The :math:`n` subscript is used to represent
the :math:`n^{th}` flight segment where :math:`n=1...N`. For
readability, these subscripts are not used in the remainder of the
manuscript, but still apply.

.. math::

   \begin{aligned}
   \sum_{n=1}^{N} R_{n} &\geq R_{req} \\
   R_{n+1} &= R_{n} \\
   R_{n} &\leq \frac{V_{\infty_{n}}}{n_{eng}c_{T_{n}} g} \frac{W_{{avg}_{n}}}{D_{n}} z_{bre_{n}}\\
   W_{fuel_{n}} &\geq \left(z_{bre_{n}} + \frac{z_{bre_{n}}^2}{2}
   + \frac{z_{bre_{n}}^{3}}{6} \right) W_{end_{n}} \\
   W_{fuel_{n}} &\geq n_{eng} {c_{T_{n}}} D_{n} t_{n} \\
   \sum_{n=1}^{N}W_{fuel_{n}} &\leq W_{f_{primary}} \\
   V_{\infty_{n}} t_{n} &= R_{n} \\
   W_{start_{n}} &\geq W_{end_{n}} + W_{fuel_{n}} \\
   W_{start_{n+1}} &= W_{end_{n}} \\
   W &\geq W_{dry} + W_{payload} + f_{fuel_{res}} W_{f_{primary}} \\
   W_{start_{0}} &= W \\
   W_{end_{N}} &\geq W_{dry} + W_{payload} + f_{fuel_{res}} W_{f_{primary}}\\
   W_{avg_{n}} &\geq \sqrt{W_{start_{n}} W_{end_{N}}} + W_{buoy_{n}} \\
   \left(\frac{L}{D}\right)_{n} &= \frac{W_{avg_{n}}}{D_{n}}
   \end{aligned}

In the remainder of this manuscript, :math:`W` refers to the
corresponding flight segment’s :math:`W_{avg}`.

The dry weight and drag of the aircraft are constrained using simple
buildups of each component’s weight and drag.

.. math::

   \begin{aligned}
   W_{dry} &\geq W_{wing} + W_{fuse} + W_{vt} + W_{ht} + W_{lg} + W_{eng} + W_{misc} \\
   D_n &\geq D_{wing_n} + D_{fuse_n} + D_{vt_n} + D_{ht_n}\end{aligned}

Mach number is constrained to be greater than a user-specified minimum
value.

.. math::

   \begin{aligned}
   M &= \frac{V_{\infty}}{a} \\
   M &\geq M_{min}\end{aligned}

The takeoff model is taken directly
from [Hoburg, 2013]. An additional constraint
on takeoff velocity is added to ensure adequate margin above stall
speed [Anderson, 2001].

.. math::

   \begin{aligned}
   {x_{TO}} &\leq {l_r} \\
   1 + {y} &\leq  2\frac{ {g} {x_{TO}}{T_e}}{{V_{TO}}^{2} {W}}  \\
   1 &\geq  0.0464\frac{{\xi}^{2.7}}{{y}^{2.9}} + \frac{{\xi}^{0.3}}{{y}^{0.049}}\\
   {\xi} &\geq \frac12 \frac{{\rho_{TO}}{V_{TO}}^{2} {S_w}{C_D}}{{T_e}} \\
   {V_{TO}} &= 1.2\sqrt{\frac{2{W}}{C_{L_{w,max}}} {S_w} {\rho_{TO}}} \end{aligned}

Atmospheric pressure, density, temperature, and speed of sound are
constrained using the atmosphere model described in
[York, 2017]. Dynamic viscosity is constrained using
the viscosity model developed in [Kirschen, 2016].

System-level Properties
~~~~~~~~~~~~~~~~~~~~~~~

The constraint for the aircraft is -compatible, and is satisfied during
each flight segment. The fuselage and payload weights are assumed to be
evenly distributed through the length of the fuselage, and the wing
weight acts directly at its area centroid, :math:`x_{wing} + \Delta
x_{ac_w}`. It is assumed that the fuel weight shifts in proportion to
the remaining fuel fraction, :math:`f_{fuel}`, and that a reserve fuel
fraction, :math:`f_{fuel_{res}}`, remains in the wing. The wingbox
forward bulkhead location, :math:`x_b`, is used as a surrogate variable
for engine .

.. math::

   \begin{aligned}
   W x_{CG_{n}} &\geq W_{wing} \left(x_{wing} + \Delta x_{ac_w}\right) 
    + W_{f_{primary}} \left(f_{fuel_{n}} + f_{fuel_{res}}\right) \left(x_{wing} +
    \Delta x_{ac_w} f_{fuel_{n}}\right)  \\
   & +\frac{1}{2} \left(W_{fuse} + W_{payload}\right) l_{fuse}
   + W_{ht} x_{CG_{ht}} + \left(W_{vt} + W_{cone} \right) x_{CG_{vt}} \nonumber \\
   & + n_{eng} W_{eng} x_b + W_{lg} x_{lg} + W_{misc} x_{misc} \nonumber\end{aligned}

In the prior constraint, :math:`f_{fuel}` is the percent of primary fuel
remaining. :math:`f_{fuel}` is represented adequately by a posynomial
inequality since it has downward pressure.

.. math:: f_{fuel_{n}} \geq \frac{\sum_{n=1}^{n}W_{fuel_{n}}}{W_{f_{primary}}}

The landing gear is constrained by the moment of each set of landing
gear about the nose of the aircraft.

.. math:: W_{lg} x_{lg} \geq W_{mg} x_m + W_{ng} x_n

The miscellaneous equipment includes only power systems in the current
model, but is defined to allow for refinements in CG modeling in future
work.

.. math::

   \begin{aligned}
   W_{misc} x_{misc} &\geq W_{hpesys} x_{hpesys}\end{aligned}

The aircraft’s moment of inertia is the sum of the inertias of its
components.

.. math::

   \label{e:Iz_sum}
   I_z \geq I_{z_{wing}} + I_{z_{fuse}} + I_{z_{tail}}

The wing moment of inertia model includes the moment of inertia of the
fuel systems and engines. It assumes that the wing and fuel weight are
evenly distributed on the planform of the wing. This is an overestimate
of the wing moment of inertia with full fuel tanks.

.. math::

   \label{e:Iz_wing}
   I_{z_{wing}} \geq \frac{n_{eng} W_{engine} y_{eng}^2}{g} + 
   \left(\frac{W_{fuel_{wing}} + W_{wing}}{g}\right) \frac{{b_{w}}^3 c_{root_{w}}}{16 S_{w}} 
   \left(\lambda_w + \frac{1}{3}\right)

The fuselage moment of inertia includes the payload moment of inertia.
It is assumed that payload and fuselage weight are evenly distributed
along the length of the fuselage. The wing root quarter-chord location
acts as a surrogate for the of the aircraft.

.. math::

   I_{z_{fuse}} \geq \left(\frac{W_{fuse} + W_{pay}}{g}\right)
   \left(\frac{x_{wing}^3 + l_{vt}^3}{3l_{fuse}}\right)

The moment of inertia of the tail is constrained by treating the tail as
a point mass.

.. math::

   \label{e:Iz_tail}
   I_{z_{tail}} \geq \left(\frac{W_{apu} + W_{tail}}{g}\right) l_{vt}^2
