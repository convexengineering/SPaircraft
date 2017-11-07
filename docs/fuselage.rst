Fuselage Model
==============

The purpose of a conventional commercial aircraft
fuselage can be decomposed into two primary functions: integrating and
connecting all of the subsystems (e.g. wing, tail, landing gear), and
carrying the payload, which typically consists of passengers, luggage,
and sometimes cargo. The design of the fuselage is therefore coupled
with virtually every aircraft subsystem.

A detailed but still approximate analysis of fuselage structure and weight
is done in TASOPT [Drela, 2011] , considering
pressure loads, torsion loads, bending loads, buoyancy weight, window
weight, payload-proportional weights, the floor, and the tail cone. The
majority of the constraints in this model are adapted directly from
these equations.



Model Assumptions
-----------------

This model assumes a single circular-cross-section fuselage. This is an
approximation, since narrow-body aircraft like the Boeing 737 and Airbus
A320 do not have perfectly circular cross sections.

The floor structural model and the horizontal bending model assume
uniform floor loading. The model leverages the analytical bending models
from [Drela, 2011], which makes assumptions
about symmetry in bending loads. Shell buckling is not explicitly
modeled while designing bending structure, but is accounted for by the
implementation of a lower yield stress for bending reinforcement
material relative to the nominal yield stress of the material.

Model Description
-----------------

Variable tables are available for download below:

* :download:`Free variables <tables/fu_freevars.pdf>`

* :download:`Fixed variables <tables/fu_fixedvars.pdf>`

Cross-sectional geometry constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The fuselage must be wide enough to accommodate the width of the seats
in a row and the width of the aisle.

.. math:: {2w_{fuse}} \geq (\mathit{SPR}) {w_{seat}} + {w_{aisle}} + 2{w_{sys}}

The cross sectional area of the fuselage skin is lower bounded using a
thin walled cylinder assumption.

.. math:: {A_{skin}} \geq 2 \pi {R_{fuse}} {t_{skin}}

The cross sectional area of the fuselage is lower bounded using the
radius of the fuselage.

.. math:: {A_{fuse}} \geq \pi {R_{fuse}}^{2}

Pressure loading constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The axial and hoop stresses in the fuselage skin are constrained by the
pressurization load due to the difference between cabin pressure and
ambient pressure at cruise altitude. The thickness of the skin is
therefore sized by the maximum allowable stress of the chosen material.

.. math::

   \begin{aligned}
   {\sigma_x} &= \frac{{\Delta P_{over}}}{2}\frac{{R_{fuse}}}{{t_{shell}}}\\
   {\sigma_{\theta}} &= {\Delta P_{over}} \frac{{R_{fuse}} }{{t_{skin}}} \\
   {\sigma_{skin}} &\geq {\sigma_x} \\
   {\sigma_{skin}} &\geq {\sigma_{\theta}}\end{aligned}

Floor loading constraints
~~~~~~~~~~~~~~~~~~~~~~~~~

The floor must be designed to withstand at least the weight of the
payload and seats multiplied by a safety factor for an emergency
landing.

.. math:: {P_{floor}} \geq {N_{land}} ({W_{payload}} + {W_{seat}})

The maximum moment and shear in the floor are determined based on this
design load and the width of the floor, assuming that the floor/wall
joints are pinned and there are no center supports.

.. math::

   \begin{aligned}
   {S_{floor}} &= \frac{P_{floor}}{2}\\
   {M_{floor}} &= \frac{{P_{floor}} {w_{floor}} }{8}\end{aligned}

The floor beam cross sectional area is constrained by the maximum
allowable cap stress and shear web stress for the beams.

.. math::

   {A_{floor}} \geq 1.5\frac{{S_{floor}}}{{\tau_{floor}}}
   + 2\frac{{M_{floor}}}{{\sigma_{floor}} {h_{floor}}}

3-dimensional geometry constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The nose must be long enough to have an aerodynamic profile and to
accommodate the cockpit. A reasonable, but arbitrary, lower bound is
employed for this work [Drela, 2011].

.. math:: {l_{nose}} \geq 5.2 \hspace{0.2cm} \rm{m}

The cylindrical shell of the fuselage sits between the nosecone and
tailcone. The variables :math:`x_{shell1}` and :math:`x_{shell2}` define
the beginning and end of the cylindrical section of the fuselage,
respectively, in the aircraft x-axis.

.. math::

   \begin{aligned}
   {x_{shell1}} &= {l_{nose}} \\
   {x_{shell2}} &\geq {l_{nose}} + {l_{shell}}\end{aligned}

The number of seats is equal to the product of the seats per row and the
number of rows. Note that non-integer numbers of rows are allowed and
necessary for GP compatibility. It is assumed that the load factor is
one, so that the number of passengers is equal to the number of seats.

.. math::

   \begin{aligned}
   {n_{seat}} &= {(\mathit{SPR})} {n_{rows}} \\
   {n_{pass}} &= n_{seat}\end{aligned}

The seat pitch and the number of rows of seats constrain the length of
the shell. The passenger floor length is lower bounded by the shell
length and twice the fuselage radius, to account for the space provided
by pressure bulkheads.

.. math::

   \begin{aligned}
   {l_{shell}} &\geq {n_{rows}} {p_s} \\
   {l_{floor}} &\geq 2{R_{fuse}} + {l_{shell}} \end{aligned}

The length of the fuselage is constrained by the sum of the nose, shell
and tail cone lengths. A signomial equality is needed, because increased
:math:`l_{fuse}` is not coupled directly to increased structural weight
although it results in improved tail control authority.

.. math::

   l_{fuse} = l_{nose} +
   l_{shell} + l_{cone}

Other locations to constrain are the wing mid-chord and the wingbox fore
and aft bulkheads, which serve as integration limits when calculating
bending loads.

.. math::

   \begin{aligned}
   x_f \leq x_{wing} + 0.5 c_0 r_{w/c}\\
   x_b  + 0.5 c_0 r_{w/c} \geq x_{wing}\end{aligned}

The skin surface area, and, in turn, skin volume for the nose, main
cabin, and rear bulkhead are constrained. The surface area of the nose,
which is approximated as an ellipse, is lower bounded using Cantrell’s
approximation [Drela, 2011].

.. math::

   \begin{aligned}
   {S_{nose}}^{\frac85} &\geq \left(2 \pi {R_{fuse}^2}\right)^{\frac85}
       \left(\frac13 + \frac23
       \left(\frac{l_{nose}}{R_{fuse}}\right)^{\frac85} \right) \\
   {S_{bulk}} &= 2 \pi {R_{fuse}}^{2} \\
   {V_{cyl}} &= {A_{skin}} {l_{shell}} \\
   {V_{nose}} &= {S_{nose}} {t_{skin}} \\
   {V_{bulk}} &= {S_{bulk}} {t_{skin}} \end{aligned}

The cabin volume is constrained assuming a cylinder with hemispherical
end caps. This is necessary for capturing buoyancy weight.

.. math::

   {V_{cabin}}\geq{A_{fuse}}\left(\frac23{l_{nose}} + {l_{shell}} +
   \frac23{R_{fuse}} \right)

Tail cone constraints
~~~~~~~~~~~~~~~~~~~~~

The tail cone needs to be able to transfer the loads exerted on the
vertical tail to the rest of the fuselage. The maximum torsion moment
imparted by the vertical tail depends on the maximum force exerted on
the tail as well as its span and taper ratio. This torsion moment, along
with the cone cross sectional area and the maximum shear stress of the
cone material, bounds the necessary cone skin thickness. The cone cross
sectional area, which varies along the cone, is coarsely approximated to
be the fuselage cross sectional area (i.e. the cross sectional area of
the cone base).

.. math::

   \begin{aligned}
   \label{eq:Qv1} {Q_v} &= \frac{{L_{vt_{max}}}
   {b_{vt}}}{3} \frac{{1 + 2{\lambda_v}}} {{1 + {\lambda_v}}} \\
   \label{eq:Qv2}
   {t_{cone}}&= \frac{Q_v}{2{A_{fuse}} {\tau_{cone}} }\end{aligned}

The volume of the cone is a definite integral from the base to the tip
of the cone. This integral is
evaluated [Drela, 2011] and combined with
Equations and to give a single signomial constraint on the cone skin
volume.

.. math::

   R_{fuse}\tau_{cone}(1+p_{\lambda_v})V_{cone} \frac{1+\lambda_{cone}}{4 l_{cone}}
   \geq L_{vt_{max}} b_{vt} \frac{p_{\lambda_v}}{3}

A change of variables is used for compatibility with the tail model,
which uses :math:`p_{\lambda_v} = 1 + 2\lambda_v` to make a structural
constraint -compatible. The same taper lower bound is introduced as in
the tail model.

.. math:: {p_{\lambda_v}} \geq 1.6

The cone skin shear stress is constrained to equal the maximum
allowable stress in the skin material.

.. math:: {\tau_{cone}} = {\sigma_{skin}}

The tail cone taper ratio constrains the length of the cone relative to
the radius of the fuselage.

.. math:: {l_{cone}} = \frac{{R_{fuse}}}{{\lambda_{cone}}}

Fuselage area moment of inertia constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The fuselage shell consists of the skin and stringers. Its area moment
of inertia determines how effectively the fuselage is able to resist
bending loads. A shell with uniform skin thickness and stringer density
has a constant area moment of inertia in both of its bending axes, shown
by the dark red line in the lower plot of
Figure [fig:fuse\_bending\_loads].

To be consistent with [Drela, 2011], the
horizontal bending moments are defined as the moments around the
aircraft’s y-axis, caused by horizontal tail loads and fuselage inertial
loads, and vertical bending moments as the moments around the aircraft’s
z-axis, caused by vertical tail loads.

.. figure:: figs/fuse_bending_loads.png
   :alt: TASOPT fuselage bending models
   (from[Drela, 2011]). The top graph shows the
   bending load distribution on the fuselage, whereas the bottom graph
   shows the area moment of inertia distribution. The horizontal bending
   loads are shown in blue, and the vertical bending loads are shown in
   red.
   :width: 100.0%

   TASOPT fuselage bending models
   (from[Drela, 2011]). The top graph shows the
   bending load distribution on the fuselage, whereas the bottom graph
   shows the area moment of inertia distribution. The horizontal bending
   loads are shown in blue, and the vertical bending loads are shown in
   red.

The effective modulus-weight shell thickness is lower bounded by
assuming that only the skin and stringers contribute to bending. This
constraint also uses an assumed fractional weight of stringers that
scales with the thickness of the skin.

.. math::

   {t_{shell}} \geq {t_{skin}}\left(1 + {f_{string}} {r_E}
   \frac{{\rho_{skin}} }{{\rho_{bend}}} \right)

It is important to consider the effects of pressurization on the yield
strength of the bending material. Since pressurization stresses the
airframe, the actual yield strength of the fuselage bending material is
lower than its nominal yield strength, an effect captured using
posynomial constraints.

.. math::

   \begin{aligned}
   \sigma_{M_h} + r_E \frac{\Delta P_{over} R_{fuse}}{2 t_{shell}}&\leq
       \sigma_{bend} \\
   \sigma_{M_v} + r_E \frac{\Delta P_{over} R_{fuse}}{2 t_{shell}}&\leq
       \sigma_{bend}\end{aligned}

The aircraft shell, which is composed of the pressurized skin and
stringers, must satisfy the following horizontal and vertical area
moment of inertia constraints.

.. math::

   \begin{aligned}
   I_{hshell} &\leq \pi R_{fuse}^3 t_{shell} \\
   I_{vshell} &\leq \pi R_{fuse}^3 t_{shell}\end{aligned}

Horizontal bending model
~~~~~~~~~~~~~~~~~~~~~~~~

There are two load cases that determine the required : maximum load
factor (MLF) at :math:`V_{ne}`, where

.. math::

   \begin{aligned}
   N &= N_{lift} \\
   L_{ht} &= L_{ht_{max}} \end{aligned}

 and emergency landing impact, where

.. math::

   \begin{aligned}
   N &= N_{land} \\
   L_{ht} &= 0. \end{aligned}

Both load cases are considered at the aircraft’s maximum takeoff weight
(MTOW). The constraints for each case are distinguished by the
subscripts :math:`MLF` and :math:`Land`. Assuming the fuselage weight is
uniformly distributed throughout the shell, the bending loads due to
fuselage inertial loads increase quadratically from the ends of the
fuselage shell to the aircraft , as shown by the blue line representing
:math:`M_h(x)` in Figure [fig:fuse\_bending\_loads]. The tail loads are
point loads at :math:`x_{tail}`, and so the horizontal tail moment
increases linearly from :math:`x_{tail}` to the aircraft’s . In the
maximum load factor case, the maximum moment exerted by the horizontal
tail is superimposed on the maximum fuselage inertial moment at load
factor :math:`N_{lift}` to size the required. For the emergency landing
impact case, only the fuselage inertial loads are considered at
:math:`N_{land}`, assuming an unloaded horizontal tail.

Several intermediate variables are introduced and used in constraints
that capture relationships. :math:`A_{0h}` represents the area that is
contributed by the aircraft shell.

.. math:: A_{0h} = \frac{I_{hshell}} {r_{E} h_{fuse}^2}

Variables :math:`A_{1h_{Land}}` and :math:`A_{1h_{MLF}}` are the lengths
that are required to sustain bending loads from the tail. Note that as
the distance from the tail increases, the moment exerted from the tail
increases linearly.

.. math::

   \begin{aligned}
   A_{1h_{Land}} &\geq N_{land} \frac{W_{tail} + W_{apu}}{h_{fuse} \sigma_{M_h}}\\
   A_{1h_{MLF}} &\geq N_{lift} \frac{W_{tail} + W_{apu} + r_{M_h}
   L_{ht_{max}}}{h_{fuse} \sigma_{M_h}}\end{aligned}

Variables :math:`A_{2h_{Land}}` and :math:`A_{2h_{MLF}}` represent the
required to sustain the distributed loads in the fuselage. As the
distance from the nose or the tail increases, the moment exerted due to
the distributed load grows with the square of length.

.. math::

   \begin{aligned}
   A_{2h_{Land}} &\geq N_{land} \frac{W_{payload} + W_{padd} + W_{shell} +
   W_{window} + W_{insul} + W_{floor} + W_{seat}} {2 l_{shell} h_{fuse}
   \sigma_{bend}} \\
   A_{2h_{MLF}} &\geq N_{lift} \frac{W_{payload} + W_{padd} + W_{shell} +
   W_{window} + W_{insul} + W_{floor }+ W_{seat}} {2 l_{shell} h_{fuse}
   \sigma_{M_h}}\end{aligned}

Bending reinforcement material in the aircraft exists where the shell
inertia is insufficient to sustain the local bending moment. Constraints
are used to determine the location over the rear fuselage
:math:`x_{hbend_\zeta}` forward of which additional is required. Some
simple constraints on geometry are added to ensure a meaningful
solution. Constraints through occur for both aforementioned load cases
in the model (with subscript :math:`\zeta` replaced by :math:`MLF` or
:math:`Land`) for worst-case fuselage sizing, but have been included
once in the paper to reduce redundancy.

.. math::

   \begin{aligned}
   \label{eq:dupBend_1} A_{0h} &= A_{2h_\zeta} (x_{shell2} - x_{hbend_\zeta}) ^ 2 +
   A_{1h_\zeta}  (x_{tail} - x_{hbend_\zeta}) \\ x_{hbend_\zeta} &\geq x_{wing}\\ x_{hbend_\zeta}
   &\leq l_{fuse}  \end{aligned}

To be able to constrain the volume of required, the area of required
must be constrained and integrated over the length of the fuselage. As
shown by [Drela, 2011], with some conservative
approximation, the volume of may be determined through the integration
of the forward and rear wingbox areas over the rear fuselage.

.. math::

   \begin{aligned}
   A_{hbendf_\zeta} &\geq A_{2h_\zeta} (x_{shell2} - x_{f})^2 + A_{1h_\zeta}
       (x_{tail} - x_{f}) - A_{0h} \\
   A_{hbendb_\zeta} &\geq A_{2h_\zeta} (x_{shell2} - x_{b})^2 + A_{1h_\zeta}
       (x_{tail} - x_{b}) - A_{0h}\end{aligned}

volumes forward, over and behind the wingbox are lower bounded by the
integration of the areas over the three fuselage sections.

.. math::

   \begin{aligned}
   V_{hbend_{f}} &\geq \frac{A_{2h_\zeta}} {3} ((x_{shell2} - x_{f})^3 -
       (x_{shell2} - x_{hbend_\zeta})^3) \\
   &+ \frac{A_{1h_\zeta}} {2} ((x_{tail} - x_{f})^2 - (x_{tail} -
       x_{hbend_\zeta})^2) - A_{0h} (x_{hbend_\zeta} - x_{f})\nonumber\\
   V_{hbend_{b}} &\geq \frac{A_{2h_\zeta}}{3} ((x_{shell2} - x_{b})^3 -
       (x_{shell2} - x_{hbend_\zeta})^3) \\
   &+ \frac{A_{1h_\zeta}}{2} ((x_{tail} - x_{b})^2 - (x_{tail} -
       x_{hbend_\zeta})^2) - A_{0h} (x_{hbend_\zeta} - x_{b}) \nonumber\\
   V_{hbend_{c}} &\geq 0.5 (A_{hbendf_\zeta} + A_{hbendb_\zeta}) c_{0} r_{w/c} 
   \label{eq:dupBend_2}\end{aligned}

The total volume is lower bounded by the sum of the volumes of required
in each fuselage section.

.. math:: V_{hbend} \geq V_{hbend_{c}} + V_{hbend_{f}} + V_{hbend_{b}}

Vertical bending model
~~~~~~~~~~~~~~~~~~~~~~

The is constrained by considering the maximum tail loads that a fuselage
must sustain. The vertical bending moment, shown in red as
:math:`M_v(x)` in Figure [fig:fuse\_bending\_loads], increases linearly
from the tail to the aircraft , since the tail lift is assumed to be a
point force.

As with horizontal bending, several intermediate variables are
introduced and used in constraints that capture relationships.
:math:`B_{1v}` is the length required to sustain the maximum vertical
tail load :math:`L_{vt_{max}}`. When multiplied by the moment arm of the
tail relative to the fuselage cross-sectional location, it gives the
local area required to sustain the loads.

.. math:: B_{1v} = \frac{r_{M_v} L_{vt_{max}}} {w_{fuse} \sigma_{M_{v}}}

:math:`B_{0v}` is the equivalent area provided by the fuselage shell.

.. math:: {B_{0v}} = \frac{{I_{vshell}}}{{r_E} {w_{fuse}}^{2}}

Since tail loads are the only vertical loads to consider, the location
forward of which additional bending material is required can be
determined. :math:`x_{vbend}` is the location where the vertical bending
moment of the inertia of the fuselage is exactly enough to sustain the
maximum vertical bending loads from the tail, expressed by a signomial
equality.

.. math::

   \begin{aligned}
    B_{0v} &= B_{1v} (x_{tail} - x_{vbend}) \\ x_{vbend}
   &\geq x_{wing} \\ x_{vbend} &\leq l_{fuse}  \end{aligned}

The area required at the rear of the wingbox is lower bounded by the
tail bending moment area minus the shell vertical bending moment area.

.. math:: A_{vbend_{b}} \geq B_{1v} (x_{tail} - x_{b}) - B_{0v}

The vertical bending volume rear of the wingbox is then constrained by
integrating :math:`A_{vbend}` over the rear fuselage, which yields the
following constraint.

.. math::

   V_{vbend_{b}} \geq 0.5 B_{1v} ((x_{tail}-x_{b})^2 - (x_{tail} - x_{vbend})^2) -
   B_{0v} (x_{vbend} - x_{b})

The vertical bending volume over the wingbox is the average of the
bending area required in the front and back of the wingbox. Since no
vertical bending reinforcement is required in the forward fuselage, the
resulting constraint is simply:

.. math:: V_{vbend_{c}} \geq 0.5 A_{vbend_{b}} c_{0} r_{w/c}

The total vertical bending reinforcement volume is the sum of the
volumes over the wingbox and the rear fuselage.

.. math:: V_{vbend} \geq V_{vbend_{b}} + V_{vbend_{c}}

Weight build-up constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The weight of the fuselage skin is the product of the skin volumes
(bulkhead, cylindrical shell, and nosecone) and the skin density.

.. math::

   {W_{skin}} \geq {\rho_{skin}} {g}  \left({V_{bulk}} + {V_{cyl}}  
   + {V_{nose}} \right)

The weight of the fuselage shell is then constrained by accounting for
the weights of the frame, stringers, and other structural components,
all of which are assumed to scale with the weight of the skin.

.. math::

   {W_{shell}} \geq {W_{skin}}\left(1 + {f_{fadd}} +  {f_{frame}}
   +  {f_{string}} \right)

The weight of the floor is lower bounded by the density of the floor
beams multiplied by the floor beam volume, in addition to an assumed
weight/area density for planking.

.. math::

   \begin{aligned}
   {V_{floor}} &\geq {A_{floor}} {w_{floor}} \\
   {W_{floor}}&\geq{V_{floor}}{\rho_{floor}}{g}+{W''_{floor}}{l_{floor}} {w_{floor}}\end{aligned}

As with the shell, the tail cone weight is bounded using assumed
proportional weights for additional structural elements, stringers, and
frames.

.. math::

   {W_{cone}}\geq{\rho_{cone}}{g}{V_{cone}}\left(1+{f_{fadd}}+{f_{frame}} +
   f_{string}\right)

The weight of the horizontal and vertical bending material is the
product of the bending material density and the and volumes required
respectively.

.. math::

   \begin{aligned}
   W_{hbend} &\geq \rho_{bend} g V_{hbend} \\
   W_{vbend} &\geq \rho_{bend} g V_{vbend}\end{aligned}

The weight of luggage is lower bounded by a buildup of 2-checked-bag
customers, 1-checked-bag customers, and average carry-on weight.

.. math::

   {W_{lugg}} \geq 2{W_{checked}} {f_{lugg,2}} {n_{pass}} +
   {W_{checked}} {f_{lugg,1}} {n_{pass}} + {W_{carry on}}

The window and insulation weight are lower bounded using assumed
weight/length and weight/area densities respectively. It is assumed that
only the passenger compartment of the the cabin is insulated and that
the passenger compartment cross sectional area is approximately 55% of
the fuselage cross sectional area.

.. math::

   \begin{aligned}
   {W_{window}} &= {W'_{window}} {l_{shell}} \\
   {W_{insul}} &\geq {W''_{insul}} \left( 0.55\left({S_{bulk}}  
   + {S_{nose}} \right) + 1.1\pi{R_{fuse}} {l_{shell}} \right) \end{aligned}

The APU and other payload proportional weights are accounted for using
weight fractions. :math:`W_{padd}` includes flight attendants, food,
galleys, toilets, furnishing, doors, lighting, air conditioning, and
in-flight entertainment systems. The total seat weight is a product of
the weight per seat and the number of seats.

.. math::

   \begin{aligned}
   {W_{apu}} &= {W_{payload}} {f_{apu}} \\
   {W_{padd}} &= {W_{payload}} {f_{padd}} \\
   {W_{seat}} &= {W'_{seat}} {n_{seat}} \end{aligned}

The effective buoyancy weight of the aircraft is constrained using a
specified cabin pressure :math:`p_{cabin}`, the ideal gas law and the
approximated cabin volume. A conservative approximation for the buoyancy
weight that does not subtract the ambient air density from the cabin air
density is used.

.. math::

   \begin{aligned}
   \rho_{cabin}&= \frac{p_{cabin}}{{R} {T_{cabin}}} \\
   {W_{buoy}} &= \rho_{cabin} {g} {V_{cabin}}\end{aligned}

There are two methods in the model that can be used to lower bound the
payload weight. The first is the sum of the cargo, luggage, and
passenger weights (Constraint ). The second is through the definition of
variable :math:`W_{avg. pass_{total}}`, which is an average payload
weight per passenger metric (Constraint ). For the purposes of this
paper, the second method is used, and as a result Constraint  is
inactive.

.. math::

   \begin{aligned}
   W_{pass} &= W_{avg. pass} n_{pass} \\
   {W_{payload}} &\geq {W_{cargo}} + {W_{lugg}} + {W_{pass}}\label{eq:payload1st} \\
   {W_{payload}} &\geq {W_{avg. pass_{total}}} {{n_{pass}}} \label{eq:payload2nd}\end{aligned}

The total weight of the fuselage is lower bounded by the sum of all of
the constituent weights. The fixed weight :math:`W_{fix}` incorporates
pilots, cockpit windows, cockpit seats, flight instrumentation,
navigation and communication equipment, which are expected to be roughly
the same for all aircraft [Drela, 2011].

.. math::

   \begin{aligned}
   {W_{fuse}} &\geq {W_{apu}} + {W_{buoy}} + {W_{cone}} + {W_{floor}} + W_{hbend}
       + W_{vbend} + {W_{insul}} \\ &+ {W_{padd}} + {W_{seat}} + {W_{shell}} +
       {W_{window}} + {W_{fix}} \nonumber\end{aligned}

Aerodynamic constraints
~~~~~~~~~~~~~~~~~~~~~~~

The drag of the fuselage is constrained using :math:`C_{D_{fuse}}` from
TASOPT, which calculates the drag using a pseudo-axisymmetric
viscous/inviscid calculation, and scaling appropriately by fuselage
dimensions and Mach number.

.. math::

   D_{fuse} = \frac{1}{2} \rho_{\infty} V_{\infty}^2 C_{D_{fuse}} \left( l_{fuse} R_{fuse}
   \frac{M^2}{M_{fuseD}^2} \right)
