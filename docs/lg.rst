Landing Gear Model
==================

The purpose of the landing gear is to support the weight of the aircraft
and allow it to manoeuvre while it is on the ground, including during
taxi, takeoff, and landing. Including the landing gear in aircraft is
important, not only because it typically weighs between three and six
percent of the maximum aircraft takeoff
weight [Chai, 1996], but also because of how
coupled its design is to other subsystems, particularly the fuselage,
wings, and engines. The landing gear geometry is constrained by wing
position, engine clearance, takeoff rotation, and tip-over criteria. In
addition to being able to withstand nominal static and dynamic loads,
the landing gear also needs to be able to absorb touchdown shock loads.
These loads and the required geometry determine the weight of the gear.
Many of the constraints imposed on landing gear design are described in
[Raymer, 1992] and
[Chai, 1996].

Model Assumptions
-----------------

The landing gear model assumes a conventional and retractable tricycle
landing gear configuration for narrowbody commercial aircraft such as a
Boeing 737-800. The nose gear consists of a single strut supported by
two wheels. The main gear consists of two struts mounted in the inboard
section of the wings, each supported by two wheels. The model only takes
one location as an input, i.e. it does not consider travel. It is also
assumed that the main landing gear retracts towards the centerline of
the aircraft, rotating about the x axis.

Model Description
-----------------

Variable tables are available for download below:

* :download:`Free variables <tables/lg_freevars.pdf>`

* :download:`Fixed variables <tables/lg_fixedvars.pdf>`

Landing Gear Position
~~~~~~~~~~~~~~~~~~~~~

The landing gear track and base are defined relative to the x- and
y-coordinates of the nose and main gear.

.. math::

   \begin{aligned}
   {T} &= 2{y_m} \\
   {x_m} &\geq {x_n} + {B}\end{aligned}

The geometric relationships between the x-coordinates of the main gear,
nose gear and the position must be enforced. These relationships are:

.. math::

   \begin{aligned}
   \label{ngdef} {x_n} + {\Delta  x_n} &= {x_{CG}} \\
   \label{mgdef} {x_{CG}} + {\Delta  x_m} &=  {x_m} \end{aligned}

Equations and must be satisfied exactly, meaning the constraints that
enforce them must be tight. As will be shown below, the load through the
nose gear and main gear is proportional to the distance from the to the
main and nose gear respectively. Because there is downward pressure on
these loads - more load generally means heavier landing gear - there is
also downward pressure on the distances :math:`{\Delta x_n}` and
:math:`{\Delta x_m}`. Therefore signomial constraints are used for both
relationships.

.. math::

   \begin{aligned}
   {x_n} + {\Delta x_n} &\geq {x_{CG}} \\
   {x_{CG}} + {\Delta x_m} &\geq {x_m}\end{aligned}

The main gear position in the spanwise (:math:`y`) direction is, on one
side, lower bounded by the length of the gear itself and, on the other
side, upper bounded by the spanwise location of the engines. Both of
these constraints are necessary to allow the landing gear to retract in
the conventional manner for typical narrowbody commercial aircraft.

.. math::

   \begin{aligned}
   {y_m} &\geq {l_m} \\
   {y_m} &\leq {y_{eng}}\end{aligned}

Wing Vertical Position and Engine Clearance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The difference between the lengths of the main gear and nose gear is
constrained by the vertical position of the wing with respect to the
bottom of the fuselage, as well as the spanwise location of the main
gear and the wing dihedral. This relationship is a signomial constraint.

.. math:: {l_n} + z_{wing} + y_m \tan(\gamma) \geq {l_m}

For aircraft with engines mounted under the wing, the length of the
main gear is also constrained by the engine diameter, because the
engines must have sufficient clearance from the ground. A signomial
constraint provides another lower bound on the length of the main gear.

.. math::

   \begin{aligned}
   \label{fan_diameter_constraint}
   {l_m} + (y_{eng} - y_m)\tan(\gamma) &\geq {d_{nacelle}} + {h_{nacelle}}
   \\
   d_{nacelle} &\geq d_{fan} + 2t_{nacelle} \end{aligned}

Takeoff Rotation
~~~~~~~~~~~~~~~~

The aircraft must be able to rotate on its main wheels at takeoff
without striking the tail of the fuselage and, similarly, must be able
to land on its main gear without striking the
tail [Raymer, 1992]. This constrains the
location of the main gear. More specifically, the horizontal distance
between the main gear and the point at which the fuselage sweeps up
towards the tail must be sufficiently small, relative to the length of
the main gear, such that the angle relative to the horizontal from the
main wheels to the upsweep point is greater than the takeoff/landing
angles. The result is a signomial constraint that imposes a lower bound
on the length of the gear and the x-location of the main gear.

.. math::

   \label{xupsweep_constraint}
   \frac{l_m}{{\tan(\theta_{max})}}\geq{x_{up}}-{x_m}

Tip-over Criteria
~~~~~~~~~~~~~~~~~

A longitudinal tip-over criterion requires that the line between the
main gear and the be at least :math:`15^\circ` relative to the vertical
such that the aircraft will not tip back on its tail at a maximum
nose-up attitude [Raymer, 1992]. This puts a
lower bound on the x-location of the main gear, as measured from the
nose of the aircraft. Note that :math:`\tan(\phi)` is a design variable
here, instead of :math:`\phi`, to make the constraint -compatible.

.. math::

   \begin{aligned}
   {x_m} &\geq \left( {l_m} + {z_{CG}} \right) {\tan(\phi)} + {x_{CG}}\\
   {\tan(\phi)} &\geq {\tan(\phi_{min})} \end{aligned}

A lateral tip-over constraint is introduced to ensure that an aircraft
does not tip over in a turn [Chai, 1996]. The
turnover angle is defined as

.. math:: \tan{\psi} = \frac{z_{CG} + l_m}{{\Delta  x_n} \sin{\delta}}

where

.. math:: \tan{\delta} = \frac{y_m}{B}.

Using the relationship

.. math:: \cos\left(\arctan\left(\frac{y_m}{B}\right)\right)=\frac{B}{\sqrt{B^2 + y_m^2}},

this constraint can be rewritten in, not only -compatible, but
-compatible form as

.. math:: 1\geq\frac{(z_{CG}+l_m)^2 ({y_m}^2 + B^2) }{ (\Delta x_n  y_m  \tan(\psi))^2}.

Typically this angle, :math:`\psi`, should be no larger than
:math:`63^\circ` [Raymer, 1992].

.. math:: {\tan(\psi)}  \leq {\tan(\psi_{max})}

Landing Gear Weight
~~~~~~~~~~~~~~~~~~~

The total landing gear system weight is lower bounded by accounting for
the weights of each assembly. An additional weight fraction is used to
account for weight that is proportional to the weight of the
wheels [Currey, 1984].

.. math::

   \begin{aligned}
   {W_{lg}} &\geq {W_{mg}} + {W_{ng}} \\
   {W_{mg}} &\geq {n_{mg}} \left(W_{ms} + {W_{mw}}(1 + f_{add_m}) \right) \\
   {W_{ng}} &\geq {W_{ns}} + {W_{nw}}(1 +  f_{add_n})\end{aligned}

The weight of each strut for both the main and nose struts is lower
bounded by simplistically assuming a thin-walled cylinder with constant
cross sectional area.

.. math::

   \begin{aligned}
   {W_{ms}} &\geq 2 \pi {r_m}{t_m} {l_m}  {\rho_{st}}g\\
   {W_{ns}} &\geq 2 \pi {r_n}{t_n} {l_n}  {\rho_{st}}g\end{aligned}

It is assumed that the strut is sized by compressive yield and, more
stringently, by buckling, again assuming a thin-walled cylinder. This
constrains the area moment of inertia of the strut cross section, which
puts upward pressure on the radius and thickness of the struts. The
buckling constraint assumes that no side force is exerted on the
cylinder, which is perhaps a weak assumption due to forces exerted in
braking, for example, and due to the fact that aircraft do not typically
land with the main gear struts perfectly normal to the runway surface.

.. math::

   \begin{aligned}
   2 \pi{r_m} {t_m}  {\sigma_{y_c}}&\geq \frac{{\lambda_{LG} L_m} {N_s}}{{n_{mg}}}
   \\
   2 \pi {r_n}  {t_n} {\sigma_{y_c}} &\geq ({L_n} + {L_{n_{dyn}}}) {N_s}\\
    {L_m} &\leq \frac{\pi^2{E}{I_m}}{{K}^{2}{l_m}^{2}}\\
   {I_m} &= \pi{r_m}^{3} {t_m} \\ % PK different
    {L_n} &\leq \frac{\pi^2{E}{I_n}}{{K}^{2}{l_n}^{2}}\\
   {I_n} &= \pi{r_n}^{3} {t_n} \end{aligned}

A machining constraint is used to ensure that the strut walls are not
too thin to be fabricated [Chai, 1996].

.. math::

   \begin{aligned}
    \frac{2 r_m}{t_m} &\leq 40 \\
    \frac{2 r_n}{t_n} &\leq 40 \end{aligned}

The wheel weights can be estimated using historical relations from
[Currey, 1984] and [Raymer, 1992], which are,
again, conveniently in monomial form.

.. math::

   \begin{aligned}
   W_{mw} &= n_{wps} W_{wa,m}\\ 
   W_{nw} &= n_{wps} W_{wa,n}\\
   W_{wa,m} &= 1.2 F_{w_m}^{0.609}\\ 
   F_{wm} &= L_{w_m} d_{t_m}\\
   L_{w_m} &= \frac{L_m}{n_{mg} n_{wps}}\\
   W_{wa,n} &= 1.2 F_{w_n}^{0.609}\\
   F_{wn} &= L_{w_n} d_{t_n} \\
   L_{w_n} &= \frac{L_n}{n_{wps}}\\
   d_{t_m} &= 1.63 L_{w_m}^{0.315} \\
   w_{t_m} &= 0.104 L_{w_m}^{0.480} \\
   d_{t_n} &= 0.8 d_{t_m} \\
   w_{t_n} &= 0.8 w_{t_m} \end{aligned}

Main gear tyre size can also be estimated using statistical relations.
The nose gear tyres are assumed to be 80% of the size of the main gear
tyres.

.. math::

   \begin{aligned}
   d_{t_m} &= 1.63 L_{w_m}^{0.315} \\
   w_{t_m} &= 0.104 L_{w_m}^{0.480} \\
   d_{t_n} &= 0.8 d_{t_m} \\
   w_{t_n} &= 0.8 w_{t_m} \end{aligned}

In addition, simple retraction space constraints are used to ensure
that the gear assemblies are not too wide to fit inside the fuselage.

.. math::

   \begin{aligned}
    2 w_{t_m} + 2 r_m &\leq h_{hold} \\
    2 w_{t_n} + 2 r_n &\leq 0.8~[\mathrm{m}] \end{aligned}

Landing Gear Loads
~~~~~~~~~~~~~~~~~~

The maximum static load through the nose and main gear is constrained by
the weight of the aircraft and the relative distances from the to the
main and nose gear, respectively.

.. math::

   \begin{aligned}
   {L_n} &= \frac{{W} {\Delta x_m}}{{B}} \\
   {L_m} &= \frac{{W} {\Delta x_n}}{{B}} \end{aligned}

 For the nose gear, there is an additional dynamic load due to the
braking condition. A typical braking deceleration of
:math:`3 \mathrm{m/s^2}` is
assumed [Raymer, 1992].

.. math:: {L_{n_{dyn}}} \geq 0.31W \frac{{l_m} + {z_{CG}}}{{B}}

The nose gear requires adequate load for satisfactory steering
performance. A typical desirable range is between 5% and 20% of the
total load [Raymer, 1992].

.. math::

   \begin{aligned}
   \frac{{L_n}}{{W}} &\geq 0.05 \\
   \frac{{L_n}}{{W}} &\leq 0.2 \end{aligned}

Shock Absorption
~~~~~~~~~~~~~~~~

Oleo-penumatic shock absorbers are common to landing gear for large
aircraft. Their purpose is to reduce the vertical load on the aircraft
at touchdown, and they are typically sized by a hard landing condition.
The maximum stroke of the shock absorber can be determined by
considering the aircraft’s kinetic energy, and the target maximum
load [Torenbeek, 1982].

.. math::

   \begin{aligned}
   E_{land} &= \frac{W}{2g} w_{ult}^2 \\
   S_{sa} &= \frac{1}{\eta_s} \frac{E_{land}}{L_m \lambda_{LG}}\end{aligned}

As a preliminary model, the oleo size can be estimated using historical
relations that are conveniently in monomial
form [Raymer, 1992]. The length of the main gear
must be greater than the length of the oleo and the radius of the tyres.

.. math::

   \begin{aligned}
   l_{oleo} &= 2.5 S_{sa} \\
   d_{oleo} &= 1.3 \sqrt{\frac{4 \lambda_{LG} L_{m}/n_{mg}}{p_{oleo} \pi}} \\
   l_{m} &\geq l_{oleo} + \frac{d_{t_m}}{2} \end{aligned}
