Horizontal Tail Model
=====================

At a conceptual design level, the purpose of the horizontal tail is
threefold: to trim the aircraft such that it can fly in steady level
flight, to provide longitudinal stability, and to give the pilot pitch
control authority over a range of flight conditions.

Model Assumptions
-----------------

The horizontal tail model assumes that the horizontal stabilizer is
mounted to the fuselage and nominally produces downforce in cruise.

Model Description
-----------------

Variable tables are available for download below:

* :download:`Free variables <tables/ht_freevars.pdf>`

* :download:`Fixed variables <tables/ht_fixedvars.pdf>`

Horizontal Tail Geometry and Structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The horizontal tail model employs many of the same geometric constraints
as the wing and vertical tail. More specifically, analogous versions of
Constraints [eq:planformarea,eq:meanaerochord,eq:spanwisemac,eq:taperratio,eq:mintaperratio]
and Constraints [eq:vtmomentarm,eq:vtleading,eq:vttrailing] enforce
planform relationships and constrain the horizontal tail moment arm,
respectively. As with the vertical tail, Constraint  needs to be
implemented as a signomial equality constraint. The horizontal tail also
reuses the same structural model
from [Hoburg, 2013].

Trim Condition
~~~~~~~~~~~~~~

The first sizing requirement is that the aircraft must satisfy the trim
condition [Burton, 2017], which implicitly requires
that the full aircraft moment coefficient be zero.

.. math::

   \frac{x_w}{\bar{c}_w} \leq \frac{x_{CG}}{\bar{c}_w} + \frac{C_{m_{ac}}}{C_{L_w}} 
   + \frac{V_{ht} C_{L_{ht}}}{C_{L_w}}

Thin airfoil theory is used to constrain the horizontal tail’s isolated
lift curve slope [Anderson, 2001].

.. math::

   \begin{aligned}
   C_{L_{ht}} &= C_{L_{\alpha,ht}} \alpha\end{aligned}

However, the horizontal tail’s lift curve slope is reduced by downwash,
:math:`\epsilon`, from the wing and
fuselage [Kroo, 2001]. Note
:math:`\eta_{h_{lift}}` is the horizontal tail sectional lift
efficiency.

.. math::

   C_{L_{\alpha,ht}} = C_{L_{\alpha,ht_0}} \left(1 - \frac{\partial \epsilon}
   {\partial \alpha}\right) \eta_{h_{lift}}

The downwash can be approximated as the downwash far behind an
elliptically loaded wing.

.. math::

   \begin{aligned}
   \epsilon &\approx \frac{2 C_{L_w}}{\pi AR_w} \\
   \implies \frac{\partial \epsilon}{\partial \alpha} &\approx
   \frac{2 C_{L_{\alpha,w}}}{\pi AR_w}\end{aligned}

Thus, an additional posynomial constraint is introduced to constrain
the corrected lift curve slope.

.. math::

   C_{L_{\alpha,ht}} + \frac{2 C_{L_{\alpha,w}} }{\pi AR_w}  \eta_{ht} C_{L_{\alpha,ht_0}}
   \leq C_{L_{\alpha,ht_0}} \eta_{ht}

Minimum Stability Margin
~~~~~~~~~~~~~~~~~~~~~~~~

The second condition is that the aircraft must maintain a minimum
stability margin at both the forward and aft
limits[Burton, 2017].

.. math::

   \begin{aligned}
   \label{e:SM_CG}
   S.M._{min} + \frac{\Delta x_{CG}}{\bar{c}_w} + \frac{C_{m_{ac}}}{C_{L_{w,max}}} 
   &\leq V_{ht} m_{ratio} + \frac{V_{ht} C_{L_{h,max}}}{C_{L_{w,max}}}\end{aligned}

The ratio of the horizontal tail and wing lift curve slopes,
:math:`m_{ratio}`, appears in Equation and is constrained using the
relationship in [Burton, 2017]. The constraint is a
signomial equality because it is not possible to know a priori whether
there will be upward or downward pressure on :math:`m_{ratio}`.

.. math:: m_{ratio} = \left(1+\frac{2}{AR_w}\right) 1 + \frac{2}{AR_{ht}}

Stability Margin
~~~~~~~~~~~~~~~~

The third condition is that the stability margin must be greater than a
minimum specified value for all intermediate locations.

.. math::

   \begin{aligned}
   S.M. &\leq \frac{x_w - x_{CG}}{\bar{c}_w}\\
   S.M. &\geq S.M._{min}\end{aligned}

Horizontal Tail Drag
~~~~~~~~~~~~~~~~~~~~

The horizontal tail employs the same drag model as the wing
(Constraints [eq:wingdrag,eq:wingdragcoeff,eq:wingpdragcoeff,eq:wingRe,eq:induceddrag]),
with the exception of the parasitic drag coefficient fit. The wing’s
parasitic drag fit  is replaced by a fit to XFOIL
[Drela, 1989] data for the
TASOPT[Drela, 2011] T-series airfoils. The TASOPT
T-series airfoils are horizontal tail airfoils intended for transonic
use. The fit considers airfoil thickness, Reynolds number, and Mach
number. The softmax affine function fit is developed with
`GPfit <https://github.com/convexengineering/gpfit>`_ and has an RMS error of 1.14%.

.. math::

   \begin{aligned}
   \label{e:HT_drag}
       {C_{D_{0_{ht}}}}^{6.49} & \geq  5.288\times10^{-20} (Re_{h})^{0.901}  
       (\tau_{h})^{0.912} (M)^{8.645}\\
       &+ 1.676\times10^{-28} (Re_{h})^{0.351} (\tau_{h})^{6.292}
       (M)^{10.256} \nonumber \\
       &+ 7.098\times10^{-25} (Re_{h})^{1.395} (\tau_{h})^{1.962} 
       (M)^{0.567} \nonumber \\
       &+ 3.731\times10^{-14} (Re_{h})^{-2.574} (\tau_{h})^{3.128} 
       (M)^{0.448} \nonumber \\
       &+ 1.443\times10^{-12} (Re_{h})^{-3.910} (\tau_{h})^{4.663} 
       (M)^{7.689} \nonumber \end{aligned}
