Fuselage Model
**************

At a high level, the purpose of a conventional commercial aircraft fuselage can
be decomposed into two primary functions: integrating and connecting all of the
subsystems (e.g. wing, tail, landing gear), and carrying the payload, which
typically consists of passengers, luggage, and sometimes cargo. The design of
the fuselage is therefore coupled with virtually every aircraft subsystem.

Prof. Mark Drela performs a detailed, but still approximate,
analysis of fuselage structure and weight in TASOPT, considering pressure loads, torsion
loads, bending loads, buoyancy weight, window weight, payload-proportional
weights, the floor, and the tail cone. The majority of the constraints in this
model are adapted directly from these equations.

Model Assumptions
=================
This model assumes a single circular-cross-section fuselage. This is an approximation,
since narrow-body aircraft like the Boeing 737 and Airbus A320 do not have perfectly
circular cross sections.

The floor structural model and the horizontal bending model assume uniform
floor loading. The model leverages the analytical bending models from TASOPT,
which makes assumptions about symmetry in bending loads.
Shell buckling is not explicitly modeled while designing bending
structure, but is accounted for by the implementation of a lower yield stress
for bending reinforcement material relative to the nominal yield stress of the
material.

Model Description
=================

.. figure:: fuselage.pdf_tex
    :align: center
    :width: 300px
    Geometric variables of the fuselage model

.. figure:: fuse_cross_section.pdf_tex
    :align: center
    :width: 300px
    Geometric variables of the fuselage model cross-section

Under construction...
