Engine Model
************



The **Engine** Model has been included in the `turbofan repository`_,
and the constraints within it are detailed in `Martin York's Master's thesis`_.

.. _turbofan repository: https://github.com/convexengineering/turbofan
.. _Martin York's Master's thesis: https://convex.mit.edu/publications/turbofanSP.pdf

Currently both of the **Engine** and **Aircraft** models take **FlightState** as their arguments. To further aid
in integrating the **Engine** model with new or existing configurations, here is a non-exhaustive list of ways
that the engine can and should be integrated into a complete aircraft model.

- Integrate the engine weight into the dry weight of the aircraft.
- Link engine thrust to aircraft drag, climb rate, and fuel burn.
- Link engine sizing variables to the landing gear model.
- Link engine fan area to vertical tail engine out sizing constraint.
- Link engine fan and LPC diameter to the nacelle drag constraints.
- Link engine weight to nacelle weight.
- Link mach numbers (M2) in the engine model to the flight segment Mach number.
- Add engine to aircraft CG model.


