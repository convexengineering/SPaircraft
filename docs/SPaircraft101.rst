SPaircraft 101
**************

These instructions will help you run SPaircraft for the first time!

SP_aircraft.py - set of aircraft models
=======================================
''SP_aircraft.py'' contains a set of methods that are used to optimize different aircraft configurations.

The process to create and run one of these models is general. The following definitions are required:
- ''Nclimb'', number of climb segments
- ''Ncruise'', number of cruise segments
- ''Nmission'', number of missions to simulate
- ''objective'', which is a string designating the objective function
- ''aircraft'', which is a string designating the type of aircraft

The following imports are also required:
.. code-block:: python
    import numpy as np

    # GPkit tools
    from gpkit import units, Model
    from gpkit import Variable, Model, units, SignomialsEnabled, SignomialEquality, Vectorize
    from gpkit.constraints.bounded import Bounded as BCS

    # Constant relaxation heuristic for SP solve
    from relaxed_constants import relaxed_constants, post_process

    # Mission model
    from D8 import Mission

For the demonstration, we define these as:
.. code-block:: python
    Nclimb = 3
    Ncruise = 2
    Nmission = 1
    objective = 'fuel' # minimizing fuel burn
    aircraft = 'D82' # optimizing the D82 configuration

Now, we have the basic building blocks of a model:

.. code-block:: python
    m = Mission(Nclimb, Ncruise, objective, aircraft, Nmission)

Since we are optimizing a D8.2 aircraft, we need to import the appropriate parameter substitutions:

.. code-block:: python
    from subsD82 import getD82subs
    substitutions = getD82subs()

And we need to specify a mission range:

.. code-block:: python
    substitutions.update({
        'ReqRng': 3000.*units('nmi'),
    })

We can update our model ''m'' with the substitutions:
.. code-block:: python
    m.substitutions.update(substitutions)
    m = Model(m.cost, BCS(m))

We use Bounded (BCS) in SPaircraft to make sure that the models are stable even if not all variables are well-constrained. U
Users are warned about unbounded variables at the end of the solve.

We implement the relaxed constants heuristic to solve the SP:
.. code-block:: python
    m_relax = relaxed_constants(m, None)

And now we solve, and print the results:
.. code-block:: python
    sol = m_relax.localsolve(verbosity=4, iteration_limit=200, reltol=0.01)
    print sol.table()

Solution visualization
======================

Model hierarchy
===============

Single-mission optimization
===========================

Multi-mission optimization
==========================
