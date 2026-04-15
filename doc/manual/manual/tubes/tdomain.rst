.. _sec-domains-tdomain:

The TDomain class
=================

  Main author: `Simon Rohou <https://www.simon-rohou.fr/research/>`_

A :class:`~codac.TDomain` is the temporal partition shared by one or more sliced tubes.
It stores an ordered list of :class:`~codac.TSlice` objects. Each stored element is either

* a non-degenerate temporal slice :math:`[t_i,t_{i+1}]`, or
* a degenerate temporal slice :math:`[t_i,t_i]`, called a *gate*.

Gates are optional. When they are explicitly represented, the partition may alternate between
``gate / slice / gate / slice / ...``. This is the structure created by
``create_tdomain(t0_tf, dt, True)``.

Although optional, gates are useful when one needs to represent tube values explicitly at specific time instants. This is particularly relevant for guaranteed integration.

Creating temporal domains
-------------------------

A temporal domain can be created in three common ways:

* from the default constructor helper ``create_tdomain()``;
* from a single time interval ``create_tdomain([t0,tf])``;
* from a sampled interval ``create_tdomain([t0,tf], dt, with_gates)``.

The following examples show the expected structures for the sampled cases above, both with and without gates.

.. tabs::

  .. group-tab:: Python

    .. code-block:: py

      td0 = create_tdomain()                  # one unbounded slice [-oo,oo]
      td1 = create_tdomain([0,1])             # one slice [0,1]
      td2 = create_tdomain([0,1], 0.5, False) # [0,0.5],[0.5,1]
      td3 = create_tdomain([0,1], 0.5, True)  # [0],[0,0.5],[0.5],[0.5,1],[1]

  .. group-tab:: C++

    .. code-block:: c++

      auto td0 = create_tdomain();                  // one unbounded slice [-oo,oo]
      auto td1 = create_tdomain({0,1});             // one slice [0,1]
      auto td2 = create_tdomain({0,1}, 0.5, false); // [0,0.5],[0.5,1]
      auto td3 = create_tdomain({0,1}, 0.5, true);  // [0],[0,0.5],[0.5],[0.5,1],[1]

Inspecting the partition
------------------------

A temporal domain exposes:

* its global time interval ``t0_tf()``,
* the number of stored temporal elements ``nb_tslices()``,
* the number of attached tubes ``nb_tubes()``,
* the temporal slice containing a given time ``tslice(t)``,
* and a vector copy of the stored temporal slices ``tslices_vector()``.

.. tabs::

  .. group-tab:: Python

    .. code-block:: py

      td = create_tdomain([0,1], 0.5, True)

      td.t0_tf()      # [0,1]
      td.nb_tslices() # 5
      td.nb_tubes()   # 0 initially
      
      td.tslice(0.0)  # [0]
      td.tslice(0.1)  # [0,0.5]
      td.tslice(0.5)  # [0.5]
      td.tslice(0.6)  # [0.5,1]
      td.tslice(1.0)  # [1]

  .. group-tab:: C++

    .. code-block:: c++

      auto td = create_tdomain({0,1}, 0.5, true); // returns a std::shared_ptr<TDomain>

      Interval dom = td->t0_tf();  // [0,1]
      size_t n = td->nb_tslices(); // 5
      size_t m = td->nb_tubes();   // 0 initially

      auto it = td->tslice(0.0);   // [0]
      auto it = td->tslice(0.1);   // [0,0.5]
      auto it = td->tslice(0.5);   // [0.5]
      auto it = td->tslice(0.6);   // [0.5,1]
      auto it = td->tslice(1.0);   // [1]


Note that ``tslice(t)`` returns the gate when ``t`` matches an explicit gate, and
otherwise returns the unique non-degenerate slice containing ``t``.

Sampling and refinement
-----------------------

The temporal partition can be refined dynamically with:

* ``sample(t, with_gate=False)`` for one time value,
* ``sample([ta,tb], dt, with_gates=False)`` for repeated sampling over an interval.

Sampling can occur:

* inside the current domain, which splits an existing slice;
* outside the current domain, which extends the partition;
* with gates, which inserts a degenerate temporal slice at the sampling time.

.. tabs::

  .. group-tab:: Python

    .. code-block:: py

      td = create_tdomain()
      td.sample(1.0, False) # [-oo,1][1,oo]
      td.sample(10.0, True) # [-oo,1],[1,10],[10],[10,oo]

  .. group-tab:: C++

    .. code-block:: c++

      auto td = create_tdomain(); // returns a std::shared_ptr<TDomain>
      td->sample(1., false); // [-oo,1],[1,oo]
      td->sample(10., true); // [-oo,1],[1,10],[10],[10,oo]

One important implementation detail is that refinement preserves the attached tubes by cloning the
corresponding slice objects onto the new temporal elements. The new slices are then reattached to
the newly inserted ``TSlice`` iterators.

Shared refinement across several tubes
--------------------------------------

One practical consequence of the shared-``TDomain`` design is that refining the time
partition from one tube automatically updates the slice structure seen by all the
other tubes attached to the same temporal domain.
The following code illustrates this effect:

.. tabs::

  .. group-tab:: Python

    .. code-block:: py

      td = create_tdomain([0,2], 1.0, False) # False: without gates

      x = SlicedTube(td, Interval(0,1))
      v = SlicedTube(td, Interval(-1,1))

      print(x)            # outputs [0,2]↦[0,1], 2 slices
      print(v)            # outputs [0,2]↦[0,1], 2 slices

      td.nb_tslices()     # 2: [0,1],[1,2]
      x.set([0.5,1], 1.3) # local update, will refine the partition at t=1.3
      td.nb_tslices()     # now 4: [0,1],[1,1.3],[1.3],[1.3,2]

      print(x)            # outputs [0,2]↦[-1,1], 4 slices
      print(v)            # outputs [0,2]↦[-1,1], 4 slices (v is also impacted by x.set(..))

  .. group-tab:: C++

    .. code-block:: c++

      auto td = create_tdomain({0,2}, 1.0, false); // false: without gates

      SlicedTube x(td, Interval(0,1));
      SlicedTube v(td, Interval(-1,1));

      cout << x << endl;   // outputs [0,2]↦[0,1], 2 slices
      cout << v << endl;   // outputs [0,2]↦[0,1], 2 slices

      td->nb_tslices();    // 2: [0,1],[1,2]
      x.set({0.5,1}, 1.3); // local update, will refine the partition at t=1.3
      td->nb_tslices();    // now 4: [0,1],[1,1.3],[1.3],[1.3,2]

      cout << x << endl;   // outputs [0,2]↦[-1,1], 4 slices
      cout << v << endl;   // outputs [0,2]↦[-1,1], 4 slices (v is also impacted by x.set(..))

.. note::

   This behavior is intentional: it keeps all related tubes on a consistent shared
   temporal support, which greatly simplifies multi-tube operations such as
   arithmetic combinations and contractions.

Gates
-----

Explicit gates are useful when you want to represent endpoint values separately from the interior
codomain of a slice. The current API provides:

* ``all_gates_defined()`` to detect a fully gate-augmented partition,
* ``delete_gates()`` to remove all explicit gates.

Therefore, a :math:`t`-domain built from ``create_tdomain([0,1], 0.5, True)`` contains five temporal elements and is a
typical example of a gate-augmented partition.
In this implementation, a gate is not a separate data type: it is simply a degenerate
:class:`~codac.TSlice`, i.e., an interval :math:`[t,t]`.

Shared temporal structure
-------------------------

A :class:`~codac.TDomain` is shared by all sliced tubes built from it.
As a consequence, structural operations such as ``truncate(...)`` affect all
attached tubes, and the number of attached tubes can be queried with
``nb_tubes()``.

.. tabs::

  .. group-tab:: Python

    .. code-block:: py

      td = create_tdomain([0,3], 1.0, False)
      x = SlicedTube(td, Interval(0,1))

      td.nb_tubes()      # 1
      td.truncate([0.5,2.5])

      print(td.t0_tf())  # [0.5,2.5]
      print(x)           # x now uses the truncated shared partition

  .. group-tab:: C++

    .. code-block:: c++

      auto td = create_tdomain({0,3}, 1.0, false);
      SlicedTube<Interval> x(td, Interval(0,1));

      std::cout << td->nb_tubes() << std::endl; // 1
      td->truncate({0.5,2.5});

      std::cout << td->t0_tf() << std::endl;    // [0.5,2.5]
      std::cout << x << std::endl;              // x now uses the truncated shared partition