.. _sec-domains-slicedtube:

The SlicedTube class
====================

  Main author: `Simon Rohou <https://www.simon-rohou.fr/research/>`_

A :class:`~codac.SlicedTube` is a tube defined over a shared :class:`~codac.TDomain`.
For each temporal slice stored in the temporal partition, the tube owns one :class:`~codac.Slice` object storing its codomain over that temporal support. The codomain type is typically :class:`~codac.Interval` or :class:`~codac.IntervalVector`, or any domain type defined by the user. The same temporal domain can be shared by many tubes.


Creating sliced tubes
---------------------

In Codac C++, :class:`~codac.SlicedTube` is a generic class template ``SlicedTube<T>``. The codomain type ``T`` is not restricted to intervals or boxes: in principle, any domain type compatible with the sliced-tube API can be used.

In practice, the most common sliced-tube types are:

* ``SlicedTube<Interval>``,
* ``SlicedTube<IntervalVector>``,
* ``SlicedTube<IntervalMatrix>``.

These standard tube types are also the ones exposed in the Python and Matlab bindings.

In C++, template argument deduction is available for the constructors based on the arguments. In the examples below, comments make the deduced type explicit when useful.

A sliced tube can be created from:

* a constant codomain,
* an analytic function of one scalar variable,
* a sampled or an analytic trajectory,
* or another sliced tube (copy constructor).

.. tabs::

  .. group-tab:: Python

    .. code-block:: py

      td = create_tdomain([0,10], 0.1, True)

      # Sliced tube from interval-type codomains
      x = SlicedTube(td, Interval(-1,1))
      y = SlicedTube(td, IntervalVector(2))
      z = SlicedTube(td, IntervalMatrix(2,2))

      # From an analytic function
      t = ScalarVar()
      f = AnalyticFunction([t], sin(t))
      xf = SlicedTube(td, f)

      # From an analytic trajectory or sampled trajectory
      at = AnalyticTraj(AnalyticFunction([t],cos(t)), [0,10])
      st = AnalyticTraj(AnalyticFunction([t],cos(t)+t/10), [0,10]).sampled(1e-2)
      xu = SlicedTube(td,at) | SlicedTube(td,st) # union (hull) of tubes

      fig.plot_tube(xf, [Color.dark_blue(),Color.light_gray()])
      fig.plot_tube(xu, [Color.dark_blue(),Color.light_gray()])

  .. group-tab:: C++

    .. code-block:: c++

      auto td = create_tdomain({0,10}, 0.1, true);

      // Sliced tube from interval-type codomains
      SlicedTube x(td, Interval(-1,1)); // x has type SlicedTube<Interval>
      SlicedTube y(td, IntervalVector(2)); // y has type SlicedTube<IntervalVector>
      SlicedTube z(td, IntervalMatrix(2,2)); // z has type SlicedTube<IntervalMatrix>

      // Explicit template notation remains possible:
      SlicedTube<Interval> x2(td, Interval(-1,1));

      // From an analytic function
      ScalarVar t;
      AnalyticFunction f({t}, sin(t));
      SlicedTube xf(td, f);
      // xf has type SlicedTube<Interval> because f outputs scalar values

      // From an analytic trajectory or sampled trajectory
      AnalyticTraj at(AnalyticFunction({t},cos(t)), {0,10});
      AnalyticTraj st(AnalyticFunction({t},cos(t)+t/10), {0,10}).sampled(1e-2);
      auto xu = SlicedTube(td,at) | SlicedTube(td,st); // union (hull) of tubes

      fig.plot_tube(xf, {Color::dark_blue(),Color::light_gray()});
      fig.plot_tube(xu, {Color::dark_blue(),Color::light_gray()});


.. figure:: interval_trajs.png
  :width: 70%

  Created tubes :math:`[x_f](\cdot)=[\sin(\cdot),\sin(\cdot)]` and :math:`[x_u](\cdot)=[\cos(\cdot),\cos(\cdot)+\frac{\cdot}{10}]`, made of 100 slices.


Basic properties
----------------

A sliced tube exposes:

* ``size()`` for the codomain dimension,
* ``shape()`` for matrix-like codomain shapes,
* ``tdomain()`` and ``t0_tf()`` from the base classes,
* ``nb_slices()`` for the number of temporal elements,
* ``codomain()`` for the global codomain hull,
* ``volume()`` for the sum of non-gate slice volumes,
* ``is_empty()`` and ``is_unbounded()``.

Accessing slices and iterating
------------------------------

You can access the first and last slices directly:

* ``first_slice()``
* ``last_slice()``

and you can retrieve a slice from a temporal-domain iterator or from a ``TSlice`` pointer.
The class also defines custom iterators so that iterating on a ``SlicedTube`` yields application-level :class:`~codac.Slice` objects rather than raw ``TSlice`` objects.

.. tabs::

  .. group-tab:: Python

    .. code-block:: py

      td = create_tdomain([0,3])
      x = SlicedTube(td, IntervalVector(2))
      x.set([[1,5],[-oo,2]], [0,1])
      x.set([[2,8],[-oo,3]], [1,2])
      x.set([[6,9],[-oo,4]], [2,3])

      s0 = x.first_slice()
      s1 = s0.next_slice()

      for s in x:
        print(s.t0_tf(), s.codomain())

  .. group-tab:: C++

    .. code-block:: c++

      auto td = create_tdomain({0,3});
      SlicedTube x(td, IntervalVector(2));
      x.set({{1,5},{-oo,2}}, {0,1});
      x.set({{2,8},{-oo,3}}, {1,2});
      x.set({{6,9},{-oo,4}}, {2,3});

      auto s0 = x.first_slice();
      auto s1 = s0->next_slice();

      for(const auto& s : x)
        std::cout << s.t0_tf() << " -> " << s.codomain() << std::endl;


Setting values and refining the partition
-----------------------------------------

A sliced tube can be updated globally or locally:

* ``set(codomain)`` sets all slices to the same codomain,
* ``set(codomain, t)`` sets the value at one time instant,
* ``set(codomain, [ta,tb])`` sets all temporal elements intersecting an interval,
* ``set_ith_slice(codomain, i)`` sets one stored slice by index.

The key point is that local assignments may refine the underlying :class:`~codac.TDomain`.
For example, setting a value at one scalar time creates an explicit gate if necessary. Similarly, setting a value over an interval may create new temporal boundaries at the interval endpoints.

.. note::

  Because the temporal domain is shared, such refinements are structural operations. If several tubes share the same ``TDomain``, all of them will observe the new partition.

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


Evaluation
----------

A sliced tube :math:`[x](\cdot)` can be evaluated over a temporal interval :math:`[t]` with ``x(t)``. The implementation walks through all relevant temporal slices, evaluates each local slice, and unions the results. If the query interval is not included in the temporal domain, the result is the unbounded value of the codomain type (for interval tubes: :math:`[-\infty,\infty]`).

When a derivative tube ``v`` is available, the overload ``x(t,v)`` uses per-slice derivative-aware evaluation and unions the resulting enclosures. This is available for ``Interval`` and ``IntervalVector`` codomain types.

.. tabs::

  .. group-tab:: Python

    .. code-block:: py

      td = create_tdomain([0,3], 1.0, False)
      x = SlicedTube(td, Interval())

      x.set([1,5], [0,1])
      x.set([2,8], [1,2])
      x.set([6,9], [2,3])

      x(0.5)      # [1,5]
      x(1.5)      # [2,8]
      x([0,3])    # [1,9]
      x(-1.0)     # [-oo,oo]

      # No explicit gates: boundary values come from adjacent-slice intersections
      x(1.0)      # [2,5]
      x(2.0)      # [6,8]
      x(3.0)      # [6,9]

  .. group-tab:: C++

    .. code-block:: c++

      auto td = create_tdomain({0,3}, 1.0, false);
      SlicedTube x(td, Interval());

      x.set({1,5}, {0,1});
      x.set({2,8}, {1,2});
      x.set({6,9}, {2,3});

      Interval y0 = x(0.5);   // [1,5]
      Interval y1 = x(1.5);   // [2,8]
      Interval y2 = x({0,3}); // [1,9]
      Interval y3 = x(-1.0);  // [-oo,oo]

      // No explicit gates: boundary values come from adjacent-slice intersections
      Interval y4 = x(1.0);   // [2,5]
      Interval y5 = x(2.0);   // [6,8]
      Interval y6 = x(3.0);   // [6,9]


Inversion
---------

The inversion of a sliced tube :math:`[x](\cdot)`, denoted :math:`[x]^{-1}([y])`, is defined by

.. math::

   [x]^{-1}([y])
   =
   \left\{\, t \mid [x](t)\cap [y]\neq\varnothing \,\right\}
   =
   \bigcup_{y\in [y]} \left\{\, t \mid y\in [x](t) \,\right\}.

It is illustrated by the figure below. Intuitively, the inversion returns the set of time values whose image under :math:`[x](\cdot)` intersects the target set :math:`[y]`.

Depending on the overload, the ``invert()`` methods

* ``invert(y, t=[...])``
* ``invert(y, v_t, t=[...])``
* ``invert(y, v, t=[...])``
* ``invert(y, v_t, v, t=[...])``

return either:

* one interval enclosing the union of all preimages,
* or the different connected components of the inversion in a vector of :class:`~codac.Interval` objects.

When a derivative tube ``v`` is provided, Codac uses derivative information slice by slice in order to compute sharper inverse images.

The following example computes the different connected components of the inversion :math:`[x]^{-1}([0,0.2])` over a temporal subdomain, and then draws their projections in red.

.. tabs::

  .. group-tab:: Python

    .. code-block:: py

      v_t = []
      y = Interval(0,0.2)
      x.invert(y,v_t)

      for t in v_t:
        z = cart_prod(t,y)
        DefaultFigure.draw_box(z, Color.red())

  .. group-tab:: C++

    .. code-block:: c++

      std::vector<Interval> v_t;
      Interval y(0,0.2);
      x.invert(y, v_t);

      for(const auto& t : v_t)
      {
        IntervalVector z = cart_prod(t,y);
        DefaultFigure::draw_box(z, Color::red());
      }

.. figure:: inversion_example.png
  :width: 70%

  Example of tube inversion for a given tube :math:`[x](\cdot)`. The red boxes correspond to the connected components of :math:`[x]^{-1}([0,0.2])`.

.. doxygengroup:: codac2_slicedtube_inversion
  :project: codac


Integration and primitives
--------------------------

Reliable integral computations are available on tubes.

.. figure:: tube_integ_inf.png
  :width: 70%
  
  Hatched part depicts :math:`\int_{a}^{b}x^-(\tau)d\tau`, the lower bound of :math:`\int_{a}^{b}[x](\tau)d\tau`.

The computation is reliable because it stands on the tube's slices. The result is an outer approximation of the integral of the tube represented by these slices:

.. figure:: tube_lb_integral_slices.png
  :width: 70%

  Outer approximation of the lower bound of :math:`\int_{a}^{b}[x](\tau)d\tau`.


.. doxygengroup:: codac2_slicedtube_integrals
   :project: codac


Inflation, extraction, and algebraic helpers
--------------------------------------------

A sliced tube can be inflated either by a constant radius or by a time-varying
sampled radius. Inflation is performed in place.


.. tabs::

  .. group-tab:: Python

    .. code-block:: py

      x.inflate(0.2) # constant inflation

      rad = SampledTraj()
      rad.set(0.1, 0.0)
      rad.set(0.3, 1.0)
      rad.set(0.2, 2.0)
      x.inflate(rad) # time-varying inflation radius

  .. group-tab:: C++

    .. code-block:: c++

      x.inflate(0.2); # constant inflation
      
      SampledTraj<double> rad;
      rad.set(0.1, 0.0);
      rad.set(0.3, 1.0);
      rad.set(0.2, 2.0);
      x.inflate(rad); # time-varying inflation radius

For vector-valued tubes, the API also provides convenient extraction operators:

* ``x[i]`` returns the :math:`i`-th scalar component as a
  :class:`~codac.SlicedTube`,
* ``x.subvector(i,j)`` returns a subvector tube.

The extracted tubes keep the same temporal partition as the original one.

.. tabs::

  .. group-tab:: Python

    .. code-block:: py

      # Component and subvector extraction
      x0 = x[0]
      x12 = x.subvector(1,2)

  .. group-tab:: C++

    .. code-block:: c++

      // Component and subvector extraction
      auto x0 = x[0];              // type: SlicedTube<Interval>
      auto x12 = x.subvector(1,2); // type: SlicedTube<IntervalVector>


Involving a tube in an analytic expression
------------------------------------------

The method ``.as_function()`` wraps a sliced tube as an analytic operator so that it can be embedded inside analytic expressions. This is a convenient bridge between the tube API and the analytic-expression API, see :ref:`the section on temporal operators <sec-functions-temporal-operator>` for its use.