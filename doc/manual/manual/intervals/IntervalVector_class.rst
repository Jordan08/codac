.. _sec-intervals-intervalvector-class:

The IntervalVector class
========================

  Main author: `Simon Rohou <https://www.simon-rohou.fr/research/>`_

An :class:`~codac.IntervalVector` represents an axis-aligned *box* (Cartesian product of intervals) in :math:`\mathbb{R}^n`:

.. math::

   [\mathbf{x}] = [x_1]\times [x_2]\times \dots \times [x_n].

In C++, ``IntervalVector`` is an alias for an Eigen dynamic column vector whose scalar type is :ref:`Interval <sec-intervals-class>` (i.e., ``Eigen::Matrix<Interval,-1,1>``). This means you can use the usual Eigen vector API (size/resize/segment/concatenation patterns) while benefiting from Codac’s interval arithmetic and utilities. In Python and Matlab, some of these functions are available.

Creating IntervalVectors
------------------------

A box can be created from:

* a dimension :math:`n` (components initialized with the default interval :math:`[-\infty,\infty]`),
* a dimension :math:`n` and a common interval (a “cube”),
* a list of intervals or bounds pairs,
* a real vector (degenerate box wrapping a point).

.. tabs::

   .. group-tab:: Python

      .. code-block:: py

         # Default box: [-oo,oo]^n
         x = IntervalVector(3)

         # Cube [-1,3]^2
         y = IntervalVector(2, Interval(-1,3))

         # From a list of bounds (each entry is [lb,ub])
         z = IntervalVector([[3,4],[4,6]]) # [3,4]×[4,6]

         # From a list of components (Intervals and/or bounds pairs)
         q = IntervalVector([y[1], z[0], [0,oo]]) # [-1,3]×[3,4]×[0,oo]

         # From a point (degenerate intervals)
         p = Vector([0.42,0.42,0.42])
         bp = IntervalVector(p) # [0.42,0.42]^3

   .. group-tab:: C++

      .. code-block:: c++

         // Default box: [-oo,oo]^n (Interval default constructor)
         IntervalVector x(3);

         // Cube [-1,3]^2
         IntervalVector y(2, Interval(-1,3));

         // From a list of bounds (initializer-list style)
         IntervalVector z{{3,4},{4,6}}; // [3,4]×[4,6]

         // From a point (degenerate intervals)
         Vector p({0.42,0.42,0.42});
         IntervalVector bp(p); // [0.42,0.42]^3

.. note::

   Codac intentionally does **not** enable Eigen’s “initialize matrices by zero” option for interval-based vectors and matrices, because default-initializing an interval should produce :math:`[-\infty,\infty]` rather than ``0``.
   As a result, resizing an ``IntervalVector`` typically creates new components initialized as :math:`[-\infty,\infty]`.

Accessing and modifying components
----------------------------------

Components are intervals; indexing is 0-based in Python/C++ and 1-based in Matlab.

.. tabs::

   .. group-tab:: Python

      .. code-block:: py

         x = IntervalVector(2, [-1,3])     # [-1,3]^2
         x[1] = Interval(0,10)             # [-1,3]×[0,10]

         # Iterating/accessing over components
         y = IntervalVector(2)
         for i, xi in enumerate(x):
            y[i] = xi

         # Unpacking (Python convenience)
         a,b = x
         assert a == x[0] and b == x[1]

         # Building a new box from existing components
         v = IntervalVector([*x, [3,6]])   # concatenation in Python
         # v == [[-1,3]×[0,10]×[3,6]]

         # Resize: new components are default-initialized ([-oo,oo])
         v.resize(4) # v == [[-1,3]×[0,10]×[3,6]×[-oo,oo]]
         s = v.subvector(1,2) # [0,10]×[3,6]

   .. group-tab:: C++

      .. code-block:: c++

        IntervalVector x(2, Interval(-1,3));  // [-1,3]^2
        x[1] = Interval(0, 10);               // [-1,3]×[0,10]

        // Accessing components
        const Interval& x0 = x[0];

        // Resize: new components are default-initialized ([-oo,oo])
        x.resize(4); // x == [-1,3]×[0,10]×[-oo,oo]×[-oo,oo]

        // Subvector / segment extraction
        IntervalVector s = x.subvector(1,2); // [0,10]×[-oo,oo]

Box properties
--------------

Most scalar properties of :ref:`Interval <sec-intervals-class>` have natural extensions to boxes:

* component-wise bounds: :math:`\mathbf{x}^-` and :math:`\mathbf{x}^+`,
* component-wise widths (diameters/radii),
* geometric metrics such as volume or max diameter.

Typical accessors you will use in practice:

.. tabs::

   .. group-tab:: Python

      .. code-block:: py

         x = IntervalVector([[0,2],[-1,3]])

         n = x.size()        # dimension
         # Common box information (component-wise):
         lo = x.lb()         # Vector of lower bounds
         hi = x.ub()         # Vector of upper bounds
         m  = x.mid()        # Vector of midpoints
         d  = x.diam()       # Vector of diameters

   .. group-tab:: C++

      .. code-block:: c++

         IntervalVector x{{0,2},{-1,3}};

         Index n   = x.size();
         Vector lo = x.lb();
         Vector hi = x.ub();
         Vector m  = x.mid();
         Vector d  = x.diam();

.. note::

   The exact set of box-wide metrics available (*e.g.*, ``max_diam()``, ``volume()``, ``argmax_diam()``) is provided by the interval-based Eigen extensions used in Codac. See the technical documentation section at the end of this page for the comprehensive list.

Testing boxes
-------------

Common predicates include:

- ``is_empty()``: tests whether the box represents the empty set (at least one empty component),
- ``is_degenerated()``: all components are degenerate,
- ``is_degenerated()``: at least one component is degenerate,
- ``is_unbounded()``: at least one bound is infinite,
- ``contains(p)``: inclusion test of a point/vector,
- ``intersects(y)``: non-empty intersection test with another box,
- ``is_subset(y)``, ``is_strict_subset(y)``, *etc.*

.. tabs::

   .. group-tab:: Python

      .. code-block:: py

         x = IntervalVector([[0,1],[2,3]])
         y = IntervalVector([[0.5,2],[1,4]])

         assert x.intersects(y)
         assert x.is_subset(y)

   .. group-tab:: C++

      .. code-block:: c++

         IntervalVector x{{0,1},{2,3}};
         IntervalVector y{{0.5,2},{1,4}};

         assert(x.intersects(y));
         assert(x.is_subset(y));

Advanced operations
-------------------

As for intervals, typical advanced operations on boxes include:

.. list-table:: Common advanced methods for a given box :math:`[\mathbf{x}]`
   :widths: 30 70
   :header-rows: 1

   * - Method
     - Description
   * - ``inflate(rad)``
     - Expands the box by ``±rad`` (scalar radius for common inflation, or vector radius for specifying different inflation on each component).
   * - ``bisect([ratio])``
     - Splits the box into two sub-boxes along the widest dimension.
   * - ``rand()``
     - Returns a random sample ``Vector`` inside the box (when non-empty).
   * - ``init()``
     - Re-initializes to :math:`[-\infty,\infty]^n`.

Arithmetic and linear algebra
-----------------------------

Because ``IntervalVector`` is a vector of ``Interval``, standard arithmetic is naturally extended component-wise:

* box–box and box–scalar operations (``+``, ``-``, ``*``, ``/``),
* interactions with real vectors/matrices,
* matrix–vector products involving :class:`~codac.IntervalMatrix`.

.. tabs::

   .. group-tab:: Python

      .. code-block:: py

         x = IntervalVector([[0,1],[2,3]])
         y = IntervalVector([[1,2],[0,1]])

         z1 = x+y
         z2 = 2*x
         z3 = x-1

   .. group-tab:: C++

      .. code-block:: c++

         IntervalVector x{{0,1},{2,3}};
         IntervalVector y{{1,2},{0,1}};

         IntervalVector z1 = x+y;
         IntervalVector z2 = 2.*x;
         IntervalVector z3 = x-1.;


.. admonition:: Technical documentation

  See the `C++ API documentation of the IntervalVector class <../../api/html/codac2___interval_vector_8h.html>`_.