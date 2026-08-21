.. _sec-affine-variables-class:

The AffineVariables class
============================

  Main author: `Jordan Ninin <https://www.ensta-bretagne.fr/ninin/>`_

``AffineVariables`` (an alias for ``AffineVarMainVector<AF_Default>``) is a
dynamic-size vector of *declared* affine variables. It is the **only** way
to introduce new, independent noise symbols :math:`\varepsilon_i` into a
computation — a plain :ref:`Affine <sec-affine-class>` built from an
``Interval`` never carries one. Every component of an ``AffineVariables``
vector shares the same noise-symbol context, and keeps its own dedicated
symbol for the lifetime of the vector.


Creating affine variables
----------------------------

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-variables-1-beg]
      :end-before: [affine-variables-1-end]
      :dedent: 4

.. list-table:: Constructors
   :widths: 40 60
   :header-rows: 1

   * - Constructor
     - Description
   * - ``AffineVariables(Index n)``
     - ``n`` unbound components, each :math:`]-\infty,+\infty[`
   * - ``AffineVariables(const IntervalVector& x)``
     - one component per entry of ``x``, each bound to its own symbol
   * - ``AffineVariables(const Vector& x)``
     - degenerate (point) components, still each with its own symbol


Component identity
----------------------

Each component exposes its own noise-symbol index through
``noise_index()``. Reassigning a component from an ``Interval`` or a
``double`` replaces its numeric value but keeps that index unchanged —
this is what lets later computations recognize that two expressions still
refer to the *same* uncertain quantity:

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-variables-2-beg]
      :end-before: [affine-variables-2-end]
      :dedent: 4

.. warning::

  Assigning a plain ``Affine`` into a component does not compile
  (``operator=(const Affine&)`` is deleted on purpose): accepting it would
  let a never-declared dependency be fabricated after the fact. Only
  ``Interval`` and ``double`` can be assigned into an existing component.


Linear combinations
------------------------

Addition, subtraction and scaling by a real constant are exact, whether the
combination mixes several distinct declared variables or reuses the same
one several times. The printed form of an ``Affine`` (``mid + c_0 eps_0 +
c_1 eps_1 + ... + err [-1,1]``) makes the effect of each operation directly
visible:

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-variables-12-beg]
      :end-before: [affine-variables-12-end]
      :dedent: 4

Coefficients scale independently of one another:

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-variables-13-beg]
      :end-before: [affine-variables-13-end]
      :dedent: 4

And reusing the same declared variable several times in a linear expression
still cancels out exactly, for any interval it was declared over:

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-variables-14-beg]
      :end-before: [affine-variables-14-end]
      :dedent: 4

.. note::

  The printed values above were computed by hand from the definition of
  addition/subtraction/scaling (which never involve any linearization) and
  are exact for these clean, exactly-representable inputs.


Nonlinear expressions
--------------------------

Applying an elementary function to a declared variable returns a plain
``Affine``, linearized around the variable's current enclosure:

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-variables-15-beg]
      :end-before: [affine-variables-15-end]
      :dedent: 4

.. note::

  ``sqr`` is a special case: it is computed by directly expanding the
  product :math:`x \times x` rather than through a generic linearization,
  so the coefficients shown above were traced exactly from the source
  (:math:`x_0^2 + \tfrac{1}{2}c^2` for the midpoint, :math:`2 x_0 c` for the
  coefficient, :math:`\tfrac{1}{2}c^2` for the remainder, where
  :math:`x_0` and :math:`c` are ``x``'s own midpoint and coefficient).

For transcendental functions, the exact printed coefficients depend on the
Chebyshev/MinRange linearization and its own floating-point error tracking
— accurate enough to trust for a *sound* enclosure, but not reproduced
here by hand:

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-variables-16-beg]
      :end-before: [affine-variables-16-end]
      :dedent: 4

.. warning::

  Only the containment check (``is_superset``) on this last example was
  verified; the exact coefficients shown are a format placeholder, not
  computed values. Run the snippet to see the actual printed output.


Resizing
------------

Two resizing operations are available, with different guarantees:

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-variables-3-beg]
      :end-before: [affine-variables-3-end]
      :dedent: 4

.. list-table::
   :widths: 30 70
   :header-rows: 1

   * - Method
     - Behavior
   * - ``resize(n)``
     - Drops all existing dependency information; every component becomes
       a fresh, independent, unbounded variable.
   * - ``conservativeResize(n)``
     - Keeps each surviving component's *interval enclosure* (not its
       correlations, which cannot survive a change of context size); newly
       added components are unbounded.


Broadcasting a single interval
----------------------------------

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-variables-4-beg]
      :end-before: [affine-variables-4-end]
      :dedent: 4

.. note::

  The generic ``Eigen::Matrix::init(const Scalar&)`` cannot be used on
  ``AffineVariables`` (it would require an implicit
  ``Interval -> AffineVarMain<T>`` conversion, which does not exist).
  ``AffineVariables`` therefore provides its own ``init(const Interval&)``
  overload, applying the same interval to every component while
  preserving each one's own symbol.


Using declared variables in a computation
----------------------------------------------

``AffineVariables`` is meant for *declaring* the uncertain quantities of a
problem, not for accumulating results: compound assignment operators
(``+=``, ``-=``, ``*=``, ``/=``) are deliberately unavailable on it (see
:ref:`Limitations <sec-affine-class>` of the ``Affine`` class). To compute
with the declared variables, convert to a plain
:ref:`AffineVector <sec-affine-vector-class>`:

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-variables-5-beg]
      :end-before: [affine-variables-5-end]
      :dedent: 4

This conversion goes through ``Affine``'s own copy constructor
component-wise, so it is always safe and preserves every noise symbol
exactly.


Whole-vector assignment
----------------------------

Assigning one ``AffineVariables`` into another (``w = v;``) is not covered
by a dedicated vector-level ``operator=``: it falls through to Eigen's
generic, component-wise matrix assignment, which calls
``AffineVarMain::operator=(const AffineVarMain&)`` per component. This
copies each component's noise-symbol identity along with its value — after
the assignment, ``w[i]`` and ``v[i]`` refer to the exact same noise symbol,
and are therefore seen as fully correlated by any later computation mixing
both:

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-variables-6-beg]
      :end-before: [affine-variables-6-end]
      :dedent: 4

.. note::

  This is by design: it lets ``w`` act as a snapshot of ``v`` that later
  code can still recognize as referring to the same underlying uncertain
  quantities, rather than as an unrelated, freshly independent copy.


Reassigning from an IntervalVector
----------------------------------------

``operator=(const IntervalVector&)`` reassigns every component from a
fresh interval vector, resizing this vector first if the sizes differ.
Unlike ``init(const Interval&)`` (same size, same symbols, new common
value), this rebuilds every component from scratch — behavior-wise it is
equivalent to ``resize()`` followed by a component-wise assignment, not to
``conservativeResize()``:

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-variables-7-beg]
      :end-before: [affine-variables-7-end]
      :dedent: 4


What is deliberately unavailable
--------------------------------------

Beyond the ``operator=(const Affine&)`` case above, a number of operations
are explicitly deleted on ``AffineVariables`` — attempting any of them is a
compile-time error, not a runtime surprise:

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-variables-8-beg]
      :end-before: [affine-variables-8-end]
      :dedent: 4

.. list-table::
   :widths: 45 55
   :header-rows: 1

   * - Deleted
     - Why
   * - ``operator+=``, ``operator-=`` (with ``Vector``, ``IntervalVector``,
       ``AffineVector`` or another ``AffineVariables``)
     - Would mutate a declared variable's value in place
   * - ``operator*=``, ``operator/=`` (with ``double``, ``Interval``,
       ``Affine`` or ``AffineVarMain``)
     - Same reason; scaling in place is not a declaration operation
   * - ``operator=(const AffineVector&)``
     - Would hand a component a value that fabricates a dependency it
       never had (mirrors the single-component ``operator=(const Affine&)``
       restriction)
   * - ``init(const Affine&)``, ``init(const AffineVarMain<T>&)``
     - Same reason as ``operator=(const Affine&)`` on a single component


Products with a matrix
---------------------------

``AffineVariables`` can appear directly as an operand of a matrix product:
it is converted internally to an ``AffineVector``, so the result always
carries ``Affine`` (not ``AffineVarMain``) coefficients — see
:ref:`the AffineMatrix page <sec-affine-matrix-class>` for more on
dependency preservation through matrix products.

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-variables-9-beg]
      :end-before: [affine-variables-9-end]
      :dedent: 4


Unary minus
---------------

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-variables-10-beg]
      :end-before: [affine-variables-10-end]
      :dedent: 4

.. note::

  Like the matrix product above, unary minus returns a plain
  ``AffineVector``, not an ``AffineVariables``: negating a declared
  variable no longer represents "the" variable itself, but a derived
  expression.


Printing
------------

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-variables-11-beg]
      :end-before: [affine-variables-11-end]
      :dedent: 4
