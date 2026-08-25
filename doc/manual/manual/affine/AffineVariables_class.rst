.. _sec-affine-variables-class:

The AffineVariables class
============================

  Main author: `Jordan Ninin <https://www.ensta-bretagne.fr/ninin/>`_

``AffineVariables`` (an alias for ``AffineVarMainVector<AF_Default>``) is a
dynamic-size vector of *declared* affine variables. It is the **only** way
to introduce new, independent noise symbols :math:`\varepsilon_i` into a
computation — a plain :ref:`Affine <sec-affine-class>` built from an
``Interval`` never carries one.

Each component of an ``AffineVariables`` vector belongs to the same
noise-symbol context and has its own dedicated symbol. This identity is
preserved when assigning a new ``Interval`` or ``double`` to an existing
component, and when copying the whole vector. Operations that rebuild the
vector, such as ``resize()``, ``conservativeResize()`` or assignment from an
``IntervalVector``, intentionally create a new context.

The role of the class can be summarized as:

.. code-block:: text

   Interval / Vector
          |
          | declaration of uncertain quantities
          v
   AffineVariables
      x[0] <-> eps_0
      x[1] <-> eps_1
      x[2] <-> eps_2
          |
          | affine / nonlinear computations
          v
   Affine / AffineVector / AffineMatrix

``AffineVariables`` is therefore a declaration container, not a general
purpose container for computed affine expressions.


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
     - Creates ``n`` unbound components, each initially enclosing
       :math:`]-\infty,+\infty[`. Every component receives a distinct
       noise symbol.
   * - ``AffineVariables(const IntervalVector& x)``
     - Creates one component per entry of ``x``. Each component is bound
       to its own noise symbol and has the corresponding interval
       enclosure.
   * - ``AffineVariables(const Vector& x)``
     - Creates one component per entry of ``x``. Each component has a
       degenerate interval enclosure and still receives its own noise
       symbol. The point value therefore does not make the component a
       constant in the dependency model.


Terminology used in this page
------------------------------

A component of ``AffineVariables`` is a **declared affine variable**. Each
such variable is associated with a dedicated **noise symbol**
:math:`\varepsilon_i`. The noise symbol, rather than the component object
itself, is the entity that is shared by the affine expressions derived from
the variable.

This terminology is used consistently throughout the affine documentation:
``AffineVariables`` declares variables and introduces noise symbols, while
``Affine`` / ``AffineVector`` / ``AffineMatrix`` carry those symbols through
computations.


A complete workflow
-------------------

The following example illustrates the typical workflow with
``AffineVariables``:

1. declare the uncertain quantities of the problem;
2. inspect their noise-symbol identities;
3. assign new interval values when needed;
4. use the declared variables in affine expressions;
5. preserve dependencies by reusing the same declared components;
6. materialize the variables as an ``AffineVector`` for vector operations.

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-variables-complete-beg]
      :end-before: [affine-variables-complete-end]
      :dedent: 4

The important point is that the declaration step and the computation step
are separate. ``AffineVariables`` introduces the independent noise symbols;
the resulting ``Affine`` and ``AffineVector`` objects carry those symbols
through subsequent computations.


Component identity
----------------------

Each component exposes its own noise-symbol index through
``noise_index()``. Reassigning a component from an ``Interval`` or a
``double`` replaces its numeric value but keeps that index unchanged.
This lets later computations recognize that two expressions still refer
to the *same* uncertain quantity:

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-variables-2-beg]
      :end-before: [affine-variables-2-end]
      :dedent: 4

Assigning a new interval is therefore **not** equivalent to declaring a new
variable: the existing noise-symbol identity is preserved.

Shared dependency versus independent affine forms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Two expressions derived from the same ``AffineVariables`` component share
the same noise symbol. Consequently, subtracting the same declared variable
from itself is exactly zero:

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-variables-dependency-beg]
      :end-before: [affine-variables-dependency-end]
      :dedent: 4

By contrast, two plain ``Affine`` objects constructed independently from the
same interval do not share a noise symbol. They therefore do not have the
same dependency information.

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

The distinction is important:

.. list-table::
   :widths: 30 30 40
   :header-rows: 1

   * - Method
     - Interval enclosure
     - Noise-symbol identity
   * - ``resize(n)``
     - Existing enclosures are discarded.
     - All components are rebuilt as fresh independent variables.
   * - ``conservativeResize(n)``
     - Surviving components keep their interval enclosures; new components
       are unbounded.
     - Surviving components are reconstructed in a new context, so their
       previous affine dependencies are discarded.

For example, after

.. code-block:: cpp

   AffineVariables v(IntervalVector({{1,2},{-1,1}}));
   v.conservativeResize(3);

the first two components still enclose ``[1,2]`` and ``[-1,1]``, while the
third component is unbounded. The original noise-symbol identities are not
preserved by the resize operation.

In contrast, assigning an ``Interval`` to an existing component or using
``init(const Interval&)`` does not rebuild the context: the component keeps
its existing noise-symbol identity.


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

For example, after

.. code-block:: cpp

   AffineVariables v(3);
   v.init(Interval(1,2));

the three components all enclose ``[1,2]`` but remain independent:
conceptually,

.. code-block:: text

   v[0] = 1.5 + 0.5 eps_0
   v[1] = 1.5 + 0.5 eps_1
   v[2] = 1.5 + 0.5 eps_2

Thus ``v[0] - v[1]`` is not zero, whereas ``v[0] - v[0]`` is exactly
zero.


Assigning values to existing variables
----------------------------------------

An existing component can be assigned from an ``Interval`` or a ``double``.
The value changes, but the component keeps its existing noise-symbol identity:

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-variables-assignment-beg]
      :end-before: [affine-variables-assignment-end]
      :dedent: 4

This is different from creating a new ``Affine`` from an ``Interval``:
the latter does not introduce a symbol that can be shared with a declared
variable.


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
and are therefore seen as fully dependent by any later computation mixing
both:

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-variables-6-beg]
      :end-before: [affine-variables-6-end]
      :dedent: 4

.. warning::

  Copying an ``AffineVariables`` object does **not** create independent
  uncertainty. The copied components refer to the same noise symbols as
  the source components. Use a newly constructed ``AffineVariables`` if
  independent variables are required.


Reassigning from an IntervalVector
----------------------------------------

``operator=(const IntervalVector&)`` reassigns every component from a
fresh interval vector, resizing this vector first if the sizes differ.
Unlike ``init(const Interval&)`` (same size, same symbols, new common
value), this rebuilds every component from scratch — behavior-wise it is
equivalent to ``resize()`` followed by a component-wise assignment, not to
``conservativeResize()``. Consequently, no dependency from before the
assignment survives:

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-variables-7-beg]
      :end-before: [affine-variables-7-end]
      :dedent: 4


What ``AffineVariables`` is — and is not
------------------------------------------

``AffineVariables`` is a declaration container for the independent
uncertain quantities of a problem. It should not be used as a general
purpose vector of computed affine expressions.

Use:

* ``AffineVariables`` to declare independent uncertain quantities;
* :ref:`Affine <sec-affine-class>` to build scalar affine expressions;
* :ref:`AffineVector <sec-affine-vector-class>` to store and manipulate
  vectors of affine expressions;
* :ref:`AffineMatrix <sec-affine-matrix-class>` to store and manipulate
  matrices of affine expressions.

A typical workflow is therefore:

.. code-block:: cpp

   AffineVariables x(...);

   Affine f = x[0] * x[1] + sin(x[2]);

   AffineVector y = A * x;

Once a value is a derived expression rather than a declared uncertain
quantity, it belongs in ``Affine`` or ``AffineVector`` rather than in
``AffineVariables``.


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
carries ``Affine`` (not ``AffineVarMain``) coefficients. This is the natural
transition from declaring uncertain quantities to performing vector and
matrix computations.

See :ref:`the AffineMatrix page <sec-affine-matrix-class>` for more on
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
