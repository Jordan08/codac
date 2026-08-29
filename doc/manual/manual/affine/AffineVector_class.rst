.. _sec-affine-vector-class:

The AffineVector class
==========================

  Main author: `Jordan Ninin <https://www.ensta-bretagne.fr/ninin/>`_

``AffineVector`` (an alias for ``AffineMainVector<AF_Default>``) is a
dynamic-size column vector of :ref:`Affine <sec-affine-class>` forms. It is
deliberately a plain ``Eigen::Matrix<Affine,-1,1>`` type alias rather than a
derived wrapper class, so every Eigen constructor, expression-template
operation and algorithm remains available without duplicating Eigen's API.

Unlike :ref:`AffineVariables <sec-affine-variables-class>`, an
``AffineVector`` does not itself introduce noise symbols: it is the
*working* type produced and consumed by computations, whether or not its
components happen to share noise symbols.

The distinction is fundamental:

.. list-table:: ``AffineVariables`` versus ``AffineVector``
   :widths: 35 32 33
   :header-rows: 1

   * - Property
     - ``AffineVariables``
     - ``AffineVector``
   * - Introduces noise symbols
     - Yes
     - No
   * - Main role
     - Declares independent uncertain quantities
     - Stores and computes with affine expressions
   * - Typical contents
     - Declared variables with dedicated noise symbols
     - Results of affine computations
   * - In-place compound assignment
     - Intentionally restricted
     - Available through the Eigen API
   * - Typical use
     - Start a dependency-aware computation
     - Continue and combine that computation


Creating and inspecting an affine vector
---------------------------------------------

The most common way to obtain an ``AffineVector`` is to convert it from
declared ``AffineVariables``, then combine it with ordinary Eigen
arithmetic:

.. tabs::

  .. group-tab:: Python

    .. literalinclude:: src.py
      :language: py
      :start-after: [affine-vector-1-beg]
      :end-before: [affine-vector-1-end]
      :dedent: 4

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-vector-1-beg]
      :end-before: [affine-vector-1-end]
      :dedent: 4

``itv()`` returns the ``IntervalVector`` enclosure of the whole vector,
component-wise — the vector equivalent of ``Affine::itv()``.


Dependency preservation
----------------------------

The dependency-preservation property demonstrated for a single ``Affine``
(see :ref:`Why affine arithmetic <sec-affine-class>`) extends directly to
vector-level linear combinations: reusing the same declared variables
several times in an expression does not widen the result the way plain
interval arithmetic would.

.. tabs::

  .. group-tab:: Python

    .. literalinclude:: src.py
      :language: py
      :start-after: [affine-vector-2-beg]
      :end-before: [affine-vector-2-end]
      :dedent: 4

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-vector-2-beg]
      :end-before: [affine-vector-2-end]
      :dedent: 4


Compound assignment
------------------------

Unlike ``AffineVariables``, a plain ``AffineVector`` supports the usual
compound assignment operators, since its scalar type (``Affine``) does not
delete them:

.. tabs::

  .. group-tab:: Python

    .. literalinclude:: src.py
      :language: py
      :start-after: [affine-vector-3-beg]
      :end-before: [affine-vector-3-end]
      :dedent: 4

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-vector-3-beg]
      :end-before: [affine-vector-3-end]
      :dedent: 4

.. note::

  This is a practical rule of thumb for choosing between the two types:
  use ``AffineVariables`` to *declare* the uncertain quantities of a
  problem, and ``AffineVector`` to *compute* with them (including
  accumulating results in place).


Block access
----------------

Being a plain Eigen column vector, the usual block accessors
(``head``, ``tail``, ``segment``) are available in C++; the Python binding
exposes the equivalent ``subvector(start_id, end_id)`` instead. They return
Eigen block expressions and can be materialized as an ``AffineVector`` when
a concrete vector is required:

.. tabs::

  .. group-tab:: Python

    .. literalinclude:: src.py
      :language: py
      :start-after: [affine-vector-4-beg]
      :end-before: [affine-vector-4-end]
      :dedent: 4

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-vector-4-beg]
      :end-before: [affine-vector-4-end]
      :dedent: 4


Transpose
-------------

.. tabs::

  .. group-tab:: Python

    .. literalinclude:: src.py
      :language: py
      :start-after: [affine-vector-5-beg]
      :end-before: [affine-vector-5-end]
      :dedent: 4

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-vector-5-beg]
      :end-before: [affine-vector-5-end]
      :dedent: 4

.. note::

  ``transpose()`` provides the corresponding Eigen row expression, which can
  be materialized as an :ref:`AffineRow <sec-affine-matrix-class>`, exactly
  like it would for any Eigen column vector.


Nonlinear functions, component by component
------------------------------------------------

There is no dedicated vector overload for the scalar nonlinear functions
listed on the :ref:`Affine class page <sec-affine-class>`. Applying
``sin``, ``sqr``, ``exp``, etc. to every component is therefore done one
component at a time:

.. tabs::

  .. group-tab:: Python

    .. literalinclude:: src.py
      :language: py
      :start-after: [affine-vector-6-beg]
      :end-before: [affine-vector-6-end]
      :dedent: 4

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-vector-6-beg]
      :end-before: [affine-vector-6-end]
      :dedent: 4


Nonlinear functions available
--------------------------------

The scalar nonlinear functions available for ``Affine`` can be applied
component by component. The complete list is given on the
:ref:`Affine class page <sec-affine-class>`. The interval-valued functions
``floor``, ``ceil`` and ``integer`` are likewise applied component by
component, but their results do not preserve affine dependencies.
The :ref:`Affine class page <sec-affine-class>` gives the complete list and
the corresponding scalar semantics.


Printing
------------

.. tabs::

  .. group-tab:: Python

    .. literalinclude:: src.py
      :language: py
      :start-after: [affine-vector-7-beg]
      :end-before: [affine-vector-7-end]
      :dedent: 4

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-vector-7-beg]
      :end-before: [affine-vector-7-end]
      :dedent: 4

.. note::

  Each component prints its own full ``Affine`` representation, joined by
  ``" ; "`` and wrapped in ``"[ ... ]"``. This gets verbose quickly, which
  is why ``itv()`` is usually the more practical way to inspect an
  ``AffineVector`` at a glance, reserving the full printout for close
  debugging of a specific component.
