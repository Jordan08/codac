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
components happen to be correlated with one another.


Creating and inspecting an affine vector
---------------------------------------------

The most common way to obtain an ``AffineVector`` is to convert it from
declared ``AffineVariables``, then combine it with ordinary Eigen
arithmetic:

.. tabs::

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
(``head``, ``tail``, ``segment``) are available and return an
``AffineVector`` over the requested range:

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-vector-4-beg]
      :end-before: [affine-vector-4-end]
      :dedent: 4


Transpose
-------------

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-vector-5-beg]
      :end-before: [affine-vector-5-end]
      :dedent: 4

.. note::

  ``transpose()`` turns an ``AffineVector`` into an
  :ref:`AffineRow <sec-affine-matrix-class>`, exactly like it would for any
  Eigen column vector.


Nonlinear functions, component by component
------------------------------------------------

There is no vectorized shortcut for elementary functions the way there is
for the linear compound assignments above (``+=``, ``*=``...): applying
``sin``, ``sqr``, etc. to every component is done one component at a time:

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-vector-6-beg]
      :end-before: [affine-vector-6-end]
      :dedent: 4


Printing
------------

.. tabs::

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
