.. _sec-affine-matrix-class:

The AffineMatrix class
==========================

  Main author: `Jordan Ninin <https://www.ensta-bretagne.fr/ninin/>`_

``AffineMatrix`` (an alias for ``AffineMainMatrix<AF_Default>``) is a
dynamic-size matrix of :ref:`Affine <sec-affine-class>` forms, a plain
``Eigen::Matrix<Affine,-1,-1>`` type alias, exactly like
:ref:`AffineVector <sec-affine-vector-class>` is for column vectors.

It is a **working container**, not a declaration mechanism: assigning an
``AffineVariables`` component to a matrix coefficient preserves that
component's noise symbol, while constructing a coefficient from an
``Interval`` does not introduce a new reusable noise symbol.


Products with a vector of affine forms
-------------------------------------------

The most common use case is multiplying a real-valued matrix by a vector
of declared affine variables — for instance evaluating a state-transition
or observation matrix on an uncertain state:

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-matrix-1-beg]
      :end-before: [affine-matrix-1-end]
      :dedent: 4

``AffineVariables`` (not just ``AffineVector``) can appear directly on
either side of such a product: it is converted internally to an
``AffineVector`` before the multiplication, so the result always has
``Affine`` (not ``AffineVarMain``) coefficients.


Supported matrix products
--------------------------

The affine matrix/vector types follow the corresponding Eigen dimensions:

.. list-table:: Common products
   :widths: 35 65
   :header-rows: 1

   * - Expression
     - Result
   * - ``RealMatrix * AffineVector``
     - ``AffineVector``
   * - ``RealMatrix * AffineVariables``
     - ``AffineVector``; the declared variables are converted to affine
       expressions without losing their noise symbols.
   * - ``AffineVector * RealRow``
     - ``AffineMatrix``
   * - ``AffineRow * RealMatrix``
     - ``AffineRow``
   * - ``AffineMatrix * AffineVector``
     - ``AffineVector``
   * - ``AffineMatrix * AffineMatrix``
     - ``AffineMatrix``

The exact scalar combinations accepted by Eigen's expression templates are
also available for affine and interval coefficients where the corresponding
Codac overloads are defined.

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-matrix-products-beg]
      :end-before: [affine-matrix-products-end]
      :dedent: 4


Dependency preservation through matrix products
-----------------------------------------------------

This is typically where affine arithmetic pays off compared to interval
arithmetic in practice: evaluating the *same* declared variables through
one or more matrices and recombining the results linearly does not lose
dependency information, e.g. in a state-space / control computation that reuses the
same uncertain state vector several times:

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-matrix-2-beg]
      :end-before: [affine-matrix-2-end]
      :dedent: 4


Building a matrix from affine forms
-----------------------------------------

``AffineMatrix`` is a general-purpose container: its components can come
from anywhere — declared variables, as below, or freestanding
``Affine(Interval(...))``/``Affine(double)`` constants (in which case they
carry no ``eps_i`` term at all, see
:ref:`the Affine class page <sec-affine-class>`). ``itv()`` gives the
enclosing ``IntervalMatrix`` either way:

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-matrix-3-beg]
      :end-before: [affine-matrix-3-end]
      :dedent: 4


Shared noise symbols in matrix coefficients
--------------------------------------------

Several matrix coefficients can deliberately refer to the same declared
variable. The resulting matrix therefore carries shared noise symbols in
multiple entries:

.. code-block:: cpp

   AffineVariables v(2);
   v[0] = Interval(-1,1);
   v[1] = Interval(2,3);

   AffineMatrix M(2,2);
   M(0,0) = v[0];
   M(0,1) = v[1];
   M(1,0) = v[0];
   M(1,1) = v[1];

Here ``M(0,0)`` and ``M(1,0)`` share the same noise symbol, as do
``M(0,1)`` and ``M(1,1)``. This is different from constructing four
independent ``Affine`` coefficients from the same intervals.

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-matrix-shared-beg]
      :end-before: [affine-matrix-shared-end]
      :dedent: 4


AffineRow
---------

``AffineRow`` is the row-vector counterpart of ``AffineVector``. It uses the
same scalar type, ``Affine``, and therefore preserves any noise symbols
carried by its components. It is primarily encountered through Eigen
operations such as ``transpose()`` and ``row(i)``.


Transpose, rows and columns
--------------------------------

Being a plain Eigen matrix alias, the usual accessors are available:
``transpose()`` swaps rows and columns, ``row(i)`` returns an
``AffineRow`` (``AffineMainRow<T>``, alias of ``Eigen::Matrix<Affine,1,-1>``),
and ``col(j)`` returns an :ref:`AffineVector <sec-affine-vector-class>`:

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-matrix-4-beg]
      :end-before: [affine-matrix-4-end]
      :dedent: 4

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-matrix-5-beg]
      :end-before: [affine-matrix-5-end]
      :dedent: 4


Printing
------------

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-matrix-6-beg]
      :end-before: [affine-matrix-6-end]
      :dedent: 4

.. note::

  As for ``AffineVector``, each component prints its own full ``Affine``
  representation — rows are separated by newlines and wrapped in
  ``"[ ... ]"``, the whole matrix in ``"[...]"``. Note that a component
  built as a freestanding ``Affine(Interval(...))`` or ``Affine(double)``
  constant (rather than derived from ``AffineVariables``) carries no
  ``eps_i`` term at all: all of its uncertainty sits in the error bound,
  as described in :ref:`the Affine class page <sec-affine-class>`.
