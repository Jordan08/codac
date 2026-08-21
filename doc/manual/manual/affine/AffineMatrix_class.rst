.. _sec-affine-matrix-class:

The AffineMatrix class
==========================

  Main author: `Jordan Ninin <https://www.ensta-bretagne.fr/ninin/>`_

``AffineMatrix`` (an alias for ``AffineMainMatrix<AF_Default>``) is a
dynamic-size matrix of :ref:`Affine <sec-affine-class>` forms, a plain
``Eigen::Matrix<Affine,-1,-1>`` type alias, exactly like
:ref:`AffineVector <sec-affine-vector-class>` is for column vectors.


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


Dependency preservation through matrix products
-----------------------------------------------------

This is typically where affine arithmetic pays off compared to interval
arithmetic in practice: evaluating the *same* declared variables through
one or more matrices and recombining the results linearly does not lose
correlation, e.g. in a state-space / control computation that reuses the
same uncertain state vector several times:

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-matrix-2-beg]
      :end-before: [affine-matrix-2-end]
      :dedent: 4


Matrix of arbitrary affine forms
-------------------------------------

``AffineMatrix`` components do not need to originate from
``AffineVariables`` at all — like ``AffineVector``, it is a general-purpose
container, and ``itv()`` gives the enclosing ``IntervalMatrix``:

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-matrix-3-beg]
      :end-before: [affine-matrix-3-end]
      :dedent: 4


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
