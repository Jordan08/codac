.. _sec-affine-class:

The Affine class
=================

  Main author: `Jordan Ninin <https://www.ensta-bretagne.fr/ninin/>`_

The ``Affine`` class represents a first-order affine approximation of a real
quantity over a set of noise symbols shared with other affine forms:

.. math::

  \hat{x} = x_0 + \sum_{i=1}^{n} x_i\,\varepsilon_i + x_{err}\,[-1,1],
  \qquad \varepsilon_i \in [-1,1].

Unlike ``Interval``, an affine form keeps track of *correlations* between
quantities that share the same noise symbols :math:`\varepsilon_i`. This
often tightens the enclosure of a composed expression compared to plain
interval arithmetic, in particular when a quantity is reused several times
in the same computation (a situation known as the *wrapping effect* in
interval analysis).

.. seealso::

  For the theoretical background, see chapter 2 of J. Ninin's PhD thesis,
  `available here <https://theses.hal.science/tel-04275036v1>`_.


Creating affine forms
----------------------

A constant, or a plain enclosure with no dependency information, can be
built directly from a real number or an ``Interval``:

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-class-1-beg]
      :end-before: [affine-class-1-end]
      :dedent: 4

.. note::

  An ``Affine`` built this way from an ``Interval`` does **not** carry a
  dedicated noise symbol: it behaves like a fresh, uncorrelated quantity.
  Two such affine forms built from the same interval will *not* cancel out
  when subtracted from one another, exactly like two independent
  ``Interval`` values would not. To create quantities that stay correlated
  through a computation, use ``AffineVariables`` (see below).


Creating affine variables
---------------------------

``AffineVariables`` (an alias for ``AffineVarMainVector``) is the only way
to introduce new, independent noise symbols. Each component of an
``AffineVariables`` vector is bound to its own symbol :math:`\varepsilon_i`,
shared by every affine form later derived from it. See
:ref:`the dedicated AffineVariables page <sec-affine-variables-class>` for
full details; a short example:

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-class-2-beg]
      :end-before: [affine-class-2-end]
      :dedent: 4

.. warning::

  ``AffineVariables`` and ``Affine`` are intentionally not
  interchangeable: an existing ``AffineVariables`` component can be
  reassigned from an ``Interval`` or a ``double`` (this replaces its value
  while preserving that component's own symbol), but assigning a plain
  ``Affine`` into it does not even compile — ``operator=(const Affine&)``
  is deleted on purpose, since silently accepting it would let a
  never-declared dependency be fabricated after the fact.


Affine form properties
------------------------

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-class-3-beg]
      :end-before: [affine-class-3-end]
      :dedent: 4

.. list-table:: Properties of an affine form :math:`\hat{x}`
   :widths: 30 70
   :header-rows: 1

   * - Method
     - Description
   * - ``itv()``
     - Interval enclosure :math:`\left[x_0-\sum|x_i|-x_{err},\,x_0+\sum|x_i|+x_{err}\right]`
   * - ``noise_count()``
     - Total number of noise symbols in the shared context (not just the
       ones with a non-zero coefficient here)
   * - ``noise(i)``
     - Coefficient :math:`x_i` of the noise symbol :math:`\varepsilon_i`
   * - ``err()``
     - Accumulated rounding/linearization error term :math:`x_{err}`
   * - ``is_empty()``, ``is_active()``, ``is_unbounded()``
     - Status predicates, mirroring the equivalent ``Interval`` predicates


Affine arithmetic
--------------------

Addition, subtraction, and scaling by a real constant are performed
*exactly* in affine arithmetic — no linearization is involved, and no
tightness is lost, whether or not the operands share noise symbols:

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-class-4-beg]
      :end-before: [affine-class-4-end]
      :dedent: 4


Elementary functions
------------------------

Nonlinear functions are supported through a first-order linearization
(Chebyshev or MinRange, see below), which guarantees a *sound* enclosure of
the true range but not necessarily the tightest possible one:

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-class-5-beg]
      :end-before: [affine-class-5-end]
      :dedent: 4

.. list-table:: Supported elementary functions
   :widths: 50 50
   :header-rows: 1

   * - Function
     - Function
   * - ``sqr``, ``sqrt``, ``inv``, ``pow``, ``root``
     - ``exp``, ``log``
   * - ``sin``, ``cos``, ``tan``
     - ``asin``, ``acos``, ``atan``, ``atan2``
   * - ``sinh``, ``cosh``, ``tanh``
     - ``asinh``, ``acosh``, ``atanh``
   * - ``abs``, ``sign``, ``chi``
     - ``floor``, ``ceil``, ``integer``

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-class-9-beg]
      :end-before: [affine-class-9-end]
      :dedent: 4

.. note::

  ``atan2`` is computed through composition with the affine division and
  ``atan`` when the box does not straddle the ``x=0`` branch cut, and
  soundly falls back to the plain ``Interval`` result otherwise.


Linearization mode
----------------------

Two linearization strategies are available for elementary functions:
Chebyshev (the default, minimizes the average approximation error) and
MinRange (minimizes the worst-case error). The mode is a global setting:

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-class-6-beg]
      :end-before: [affine-class-6-end]
      :dedent: 4

.. warning::

  ``change_mode`` affects every ``Affine`` computation from the point it is
  called, not just the affine forms created afterwards. Restore the default
  mode once you are done, especially if the change was meant to be local to
  a specific computation.


Comparisons and equality
----------------------------

``operator<`` and ``operator>`` behave like their ``Interval`` counterparts
and return a :ref:`BoolInterval <sec-intervals-boolinterval-class>`. Equality,
however, is defined purely in terms of the interval enclosure:

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-class-7-beg]
      :end-before: [affine-class-7-end]
      :dedent: 4


Why affine arithmetic: dependency preservation
---------------------------------------------------

The reason to reach for ``AffineVariables`` instead of plain ``Interval``
computations is that a quantity built from a shared noise symbol keeps its
identity through a chain of operations. Subtracting a variable from itself
collapses to zero — something plain interval arithmetic cannot do, since
:math:`[-2,3]-[-2,3]=[-5,5]` under the pessimistic assumption that both
occurrences are independent:

.. tabs::

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-class-8-beg]
      :end-before: [affine-class-8-end]
      :dedent: 4

This property extends to any linear combination reusing the same
``AffineVariables`` components, including through matrix products — which
is typically where affine arithmetic pays off compared to interval
arithmetic in practice (e.g. repeated evaluation of the same uncertain
state vector through several matrices).


Limitations
--------------

- Union and intersection of two affine forms (``operator|=``,
  ``operator&=``) are intentionally not implemented: a first-order affine
  representation cannot in general represent the union or intersection of
  two such sets without losing the affine structure, so these operators
  are deleted rather than silently degrading to an interval result.
- ``min`` and ``max`` return a plain ``Interval``, not an ``Affine``: since
  these functions are not differentiable, composing them further loses any
  affine correlation with the rest of the computation from that point on.
- Compound assignment operators (``+=``, ``-=``, ``*=``, ``/=`` with a
  scalar) are not available on ``AffineVariables`` components, to avoid
  accidentally mutating a bound noise symbol in place.
- This module is currently C++ only; Python bindings are not yet
  available.
