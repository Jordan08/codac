.. _sec-affine-class:

The Affine class
=================

  Main author: `Jordan Ninin <https://www.ensta-bretagne.fr/ninin/>`_

The ``Affine`` class represents a first-order affine approximation of a real
quantity over a set of noise symbols shared with other affine forms:

.. math::

  \hat{x} = x_0 + \sum_{i=1}^{n} x_i\,\varepsilon_i + x_{err}\,[-1,1],
  \qquad \varepsilon_i \in [-1,1].

Unlike ``Interval``, an affine form keeps track of first-order *dependencies*
between quantities that share the same noise symbols :math:`\varepsilon_i`.
This often tightens the enclosure of a composed expression compared to
plain interval arithmetic, in particular when a quantity is reused several
times in the same computation (a situation known as the *wrapping effect*
in interval analysis).

Understanding an affine form
----------------------------

An affine form contains three distinct kinds of information:

.. list-table:: Components of an affine form
   :widths: 25 75
   :header-rows: 1

   * - Term
     - Meaning
   * - ``x0``
     - Central value of the affine approximation.
   * - ``xi * eps_i``
     - A reusable dependency term. The same noise symbol in another affine
       form represents the same underlying uncertainty and therefore carries
       dependency information between the two forms.
   * - ``xerr * [-1,1]``
     - Residual, non-reusable uncertainty, including linearization and
       floating-point errors. It is deliberately not associated with a new
       noise symbol.

The practical benefit of affine arithmetic comes from preserving the
``eps_i`` terms through subsequent computations rather than replacing every
uncertainty by an unrelated interval at each operation. A noise-symbol
coefficient represents uncertainty whose identity can be reused by later
affine expressions; the error term represents uncertainty for which no
reusable dependency is maintained.

.. seealso::

  For the theoretical background, see chapter 2 of J. Ninin's PhD thesis,
  `available here <https://theses.hal.science/tel-04275036v1>`_.


Creating affine forms
----------------------

A constant, or a plain enclosure with no dependency information, can be
built directly from a real number or an ``Interval``. Contrast this with a
component of ``AffineVariables`` (introduced below), or a plain ``Affine``
copied from one, which *do* carry a dedicated noise symbol:

.. tabs::

  .. group-tab:: Python

    .. literalinclude:: src.py
      :language: py
      :start-after: [affine-class-1-beg]
      :end-before: [affine-class-1-end]
      :dedent: 4

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-class-1-beg]
      :end-before: [affine-class-1-end]
      :dedent: 4

.. note::

  An ``Affine`` built directly from an ``Interval`` or a ``double`` does
  **not** carry a dedicated noise symbol (no ``eps_i`` term at all, as the
  printed values above show for ``y`` and ``z``): it behaves like a fresh,
  uncorrelated quantity, with all of its uncertainty folded into the error
  bound. Two such affine forms built from the same interval will *not*
  cancel out when subtracted from one another, exactly like two
  independent ``Interval`` values would not. Copying an existing
  ``AffineVariables`` component into a plain ``Affine`` (like ``x`` above),
  by contrast, *does* preserve its noise symbol — see
  :ref:`AffineVariables <sec-affine-variables-class>` for full details.


Creating affine variables
---------------------------

``AffineVariables`` (an alias for ``AffineVarMainVector``) is the only way
to introduce new, independent noise symbols. Each component of an
``AffineVariables`` vector is bound to its own symbol :math:`\varepsilon_i`,
shared by every affine form later derived from it. See
:ref:`the dedicated AffineVariables page <sec-affine-variables-class>` for
full details; a short example:

.. tabs::

  .. group-tab:: Python

    .. literalinclude:: src.py
      :language: py
      :start-after: [affine-class-2-beg]
      :end-before: [affine-class-2-end]
      :dedent: 4

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

  .. group-tab:: Python

    .. literalinclude:: src.py
      :language: py
      :start-after: [affine-class-3-beg]
      :end-before: [affine-class-3-end]
      :dedent: 4

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

``mig()``, ``mag()``, ``smag()`` and ``smig()`` mirror the equivalent
``Interval`` methods, applied to the affine form's own interval enclosure;
``volume()`` is an alias for ``diam()`` in this scalar, 1-dimensional case:

.. tabs::

  .. group-tab:: Python

    .. literalinclude:: src.py
      :language: py
      :start-after: [affine-class-10-beg]
      :end-before: [affine-class-10-end]
      :dedent: 4

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-class-10-beg]
      :end-before: [affine-class-10-end]
      :dedent: 4

.. list-table:: Magnitude of an affine form :math:`\hat{x}`
   :widths: 30 70
   :header-rows: 1

   * - Method
     - Description
   * - ``mig()``
     - Mignitude: :math:`\mathrm{lb}(\hat{x})` if :math:`\hat{x}>0`,
       :math:`-\mathrm{ub}(\hat{x})` if :math:`\hat{x}<0`, :math:`0` otherwise
   * - ``mag()``
     - Magnitude: :math:`\max(|\mathrm{lb}(\hat{x})|,|\mathrm{ub}(\hat{x})|)`
   * - ``smag()``
     - Signed magnitude: whichever bound has the larger absolute value,
       kept with its sign
   * - ``smig()``
     - Signed mignitude: :math:`\mathrm{lb}(\hat{x})` if
       :math:`\mathrm{lb}(\hat{x})>0`, :math:`\mathrm{ub}(\hat{x})` if
       :math:`\mathrm{ub}(\hat{x})<0`, :math:`0` otherwise
   * - ``volume()``
     - Alias for ``diam()``


Enlarging an affine form
----------------------------

``inflate(r)`` adds :math:`[-r,+r]` to the affine form's remainder error
term, in place, widening its interval enclosure symmetrically by
:math:`r` on each side (the result is outward-rounded, so it may be a
tiny bit wider than the mathematical :math:`r`-enlargement):

.. tabs::

  .. group-tab:: Python

    .. literalinclude:: src.py
      :language: py
      :start-after: [affine-class-11-beg]
      :end-before: [affine-class-11-end]
      :dedent: 4

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-class-11-beg]
      :end-before: [affine-class-11-end]
      :dedent: 4


Affine arithmetic
--------------------

Addition, subtraction, and scaling by a real constant are performed
*exactly* in affine arithmetic — no linearization is involved, and no
tightness is lost, whether or not the operands share noise symbols. Here
``x`` and ``y`` are copies of two declared variables, so the printed
result keeps both ``eps_0`` and ``eps_1`` visible throughout — including
in the last line, where ``z - x`` cancels ``x``'s own contribution exactly
and recovers ``2*x - x = x``:

.. tabs::

  .. group-tab:: Python

    .. literalinclude:: src.py
      :language: py
      :start-after: [affine-class-4-beg]
      :end-before: [affine-class-4-end]
      :dedent: 4

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-class-4-beg]
      :end-before: [affine-class-4-end]
      :dedent: 4


Nonlinear operations and elementary functions
-------------------------------------------------

Nonlinear operations are supported through first-order affine
linearization (Chebyshev or MinRange, see below). The resulting affine form
is a sound enclosure of the true range, but the enclosure is not
necessarily the tightest possible one.

For nonlinear unary functions, the nonlinear remainder is accumulated in
the error term. Existing noise symbols from the input are retained in the
linear part of the result. Thus nonlinear operations preserve useful
dependency information while accounting conservatively for the terms that
cannot be represented linearly. Some operations are deliberately
interval-valued instead; these are identified explicitly in the list below.

The following nonlinear operations are available on ``Affine``:

.. list-table:: Nonlinear operations
   :widths: 28 72
   :header-rows: 1

   * - Operation
     - Description
   * - ``inv(x)``
     - Reciprocal :math:`1/x`.
   * - ``sqr(x)``
     - Square :math:`x^2`.
   * - ``sqrt(x)``
     - Square root.
   * - ``exp(x)``
     - Exponential.
   * - ``log(x)``
     - Natural logarithm.
   * - ``pow(x, n)``
     - Integer power, with ``n`` an ``int``.
   * - ``pow(x, d)``
     - Real power, with ``d`` a ``double``.
   * - ``pow(x, y)``
     - Power with an ``Interval`` exponent.
   * - ``pow(x, y)``
     - Power with an ``Affine`` exponent. The exponent's affine dependency
       structure is evaluated through its interval enclosure.
   * - ``pow(c, y)``
     - Power with a real ``double`` base ``c`` and an ``Affine`` exponent.
       The exponent's affine dependency structure is evaluated through its
       interval enclosure.
   * - ``pow(I, y)``
     - Power with an ``Interval`` base ``I`` and an ``Affine`` exponent. The
       exponent's affine dependency structure is evaluated through its
       interval enclosure.
   * - ``x * y``
     - Product of two affine forms; nonlinear terms are enclosed in the
       error term.
   * - ``x / y``
     - Quotient of two affine forms; nonlinear terms are conservatively
       enclosed in the error term.
   * - ``root(x, n)``
     - ``n``-th root.
   * - ``sin(x)``, ``cos(x)``, ``tan(x)``
     - Trigonometric functions.
   * - ``asin(x)``, ``acos(x)``, ``atan(x)``
     - Inverse trigonometric functions.
   * - ``atan2(y, x)``
     - Two-argument arctangent.
   * - ``sinh(x)``, ``cosh(x)``, ``tanh(x)``
     - Hyperbolic functions.
   * - ``asinh(x)``, ``acosh(x)``, ``atanh(x)``
     - Inverse hyperbolic functions.
   * - ``abs(x)``
     - Absolute value.
   * - ``floor(x)``, ``ceil(x)``
     - Floor and ceiling. These operations return an ``Interval`` and do not
       preserve affine dependencies.
   * - ``integer(x)``
     - Largest integer interval contained in the interval represented by
       ``x``. This operation returns an ``Interval`` and does not preserve
       affine dependencies.
   * - ``sign(x)``
     - Enclosure of the sign of ``x``; if zero belongs to the enclosure,
       the result is ``[-1,1]`` by convention.
   * - ``chi(a,b,c)``
     - Conditional selection based on the sign of ``a``: selects ``b`` if
       ``a <= 0``, ``c`` if ``a > 0``, and encloses both otherwise.

The functions ``min``, ``max``, ``floor``, ``ceil`` and ``integer`` are
interval-valued operations. They intentionally leave the affine dependency
model because there is no linear dependency between the input and output.
``chi`` is different: it returns an ``Affine`` and is intended for
conditional selection between affine expressions.

.. tabs::

  .. group-tab:: Python

    .. literalinclude:: src.py
      :language: py
      :start-after: [affine-class-nonlinear-list-beg]
      :end-before: [affine-class-nonlinear-list-end]
      :dedent: 4

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-class-nonlinear-list-beg]
      :end-before: [affine-class-nonlinear-list-end]
      :dedent: 4

.. tabs::

  .. group-tab:: Python

    .. literalinclude:: src.py
      :language: py
      :start-after: [affine-class-5-beg]
      :end-before: [affine-class-5-end]
      :dedent: 4

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-class-5-beg]
      :end-before: [affine-class-5-end]
      :dedent: 4


.. tabs::

  .. group-tab:: Python

    .. literalinclude:: src.py
      :language: py
      :start-after: [affine-class-9-beg]
      :end-before: [affine-class-9-end]
      :dedent: 4

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-class-9-beg]
      :end-before: [affine-class-9-end]
      :dedent: 4

.. note::

  ``atan2(y,x)`` is handled using the affine division/``atan`` construction
  when the required domain conditions are satisfied. In cases where that
  construction cannot provide a sound affine representation, the
  implementation falls back to the corresponding interval enclosure.


Linearization mode
----------------------

Two linearization strategies are available for nonlinear functions:

* **Chebyshev** (the default): uses the minimax/Chebyshev affine
  approximation, minimizing the maximum deviation of the linearization
  from the nonlinear function over the input interval.
* **MinRange**: chooses the slope from the endpoint derivative so that the
  resulting affine range is minimized for the monotone function on the
  interval.

The implementation follows the classical affine-arithmetic construction.
For an interval :math:`x=[a,b]` and a differentiable function :math:`f`,
the Chebyshev construction used in ``codac2_AffineMain.h`` computes the
secant slope

.. math::

   \alpha = \frac{f(b)-f(a)}{b-a},

then finds a point :math:`u` satisfying

.. math::

   f'(u)=\alpha.

The two endpoint residuals are

.. math::

   d_a = f(a)-\alpha a, \qquad
   d_b = f(b)-\alpha b,

while the interior extremum is evaluated as

.. math::

   d_u = f(u)-\alpha u.

The code encloses the residual extrema with
:math:`[d_{min},d_{max}]`, where :math:`d_{min}` includes the endpoint
minimum and :math:`d_{max}` includes the interior maximum. It then uses

.. math::

   \beta = \operatorname{mid}([d_{min},d_{max}]),\qquad
   \delta = \operatorname{rad}([d_{min},d_{max}]),

so that the nonlinear function is enclosed by

.. math::

   f(x) \approx \alpha x + \beta + \delta[-1,1].

The implementation applies this construction to the existing affine form:
the linear part reuses its noise symbols, while the residual radius is
added to the affine error term. The code also contains dedicated domain
checks and special cases for functions such as ``inv``, ``sqrt`` and
``exp``.

The secant-slope construction above assumes :math:`f` is monotonic with a
single sign of curvature over :math:`x`, so that :math:`g(t)=f(t)-\alpha t`
has exactly one interior extremum. That assumption fails for functions
that can have an inflection point inside :math:`x`. For those cases —
the trigonometric functions (``sin``, ``cos``, ``tan``), their inverses
(``asin``, ``acos``, ``atan``), the hyperbolic functions (``sinh``,
``tanh``) and their inverses (``asinh``, ``atanh``), and ``pow(x, n)`` for
odd integer :math:`n` — ``codac2_AffineMain.h`` uses a different,
two-point Chebyshev construction to compute :math:`\alpha`. For
:math:`x=[a,b]`, :math:`f` is sampled at the Chebyshev (Gauss) nodes
:math:`\pm x_0` of the interval, with :math:`x_0=1/\sqrt2`:

.. math::

   x_{b0} = \frac{(b-a)x_0+(a+b)}{2}, \qquad
   x_{b1} = \frac{-(b-a)x_0+(a+b)}{2},

.. math::

   c_1 = x_0\big(f(x_{b0})-f(x_{b1})\big), \qquad
   \alpha = \frac{2c_1}{b-a}.

This :math:`\alpha` is the degree-1 coefficient of the two-point discrete
Chebyshev expansion of :math:`f` on :math:`x`, and is generally different
from the secant slope. Once :math:`\alpha` is obtained this way,
:math:`\beta` and :math:`\delta` are still built as above: the residual
:math:`f(t)-\alpha t` is evaluated at the endpoints of :math:`x` and at
the interior point(s) :math:`u` solving :math:`f'(u)=\alpha` that fall
inside :math:`x` — a single pair :math:`\pm u` for the functions above
other than ``sin``, ``cos`` and ``tan``, since their derivative is even,
and one candidate per period for the periodic ``sin``, ``cos`` and
``tan`` — and :math:`\beta`, :math:`\delta` are set to the midpoint and
radius of the union of these residual intervals.

For the MinRange construction, the implementation uses the derivative at
the endpoint having the smallest derivative magnitude (for the monotone
case), builds the endpoint residual band

.. math::

   [f(a)-\alpha a,\; f(b)-\alpha b],

with the corresponding orientation for decreasing functions, and then sets
:math:`\beta` to the midpoint and :math:`\delta` to the radius of that
band.

These constructions are described in the affine-arithmetic literature. See
chapter 2 of J. Ninin's PhD thesis and the foundational paper by Stolfi and
de Figueiredo:

* `J. Ninin, PhD thesis <https://theses.hal.science/tel-04275036v1>`_
* `J. Stolfi and L. H. de Figueiredo, An Introduction to Affine Arithmetic,
  TEMA 4(3), 2003 <https://doi.org/10.5540/tema.2003.04.03.0297>`_

The mode is stored as a thread-local setting:

.. tabs::

  .. group-tab:: Python

    .. literalinclude:: src.py
      :language: py
      :start-after: [affine-class-6-beg]
      :end-before: [affine-class-6-end]
      :dedent: 4

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

``operator<`` and ``operator>`` return a :ref:`BoolInterval
<sec-intervals-boolinterval-class>` rather than a C++ ``bool``. This is
necessary because an affine enclosure may not contain enough information to
decide the truth value uniquely:

* ``True`` means the relation is certainly true;
* ``False`` means it is certainly false;
* ``Maybe`` means that both truth values remain possible from the available
  enclosure.

Equality is different: ``operator==`` is defined from the interval
enclosures rather than from the affine coefficients or noise-symbol
identities. Consequently, two independent affine forms with the same
interval enclosure can compare equal.

.. tabs::

  .. group-tab:: Python

    .. literalinclude:: src.py
      :language: py
      :start-after: [affine-class-7-beg]
      :end-before: [affine-class-7-end]
      :dedent: 4

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-class-7-beg]
      :end-before: [affine-class-7-end]
      :dedent: 4


Set-relation predicates
----------------------------

These predicates all mirror the equivalent ``Interval`` predicate, applied
to the affine form's own interval enclosure — they return a plain C++
``bool``, unlike ``operator<``/``operator>`` above. Each accepts either an
``Interval`` or another ``Affine`` (``contains()`` takes a plain
``double``):

.. tabs::

  .. group-tab:: Python

    .. literalinclude:: src.py
      :language: py
      :start-after: [affine-class-12-beg]
      :end-before: [affine-class-12-end]
      :dedent: 4

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [affine-class-12-beg]
      :end-before: [affine-class-12-end]
      :dedent: 4

.. list-table:: Set-relation predicates on an affine form :math:`\hat{x}`
   :widths: 30 70
   :header-rows: 1

   * - Method
     - Description
   * - ``is_subset(x)``
     - ``true`` iff this form's enclosure is a subset of ``x``
   * - ``is_superset(x)``
     - ``true`` iff this form's enclosure is a superset of ``x``
   * - ``is_disjoint(x)``
     - ``true`` iff this form's enclosure and ``x`` share no point at all
   * - ``intersects(x)``
     - ``true`` iff this form's enclosure and ``x`` share at least one point
   * - ``overlaps(x)``
     - ``true`` iff their intersection has non-zero volume (a single
       touching point does not count)
   * - ``contains(d)``
     - ``true`` iff the real value ``d`` belongs to this form's enclosure


Dependency information by operation
-------------------------------------

A useful way to understand the API is to distinguish operations that keep
the affine representation from operations that return an interval:

.. list-table:: Dependency behavior
   :widths: 30 35 35
   :header-rows: 1

   * - Operation
     - Result
     - Dependency information
   * - ``x + y``, ``x - y``
     - ``Affine``
     - Preserved exactly.
   * - ``a * x``, ``x / a`` for a non-zero real constant ``a``
     - ``Affine``
     - Preserved exactly.
   * - ``x * y``, ``x / y``
     - ``Affine``
     - Existing linear dependencies are retained; nonlinear terms are
       enclosed in the error term.
   * - ``f(x)`` for an elementary nonlinear function
     - ``Affine``
     - Existing noise symbols are retained through the linearization; the
       nonlinear remainder is added to the error term.
   * - ``min(x,y)``, ``max(x,y)``
     - ``Interval``
     - The affine dependency structure is no longer available in the
       result.
   * - ``floor(x)``, ``ceil(x)``, ``integer(x)``
     - ``Interval``
     - The result is interval-valued; affine dependencies are not preserved.
   * - ``x & y``, ``x | y``
     - ``Interval``
     - Interval intersection/hull operations are available as free
       functions, but are not affine-preserving operations.

This distinction is why ``Affine`` can be used throughout a computation
without repeatedly falling back to interval arithmetic, while operations
whose natural result is only an interval are explicitly documented as such.


Why affine arithmetic: dependency preservation
---------------------------------------------------

The reason to reach for ``AffineVariables`` instead of plain ``Interval``
computations is that a quantity built from a shared noise symbol keeps its
identity through a chain of operations. Subtracting a variable from itself
collapses to zero — something plain interval arithmetic cannot do, since
:math:`[-2,3]-[-2,3]=[-5,5]` under the pessimistic assumption that both
occurrences are independent. The contrast below is worth noting: ``v[0]``
and ``v[1]`` are broadcast to the *same* enclosure via ``init``, yet only
subtracting a variable from *itself* (``x``, ``y``) cancels — subtracting
two genuinely distinct variables (``z``) does not, even with an identical
range. The property also holds for a *derived* expression subtracted from
itself (``zz``), not just a raw declared variable:

.. tabs::

  .. group-tab:: Python

    .. literalinclude:: src.py
      :language: py
      :start-after: [affine-class-8-beg]
      :end-before: [affine-class-8-end]
      :dedent: 4

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


Limitations and interval-valued operations
-----------------------------------------------

- Union and intersection of two affine forms (``operator|=``,
  ``operator&=``) are intentionally not implemented: a first-order affine
  representation cannot in general represent the union or intersection of
  two such sets without losing the affine structure, so these operators
  are deleted rather than silently degrading to an interval result.
- ``min`` and ``max`` with affine operands return a plain ``Interval``, not
  an ``Affine``. The branch selected by ``min``/``max`` depends on the
  actual values of the operands, so a single affine form cannot in general
  represent the result while preserving the required dependency
  information. Converting the result to an interval also means that
  subsequent affine computations cannot reuse the previous noise-symbol
  structure through this operation. They should not be confused with
  ``chi``, which provides an affine result for conditional selection.
- Compound assignment operators (``+=``, ``-=``, ``*=``, ``/=`` with a
  scalar) are not available on ``AffineVariables`` components, to avoid
  accidentally mutating a bound noise symbol in place.

