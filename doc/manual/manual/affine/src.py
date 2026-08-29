#!/usr/bin/env python

#  Codac tests
# ----------------------------------------------------------------------------
#  \date       2026
#  \author     Jordan Ninin
#  \copyright  Copyright 2026 Codac Team
#  \license    GNU Lesser General Public License (LGPL)

import sys, os
import unittest
import math
from codac import *

class TestAffineManual(unittest.TestCase):

  def tests_Affine_manual(test):

    # [affine-class-1-beg]
    v = AffineVariables(1)
    v[0] = Interval(-1,5)      # a variable affine form, enclosure [-1,5]
                                # v[0] = [-1, 5] : 2 + 3 eps_0 + 0 [-1,1]
    x = Affine(v[0])           # a constant affine form, enclosure [-1,5]
                                # x = [-1, 5] : 2 + 3 eps_0 + 0 [-1,1]
    y = Affine(Interval(-1,5)) # an affine form, enclosure [-1,5], but no dedicated noise symbol
                                # y = [-1, 5] : 2  + 3 [-1,1]
    z = Affine(0.)              # a constant affine form, enclosure [0,0]
                                # z = [0, 0] : 0 + 0 [-1,1]
    v[0] = Interval(4.)        # v[0] = [4, 4] : 4 + 0 eps_0 + 0 [-1,1]
    v[0] = Interval(2,3)       # v[0] = [2, 3] : 2.5 + 0.5 eps_0 + 0 [-1,1]
    zz = Affine()                # unbounded, ]-oo,oo[, like Interval()
    # [affine-class-1-end]
    test.assertTrue(x.itv() == Interval(-1,5))
    test.assertTrue(y.itv() == Interval(-1,5))
    test.assertTrue(zz.itv() == Interval())
    test.assertTrue(z.itv() == Interval(0))

    # [affine-class-2-beg]
    # Only AffineVariables introduces new, independent noise symbols.
    # A plain Affine built from an Interval never carries a dedicated
    # symbol: it cannot be correlated with anything else later on.
    v = AffineVariables([[1,2],[-1,1],[3,4]])
    # v[0] = [1,2] : 0.5 + 0.5 eps_0 + 0 eps_1 + 0 eps_2 + 0 [-1,1],
    # v[1] = [-1,1] : 0 + 0 eps_0 + 1 eps_1 + 0 eps_2 + 0 [-1,1],
    # v[2] = [3,4] : 3.5 + 0 eps_0 + 0 eps_1 + 0.5 eps_2 + 0 [-1,1]
    # [affine-class-2-end]
    test.assertTrue(v[0].itv() == Interval(1,2))
    test.assertTrue(v[1].itv() == Interval(-1,1))
    test.assertTrue(v[2].itv() == Interval(3,4))

    # [affine-class-3-beg]
    v = AffineVariables([[1,2],[-1,1],[3,4]])

    v[0].itv()           # interval enclosure: [1,2]
    v[0].noise_count()   # total number of noise symbols in the shared context: 3
    v[0].noise(0)        # coefficient of eps_0 : 0.5 (the radius)
    v[0].noise(1)        # coefficient of eps_1 : 0
    v[0].noise(2)        # coefficient of eps_2 : 0
    v[0].err()           # accumulated rounding-error term: 0 (nothing computed yet)
    v.noise_index(0)     # index of the noise symbol associated with v[0] : 0
    v.noise_index(2)     # index of the noise symbol associated with v[2] : 2
    v[1]                    # v[1] = [-1,1] : 0 + 0 eps_0 + 1 eps_1 + 0 eps_2 + 0 [-1,1]
    v[1] = Interval(10,20)  # v[1] = [10,20] : 15 + 0 eps_0 + 5 eps_1 + 0 eps_2 + 0 [-1,1]
    # [affine-class-3-end]
    test.assertTrue(v[0].noise_count() == 3)
    test.assertTrue(v[0].noise(0) == 0.5)
    test.assertTrue(v[0].noise(1) == 0.0)
    test.assertTrue(v[0].noise(2) == 0.0)
    test.assertTrue(v[0].err() == 0.0)
    test.assertTrue(v[2].noise_count() == 3)
    test.assertTrue(v[2].noise(0) == 0.0)
    test.assertTrue(v[2].noise(1) == 0.0)
    test.assertTrue(v[2].noise(2) == 0.5)
    test.assertTrue(v[2].err() == 0.0)

    # [affine-class-4-beg]
    # Addition, subtraction and scaling by a constant are exact in affine
    # arithmetic: no linearization, no loss of tightness, whether or not
    # the operands are correlated.
    var = AffineVariables(2)
    var[0] = Interval(2,3)
    x = var[0]    # [2, 3]: 2.5 + 0.5 eps_0 + 0 eps_1 + 0 [-1,1]
    var[1] = Interval(1,3)
    y = var[1]    # [1, 3]: 2 + 0 eps_0 + 1 eps_1 + 0 [-1,1]
    z = x + y     # [3, 6]: 4.5 + 0.5 eps_0 + 1 eps_1 + 0 [-1,1]
    z = x - 1     # [1, 2]: 1.5 + 0.5 eps_0 + 0 eps_1 + 0 [-1,1]
    z = 2 * x     # [4, 6]: 5 + 1 eps_0 + 0 eps_1 + 0 [-1,1]
    z = z - x     # [2, 3]: 2.5 + 0.5 eps_0 + 0 eps_1 + 0 [-1,1]
    # [affine-class-4-end]
    test.assertTrue((x + y).itv() == Interval(3,6))
    test.assertTrue((x - 1).itv() == Interval(1,2))
    test.assertTrue((2 * x).itv() == Interval(4,6))

    # [affine-class-5-beg]
    # Nonlinear functions rely on a first-order (Chebyshev or MinRange)
    # linearization: the result is a sound, but not necessarily minimal,
    # enclosure of the true range.
    v = AffineVariables([[0, 1.5707963267948966]]) # [0, pi/2]
    y = sin(v[0]) #[-0.0544991, 1.25618] : 0.600837 + 0.52725 eps_0 + 0.128086 [-1,1]
    z = exp(v[0]) #[0.276127, 4.81048] : 2.54331 + 1.90524 eps_0 + 0.361937 [-1,1]
    # [affine-class-5-end]
    test.assertTrue(y.itv().is_superset(Interval(0,1)))
    test.assertTrue(z.itv().is_superset(Interval(1,math.exp(1.5707963267948966))))

    # [affine-class-nonlinear-list-beg]
    v = AffineVariables([[1,2],[2,3],[0.5,1.0]])
    x = v[0]
    y = v[1]
    z = v[2]

    a1 = inv(x)
    a2 = sqr(x)
    a3 = sqrt(x)
    a4 = exp(x)
    a5 = log(x)
    a6 = pow(x, 2)
    a7 = pow(x, 0.5)
    a8 = pow(x, Interval(1,2))
    a9 = pow(x, y)
    a10 = root(x, 2)
    a11 = sin(z)
    a12 = cos(z)
    a13 = tan(z)
    a14 = asin(z)
    a15 = acos(z)
    a16 = atan(z)
    a17 = sinh(z)
    a18 = cosh(z)
    a19 = tanh(z)
    a20 = asinh(z)
    a21 = acosh(x)
    a22 = atanh(z)
    a23 = atan2(z, x)
    a24 = abs(x)
    a25 = floor(x)
    a26 = ceil(x)
    a27 = integer(x)
    a28 = sign(z)
    a29 = chi(z, x, y)
    # [affine-class-nonlinear-list-end]

    # [affine-class-6-beg]
    # The linearization mode is a global (thread-local) setting shared by
    # every Affine instance; restore the default after use so that later
    # computations are not silently affected.
    Affine.change_mode(Affine.AF_Lin_MinRange)
    # ... nonlinear computations using the MinRange linearization ...
    Affine.change_mode(Affine.AF_Lin_Chebyshev) # back to the default
    # [affine-class-6-end]
    test.assertTrue(Affine.get_mode() == Affine.AF_Lin_Chebyshev)

    # [affine-class-9-beg]
    # atan2(y,x) is available for affine forms. When the affine
    # division/atan construction is not applicable, the implementation
    # falls back to a sound interval enclosure.
    var = AffineVariables([[1,2],[0.5,1]])
    x = var[0]
    y = var[1]
    z = atan2(y,x)
    # [affine-class-9-end]
    test.assertTrue(z.itv().is_subset(Interval(0, Interval.half_pi().ub())))

    # [affine-class-7-beg]
    # operator== only compares interval enclosures, not the underlying
    # coefficients or dependency structure: two affine forms with an
    # identical enclosure compare equal even if one of them originates
    # from a shared noise symbol and the other does not.
    var = AffineVariables(2)
    var.init(Interval(1,2))
    a = var[0]
    b = var[1]
    same_enclosure = (a == b) # true
    # [affine-class-7-end]
    test.assertTrue(same_enclosure)

    # [affine-class-12-beg]
    # These predicates mirror the equivalent Interval methods, applied to
    # this affine form's own interval enclosure; each accepts either an
    # Interval or another Affine (contains() takes a plain double).
    v = AffineVariables([[1,2],[0,3],[5,6]])
    x = v[0] # [1,2]
    y = v[1] # [0,3]
    z = v[2] # [5,6]

    x.is_subset(y)    # true:  [1,2] is a subset of [0,3]
    y.is_superset(x)  # true:  [0,3] is a superset of [1,2]
    x.is_disjoint(z)  # true:  [1,2] and [5,6] share no point at all
    x.intersects(y)   # true:  [1,2] and [0,3] share at least one point
    x.overlaps(y)     # true:  their intersection has non-zero volume
    x.contains(1.5)   # true:  1.5 belongs to [1,2]
    # [affine-class-12-end]
    test.assertTrue(x.is_subset(y))
    test.assertTrue(y.is_superset(x))
    test.assertTrue(x.is_disjoint(z))
    test.assertTrue(x.intersects(y))
    test.assertTrue(x.overlaps(y))
    test.assertTrue(x.contains(1.5))

    # [affine-class-8-beg]
    # Because AffineVariables preserves the noise symbol across
    # computations, subtracting a variable from itself is exact: the
    # dependency is not lost the way it would be with plain intervals
    # ([-2,3] - [-2,3] would give [-5,5], not {0}).
    v = AffineVariables(2)
    v.init(Interval(-2,3))
    x = v[0] - v[0] # [0]
    y = v[1] - v[1] # [0]
    z = v[0] - v[1] # [-5,5]
    zz = z  - z     # [0]
    # [affine-class-8-end]
    test.assertTrue(x.itv() == Interval(0.0))
    test.assertTrue(y.itv() == Interval(0.0))
    test.assertTrue(z.itv() == Interval(-5,5))
    test.assertTrue(zz.itv() == Interval(0.0))

    # [affine-class-10-beg]
    # mig(), mag(), smag() and smig() mirror the equivalent Interval
    # methods, applied to this affine form's own interval enclosure;
    # volume() is an alias for diam() in this scalar, 1-dimensional case.
    v = AffineVariables([[-3,5],[2,4]])
    x = v[0] # [-3,5]: straddles zero
    y = v[1] # [2,4]: strictly positive

    x.mig()     # mignitude: 0 (0 belongs to [-3,5])
    x.mag()     # magnitude: 5 (max(|-3|,|5|))
    x.smag()    # signed magnitude: 5 (the bound with the largest |.|)
    x.smig()    # signed mignitude: 0 (the enclosure straddles zero)
    x.volume()  # == x.diam(): 8

    y.mig()   # mignitude: 2 (y is strictly positive: its lower bound)
    y.mag()   # magnitude: 4
    y.smag()  # signed magnitude: 4
    y.smig()  # signed mignitude: 2 (y is strictly positive: its lower bound)
    # [affine-class-10-end]
    test.assertTrue(x.mig() == 0.0)
    test.assertTrue(x.mag() == 5.0)
    test.assertTrue(x.smag() == 5.0)
    test.assertTrue(x.smig() == 0.0)
    test.assertTrue(x.volume() == x.diam())
    test.assertTrue(x.volume() == 8.0)
    test.assertTrue(y.mig() == 2.0)
    test.assertTrue(y.mag() == 4.0)
    test.assertTrue(y.smag() == 4.0)
    test.assertTrue(y.smig() == 2.0)

    # [affine-class-11-beg]
    # inflate(r) widens the enclosure by [-r,+r] on each side, in place
    # (outward-rounded, so the result may be a hair wider than [0.5,2.5]).
    # A genuine copy is taken explicitly here: v[0] itself is a reference
    # into v's own storage (see the note on AffineVariables indexing
    # below), so inflating it in place would also widen v[0].
    v = AffineVariables(1)
    v[0] = Interval(1,2)
    x = Affine(v[0]) # [1,2]
    x.inflate(0.5)   # [0.5,2.5]
    # [affine-class-11-end]
    test.assertTrue(x.itv().is_superset(Interval(0.5,2.5)))


  def tests_AffineVariables_manual(test):

    # [affine-variables-complete-beg]
    # 1. Declare uncertain quantities.
    x = AffineVariables([
      [0., 1.],
      [2., 4.],
      [-1., 1.]
    ])

    # 2. Each component receives a distinct noise symbol.
    i0 = x.noise_index(0)
    i1 = x.noise_index(1)
    i2 = x.noise_index(2)

    # 3. Existing variables can be assigned a new interval without
    # changing their noise-symbol identity.
    x[1] = Interval(5., 6.)

    # 4. Use the declared variables in affine expressions.
    a = x[0] + 2. * x[1]
    b = x[1] - x[2]

    # 5. Reusing the same declared variable preserves its dependency.
    c = x[0] - x[0]

    # 6. Materialize the declared variables as an AffineVector when vector
    # operations are required.
    v = AffineVector(x)
    y = 2. * v
    # [affine-variables-complete-end]

    test.assertTrue(i0 == 0)
    test.assertTrue(i1 == 1)
    test.assertTrue(i2 == 2)
    test.assertTrue(x.noise_index(1) == i1)
    test.assertTrue(c.itv() == Interval(0.0))
    test.assertTrue(y.itv() == IntervalVector([[0.,2.],[10.,12.],[-2.,2.]]))

    # [affine-variables-1-beg]
    u = AffineVariables(3)                        # 3 unbound components, ]-oo,+oo[ each
    v = AffineVariables([[1,2],[-1,1],[3,4]])
    w = AffineVariables(Vector([1.5,0.,3.5]))     # degenerate (point) components
    # [affine-variables-1-end]
    test.assertTrue(u.size() == 3)
    test.assertTrue(u[0].itv() == Interval())
    test.assertTrue(v[0].itv() == Interval(1,2))
    test.assertTrue(w[0].itv() == Interval(1.5,1.5))

    # [affine-variables-2-beg]
    # Each component owns a dedicated, distinct noise symbol, exposed by
    # noise_index(i) on the AffineVariables container -- Affine itself has
    # no noise_index() of its own, since a plain component reference does
    # not know its own position in the vector.
    v = AffineVariables([[1,2],[-1,1]])
    i0 = v.noise_index(0) # 0
    i1 = v.noise_index(1) # 1
    v[0] = Interval(5,6)  # v[0] now encloses [5,6]...
    # [affine-variables-2-end]
    test.assertTrue(i0 == 0)
    test.assertTrue(i1 == 1)
    test.assertTrue(v.noise_index(0) == i0) # ...but keeps its own symbol
    test.assertTrue(v[0].itv() == Interval(5,6))

    # [affine-variables-assignment-beg]
    v = AffineVariables([[1.,2.],[-1.,1.]])

    i0 = v.noise_index(0)
    v[0] = Interval(5.,6.)
    v[1] = 3.

    # The numeric values change, but the declared identities remain.
    # [affine-variables-assignment-end]
    test.assertTrue(v.noise_index(0) == i0)
    test.assertTrue(v.noise_index(1) == 1)
    test.assertTrue(v[0].itv() == Interval(5.,6.))
    test.assertTrue(v[1].itv() == Interval(3.))

    # [affine-variables-dependency-beg]
    v = AffineVariables([[1.,2.]])
    a = Affine(v[0])
    b = Affine(v[0])
    correlated = a - b

    # These two plain Affine objects are created independently.
    p = Affine(Interval(1.,2.))
    q = Affine(Interval(1.,2.))
    independent = p - q
    # [affine-variables-dependency-end]
    test.assertTrue(correlated.itv() == Interval(0.0))
    test.assertTrue(independent.itv() == Approx_Interval(Interval(-1.,1.)))

    # [affine-variables-3-beg]
    # resize() drops any existing dependency information and reconstructs
    # every component as a fresh, independent, unbounded variable.
    # conservativeResize() instead keeps each surviving component's
    # interval enclosure (not its correlations), and initializes any newly
    # added component to ]-oo,+oo[.
    v = AffineVariables([[1,2],[-1,1],[3,4]])
    old_i0 = v.noise_index(0)

    v.conservativeResize(4)
    # v[0..2] still enclose their original intervals, v[3] is unbounded.
    # Their previous noise-symbol identities are intentionally not kept.
    new_i0 = v.noise_index(0)

    # resize() discards even the surviving interval enclosures and creates
    # a completely fresh set of unbounded variables.
    v.resize(2)
    # [affine-variables-3-end]
    test.assertTrue(v.size() == 2)
    test.assertTrue(v[0].itv() == Interval())
    test.assertTrue(old_i0 == new_i0)
    test.assertTrue(v.noise_index(0) == new_i0)
    test.assertTrue(v[1].itv() == Interval())

    # [affine-variables-4-beg]
    # init(Interval) broadcasts the same interval to every component,
    # while preserving each component's own dedicated noise symbol.
    v = AffineVariables(3)
    v.init(Interval(1,2))
    # [affine-variables-4-end]
    for i in range(v.size()):
      test.assertTrue(v[i].itv() == Interval(1,2))
      test.assertTrue(v.noise_index(i) == i)

    # [affine-variables-5-beg]
    # To use the declared variables in a computation, convert them to a
    # plain AffineVector (or let arithmetic operators do it implicitly):
    # the conversion goes through Affine's own copy constructor and
    # therefore preserves every noise symbol exactly.
    v = AffineVariables([[1,2],[-1,1],[3,4]])
    x = AffineVector(v) # safe: element-wise Affine copy construction
    # [affine-variables-5-end]
    test.assertTrue(x[0].itv() == v[0].itv())

    # [affine-variables-6-beg]
    # C++'s "w = v;" performs a whole-vector copy assignment, which calls
    # AffineVarMain's own assignment per component -- this copies each
    # component's noise-symbol identity along with its value, so w[i] and
    # v[i] end up sharing the exact same noise symbol. Python's "=" only
    # rebinds a name (it never calls a C++ assignment operator), so the
    # same observable effect is obtained here through the copy
    # constructor instead.
    v = AffineVariables([[1,2],[-1,1]])
    w = AffineVariables(2)
    w = AffineVariables(v)
    d = w[0] - v[0]
    # [affine-variables-6-end]
    test.assertTrue(w.noise_index(0) == v.noise_index(0))
    test.assertTrue(d.itv() == Interval(0.0))

    # A copy is therefore not a declaration of independent uncertainty:
    # w[0] and v[0] remain fully correlated.

    # [affine-variables-7-beg]
    # C++'s operator=(const IntervalVector&) reassigns every component of
    # the existing object in place, resizing it first if the sizes differ.
    # Python's "=" cannot invoke that operator either (see above): rebuilding
    # v from a new AffineVariables produces the same observable result --
    # every component recreated fresh from the given interval vector, no
    # dependency from before surviving.
    v = AffineVariables(2)
    v = AffineVariables([[1,2],[-1,1],[3,4]]) # v is (re)built with 3 components
    # [affine-variables-7-end]
    test.assertTrue(v.size() == 3)
    test.assertTrue(v[2].itv() == Interval(3,4))
    test.assertTrue(v.noise_index(2) == 2)

    # [affine-variables-8-beg]
    # Compound assignment, and re-binding a component or the whole vector
    # to a foreign Affine/AffineVector, are deliberately not exposed on
    # AffineVariables: they would either mutate a declared variable in
    # place, or hand it a value that fabricates a dependency it never had.
    # Attempting either of them raises a TypeError, since no matching
    # overload is bound for these operations.
    v = AffineVariables(1)
    with test.assertRaises(TypeError):
      v += 1.
    with test.assertRaises(TypeError):
      v.init(Affine(1.))
    # [affine-variables-8-end]

    # [affine-variables-9-beg]
    # AffineVariables can appear directly as an operand of a matrix
    # product: it is converted internally to an AffineVector, so the
    # result always carries Affine (not a bound-variable) coefficients.
    A = Matrix(2,2)
    A[0,0] = 1.; A[0,1] = 1.
    A[1,0] = 0.; A[1,1] = 1.
    v = AffineVariables([[-1,1],[-1,1]])
    y = A @ v
    # [affine-variables-9-end]
    test.assertTrue(y.itv() == IntervalVector([[-2,2],[-1,1]]))

    # [affine-variables-10-beg]
    # Unary minus also converts to a plain AffineVector: negating a
    # declared variable no longer represents "the" variable itself, but a
    # derived expression.
    v = AffineVariables([[1,2],[-1,1]])
    m = -v
    # [affine-variables-10-end]
    test.assertTrue(m.itv() == IntervalVector([[-2,-1],[-1,1]]))

    # [affine-variables-11-beg]
    v = AffineVariables([[1,2],[-1,1]])
    # Expected `print(v)` output:
    # [ [ 1, 2] : 1.5 + 0.5 eps_0 + 0 eps_1 + 0 [-1,1] ; 
    #   [-1, 1] : 0   + 0 eps_0   + 1 eps_1 + 0 [-1,1]  ]

    # [affine-variables-11-end]

    # [affine-variables-12-beg]
    # Linear combination of two DIFFERENT declared variables: both noise
    # symbols show up in the result, each with its own coefficient.
    v = AffineVariables([[1,2],[3,4]])
    z = v[0] + v[1]
    # Expected `print(z)` output:
    #   [4, 6] : 5 + 0.5 eps_0 + 0.5 eps_1 + 0 [-1,1]
    # [affine-variables-12-end]
    test.assertTrue(z.itv() == Interval(4,6))
    test.assertTrue(z.noise(0) == 0.5)
    test.assertTrue(z.noise(1) == 0.5)
    test.assertTrue(z.err() == 0.0)

    # [affine-variables-13-beg]
    # Linear combination with real coefficients and a constant term: the
    # constant only shifts the midpoint, each variable's own coefficient
    # is scaled independently.
    v = AffineVariables([[0,2],[0,4]])
    z = 2.*v[0] - 0.5*v[1] + 1.
    # Expected `print(z)` output:
    #   [-1, 5] : 2 + 2 eps_0 + -1 eps_1 + 0 [-1,1]
    # [affine-variables-13-end]
    test.assertTrue(z.itv() == Interval(-1,5))
    test.assertTrue(z.noise(0) == 2.0)
    test.assertTrue(z.noise(1) == -1.0)

    # [affine-variables-14-beg]
    # Reusing the same declared variable several times still collapses
    # exactly, whatever the linear combination (see also the AffineVector
    # page for the same property at vector level):
    # 4x - 2(x+x) + (x-x) = 4x - 4x + 0 = 0 for any x.
    v = AffineVariables([[-2,3]])
    x = v[0]
    z = 4.*x - 2.*(x + x) + (x - x)
    # Expected `print(z)` output:
    #   [0, 0] : 0 + 0 eps_0 + 0 [-1,1]
    # [affine-variables-14-end]
    test.assertTrue(z.itv() == Interval(0.0))
    test.assertTrue(z.noise(0) == 0.0)
    test.assertTrue(z.err() == 0.0)

    # [affine-variables-15-beg]
    # Nonlinear function of a single declared variable: sqr(x) = x*x is
    # computed by directly expanding the product rather than through a
    # generic Chebyshev linearization, so it stays exact here (the
    # remainder term only has to account for the eps_0^2 <= 1 cross term,
    # not for a derivative approximation).
    v = AffineVariables([[1,3]]) # x0 = 2, radius = 1
    x = v[0]
    y = sqr(x)
    # Expected `print(y)` output (hand-traced from the source: new
    # midpoint x0^2 + 0.5*coeff^2 = 4 + 0.5, new coefficient 2*x0*coeff = 4,
    # remainder 0.5 plus a floating-point safety margin far below display
    # precision):
    #   [0, 9] : 4.5 + 4 eps_0 + 0.5 [-1,1]
    # [affine-variables-15-end]
    test.assertTrue(y.itv().is_superset(Interval(1,9))) # true range of x^2 on [1,3]
    test.assertTrue(y == Approx_Affine(Interval(0,9),1.e-6))

    # [affine-variables-16-beg]
    # Transcendental functions (sin, cos, exp, ...) go through a
    # Chebyshev/MinRange linearization with its own floating-point error
    # tracking; unlike the linear and sqr() cases above, the exact printed
    # coefficients are not reproduced here (they were not hand-verified
    # against the source) -- only the guaranteed containment is checked.
    # Run this snippet to see the actual printed values.
    v = AffineVariables([[0, 1.5707963267948966]]) # [0, pi/2]
    y = sin(v[0])  #  [-0.0544991, 1.25618] : 0.600837 + 0.52725 eps_0 + 0.128086 [-1,1]
    # The exact coefficients depend on the selected linearization mode.
    # [affine-variables-16-end]
    test.assertTrue(y.itv().is_superset(Interval(0,1)))


  def tests_AffineVector_manual(test):

    # [affine-vector-1-beg]
    v = AffineVariables([[1,2],[-1,1],[3,4]])
    x = AffineVector(v)          # materialized from declared variables
    y = x + x                    # ordinary vector expression
    box = x.itv()
    # [affine-vector-1-end]
    test.assertTrue(box == IntervalVector([[1,2],[-1,1],[3,4]]))
    test.assertTrue(y.itv() == IntervalVector([[2,4],[-2,2],[6,8]]))

    # [affine-vector-2-beg]
    # Dependency is preserved through vector-level linear combinations
    # exactly as it is for a single Affine (see "Why affine arithmetic" in
    # the Affine class page): reusing the same declared variables several
    # times does not widen the result the way plain interval arithmetic
    # would.
    v = AffineVariables([[-2,3],[0,1]])
    x = AffineVector(v)
    z = 4.*x - 2.*(x + x) + (x - x)
    # [affine-vector-2-end]
    test.assertTrue(z.itv() == IntervalVector([[0,0],[0,0]]))

    # [affine-vector-3-beg]
    # Compound assignment operators are available on a plain AffineVector.
    v = AffineVariables([[1,2],[3,4]])
    x = AffineVector(v)
    x += 1.
    x *= 2.
    # [affine-vector-3-end]
    test.assertTrue(x.itv() == IntervalVector([[4,6],[8,10]]))

    # [affine-vector-4-beg]
    # AffineVector offers subvector(start_id, end_id), extracting the
    # (inclusive) range of components between the two given 0-based
    # indices -- the Python binding does not expose Eigen's own head() /
    # tail() / segment() names.
    v = AffineVariables([[1,2],[-1,1],[3,4],[5,6]])
    x = AffineVector(v)
    first_two = x.subvector(0,1)
    last_two  = x.subvector(2,3)
    middle    = x.subvector(1,2)
    # [affine-vector-4-end]
    test.assertTrue(first_two.itv() == IntervalVector([[1,2],[-1,1]]))
    test.assertTrue(last_two.itv()  == IntervalVector([[3,4],[5,6]]))
    test.assertTrue(middle.itv()    == IntervalVector([[-1,1],[3,4]]))

    # [affine-vector-5-beg]
    # transpose() turns an AffineVector into a plain AffineMatrix (here of
    # shape 1x2): unlike in C++, there is no implicit conversion to
    # AffineRow on assignment, since Python's "=" never triggers one.
    v = AffineVariables([[1,2],[-1,1]])
    x = AffineVector(v)
    r = x.transpose()
    # [affine-vector-5-end]
    test.assertTrue(r[0,0].itv() == Interval(1,2))
    test.assertTrue(r[0,1].itv() == Interval(-1,1))

    # [affine-vector-6-beg]
    # Applying an elementary function component-wise: there is no
    # vectorized shortcut for nonlinear functions (unlike the linear
    # compound assignments above), so it is done component by component.
    v = AffineVariables([[0,1],[1,2]])
    x = AffineVector(v)
    y = AffineVector(x.size())
    for i in range(x.size()):
      y[i] = sin(x[i])
    # [affine-vector-6-end]
    test.assertTrue(y[0].itv().is_superset(Interval(0, 0.8)))     # true range: [0, sin(1)]
    test.assertTrue(y[1].itv().is_superset(Interval(0.85, 0.99))) # true range: [sin(1), 1] (pi/2 in [1,2])

    # [affine-vector-7-beg]
    v = AffineVariables([[1,2],[3,4]])
    x = AffineVector(v)
    # Expected `print(x)` output (each component prints its own full
    # Affine representation, joined by " ; " and wrapped in "[ ... ]" --
    # this is why itv() is usually preferred for a quick summary):
    #   [ [1, 2] : 1.5 + 0.5 eps_0 + 0 eps_1 + 0 [-1,1]  ;
    #     [3, 4] : 3.5 + 0 eps_0 + 0.5 eps_1 + 0 [-1,1]  ]
    # [affine-vector-7-end]


  def tests_AffineMatrix_manual(test):

    # [affine-matrix-1-beg]
    A = Matrix(2,2)
    A[0,0] = 1.; A[0,1] = 1.
    A[1,0] = 0.; A[1,1] = 1.

    v = AffineVariables([[-1,1],[-1,1]])
    x = AffineVector(v)

    y = A @ x    # real matrix times a vector of affine forms
    box = y.itv()
    # [affine-matrix-1-end]
    test.assertTrue(box == IntervalVector([[-2,2],[-1,1]]))

    # [affine-matrix-2-beg]
    # Dependency preservation extends to matrix products: evaluating the
    # same declared variables through two matrices and recombining them
    # linearly does not lose correlation, which is typically where affine
    # arithmetic pays off compared to interval arithmetic in practice
    # (e.g. repeated evaluation of the same uncertain state vector through
    # several matrices, as in a state-space / control computation).
    A = Matrix(2,2)
    A[0,0] = 1.; A[0,1] = 0.
    A[1,0] = 0.; A[1,1] = 1.

    v = AffineVariables([[-1,1],[-1,1]])
    x = AffineVector(v)

    diff = (A @ x) - (A @ x)
    # [affine-matrix-2-end]
    test.assertTrue(diff.itv() == IntervalVector([[0,0],[0,0]]))

    # [affine-matrix-3-beg]
    # AffineMatrix is a plain Eigen matrix alias, like AffineVector;
    # itv() converts it to an IntervalMatrix.
    M = AffineMatrix(2,2)
    var = AffineVariables([[1,2],[0],[0],[1,2]])
    M[0,0] = var[0]
    M[0,1] = var[1]
    M[1,0] = var[2]
    M[1,1] = var[3]
    box = M.itv()
    # [affine-matrix-3-end]
    test.assertTrue(box[0,0] == Interval(1,2))
    test.assertTrue(box[0,1] == Interval(0,0))

    # [affine-matrix-shared-beg]
    v = AffineVariables(2)
    v[0] = Interval(-1,1)
    v[1] = Interval(2,3)

    M = AffineMatrix(2,2)
    M[0,0] = v[0]
    M[0,1] = v[1]
    M[1,0] = v[0]
    M[1,1] = v[1]
    # [affine-matrix-shared-end]

    # [affine-matrix-4-beg]
    # transpose(), like on any Eigen matrix, swaps rows and columns.
    M = AffineMatrix(2,2)
    M[0,0] = Affine(Interval(1,2)); M[0,1] = Affine(3.)
    M[1,0] = Affine(-1.);           M[1,1] = Affine(Interval(0,1))
    Mt = M.transpose()
    # [affine-matrix-4-end]
    test.assertTrue(Mt[0,1].itv() == M[1,0].itv())
    test.assertTrue(Mt[1,0].itv() == M[0,1].itv())

    # [affine-matrix-5-beg]
    # Row and column extraction return an AffineRow / AffineVector.
    M = AffineMatrix(2,2)
    M[0,0] = Affine(Interval(1,2)); M[0,1] = Affine(3.)
    M[1,0] = Affine(-1.);           M[1,1] = Affine(Interval(0,1))
    row0 = M.row(0)
    col1 = M.col(1)
    # [affine-matrix-5-end]
    test.assertTrue(row0[0].itv() == Interval(1,2))
    test.assertTrue(row0[1].itv() == Interval(3,3))
    test.assertTrue(col1[0].itv() == Interval(3,3))
    test.assertTrue(col1[1].itv() == Interval(0,1))

    # [affine-matrix-products-beg]
    v = AffineVariables([[-1,1],[2,3]])

    M = AffineMatrix(2,2)
    M[0,0] = v[0]
    M[0,1] = v[1]
    M[1,0] = v[0]
    M[1,1] = v[1]

    x = AffineVector(v)
    y1 = M @ x

    N = AffineMatrix.eye(2,2)
    y2 = N @ x
    P = M @ N
    # [affine-matrix-products-end]
    test.assertTrue(y1.size() == 2)
    test.assertTrue(y2.size() == 2)
    test.assertTrue(P.rows() == 2)
    test.assertTrue(P.cols() == 2)

    # [affine-matrix-6-beg]
    M = AffineMatrix(1,2)
    M[0,0] = Affine(Interval(1,2))
    M[0,1] = Affine(3.)
    # Expected `print(M)` format (rows separated by newlines, each
    # component printing its own full Affine representation -- verbose,
    # as for AffineVector; itv() is usually preferred for a quick view).
    # Both M[0,0] and M[0,1] were built as freestanding Affine(Interval)/
    # Affine(double) constants (not through AffineVariables), so neither
    # carries an eps_i term: all uncertainty sits in the error bound.
    #   [ [1, 2] : 1.5 + 0.5 [-1,1]  , <3, 3> : 3 + 0 [-1,1]  ]
    # [affine-matrix-6-end]


if __name__ ==  '__main__':
  unittest.main()
