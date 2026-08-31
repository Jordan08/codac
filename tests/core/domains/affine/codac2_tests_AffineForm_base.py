#!/usr/bin/env python

# Codac tests
#
# Base tests for Affine (the Python binding of AffineMain<T>): construction,
# copy/assignment-like semantics, comparisons and conversions to/from
# Interval, and the numerical robustness of the underlying affine
# representation (coefficient storage, compact(), products/sums under
# extreme scales).
#
# Python has no assignment operator distinct from construction: wherever
# the C++ test reassigns an *existing* Affine from an Interval or a double
# (relying on AffineMain::operator=, e.g. to check that noise_count() is
# preserved), this file uses x.init(...) instead, since
# AffineMain::init(const Interval&) is implemented as "*this = x;" and is
# therefore the exact same operation reachable from Python. Reassigning an
# existing Affine from another Affine has no such equivalent; those spots
# fall back to copy-construction (Affine(source)).
#
# TODO: check this with Simon Rohou
# AffineVariables.__getitem__ returns a live reference to the container's
# element (like Affine/AffineVector/AffineRow/AffineMatrix indexing
# already did), not a copy: wherever this file names a single component
# (e.g. "value = variables[0]") and then mutates that name in place
# (+=, .compact(), .inflate(), ...), it wraps the read in Affine(...) to
# get an independent copy first -- otherwise the mutation would silently
# reach back into the container, exactly as it would with a live C++
# reference, and corrupt any later read of that same component.
#
# ----------------------------------------------------------------------------
#  \date       2026
#  \author     Jordan Ninin
#  \copyright  Copyright 2026 Codac Team
#  \license    GNU Lesser General Public License (LGPL)

import unittest
import sys
import math
from codac import *

ERROR = sys.float_info.epsilon*10


class TestAffineFormBase(unittest.TestCase):

  def check_change_mode_minrange(self):
    Affine.change_mode(Affine.Affine_Mode.AF_Lin_MinRange)
    self.assertTrue(Affine.Affine_Mode.AF_Lin_MinRange == Affine.get_mode())
    Affine.change_mode(Affine.Affine_Mode.AF_Lin_Chebyshev)
    self.assertTrue(Affine.Affine_Mode.AF_Lin_Chebyshev == Affine.get_mode())

  def check_affine_inclu(self, y_actual, y_expected):
    if y_expected.is_empty():
      self.assertTrue(y_actual.is_empty())
      return
    self.assertFalse(y_actual.is_empty())
    self.assertTrue(y_expected.lb() >= y_actual.itv().lb())
    self.assertTrue(y_expected.ub() <= y_actual.itv().ub())

  def check_interval_enclosure(self, actual, expected):
    self.assertTrue(actual.itv().is_superset(expected))

  def test_approx_affinemain(self):

    # empty
    ax = AffineVariables(1)
    ax[0] = Interval.empty()
    self.assertTrue(ax[0] == Approx_Affine(Interval.empty()))

    # unbounded
    ax = AffineVariables(1)
    self.assertTrue(ax[0] == Approx_Affine(Interval()))

    ax = AffineVariables(IntervalVector([[0,1],[0,1],[0,1]]))
    self.assertTrue(ax[0] == Approx_Affine(Interval(0,1)))
    self.assertTrue(ax[1] == Approx_Affine(Interval(0,1)))
    self.assertTrue(ax[2] == Approx_Affine(Interval(0,1)))

    # simple linear cases: x in [0,1], 1 noise var
    ax = AffineVariables(IntervalVector([[0,1]]))
    self.assertTrue(ax[0] == Approx_Affine(Interval(0,1)))
    self.assertTrue((ax[0]+1.0) == Approx_Affine(Interval(1,2)))
    self.assertTrue((-ax[0]) == Approx_Affine(Interval(-1,0)))
    self.assertTrue((2.0*ax[0]) == Approx_Affine(Interval(0,2)))

    # two affine variables summed together (2 noise variables)
    ax = AffineVariables(IntervalVector([[0,1],[0,1]]))
    y = ax[0] + ax[1]
    self.assertTrue(y == Approx_Affine(Interval(0,2)))

    # must not match a wrong expected interval
    bx = AffineVariables(2)
    bx[0] = Interval(0,1)
    self.assertFalse(bx[0] == Approx_Affine(Interval(5,6)))
    bx[0] = Interval(5,6)
    self.assertTrue(bx[0] == Approx_Affine(Interval(5,6)))

  def test_affineform_operations(self):

    self.check_change_mode_minrange()

    # static_assert-only checks (is_constructible_v, is_convertible_v, ...)
    # are not translated: they check C++ compile-time properties.

    # operator-()
    ax = AffineVariables(1)
    ax[0] = Interval(0,1)
    self.assertTrue(-(ax[0]) == Approx_Affine(Interval(-1,0), ERROR))

    ax = AffineVariables(1)
    self.assertTrue(-(ax[0]) == Approx_Affine(Interval(), ERROR))

    ax = AffineVariables(1)
    ax[0] = Interval(-oo,0)
    self.assertTrue(-(ax[0]) == Approx_Affine(Interval(0,oo), ERROR))

    ax = AffineVariables(4)
    ax[0] = Interval(-oo,1)
    self.assertTrue(-(ax[0]) == Approx_Affine(Interval(-1,oo), ERROR))
    self.assertTrue(ax[2] == Approx_Affine(Interval(), ERROR))

    x = Affine(3.)
    y = Affine(x)
    self.assertTrue(y.itv() == Interval(3.))

    z = Affine(x) # Python has no operator=(const Affine&) distinct from
                   # construction; copy-constructing is the closest
                   # equivalent to "z = x;" reachable from the bindings.
    self.assertTrue(z.itv() == Interval(3.))
    self.assertTrue((-z).itv() == Interval(-3.))

    xv = AffineVariables(IntervalVector([[1.,2.],[3.,4.],[5.,6.]]))
    source = Affine(xv[0]) # a copy: set_empty() below must not empty xv[0]
    source.set_empty()
    destination = Affine(source) # see note above about operator=
    self.assertTrue(destination.is_empty())
    self.assertTrue(destination.noise_count() == source.noise_count())

    # AF_fAF2 (the raw representation struct) type-trait checks are not
    # translated: AF_fAF2 is not exposed to Python as its own type.

    xv = AffineVariables(2)
    xv[0] = Interval(5.0)
    self.assertTrue(xv[0].itv() == Interval(5.))

    variables = AffineVariables(IntervalVector([[1., 2.]]))
    x = variables[0]
    residual = floor(x) - x
    self.assertTrue(residual.itv() != Interval(0.))
    self.assertTrue(Interval(-1., 0.).is_subset(residual.itv()))

  def test_af_inactive_assignment_keeps_storage_consistent(self):

    # destination initially owns storage for one affine variable
    small = AffineVariables(IntervalVector([[0.0, 1.0]]))
    destination = small[0]

    # source is inactive and has a larger logical dimension. Python has no
    # operator= reusing destination's existing storage (only construction
    # is reachable), so this only exercises the resulting value, not the
    # "does not corrupt stale storage" regression itself.
    large_inactive = AffineVariables(16)
    self.assertFalse(large_inactive[0].is_active())

    destination = Affine(large_inactive[0])
    self.assertTrue(destination.noise_count() == 16)
    self.assertTrue(destination.is_unbounded())

    # Re-activating the destination writes all coefficients from 0 through
    # noise_count(). init(x) is AffineMain::operator=(const Interval&).
    destination.init(Interval(-2.0, 3.0))

    self.assertTrue(destination.noise_count() == 16)
    self.assertTrue(destination.is_active())
    self.check_affine_inclu(destination, Interval(-2.0, 3.0))

    for i in range(destination.noise_count()):
      self.assertTrue(destination.noise(i) == 0.0)

  def test_af_compact_safe_for_every_inactive_status(self):

    variables = AffineVariables(12)
    value = Affine(variables[0])

    self.assertFalse(value.is_active())
    value.compact()

    value.init(Interval.empty())
    self.assertTrue(value.is_empty())
    value.compact()

    value.init(Interval())
    self.assertTrue(value.is_unbounded())
    value.compact()

    value.init(Interval(-oo, 2.0))
    self.assertTrue(value.is_unbounded())
    value.compact()

    value.init(Interval(-2.0, oo))
    self.assertTrue(value.is_unbounded())
    value.compact()

  def test_af_dimension_changes_preserve_coefficients_and_inclusion(self):

    # Public arithmetic currently exposes dimension growth, not
    # contraction. This exercises the observable grow path used by
    # operator+= and operator*=.
    small = AffineVariables(IntervalVector([[1.0, 2.0]]))
    large = AffineVariables(IntervalVector([
      [-1.0, 1.0], [-2.0, 2.0], [-3.0, 3.0], [-4.0, 4.0],
      [-5.0, 5.0], [-6.0, 6.0], [-7.0, 7.0], [-8.0, 8.0]
    ]))

    # Affine(...) copies: small[0]/large[7] are live references into their
    # container (see codac2_py_AffineVariables.cpp), and small[0] is read
    # again below, so mutating it in place here would corrupt that read.
    total = Affine(small[0])
    total += large[7]
    self.assertTrue(total.noise_count() == 8)
    self.check_affine_inclu(total, small[0].itv() + large[7].itv())

    product = Affine(small[0])
    product *= large[7]
    self.assertTrue(product.noise_count() == 8)
    self.check_affine_inclu(product, small[0].itv() * large[7].itv())

  def test_af_multiplication_handles_opposite_finite_scales(self):

    inputs = [
      (Interval(1.e299, 1.e300), Interval(1.e-300, 1.e-299)),
      (Interval(-1.e300, -1.e299), Interval(1.e-300, 1.e-299)),
      (Interval(1.e200, 2.e200), Interval(1.e-200, 2.e-200)),
      (Interval(-1.e200, 2.e200), Interval(-2.e-200, 1.e-200))
    ]

    for first, second in inputs:
      variables = AffineVariables(IntervalVector([first, second]))

      reference = first * second
      result = variables[0] * variables[1]

      self.assertFalse(result.is_empty())
      self.check_affine_inclu(result, reference)

      if result.is_active():
        self.assertTrue(math.isfinite(result.err()))
        self.assertTrue(result.err() >= 0.0)

  def test_af_scalar_products_and_sums_handle_extreme_cancellation(self):

    large = 1.e300
    small = 1.e-300

    variables = AffineVariables(IntervalVector([[large, next_float(large)]]))

    scaled = variables[0] * small
    self.assertFalse(scaled.is_empty())
    self.check_affine_inclu(scaled, variables[0].itv() * small)

    cancelled = variables[0] + (-large)
    self.assertFalse(cancelled.is_empty())
    self.check_affine_inclu(cancelled, variables[0].itv() - large)

    if scaled.is_active():
      self.assertTrue(math.isfinite(scaled.err()))
      self.assertTrue(scaled.err() >= 0.0)
    if cancelled.is_active():
      self.assertTrue(math.isfinite(cancelled.err()))
      self.assertTrue(cancelled.err() >= 0.0)

  def test_af_active_remainder_stays_finite_and_nonnegative(self):

    variables = AffineVariables(IntervalVector([[-2.0, 3.0], [0.5, 2.0], [-1.0, 1.0]]))

    results = [
      variables[0] + variables[1],
      variables[0] * variables[1],
      sqr(variables[0]),
      exp(variables[2]),
      sin(variables[0]),
      asinh(variables[0]),
      inv(variables[1])
    ]

    for result in results:
      if result.is_active():
        self.assertTrue(math.isfinite(result.err()))
        self.assertTrue(result.err() >= 0.0)

  def test_af_repeated_multiplication_and_status_transitions_remain_stable(self):

    # Stress/regression test for repeated active/inactive transitions.
    for iteration in range(256):
      variables = AffineVariables(IntervalVector([[-1.0, 2.0], [2.0, 3.0], [-0.5, 0.5]]))

      value = Affine(variables[0]) # a copy: variables[0] is read again below
      value *= variables[1]
      value += variables[2]
      self.check_affine_inclu(
        value,
        variables[0].itv() * variables[1].itv() + variables[2].itv()
      )

      value.init(Interval())
      self.assertTrue(value.is_unbounded())
      value.init(Interval(-1.0, 1.0))
      self.check_affine_inclu(value, Interval(-1.0, 1.0))
      value.compact()

  def test_af_interval_extraction_matches_all_status_transitions(self):

    value = Affine()

    value.init(Interval.empty())
    self.assertTrue(value.itv().is_empty())

    value.init(Interval())
    self.assertTrue(value.itv() == Interval())

    value.init(Interval(-oo, 4.0))
    self.assertTrue(value.itv() == Interval(-oo, 4.0))

    value.init(Interval(-3.0, oo))
    self.assertTrue(value.itv() == Interval(-3.0, oo))

    value.init(Interval(2.0))
    self.assertTrue(value.itv() == Interval(2.0))

    value.init(Interval(-2.0, 5.0))
    self.check_affine_inclu(value, Interval(-2.0, 5.0))

  def test_af_square_preserves_small_coefficients_amplified_by_large_remainder(self):

    cases = [
      (Interval(-1.e-20, 1.e-20), Interval(-1.e20, 1.e20)),
      (Interval(-1.e-30, 1.e-30), Interval(-1.e30, 1.e30)),
      (Interval(-1.e-50, 1.e-50), Interval(-1.e50, 1.e50))
    ]

    for variable, uncertainty in cases:
      variables = AffineVariables(IntervalVector([variable]))

      value = Affine(variables[0])
      value += uncertainty

      before = value.itv()
      result = sqr(value)
      reference = sqr(before)

      self.assertFalse(result.is_empty())
      self.check_affine_inclu(result, reference)

      if result.is_active():
        self.assertTrue(math.isfinite(result.err()))
        self.assertTrue(result.err() >= 0.0)

  def test_af_normalized_square_preserves_mixed_scale_contributions(self):

    scales = [1.e20, 1.e30, 1.e50, 1.e100]

    for scale in scales:
      inverse_scale = 1.0/scale

      variables = AffineVariables(IntervalVector([[-inverse_scale, inverse_scale]]))

      value = Affine(variables[0])
      value += Interval(-scale, scale)

      reference = sqr(value.itv()) / (scale*scale)
      normalized = sqr(value) / (scale*scale)

      self.assertFalse(normalized.is_empty())
      self.check_affine_inclu(normalized, reference)

      if normalized.is_active():
        self.assertTrue(math.isfinite(normalized.err()))
        self.assertTrue(normalized.err() >= 0.0)

  def test_af_scalar_multiplication_transfers_discarded_coefficients(self):

    variables = AffineVariables(IntervalVector([[-1.e-20, 1.e-20]]))

    result = variables[0] * 0.5

    self.check_affine_inclu(result, variables[0].itv() * 0.5)

  def test_af_addition_transfers_cancelled_coefficients(self):

    variables = AffineVariables(IntervalVector([[-1.e-20, 1.e-20]]))

    x = variables[0]

    # Introduces a nearly complete cancellation while retaining a
    # representable small residual coefficient.
    result = x + (-x * prev_float(1.0))
    reference = x.itv() + (-x.itv() * prev_float(1.0))

    self.check_affine_inclu(result, reference)

  def test_af_compact_transfers_removed_coefficients(self):

    variables = AffineVariables(IntervalVector([[-1.e-8, 1.e-8]]))

    value = Affine(variables[0])
    before = value.itv()

    value.compact(1.e-6)

    self.check_affine_inclu(value, before)

    if value.is_active():
      self.assertTrue(value.err() >= 1.e-8)

  def test_af_square_handles_small_coefficient_with_large_remainder(self):

    small = 0.5 * 2**-55 # inférieur à AF_EC
    remainders = [1.0, 2.0, 10.0, 1.e4, 1.e8, 1.e16]

    for remainder in remainders:
      variables = AffineVariables(IntervalVector([[-small, small]]))

      value = Affine(variables[0])

      # Le petit coefficient reste dans _val[1]. L'incertitude ajoutée est
      # placée dans _err.
      value += Interval(-remainder, remainder)

      self.assertTrue(value.is_active())
      self.assertTrue(value.noise_count() == 1)

      self.assertTrue(math.fabs(value.noise(0)) > 0.0)
      self.assertTrue(math.fabs(value.noise(0)) < 2**-55)
      self.assertTrue(value.err() >= remainder)

      input_enclosure = value.itv()
      expected = sqr(input_enclosure)

      result = sqr(value)

      self.assertFalse(result.is_empty())
      self.check_affine_inclu(result, expected)

      if result.is_active():
        self.assertTrue(math.isfinite(result.err()))
        self.assertTrue(result.err() >= 0.0)

  def test_af_square_handles_several_small_coefficients_with_large_remainder(self):

    small = 0.25 * 2**-55

    variables = AffineVariables(IntervalVector([
      [-small, small], [-small, small], [-small, small], [-small, small]
    ]))

    value = variables[0] + variables[1] + variables[2] + variables[3]
    value += Interval(-1.e8, 1.e8)

    self.assertTrue(value.is_active())

    before = value.itv()
    expected = sqr(before)
    result = sqr(value)

    self.assertFalse(result.is_empty())
    self.check_affine_inclu(result, expected)

    if result.is_active():
      self.assertTrue(math.isfinite(result.err()))
      self.assertTrue(result.err() >= 0.0)

  def test_af_normalized_square_remains_inclusive_for_mixed_scales(self):

    cases = [
      (0.5 * 2**-55, 2.0),
      (0.5 * 2**-55, 10.0),
      (0.25 * 2**-55, 1.e4),
      (0.125 * 2**-55, 1.e8)
    ]

    for small, remainder in cases:
      variables = AffineVariables(IntervalVector([[-small, small]]))

      value = Affine(variables[0])
      value += Interval(-remainder, remainder)

      normalization = remainder * remainder

      self.assertTrue(math.isfinite(normalization))
      self.assertTrue(normalization > 0.0)

      expected = sqr(value.itv()) / normalization
      result = sqr(value) / normalization

      self.assertFalse(result.is_empty())
      self.check_affine_inclu(result, expected)

      if result.is_active():
        self.assertTrue(math.isfinite(result.err()))
        self.assertTrue(result.err() >= 0.0)

  def test_affinemain_empty_factory_matches_interval_empty_assignment(self):

    from_factory = Affine.empty()
    from_assignment = Affine(Interval.empty())

    self.assertTrue(from_factory.is_empty())
    self.assertTrue(from_assignment.is_empty())
    self.assertTrue(from_factory.itv().is_empty())
    self.assertTrue(from_factory.noise_count() == from_assignment.noise_count())
    self.assertTrue(from_factory == Approx_Affine(Interval.empty()))

  def test_strict_comparison_operators_return_documented_boolinterval(self):

    ax = AffineVariables(2)
    ax[0] = Interval(0, 1)
    ax[1] = Interval(2, 3)

    # ub(ax[0]) < lb(ax[1]): certainly true, both ways and for every
    # Interval/Affine combination of the operands.
    self.assertTrue((ax[0] < ax[1]) == BoolInterval.TRUE)
    self.assertTrue((ax[1] > ax[0]) == BoolInterval.TRUE)
    self.assertTrue((ax[0] < Interval(2, 3)) == BoolInterval.TRUE)
    self.assertTrue((Interval(2, 3) > ax[0]) == BoolInterval.TRUE)

    # The reverse comparisons are certainly false.
    self.assertTrue((ax[1] < ax[0]) == BoolInterval.FALSE)
    self.assertTrue((ax[0] > ax[1]) == BoolInterval.FALSE)
    self.assertTrue((Interval(2, 3) < ax[0]) == BoolInterval.FALSE)

    # Overlapping ranges: the sign is undetermined.
    bx = AffineVariables(2)
    bx[0] = Interval(0, 2)
    bx[1] = Interval(1, 3)
    self.assertTrue((bx[0] < bx[1]) == BoolInterval.UNKNOWN)
    self.assertTrue((bx[0] > bx[1]) == BoolInterval.UNKNOWN)

    # An empty operand makes the comparison empty.
    ex = AffineVariables(1)
    ex[0] = Interval.empty()
    self.assertTrue((ex[0] < ax[1]) == BoolInterval.EMPTY)
    self.assertTrue((ex[0] > ax[1]) == BoolInterval.EMPTY)
    self.assertTrue((ex[0] < Interval(2, 3)) == BoolInterval.EMPTY)
    self.assertTrue((Interval(2, 3) < ex[0]) == BoolInterval.EMPTY)

  def test_set_relation_predicates_delegate_to_interval_enclosure(self):

    ax = AffineVariables(4)
    ax[0] = Interval(1, 2)   # subset candidate
    ax[1] = Interval(0, 3)   # superset candidate
    ax[2] = Interval(5, 6)   # disjoint from ax[0]/ax[1]
    ax[3] = Interval(2, 4)   # overlaps ax[1] without being a subset

    # is_subset / is_superset, against both an Interval and an Affine.
    self.assertTrue(ax[0].is_subset(Interval(0, 3)))
    self.assertTrue(ax[0].is_subset(ax[1]))
    self.assertFalse(ax[1].is_subset(ax[0]))
    self.assertTrue(ax[1].is_superset(ax[0]))
    self.assertTrue(ax[1].is_superset(Interval(1, 2)))
    self.assertFalse(ax[0].is_superset(ax[1]))

    # Strict variants fail on touching boundaries, succeed strictly inside.
    self.assertTrue(ax[0].is_strict_subset(ax[1]))
    self.assertFalse(ax[1].is_strict_subset(ax[1]))
    self.assertTrue(ax[1].is_strict_superset(ax[0]))
    self.assertFalse(ax[1].is_strict_superset(ax[1]))

    # Interior variant: ax[0]=[1,2] lies in the interior of [0,3], but its
    # own lower bound touches the boundary of [1,3].
    self.assertTrue(ax[0].is_interior_subset(Interval(0, 3)))
    self.assertFalse(ax[0].is_interior_subset(Interval(1, 3)))
    self.assertTrue(ax[0].is_strict_interior_subset(Interval(0, 3)))

    # contains / interior_contains.
    self.assertTrue(ax[0].contains(1.0))
    self.assertTrue(ax[0].contains(2.0))
    self.assertFalse(ax[0].interior_contains(1.0))
    self.assertTrue(ax[0].interior_contains(1.5))

    # intersects / overlaps / is_disjoint, on overlapping and disjoint pairs.
    self.assertTrue(ax[0].intersects(ax[1]))
    self.assertTrue(ax[0].overlaps(ax[1]))
    self.assertFalse(ax[0].is_disjoint(ax[1]))

    self.assertFalse(ax[0].intersects(ax[2]))
    self.assertFalse(ax[0].overlaps(ax[2]))
    self.assertTrue(ax[0].is_disjoint(ax[2]))

    # A boundary-touching pair intersects (non-empty intersection) but does
    # not overlap (zero-volume intersection), and is not disjoint.
    touching = AffineVariables(2)
    touching[0] = Interval(0, 2)
    touching[1] = Interval(2, 4)
    self.assertTrue(touching[0].intersects(touching[1]))
    self.assertFalse(touching[0].overlaps(touching[1]))
    self.assertFalse(touching[0].is_disjoint(touching[1]))

    # Every predicate must agree exactly with the equivalent Interval query
    # on the enclosure, since Affine simply forwards to itv().
    self.assertTrue(ax[3].is_subset(ax[1]) == ax[3].itv().is_subset(ax[1].itv()))
    self.assertTrue(ax[3].overlaps(ax[1]) == ax[3].itv().overlaps(ax[1].itv()))

  def test_repr_streams_interval_enclosure_and_coefficients(self):

    ax = AffineVariables(1)
    ax[0] = Interval(1.0, 2.0)
    s = str(ax[0])
    expected_prefix = str(ax[0].itv())

    self.assertTrue(s[:len(expected_prefix)] == expected_prefix)
    self.assertTrue(s.find(" : ") != -1)
    self.assertTrue(s.find("eps_0") != -1)

    empty_form = Affine.empty()
    self.assertTrue(str(empty_form).find("not enabled") != -1)

  def test_size_init_and_init_from_list_reinitialize_in_place(self):

    # size() is always 1, matching Interval::size()'s template-uniformity role.
    self.assertTrue(Affine().size() == 1)
    self.assertTrue(Affine.empty().size() == 1)
    self.assertTrue(Affine(Interval(1.0, 2.0)).size() == 1)

    ax = AffineVariables(1)
    ax[0] = Interval(1.0, 2.0)
    x = Affine(ax[0])
    self.assertTrue(x.size() == 1)

    # init() resets to [-oo,oo].
    x.init()
    self.assertTrue(x.itv() == Interval())
    self.assertTrue(x.is_unbounded())

    # init(const Interval&) assigns like operator=(const Interval&): for a
    # plain Affine (with no dedicated noise symbol of its own) the whole
    # uncertainty is carried by the remainder error term, and every noise
    # coefficient is left at zero.
    x.init(Interval(3.0, 7.0))
    self.assertTrue(x.itv() == Interval(3.0, 7.0))
    self.assertTrue(x.mid() == 5.0)
    self.assertTrue(x.err() == 2.0)
    for i in range(x.noise_count()):
      self.assertTrue(x.noise(i) == 0.0)

    # A degenerate interval carries no uncertainty at all.
    y = Affine()
    y.init(Interval(4.0))
    self.assertTrue(y.itv() == Interval(4.0))
    self.assertTrue(y.err() == 0.0)
    self.assertTrue(y.is_degenerated())

    # init_from_list mirrors Interval::init_from_list: one value gives a
    # degenerate point, two values give the bounds of the resulting hull.
    single = Affine()
    single.init_from_list([4.0])
    self.assertTrue(single.itv() == Interval(4.0))
    self.assertTrue(single.is_degenerated())

    pair = Affine()
    pair.init_from_list([1.0, 3.0])
    self.assertTrue(pair.itv() == Interval(1.0, 3.0))
    self.assertTrue(pair.mid() == 2.0)
    self.assertTrue(pair.err() == 1.0)

  def test_equality_with_double(self):

    a = Affine(4.0)

    self.assertTrue(a == 4.0)
    self.assertFalse(a == 5.0)
    self.assertTrue(a != 5.0)
    self.assertFalse(a != 4.0)

  def test_implicit_conversion_to_and_from_interval(self):

    a = Affine(4.0)

    # Python has no implicit "operator Interval()"; itv() is the explicit
    # equivalent (and is, in fact, exactly what that operator forwards to).
    as_interval = a.itv()
    self.assertTrue(as_interval == Interval(4.0))

    from_interval = Affine(Interval(2.0, 3.0)) # implicit Affine(const Interval&)
    self.check_interval_enclosure(from_interval, Interval(2.0, 3.0))

    # Plain arithmetic between an Affine and an Interval still resolves
    # cleanly to an Affine.
    itv = Interval(1.0, 2.0)
    total = a + itv
    self.check_interval_enclosure(total, Interval(4.0) + itv)


if __name__ ==  '__main__':
  unittest.main()
