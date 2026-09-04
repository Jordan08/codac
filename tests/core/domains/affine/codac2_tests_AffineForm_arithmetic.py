#!/usr/bin/env python

# Codac tests
#
# Tests for Affine (the Python binding of AffineMain<T>) arithmetic
# operators (+,-,*,/) and elementary/transcendental functions (sqrt, log,
# exp, pow, root, trigonometric and hyperbolic functions, floor/ceil/sign,
# min/max/intersection/union), including a sampling-based regression suite
# checking that every nonlinear function returns a sound affine enclosure
# of its reference interval image.
#
# Python has no operator overloading distinct from construction, and
# "from codac import *" shadows several Python builtins (min, max, abs,
# pow, ...) with their codac counterparts: this file captures the plain
# Python min/max under different names before the star-import, and uses
# math.fabs instead of the (now Affine/Interval-only) abs() whenever a
# plain float is involved. It also explicitly wraps floor()/ceil()/
# integer() results with Affine(...): these free functions return an
# Interval even when given an Affine argument (mirroring the C++ API),
# and C++ silently promotes that back through the implicit
# AffineMain(const Interval&) constructor wherever the result is bound to
# an AffineT-typed variable -- a promotion Python has no equivalent for.
#
# AffineVariables.__getitem__ returns a live reference to the container's
# element, not a copy: wherever this file names a single component (e.g.
# "x_aff = ax[0]") and later mutates that name in place via a compound
# operator, it wraps the read in Affine(...) first to get an independent
# copy -- otherwise the mutation would reach back into the container and
# corrupt any later read of that same component (this notably applies to
# the check_add/check_mul/check_div helpers below, which reuse their
# operands across many assertions).
#
# ----------------------------------------------------------------------------
#  \date       2026
#  \author     Jordan Ninin
#  \copyright  Copyright 2026 Codac Team
#  \license    GNU Lesser General Public License (LGPL)

import unittest
import sys
import math
import builtins
_min = builtins.min
_max = builtins.max
from codac import *

MAX_DOUBLE = sys.float_info.max
SAMPLE_SIZE = 300
ERROR = sys.float_info.epsilon*10

piL = Interval.pi().lb()
piU = Interval.pi().ub()


def safe_math(f, x):
  """Calls a math.* function and returns nan instead of raising, mirroring
  how the C++ standard library reports a domain/range error (e.g.
  std::asin outside [-1,1], std::cosh on a huge input) by returning NaN
  rather than throwing."""
  try:
    return f(x)
  except (ValueError, OverflowError):
    return float('nan')

def odd_root_value(x, n):
  return math.copysign(math.fabs(x) ** (1.0/n), x)


class TestAffineFormArithmetic(unittest.TestCase):

  def check_affine_inclu(self, y_actual, y_expected):
    if y_expected.is_empty():
      self.assertTrue(y_actual.is_empty())
      return
    self.assertFalse(y_actual.is_empty())
    self.assertTrue(y_expected.lb() >= y_actual.itv().lb())
    self.assertTrue(y_expected.ub() <= y_actual.itv().ub())

  def check_add(self, x, z, y_expected):
    # Mirrors both the (Interval,double) and the (Interval,Interval)
    # overloads of CHECK_add<T> in the C++ source: when z is itself an
    # Interval it is given its own noise symbol in a shared 2-variable
    # context (as the C++ overload that forwards through two Affine
    # operands does); otherwise it stays a plain double.
    if isinstance(z, Interval):
      ax = AffineVariables(2)
      ax[0] = x
      ax[1] = z
      x_aff = ax[0]
      z_op = ax[1]
      active_guard = x_aff.is_active() and z_op.is_active()
    else:
      ax = AffineVariables(1)
      ax[0] = x
      x_aff = ax[0]
      z_op = z
      active_guard = x_aff.is_active() and math.fabs(z) < 1.e100

    y_actual = x_aff + z_op
    self.assertTrue(y_actual == Approx_Affine(y_expected, ERROR))

    # symmetrical case
    y_actual = z_op + x_aff
    self.assertTrue(y_actual == Approx_Affine(y_expected, ERROR))

    # +=operator: Affine(...) makes an independent copy to mutate here --
    # plain "y_actual = x_aff" would alias x_aff itself (Python assignment
    # never copies), corrupting it for the many checks reusing it below.
    y_actual = Affine(x_aff)
    y_actual += z_op
    self.assertTrue(y_actual == Approx_Affine(y_expected, ERROR))

    # +=operator in the other direction
    y_actual = Affine(z_op)
    y_actual += x_aff
    self.assertTrue(y_actual == Approx_Affine(y_expected, ERROR))

    # subtraction
    y_actual = (-x_aff) - z_op
    self.assertTrue(y_actual == Approx_Affine(-y_expected, ERROR))

    # symmetrical case
    y_actual = (-z_op) - x_aff
    self.assertTrue(y_actual == Approx_Affine(-y_expected, ERROR))

    # -=operator
    y_actual = -x_aff
    y_actual -= z_op
    self.assertTrue(y_actual == Approx_Affine(-y_expected, ERROR))

    # -=operator in the other direction
    y_actual = -z_op
    y_actual -= x_aff
    self.assertTrue(y_actual == Approx_Affine(-y_expected, ERROR))

    # linear simplification
    if active_guard:
      y_actual = Affine(z_op)
      y_actual += x_aff
      y_actual -= x_aff
      y_actual += x_aff
      self.assertTrue(y_actual == Approx_Affine(y_expected, ERROR))

      y_actual = Affine(x_aff)
      y_actual += z_op
      y_actual -= z_op
      y_actual += z_op
      self.assertTrue(y_actual == Approx_Affine(y_expected, ERROR))

  def check_add_scalar(self, x, z, y_expected):
    self.check_add(x, z, y_expected)
    self.check_add(-x, -z, -y_expected)

  def check_mul(self, x, z, y_expected):
    ax = AffineVariables(2)
    ax[0] = x
    ax[1] = z
    x_aff = ax[0]
    z_aff = ax[1]

    y_actual = x_aff * z_aff
    self.check_affine_inclu(y_actual, y_expected)

    y_actual = z_aff * x_aff # symmetrical case
    self.check_affine_inclu(y_actual, y_expected)

    y_actual = Affine(x_aff) # a copy: mutating y_actual below must not
                              # mutate the shared x_aff (Python "=" aliases)
    y_actual *= z_aff # *=operator
    self.check_affine_inclu(y_actual, y_expected)

  def check_mul_scal(self, x, z, y_expected):
    xa = AffineVariables(1)
    xa[0] = x
    self.assertTrue(xa[0]*z == Approx_Affine(y_expected, ERROR))
    self.assertTrue(z*xa[0] == Approx_Affine(y_expected, ERROR))

  def check_div(self, x, z, y_expected):
    ax = AffineVariables(2)
    ax[0] = x
    ax[1] = z
    x_aff = ax[0]
    z_aff = ax[1]

    y_actual = x_aff / z_aff
    self.check_affine_inclu(y_actual, y_expected)

    y_actual = Affine(x_aff)
    y_actual /= z_aff # /=operator
    self.check_affine_inclu(y_actual, y_expected)

  def check_div_scal(self, x, z, y_expected):
    xa = AffineVariables(1)
    xa[0] = x
    y_actual = xa[0] / z
    self.check_affine_inclu(y_actual, y_expected)

    # a copy: mutating y_actual below must not mutate xa[0] itself.
    y_actual = Affine(xa[0])
    y_actual /= z # /=operator
    self.check_affine_inclu(y_actual, y_expected)

  def check_sqrt(self, x, y_expected):
    ax = AffineVariables(1)
    ax[0] = x
    self.check_affine_inclu(sqrt(ax[0]), y_expected)

  def check_log(self, x, y_expected):
    ax = AffineVariables(1)
    ax[0] = x
    self.check_affine_inclu(log(ax[0]), y_expected)

  def check_exp(self, x, y_expected):
    ax = AffineVariables(1)
    ax[0] = x
    self.check_affine_inclu(exp(ax[0]), y_expected)

  def check_trigo(self, x, expected):
    ax = AffineVariables(1)
    ax[0] = x
    self.check_affine_inclu(sin(ax[0]), expected)
    self.check_affine_inclu(sin(Interval.pi()-ax[0]), expected)
    self.check_affine_inclu(sin(ax[0]+Interval.two_pi()), expected)
    self.check_affine_inclu(sin(-ax[0]), -expected)
    self.check_affine_inclu(cos(ax[0]-Interval.half_pi()), expected)
    self.check_affine_inclu(cos(Interval.half_pi()-ax[0]), expected)
    self.check_affine_inclu(cos(ax[0]+Interval.half_pi()), -expected)
    self.check_affine_inclu(cos(ax[0]+Interval.two_pi()-Interval.half_pi()), expected)

  def check_pow(self, x, p, y_expected):
    ax = AffineVariables(1)
    ax[0] = x
    self.check_affine_inclu(pow(ax[0], p), y_expected)
    self.check_affine_inclu(pow(-ax[0], p), y_expected if p % 2 == 0 else -y_expected)

  def check_root(self, x, p, y_expected):
    ax = AffineVariables(1)
    ax[0] = x
    self.check_affine_inclu(root(ax[0], p), y_expected)

  def check_sqr(self, x, y_expected):
    ax = AffineVariables(1)
    ax[0] = x
    self.check_affine_inclu(sqr(ax[0]), y_expected)

  def check_asin(self, x, y_expected):
    ax = AffineVariables(1)
    ax[0] = x
    self.check_affine_inclu(asin(ax[0]), y_expected)

    if not x.is_empty() and not x.is_unbounded() and not y_expected.is_empty() and not y_expected.is_unbounded():
      variables = AffineVariables(IntervalVector([x]))
      y = asin(variables[0])
      for k in range(10000):
        d = x.lb() + (x.ub() - x.lb()) * k / 10000.0
        dd = safe_math(math.asin, d)
        if not math.isnan(dd):
          self.assertTrue(y.contains(dd))

  def check_acos(self, x, y_expected):
    ax = AffineVariables(1)
    ax[0] = x
    self.check_affine_inclu(acos(ax[0]), y_expected)

    if not x.is_empty() and not x.is_unbounded() and not y_expected.is_empty() and not y_expected.is_unbounded():
      variables = AffineVariables(IntervalVector([x]))
      y = acos(variables[0])
      for k in range(10000):
        d = x.lb() + (x.ub() - x.lb()) * k / 10000.0
        dd = safe_math(math.acos, d)
        if not math.isnan(dd):
          self.assertTrue(y.contains(dd))

  def check_atan(self, x, y_expected):
    ax = AffineVariables(1)
    ax[0] = x
    self.check_affine_inclu(atan(ax[0]), y_expected)

    if not x.is_empty() and not x.is_unbounded() and not y_expected.is_empty() and not y_expected.is_unbounded():
      variables = AffineVariables(IntervalVector([x]))
      y = atan(variables[0])
      for k in range(10000):
        d = x.lb() + (x.ub() - x.lb()) * k / 10000.0
        dd = safe_math(math.atan, d)
        if not math.isnan(dd):
          self.assertTrue(y.contains(dd))

  def check_tanh(self, x):
    y_expected = tanh(x)
    ax = AffineVariables(1)
    ax[0] = x
    self.check_affine_inclu(tanh(ax[0]), y_expected)

    if not x.is_empty() and not x.is_unbounded() and not y_expected.is_empty() and not y_expected.is_unbounded():
      variables = AffineVariables(IntervalVector([x]))
      y = tanh(variables[0])
      for k in range(10000):
        d = x.lb() + (x.ub() - x.lb()) * k / 10000.0
        dd = safe_math(math.tanh, d)
        if not math.isnan(dd):
          self.assertTrue(y.contains(dd))

  def check_sinh(self, x):
    y_expected = sinh(x)
    ax = AffineVariables(1)
    ax[0] = x
    self.check_affine_inclu(sinh(ax[0]), y_expected)

    if not x.is_empty() and not x.is_unbounded() and not y_expected.is_empty() and not y_expected.is_unbounded():
      variables = AffineVariables(IntervalVector([x]))
      y = sinh(variables[0])
      for k in range(10000):
        d = x.lb() + (x.ub() - x.lb()) * k / 10000.0
        dd = safe_math(math.sinh, d)
        if not math.isnan(dd):
          self.assertTrue(y.contains(dd))

  def check_cosh(self, x):
    y_expected = cosh(x)
    ax = AffineVariables(1)
    ax[0] = x
    self.check_affine_inclu(cosh(ax[0]), y_expected)

    if not x.is_empty() and not x.is_unbounded() and not y_expected.is_empty() and not y_expected.is_unbounded():
      variables = AffineVariables(IntervalVector([x]))
      y = cosh(variables[0])
      for k in range(10000):
        d = x.lb() + (x.ub() - x.lb()) * k / 10000.0
        dd = safe_math(math.cosh, d)
        if not math.isnan(dd):
          self.assertTrue(y.contains(dd))

  def check_sampled_enclosure(self, input, output, function, sample_count=SAMPLE_SIZE):
    self.assertFalse(output.is_empty())
    self.assertTrue(sample_count > 100)

    for k in range(sample_count+1):
      if k == 0:
        x = input.lb()
      elif k == sample_count:
        x = input.ub()
      else:
        ratio = k / sample_count
        x = (1.0-ratio)*input.lb() + ratio*input.ub()
        x = _max(input.lb(), _min(input.ub(), x))

      self.assertTrue(input.contains(x))
      expected = function(x)
      self.assertTrue(output.contains(expected))

  def check_sampled_interval_enclosure(self, input, output, global_reference, interval_function, sample_count=SAMPLE_SIZE):
    self.assertFalse(output.is_empty())
    self.assertTrue(sample_count > 100)

    output_interval = output.itv()
    tested_point_count = 0

    for k in range(sample_count+1):
      if k == 0:
        x = input.lb()
      elif k == sample_count:
        x = input.ub()
      else:
        ratio = k / sample_count
        x = (1.0-ratio)*input.lb() + ratio*input.ub()
        x = _max(input.lb(), _min(input.ub(), x))

      self.assertTrue(input.contains(x))

      # Both operands are certified enclosures of the same point image:
      # interval_function([x]) is the local singleton evaluation, while
      # global_reference encloses the image of the complete valid domain.
      # Their intersection remains a certified enclosure and removes
      # harmless overshoots such as cos([0]) extending by a few ULP
      # beyond 1.
      local_reference = interval_function(Interval(x))

      # An empty singleton image means x itself is outside the real
      # pointwise domain (e.g. x=0 for log). Such a point imposes no
      # pointwise inclusion condition.
      if local_reference.is_empty():
        continue

      expected = local_reference & global_reference
      self.assertFalse(expected.is_empty())
      self.assertTrue(expected.is_subset(output_interval))
      tested_point_count += 1

    # Prevent a wrongly specified domain or interval function from making
    # the whole sampling test vacuous.
    self.assertTrue(tested_point_count > 100)

  def check_nonlinear_operation(self, operation_name, input, valid_domain, affine_function, interval_function, sample_count=SAMPLE_SIZE):
    variables = AffineVariables(IntervalVector([input]))
    result = affine_function(variables[0])
    reference = interval_function(valid_domain)

    self.check_affine_inclu(result, reference)
    if not valid_domain.is_empty() and not valid_domain.is_degenerated() and not valid_domain.is_unbounded():
      # Use the certified interval image of each singleton as the oracle:
      # a scalar math.* result may differ by a few ULP from the
      # outward-rounded CODAC interval implementation.
      self.check_sampled_interval_enclosure(valid_domain, result, reference, interval_function, sample_count)

  def check_nonlinear_interval_image(self, operation_name, input, affine_function, interval_function):
    variables = AffineVariables(IntervalVector([input]))
    result = affine_function(variables[0])
    self.check_affine_inclu(result, interval_function(input))

  def test_arithmetic_operators_and_elementary_functions(self):

    self.check_add(Interval.empty(), Interval(0,1), Interval.empty())
    self.check_add(Interval(0,1), Interval.empty(), Interval.empty())
    self.check_add(Interval(-oo,1), Interval(0,1), Interval(-oo, 2))
    self.check_add(Interval(1,oo), Interval(0,1), Interval(1,oo))
    self.check_add(Interval(-oo,oo), Interval(0,1), Interval(-oo,oo))
    self.check_add(Interval(MAX_DOUBLE,oo), 1, Interval(MAX_DOUBLE,oo))
    self.check_add(Interval(MAX_DOUBLE,oo), -1, Interval(prev_float(MAX_DOUBLE),oo))
    self.check_add(Interval(MAX_DOUBLE,oo), Interval(MAX_DOUBLE,oo), Interval(MAX_DOUBLE,oo))
    self.check_add(Interval(MAX_DOUBLE,oo), -oo, Interval.empty())
    self.check_add(Interval(MAX_DOUBLE,oo), oo, Interval.empty())
    self.check_add(Interval(-oo,-MAX_DOUBLE), 1, Interval(-oo,next_float(-MAX_DOUBLE)))
    self.check_add(Interval(-oo,-MAX_DOUBLE), -1, Interval(-oo,-MAX_DOUBLE))
    self.check_add(Interval(-oo,-MAX_DOUBLE), Interval(-oo,-MAX_DOUBLE), Interval(-oo,-MAX_DOUBLE))

    self.check_add_scalar(Interval.empty(), oo, Interval.empty())
    self.check_add_scalar(Interval.empty(), 0, Interval.empty())
    self.check_add_scalar(Interval(0,1), 1, Interval(1,2))
    self.check_add_scalar(Interval(0,1), -oo, Interval.empty())
    self.check_add_scalar(Interval(0,1), oo, Interval.empty())
    self.check_add_scalar(Interval(-oo,1), 1, Interval(-oo,2))

    self.check_mul(Interval.empty(), Interval(0,1), Interval.empty())
    self.check_mul(Interval.zero(), Interval(), Interval.zero())
    self.check_mul(Interval(-1,1), Interval(-oo,0), Interval())
    self.check_mul(Interval(-oo,-1), Interval(-1,0), Interval(0,oo))
    self.check_mul(Interval(-oo, 1), Interval(-1,0), Interval(-1,oo))
    self.check_mul(Interval(0, 1), Interval(1,oo), Interval(0,oo))
    self.check_mul(Interval(0, 1), Interval(-1,oo), Interval(-1,oo))
    self.check_mul(Interval(-oo,-1), Interval(0,1), Interval(-oo,0))
    self.check_mul(Interval(-oo, 1), Interval(0,1), Interval(-oo,1))
    self.check_mul(Interval(0, 1), Interval(-oo,-1), Interval(-oo,0))
    self.check_mul(Interval(0, 1), Interval(-oo,1), Interval(-oo,1))
    self.check_mul(Interval(1,oo), Interval(0,1), Interval(0,oo))
    self.check_mul(Interval(-1,oo), Interval(0,1), Interval(-1,oo))
    self.check_mul(Interval(1,2), Interval(1,2), Interval(1,4))
    self.check_mul(Interval(1,2), Interval(-2,3), Interval(-4,6))
    self.check_mul(Interval(-1,1), Interval(-2,3), Interval(-3,3))
    self.check_mul_scal(Interval(1,2), -oo, Interval.empty())
    self.check_mul_scal(Interval(1,2), oo, Interval.empty())
    self.check_mul_scal(Interval(1,2), -1, Interval(-2,-1))

    a = AffineVariables(2)   # contexte à 2 variables
    a[0] = Interval(1.,2.)
    b = AffineVariables(4)   # contexte à 4 variables
    b[3] = Interval(3.,4.)

    p = a[0] * b[3]
    self.assertTrue(p.itv().is_superset(Interval(1.,2.)*Interval(3.,4.)))

    self.check_div(Interval.empty(), Interval(0,1), Interval.empty())
    self.check_div(Interval.zero(), Interval.zero(), Interval.empty())
    self.check_div(Interval(1,2), Interval.zero(), Interval.empty())
    self.check_div(Interval(), Interval.zero(), Interval.empty())
    self.check_div(Interval.zero(), Interval(0,1), Interval.zero())
    self.check_div(Interval.zero(), Interval(), Interval.zero())
    self.check_div(Interval(6,12), Interval(2,3), Interval(2,6))
    self.check_div(Interval(6,12), Interval(-3,-2), Interval(-6,-2))
    self.check_div(Interval(6,12), Interval(-2,3), Interval())
    self.check_div(Interval(-oo,-1), Interval(-1,0), Interval(1,oo))
    self.check_div(Interval(-oo,-1), Interval(0,1), Interval(-oo,-1))
    self.check_div(Interval(1,oo), Interval(-1,0), Interval(-oo,-1))
    self.check_div(Interval(1,oo), Interval(0,1), Interval(1,oo))
    self.check_div(Interval(-1,1), Interval(-1,1), Interval())
    self.check_div_scal(Interval(1,2), -oo, Interval.empty())
    self.check_div_scal(Interval(1,2), oo, Interval.empty())
    self.check_div_scal(Interval(1,2), -1, Interval(-2,-1))

    self.check_log(Interval.empty(), Interval.empty())
    self.check_log(Interval(), Interval())
    self.check_log(Interval(0,oo), Interval())
    self.check_log(Interval(-oo,0), Interval.empty())
    self.check_log(Interval(1,2), Interval(0,math.log(2)))
    self.check_log(Interval(-1,1), Interval(-oo,0))

    ax = AffineVariables(1)
    ax[0] = Interval(0,next_float(0))
    self.assertTrue(log(ax[0]).itv().ub() > -744.5)

    self.check_log(Interval(0,1), Interval(-oo,0))
    self.check_log(Interval(1,oo), Interval(0,oo))
    self.check_log(Interval(0), Interval.empty())
    self.check_log(Interval(-2,-1), Interval.empty())

    self.check_exp(Interval.empty(), Interval.empty())
    self.check_exp(Interval(), Interval(0,oo))
    self.check_exp(Interval(0,oo), Interval(1,oo))
    self.check_exp(Interval(-oo,0), Interval(0,1))
    self.check_exp(Interval(0,2), Interval(1,math.exp(2)))
    self.check_exp(Interval(-1,1), Interval(math.exp(-1),math.exp(1)))
    self.check_exp(Interval(1.e100,1.e111), Interval(MAX_DOUBLE,oo))
    self.check_exp(Interval(MAX_DOUBLE,oo), Interval(MAX_DOUBLE,oo))
    self.check_exp(Interval(0, MAX_DOUBLE), Interval(1,oo))

    av = AffineVariables(IntervalVector([[-8.]]))
    self.assertTrue(root(av[0], 3).itv() == Approx(Interval(-2.), ERROR))
    self.assertTrue(root(av[0], 5).itv() == Approx(Interval(-(8.**0.2)), ERROR))

    av[0] = Interval(1., 5.)
    self.assertTrue(root(av[0], 0).itv() == Interval.empty())
    self.assertTrue(root(av[0].itv(), 0) == Interval.empty())

    ax = AffineVariables(1)
    ax[0] = Interval(-8., 8.)
    y = root(ax[0], 3)
    self.assertTrue(y.itv().is_superset(Interval(-2., 2.))) # root(-8,3)=-2, root(8,3)=2
    # vérifie que la forme reste affine (pas juste un intervalle plat) :
    self.assertFalse(y.is_degenerated())

    self.check_trigo(Interval(), Interval(-1,1))
    self.check_trigo(Interval.empty(), Interval.empty())
    self.check_trigo(Interval(0,piU/2.0), Interval(0,1))
    self.check_trigo(Interval(0,piU), Interval(0,1))
    self.check_trigo(Interval(0,3*piU/2.0), Interval(-1,1))
    self.check_trigo(Interval(piL,3*piU/2.0), Interval(-1,0))
    self.check_trigo(Interval(0.5,1.5), Interval(math.sin(0.5),math.sin(1.5)))
    self.check_trigo(Interval(1.5,3), Interval(math.sin(3.0),1))
    self.check_trigo(Interval(3,4), Interval(math.sin(4.0),math.sin(3.0)))
    self.check_trigo(Interval(3,5), Interval(-1,math.sin(3.0)))
    self.check_trigo(Interval(3,2*piU+1.5), Interval(-1,math.sin(1.5)))
    self.check_trigo(Interval(5,2*piU+1.5), Interval(math.sin(5.0),math.sin(1.5)))
    self.check_trigo(Interval(5,2*piU+3), Interval(math.sin(5.0),1))

    ax = AffineVariables(1)
    self.check_affine_inclu(tan(ax[0]), Interval())

    ax = AffineVariables(1)
    ax[0] = (-Interval.pi()/4.0 | Interval.pi()/4.0)
    self.check_affine_inclu(tan(ax[0]), Interval(-1,1))

    # tan(pi/4,pi/2)=[1,+oo)
    x = AffineVariables(1) # upper bound of x is close to pi/2
    x[0] = Interval(piL/4.0,(1-1e-10)*piL/2.0)
    y = tan(x[0])
    self.assertTrue(y.itv().lb() <= 1.0)
    self.assertTrue(y.itv().ub() > 1.e8) # upper bound of tan(x) is close to +oo

    # tan(-pi/2,pi/4)=(-oo,1]
    ax = AffineVariables(1)
    ax[0] = Interval(-(1-1e-10)*piL/2.0,piL/4.0)
    y = (tan(ax[0])).itv()
    self.assertTrue(y.lb() <= -1.e8) # lower bound is close to -oo
    self.assertTrue(y.ub() >= 1.0)

    ax = AffineVariables(1)
    ax[0] = Interval.pi()/2.0
    self.check_affine_inclu(tan(ax[0]), Interval())

    ax = AffineVariables(1)
    ax[0] = Interval.pi()
    y = tan(-ax[0])
    self.assertTrue(y == Approx_Affine(Interval(0), ERROR))
    self.assertTrue(y.contains(0.))
    self.assertTrue(y.diam() < 1e-8)

    ax = AffineVariables(1)
    ax[0] = Interval(3*piL/4.0, 5*piU/4.0)
    y = tan(ax[0])
    self.check_affine_inclu(y, Interval(-1,1))
    self.assertTrue(y.lb() == Approx(-1.0,1.e-5))
    self.assertTrue(y.ub() == Approx(1.0,1.e-5))

    ax[0] = Interval(-oo,oo)
    y = tan(ax[0])
    self.assertTrue(y == Approx_Affine(Interval(), ERROR))

    ax[0] = ((-Interval.pi()/4.0)|(Interval.pi()/4.0))
    y = tan(ax[0])
    self.check_affine_inclu(y, Interval(-1,1))
    self.assertTrue(y.lb() == Approx(-1.0,1.e-5))
    self.assertTrue(y.ub() == Approx(1.0,1.e-5))

    # tan(pi/4,pi/2)=[1,+oo)
    ax = AffineVariables(1)
    ax[0] = Interval(piL/4.0,(1-1e-10)*piL/2.0)
    y = tan(ax[0])
    self.assertTrue(y.lb() <= 1.0)
    self.assertTrue(y.ub() >= 1.e8) # upper bound of tan(x) is close to +oo

    # tan(-pi/2,pi/4)=(-oo,1]
    ax = AffineVariables(1)
    ax[0] = Interval(-(1-1e-10)*piL/2.0,piL/4.0)
    y = tan(ax[0])
    self.check_affine_inclu(y, Interval(-1.e8,1))

    ax = AffineVariables(1)
    ax[0] = (3*Interval.pi()/4.0 | 5*Interval.pi()/4.0)
    self.check_affine_inclu(tan(ax[0]), Interval(-1,1))

    self.check_pow(Interval(), 4, Interval(0,oo))
    self.check_pow(Interval.empty(), 4, Interval.empty())
    self.check_pow(Interval(2,3), 4, Interval(16,81))
    self.check_pow(Interval(-2,3), 4, Interval(0,81))
    self.check_pow(Interval(-3,2), 4, Interval(0,81))
    self.check_pow(Interval(2,oo), 4, Interval(16,oo))
    self.check_pow(Interval(), 3, Interval())
    self.check_pow(Interval.empty(), 3, Interval.empty())
    self.check_pow(Interval(2,3), 3, Interval(8,27))
    self.check_pow(Interval(-2,3), 3, Interval(-8,27))
    self.check_pow(Interval(-3,2), 3, Interval(-27,8))
    self.check_pow(Interval(2,oo), 3, Interval(8,oo))
    self.check_pow(Interval(-10,10), -2, Interval(1.0/100,oo))

    self.check_root(Interval(0,1), -1, Interval(1.0,oo))
    self.check_root(Interval(-27,-8), 3, Interval(-3,-2))
    self.check_root(Interval(-4,1), 2, Interval(0,1))
    self.check_root(Interval(-8,1), 3, Interval(-2,1))

    self.check_sqrt(Interval(), Interval(0,oo))
    self.check_sqrt(Interval(-oo,0), Interval.zero())
    self.check_sqrt(Interval(-9,4), Interval(0,2))
    self.check_sqrt(Interval(4,9), Interval(2,3))
    self.check_sqrt(Interval(-9,-4), Interval.empty())
    self.check_sqrt(Interval(-9,oo), Interval(0,oo))

    self.check_sinh(Interval())
    self.check_sinh(Interval(0,oo))
    self.check_sinh(Interval(0,1))
    self.check_sinh(Interval(1,oo))
    self.check_sinh(Interval(1,1))
    self.check_sinh(Interval(2,3))
    self.check_sinh(Interval(4,5))
    self.check_sinh(Interval(0,0.5))
    self.check_sinh(Interval(-0.2,0.5))
    self.check_sinh(Interval(0.6,0.9))
    self.check_sinh(Interval(-0.5,-0.3))

    self.check_cosh(Interval())
    self.check_cosh(Interval(0,oo))
    self.check_cosh(Interval(0,1))
    self.check_cosh(Interval(1,oo))
    self.check_cosh(Interval(1,1))
    self.check_cosh(Interval(2,3))
    self.check_cosh(Interval(4,5))
    self.check_cosh(Interval(0,0.5))
    self.check_cosh(Interval(-0.2,0.5))
    self.check_cosh(Interval(0.6,0.9))
    self.check_cosh(Interval(-0.5,-0.3))

    self.check_tanh(Interval.empty())
    self.check_tanh(Interval())
    self.check_tanh(Interval(0,oo))
    self.check_tanh(Interval(-oo,0))
    self.check_tanh(Interval(0,0.5))
    self.check_tanh(Interval(-0.2,0.5))
    self.check_tanh(Interval(0.6,0.9))
    self.check_tanh(Interval(-0.5,-0.3))

    # issue 248
    aff = AffineVariables(1)
    aff[0] = Interval(-1.57079632679489678, 1.1780972450961728626)
    self.assertFalse(tan(aff[0]).is_empty())

  def test_additional_operations(self):

    self.check_sqr(Interval.empty(), Interval.empty())
    self.check_sqr(Interval(), Interval(0,oo))
    self.check_sqr(Interval(2,3), Interval(4,9))
    self.check_sqr(Interval(-2,3), Interval(0,9))
    self.check_sqr(Interval(-3,2), Interval(0,9))
    self.check_sqr(Interval(2,oo), Interval(4,oo))

    self.check_asin(Interval.empty(), Interval.empty())
    self.check_asin(Interval(2,3), Interval.empty()) # entirely out of domain
    self.check_asin(Interval(0,1), Interval(0,Interval.half_pi().ub()))
    self.check_asin(Interval(-1,0), Interval(-Interval.half_pi().ub(),0))
    self.check_asin(Interval(-2,0.5), Interval(-Interval.half_pi().ub(),math.asin(0.5))) # truncated on the left

    for x_itv in [Interval(-0.95, 0.20), Interval(-0.75, 0.85), Interval(0.10, 0.99)]:
      variables = AffineVariables(IntervalVector([x_itv]))
      y = asin(variables[0])

      for k in range(10000):
        x = x_itv.lb() + (x_itv.ub() - x_itv.lb()) * k / 10000.0
        self.assertTrue(y.contains(math.asin(x)))

    self.check_acos(Interval.empty(), Interval.empty())
    self.check_acos(Interval(2,3), Interval.empty()) # entirely out of domain
    self.check_acos(Interval(-1,1), Interval(0,Interval.pi().ub()))
    self.check_acos(Interval(0,1), Interval(0,Interval.half_pi().ub()))

    self.check_atan(Interval.empty(), Interval.empty())
    self.check_atan(Interval(), Interval(-Interval.half_pi().ub(),Interval.half_pi().ub()))
    self.check_atan(Interval(0,1), Interval(0,math.atan(1.)))
    self.check_atan(Interval(-1,0), Interval(-math.atan(1.),0))

    ax = AffineVariables(1)
    ax[0] = Interval(2., 3.)
    self.assertTrue(acosh(ax[0]).itv().is_superset(Interval(math.acosh(2.), math.acosh(3.))))

    ax = AffineVariables(1)
    ax[0] = Interval(-1., 1.)
    self.assertTrue(asinh(ax[0]).itv().is_superset(Interval(math.asinh(-1.), math.asinh(1.))))

    ax = AffineVariables(1)
    ax[0] = Interval(-0.5, 0.5)
    self.assertTrue(atanh(ax[0]).itv().is_superset(Interval(math.atanh(-0.5), math.atanh(0.5))))

    # min, max, operator&, operator| between two affine forms (and between
    # an affine form and a plain Interval)
    ax = AffineVariables(2)
    ax[0] = Interval(1,3)
    ax[1] = Interval(2,5)

    self.assertTrue(min(ax[0],ax[1]) == Interval(1,3))
    self.assertTrue(max(ax[0],ax[1]) == Interval(2,5))
    self.assertTrue((ax[0] & ax[1]) == Interval(2,3))
    self.assertTrue((ax[0] | ax[1]) == Interval(1,5))

    # disjoint intervals: intersection is empty, union is the hull
    ax = AffineVariables(2)
    ax[0] = Interval(0,1)
    ax[1] = Interval(2,3)

    self.assertTrue((ax[0] & ax[1]).is_empty())
    self.assertTrue((ax[0] | ax[1]) == Interval(0,3))

    # mixed operand: Affine and plain Interval, both operand orders
    ax = AffineVariables(1)
    ax[0] = Interval(1,4)
    b = Interval(2,6)

    self.assertTrue((ax[0] & b) == Interval(2,4))
    self.assertTrue((b & ax[0]) == Interval(2,4))
    self.assertTrue((ax[0] | b) == Interval(1,6))
    self.assertTrue((b | ax[0]) == Interval(1,6))
    self.assertTrue(min(ax[0],b) == Interval(1,4))
    self.assertTrue(max(ax[0],b) == Interval(2,6))
    self.assertTrue(min(b,ax[0]) == Interval(1,4))
    self.assertTrue(max(b,ax[0]) == Interval(2,6))

  def test_critical_nonlinear_regressions(self):

    # Racine impaire d'un singleton négatif
    x = AffineVariables(IntervalVector([[-8.]]))
    self.assertTrue(root(x[0], 3).itv() == Approx(Interval(-2.), ERROR))

    # Indice zéro
    x = AffineVariables(IntervalVector([[1., 2.]]))
    self.assertTrue(root(x[0], 0).is_empty())

    # floor ne doit pas conserver une fausse dépendance
    variables = AffineVariables(IntervalVector([[1., 2.]]))
    x = variables[0]
    residual = floor(x) - x

    self.assertFalse(residual.itv() == Interval(0.))
    self.assertTrue(Interval(-1., 0.).is_subset(residual.itv()))

  def test_large_and_exceptional_arguments(self):

    int_max = 2147483647.0 # std::numeric_limits<int>::max(), 32-bit int

    inputs = [
      Interval(1e100, 1e100 + 1e90),
      Interval(-1e100, -1e100 + 1e90),
      Interval(int_max, int_max + 2.0),
      Interval(-oo, oo),
      Interval.empty()
    ]

    for input in inputs:
      x = AffineVariables(IntervalVector([input]))

      sin(x[0])
      cos(x[0])
      floor(x[0])
      ceil(x[0])

  def test_floor_and_ceil_of_noninteger_singleton(self):

    variables = AffineVariables(IntervalVector([[1.5], [-1.5], [2.0]]))

    self.assertTrue(floor(variables[0]) == Interval(1.0))
    self.assertTrue(ceil(variables[0]) == Interval(2.0))

    self.assertTrue(floor(variables[1]) == Interval(-2.0))
    self.assertTrue(ceil(variables[1]) == Interval(-1.0))

    self.assertTrue(floor(variables[2]) == Interval(2.0))
    self.assertTrue(ceil(variables[2]) == Interval(2.0))

  def test_odd_affine_roots_crossing_zero_enclose_sampled_values(self):

    inputs = [
      Interval(-8.0, 27.0),
      Interval(-1.e-12, 1.e6),
      Interval(-1.e6, 1.e-12),
      Interval(-1.e-300, 1.e-300),
      Interval(-1.0, 1024.0),
      Interval(-125.0, 32.0)
    ]
    orders = [3, 5, 7, 9, 101]

    for input in inputs:
      for order in orders:
        variables = AffineVariables(IntervalVector([input]))
        result = root(variables[0], order)
        reference = root(input, order)

        self.check_affine_inclu(result, reference)
        self.check_sampled_enclosure(
          input, result, lambda x, order=order: odd_root_value(x, order))

  def test_odd_affine_roots_crossing_zero_preserve_symmetry(self):

    orders = [3, 5, 7, 9]

    for order in orders:
      input = Interval(-64.0, 64.0)
      variables = AffineVariables(IntervalVector([input]))
      result = root(variables[0], order)

      self.assertTrue(result.contains(0.0))
      self.assertTrue(result.lb() == Approx(-result.ub(), 1.e-10))
      self.check_sampled_enclosure(
        input, result, lambda x, order=order: odd_root_value(x, order))

  def test_partially_valid_domains_keep_sound_affine_enclosure(self):

    # sqrt with a negative invalid part
    input = Interval(-9.0, 4.0)
    valid_domain = input & Interval(0.0, oo)
    variables = AffineVariables(IntervalVector([input]))
    result = sqrt(variables[0])

    self.check_affine_inclu(result, sqrt(valid_domain))
    self.check_sampled_enclosure(valid_domain, result, math.sqrt)

    # asin with values below its domain
    input = Interval(-2.0, 0.75)
    valid_domain = input & Interval(-1.0, 1.0)
    variables = AffineVariables(IntervalVector([input]))
    result = asin(variables[0])

    self.check_affine_inclu(result, asin(valid_domain))
    self.check_sampled_enclosure(valid_domain, result, math.asin)

    # acos with values above its domain
    input = Interval(-0.75, 2.0)
    valid_domain = input & Interval(-1.0, 1.0)
    variables = AffineVariables(IntervalVector([input]))
    result = acos(variables[0])

    self.check_affine_inclu(result, acos(valid_domain))
    self.check_sampled_enclosure(valid_domain, result, math.acos)

    # acosh with values below its domain
    input = Interval(0.0, 4.0)
    valid_domain = input & Interval(1.0, oo)
    variables = AffineVariables(IntervalVector([input]))
    result = acosh(variables[0])

    self.check_affine_inclu(result, acosh(valid_domain))
    self.check_sampled_enclosure(valid_domain, result, math.acosh)

  def test_entirely_invalid_nonlinear_domains_remain_empty(self):

    variables = AffineVariables(IntervalVector([[-9.0, -4.0], [2.0, 3.0], [-3.0, 0.5]]))

    self.assertTrue(sqrt(variables[0]).is_empty())
    self.assertTrue(asin(variables[1]).is_empty())
    self.assertTrue(acos(variables[1]).is_empty())
    self.assertTrue(acosh(variables[2]).is_empty())

  def test_all_nonlinear_affine_functions_enclose_reference_values(self):

    modes = [Affine.Affine_Mode.AF_Lin_Chebyshev, Affine.Affine_Mode.AF_Lin_MinRange]

    for mode in modes:
      Affine.change_mode(mode)

      # Reciprocal: both connected components of R\{0}, including points
      # close to the singularity.
      for input in [Interval(1.e-8, 3.0), Interval(-4.0, -1.e-8),
                    Interval(1.e-150, 1.e-140), Interval(-1.e8, -1.e-4)]:
        self.check_nonlinear_operation("inv", input, input, inv, lambda x: 1.0/x)

      # Square and absolute value, including intervals crossing zero.
      for input in [Interval(-4.0, 3.0), Interval(-1.e-12, 1.e12),
                    Interval(-3.0, -0.25), Interval(0.25, 5.0)]:
        self.check_nonlinear_operation("sqr", input, input, sqr, sqr)
        self.check_nonlinear_operation("abs", input, input, abs, abs)

      # Square root, including a large invalid negative part.
      for input in [Interval(0.0, 4.0), Interval(1.e-24, 1.e6),
                    Interval(-1.e12, 1.e-12), Interval(-1.e300, 4.0)]:
        domain = input & Interval(0.0, oo)
        self.check_nonlinear_operation("sqrt", input, domain, sqrt, sqrt, 8192)

      # Exponential and logarithm. The chosen exponential intervals avoid a
      # floating overflow in the scalar reference while still spanning very
      # different scales.
      for input in [Interval(-20.0, 20.0), Interval(-700.0, -600.0),
                    Interval(-1.e-12, 1.e-12), Interval(1.0, 10.0)]:
        self.check_nonlinear_operation("exp", input, input, exp, exp)
      for input in [Interval(1.e-300, 1.e-12), Interval(1.e-12, 1.e12),
                    Interval(0.25, 4.0), Interval(-1.e6, 3.0)]:
        domain = input & Interval(0.0, oo)
        self.check_nonlinear_operation("log", input, domain, log, log, 8192)

      # Integer and real powers, plus odd roots crossing zero.
      for input in [Interval(-8.0, 27.0), Interval(-1.e6, 1.e-12), Interval(-3.0, 2.0)]:
        self.check_nonlinear_operation(
          "pow(x,3)", input, input, lambda x: pow(x,3), lambda x: pow(x,3))
        self.check_nonlinear_operation(
          "root(x,3)", input, input, lambda x: root(x,3), lambda x: root(x,3), 8192)
      for input in [Interval(1.e-12, 4.0), Interval(0.25, 1.e6)]:
        self.check_nonlinear_operation(
          "pow(x,1.5)", input, input, lambda x: pow(x,1.5), lambda x: pow(x,1.5))

      # Interval and affine exponents deliberately lose the exponent's
      # affine dependency, but their interval image must still be enclosed.
      for input in [Interval(0.25, 4.0), Interval(1.e-6, 1.e3)]:
        exponent = Interval(0.5, 2.0)
        self.check_nonlinear_interval_image(
          "pow(x,Interval)", input,
          lambda x: pow(x, exponent), lambda x: pow(x, exponent))

        variables = AffineVariables(IntervalVector([input, exponent]))
        result = pow(variables[0], variables[1])
        self.check_affine_inclu(result, pow(input, exponent))

      # Trigonometric functions. Tangent intervals stay inside one
      # continuity branch, while sine and cosine additionally cover
      # several extrema.
      for input in [Interval(-20.0, 20.0), Interval(-1.e-8, 1.e-8), Interval(0.25, 5.75)]:
        self.check_nonlinear_operation("sin", input, input, sin, sin, 8192)
        self.check_nonlinear_operation("cos", input, input, cos, cos, 8192)
      for input in [Interval(-1.4, 1.4), Interval(1.7, 4.5), Interval(-1.e-8, 1.e-8)]:
        self.check_nonlinear_operation("tan", input, input, tan, tan, 8192)

      # Inverse trigonometric functions, including largely invalid input
      # intervals and narrow valid parts near the singular endpoints.
      for input in [Interval(-1.0, 1.0), Interval(-1.e12, -1.0 + 1.e-12),
                    Interval(1.0 - 1.e-12, 1.e12), Interval(-1.0 - 1.e-12, 1.0 + 1.e-12)]:
        domain = input & Interval(-1.0, 1.0)
        self.check_nonlinear_operation("asin", input, domain, asin, asin, 8192)
        self.check_nonlinear_operation("acos", input, domain, acos, acos, 8192)
      for input in [Interval(-1.e12, 1.e12), Interval(-1.e-12, 1.e-12), Interval(-10.0, 3.0)]:
        self.check_nonlinear_operation("atan", input, input, atan, atan, 8192)

      # Hyperbolic functions. Inputs stay below scalar overflow for
      # sinh/cosh.
      for input in [Interval(-20.0, 20.0), Interval(-1.e-8, 1.e-8), Interval(-5.0, 2.0)]:
        self.check_nonlinear_operation("sinh", input, input, sinh, sinh)
        self.check_nonlinear_operation("cosh", input, input, cosh, cosh)
        self.check_nonlinear_operation("tanh", input, input, tanh, tanh)
        self.check_nonlinear_operation("asinh", input, input, asinh, asinh)

      # acosh and atanh with their restricted real domains.
      for input in [Interval(1.0, 10.0), Interval(-1.e12, 1.0 + 1.e-12),
                    Interval(1.0 - 1.e-12, 1.0 + 1.e-8), Interval(-1.e12, 2.0)]:
        domain = input & Interval(1.0, oo)
        self.check_nonlinear_operation("acosh", input, domain, acosh, acosh, 8192)
      for input in [Interval(-0.99, 0.99), Interval(-1.e-12, 1.e-12),
                    Interval(-0.999999, 0.5), Interval(-0.5, 0.999999)]:
        self.check_nonlinear_operation("atanh", input, input, atanh, atanh, 8192)

    Affine.change_mode(Affine.Affine_Mode.AF_Lin_Chebyshev)

  def test_all_nonlinear_functions_accept_correlated_affine_expressions(self):

    modes = [Affine.Affine_Mode.AF_Lin_Chebyshev, Affine.Affine_Mode.AF_Lin_MinRange]

    for mode in modes:
      Affine.change_mode(mode)

      variables = AffineVariables(IntervalVector([[-1.0, 1.0]]))
      t = variables[0]

      general = 3.0*t - 0.5      # [-3.5,2.5]
      positive = 1.5*t + 2.0     # [0.5,3.5]
      nonzero = 0.25*t + 1.0     # [0.75,1.25]
      unit = 0.75*t              # [-0.75,0.75]
      sqrt_arg = 5.0*t - 4.0     # [-9,1]
      acosh_arg = 5.0*t - 1.0    # [-6,4]

      self.check_affine_inclu(inv(nonzero), 1.0/nonzero.itv())
      self.check_affine_inclu(sqr(general), sqr(general.itv()))
      self.check_affine_inclu(abs(general), abs(general.itv()))
      self.check_affine_inclu(sqrt(sqrt_arg), sqrt(sqrt_arg.itv()))
      self.check_affine_inclu(exp(general), exp(general.itv()))
      self.check_affine_inclu(log(positive), log(positive.itv()))
      self.check_affine_inclu(pow(general, 3), pow(general.itv(), 3))
      self.check_affine_inclu(pow(positive, 1.5), pow(positive.itv(), 1.5))
      self.check_affine_inclu(root(general, 3), root(general.itv(), 3))

      self.check_affine_inclu(sin(general), sin(general.itv()))
      self.check_affine_inclu(cos(general), cos(general.itv()))
      self.check_affine_inclu(tan(unit), tan(unit.itv()))
      self.check_affine_inclu(asin(general), asin(general.itv()))
      self.check_affine_inclu(acos(general), acos(general.itv()))
      self.check_affine_inclu(atan(general), atan(general.itv()))

      self.check_affine_inclu(sinh(general), sinh(general.itv()))
      self.check_affine_inclu(cosh(general), cosh(general.itv()))
      self.check_affine_inclu(tanh(general), tanh(general.itv()))
      self.check_affine_inclu(asinh(general), asinh(general.itv()))
      self.check_affine_inclu(acosh(acosh_arg), acosh(acosh_arg.itv()))
      self.check_affine_inclu(atanh(unit), atanh(unit.itv()))

      exponent = Interval(0.5, 2.0)
      self.check_affine_inclu(pow(positive, exponent), pow(positive.itv(), exponent))

      exponent_variable = AffineVariables(IntervalVector([exponent]))
      self.check_affine_inclu(
        pow(positive, exponent_variable[0]), pow(positive.itv(), exponent))

    Affine.change_mode(Affine.Affine_Mode.AF_Lin_Chebyshev)

  def test_nonsmooth_nonlinear_affine_functions_enclose_reference_values(self):

    modes = [Affine.Affine_Mode.AF_Lin_Chebyshev, Affine.Affine_Mode.AF_Lin_MinRange]

    for mode in modes:
      Affine.change_mode(mode)

      for input in [Interval(-4.75, 3.25), Interval(-1.0, 1.0),
                    Interval(1.25, 1.75), Interval(-2.0, -2.0)]:
        # floor()/ceil() return an Interval even for an Affine argument
        # (like their Interval overload); C++ promotes it back to AffineT
        # implicitly where the result is bound to an AffineT-typed
        # variable. Affine(...) makes that promotion explicit here.
        self.check_nonlinear_operation(
          "floor", input, input, lambda x: Affine(floor(x)), floor)
        self.check_nonlinear_operation(
          "ceil", input, input, lambda x: Affine(ceil(x)), ceil)
        self.check_nonlinear_interval_image(
          "integer", input, lambda x: Affine(integer(x)), integer)

      for input in [Interval(-4.0, 3.0), Interval(-3.0, -0.25),
                    Interval(0.25, 5.0), Interval(0.0, 0.0)]:
        self.check_nonlinear_interval_image("sign", input, sign, sign)

    Affine.change_mode(Affine.Affine_Mode.AF_Lin_Chebyshev)

  def test_atan2_encloses_reference_two_argument_arctangent(self):

    cases = [
      (Interval(1.0, 2.0), Interval(1.0, 2.0)),     # first quadrant
      (Interval(1.0, 2.0), Interval(-2.0, -1.0)),   # second quadrant
      (Interval(-2.0, -1.0), Interval(-2.0, -1.0)), # third quadrant
      (Interval(-2.0, -1.0), Interval(1.0, 2.0)),   # fourth quadrant
      (Interval(-0.1, 0.1), Interval(1.0, 2.0))     # straddles y=0, x>0
    ]

    for y_itv, x_itv in cases:
      variables = AffineVariables(IntervalVector([y_itv, x_itv]))
      result = atan2(variables[0], variables[1])
      reference = atan2(y_itv, x_itv)

      self.check_affine_inclu(result, reference)

      for k in range(51):
        y = y_itv.lb() + (y_itv.ub()-y_itv.lb())*k/50.0
        x = x_itv.lb() + (x_itv.ub()-x_itv.lb())*k/50.0
        self.assertTrue(result.contains(math.atan2(y, x)))

  def test_chi_selects_branch_according_to_condition_sign(self):

    bx = AffineVariables(2)
    bx[0] = Interval(10.0, 11.0)  # branch b
    bx[1] = Interval(20.0, 21.0)  # branch c

    negative = Affine(Interval(-2.0, -1.0))
    positive = Affine(Interval(1.0, 2.0))

    # Condition strictly negative: selects b, for every Interval/Affine
    # combination of the operands that the API exposes.
    self.assertTrue(chi(negative, bx[0], bx[1]).itv() == bx[0].itv())
    self.assertTrue(chi(negative.itv(), bx[0], bx[1]).itv() == bx[0].itv())
    self.assertTrue(chi(negative.itv(), bx[0].itv(), bx[1]).itv() == bx[0].itv())
    self.assertTrue(chi(negative.itv(), bx[0], bx[1].itv()).itv() == bx[0].itv())
    self.assertTrue(chi(negative, Interval(30.0, 31.0), bx[1]).itv() == Interval(30.0, 31.0))

    # Condition strictly positive: selects c.
    self.assertTrue(chi(positive, bx[0], bx[1]).itv() == bx[1].itv())
    self.assertTrue(chi(positive.itv(), bx[0], bx[1]).itv() == bx[1].itv())
    self.assertTrue(chi(positive, bx[0], Interval(30.0, 31.0)).itv() == Interval(30.0, 31.0))

    # Condition straddling zero: both branches remain possible, so the
    # result must enclose their union.
    straddling = chi(Interval(-1.0, 1.0), bx[0], bx[1])
    self.assertTrue(straddling.itv().is_superset(bx[0].itv() | bx[1].itv()))

    # An empty condition yields an empty result.
    self.assertTrue(chi(Interval.empty(), bx[0], bx[1]).is_empty())

  # NOTE: distance(x1,x2) is declared and documented on AffineMain<T>
  # (codac2_AffineMain.h) for 3 overloads, and each one forwards to a free
  # function codac2::distance(const Interval&, const Interval&) that does
  # not exist anywhere in the codebase, so instantiating any of them fails
  # to compile in C++ itself. distance() is therefore not exposed to
  # Python either, and is left untested here, matching the C++ test file.

  def test_inflate_widens_enclosure_by_given_radius_on_each_side(self):

    ax = AffineVariables(1)
    ax[0] = Interval(1.0, 2.0)
    before = ax[0].itv()

    # Affine(...) copies: inflate() mutates in place, and ax[0] is read
    # again below, so aliasing it here would corrupt that later read.
    inflated = Affine(ax[0])
    inflated.inflate(0.5)

    self.check_affine_inclu(inflated, before + Interval(-0.5, 0.5))
    self.assertTrue(inflated.noise_count() == ax[0].noise_count())

    # Inflating by zero must not change the enclosure.
    unchanged = Affine(ax[0])
    unchanged.inflate(0.0)
    self.assertTrue(unchanged.itv() == before)

  def test_mig_mag_smag_smig_and_volume_delegate_to_interval_enclosure(self):

    ax = AffineVariables(3)
    ax[0] = Interval(-3.0, 5.0)   # straddles zero
    ax[1] = Interval(2.0, 4.0)    # strictly positive
    ax[2] = Interval(-5.0, -1.0)  # strictly negative

    for i in range(ax.size()):
      self.assertTrue(ax[i].mig() == ax[i].itv().mig())
      self.assertTrue(ax[i].mag() == ax[i].itv().mag())
      self.assertTrue(ax[i].smag() == ax[i].itv().smag())
      self.assertTrue(ax[i].smig() == ax[i].itv().smig())
      self.assertTrue(ax[i].volume() == ax[i].itv().volume())
      self.assertTrue(ax[i].volume() == ax[i].diam())

    # The documented formula: mig = 0 when the interval straddles zero,
    # +lb when strictly positive, -ub when strictly negative.
    self.assertTrue(ax[0].mig() == 0.0)
    self.assertTrue(ax[1].mig() == 2.0)
    self.assertTrue(ax[2].mig() == 1.0)
    self.assertTrue(ax[0].mag() == 5.0)


  def affine_eval_at(self, y, eps):
    # Evaluates a one-variable affine form at a given noise value eps in
    # [-1,1]. This is the linearization itself, as opposed to itv() which is
    # its range over the whole box.
    # A form that fell back to its interval enclosure keeps no noise symbol,
    # and noise(0) asserts in that case.
    if not y.is_active() or y.noise_count() < 1:
      return y.itv()
    return Interval(y.mid()) + Interval(y.noise(0))*Interval(eps) \
         + Interval(-y.err(), y.err())

  def check_pointwise_linearization(self, name, input, domain,
                                    affine_function, interval_function,
                                    sample_count = 50):
    # At every point of the box, the linearization must enclose the function.
    # This is strictly stronger than checking the image, which only constrains
    # the hull and would accept a linearization with a wrong slope.
    variables = AffineVariables(IntervalVector([input]))
    y = affine_function(variables[0])
    reference = interval_function(input & domain)

    if reference.is_empty():
      self.assertTrue(y.is_empty(), name)
      return

    self.assertFalse(y.is_empty(), name)
    self.assertTrue(y.itv().is_superset(reference), name)

    if input.is_unbounded() or y.itv().is_unbounded():
      return

    for k in range(sample_count+1):
      eps = -1.0 + (2.0*k)/sample_count
      x = (Interval(input.mid()) + Interval(input.rad())*Interval(eps)) \
          & input & domain
      if x.is_empty():
        continue
      fx = interval_function(x)
      if fx.is_empty() or fx.is_unbounded():
        continue
      y_at_eps = self.affine_eval_at(y, eps)
      if y_at_eps.is_unbounded():
        continue
      self.assertTrue(y_at_eps.is_superset(fx),
                      "%s input=%s eps=%g x=%s f(x)=%s y=%s"
                      % (name, input, eps, x, fx, y_at_eps))

  def test_every_linearization_encloses_the_function_at_each_point(self):

    modes = [Affine.Affine_Mode.AF_Lin_Chebyshev, Affine.Affine_Mode.AF_Lin_MinRange]

    for mode in modes:
      Affine.change_mode(mode)

      for input in [Interval(0.5,3.0), Interval(1.e-3,1.e3), Interval(-4.0,-0.25)]:
        self.check_pointwise_linearization("inv", input, Interval(-oo,oo), inv, lambda x: 1.0/x)

      for input in [Interval(-4.0,3.0), Interval(0.25,5.0), Interval(-1.e6,1.e6)]:
        self.check_pointwise_linearization("sqr", input, Interval(-oo,oo), sqr, sqr)

      for input in [Interval(0.0,4.0), Interval(1.e-6,1.e6), Interval(-1.0,4.0)]:
        self.check_pointwise_linearization("sqrt", input, Interval(0.0,oo), sqrt, sqrt)

      for input in [Interval(-20.0,20.0), Interval(-1.0,1.0), Interval(-700.0,-600.0)]:
        self.check_pointwise_linearization("exp", input, Interval(-oo,oo), exp, exp)

      for input in [Interval(1.e-6,1.e6), Interval(0.25,4.0), Interval(-1.0,4.0)]:
        self.check_pointwise_linearization("log", input, Interval(0.0,oo), log, log)

      for input in [Interval(-8.0,27.0), Interval(-3.0,2.0), Interval(0.5,4.0)]:
        self.check_pointwise_linearization("pow(x,3)", input, Interval(-oo,oo), lambda x: pow(x,3), lambda x: pow(x,3))

      for input in [Interval(-2.0,3.0), Interval(1.0,5.0), Interval(-6.0,-2.0)]:
        self.check_pointwise_linearization("pow(x,4)", input, Interval(-oo,oo), lambda x: pow(x,4), lambda x: pow(x,4))

      for input in [Interval(0.5,4.0), Interval(-4.0,-0.5)]:
        self.check_pointwise_linearization("pow(x,-2)", input, Interval(-oo,oo), lambda x: pow(x,-2), lambda x: pow(x,-2))

      for input in [Interval(0.25,9.0), Interval(1.e-3,1.e3)]:
        self.check_pointwise_linearization("pow(x,1.5)", input, Interval(0.0,oo), lambda x: pow(x,1.5), lambda x: pow(x,1.5))

      for input in [Interval(0.0,9.0), Interval(-1.0,9.0)]:
        self.check_pointwise_linearization("root(x,2)", input, Interval(0.0,oo), lambda x: root(x,2), lambda x: root(x,2))

      for input in [Interval(-8.0,27.0), Interval(-3.0,2.0), Interval(1.0,8.0)]:
        self.check_pointwise_linearization("root(x,3)", input, Interval(-oo,oo), lambda x: root(x,3), lambda x: root(x,3))

      for input in [Interval(-1.0,1.0), Interval(2.0,5.0), Interval(-12.0,-8.0)]:
        self.check_pointwise_linearization("cos", input, Interval(-oo,oo), cos, cos)

      for input in [Interval(-1.0,1.0), Interval(2.0,5.0), Interval(-12.0,-8.0)]:
        self.check_pointwise_linearization("sin", input, Interval(-oo,oo), sin, sin)

      for input in [Interval(-1.4,1.4), Interval(0.1,0.9), Interval(-0.5,1.2)]:
        self.check_pointwise_linearization("tan", input, Interval(-oo,oo), tan, tan)

      for input in [Interval(-0.9,0.9), Interval(-1.0,1.0), Interval(-2.0,0.5)]:
        self.check_pointwise_linearization("acos", input, Interval(-1.0,1.0), acos, acos)

      for input in [Interval(-0.9,0.9), Interval(-1.0,1.0), Interval(-0.5,2.0)]:
        self.check_pointwise_linearization("asin", input, Interval(-1.0,1.0), asin, asin)

      for input in [Interval(-20.0,20.0), Interval(-1.0,1.0), Interval(3.0,9.0)]:
        self.check_pointwise_linearization("atan", input, Interval(-oo,oo), atan, atan)

      for input in [Interval(-3.0,5.0), Interval(1.0,4.0), Interval(-6.0,-2.0)]:
        self.check_pointwise_linearization("cosh", input, Interval(-oo,oo), cosh, cosh)

      for input in [Interval(-3.0,5.0), Interval(1.0,4.0)]:
        self.check_pointwise_linearization("sinh", input, Interval(-oo,oo), sinh, sinh)

      for input in [Interval(-3.0,5.0), Interval(-8.0,8.0)]:
        self.check_pointwise_linearization("tanh", input, Interval(-oo,oo), tanh, tanh)

      for input in [Interval(1.0,9.0), Interval(2.0,20.0), Interval(-1.0,3.0)]:
        self.check_pointwise_linearization("acosh", input, Interval(1.0,oo), acosh, acosh)

      for input in [Interval(-8.0,8.0), Interval(1.0,20.0)]:
        self.check_pointwise_linearization("asinh", input, Interval(-oo,oo), asinh, asinh)

      for input in [Interval(-0.9,0.9), Interval(-0.5,0.999), Interval(-2.0,0.5)]:
        self.check_pointwise_linearization("atanh", input, Interval(-1.0,1.0), atanh, atanh)

      for input in [Interval(-4.0,3.0), Interval(1.0,5.0), Interval(-5.0,-1.0)]:
        self.check_pointwise_linearization("abs", input, Interval(-oo,oo), abs, abs)

    Affine.change_mode(Affine.Affine_Mode.AF_Lin_Chebyshev)


if __name__ ==  '__main__':
  unittest.main()
