#!/usr/bin/env python

#  Codac tests
#
#  Rounding of the Interval operations, which Codac delegates to GAOL. The
#  bounds are checked exactly, with rational arithmetic: they have to enclose
#  the exact result of each operation, as tightly as possible for the
#  arithmetic operations, and they have to go on doing so after any operation,
#  GAOL computing with the rounding direction set upward once for all.
#
# ----------------------------------------------------------------------------
#  \date       2026
#  \author     Jordan Ninin
#  \copyright  Copyright 2026 Codac Team
#  \license    GNU Lesser General Public License (LGPL)

import unittest
from codac import *
from fractions import Fraction
import random
import struct

# Note: "from codac import *" replaces the built-in abs, min, max and pow with
# the interval functions of Codac, which is why they are not used on floats
# below. Only exact operations are done on floats: negation and comparison.

NB_RANDOM_VALUES = 1000

def double_from_bits(bits):
  return struct.unpack('<d', struct.pack('<Q', bits))[0]

def bits_from_double(x):
  return struct.unpack('<Q', struct.pack('<d', x))[0]

def random_double(rng, emin, emax):
  # A random sign and mantissa, and an exponent between emin and emax
  return double_from_bits((rng.getrandbits(1) << 63) | ((1023 + rng.randint(emin, emax)) << 52) | rng.getrandbits(52))

def next_double(x):
  # The double immediately above x, a finite double
  if x == 0.0:
    return double_from_bits(1)
  bits = bits_from_double(x)
  return double_from_bits(bits + 1 if x > 0.0 else bits - 1)

def previous_double(x):
  return -next_double(-x)

def ordered(a, b):
  return (a, b) if a <= b else (b, a)

# The doubles immediately below and above the exact value of each function at
# x, computed with 400 bits of precision (mpmath): the enclosure computed by
# Codac has to contain both.
ELEMENTARY_VALUES = [
  ("exp", lambda x: exp(x),  "0x1.0000000000000p+0",  "0x1.5bf0a8b145769p+1",  "0x1.5bf0a8b14576ap+1"),
  ("exp", lambda x: exp(x),  "0x1.0000000000000p-1",  "0x1.a61298e1e069bp+0",  "0x1.a61298e1e069cp+0"),
  ("exp", lambda x: exp(x), "-0x1.0000000000000p+0",  "0x1.78b56362cef37p-2",  "0x1.78b56362cef38p-2"),
  ("exp", lambda x: exp(x),  "0x1.4000000000000p+3",  "0x1.5829dcf95055fp+14", "0x1.5829dcf950560p+14"),
  ("exp", lambda x: exp(x),  "0x1.0624dd2f1a9fcp-10", "0x1.0041919b7ee33p+0",  "0x1.0041919b7ee34p+0"),
  ("exp", lambda x: exp(x), "-0x1.4000000000000p+4",  "0x1.1b48655f37266p-29", "0x1.1b48655f37267p-29"),
  ("log", lambda x: log(x),  "0x1.0000000000000p+1",  "0x1.62e42fefa39efp-1",  "0x1.62e42fefa39f0p-1"),
  ("log", lambda x: log(x),  "0x1.4000000000000p+3",  "0x1.26bb1bbb55515p+1",  "0x1.26bb1bbb55516p+1"),
  ("log", lambda x: log(x),  "0x1.999999999999ap-4", "-0x1.26bb1bbb55516p+1", "-0x1.26bb1bbb55515p+1"),
  ("log", lambda x: log(x),  "0x1.8000000000000p+0",  "0x1.9f323ecbf984bp-2",  "0x1.9f323ecbf984cp-2"),
  ("log", lambda x: log(x),  "0x1.4f8b588e368f1p-17", "-0x1.7069e2aa2aa5bp+3", "-0x1.7069e2aa2aa5ap+3"),
  ("log", lambda x: log(x),  "0x1.2a05f20000000p+33", "0x1.7069e2aa2aa5ap+4",  "0x1.7069e2aa2aa5bp+4"),
  ("sin", lambda x: sin(x),  "0x1.0000000000000p+0",  "0x1.aed548f090ceep-1",  "0x1.aed548f090cefp-1"),
  ("sin", lambda x: sin(x),  "0x1.0000000000000p-1",  "0x1.eaee8744b05efp-2",  "0x1.eaee8744b05f0p-2"),
  ("sin", lambda x: sin(x),  "0x1.8000000000000p+1",  "0x1.210386db6d55bp-3",  "0x1.210386db6d55cp-3"),
  ("sin", lambda x: sin(x), "-0x1.0000000000000p+1", "-0x1.d18f6ead1b446p-1", "-0x1.d18f6ead1b445p-1"),
  ("sin", lambda x: sin(x),  "0x1.9000000000000p+6", "-0x1.03425b78c4db9p-1", "-0x1.03425b78c4db8p-1"),
  ("sin", lambda x: sin(x),  "0x1.0624dd2f1a9fcp-10", "0x1.0624da5218a62p-10", "0x1.0624da5218a63p-10"),
  ("cos", lambda x: cos(x),  "0x1.0000000000000p+0",  "0x1.14a280fb5068bp-1",  "0x1.14a280fb5068cp-1"),
  ("cos", lambda x: cos(x),  "0x1.0000000000000p-1",  "0x1.c1528065b7d4fp-1",  "0x1.c1528065b7d50p-1"),
  ("cos", lambda x: cos(x),  "0x1.8000000000000p+1", "-0x1.fae04be85e5d3p-1", "-0x1.fae04be85e5d2p-1"),
  ("cos", lambda x: cos(x), "-0x1.0000000000000p+1", "-0x1.aa22657537205p-2", "-0x1.aa22657537204p-2"),
  ("cos", lambda x: cos(x),  "0x1.9000000000000p+6",  "0x1.b981dbf665fdfp-1",  "0x1.b981dbf665fe0p-1"),
  ("cos", lambda x: cos(x),  "0x1.8000000000000p+0",  "0x1.21bd54fc5f9a7p-4",  "0x1.21bd54fc5f9a8p-4"),
  ("tan", lambda x: tan(x),  "0x1.0000000000000p+0",  "0x1.8eb245cbee3a5p+0",  "0x1.8eb245cbee3a6p+0"),
  ("tan", lambda x: tan(x),  "0x1.0000000000000p-1",  "0x1.17b4f5bf3474ap-1",  "0x1.17b4f5bf3474bp-1"),
  ("tan", lambda x: tan(x), "-0x1.3333333333333p+0", "-0x1.493c43acb164dp+1", "-0x1.493c43acb164cp+1"),
  ("tan", lambda x: tan(x),  "0x1.8000000000000p+1", "-0x1.23ef71254b870p-3", "-0x1.23ef71254b86fp-3"),
  ("atan", lambda x: atan(x),  "0x1.0000000000000p+0",  "0x1.921fb54442d18p-1",  "0x1.921fb54442d19p-1"),
  ("atan", lambda x: atan(x),  "0x1.0000000000000p-1",  "0x1.dac670561bb4fp-2",  "0x1.dac670561bb50p-2"),
  ("atan", lambda x: atan(x), "-0x1.0000000000000p+1", "-0x1.1b6e192ebbe45p+0", "-0x1.1b6e192ebbe44p+0"),
  ("atan", lambda x: atan(x),  "0x1.f400000000000p+9",  "0x1.91de2c0e658bcp+0",  "0x1.91de2c0e658bdp+0"),
  ("sqrt", lambda x: sqrt(x),  "0x1.0000000000000p+1",  "0x1.6a09e667f3bccp+0",  "0x1.6a09e667f3bcdp+0"),
  ("sqrt", lambda x: sqrt(x),  "0x1.8000000000000p+1",  "0x1.bb67ae8584caap+0",  "0x1.bb67ae8584cabp+0"),
  ("sqrt", lambda x: sqrt(x),  "0x1.999999999999ap-4",  "0x1.43d136248490fp-2",  "0x1.43d1362484910p-2"),
  ("sqrt", lambda x: sqrt(x),  "0x1.ad7f29abcaf48p-24", "0x1.4b96be9c2da2bp-12", "0x1.4b96be9c2da2cp-12"),
]

class TestInterval_rounding(unittest.TestCase):

  def assert_tightest(self, x, exact_values, description):
    # [x] is the tightest interval of doubles enclosing the exact values
    exact_values = sorted(exact_values)
    lo, hi = exact_values[0], exact_values[-1]
    lb, ub = x.lb(), x.ub()
    self.assertTrue(Fraction(lb) <= lo and Fraction(next_double(lb)) > lo,
      "lower bound of " + description + ": " + lb.hex())
    self.assertTrue(Fraction(ub) >= hi and Fraction(previous_double(ub)) < hi,
      "upper bound of " + description + ": " + ub.hex())

  def test_arithmetic_operations(self):

    rng = random.Random(2026)

    for i in range(NB_RANDOM_VALUES):

      a, b = random_double(rng, -30, 30), random_double(rng, -30, 30)
      A, B = Fraction(a), Fraction(b)
      values = "a=" + a.hex() + " b=" + b.hex()

      self.assert_tightest(Interval(a) + Interval(b), [A + B], "a+b, " + values)
      self.assert_tightest(Interval(a) - Interval(b), [A - B], "a-b, " + values)
      self.assert_tightest(Interval(a) * Interval(b), [A * B], "a*b, " + values)
      self.assert_tightest(Interval(a) * b, [A * B], "[a]*b, " + values)
      self.assert_tightest(Interval(a) / Interval(b), [A / B], "a/b, " + values)
      self.assert_tightest(sqr(Interval(a)), [A * A], "sqr(a), " + values)

      # GAOL encloses square roots, not always as tightly as possible:
      # lb^2 <= |a| <= ub^2 only
      m = -a if a < 0.0 else a
      r = sqrt(Interval(m))
      self.assertTrue(r.lb() >= 0.0 and Fraction(r.lb())**2 <= Fraction(m) <= Fraction(r.ub())**2,
        "sqrt(|a|), |a|=" + m.hex() + ": [" + r.lb().hex() + ", " + r.ub().hex() + "]")

      # Products of intervals, whose bounds have all the signs GAOL distinguishes:
      # the tightest enclosure of the four products of bounds
      c, e = random_double(rng, -30, 30), random_double(rng, -30, 30)
      x, y = Interval(*ordered(a, c)), Interval(*ordered(b, e))
      corners = [Fraction(u) * Fraction(v) for u in (x.lb(), x.ub()) for v in (y.lb(), y.ub())]
      self.assert_tightest(x * y, corners,
        "[x]*[y], x=[" + x.lb().hex() + ", " + x.ub().hex() + "] y=[" + y.lb().hex() + ", " + y.ub().hex() + "]")

  def test_elementary_functions(self):

    for name, f, x, below, above in ELEMENTARY_VALUES:
      y = f(Interval(float.fromhex(x)))
      description = name + "(" + x + ") = [" + y.lb().hex() + ", " + y.ub().hex() + "]"
      self.assertTrue(y.lb() <= float.fromhex(below), description)
      self.assertTrue(y.ub() >= float.fromhex(above), description)

    # pi, whose neighbours scaled by a power of two are those of pi/2 and 2pi
    pi_below, pi_above = float.fromhex("0x1.921fb54442d18p+1"), float.fromhex("0x1.921fb54442d19p+1")
    self.assertTrue(Interval.pi().lb() <= pi_below and Interval.pi().ub() >= pi_above)
    self.assertTrue(Interval.half_pi().lb() <= float.fromhex("0x1.921fb54442d18p+0")
      and Interval.half_pi().ub() >= float.fromhex("0x1.921fb54442d19p+0"))
    self.assertTrue(Interval.two_pi().lb() <= float.fromhex("0x1.921fb54442d18p+2")
      and Interval.two_pi().ub() >= float.fromhex("0x1.921fb54442d19p+2"))

  def test_after_any_operation(self):

    # Codac and GAOL compute with the rounding direction set upward, once for
    # all: every operation has to leave it so, or the operations after it return
    # wrong bounds. x is within [0,1] and y above 1, which is the domain of every
    # function below.
    x, y = Interval(0.1, 0.3), Interval(1.5, 2.5)

    operations = [
      ("Interval.mid()", lambda: x.mid()),
      ("Interval.rad()", lambda: x.rad()),
      ("Interval.diam()", lambda: x.diam()),
      ("Interval.mag()", lambda: x.mag()),
      ("Interval.mig()", lambda: x.mig()),
      ("Interval.rand()", lambda: x.rand()),
      ("Interval.inflate()", lambda: Interval(x).inflate(0.1)),
      ("str(Interval)", lambda: str(x)),
      ("sqr", lambda: sqr(x)),
      ("sqrt", lambda: sqrt(x)),
      ("pow(x,int)", lambda: pow(y, 3)),
      ("pow(x,float)", lambda: pow(y, 2.5)),
      ("pow(x,Interval)", lambda: pow(y, x)),
      ("root", lambda: root(y, 3)),
      ("exp", lambda: exp(x)),
      ("log", lambda: log(x)),
      ("cos", lambda: cos(x)),
      ("sin", lambda: sin(x)),
      ("tan", lambda: tan(x)),
      ("acos", lambda: acos(x)),
      ("asin", lambda: asin(x)),
      ("atan", lambda: atan(x)),
      ("atan2", lambda: atan2(y, x)),
      ("cosh", lambda: cosh(x)),
      ("sinh", lambda: sinh(x)),
      ("tanh", lambda: tanh(x)),
      ("acosh", lambda: acosh(y)),
      ("asinh", lambda: asinh(x)),
      ("atanh", lambda: atanh(x)),
      ("abs", lambda: abs(x)),
      ("min", lambda: min(x, y)),
      ("max", lambda: max(x, y)),
      ("sign", lambda: sign(x)),
      ("integer", lambda: integer(y)),
      ("floor", lambda: floor(y)),
      ("ceil", lambda: ceil(y)),
      ("Interval.pi()", lambda: Interval.pi()),
    ]

    rng = random.Random(2026)

    for name, operation in operations:
      operation()
      # A product and a sum whose exact results are not doubles, as the random
      # mantissas make almost certain: their bounds are the tightest ones only
      # if the rounding direction is still upward
      a, b = random_double(rng, -30, 30), random_double(rng, -30, 30)
      values = " after " + name + ", a=" + a.hex() + " b=" + b.hex()
      self.assert_tightest(Interval(a) * Interval(b), [Fraction(a) * Fraction(b)], "a*b" + values)
      self.assert_tightest(Interval(a) + Interval(b), [Fraction(a) + Fraction(b)], "a+b" + values)

if __name__ ==  '__main__':
  unittest.main()
