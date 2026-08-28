#!/usr/bin/env python

# Codac tests
#
# Unit tests for the "value" containers AffineVector and AffineRow
# (the Python bindings of AffineMainVector<T> and AffineMainRow<T>):
# construction, indexing, segment views, transpose, dot product, sum/
# reverse, cwiseProduct/cwiseQuotient, compound assignments (scalar,
# IntervalVector, AffineVector, AffineVariables), and the arithmetic of
# AffineVector (linear combinations, empty propagation, component-wise
# nonlinear expressions).
#
# codac2_py_AffineVector_templ.h / codac2_py_AffineMatrixBase.h expose a
# hand-picked, reduced subset of Eigen's API (documented as such in the
# latter's header comment) rather than the full generic MatrixBase
# binding. In particular there is no dot()/sum()/reverse()/cwiseProduct()/
# cwiseQuotient()/coeffRef()/head()/tail()/array()/swap()/eval(), no
# lvalue segment() (only the value-returning subvector()), and no
# fixed-size Eigen container (Eigen::Matrix<Affine,N,M>) at all. Where a
# C++ test relies on one of these, this file either reconstructs the same
# observable behavior from the bound primitives (indexing, +, *, /,
# AffineVector(list)) or, for fixed-size containers, is not translated
# (there is nothing in the Python API to translate it to).
#
# Python also has no assignment operator distinct from construction:
# wherever the C++ source reassigns an *existing* AffineVariables from an
# IntervalVector (relying on AffineVarMainVector::operator=(const
# IntervalVector&), which also resizes), this file constructs a fresh
# AffineVariables from that IntervalVector directly. And wherever it
# copies an existing AffineVector/AffineVariables into a new C++ variable
# (a real, independent copy), this file uses AffineVector(source)
# explicitly: a bare Python "y = x" would alias the same object rather
# than copy it, which would silently defeat every test that checks that a
# copy is independent from its source. The same applies to a single
# component read through AffineVariables/AffineVector/AffineRow
# indexing, which -- like the C++ operator[] it mirrors -- returns a live
# reference, not a copy: naming it (e.g. "x = variables[0]") and then
# mutating that name in place would reach back into the container.
#
# codac2::Approx<AffineMain<T>>/Approx<AffineMainVector<T>> are exposed
# to Python as Approx_Affine/Approx_AffineVector (mirroring
# Approx_double/Approx_Interval/Approx_Vector/...), used here directly
# rather than through the generic Approx() dispatcher: their constructor
# takes an Interval/IntervalVector, exactly like Approx_Interval/
# Approx_IntervalVector's, so Approx() cannot tell them apart from the
# argument's type alone.
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

ERROR = sys.float_info.epsilon*100


def affine_dot(x, y):
  """x.dot(y): dot() is not exposed to Python, only reconstructed here
  from indexing and the bound +/* operators."""
  result = x[0] * y[0]
  for i in range(1, len(x)):
    result = result + x[i] * y[i]
  return result

def affine_sum(x):
  result = x[0]
  for i in range(1, len(x)):
    result = result + x[i]
  return result

def affine_reverse(x):
  return AffineVector([x[len(x)-1-i] for i in range(len(x))])

def affine_cwise_product(x, y):
  return AffineVector([x[i]*y[i] for i in range(len(x))])

def affine_cwise_quotient(x, y):
  return AffineVector([x[i]/y[i] for i in range(len(x))])

def make_point_vector_5():
  variables = AffineVariables(IntervalVector([[1],[2],[3],[4],[5]]))
  return AffineVector(variables)


class TestAffineVector(unittest.TestCase):

  def check_point(self, actual, expected):
    self.assertTrue(actual == Approx_Affine(Interval(expected), ERROR))

  def check_interval(self, actual, expected):
    self.assertTrue(actual == Approx_Affine(expected, ERROR))

  def check_interval_enclosure(self, actual, expected):
    self.assertTrue(actual.itv().is_superset(expected))

  def test_construction_copy_and_access(self):

    # Construction from an AffineVariables and copy
    x = AffineVariables(IntervalVector([[1.,2.],[3.,4.],[5.,6.]]))
    y = AffineVector(x)

    self.assertTrue(y.size() == 3)
    self.assertTrue(y.itv() == IntervalVector([[1.,2.],[3.,4.],[5.,6.]]))
    self.assertTrue(y[0].itv() == Interval(1.,2.))
    self.assertTrue(y[1].itv() == Interval(3.,4.))
    self.assertTrue(y[2].itv() == Interval(5.,6.))

    # Copy of a copy
    x = AffineVariables(IntervalVector([[-2.,1.],[4.,7.]]))
    y = AffineVector(x)
    z = AffineVector(y)

    self.assertTrue(z.size() == 2)
    self.assertTrue(z.itv() == y.itv())
    self.assertTrue(z.itv() == IntervalVector([[-2.,1.],[4.,7.]]))

    # Component assignment on the destination leaves the source vector
    # unchanged
    x = AffineVariables(IntervalVector([[1.,2.],[3.,4.]]))
    y = AffineVector(x)
    y[0] = Interval(-1.,1.)

    self.assertTrue(y.itv() == IntervalVector([[-1.,1.],[3.,4.]]))
    self.assertTrue(x.itv() == IntervalVector([[1.,2.],[3.,4.]]))

    # Direct component assignment via operator[]
    x = AffineVariables(3)
    x[2] = Interval(-1., 1.)
    self.assertTrue(x[2].itv() == Interval(-1., 1.))

    # conservativeResize grows and then shrinks the vector
    x = AffineVariables(IntervalVector([[-1.,2.],[3.,4.]]))
    x.conservativeResize(3)
    self.assertTrue(x.size() == 3)
    self.assertTrue(x[0].itv() == Interval(-1.,2.))
    self.assertTrue(x[1].itv() == Interval(3.,4.))
    self.assertTrue(x[2].itv() == Interval())
    x.conservativeResize(1)

    self.assertTrue(x.size() == 1)
    self.assertTrue(x[0].itv() == Interval(-1.,2.))

    # conservativeResize to the same size is a no-op
    x = AffineVariables(IntervalVector([[1.,2.],[3.,4.]]))
    x.conservativeResize(2) # même taille
    self.assertTrue(x[0].itv() == Interval(1.,2.))
    self.assertTrue(x[1].itv() == Interval(3.,4.))

    # Construction through AffineVariables from an IntervalVector
    box = IntervalVector([[-5.,-2.],[0.],[7.,9.],[-oo,oo]])
    x = AffineVariables(box)
    y = AffineVector(x)

    self.assertTrue(y.size() == 4)
    self.assertTrue(y.rows() == 4)
    self.assertTrue(y.cols() == 1)
    self.assertTrue(y.itv() == box)

    # Construction through AffineVariables from a real Vector
    values = Vector([-2.,0.,3.5,8.])
    x = AffineVariables(values)
    y = AffineVector(x)

    self.assertTrue(y.size() == 4)
    self.assertTrue(y.itv() == IntervalVector([[-2.],[0.],[3.5],[8.]]))

    # The default affine variables are unbounded
    x = AffineVariables(3)
    y = AffineVector(x)

    self.assertTrue(y.size() == 3)
    self.assertTrue(y.itv() == IntervalVector(3))
    self.assertFalse(y.is_empty())

    # Parentheses, brackets and Eigen's coeff() address the same elements
    # (coeff() is not exposed: [] already covers the same access)
    x = AffineVariables(IntervalVector([[-1.,1.],[2.,4.],[6.]]))
    y = AffineVector(x)

    self.assertTrue(y[2].itv() == Interval(6.))

    # Component assignment does not change the vector dimension
    x = AffineVariables(IntervalVector([[1.,2.],[3.,4.],[5.,6.],[7.,8.]]))
    y = AffineVector(x)

    y[1] = Interval(-10.,-8.)
    y[3] = Interval(0.)

    self.assertTrue(y.size() == 4)
    self.assertTrue(y.itv() == IntervalVector([[1.,2.],[-10.,-8.],[5.,6.],[0.]]))

    # Unbounded components are preserved by itv()
    x = AffineVariables(IntervalVector([[-oo,0.],[1.,oo],[-oo,oo]]))
    y = AffineVector(x)

    self.assertTrue(y.itv() == IntervalVector([[-oo,0.],[1.,oo],[-oo,oo]]))
    self.assertTrue(y.itv().is_unbounded())

    # AffineVariables must be built directly from the IntervalVector
    variables = AffineVariables(IntervalVector([[1.,2.],[3.,4.],[5.,6.]]))

    self.assertTrue(variables.size() == 3)
    self.assertTrue(AffineVector(variables) == Approx_AffineVector(
      IntervalVector([[1.,2.],[3.,4.],[5.,6.]]), ERROR))

    # Explicit conversion from affine variables to ordinary affine forms
    variables = AffineVariables(IntervalVector([[-2.,1.],[4.,7.]]))
    x = AffineVector(variables)

    self.assertTrue(x.size() == 2)
    self.assertTrue(x.itv() == IntervalVector([[-2.,1.],[4.,7.]]))

    # Copy construction of an AffineVector
    variables = AffineVariables(IntervalVector([[1.,2.],[3.,4.],[5.,6.]]))
    x = AffineVector(variables)
    y = AffineVector(x)

    self.assertTrue(y.size() == 3)
    self.assertTrue(y.itv() == x.itv())

    # Modification of components and independence after copy
    variables = AffineVariables(IntervalVector([[1.,2.],[3.,4.],[5.,6.]]))
    x = AffineVector(variables)
    y = AffineVector(x)

    y[0] = Interval(-1.,1.)
    y[2] = Interval(0.)

    self.assertTrue(x.itv() == IntervalVector([[1.,2.],[3.,4.],[5.,6.]]))
    self.assertTrue(y.itv() == IntervalVector([[-1.,1.],[3.,4.],[0.]]))

  def test_empty_propagation(self):

    # set_empty() on a copy leaves the source unaffected
    x = AffineVariables(IntervalVector([[1.,2.],[3.,4.],[5.,6.]]))
    y = AffineVector(x)

    self.assertFalse(y.is_empty())
    y.set_empty()
    self.assertTrue(y.is_empty())
    self.assertTrue(y.size() == 3)

    # A single empty component makes the whole vector empty
    x = AffineVariables(IntervalVector([[1.,2.],Interval.empty(),[5.,6.]]))
    y = AffineVector(x)

    self.assertTrue(y.is_empty())
    self.assertTrue(y.size() == 3)
    self.assertTrue(y.itv().is_empty())

    # An empty component survives being copied twice
    x = AffineVariables(IntervalVector([[1.,2.],Interval.empty(),[3.,4.]]))
    y = AffineVector(x)
    z = AffineVector(y)

    self.assertTrue(y.size() == 3)
    self.assertTrue(y.is_empty())
    self.assertTrue(z.size() == 3)
    self.assertTrue(z.is_empty())

    # Setting one component to empty makes the whole vector empty
    x = AffineVariables(IntervalVector([[1.,2.],[3.,4.],[5.,6.]]))
    y = AffineVector(x)
    y[1] = Interval.empty()

    self.assertTrue(y.is_empty())
    self.assertTrue(y.itv().is_empty())
    self.assertTrue(y.size() == 3)

    # Empty vectors keep their dimension through a copy
    x = AffineVariables(IntervalVector([[1.,2.],[3.,4.]]))
    y = AffineVector(x)
    y.set_empty()
    z = AffineVector(y)

    self.assertTrue(z.size() == 2)
    self.assertTrue(z.is_empty())
    self.assertTrue(z.itv().is_empty())

    # Empty values propagate through vector addition and subtraction
    x = AffineVariables(IntervalVector([[1.,2.],[3.,4.]]))
    e = AffineVariables(IntervalVector([Interval.empty(),Interval.empty()]))

    total = AffineVector(x) + AffineVector(e)
    difference = AffineVector(e) - AffineVector(x)

    self.assertTrue(total.is_empty())
    self.assertTrue(difference.is_empty())
    self.assertTrue(total.size() == 2)
    self.assertTrue(difference.size() == 2)

    # Empty vector propagation through unary minus and scalar operations
    e = AffineVariables(IntervalVector([Interval.empty(),Interval.empty(),Interval.empty()]))
    ev = AffineVector(e)

    opposite = -ev
    product_left = 2. * ev
    product_zero = 0. * ev
    quotient = ev / 2.

    self.assertTrue(opposite.is_empty())
    self.assertTrue(product_left.is_empty())
    self.assertTrue(product_zero.is_empty())
    self.assertTrue(quotient.is_empty())
    self.assertTrue(opposite.size() == 3)
    self.assertTrue(product_left.size() == 3)
    self.assertTrue(product_zero.size() == 3)
    self.assertTrue(quotient.size() == 3)

    # Empty vector propagation for both operand orders
    x = AffineVector(AffineVariables(IntervalVector([[1.,2.],[-3.,4.],[5.,6.]])))
    e = AffineVector(AffineVariables(IntervalVector([Interval.empty(),Interval.empty(),Interval.empty()])))

    sum_left = e + x
    sum_right = x + e
    difference_left = e - x
    difference_right = x - e

    self.assertTrue(sum_left.is_empty())
    self.assertTrue(sum_right.is_empty())
    self.assertTrue(difference_left.is_empty())
    self.assertTrue(difference_right.is_empty())
    self.assertTrue(sum_left.size() == 3)
    self.assertTrue(sum_right.size() == 3)
    self.assertTrue(difference_left.size() == 3)
    self.assertTrue(difference_right.size() == 3)

    # Empty propagation through a nested expression
    x = AffineVector(AffineVariables(IntervalVector([[1.,2.],[3.,4.]])))
    y = AffineVector(AffineVariables(IntervalVector([[-2.,1.],[5.,7.]])))
    e = AffineVector(AffineVariables(IntervalVector([Interval.empty(),Interval.empty()])))

    result = 2. * (x + y) - (3. * e - x) / 2.

    self.assertTrue(result.is_empty())
    self.assertTrue(result.size() == 2)
    self.assertTrue(result.itv().is_empty())

    # Empty propagation through intermediate AffineVector results
    x = AffineVector(AffineVariables(IntervalVector([[1.,3.],[-2.,4.]])))
    e = AffineVector(AffineVariables(IntervalVector([Interval.empty(),Interval.empty()])))

    first = x + e
    second = 3. * first - x
    third = (second + x) / 4.

    self.assertTrue(first.is_empty())
    self.assertTrue(second.is_empty())
    self.assertTrue(third.is_empty())
    self.assertTrue(first.size() == 2)
    self.assertTrue(second.size() == 2)
    self.assertTrue(third.size() == 2)

    # A single empty component propagates through every resulting vector
    x = AffineVector(AffineVariables(IntervalVector([[1.,2.],Interval.empty(),[5.,6.]])))
    y = AffineVector(AffineVariables(IntervalVector([[-1.,1.],[3.,4.],[2.,3.]])))

    total = x + y
    difference = y - x
    expression = 2. * x - 3. * y

    self.assertTrue(total.is_empty())
    self.assertTrue(difference.is_empty())
    self.assertTrue(expression.is_empty())
    self.assertTrue(total.size() == 3)
    self.assertTrue(difference.size() == 3)
    self.assertTrue(expression.size() == 3)
    self.assertTrue(total.itv().is_empty())
    self.assertTrue(difference.itv().is_empty())
    self.assertTrue(expression.itv().is_empty())

    # Assigning an empty expression to a non-empty vector propagates emptiness
    x = AffineVector(AffineVariables(IntervalVector([[1.,2.],[3.,4.]])))
    e = AffineVector(AffineVariables(IntervalVector([Interval.empty(),Interval.empty()])))

    result = AffineVector(x)
    self.assertFalse(result.is_empty())

    result = x + e

    self.assertTrue(result.is_empty())
    self.assertTrue(result.size() == 2)
    self.assertTrue(result.itv().is_empty())

    # Empty vector and dimension preservation
    variables = AffineVariables(IntervalVector([[1.,2.],[3.,4.],[5.,6.]]))
    x = AffineVector(variables)

    self.assertFalse(x.is_empty())
    x.set_empty()
    self.assertTrue(x.is_empty())
    self.assertTrue(x.size() == 3)
    self.assertTrue(x.itv().is_empty())

    # Propagation of one empty component during conversion
    variables = AffineVariables(IntervalVector([[1.,2.],Interval.empty(),[5.,6.]]))
    x = AffineVector(variables)

    self.assertTrue(x.is_empty())
    self.assertTrue(x.size() == 3)
    self.assertTrue(x.itv().is_empty())

    # Propagation of emptiness through addition, subtraction and unary minus
    x = AffineVector(AffineVariables(IntervalVector([[1.,2.],[3.,4.]])))
    empty = AffineVector(AffineVariables(IntervalVector([Interval.empty(),Interval.empty()])))

    sum_left = empty + x
    sum_right = x + empty
    difference_left = empty - x
    difference_right = x - empty
    opposite = -empty

    self.assertTrue(sum_left.is_empty())
    self.assertTrue(sum_right.is_empty())
    self.assertTrue(difference_left.is_empty())
    self.assertTrue(difference_right.is_empty())
    self.assertTrue(opposite.is_empty())
    self.assertTrue(sum_left.size() == 2)
    self.assertTrue(difference_right.size() == 2)

    # Propagation of emptiness through scalar and nested operations
    x = AffineVector(AffineVariables(IntervalVector([[1.,2.],[3.,4.]])))
    empty = AffineVector(AffineVariables(IntervalVector([Interval.empty(),Interval.empty()])))

    product = 2. * empty
    zero_product = 0. * empty
    nested = 2. * (x + empty) - 0.5 * x

    self.assertTrue(product.is_empty())
    self.assertTrue(zero_product.is_empty())
    self.assertTrue(nested.is_empty())
    self.assertTrue(product.size() == 2)
    self.assertTrue(nested.itv().is_empty())

  def test_linear_arithmetic(self):

    # Addition of two affine-variable vectors
    x = AffineVariables(IntervalVector([[1.,2.],[-3.,-1.]]))
    y = AffineVariables(IntervalVector([[4.,5.],[2.,6.]]))
    z = AffineVector(x) + AffineVector(y)

    self.assertTrue(z.size() == 2)
    self.assertTrue(z.itv() == IntervalVector([[5.,7.],[-1.,5.]]))

    # Subtraction of an affine-variable vector from an equal copy is
    # exactly zero
    x = AffineVariables(IntervalVector([[1.,3.],[-3.,1.]]))
    z = AffineVector(x) - AffineVector(x)

    self.assertTrue(z.size() == 2)
    self.assertTrue(z.itv() == IntervalVector([[0.,0.],[0.,0.]]))

    # Subtraction of two independent affine-variable vectors
    x = AffineVariables(IntervalVector([[1.,3.],[-3.,-1.]]))
    y = AffineVariables(IntervalVector([[4.,6.],[2.,6.]]))
    z = AffineVector(x) - AffineVector(y)

    self.assertTrue(z.size() == 2)
    self.assertTrue(z.itv() == IntervalVector([[-3.,-3.],[-7.,-5.]]))

    # Unary minus
    x = AffineVariables(IntervalVector([[1.,2.],[-3.,-1.],[0.,4.]]))
    y = -AffineVector(x)

    self.assertTrue(y.size() == 3)
    self.assertTrue(y.itv() == IntervalVector([[-2.,-1.],[1.,3.],[-4.,0.]]))

    # Scalar multiplication by a positive constant
    x = AffineVariables(IntervalVector([[1.,2.],[-3.,-1.],[0.,4.]]))
    y = 2. * AffineVector(x)

    self.assertTrue(y.size() == 3)
    self.assertTrue(y.itv() == IntervalVector([[2.,4.],[-6.,-2.],[0.,8.]]))

    # Scalar division by a positive constant
    x = AffineVariables(IntervalVector([[1.,2.],[-3.,-1.],[0.,4.]]))
    y = AffineVector(x) / 2.

    self.assertTrue(y.size() == 3)
    self.assertTrue(y.itv() == IntervalVector([[0.5,1.],[-1.5,-0.5],[0.,2.]]))

    # Addition with a zero vector
    x = AffineVector(AffineVariables(IntervalVector([[-2.,3.],[4.,5.],[-7.,-1.]])))
    zero = AffineVector(AffineVariables(IntervalVector([[0.],[0.],[0.]])))

    y = x + zero
    z = zero + x

    self.assertTrue(y.itv() == x.itv())
    self.assertTrue(z.itv() == x.itv())

    # Subtraction from a zero vector
    x = AffineVector(AffineVariables(IntervalVector([[-2.,3.],[4.,5.],[-7.,-1.]])))
    zero = AffineVector(AffineVariables(IntervalVector([[0.],[0.],[0.]])))

    y = zero - x
    self.assertTrue(y.itv() == IntervalVector([[-3.,2.],[-5.,-4.],[1.,7.]]))

    # Dependency is preserved: an affine variable minus itself is exactly zero
    x = AffineVariables(IntervalVector([[-2.,3.],[4.,5.],[-7.,-1.]]))
    xv = AffineVector(x)
    y = xv - xv

    self.assertTrue(y.itv() == IntervalVector([[0.],[0.],[0.]]))

    # Dependency is preserved in a linear expression
    x = AffineVariables(IntervalVector([[1.,2.],[-3.,4.]]))
    xv = AffineVector(x)
    y = xv + 2. * xv

    self.assertTrue(y.itv() == IntervalVector([[3.,6.],[-9.,12.]]))

    # Difference of two independent affine-variable vectors
    x = AffineVector(AffineVariables(IntervalVector([[0.,2.],[4.,6.]])))
    y = AffineVector(AffineVariables(IntervalVector([[-1.,1.],[2.,3.]])))

    z1 = x - y
    z2 = y - x

    self.assertTrue(z1.itv() == IntervalVector([[1.,1.],[2.,3.]]))
    self.assertTrue(z2.itv() == -z1.itv())

    # Multiplication by zero
    x = AffineVariables(IntervalVector([[-8.,10.],[-3.,4.],[5.]]))
    y = 0. * AffineVector(x)

    self.assertTrue(y.itv() == IntervalVector([[0.],[0.],[0.]]))

    # Multiplication by a negative scalar reverses interval bounds
    x = AffineVariables(IntervalVector([[1.,2.],[-3.,-1.],[0.,4.]]))
    y = -2. * AffineVector(x)

    self.assertTrue(y.itv() == IntervalVector([[-4.,-2.],[2.,6.],[-8.,0.]]))

    # Division by a negative scalar
    x = AffineVariables(IntervalVector([[2.,4.],[-6.,-2.],[0.,8.]]))
    y = AffineVector(x) / -2.

    self.assertTrue(y.itv() == IntervalVector([[-2.,-1.],[1.,3.],[-4.,0.]]))

    # Chained linear expression inspired by ex_affineform.cpp
    x = AffineVector(AffineVariables(IntervalVector([[1.,2.],[1.,3.],[1.,2.]])))
    y = x + 2. * x
    z = y - x

    self.assertTrue(y.itv() == IntervalVector([[3.,6.],[3.,9.],[3.,6.]]))
    self.assertTrue(z.itv() == IntervalVector([[2.,4.],[2.,6.],[2.,4.]]))

    # Assignment of an expression to an existing AffineVector
    x = AffineVariables(IntervalVector([[1.,2.],[3.,4.]]))
    y = AffineVariables(IntervalVector([[-1.,1.],[5.,6.]]))
    xv, yv = AffineVector(x), AffineVector(y)

    z = AffineVector(xv)
    z = xv + yv

    self.assertTrue(z.size() == 2)
    self.assertTrue(z.itv() == IntervalVector([[0.,3.],[8.,10.]]))

    # Aliasing during assignment must not corrupt the expression
    x = AffineVariables(IntervalVector([[1.,2.],[-3.,4.]]))
    xv = AffineVector(x)
    y = AffineVector(xv)
    y = y + xv

    self.assertTrue(y.itv() == IntervalVector([[2.,4.],[-6.,8.]]))

    # A vector of degenerate intervals behaves as a real vector
    x = AffineVector(AffineVariables(Vector([1.,-2.,3.])))
    y = AffineVector(AffineVariables(Vector([4.,5.,-6.])))

    total = x + y
    difference = x - y

    self.assertTrue(total.itv() == IntervalVector([[5.],[3.],[-3.]]))
    self.assertTrue(difference.itv() == IntervalVector([[-3.],[-7.],[9.]]))

    # Complex linear expression with three independent affine vectors
    x = AffineVector(AffineVariables(IntervalVector([[1.,2.],[-2.,1.],[4.,5.]])))
    y = AffineVector(AffineVariables(IntervalVector([[-1.,1.],[3.,4.],[-2.,0.]])))
    z = AffineVector(AffineVariables(IntervalVector([[2.,6.],[-4.,2.],[1.,3.]])))

    result = 2. * x - 3. * y + z / 2.

    self.assertTrue(result.size() == 3)
    self.assertTrue(result.itv() == IntervalVector([[4.,6.],[-15.,-9.],[11.5,14.5]]))

    # Nested expression combining sums, differences and scalar operations
    x = AffineVector(AffineVariables(IntervalVector([[-2.,4.],[1.,3.]])))
    y = AffineVector(AffineVariables(IntervalVector([[5.,7.],[-4.,2.]])))

    result = (2. * (x + y) - (x - y)) / 3.

    # Algebraically, the expression is (x + 3*y) / 3.
    self.assertTrue(result == Approx_AffineVector(
      IntervalVector([[13./3.,25./3.],[-11./3.,3.]]), ERROR))

    # A complex expression must preserve affine dependencies
    x = AffineVector(AffineVariables(IntervalVector([[-5.,2.],[1.,7.],[-3.,4.]])))
    result = 4. * x - 2. * (x + x) + (x - x)

    self.assertTrue(result.itv() == IntervalVector([[0.],[0.],[0.]]))

    # Several intermediate AffineVector expressions can be chained
    x = AffineVector(AffineVariables(IntervalVector([[1.,2.],[-1.,3.],[4.,6.]])))
    y = AffineVector(AffineVariables(IntervalVector([[-2.,0.],[2.,5.],[-3.,1.]])))

    u = 2. * x - y
    v = x + 3. * y
    result = (u + v) / 3.

    # (2*x-y+x+3*y)/3 = x + 2*y/3.
    self.assertTrue(result == Approx_AffineVector(
      IntervalVector([[-1./3.,2.],[1./3.,19./3.],[2.,20./3.]]), ERROR))

  def test_product_by_interval_scalar(self):

    # Product of an affine vector by a degenerate interval
    x = AffineVector(AffineVariables(IntervalVector([[1.,2.],[-3.,-1.],[0.,4.]])))
    a = Interval(2.)
    result = x * a

    self.assertTrue(result.size() == 3)
    self.assertTrue(result.itv() == IntervalVector([[2.,4.],[-6.,-2.],[0.,8.]]))

    # Product with a degenerate interval in the opposite operand order
    x = AffineVector(AffineVariables(IntervalVector([[1.,2.],[-3.,-1.],[0.,4.]])))
    a = Interval(-2.)
    result = a * x

    self.assertTrue(result.size() == 3)
    self.assertTrue(result.itv() == IntervalVector([[-4.,-2.],[2.,6.],[-8.,0.]]))

    # Product by a non-degenerate positive interval encloses interval
    # arithmetic
    x = AffineVector(AffineVariables(IntervalVector([[1.,2.],[-3.,-1.],[0.,4.]])))
    a = Interval(2.,3.)
    result = x * a
    expected = IntervalVector([[2.,6.],[-9.,-2.],[0.,12.]])

    self.assertTrue(result.size() == 3)
    for i in range(result.size()):
      self.assertTrue(result[i].itv().is_superset(expected[i]))

    # Product by an interval containing zero
    x = AffineVector(AffineVariables(IntervalVector([[1.,2.],[-3.,-1.],[0.,4.]])))
    a = Interval(-1.,2.)
    result = a * x
    expected = IntervalVector([[-2.,4.],[-6.,3.],[-4.,8.]])

    for i in range(result.size()):
      self.assertTrue(result[i].itv().is_superset(expected[i]))

    # Empty interval propagation through a vector product
    x = AffineVector(AffineVariables(IntervalVector([[1.,2.],[-3.,-1.],[0.,4.]])))
    result = x * Interval.empty()

    self.assertTrue(result.size() == 3)
    self.assertTrue(result.is_empty())
    self.assertTrue(result.itv().is_empty())

  def test_dot_product(self):

    # Dot product with a point IntervalVector
    x = AffineVector(AffineVariables(IntervalVector([[1.,2.],[-3.,-1.],[4.,5.]])))
    weights = IntervalVector([[2.],[-1.],[3.]])

    result = affine_dot(x, weights)
    self.assertTrue(result.itv() == Interval(15.,22.))

    # Dot product with the IntervalVector on the left (same value, dot()
    # is symmetric in the reconstruction used here)
    x = AffineVector(AffineVariables(Vector([2.,-3.,4.])))
    weights = IntervalVector([[5.],[2.],[-1.]])

    result = affine_dot(weights, x)
    self.assertTrue(result.itv() == Interval(0.))

    # Dot product with non-degenerate interval coefficients
    x = AffineVector(AffineVariables(IntervalVector([[1.,2.],[-3.,-1.]])))
    weights = IntervalVector([[2.,3.],[-1.,2.]])

    result = affine_dot(x, weights)
    expected = Interval(-4.,9.)
    self.assertTrue(result.itv().is_superset(expected))

    # Dependency in a dot product: [1,-1] dot [x,x] is exactly zero
    variable = AffineVariables(IntervalVector([[-2.,5.]]))
    repeated = AffineVector([variable[0], variable[0]])
    weights = IntervalVector([[1.],[-1.]])

    result = affine_dot(repeated, weights)
    self.assertTrue(result.itv() == Interval(0.))

    # Empty coefficient propagation through a dot product
    x = AffineVector(AffineVariables(Vector([1.,2.,3.])))
    weights = IntervalVector([Interval(1.),Interval.empty(),Interval(2.)])

    result = affine_dot(x, weights)
    self.assertTrue(result.is_empty())
    self.assertTrue(result.itv().is_empty())

    # Dot product between two affine vectors, and in reverse order
    x = AffineVector(AffineVariables(Vector([2.,-3.,4.])))
    y = AffineVector(AffineVariables(Vector([5.,6.,-2.])))

    direct = affine_dot(x, y)
    reverse = affine_dot(y, x)

    # 2*5 + (-3)*6 + 4*(-2) = -16.
    self.assertTrue(direct.itv() == Interval(-16.))
    self.assertTrue(reverse.itv() == Interval(-16.))
    self.assertTrue(reverse.itv() == direct.itv())

    # Dot product between affine-variable vectors with disjoint
    # per-component ranges
    x = AffineVector(AffineVariables(IntervalVector([[1.,2.],[-3.,-1.],[4.,5.]])))
    y = AffineVector(AffineVariables(IntervalVector([[2.,3.],[-2.,1.],[-1.,2.]])))

    direct = affine_dot(x, y)
    reverse = affine_dot(y, x)
    expected = Interval(-5.5,22.)
    self.assertTrue(direct.itv().is_superset(expected))
    self.assertTrue(reverse.itv().is_superset(expected))
    self.assertTrue(direct == reverse)

    # Dot product between affine-variable vectors with overlapping
    # per-component ranges
    x = AffineVector(AffineVariables(IntervalVector([[1.,2.],[3.,10.],[4.,5.]])))
    y = AffineVector(AffineVariables(IntervalVector([[2.,3.],[2.,10.],[1.,2.]])))

    direct = affine_dot(x, y)
    reverse = affine_dot(y, x)
    expected = Interval(12,116.)
    self.assertTrue(direct.itv().is_superset(expected))
    self.assertTrue(reverse.itv().is_superset(expected))
    self.assertTrue(reverse == direct)

    # Dependency is preserved in both orders for IntervalVector
    # initialization
    x = AffineVariables(IntervalVector([[-2.,5.],[-2.,5.]]))
    repeated = AffineVector([x[0], x[0]])
    affine_weights = AffineVector(AffineVariables(IntervalVector([[1.],[-1.]])))

    direct = affine_dot(repeated, affine_weights)
    reverse = affine_dot(affine_weights, repeated)

    self.assertTrue(direct.itv() == Interval(0.))
    self.assertTrue(reverse.itv() == Interval(0.))

    # Eigen AffineT dot product (plain point AffineVector)
    lhs = AffineVector([Interval(1.0), Interval(2.0), Interval(3.0)])
    rhs = AffineVector([Interval(4.0), Interval(-1.0), Interval(2.0)])

    result = affine_dot(lhs, rhs)
    self.check_point(result, 8.0)

  def test_componentwise_nonlinear_expressions(self):

    # Component-wise nonlinear expressions on point affine variables
    x = AffineVariables(Vector([2.,3.,4.]))
    y = AffineVariables(Vector([5.,6.,8.]))
    result = AffineVector(x)

    result[0] = (x[0] + y[0]) * (y[0] - x[0])
    result[1] = x[1] * y[1] + x[1] / y[1]
    result[2] = (x[2] * x[2] + y[2]) / x[2]

    self.assertTrue(result == Approx_AffineVector(IntervalVector([[21.],[18.5],[6.]]), ERROR))

    # Quadratic polynomial (implicit conversion from AffineVariables)
    x = AffineVariables(Vector([-2.,3.,4.]))
    result = AffineVector(x)

    for i in range(x.size()):
      result[i] = x[i] * x[i] + 2. * x[i] + 1.

    self.assertTrue(result.size() == 3)
    self.assertTrue(result.itv() == IntervalVector([[1.],[16.],[25.]]))

    # Products between vector components
    x = AffineVariables(Vector([2.,-3.,4.]))
    y = AffineVariables(Vector([5.,6.,-2.]))
    result = AffineVector(x)

    result[0] = x[0] * y[1] + x[1] * y[0]
    result[1] = (x[0] + y[0]) * (x[1] - y[1])
    result[2] = x[2] * y[2] - x[0] * x[1]

    self.assertTrue(result.size() == 3)
    self.assertTrue(result.itv() == IntervalVector([[-3.],[-63.],[-2.]]))

    # Rational expressions with point denominators
    x = AffineVariables(Vector([2.,3.,4.]))
    y = AffineVariables(Vector([5.,6.,8.]))
    result = AffineVector(x)

    result[0] = (x[0] * x[0] + y[0]) / x[0]
    result[1] = (x[1] + y[1]) / (y[1] - x[1])
    result[2] = (x[2] * y[2]) / (x[2] + y[2])

    self.assertTrue(result.size() == 3)
    self.assertTrue(result == Approx_AffineVector(IntervalVector([[4.5],[3.],[8./3.]]), ERROR))

    # Cubic and mixed polynomial terms
    x = AffineVariables(Vector([-2.,3.,4.]))
    y = AffineVariables(Vector([1.,-1.,2.]))
    result = AffineVector(x)

    result[0] = x[0] * x[0] * x[0] - y[0]
    result[1] = x[1] * x[1] * y[1] + 2. * x[1]
    result[2] = (x[2] + y[2]) * (x[2] - y[2]) + y[2] * y[2]

    self.assertTrue(result.size() == 3)
    self.assertTrue(result == Approx_AffineVector(IntervalVector([[-9.],[-3.],[16.]]), ERROR))

  def test_componentwise_elementary_functions(self):

    # Component-wise sqr() and abs() on point (degenerate) affine variables
    x = AffineVector(AffineVariables(Vector([2.,-3.,-0.5])))
    result = AffineVector(x)

    result[0] = sqr(x[0])
    result[1] = sqr(x[1])
    result[2] = abs(x[2])

    self.assertTrue(result == Approx_AffineVector(IntervalVector([[4.],[9.],[0.5]]), ERROR))

    # Component-wise sqrt() encloses the exact interval-arithmetic result
    x = AffineVector(AffineVariables(IntervalVector([[4.,9.],[1.,4.]])))
    result = AffineVector(x)

    result[0] = sqrt(x[0])
    result[1] = sqrt(x[1])

    self.assertTrue(result[0].itv().is_superset(Interval(2.,3.)))
    self.assertTrue(result[1].itv().is_superset(Interval(1.,2.)))

    # Component-wise exp() and log() enclose the exact interval-arithmetic
    # result
    x = AffineVector(AffineVariables(IntervalVector([[0.,1.],[1.,2.]])))
    result = AffineVector(x)

    result[0] = exp(x[0])
    result[1] = log(x[1])

    self.assertTrue(result[0].itv().is_superset(Interval(math.exp(0.),math.exp(1.))))
    self.assertTrue(result[1].itv().is_superset(Interval(math.log(1.),math.log(2.))))

    # Component-wise sin() and cos(), on a range where both are monotonic
    x = AffineVector(AffineVariables(IntervalVector([[0.,0.5],[0.,0.5]])))
    result = AffineVector(x)

    result[0] = sin(x[0])
    result[1] = cos(x[1])

    self.assertTrue(result[0].itv().is_superset(Interval(math.sin(0.),math.sin(0.5))))
    self.assertTrue(result[1].itv().is_superset(Interval(math.cos(0.5),math.cos(0.))))

  def test_componentwise_min_max_set_operations(self):

    a = AffineVector(AffineVariables(Vector([1.,5.,3.])))
    b = AffineVector(AffineVariables(Vector([4.,2.,3.])))

    result_min = IntervalVector(3)
    result_max = IntervalVector(3)
    result_union = IntervalVector(3)
    for i in range(a.size()):
      result_min[i] = min(a[i], b[i])
      result_max[i] = max(a[i], b[i])
      result_union[i] = a[i] | b[i]
    inter_0 = a[0] & b[0] # disjoint point values
    inter_2 = a[2] & b[2] # equal point values

    self.assertTrue(result_min == IntervalVector([[1.],[2.],[3.]]))
    self.assertTrue(result_max == IntervalVector([[4.],[5.],[3.]]))
    self.assertTrue(result_union == IntervalVector([[1.,4.],[2.,5.],[3.]]))
    self.assertTrue(inter_0.is_empty())
    self.assertTrue(inter_2 == Interval(3.))

  def test_addon_lb_ub_mid(self):

    values = IntervalVector([[-1.0, 1.0], [2.0, 3.0], [4.0, 5.0]])
    x = AffineVariables(values)

    self.assertTrue(x.lb() == Vector([-1.0, 2.0, 4.0]))
    self.assertTrue(x.ub() == Vector([1.0, 3.0, 5.0]))
    self.assertTrue(x.size() == 3)
    self.assertTrue(x.mid() == Vector([0.0, 2.5, 4.5]))

    y = AffineVector(x)
    self.assertTrue(y.lb() == Vector([-1.0, 2.0, 4.0]))
    self.assertTrue(y.ub() == Vector([1.0, 3.0, 5.0]))
    self.assertTrue(y.size() == 3)
    self.assertTrue(y.mid() == Vector([0.0, 2.5, 4.5]))

    y.set_empty()
    self.assertTrue(y.is_empty())

    x.init(Interval.empty()) # broadcasts, but every component ends up
                              # empty either way; matches x.set_empty()
    self.assertTrue(x.is_empty())

    x = AffineVariables(values)
    self.assertTrue(x.itv() == values)

  def test_vector_indexing(self):

    vector = make_point_vector_5()

    self.assertTrue(vector.rows() == 5)
    self.assertTrue(vector.cols() == 1)
    self.assertTrue(vector.size() == 5)

    self.check_point(vector[0], 1.0)
    self.check_point(vector[2], 3.0)
    self.check_point(vector[4], 5.0)

    # coeffRef() is not exposed: __setitem__ is the equivalent in-place
    # single-component mutation.
    vector[3] = Interval(-4.0)
    self.check_point(vector[3], -4.0)

  def test_vector_segment_head_and_tail(self):

    vector = make_point_vector_5()

    # head()/tail() are not exposed; subvector(start,end) (inclusive
    # bounds) is the equivalent of segment(start,len) here.
    segment = vector.subvector(1, 3)
    head = vector.subvector(0, 1)
    tail = vector.subvector(3, 4)

    self.assertTrue(segment.size() == 3)
    self.check_point(segment[0], 2.0)
    self.check_point(segment[1], 3.0)
    self.check_point(segment[2], 4.0)

    self.check_point(head[0], 1.0)
    self.check_point(head[1], 2.0)
    self.check_point(tail[0], 4.0)
    self.check_point(tail[1], 5.0)

  def test_segment_assignment(self):

    vector = AffineVector.zero(6)

    source = AffineVector(AffineVariables(IntervalVector([[2],[4],[6]])))
    # There is no lvalue segment(): put(start_id, x) is the exposed
    # equivalent of "vector.segment(2, 3) = source;".
    vector.put(2, source)

    self.check_point(vector[0], 0.0)
    self.check_point(vector[1], 0.0)
    self.check_point(vector[2], 2.0)
    self.check_point(vector[3], 4.0)
    self.check_point(vector[4], 6.0)
    self.check_point(vector[5], 0.0)

  def test_vector_transpose_to_row(self):

    vector = make_point_vector_5()
    # transpose() returns a generic AffineMatrix (1 x n), not an AffineRow.
    row = vector.transpose()

    self.assertTrue(row.rows() == 1)
    self.assertTrue(row.cols() == 5)
    self.check_point(row[(0,0)], 1.0)
    self.check_point(row[(0,4)], 5.0)

  def test_swap(self):

    lhs = AffineVector(AffineVariables(IntervalVector([[1],[2],[3]])))
    rhs = AffineVector(AffineVariables(IntervalVector([[4],[5],[6]])))

    # swap() is not exposed; Python has no reference semantics to swap
    # storage in place, so this rebinds the two names instead, which is
    # observationally the same for this test (only the resulting values
    # are checked, not object identity).
    lhs, rhs = rhs, lhs

    self.check_point(lhs[0], 4.0)
    self.check_point(lhs[2], 6.0)
    self.check_point(rhs[0], 1.0)
    self.check_point(rhs[2], 3.0)

  def test_compound_assignment(self):

    value = AffineVector(AffineVariables(IntervalVector([[1],[2],[3]])))
    other = AffineVector(AffineVariables(IntervalVector([[4],[5],[6]])))

    value += other
    self.check_point(value[0], 5.0)
    self.check_point(value[2], 9.0)

    value -= other
    self.check_point(value[0], 1.0)
    self.check_point(value[2], 3.0)

    # array() is not exposed: +=/-= only accept another AffineVector (see
    # the module docstring), so a same-size constant vector stands in for
    # the scalar broadcast; *= and /= already accept a plain scalar.
    value += AffineVector.constant(3, Affine(3.0))
    self.check_point(value[0], 4.0)
    self.check_point(value[2], 6.0)

    value -= AffineVector.constant(3, Affine(5.0))
    self.check_point(value[0], -1.0)
    self.check_point(value[2], 1.0)

    value *= 2.0
    self.check_point(value[0], -2.0)
    self.check_point(value[2], 2.0)

    value /= 2.0
    self.check_point(value[0], -1.0)
    self.check_point(value[2], 1.0)

  def test_compound_assignment_with_affine_variables(self):

    value_variables = AffineVariables(IntervalVector([[1,2], [-2,5], [3,6]]))
    value = AffineVector(value_variables)

    value -= value_variables
    self.check_point(value[0], 0)
    self.check_point(value[1], 0)
    self.check_point(value[2], 0)

    value += value_variables
    self.check_interval(value[0], value_variables[0].itv())
    self.check_interval(value[1], value_variables[1].itv())
    self.check_interval(value[2], value_variables[2].itv())

  def test_scalar_compound_assignments_intervalvector(self):

    value = IntervalVector([[1.0], [-2.0], [4.0]])

    value += IntervalVector.constant(3, Interval(2.0))
    self.assertTrue(value[0] == Interval(3.0))
    self.assertTrue(value[1] == Interval(0.0))
    self.assertTrue(value[2] == Interval(6.0))

    value -= IntervalVector.constant(3, Interval(2.0))
    self.assertTrue(value[0] == Interval(1.0))
    self.assertTrue(value[1] == Interval(-2.0))
    self.assertTrue(value[2] == Interval(4.0))

    value *= 2.0
    self.assertTrue(value[0] == Interval(2.0))
    self.assertTrue(value[1] == Interval(-4.0))
    self.assertTrue(value[2] == Interval(8.0))

    value /= 2.0
    self.assertTrue(value[0] == Interval(1.0))
    self.assertTrue(value[1] == Interval(-2.0))
    self.assertTrue(value[2] == Interval(4.0))

  def test_scalar_compound_assignments_affinevector(self):

    value = AffineVector(3)
    value[0] = Interval(1.0)
    value[1] = Interval(-2.0)
    value[2] = Interval(4.0)

    value += AffineVector.constant(3, Affine(2.0))
    self.check_point(value[0], 3.0)
    self.check_point(value[1], 0.0)
    self.check_point(value[2], 6.0)

    value -= AffineVector.constant(3, Affine(2.0))
    self.check_point(value[0], 1.0)
    self.check_point(value[1], -2.0)
    self.check_point(value[2], 4.0)

    value *= 2.0
    self.check_point(value[0], 2.0)
    self.check_point(value[1], -4.0)
    self.check_point(value[2], 8.0)

    value /= 2.0
    self.check_point(value[0], 1.0)
    self.check_point(value[1], -2.0)
    self.check_point(value[2], 4.0)

  def test_unary_minus(self):

    vector = make_point_vector_5()
    negative = -vector

    self.check_point(negative[0], -1.0)
    self.check_point(negative[2], -3.0)
    self.check_point(negative[4], -5.0)

  def test_cwise_product(self):

    lhs = AffineVector(AffineVariables(IntervalVector([[1],[-2],[3]])))
    rhs = AffineVector(AffineVariables(IntervalVector([[4],[5],[-6]])))

    result = affine_cwise_product(lhs, rhs)
    self.check_point(result[0], 4.0)
    self.check_point(result[1], -10.0)
    self.check_point(result[2], -18.0)

  def test_cwise_quotient(self):

    numerator = AffineVector(AffineVariables(IntervalVector([[4],[-10],[18]])))
    denominator = AffineVector(AffineVariables(IntervalVector([[2],[5],[-3]])))

    result = affine_cwise_quotient(numerator, denominator)
    self.check_point(result[0], 2.0)
    self.check_point(result[1], -2.0)
    self.check_point(result[2], -6.0)

  def test_sum(self):

    vector = make_point_vector_5()
    result = affine_sum(vector)
    self.check_point(result, 15.0)

  def test_reverse(self):

    vector = make_point_vector_5()
    reversed_vector = affine_reverse(vector)

    self.check_point(reversed_vector[0], 5.0)
    self.check_point(reversed_vector[1], 4.0)
    self.check_point(reversed_vector[4], 1.0)

  def test_interval_enclosure_through_linear_expression(self):

    variables = AffineVariables(IntervalVector([[-2.0, 1.0], [3.0, 5.0], [-1.0, 4.0]]))
    x = AffineVector(variables)

    result = 2.0 * x - x

    for i in range(result.size()):
      self.check_interval_enclosure(result[i], variables[i].itv())

  def test_eval_materializes_expressions(self):

    # Python operators already return concrete AffineVector objects
    # immediately (there is no separate lazy expression type to
    # materialize with an .eval() call).
    lhs = AffineVector(AffineVariables(IntervalVector([[1],[2],[3]])))
    rhs = AffineVector(AffineVariables(IntervalVector([[4],[5],[6]])))

    evaluated = lhs + rhs

    self.check_point(evaluated[0], 5.0)
    self.check_point(evaluated[1], 7.0)
    self.check_point(evaluated[2], 9.0)

  def test_stream_output(self):

    # AffineVector, non-empty
    v = make_point_vector_5()
    s = str(v)
    self.assertTrue(len(s) != 0)
    self.assertTrue(s.find("empty") == -1)

    # AffineVector, empty
    v = AffineVector(2)
    v[0] = Interval.empty()
    self.assertTrue(str(v) == "[ empty 2d box ]")

    # AffineRow, non-empty
    row = AffineRow(3)
    row[0] = Interval(1.0)
    row[1] = Interval(2.0)
    row[2] = Interval(3.0)
    s = str(row)
    self.assertTrue(len(s) != 0)
    self.assertTrue(s.find("empty") == -1)

    # AffineRow, empty
    row = AffineRow(2)
    row[0] = Interval.empty()
    self.assertTrue(str(row) == "[ empty row ]")


if __name__ ==  '__main__':
  unittest.main()
