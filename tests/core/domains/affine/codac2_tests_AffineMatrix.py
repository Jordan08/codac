#!/usr/bin/env python

# Codac tests
#
# Unit tests for AffineMatrix (the Python binding of AffineMainMatrix<T>):
# construction/access/resize, block/row/column/diagonal views and
# assignment, transpose, factories (zero/ones/constant/eye), compound
# assignments, and every product involving a matrix or a vector-by-row
# outer product (real, Interval and Affine operands).
#
# Matching Python/NumPy convention, every matrix/vector product (matrix *
# matrix, matrix * vector, row * matrix, row . vector) is bound on "@"
# (__matmul__/__rmatmul__), not "*": "*" is reserved for scalar/elementwise
# multiplication (see codac2_py_AffineMatrixBase.h).
#
# codac2_py_AffineMatrixBase.h exposes a hand-picked, reduced subset of
# Eigen's API (documented as such in its own header comment) rather than
# the full generic MatrixBase binding: there is no coeffRef()/
# topLeftCorner()/bottomRightCorner()/middleRows()/middleCols()/
# transposeInPlace()/diagonal()/trace()/replicate()/dot()/array(), no
# lvalue block()/row()/col() (only the value-returning forms, paired with
# set_block()/set_row()/set_col()), and no fixed-size Eigen container
# (Eigen::Matrix<Affine,N,M>) at all. Where a C++ test relies on one of
# these, this file either reconstructs the same observable behavior from
# the bound primitives (indexing, block(), +, @) or, for fixed-size
# containers and Eigen lazy-expression-template mechanics (which have no
# meaning in Python: every Python expression is already evaluated
# immediately), is not translated.
#
# A read of an Affine matrix entry, matrix[(i,j)], uses a tuple; a write,
# matrix[[i,j]] = value, uses a list -- these are genuinely different
# __getitem__/__setitem__ overloads in the binding, not a typo. That read
# returns a live reference to the entry, not a copy (like Affine/
# AffineVector/AffineVariables indexing): naming it and then mutating
# that name in place (e.g. via .set_empty()) reaches back into the
# matrix, exactly as the equivalent C++ reference would.
#
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
from codac import *

ERROR = sys.float_info.epsilon*100


def make_point_matrix_3x4():
  matrix = AffineMatrix(3, 4)
  matrix[[0,0]] = Interval(1.0); matrix[[0,1]] = Interval(2.0)
  matrix[[0,2]] = Interval(3.0); matrix[[0,3]] = Interval(4.0)
  matrix[[1,0]] = Interval(5.0); matrix[[1,1]] = Interval(6.0)
  matrix[[1,2]] = Interval(7.0); matrix[[1,3]] = Interval(8.0)
  matrix[[2,0]] = Interval(9.0); matrix[[2,1]] = Interval(10.0)
  matrix[[2,2]] = Interval(11.0); matrix[[2,3]] = Interval(12.0)
  return matrix

def interval_hull_matrix(value):
  result = IntervalMatrix(value.rows(), value.cols())
  for i in range(value.rows()):
    for j in range(value.cols()):
      result[[i,j]] = value[(i,j)].itv()
  return result

def interval_hull_row_matrix(row):
  # 1 x n IntervalMatrix hull of a 1-D indexable AffineRow: interval_hull_matrix()
  # cannot be reused directly since AffineRow only supports 1-D indexing (row[j]),
  # not the (i,j) tuple indexing of a 2-D AffineMatrix.
  n = len(row)
  result = IntervalMatrix(1, n)
  for j in range(n):
    result[[0,j]] = row[j].itv()
  return result

def interval_hull_vector(value):
  # Works for AffineVector, AffineRow or AffineVariables alike: all three
  # use the same 1-D [i] indexing.
  result = IntervalVector(len(value))
  for i in range(len(value)):
    result[i] = value[i].itv()
  return result

def multiply_intervals(lhs, rhs):
  return lhs * rhs

def make_affine_matrix_2x3():
  value = AffineMatrix(2, 3)
  value[[0,0]] = Interval(1.0);  value[[0,1]] = Interval(-2.0); value[[0,2]] = Interval(3.0)
  value[[1,0]] = Interval(4.0);  value[[1,1]] = Interval(0.5);  value[[1,2]] = Interval(-1.0)
  return value

def make_real_matrix_2x3():
  return Matrix([[1.0,-2.0,3.0],[4.0,0.5,-1.0]])

def make_interval_matrix_2x3():
  return IntervalMatrix([
    [Interval(0.9,1.1), Interval(-2.1,-1.9), Interval(2.9,3.1)],
    [Interval(3.9,4.1), Interval(0.4,0.6), Interval(-1.1,-0.9)]
  ])

def make_variable_vector_3():
  return AffineVariables(IntervalVector([[1.0,2.0],[-1.0,1.0],[3.0,4.0]]))

def make_affine_vector_3():
  return AffineVector(make_variable_vector_3())

def make_affine_row_3():
  return AffineRow([Interval(1.0,2.0), Interval(-1.0,1.0), Interval(3.0,4.0)])

def affine_diagonal(matrix):
  # _min, not min: "from codac import *" shadows the plain Python min()
  # with codac's Interval/Affine-valued min(), which would turn n into an
  # Interval and break range(n) below.
  n = _min(matrix.rows(), matrix.cols())
  return AffineVector([matrix[(i,i)] for i in range(n)])

def affine_set_diagonal(matrix, values):
  n = _min(matrix.rows(), matrix.cols())
  for i in range(n):
    matrix[[i,i]] = values[i]

def affine_trace(matrix):
  n = _min(matrix.rows(), matrix.cols())
  result = matrix[(0,0)]
  for i in range(1, n):
    result = result + matrix[(i,i)]
  return result

def affine_replicate(row, row_times, col_times):
  # AffineMatrix's list constructor takes a Sequence[AffineVector], one
  # AffineVector per row (see codac2_py_AffineMatrix.cpp): an AffineRow is
  # not accepted there, hence AffineVector rather than AffineRow here.
  n = len(row)
  return AffineMatrix([
    AffineVector([row[j % n] for j in range(n*col_times)]) for _ in range(row_times)
  ])

def affine_outer_product(col, row):
  return AffineMatrix([
    AffineVector([col[i]*row[j] for j in range(len(row))]) for i in range(len(col))
  ])

def _as_interval(v):
  return v.itv() if hasattr(v, 'itv') else Interval(v)

def interval_outer_product(col, row):
  # The n x m IntervalMatrix outer product of two 1-D indexable
  # containers (Row, IntervalRow, AffineRow, Vector, IntervalVector,
  # AffineVector or AffineVariables all qualify): built from itv()-wrapped
  # col[i]*row[j], so it never depends on a Matrix/IntervalMatrix
  # constructor accepting a Row/IntervalRow/AffineRow argument.
  n, m = len(col), len(row)
  result = IntervalMatrix(n, m)
  for i in range(n):
    for j in range(m):
      result[[i,j]] = _as_interval(col[i]) * _as_interval(row[j])
  return result


class TestAffineMatrix(unittest.TestCase):

  def check_point(self, actual, expected):
    self.assertTrue(actual == Approx_Affine(Interval(expected), ERROR))

  def check_encloses_matrix(self, actual, expected):
    self.assertTrue(actual.rows() == expected.rows())
    self.assertTrue(actual.cols() == expected.cols())
    for i in range(actual.rows()):
      for j in range(actual.cols()):
        self.assertTrue(actual[(i,j)].itv().is_superset(expected[(i,j)]))

  def check_encloses_vector(self, actual, expected):
    self.assertTrue(len(actual) == expected.size())
    for i in range(len(actual)):
      self.assertTrue(actual[i].itv().is_superset(expected[i]))

  def test_construction_access_and_assignment(self):

    a = AffineMatrix(2, 3)
    self.assertTrue(a.rows() == 2)
    self.assertTrue(a.cols() == 3)

    a[[0,0]] = Interval(1.,2.); a[[0,1]] = Interval(3.); a[[0,2]] = Interval(-1.,1.)
    a[[1,0]] = Interval(-2.);   a[[1,1]] = Interval(4.,5.); a[[1,2]] = Interval(0.)

    self.assertTrue(a[(0,0)].itv() == Interval(1.,2.))
    self.assertTrue(a[(0,1)].itv() == Interval(3.))
    self.assertTrue(a[(1,1)].itv() == Interval(4.,5.))

    b = AffineMatrix(a)
    self.assertTrue(b.rows() == 2)
    self.assertTrue(b.cols() == 3)
    self.assertTrue(b.itv() == a.itv())

    b[[0,0]] = Interval(10.,11.)
    self.assertTrue(a[(0,0)].itv() == Interval(1.,2.))
    self.assertTrue(b[(0,0)].itv() == Interval(10.,11.))

    c = AffineMatrix(a) # "c(1,1); c = a;" reassignment: Python has no
                         # assignment operator distinct from construction.
    self.assertTrue(c.itv() == a.itv())

  def test_resize_preserves_coefficients(self):

    a = AffineMatrix(2,2)
    a[[0,0]]=Interval(1.); a[[0,1]]=Interval(2.); a[[1,0]]=Interval(3.); a[[1,1]]=Interval(4.)

    a.conservativeResize(3,3)
    self.assertTrue(a.rows() == 3)
    self.assertTrue(a.cols() == 3)
    self.assertTrue(a[(0,0)].itv() == Interval(1.))
    self.assertTrue(a[(0,1)].itv() == Interval(2.))
    self.assertTrue(a[(1,0)].itv() == Interval(3.))
    self.assertTrue(a[(1,1)].itv() == Interval(4.))

    a[[2,0]]=Interval(5.); a[[2,1]]=Interval(6.); a[[2,2]]=Interval(7.)
    a[[0,2]]=Interval(8.); a[[1,2]]=Interval(9.)
    a.conservativeResize(2,1)
    self.assertTrue(a.rows() == 2)
    self.assertTrue(a.cols() == 1)
    self.assertTrue(a[(0,0)].itv() == Interval(1.))
    self.assertTrue(a[(1,0)].itv() == Interval(3.))

    a.resize(4,2)
    self.assertTrue(a.rows() == 4)
    self.assertTrue(a.cols() == 2)

  def test_empty_coefficients_and_stream_output(self):

    a = AffineMatrix(2,2)
    a[[0,0]]=Interval(1.); a[[0,1]]=Interval(2.); a[[1,0]]=Interval(3.); a[[1,1]]=Interval(4.)
    # __getitem__ returns a live reference (reference_internal), just
    # like C++'s a(0,1): mutating it mutates the matrix in place.
    a[(0,1)].set_empty()

    self.assertTrue(a.is_empty())
    self.assertTrue(a.itv()[(0,0)] == Interval(1.))
    self.assertTrue(a.itv()[(0,1)].is_empty())

    self.assertTrue(str(a) == "[ empty 2x2 mat ]")

    a.set_empty()
    for i in range(a.rows()):
      for j in range(a.cols()):
        self.assertTrue(a[(i,j)].is_empty())

  def test_linear_addition_and_subtraction(self):

    a = AffineMatrix(2,2); b = AffineMatrix(2,2)
    a[[0,0]]=Interval(1.,2.); a[[0,1]]=Interval(-1.,1.)
    a[[1,0]]=Interval(3.);    a[[1,1]]=Interval(4.,5.)
    b[[0,0]]=Interval(2.,3.); b[[0,1]]=Interval(2.)
    b[[1,0]]=Interval(-2.,0.); b[[1,1]]=Interval(1.,2.)

    total = a + b
    difference = a - b
    for i in range(2):
      for j in range(2):
        self.assertTrue(total[(i,j)].itv().is_superset(a[(i,j)].itv()+b[(i,j)].itv()))
        self.assertTrue(difference[(i,j)].itv().is_superset(a[(i,j)].itv()-b[(i,j)].itv()))

    variables = AffineVariables(IntervalVector([[1.,2.],[-1.,1.]]))
    a.set_col(0, variables)
    a.set_col(1, variables)
    zero = a - a

    for i in range(2):
      for j in range(2):
        self.assertTrue(zero[(i,j)].itv() == Interval(0.))

  def test_linear_addition_with_row_and_segment_assignment(self):

    a = AffineMatrix(2,2); b = AffineMatrix(2,2)
    v = AffineVariables(IntervalVector([[1.,2.],[-1.,1.],[3.,4.],[5.,6.]]))
    a.set_row(0, AffineRow([v[0], v[1]]))
    a.set_row(1, AffineRow([v[2], v[3]]))
    b[[0,0]] = v[0]
    b[[0,1]] = v[1]
    b[[1,0]] = v[2]
    b[[1,1]] = v[3]

    zero = a - b

    for i in range(2):
      for j in range(2):
        self.assertTrue(zero[(i,j)].itv() == Interval(0.))

  def test_scalar_operations(self):

    a = AffineMatrix(2,2)
    a[[0,0]]=Interval(1.,2.); a[[0,1]]=Interval(-3.,-1.)
    a[[1,0]]=Interval(0.,4.); a[[1,1]]=Interval(5.)

    twice = 2. * a
    negative = a * -3.
    half = a / 2.
    for i in range(2):
      for j in range(2):
        self.assertTrue(twice[(i,j)].itv().is_superset(2.*a[(i,j)].itv()))
        self.assertTrue(negative[(i,j)].itv().is_superset(-3.*a[(i,j)].itv()))
        self.assertTrue(half[(i,j)].itv().is_superset(a[(i,j)].itv()/2.))

  def test_matrix_product(self):

    a = AffineMatrix(2,3); b = AffineMatrix(3,2)
    av = [[1.,2.,3.],[-1.,4.,2.]]
    bv = [[2.,1.],[0.,-1.],[3.,2.]]
    for i in range(2):
      for j in range(3):
        a[[i,j]] = Interval(av[i][j])
    for i in range(3):
      for j in range(2):
        b[[i,j]] = Interval(bv[i][j])

    result = a @ b
    self.assertTrue(result.rows()==2)
    self.assertTrue(result.cols()==2)
    self.assertTrue(result[(0,0)].itv() == Interval(11.))
    self.assertTrue(result[(0,1)].itv() == Interval(5.))
    self.assertTrue(result[(1,0)].itv() == Interval(4.))
    self.assertTrue(result[(1,1)].itv() == Interval(-1.))

  def test_matrix_times_affinevector(self):

    a = AffineMatrix(2,3)
    a[[0,0]]=Interval(1.); a[[0,1]]=Interval(2.);  a[[0,2]]=Interval(-1.)
    a[[1,0]]=Interval(3.); a[[1,1]]=Interval(-2.); a[[1,2]]=Interval(4.)

    x = AffineVector(AffineVariables(IntervalVector([[1.,2.],[-1.,1.],[3.,4.]])))
    y = a @ x

    self.assertTrue(y.size()==2)
    expected0 = x[0].itv()+2.*x[1].itv()-x[2].itv()
    expected1 = 3.*x[0].itv()-2.*x[1].itv()+4.*x[2].itv()
    self.assertTrue(y[0].itv().is_superset(expected0))
    self.assertTrue(y[1].itv().is_superset(expected1))

  def test_matrix_times_affinevariables(self):

    a = AffineMatrix(2,2)
    a[[0,0]]=Interval(1.); a[[0,1]]=Interval(-1.)
    a[[1,0]]=Interval(2.); a[[1,1]]=Interval(3.)

    x = AffineVariables(IntervalVector([[-2.,5.],[1.,4.]]))
    y = a @ x

    self.assertTrue(y.size()==2)
    self.assertTrue(y[0].itv().is_superset(x[0].itv()-x[1].itv()))
    self.assertTrue(y[1].itv().is_superset(2.*x[0].itv()+3.*x[1].itv()))

  def test_vector_row_times_matrix(self):

    a = AffineMatrix(2,3)
    a[[0,0]]=Interval(1.); a[[0,1]]=Interval(2.); a[[0,2]]=Interval(3.)
    a[[1,0]]=Interval(4.); a[[1,1]]=Interval(5.); a[[1,2]]=Interval(6.)

    x = AffineVector(AffineVariables(Vector([2.,-1.])))
    y = x.transpose() @ a

    self.assertTrue(y[(0,0)].itv() == Interval(-2.))
    self.assertTrue(y[(0,1)].itv() == Interval(-1.))
    self.assertTrue(y[(0,2)].itv() == Interval(0.))

  def test_chained_linear_expression(self):

    a = AffineMatrix(2,2); b = AffineMatrix(2,2); c = AffineMatrix(2,2)
    a[[0,0]]=Interval(1.); a[[0,1]]=Interval(2.); a[[1,0]]=Interval(3.); a[[1,1]]=Interval(4.)
    b[[0,0]]=Interval(2.); b[[0,1]]=Interval(-1.); b[[1,0]]=Interval(0.); b[[1,1]]=Interval(5.)
    c[[0,0]]=Interval(-1.); c[[0,1]]=Interval(3.); c[[1,0]]=Interval(2.); c[[1,1]]=Interval(1.)

    result = 2.*a - 3.*b + c
    self.assertTrue(result[(0,0)].itv() == Interval(-5.))
    self.assertTrue(result[(0,1)].itv() == Interval(10.))
    self.assertTrue(result[(1,0)].itv() == Interval(8.))
    self.assertTrue(result[(1,1)].itv() == Interval(-6.))

  def test_matrix_dimensions_and_coefficient_access(self):

    matrix = make_point_matrix_3x4()

    self.assertTrue(matrix.rows() == 3)
    self.assertTrue(matrix.cols() == 4)
    self.assertTrue(matrix.size() == 12)

    self.check_point(matrix[(0,0)], 1.0)
    self.check_point(matrix[(1,2)], 7.0)
    self.check_point(matrix[(2,3)], 12.0)

    # coeffRef() is not exposed: __setitem__ is the equivalent in-place
    # single-component mutation.
    matrix[[1,1]] = Interval(-6.0)
    self.check_point(matrix[(1,1)], -6.0)

  def test_matrix_row_and_column_views(self):

    matrix = make_point_matrix_3x4()

    column = matrix.col(2)
    self.assertTrue(column.size() == 3)
    self.check_point(column[0], 3.0)
    self.check_point(column[1], 7.0)
    self.check_point(column[2], 11.0)

    row = matrix.row(1)
    self.assertTrue(row.size() == 4)
    self.check_point(row[0], 5.0)
    self.check_point(row[1], 6.0)
    self.check_point(row[2], 7.0)
    self.check_point(row[3], 8.0)

  def test_row_and_column_assignments(self):

    matrix = AffineMatrix(3, 3)
    matrix.init(Affine(0.0))

    column_variables = AffineVariables(IntervalVector([[1],[2],[3]]))
    matrix.set_col(1, column_variables)

    self.check_point(matrix[(0,1)], 1.0)
    self.check_point(matrix[(1,1)], 2.0)
    self.check_point(matrix[(2,1)], 3.0)

    row_variables = AffineVariables(IntervalVector([[4],[5],[6]]))
    row = AffineRow([row_variables[0], row_variables[1], row_variables[2]])
    matrix.set_row(2, row)

    self.check_point(matrix[(2,0)], 4.0)
    self.check_point(matrix[(2,1)], 5.0)
    self.check_point(matrix[(2,2)], 6.0)

  def test_affinevariables_assignment_to_matrix_column(self):

    matrix = AffineMatrix(3, 2)
    matrix.init(Affine(0.0))

    variables = AffineVariables(IntervalVector([[-2.0, 1.0], [3.0, 4.0], [5.0, 7.0]]))
    matrix.set_col(0, variables)

    for i in range(variables.size()):
      self.assertTrue(matrix[(i,0)].itv() == variables[i].itv())

  def test_block_extraction(self):

    matrix = make_point_matrix_3x4()
    block = matrix.block(1, 1, 2, 3)

    self.assertTrue(block.rows() == 2)
    self.assertTrue(block.cols() == 3)

    self.check_point(block[(0,0)], 6.0)
    self.check_point(block[(0,2)], 8.0)
    self.check_point(block[(1,0)], 10.0)
    self.check_point(block[(1,2)], 12.0)

  def test_block_assignment(self):

    matrix = AffineMatrix(4, 5)
    matrix.init(Affine(0.0))

    block = AffineMatrix(2, 3)
    block[[0,0]] = Interval(1.0); block[[0,1]] = Interval(2.0); block[[0,2]] = Interval(3.0)
    block[[1,0]] = Interval(4.0); block[[1,1]] = Interval(5.0); block[[1,2]] = Interval(6.0)

    matrix.set_block(1, 1, 2, 3, block)

    self.check_point(matrix[(1,1)], 1.0)
    self.check_point(matrix[(1,3)], 3.0)
    self.check_point(matrix[(2,1)], 4.0)
    self.check_point(matrix[(2,3)], 6.0)
    self.check_point(matrix[(0,0)], 0.0)
    self.check_point(matrix[(3,4)], 0.0)

  def test_corners(self):

    # topLeftCorner()/bottomRightCorner() are not exposed: block() is
    # the equivalent, at the appropriate offset.
    matrix = make_point_matrix_3x4()

    top_left = matrix.block(0, 0, 2, 2)
    bottom_right = matrix.block(matrix.rows()-2, matrix.cols()-2, 2, 2)

    self.check_point(top_left[(0,0)], 1.0)
    self.check_point(top_left[(1,1)], 6.0)
    self.check_point(bottom_right[(0,0)], 7.0)
    self.check_point(bottom_right[(1,1)], 12.0)

  def test_middle_rows_and_columns(self):

    # middleRows()/middleCols() are not exposed: block() is the
    # equivalent, spanning the full orthogonal dimension.
    matrix = make_point_matrix_3x4()

    middle_columns = matrix.block(0, 1, matrix.rows(), 2)
    self.assertTrue(middle_columns.rows() == 3)
    self.assertTrue(middle_columns.cols() == 2)
    self.check_point(middle_columns[(0,0)], 2.0)
    self.check_point(middle_columns[(2,1)], 11.0)

    middle_rows = matrix.block(1, 0, 2, matrix.cols())
    self.assertTrue(middle_rows.rows() == 2)
    self.assertTrue(middle_rows.cols() == 4)
    self.check_point(middle_rows[(0,0)], 5.0)
    self.check_point(middle_rows[(1,3)], 12.0)

  def test_transpose(self):

    matrix = make_point_matrix_3x4()
    transpose = matrix.transpose()

    self.assertTrue(transpose.rows() == 4)
    self.assertTrue(transpose.cols() == 3)

    self.check_point(transpose[(0,0)], 1.0)
    self.check_point(transpose[(2,1)], 7.0)
    self.check_point(transpose[(3,2)], 12.0)

  def test_in_place_transpose_for_square_matrix(self):

    matrix = AffineMatrix(3, 3)
    matrix[[0,0]] = Interval(1.0); matrix[[0,1]] = Interval(2.0); matrix[[0,2]] = Interval(3.0)
    matrix[[1,0]] = Interval(4.0); matrix[[1,1]] = Interval(5.0); matrix[[1,2]] = Interval(6.0)
    matrix[[2,0]] = Interval(7.0); matrix[[2,1]] = Interval(8.0); matrix[[2,2]] = Interval(9.0)

    # transposeInPlace() is not exposed: reassigning from transpose() is
    # the equivalent (Python has no in-place mutation of the storage
    # layout, but nothing else references this matrix at this point).
    matrix = matrix.transpose()

    self.check_point(matrix[(0,1)], 4.0)
    self.check_point(matrix[(1,0)], 2.0)
    self.check_point(matrix[(2,1)], 6.0)

  def test_diagonal_extraction(self):

    matrix = AffineMatrix(3, 3)
    matrix[[0,0]] = Interval(1.0); matrix[[0,1]] = Interval(2.0); matrix[[0,2]] = Interval(3.0)
    matrix[[1,0]] = Interval(4.0); matrix[[1,1]] = Interval(5.0); matrix[[1,2]] = Interval(6.0)
    matrix[[2,0]] = Interval(7.0); matrix[[2,1]] = Interval(8.0); matrix[[2,2]] = Interval(9.0)

    diagonal = affine_diagonal(matrix)
    self.assertTrue(diagonal.size() == 3)
    self.check_point(diagonal[0], 1.0)
    self.check_point(diagonal[1], 5.0)
    self.check_point(diagonal[2], 9.0)

  def test_diagonal_assignment(self):

    matrix = AffineMatrix(3, 3)
    matrix.init(Affine(0.0))

    diagonal_variables = AffineVariables(IntervalVector([[2],[4],[6]]))
    diagonal = AffineVector(diagonal_variables)
    affine_set_diagonal(matrix, diagonal)

    self.check_point(matrix[(0,0)], 2.0)
    self.check_point(matrix[(1,1)], 4.0)
    self.check_point(matrix[(2,2)], 6.0)
    self.check_point(matrix[(0,1)], 0.0)

  def test_setzero_setones_setconstant(self):

    matrix = AffineMatrix(2, 3)

    # setZero()/setOnes()/setConstant() are not exposed: init(x)
    # broadcasts a single Affine value to every entry, which is exactly
    # their common implementation.
    matrix.init(Affine(0.0))
    for i in range(matrix.rows()):
      for j in range(matrix.cols()):
        self.check_point(matrix[(i,j)], 0.0)

    matrix.init(Affine(1.0))
    for i in range(matrix.rows()):
      for j in range(matrix.cols()):
        self.check_point(matrix[(i,j)], 1.0)

    matrix.init(Affine(3.5))
    for i in range(matrix.rows()):
      for j in range(matrix.cols()):
        self.check_point(matrix[(i,j)], 3.5)

  def test_zero_ones_constant_factories(self):

    zero = AffineMatrix.zero(2, 2)
    ones = AffineMatrix.ones(2, 2)
    constant = AffineMatrix.constant(2, 2, Affine(-3.0))

    for i in range(2):
      for j in range(2):
        self.check_point(zero[(i,j)], 0.0)
        self.check_point(ones[(i,j)], 1.0)
        self.check_point(constant[(i,j)], -3.0)

  def test_identity_factory(self):

    identity = AffineMatrix.eye(3, 3)

    for i in range(identity.rows()):
      for j in range(identity.cols()):
        self.check_point(identity[(i,j)], 1.0 if i == j else 0.0)

  def test_resize_and_conservativeresize(self):

    matrix = AffineMatrix(2, 2)
    matrix[[0,0]] = Interval(1.0); matrix[[0,1]] = Interval(2.0)
    matrix[[1,0]] = Interval(3.0); matrix[[1,1]] = Interval(4.0)

    matrix.conservativeResize(3, 3)
    self.assertTrue(matrix.rows() == 3)
    self.assertTrue(matrix.cols() == 3)
    self.check_point(matrix[(0,0)], 1.0)
    self.check_point(matrix[(0,1)], 2.0)
    self.check_point(matrix[(1,0)], 3.0)
    self.check_point(matrix[(1,1)], 4.0)

    matrix.resize(4, 2)
    self.assertTrue(matrix.rows() == 4)
    self.assertTrue(matrix.cols() == 2)

  def test_arithmetic_expressions(self):

    lhs = AffineMatrix(2, 2); rhs = AffineMatrix(2, 2)
    lhs[[0,0]] = Interval(1.0); lhs[[0,1]] = Interval(2.0)
    lhs[[1,0]] = Interval(3.0); lhs[[1,1]] = Interval(4.0)
    rhs[[0,0]] = Interval(5.0); rhs[[0,1]] = Interval(6.0)
    rhs[[1,0]] = Interval(7.0); rhs[[1,1]] = Interval(8.0)

    total = lhs + rhs
    difference = rhs - lhs
    expression = 2.0 * lhs - rhs

    self.check_point(total[(0,0)], 6.0)
    self.check_point(total[(1,1)], 12.0)
    self.check_point(difference[(0,1)], 4.0)
    self.check_point(difference[(1,0)], 4.0)
    self.check_point(expression[(0,0)], -3.0)
    self.check_point(expression[(1,1)], 0.0)

  def test_scalar_compound_assignments_intervalmatrix(self):

    value = IntervalMatrix([
      [Interval(1.0), Interval(-2.0), Interval(3.0)],
      [Interval(4.0), Interval(5.0), Interval(-6.0)]
    ])

    # array() is not exposed: += and -= only accept another matrix (see
    # the module docstring), so a same-size constant matrix stands in for
    # the scalar broadcast; *= and /= already accept a plain scalar.
    value += IntervalMatrix.constant(2, 3, Interval(2.0))
    self.assertTrue(value[(0,0)] == Interval(3.0))
    self.assertTrue(value[(0,1)] == Interval(0.0))
    self.assertTrue(value[(1,2)] == Interval(-4.0))

    value -= IntervalMatrix.constant(2, 3, Interval(2.0))
    self.assertTrue(value[(0,0)] == Interval(1.0))
    self.assertTrue(value[(0,1)] == Interval(-2.0))
    self.assertTrue(value[(1,2)] == Interval(-6.0))

    value *= 2.0
    self.assertTrue(value[(0,0)] == Interval(2.0))
    self.assertTrue(value[(0,1)] == Interval(-4.0))
    self.assertTrue(value[(1,2)] == Interval(-12.0))

    value /= 2.0
    self.assertTrue(value[(0,0)] == Interval(1.0))
    self.assertTrue(value[(0,1)] == Interval(-2.0))
    self.assertTrue(value[(1,2)] == Interval(-6.0))

  def test_scalar_compound_assignments_affinematrix(self):

    value = AffineMatrix(make_interval_matrix_2x3())
    value.init(Affine(0.0)) # start from a clean point matrix, see below
    value = AffineMatrix(IntervalMatrix([
      [Interval(1.0), Interval(-2.0), Interval(3.0)],
      [Interval(4.0), Interval(5.0), Interval(-6.0)]
    ]))

    value += AffineMatrix.constant(2, 3, Affine(2.0))
    self.check_point(value[(0,0)], 3.0)
    self.check_point(value[(0,1)], 0.0)
    self.check_point(value[(1,2)], -4.0)

    value -= AffineMatrix.constant(2, 3, Affine(2.0))
    self.check_point(value[(0,0)], 1.0)
    self.check_point(value[(0,1)], -2.0)
    self.check_point(value[(1,2)], -6.0)

    value *= 2.0
    self.check_point(value[(0,0)], 2.0)
    self.check_point(value[(0,1)], -4.0)
    self.check_point(value[(1,2)], -12.0)

    value /= 2.0
    self.check_point(value[(0,0)], 1.0)
    self.check_point(value[(0,1)], -2.0)
    self.check_point(value[(1,2)], -6.0)

  def test_matrix_product_2(self):

    lhs = AffineMatrix(2, 3); rhs = AffineMatrix(3, 2)
    lhs[[0,0]] = Interval(1.0); lhs[[0,1]] = Interval(2.0); lhs[[0,2]] = Interval(3.0)
    lhs[[1,0]] = Interval(4.0); lhs[[1,1]] = Interval(5.0); lhs[[1,2]] = Interval(6.0)
    rhs[[0,0]] = Interval(7.0);  rhs[[0,1]] = Interval(8.0)
    rhs[[1,0]] = Interval(9.0);  rhs[[1,1]] = Interval(10.0)
    rhs[[2,0]] = Interval(11.0); rhs[[2,1]] = Interval(12.0)

    result = lhs @ rhs
    self.check_point(result[(0,0)], 58.0)
    self.check_point(result[(0,1)], 64.0)
    self.check_point(result[(1,0)], 139.0)
    self.check_point(result[(1,1)], 154.0)

  def test_matrix_vector_product(self):

    matrix = AffineMatrix(2, 3)
    matrix[[0,0]] = Interval(1.0); matrix[[0,1]] = Interval(2.0); matrix[[0,2]] = Interval(3.0)
    matrix[[1,0]] = Interval(4.0); matrix[[1,1]] = Interval(5.0); matrix[[1,2]] = Interval(6.0)

    vector = AffineVector(AffineVariables(IntervalVector([[2],[-1],[3]])))

    result = matrix @ vector
    self.check_point(result[0], 9.0)
    self.check_point(result[1], 21.0)

  def test_real_matrix_times_affinevariables(self):

    matrix = IntervalMatrix([
      [Interval(1.0), Interval(2.0), Interval(3.0)],
      [Interval(4.0), Interval(5.0), Interval(6.0)]
    ])

    variables = AffineVariables(Vector([2.0, -1.0, 3.0]))
    result = matrix @ variables

    self.check_point(result[0], 9.0)
    self.check_point(result[1], 21.0)

  def test_affinevariables_outer_product(self):

    variables = AffineVariables(Vector([2.0, -1.0, 3.0]))
    row = IntervalRow([Interval(4.0), Interval(-2.0)])

    result = affine_outer_product(AffineVector(variables), row)

    self.assertTrue(result.rows() == 3)
    self.assertTrue(result.cols() == 2)
    self.check_point(result[(0,0)], 8.0)
    self.check_point(result[(0,1)], -4.0)
    self.check_point(result[(2,0)], 12.0)
    self.check_point(result[(2,1)], -6.0)

  def test_trace(self):

    matrix = AffineMatrix(3, 3)
    matrix[[0,0]] = Interval(1.0); matrix[[0,1]] = Interval(2.0); matrix[[0,2]] = Interval(3.0)
    matrix[[1,0]] = Interval(4.0); matrix[[1,1]] = Interval(5.0); matrix[[1,2]] = Interval(6.0)
    matrix[[2,0]] = Interval(7.0); matrix[[2,1]] = Interval(8.0); matrix[[2,2]] = Interval(9.0)

    result = affine_trace(matrix)
    self.check_point(result, 15.0)

  def test_replicate(self):

    row_variables = AffineVariables(IntervalVector([[1],[2]]))
    row = AffineRow([row_variables[0], row_variables[1]])

    replicated = affine_replicate(row, 3, 2)
    self.assertTrue(replicated.rows() == 3)
    self.assertTrue(replicated.cols() == 4)

    self.check_point(replicated[(0,0)], 1.0)
    self.check_point(replicated[(0,3)], 2.0)
    self.check_point(replicated[(2,0)], 1.0)
    self.check_point(replicated[(2,3)], 2.0)

  def test_interval_valued_block_operations(self):

    variables = AffineVariables(IntervalVector([[-2.0, 1.0], [3.0, 5.0], [-1.0, 4.0]]))

    matrix = AffineMatrix(3, 2)
    matrix.set_col(0, variables)
    matrix.set_col(1, -AffineVector(variables))

    block = matrix.block(0, 0, 2, 2)

    self.assertTrue(block[(0,0)].itv() == variables[0].itv())
    self.assertTrue(block[(1,0)].itv() == variables[1].itv())
    self.assertTrue(block[(0,1)].itv() == (-variables[0]).itv())
    self.assertTrue(block[(1,1)].itv() == (-variables[1]).itv())

  def test_compound_assignments_with_interval_scalar(self):

    initial = AffineMatrix(2, 2)
    initial[[0,0]] = Interval(1.0); initial[[0,1]] = Interval(2.0)
    initial[[1,0]] = Interval(3.0); initial[[1,1]] = Interval(4.0)

    value = AffineMatrix(initial)
    value += Interval(2.0)
    self.check_point(value[(0,0)], 3.0)
    self.check_point(value[(0,1)], 4.0)
    self.check_point(value[(1,0)], 5.0)
    self.check_point(value[(1,1)], 6.0)

    value = AffineMatrix(initial)
    value -= Interval(2.0)
    self.check_point(value[(0,0)], -1.0)
    self.check_point(value[(0,1)], 0.0)
    self.check_point(value[(1,0)], 1.0)
    self.check_point(value[(1,1)], 2.0)

    value = AffineMatrix(initial)
    value *= Interval(2.0)
    self.check_point(value[(0,0)], 2.0)
    self.check_point(value[(0,1)], 4.0)
    self.check_point(value[(1,0)], 6.0)
    self.check_point(value[(1,1)], 8.0)

    value = AffineMatrix(initial)
    value /= Interval(2.0)
    self.check_point(value[(0,0)], 0.5)
    self.check_point(value[(0,1)], 1.0)
    self.check_point(value[(1,0)], 1.5)
    self.check_point(value[(1,1)], 2.0)

  def test_compound_assignments_with_affine_scalar(self):

    initial = AffineMatrix(2, 2)
    initial[[0,0]] = Interval(1.0); initial[[0,1]] = Interval(2.0)
    initial[[1,0]] = Interval(3.0); initial[[1,1]] = Interval(4.0)

    scalar = Affine(Interval(2.0))

    value = AffineMatrix(initial)
    value += scalar
    self.check_point(value[(0,0)], 3.0)
    self.check_point(value[(0,1)], 4.0)
    self.check_point(value[(1,0)], 5.0)
    self.check_point(value[(1,1)], 6.0)

    value = AffineMatrix(initial)
    value -= scalar
    self.check_point(value[(0,0)], -1.0)
    self.check_point(value[(0,1)], 0.0)
    self.check_point(value[(1,0)], 1.0)
    self.check_point(value[(1,1)], 2.0)

    value = AffineMatrix(initial)
    value *= scalar
    self.check_point(value[(0,0)], 2.0)
    self.check_point(value[(0,1)], 4.0)
    self.check_point(value[(1,0)], 6.0)
    self.check_point(value[(1,1)], 8.0)

    value = AffineMatrix(initial)
    value /= scalar
    self.check_point(value[(0,0)], 0.5)
    self.check_point(value[(0,1)], 1.0)
    self.check_point(value[(1,0)], 1.5)
    self.check_point(value[(1,1)], 2.0)

  def test_compound_assignments_with_intervalmatrix(self):

    initial = AffineMatrix(2, 2)
    initial[[0,0]] = Interval(1.0); initial[[0,1]] = Interval(2.0)
    initial[[1,0]] = Interval(3.0); initial[[1,1]] = Interval(4.0)

    other = IntervalMatrix(2, 2)
    other[[0,0]] = Interval(2.0); other[[0,1]] = Interval(3.0)
    other[[1,0]] = Interval(4.0); other[[1,1]] = Interval(5.0)

    value = AffineMatrix(initial)
    value += other
    self.check_point(value[(0,0)], 3.0)
    self.check_point(value[(0,1)], 5.0)
    self.check_point(value[(1,0)], 7.0)
    self.check_point(value[(1,1)], 9.0)

    value = AffineMatrix(initial)
    value -= other
    self.check_point(value[(0,0)], -1.0)
    self.check_point(value[(0,1)], -1.0)
    self.check_point(value[(1,0)], -1.0)
    self.check_point(value[(1,1)], -1.0)

    # There is no compound @= for a full matrix-matrix product (only the
    # scalar overloads are exposed for +=/-=/*=//=, see the module
    # docstring): reassigning from the non-compound @ is the equivalent.
    value = AffineMatrix(initial)
    value = value @ other
    self.check_point(value[(0,0)], 10.0)
    self.check_point(value[(0,1)], 13.0)
    self.check_point(value[(1,0)], 22.0)
    self.check_point(value[(1,1)], 29.0)

  def test_intervalmatrix_compound_assignments_with_intervalmatrix(self):

    initial = IntervalMatrix(2, 2)
    initial[[0,0]] = Interval(1.0); initial[[0,1]] = Interval(2.0)
    initial[[1,0]] = Interval(3.0); initial[[1,1]] = Interval(4.0)

    other = IntervalMatrix(2, 2)
    other[[0,0]] = Interval(2.0); other[[0,1]] = Interval(3.0)
    other[[1,0]] = Interval(4.0); other[[1,1]] = Interval(5.0)

    value = IntervalMatrix(initial)
    value += other
    self.assertTrue(value[(0,0)]==Interval(3.0))
    self.assertTrue(value[(0,1)]==Interval(5.0))
    self.assertTrue(value[(1,0)]==Interval(7.0))
    self.assertTrue(value[(1,1)]==Interval(9.0))

    value = IntervalMatrix(initial)
    value -= other
    self.assertTrue(value[(0,0)]==Interval(-1.0))
    self.assertTrue(value[(0,1)]==Interval(-1.0))
    self.assertTrue(value[(1,0)]==Interval(-1.0))
    self.assertTrue(value[(1,1)]==Interval(-1.0))

    value = IntervalMatrix(initial)
    value = value * other
    self.assertTrue(value[(0,0)]==Interval(10.0))
    self.assertTrue(value[(0,1)]==Interval(13.0))
    self.assertTrue(value[(1,0)]==Interval(22.0))
    self.assertTrue(value[(1,1)]==Interval(29.0))

  def test_compound_assignments_with_affinematrix(self):

    initial = AffineMatrix(2, 2)
    initial[[0,0]] = Interval(1.0); initial[[0,1]] = Interval(2.0)
    initial[[1,0]] = Interval(3.0); initial[[1,1]] = Interval(4.0)

    rhs = AffineMatrix(2, 2)
    rhs[[0,0]] = Interval(2.0); rhs[[0,1]] = Interval(3.0)
    rhs[[1,0]] = Interval(4.0); rhs[[1,1]] = Interval(5.0)

    value = AffineMatrix(initial)
    value += rhs
    self.check_point(value[(0,0)], 3.0)
    self.check_point(value[(0,1)], 5.0)
    self.check_point(value[(1,0)], 7.0)
    self.check_point(value[(1,1)], 9.0)

    value = AffineMatrix(initial)
    value -= rhs
    self.check_point(value[(0,0)], -1.0)
    self.check_point(value[(0,1)], -1.0)
    self.check_point(value[(1,0)], -1.0)
    self.check_point(value[(1,1)], -1.0)

    value = AffineMatrix(initial)
    value = value @ rhs
    self.check_point(value[(0,0)], 10.0)
    self.check_point(value[(0,1)], 13.0)
    self.check_point(value[(1,0)], 22.0)
    self.check_point(value[(1,1)], 29.0)

  def test_real_matrix_times_affinevector(self):

    matrix = make_real_matrix_2x3()
    vector = make_affine_vector_3()

    result = matrix @ vector
    expected = multiply_intervals(IntervalMatrix(matrix), interval_hull_vector(vector))

    self.check_encloses_vector(result, expected)

  def test_real_matrix_times_affinevariables_2(self):

    matrix = make_real_matrix_2x3()
    vector = make_variable_vector_3()

    result = matrix @ vector
    self.assertTrue(result.size() == 2)

    expected = multiply_intervals(IntervalMatrix(matrix), interval_hull_vector(vector))
    self.check_encloses_vector(AffineVector(result), expected)

  def test_interval_matrix_times_affinevector(self):

    matrix = make_interval_matrix_2x3()
    vector = make_affine_vector_3()

    result = matrix @ vector
    expected = multiply_intervals(matrix, interval_hull_vector(vector))

    self.check_encloses_vector(result, expected)

  def test_interval_matrix_times_affinevariables(self):

    matrix = make_interval_matrix_2x3()
    vector = make_variable_vector_3()

    result = matrix @ vector
    self.assertTrue(result.size() == 2)

    expected = multiply_intervals(matrix, interval_hull_vector(vector))
    self.check_encloses_vector(AffineVector(result), expected)

  def test_affine_matrix_times_affinevector(self):

    matrix = make_affine_matrix_2x3()
    vector = make_affine_vector_3()

    result = matrix @ vector
    expected = multiply_intervals(interval_hull_matrix(matrix), interval_hull_vector(vector))

    self.check_encloses_vector(result, expected)

  def test_affine_matrix_times_affinevariables(self):

    matrix = make_affine_matrix_2x3()
    vector = make_variable_vector_3()

    result = matrix @ vector
    self.assertTrue(result.size() == 2)

    expected = multiply_intervals(interval_hull_matrix(matrix), interval_hull_vector(vector))
    self.check_encloses_vector(AffineVector(result), expected)

  def test_real_matrix_expression_times_affinevariables(self):

    # Eigen lazy expression templates ((a+b), 3.0*a, .block(...) left
    # unevaluated) have no Python equivalent: every Python expression is
    # already a concrete value, so the three C++ SECTIONs collapse to
    # plain direct computations here.
    a = make_real_matrix_2x3()
    b = 2.0 * make_real_matrix_2x3()
    vector = make_variable_vector_3()

    # sum expression
    result = (a + b) @ vector
    expected = multiply_intervals(IntervalMatrix(a + b), interval_hull_vector(vector))
    self.check_encloses_vector(AffineVector(result), expected)

    # scaled expression
    result = (3.0 * a) @ vector
    expected = multiply_intervals(IntervalMatrix(3.0 * a), interval_hull_vector(vector))
    self.check_encloses_vector(AffineVector(result), expected)

    # block expression
    large = Matrix([[1.0,-2.0,3.0,8.0],[4.0,0.5,-1.0,9.0],[7.0,6.0,5.0,4.0]])
    block = large.block(0, 0, 2, 3)
    result = block @ vector
    expected = multiply_intervals(IntervalMatrix(block), interval_hull_vector(vector))
    self.check_encloses_vector(AffineVector(result), expected)

  def test_interval_matrix_expression_times_affinevariables(self):

    a = make_interval_matrix_2x3()
    b = make_interval_matrix_2x3()
    vector = make_variable_vector_3()

    result = (a + b) @ vector
    expected = multiply_intervals(a + b, interval_hull_vector(vector))
    self.check_encloses_vector(AffineVector(result), expected)

  def test_affine_matrix_expression_times_affinevariables(self):

    a = make_affine_matrix_2x3()
    b = make_affine_matrix_2x3()
    vector = make_variable_vector_3()

    result = (a + b) @ vector
    evaluated = a + b
    expected = multiply_intervals(interval_hull_matrix(evaluated), interval_hull_vector(vector))
    self.check_encloses_vector(AffineVector(result), expected)

  def test_affinevariables_outer_product_with_real_row(self):

    vector = make_variable_vector_3()
    # Row has no Python constructor (it exists only as a return type, e.g.
    # from Vector.transpose()): a plain list of floats stands in for a
    # "real" row here -- _as_interval()/Affine's scalar multiplication
    # both already accept a plain float, so the helpers below need no Row.
    row = [2.0, -1.0, 0.5, 3.0]

    result = affine_outer_product(AffineVector(vector), row)
    self.assertTrue(result.rows() == 3)
    self.assertTrue(result.cols() == 4)

    expected = interval_outer_product(vector, row)
    self.check_encloses_matrix(result, expected)

  def test_affinevariables_outer_product_with_interval_row(self):

    vector = make_variable_vector_3()
    row = IntervalRow([Interval(1.0,2.0), Interval(-1.0,1.0), Interval(0.5), Interval(2.0,3.0)])

    result = affine_outer_product(AffineVector(vector), row)
    self.assertTrue(result.rows() == 3)
    self.assertTrue(result.cols() == 4)

    expected = interval_outer_product(vector, row)
    self.check_encloses_matrix(result, expected)

  def test_affinevariables_outer_product_with_affinerow(self):

    vector = make_variable_vector_3()
    row = AffineRow([Interval(1.0,2.0), Interval(-1.0,1.0), Interval(0.5), Interval(2.0,3.0)])

    result = affine_outer_product(AffineVector(vector), row)
    self.assertTrue(result.rows() == 3)
    self.assertTrue(result.cols() == 4)

    expected = interval_outer_product(vector, row)
    self.check_encloses_matrix(result, expected)

  def test_affinevector_outer_product_with_real_row(self):

    vector = make_affine_vector_3()
    # See test_affinevariables_outer_product_with_real_row: Row has no
    # Python constructor, so a plain list of floats stands in for it.
    row = [-2.0, 4.0]

    result = affine_outer_product(vector, row)
    expected = interval_outer_product(vector, row)
    self.check_encloses_matrix(result, expected)

  def test_affinevector_outer_product_with_interval_row(self):

    vector = make_affine_vector_3()
    row = IntervalRow([Interval(-2.0,-1.0), Interval(3.0,4.0)])

    result = affine_outer_product(vector, row)
    expected = interval_outer_product(vector, row)
    self.check_encloses_matrix(result, expected)

  def test_affinevector_outer_product_with_affinerow(self):

    vector = make_affine_vector_3()
    row = AffineRow([Interval(-2.0,-1.0), Interval(3.0,4.0)])

    result = affine_outer_product(vector, row)
    expected = interval_outer_product(vector, row)
    self.check_encloses_matrix(result, expected)

  def test_affinerow_times_real_matrix(self):

    row = make_affine_row_3()
    matrix = Matrix([[1.0,-1.0],[2.0,3.0],[-2.0,4.0]])

    # AffineRow @ Matrix returns an AffineRow (only 1-D indexable, see the
    # module docstring): wrap it as a 1 x n AffineMatrix so
    # check_encloses_matrix can use its (i,j) tuple indexing.
    result = AffineMatrix(row @ matrix)
    expected = multiply_intervals(interval_hull_row_matrix(row), IntervalMatrix(matrix))
    self.check_encloses_matrix(result, expected)

  def test_affinerow_times_interval_matrix(self):

    row = make_affine_row_3()
    matrix = IntervalMatrix([
      [Interval(1.0,1.1), Interval(-1.0,-0.9)],
      [Interval(2.0,2.1), Interval(3.0,3.1)],
      [Interval(-2.0,-1.9), Interval(4.0,4.1)]
    ])

    # See test_affinerow_times_real_matrix: wrap the AffineRow result as a
    # 1 x n AffineMatrix for check_encloses_matrix's (i,j) indexing.
    result = AffineMatrix(row @ matrix)
    expected = multiply_intervals(interval_hull_row_matrix(row), matrix)
    self.check_encloses_matrix(result, expected)

  def test_affinerow_times_affine_matrix(self):

    row = make_affine_row_3()
    matrix = AffineMatrix(3, 2)
    matrix[[0,0]] = Interval(1.0); matrix[[0,1]] = Interval(-1.0)
    matrix[[1,0]] = Interval(2.0); matrix[[1,1]] = Interval(3.0)
    matrix[[2,0]] = Interval(-2.0); matrix[[2,1]] = Interval(4.0)

    # See test_affinerow_times_real_matrix: wrap the AffineRow result as a
    # 1 x n AffineMatrix for check_encloses_matrix's (i,j) indexing.
    result = AffineMatrix(row @ matrix)
    expected = multiply_intervals(interval_hull_row_matrix(row), interval_hull_matrix(matrix))
    self.check_encloses_matrix(result, expected)

  def test_real_column_vector_times_affinerow(self):

    vector = Vector([1.0, -2.0, 4.0])
    row = make_affine_row_3()

    result = affine_outer_product(vector, row)
    expected = interval_outer_product(vector, row)
    self.check_encloses_matrix(result, expected)

  def test_matrix_products_preserve_dynamic_dimensions(self):

    lhs = AffineMatrix(4, 3)
    rhs = AffineMatrix(3, 5)

    for i in range(4):
      for j in range(3):
        lhs[[i,j]] = Interval(float(1 + i + 2*j))

    for i in range(3):
      for j in range(5):
        rhs[[i,j]] = Interval(float(2 - i + j))

    result = lhs @ rhs
    self.assertTrue(result.rows() == 4)
    self.assertTrue(result.cols() == 5)

    expected = multiply_intervals(interval_hull_matrix(lhs), interval_hull_matrix(rhs))
    self.check_encloses_matrix(result, expected)

  def test_degenerate_products_produce_degenerate_results(self):

    matrix = Matrix([[1.0,2.0,3.0],[-1.0,4.0,0.5]])
    variables = AffineVariables(Vector([2.0, -1.0, 3.0]))

    result = matrix @ variables
    self.assertTrue(result.size() == 2)
    self.assertTrue(result[0].itv() == Interval(9.0))
    self.assertTrue(result[1].itv() == Interval(-4.5))

  def test_identity_matrix_preserves_hulls(self):

    variables = make_variable_vector_3()
    identity = Matrix.eye(3, 3)

    result = identity @ variables
    self.assertTrue(result.size() == variables.size())

    for i in range(variables.size()):
      self.assertTrue(result[i].itv().is_superset(variables[i].itv()))

  def test_zero_real_matrix_annihilates(self):

    variables = make_variable_vector_3()
    zero = Matrix.zero(2, 3)

    result = zero @ variables
    self.assertTrue(result.size() == 2)
    self.assertTrue(result[0].itv() == Interval(0.0))
    self.assertTrue(result[1].itv() == Interval(0.0))

  def test_one_by_one_products_retain_scalar(self):

    variables = AffineVariables(IntervalVector([[-2.0, 3.0]]))

    # real matrix
    matrix = Matrix([[-2.0]])
    result = matrix @ variables
    self.assertTrue(result[0].itv().is_superset(Interval(-6.0, 4.0)))

    # interval matrix
    matrix = IntervalMatrix([[Interval(-2.0, -1.0)]])
    result = matrix @ variables
    self.assertTrue(result[0].itv().is_superset(matrix[(0,0)] * variables[0].itv()))

    # Affine matrix
    matrix = AffineMatrix(1, 1)
    matrix[[0,0]] = Interval(-2.0, -1.0)
    result = matrix @ variables
    self.assertTrue(result[0].itv().is_superset(matrix[(0,0)].itv() * variables[0].itv()))

  def test_empty_coefficient_propagation_through_products(self):

    matrix = AffineMatrix(1, 2)
    matrix[[0,0]] = Interval(1.0)
    matrix[[0,1]] = Interval.empty()

    variables = AffineVariables(IntervalVector([[1.0, 2.0], [3.0, 4.0]]))
    result = matrix @ variables

    self.assertTrue(result.size() == 1)
    self.assertTrue(result[0].is_empty())

  def test_matrix_product_preserves_dependency(self):

    matrix = make_real_matrix_2x3()
    x = make_variable_vector_3()

    result = (matrix @ x) - (matrix @ x)

    for i in range(result.size()):
      self.assertTrue(result[i].itv() == Interval(0.0))

  def test_stream_output_for_nonempty_matrix(self):

    matrix = make_point_matrix_3x4()
    s = str(matrix)
    self.assertTrue(len(s) != 0)
    self.assertTrue(s.find("empty") == -1)


if __name__ ==  '__main__':
  unittest.main()
