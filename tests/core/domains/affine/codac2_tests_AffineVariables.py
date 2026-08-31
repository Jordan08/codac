#!/usr/bin/env python

# Codac tests
#
# Tests specific to AffineVariables, the Python binding of
# AffineVarMainVector<T>, the "affine variable" container that introduces
# fresh, uniquely-indexed noise symbols: distinct noise symbols per
# component, resize()/conservativeResize() semantics, streaming, and
# transparent conversion to AffineVector.
#
# Note: a few C++ tests from codac2_tests_AffineVariables.cpp are not
# translated here because they exercise capabilities that are deliberately
# not exposed to Python:
#  - index-to-index assignment between two elements of an AffineVariables
#    (AffineVarMain::operator=(const AffineVarMain&), which copies the
#    whole noise representation, including the source's own noise index)
#    has no Python binding: __setitem__ only accepts an Interval, which
#    reaches AffineVarMain::operator=(const Interval&) instead (see the
#    comment in codac2_py_AffineVariables.cpp);
#  - there is no standalone Python class for a single AffineVarMain (only
#    AffineVariables, the vector container).
#
# ----------------------------------------------------------------------------
#  \date       2026
#  \author     Jordan Ninin
#  \copyright  Copyright 2026 Codac Team
#  \license    GNU Lesser General Public License (LGPL)

import unittest
from codac import *


def make_variable_vector_5():
  return AffineVariables(IntervalVector([
    [-2.0,-1.0], [0.0,2.0], [3.0,4.0], [-1.0,1.0], [5.0,7.0]
  ]))

def make_variable_vector_3():
  return AffineVariables(IntervalVector([[1.0,2.0],[-1.0,1.0],[3.0,4.0]]))


class TestAffineVariables(unittest.TestCase):

  def test_itv(self):

    variables = make_variable_vector_5()
    boxes = variables.itv()

    self.assertTrue(len(boxes) == 5)
    self.assertTrue(boxes[0] == Interval(-2.0,-1.0))
    self.assertTrue(boxes[1] == Interval(0.0,2.0))
    self.assertTrue(boxes[2] == Interval(3.0,4.0))
    self.assertTrue(boxes[3] == Interval(-1.0,1.0))
    self.assertTrue(boxes[4] == Interval(5.0,7.0))

  def test_resize_discards_previous_values_unlike_conservativeresize(self):

    box = IntervalVector([[1.0,2.0],[3.0,4.0]])

    preserved = AffineVariables(box)
    preserved.conservativeResize(3)
    self.assertTrue(preserved[0].itv() == Interval(1.0,2.0))
    self.assertTrue(preserved[1].itv() == Interval(3.0,4.0))
    self.assertTrue(preserved[2].itv() == Interval())

    reset = AffineVariables(box)
    reset.resize(3)
    self.assertTrue(reset[0].itv() == Interval())
    self.assertTrue(reset[1].itv() == Interval())
    self.assertTrue(reset[2].itv() == Interval())

    # resize() also reassigns fresh, sequential noise symbols.
    for i in range(reset.size()):
      self.assertTrue(reset.noise_index(i) == i)

  def test_repr_streams_the_interval_box(self):

    x = AffineVariables(IntervalVector([[1.0,2.0],[3.0,4.0]]))
    s = str(x)
    self.assertTrue(len(s) != 0)
    self.assertTrue(s.find("empty") == -1)

    e = AffineVariables(1)
    e[0] = Interval.empty()
    self.assertTrue(str(e) == "[ empty 1d box ]")

  def test_init_preserves_distinct_noise_symbols(self):

    v = AffineVariables(3)
    v.init(Interval(1.0, 2.0))

    self.assertTrue(v.size() == 3)
    rad = Interval(1.0, 2.0).rad() # 0.5

    for i in range(3):
      # Every component carries the same interval enclosure...
      self.assertTrue(v[i].itv() == Interval(1.0, 2.0))
      self.assertTrue(v[i].rad() == rad)
      self.assertTrue(v[i].noise_count() == 3)
      # ...but each keeps its own, distinct noise symbol: component i has
      # a coefficient equal to the interval's radius on eps_i, and zero
      # on every other eps_j.
      for j in range(3):
        self.assertTrue(v[i].noise(j) == (rad if i == j else 0.0))

    v.conservativeResize(5)

    for i in range(3):
      self.assertTrue(v[i].itv() == Interval(1.0, 2.0))
      self.assertTrue(v[i].rad() == rad)
      self.assertTrue(v[i].noise_count() == 5)
      for j in range(5):
        self.assertTrue(v[i].noise(j) == (rad if i == j else 0.0))

    for i in range(3, 5):
      self.assertTrue(v[i].itv() == Interval())
      self.assertTrue(v[i].noise_count() == 5)

  def test_segment_conversion(self):

    variables = make_variable_vector_5()
    # AffineVariables has no segment(): going through AffineVector, whose
    # subvector(start,end) uses the same inclusive-bounds convention,
    # gives the same three components (indices 1, 2, 3).
    segment = AffineVector(variables).subvector(1, 3)

    self.assertTrue(segment.size() == 3)
    self.assertTrue(segment[0].itv() == variables[1].itv())
    self.assertTrue(segment[1].itv() == variables[2].itv())
    self.assertTrue(segment[2].itv() == variables[3].itv())

  def test_converts_transparently_to_affinevector(self):

    variables = make_variable_vector_3()

    # construction
    values = AffineVector(variables)
    self.assertTrue(values.size() == variables.size())
    for i in range(values.size()):
      self.assertTrue(values[i].itv() == variables[i].itv())

    # assignment: Python has no separate assignment operator distinct from
    # construction, so re-building from the same implicit conversion is
    # the faithful equivalent.
    values = AffineVector(variables)
    self.assertTrue(values.size() == variables.size())
    for i in range(values.size()):
      self.assertTrue(values[i].itv() == variables[i].itv())

  def test_construction_and_assignment_from_intervalvector(self):

    values = IntervalVector([[-1.0,1.0],[2.0,3.0],[4.0,5.0]])

    # Python has no operator=(const IntervalVector&): re-building from the
    # same values exercises the same resize-and-fresh-noise-symbols path.
    x = AffineVariables(values)

    self.assertTrue(x.size() == 3)
    for i in range(x.size()):
      self.assertTrue(x[i].noise_count() == 3)
      self.assertTrue(x[i].itv() == values[i])

  def test_does_not_inherit_eigen_constructors(self):

    # The C++ test parametrizes this over AF_fAF2 specifically (as opposed
    # to AF_Default elsewhere in this file): only one representation is
    # exported to Python, but the sequential-noise-symbol behavior checked
    # here does not depend on that choice.

    x = AffineVariables(3)
    self.assertTrue(x.noise_index(0) == 0)
    self.assertTrue(x.noise_index(1) == 1)
    self.assertTrue(x.noise_index(2) == 2)

    y = AffineVariables(IntervalVector([[-1.,1.],[2.,3.],[4.,5.]]))
    self.assertTrue(y.noise_index(0) == 0)
    self.assertTrue(y.noise_index(1) == 1)
    self.assertTrue(y.noise_index(2) == 2)

    z = AffineVariables(Vector([-1., 2., 3.]))
    self.assertTrue(z.noise_index(0) == 0)
    self.assertTrue(z.noise_index(1) == 1)
    self.assertTrue(z.noise_index(2) == 2)

    w = AffineVariables(0)
    w.resize(3)
    self.assertTrue(w.noise_index(0) == 0)
    self.assertTrue(w.noise_index(1) == 1)
    self.assertTrue(w.noise_index(2) == 2)

  def test_init_places_radius_on_own_noise_symbol(self):

    variables = AffineVariables(3)
    # Give every component a distinct, non-trivial value first, so that
    # the re-init below is genuinely observed to reset the state rather
    # than leaving stale coefficients around.
    variables[0] = Interval(10.0, 20.0)
    variables[1] = Interval(-5.0, 5.0)
    variables[2] = Interval(100.0, 200.0)

    # __getitem__ returns a live reference to the element (matching C++'s
    # own AffineVarMain& here), so mutating it in place reaches the real
    # component directly, exactly like the C++ variables[1].init(x).
    variables[1].init(Interval(3.0, 7.0))

    self.assertTrue(variables[1].noise_count() == 3)
    self.assertTrue(variables.noise_index(1) == 1)
    self.assertTrue(variables[1].itv() == Interval(3.0, 7.0))
    self.assertTrue(variables[1].mid() == 5.0)

    # Unlike a plain Affine (where the whole radius goes to the remainder
    # error term), an affine variable carries it on its own dedicated
    # noise symbol: every other component's contribution to variables[1],
    # and its remainder error term, stay at zero.
    self.assertTrue(variables[1].noise(1) == 2.0)
    self.assertTrue(variables[1].noise(0) == 0.0)
    self.assertTrue(variables[1].noise(2) == 0.0)
    self.assertTrue(variables[1].err() == 0.0)

    # The other components are untouched.
    self.assertTrue(variables[0].itv() == Interval(10.0, 20.0))
    self.assertTrue(variables[2].itv() == Interval(100.0, 200.0))

    # init() with no argument resets to an unbounded placeholder, still
    # tied to the same noise index.
    variables[1].init()
    self.assertTrue(variables[1].is_unbounded())
    self.assertTrue(variables.noise_index(1) == 1)

    # A degenerate interval carries no uncertainty at all, on the noise
    # symbol or anywhere else.
    variables[2].init(Interval(50.0))
    self.assertTrue(variables[2].itv() == Interval(50.0))
    self.assertTrue(variables[2].is_degenerated())
    self.assertTrue(variables[2].noise(2) == 0.0)
    self.assertTrue(variables[2].err() == 0.0)


if __name__ ==  '__main__':
  unittest.main()
