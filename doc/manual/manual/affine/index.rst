.. _sec-affine:

Affine forms
============

Codac provides an affine-arithmetic domain, complementary to the
:ref:`Interval structures <sec-intervals-class>`. An affine form tracks
first-order linear correlations between quantities that share common noise
symbols, which can significantly reduce the *wrapping effect* compared to
plain interval arithmetic whenever a computation reuses the same uncertain
quantity several times.

- :ref:`The Affine class <sec-affine-class>`: represents a single affine
  form, and the closely related ``AffineVariables`` type used to introduce
  new, independent noise symbols.
- :ref:`The AffineVariables class <sec-affine-variables-class>`: a vector
  of *declared*, independent affine variables — the only way to introduce
  new noise symbols.
- :ref:`The AffineVector class <sec-affine-vector-class>`: a general-purpose
  column vector of affine forms, used to compute with declared variables.
- :ref:`The AffineMatrix class <sec-affine-matrix-class>`: a general-purpose
  matrix of affine forms, notably for products with a vector of affine
  variables.

.. toctree::
   :hidden:

   AffineForm_class.rst
   AffineVariables_class.rst
   AffineVector_class.rst
   AffineMatrix_class.rst

.. note::

   This module is currently C++ only; Python bindings are not yet available.
