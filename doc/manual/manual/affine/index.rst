.. _sec-affine:

Affine forms
============

Codac provides an affine-arithmetic domain, complementary to the
:ref:`Interval structures <sec-intervals-class>`. An affine form tracks first-order dependencies between quantities that
share common noise symbols. These shared symbols preserve the identity of
uncertain quantities through a computation and can significantly reduce the
*wrapping effect* compared to plain interval arithmetic whenever a
computation reuses the same uncertain quantity several times.

The main affine types have distinct roles:

.. list-table:: Main affine types
   :widths: 30 70
   :header-rows: 1

   * - Type
     - Role
   * - ``Affine``
     - Scalar affine form used for arithmetic and nonlinear computations.
   * - ``AffineVariables``
     - Declares independent uncertain quantities and introduces their noise
       symbols.
   * - ``AffineVector``
     - Working column vector of affine forms used for vector computations.
   * - ``AffineRow``
     - Working row vector of affine forms, typically obtained from an
       ``AffineVector`` transpose or a matrix row.
   * - ``AffineMatrix``
     - Working matrix of affine forms used for matrix and matrix-vector
       computations.

See the dedicated pages:

- :ref:`The Affine class <sec-affine-class>`: scalar affine forms.
- :ref:`The AffineVariables class <sec-affine-variables-class>`: declared
  affine variables and noise-symbol management.
- :ref:`The AffineVector class <sec-affine-vector-class>`: vector
  computations with affine forms.
- :ref:`The AffineMatrix class <sec-affine-matrix-class>`: matrix and
  matrix-vector computations while preserving dependencies.

.. toctree::
   :hidden:

   AffineForm_class.rst
   AffineVariables_class.rst
   AffineVector_class.rst
   AffineMatrix_class.rst

.. note::

   This module is currently C++ only; Python bindings are not yet available.


Conceptual workflow
-------------------

The affine domain separates the **declaration** of independent uncertain
quantities from the **computation** performed with the resulting affine
forms:

.. code-block:: text

   AffineVariables
       x[0] <-> eps_0
       x[1] <-> eps_1
       ...
           |
           | arithmetic / nonlinear functions
           v
   Affine / AffineVector / AffineMatrix

``AffineVariables`` is the only type in this module that introduces new
noise symbols. ``Affine``, ``AffineVector`` and ``AffineMatrix`` are the
working types used to carry those symbols through computations.

The main pages should therefore be read in the following order:

#. :ref:`The Affine class <sec-affine-class>` — the scalar affine form and
   its arithmetic.
#. :ref:`The AffineVariables class <sec-affine-variables-class>` — how
   independent noise symbols are introduced and managed.
#. :ref:`The AffineVector class <sec-affine-vector-class>` — vector
   computations with affine forms.
#. :ref:`The AffineMatrix class <sec-affine-matrix-class>` — matrix and
   matrix-vector computations while preserving dependencies.

.. seealso::

   The affine-arithmetic model and the Chebyshev and MinRange linearization
   techniques are described in J. Ninin's PhD thesis, chapter 2, and in the
   foundational work of J. Stolfi and L. H. de Figueiredo, *An Introduction
   to Affine Arithmetic*, TEMA 4(3), 2003.
