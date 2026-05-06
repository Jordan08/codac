.. _sec-ctc-geom-ctcparallelepiped:

The CtcParallelpiped contractor
===============================

  Main author: `Maël Godard <https://godardma.github.io/>`_

.. doxygenclass:: codac2::CtcParallelepiped
  :project: codac

Methods
-------

.. doxygenfunction:: codac2::CtcParallelepiped::contract(IntervalVector&) const
  :project: codac

.. tabs::

  .. group-tab:: Python

    .. literalinclude:: src.py
      :language: py
      :start-after: [ctcparallelepiped-1-beg]
      :end-before: [ctcparallelepiped-1-end]
      :dedent: 4

  .. group-tab:: C++

    .. literalinclude:: src.cpp
      :language: c++
      :start-after: [ctcparallelepiped-1-beg]
      :end-before: [ctcparallelepiped-1-end]
      :dedent: 4

.. admonition:: Technical documentation

  See the `C++ API documentation of this class <../../api/html/classcodac2_1_1_ctc_parallelepiped.html>`_.