.. _sec-graphics-functions:

2D drawing functions
====================

  Main author: `Maël Godard <https://godardma.github.io>`_

This page lists all the 2D drawing functions available in Codac. Note that the functions described here on a
:ref:`Figure2D <subsec-graphics-2d-figures-figure2d>` object also work on the :ref:`DefaultFigure <subsec-graphics-2d-figures-defaultfigure>`.

For example : 

.. tabs::

  .. code-tab:: py
    
    # On a figure
    fig.draw_box([[2.2,2.5],[2.2,2.5]])

    # On the DefaultFigure
    DefaultFigure.draw_box([[2.2,2.5],[2.2,2.5]])

  .. code-tab:: c++

    // On a figure
    fig.draw_box({{2.2,2.5},{2.2,2.5}});

    // On the DefaultFigure
    DefaultFigure::draw_box({{2.2,2.5},{2.2,2.5}});

  .. code-tab:: matlab
    
    % On a figure
    fig.draw_box(IntervalVector({{2.2,2.5},{2.2,2.5}}));

    % On the DefaultFigure
    DefaultFigure().draw_box(IntervalVector({{2.2,2.5},{2.2,2.5}}));

In addition to the arguments described here for every function, an optionnal argument of type StyleProperties can be added to choose the edge
and fill color (by default black edge no fill). This object can be deduced from one or two Color objects, see :ref:`this page <sec-graphics-colors>` 
for more information.

.. tabs::

  .. code-tab:: py
    
    fig.draw_box([[2.2,2.5],[2.2,2.5]],Color.red()) # Red edge
    fig.draw_box([[2.2,2.5],[2.2,2.5]],[Color.red()]) # Red edge
    fig.draw_box([[2.2,2.5],[2.2,2.5]],[Color.red(),Color.blue()]) # Red edge and blue fill

  .. code-tab:: c++

    fig.draw_box({{2.2,2.5},{2.2,2.5}},Color.red()); // Red edge
    fig.draw_box({{2.2,2.5},{2.2,2.5}},{Color.red()}); // Red edge
    fig.draw_box({{2.2,2.5},{2.2,2.5}},{Color.red(),Color.blue()}); // Red edge and blue fill

  .. code-tab:: matlab

    fig.draw_box(IntervalVector({{2.2,2.5},{2.2,2.5}}),Color().red()); % Red edge
    fig.draw_box(IntervalVector({{2.2,2.5},{2.2,2.5}}),StyleProperties({Color().red()})); % Red edge
    fig.draw_box(IntervalVector({{2.2,2.5},{2.2,2.5}}),StyleProperties({Color().red(),Color().blue()})); % Red edge and blue fill


A line style and a layer can be added to the StyleProperties object, see :ref:`the dedicated page <sec-graphics-colors>` for more information.


Overview of the functions
-------------------------

Below are the detailled available drawing functions. The shapes that can be drawn are:

Geometric shapes
  - Point
  - Box
  - Circle
  - Ring
  - Line
  - Arrow
  - Polyline
  - Polygone
  - Parallelepiped
  - Zonotope
  - Pie
  - Ellipse
  - Ellipsoid

Trajectories
  - Sampled Trajectory
  - Analytic Trajectory

Vehicles
  - Tank
  - AUV
  - Motor Boat

Paving
  - PavingOut (Paving with contractors)
  - PavingInOut (Paving with separators)
  - Subpaving

Geometric shapes
----------------

All the drawable geometric objects can take a last optionnal argument to set up their stroke color, fill color, line style etc.
For further details, refer to :ref:`subsec-graphics-colors-styleproperties`. 

The geometric shapes that can be drawn and their arguments are listed below :

- draw_point

  - Vector : 2D coordinate of the point to draw

- draw_box

  - IntervalVector : the 2D box to draw

- draw_circle

  - Vector : the center of the circle
  - double : the radius of the circle
  
- draw_ring

  - Vector : the center of the ring
  - Interval : the ring range : inner and outer radius

- draw_pie

  - Vector : the center of the pie
  - Interval : the pie range : inner and outer radius
  - Interval : the pie angle : begin and end angle (in radian)

- draw_line

  - Vector : the first point of the line
  - vector : the second point of the line

- draw_line

  - Segment : the segment to draw

- draw_arrow

  - Vector : the first point of the line
  - vector : the second point of the line
  - float : the length of the tip of the arrow

- draw_polyline

  - vector<Vector> : vector where each element is a point of the polyline to draw

- draw_polyline

  - vector<Vector> : vector where each element is a point of the polyline to draw
  - float : the length of the tip of the arrow

- draw_polygone

  - Polygon : the polygon to draw

- draw_parallelepiped

  - Parallelepiped : the parallelepiped to draw

- draw_zonotope

  - Zonotope : the zonotope to draw

- draw_ellipse

  - Vector : center of the ellipse
  - Vector : Half-lengths of the ellipse
  - double : rotation angle of the ellipse (in radian)

- draw_ellipsoid

  - Ellipsoid : the ellipsoid to draw

Trajectories and tubes
----------------------


All the drawable trajectories and tubes can take a last optionnal argument to set up their stroke color, fill color, line style etc.
For further details, refer to :ref:`subsec-graphics-colors-styleproperties` for constant fill and edge color, 
and :ref:`subsec-graphics-colors-stylegradientproperties` for the use of a colormap.

The trajectories and tubes that can be drawn are listed below:

- plot_trajectory

  - SampledTraj | AnalyticTraj : the 1D trajectory to draw (x-axis is the time)

- plot_trajectories

  - SampledTraj | AnalyticTraj : the nD trajectory to draw component by component (x-axis is the time)

- plot_tube

  - SlicedTube : the 1D tube to draw (x-axis is the time)

- plot_tube

  - SlicedTube : the 1D tube to draw (x-axis is the time)
  - SlicedTube : the derivative of the tube to draw

- draw_trajectory

  - SampledTraj | AnalyticTraj : the 2D trajectory to draw

- draw_tube

  - SlicedTube : the 2D tube to draw

Vehicles
--------

All the drawable vehicles can take a last optionnal argument to set up their stroke color, fill color, line style etc.
For further details, refer to :ref:`subsec-graphics-colors-styleproperties`. 

The vehicles that can be drawn and their arguments are listed below :

- draw_tank

  - Vector : the 3D state of the tank to draw (x,y,theta in radian)
  - float : the size of the tank

- draw_AUV

  - Vector : the 3D state of the AUV to draw (x,y,theta in radian)
  - float : the size of the AUV

- draw_motor_boat

  - Vector : the 3D state of the motor boat to draw (x,y,theta in radian)
  - float : the size of the motor boat

Paving
------

All the drawable pavings can take a last optionnal argument to set up their style for the insider, boundary and outside boxes.
For further details, refer to :ref:`subsec-graphics-colors-styleproperties`. 


Miscellaneous
-------------

.. doxygenfunction:: codac2::Figure2D::draw_text(const std::string&, const Vector&, double, const StyleProperties&)
  :project: codac

.. doxygenfunction:: codac2::Figure2D::draw_raster(const std::string&, const IntervalVector&, const StyleProperties&)
  :project: codac

Paving
------

.. doxygenfunction:: codac2::Figure2D::draw_paving(const PavingOut&, const PavingStyle&)
  :project: codac

.. doxygenfunction:: codac2::Figure2D::draw_paving(const PavingOut&, const std::function<void(Figure2D&,const IntervalVector&,const StyleProperties&)>&, const PavingStyle&)
  :project: codac

.. doxygenfunction:: codac2::Figure2D::draw_paving(const PavingInOut&, const PavingStyle&)
  :project: codac

.. doxygenfunction:: codac2::Figure2D::draw_paving(const PavingInOut&, const std::function<void(Figure2D&,const IntervalVector&,const StyleProperties&)>&, const PavingStyle&)
  :project: codac

.. doxygenfunction:: codac2::Figure2D::draw_subpaving(const Subpaving<P>&, const StyleProperties&)
  :project: codac
