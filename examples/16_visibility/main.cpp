#include <codac>
using namespace codac2;

int main()
{
  // Observation Point and Obstacle Segment
  Vector a({1, 1});
  Segment s({{1,1}, {4,4}}, {{3,3}, {2,2}});

  // Set up the figure
  DefaultFigure::set_axes(axis(0,{-1,6}), axis(1,{-1,6}));
  
  // Show the observation point and the segment
  DefaultFigure::draw_circle(a, 0.05, StyleProperties({Color::dark_green(), Color::green()}, "w:0.025", "z:5"));
  DefaultFigure::draw_line(s[0].mid(), s[1].mid(), StyleProperties(Color::red(), "w:0.05", "z:5"));

  // Paving of the visibility separator
  DefaultFigure::pave(
    {{-1,6},{-1,6}},
    CtcVisible(a, s),
    1e-1
  );
}