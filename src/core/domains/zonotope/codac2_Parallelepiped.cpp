/** 
 *  codac2_Parallelepiped.cpp
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Maël Godard
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include "codac2_Parallelepiped.h"
#include "codac2_inversion.h"

using namespace codac2;

Parallelepiped::Parallelepiped(const Vector& c_, const Matrix& A_)
  : Zonotope(c_, A_)
{
  assert_release(A.cols() <= c.size() && "too many vectors, you are describing a zonotope");
}

void generate_vertices(Index i, Index n, const Vector& c, const Matrix& A, std::vector<Vector>& L_v)
{
  if (i == n)
  {
    L_v.push_back(c);
  }
  else if (i<n)
  {
    generate_vertices(i+1, n, c + A.col(i), A, L_v);
    generate_vertices(i+1, n, c - A.col(i), A, L_v);
  }
}

std::vector<Vector> Parallelepiped::vertices() const
{
  std::vector<Vector> L_v;
  generate_vertices(0, c.size(),c,A,L_v);
  return L_v;
}

BoolInterval Parallelepiped::contains(const Vector& v) const
{
  return is_superset(v.template cast<Interval>());
}

BoolInterval Parallelepiped::is_superset(const IntervalVector& x) const
{
  assert_release(A.rows() == A.cols() && "Matrix A must be square to check containment.");
  assert_release(x.size() == c.size() && "Point dimension must match parallelepiped dimension.");

  IntervalVector B = inverse_enclosure(A)*(x - c);
  IntervalVector IV = IntervalVector::constant(A.cols(),{-1,1});

  if (!(B.intersects(IV)))
    return BoolInterval::FALSE;

  else
  {
    if (B.is_subset(IV))
      return BoolInterval::TRUE;
    else
      return BoolInterval::UNKNOWN;
  }
}