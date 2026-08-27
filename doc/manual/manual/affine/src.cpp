/**
 *  Codac tests
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Jordan Ninin
 *  \copyright  Copyright 2026 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include <catch2/catch_test_macros.hpp>
#include <sstream>
#include <concepts>
#include "codac2_Affine.h"
#include "codac2_Matrix.h"
#include "codac2_IntervalMatrix.h"
#include "codac2_Approx.h"
#include <iostream>

using namespace std;
using namespace codac2;

namespace {
  // Small helper concepts used to demonstrate, at compile time, which
  // operations are deliberately unavailable on AffineVariables.
  template<class V>
  concept has_plus_eq_double = requires(V& v, double d) { v += d; };

  template<class V>
  concept has_init_from_affine = requires(V& v, const Affine& a) { v.init(a); };
}

TEST_CASE("Affine class - manual")
{
  {
    // [affine-class-1-beg]
    AffineVariables v(1);
    v[0] = Interval(-1,5);    // a variable affine form, enclosure [-1,5]
                              // v[0] = [-1, 5] : 2 + 3 eps_0 + 0 [-1,1]
    Affine x = v[0];          // a constant affine form, enclosure [-1,5]
                              // x = [-1, 5] : 2 + 3 eps_0 + 0 [-1,1]
    Affine y =Interval(-1,5); // a affine form, enclosure [-1,5], but no dedicated noise symbol
                              // y = [-1, 5] : 2  + 3 [-1,1]
    Affine z(0.);             // a constant affine form, enclosure [0,0]
                              // z = [0, 0] : 0 + 0 [-1,1]
    v[0] = Interval(4.);      // v[0] = [4, 4] : 4 + 0 eps_0 + 0 [-1,1]
    v[0] = Interval(2,3);    // v[0] = [2, 3] : 2.5 + 0.5 eps_0 + 0 [-1,1]
    Affine zz;                 // unbounded, ]-oo,oo[, like Interval()
    // [affine-class-1-end]
    CHECK(x.itv() == Interval(-1,5));
    CHECK(y.itv() == Interval(-1,5));
    CHECK(zz.itv() == Interval());
    CHECK(z.itv() == Interval(0));
  }

  {
    // [affine-class-2-beg]
    // Only AffineVariables introduces new, independent noise symbols.
    // A plain Affine built from an Interval never carries a dedicated
    // symbol: it cannot be correlated with anything else later on.
    AffineVariables v(IntervalVector({{1,2},{-1,1},{3,4}}));
    // v[0] = [1,2] : 0.5 + 0.5 eps_0 + 0 eps_1 + 0 eps_2 + 0 [-1,1],
    // v[1] = [-1,1] : 0 + 0 eps_0 + 1 eps_1 + 0 eps_2 + 0 [-1,1],
    // v[2] = [3,4] : 3.5 + 0 eps_0 + 0 eps_1 + 0.5 eps_2 + 0 [-1,1]
    // [affine-class-2-end]
    CHECK(v[0].itv() == Interval(1,2));
    CHECK(v[1].itv() == Interval(-1,1));
    CHECK(v[2].itv() == Interval(3,4));
  }

  {
    // [affine-class-3-beg]
    AffineVariables v(IntervalVector({{1,2},{-1,1},{3,4}}));

    v[0].itv();          // interval enclosure: [1,2]
    v[0].noise_count();  // total number of noise symbols in the shared context: 3
    v[0].noise(0);       // coefficient of eps_0 : 0.5 (the radius)
    v[0].noise(1);       // coefficient of eps_1 : 0
    v[0].noise(2);       // coefficient of eps_2 : 0
    v[0].err();          // accumulated rounding-error term: 0 (nothing computed yet)
    v[0].noise_index();  // index of the noise symbol associated with v[0] : 0
    v[2].noise_index();  // index of the noise symbol associated with v[2] : 2
    v[1];                   // v[1] = [-1,1] : 0 + 0 eps_0 + 1 eps_1 + 0 eps_2 + 0 [-1,1]
    v[1] = Interval(10,20); // v[1] = [10,20] : 15 + 0 eps_0 + 5 eps_1 + 0 eps_2 + 0 [-1,1]
    // [affine-class-3-end]
    CHECK(v[0].noise_count() == 3);
    CHECK(v[0].noise(0) == 0.5);
    CHECK(v[0].noise(1) == 0.0);
    CHECK(v[0].noise(2) == 0.0);
    CHECK(v[0].err() == 0.0);
    CHECK(v[2].noise_count() == 3);
    CHECK(v[2].noise(0) == 0.0);
    CHECK(v[2].noise(1) == 0.0);
    CHECK(v[2].noise(2) == 0.5);
    CHECK(v[2].err() == 0.0); 
  }

  {
    // [affine-class-4-beg]
    // Addition, subtraction and scaling by a constant are exact in affine
    // arithmetic: no linearization, no loss of tightness, whether or not
    // the operands are correlated.
    AffineVariables var(2);
    var[0] = Interval(2,3);
    Affine x = var[0];   // [2, 3]: 2.5 + 0.5 eps_0 + 0 eps_1 + 0 [-1,1]
    var[1] = Interval(1,3);
    Affine y = var[1];   // [1, 3]: 2 + 0 eps_0 + 1 eps_1 + 0 [-1,1]
    Affine z = x + y;    // [3, 6]: 4.5 + 0.5 eps_0 + 1 eps_1 + 0 [-1,1]
    z = x - 1;           // [1, 2]: 1.5 + 0.5 eps_0 + 0 eps_1 + 0 [-1,1]
    z = 2 * x;           // [4, 6]: 5 + 1 eps_0 + 0 eps_1 + 0 [-1,1]
    z = z - x;           // [2, 3]: 2.5 + 0.5 eps_0 + 0 eps_1 + 0 [-1,1]
    // [affine-class-4-end]
    CHECK((x + y).itv() == Interval(3,6));
    CHECK((x - 1).itv() == Interval(1,2));
    CHECK((2 * x).itv() == Interval(4,6));
  }

  {
    // [affine-class-5-beg]
    // Nonlinear functions rely on a first-order (Chebyshev or MinRange)
    // linearization: the result is a sound, but not necessarily minimal,
    // enclosure of the true range.
    AffineVariables x(IntervalVector({{0, 1.5707963267948966}})); // [0, pi/2]
    Affine y = sin(x[0]); //[-0.0544991, 1.25618] : 0.600837 + 0.52725 eps_0 + 0.128086 [-1,1] 
    Affine z = exp(x[0]); //[0.276127, 4.81048] : 2.54331 + 1.90524 eps_0 + 0.361937 [-1,1] 
    // [affine-class-5-end]
    CHECK(y.itv().is_superset(Interval(0,1)));
    CHECK(z.itv().is_superset(Interval(1,std::exp(1.5707963267948966))));
  }

  {
    // [affine-class-nonlinear-list-beg]
    AffineVariables v(IntervalVector({{1,2},{2,3},{0.5,1.0}}));
    Affine x = v[0];
    Affine y = v[1];
    Affine z = v[2];

    Affine a1 = inv(x);
    Affine a2 = sqr(x);
    Affine a3 = sqrt(x);
    Affine a4 = exp(x);
    Affine a5 = log(x);
    Affine a6 = pow(x, 2);
    Affine a7 = pow(x, 0.5);
    Affine a8 = pow(x, Interval(1,2));
    Affine a9 = pow(x, y);
    Affine a10 = root(x, 2);
    Affine a11 = sin(z);
    Affine a12 = cos(z);
    Affine a13 = tan(z);
    Affine a14 = asin(z);
    Affine a15 = acos(z);
    Affine a16 = atan(z);
    Affine a17 = sinh(z);
    Affine a18 = cosh(z);
    Affine a19 = tanh(z);
    Affine a20 = asinh(z);
    Affine a21 = acosh(x);
    Affine a22 = atanh(z);
    Affine a23 = atan2(z, x);
    Affine a24 = abs(x);
    Interval a25 = floor(x);
    Interval a26 = ceil(x);
    Interval a27 = integer(x);
    Affine a28 = sign(z);
    Affine a29 = chi(z, x, y);
    // [affine-class-nonlinear-list-end]
    (void)a1; (void)a2; (void)a3; (void)a4; (void)a5; (void)a6; (void)a7;
    (void)a8; (void)a9; (void)a10; (void)a11; (void)a12; (void)a13;
    (void)a14; (void)a15; (void)a16; (void)a17; (void)a18; (void)a19;
    (void)a20; (void)a21; (void)a22; (void)a23; (void)a24; (void)a25;
    (void)a26; (void)a27; (void)a28; (void)a29;
  }

  {
    // [affine-class-6-beg]
    // The linearization mode is a global (thread-local) setting shared by
    // every AffineMain instantiation; restore the default after use so
    // that later computations are not silently affected.
    Affine::change_mode(Affine::AF_Lin_MinRange);
    // ... nonlinear computations using the MinRange linearization ...
    Affine::change_mode(Affine::AF_Lin_Chebyshev); // back to the default
    // [affine-class-6-end]
    CHECK(Affine::get_mode() == Affine::AF_Lin_Chebyshev);
  }

  {
    // [affine-class-7-beg]
    // operator== only compares interval enclosures, not the underlying
    // coefficients or dependency structure: two affine forms with an
    // identical enclosure compare equal even if one of them originates
    // from a shared noise symbol and the other does not.
    AffineVariables var(2);
    var.init(Interval(1,2));
    Affine a = var[0];
    Affine b = var[1];
    bool same_enclosure = (a == b); // true
    // [affine-class-7-end]
    CHECK(same_enclosure);
  }

  {
    // [affine-class-8-beg]
    // Because AffineVariables preserves the noise symbol across
    // computations, subtracting a variable from itself is exact: the
    // dependency is not lost the way it would be with plain intervals
    // ([-2,3] - [-2,3] would give [-5,5], not {0}).
    AffineVariables v(2);
    v.init(Interval(-2,3));
    Affine x = v[0] - v[0]; // [0]
    Affine y = v[1] - v[1]; // [0]
    Affine z = v[0] - v[1]; // [-5,5]
    Affine zz= z  - z;      // [0]
    // [affine-class-8-end]
    CHECK(x.itv() == Interval(0.0));
    CHECK(y.itv() == Interval(0.0));
    CHECK(z.itv() == Interval(-5,5));
    CHECK(zz.itv() == Interval(0.0));
  }

  {
    // [affine-class-9-beg]
    // atan2(y,x) is available for affine forms. When the affine
    // division/atan construction is not applicable, the implementation
    // falls back to a sound interval enclosure.
    AffineVariables var(IntervalVector({{1,2},{0.5,1}}));
    Affine x = var[0];
    Affine y = var[1];
    Affine z = atan2(y,x);
    // [affine-class-9-end]
    CHECK(z.itv().is_subset(Interval(0, Interval::half_pi().ub())));
  }

  {
    // [affine-class-10-beg]
    // mig(), mag(), smag() and smig() mirror the equivalent Interval
    // methods, applied to this affine form's own interval enclosure;
    // volume() is an alias for diam() in this scalar, 1-dimensional case.
    AffineVariables v(IntervalVector({{-3,5},{2,4}}));
    Affine x = v[0]; // [-3,5]: straddles zero
    Affine y = v[1]; // [2,4]: strictly positive

    x.mig();    // mignitude: 0 (0 belongs to [-3,5])
    x.mag();    // magnitude: 5 (max(|-3|,|5|))
    x.smag();   // signed magnitude: 5 (the bound with the largest |.|)
    x.smig();   // signed mignitude: 0 (the enclosure straddles zero)
    x.volume(); // == x.diam(): 8

    y.mig();  // mignitude: 2 (y is strictly positive: its lower bound)
    y.mag();  // magnitude: 4
    y.smag(); // signed magnitude: 4
    y.smig(); // signed mignitude: 2 (y is strictly positive: its lower bound)
    // [affine-class-10-end]
    CHECK(x.mig() == 0.0);
    CHECK(x.mag() == 5.0);
    CHECK(x.smag() == 5.0);
    CHECK(x.smig() == 0.0);
    CHECK(x.volume() == x.diam());
    CHECK(x.volume() == 8.0);
    CHECK(y.mig() == 2.0);
    CHECK(y.mag() == 4.0);
    CHECK(y.smag() == 4.0);
    CHECK(y.smig() == 2.0);
  }

  {
    // [affine-class-11-beg]
    // inflate(r) widens the enclosure by [-r,+r] on each side, in place
    // (outward-rounded, so the result may be a hair wider than [0.5,2.5]).
    AffineVariables v(1);
    v[0] = Interval(1,2);
    Affine x = v[0]; // [1,2]
    x.inflate(0.5);  // [0.5,2.5]
    // [affine-class-11-end]
    CHECK(x.itv().is_superset(Interval(0.5,2.5)));
  }

  {
    // [affine-class-12-beg]
    // These predicates mirror the equivalent Interval methods, applied to
    // this affine form's own interval enclosure; each accepts either an
    // Interval or another Affine (contains() takes a plain double).
    AffineVariables v(IntervalVector({{1,2},{0,3},{5,6}}));
    Affine x = v[0]; // [1,2]
    Affine y = v[1]; // [0,3]
    Affine z = v[2]; // [5,6]

    x.is_subset(y);   // true:  [1,2] is a subset of [0,3]
    y.is_superset(x); // true:  [0,3] is a superset of [1,2]
    x.is_disjoint(z); // true:  [1,2] and [5,6] share no point at all
    x.intersects(y);  // true:  [1,2] and [0,3] share at least one point
    x.overlaps(y);    // true:  their intersection has non-zero volume
    x.contains(1.5);  // true:  1.5 belongs to [1,2]
    // [affine-class-12-end]
    CHECK(x.is_subset(y));
    CHECK(y.is_superset(x));
    CHECK(x.is_disjoint(z));
    CHECK(x.intersects(y));
    CHECK(x.overlaps(y));
    CHECK(x.contains(1.5));
  }
}

TEST_CASE("AffineVariables class - manual")
{
  {
    // [affine-variables-complete-beg]
    // 1. Declare uncertain quantities.
    AffineVariables x(IntervalVector({
      {0., 1.},
      {2., 4.},
      {-1., 1.}
    }));

    // 2. Each component receives a distinct noise symbol.
    Index i0 = x[0].noise_index();
    Index i1 = x[1].noise_index();
    Index i2 = x[2].noise_index();

    // 3. Existing variables can be assigned a new interval without
    // changing their noise-symbol identity.
    x[1] = Interval(5., 6.);

    // 4. Use the declared variables in affine expressions.
    Affine a = x[0] + 2. * x[1];
    Affine b = x[1] - x[2];

    // 5. Reusing the same declared variable preserves its dependency.
    Affine c = x[0] - x[0];

    // 6. Materialize the declared variables as an AffineVector when vector
    // operations are required.
    AffineVector v(x);
    AffineVector y = 2. * v;
    // [affine-variables-complete-end]

    CHECK(i0 == 0);
    CHECK(i1 == 1);
    CHECK(i2 == 2);
    CHECK(x[1].noise_index() == i1);
    CHECK(c.itv() == Interval(0.0));
    CHECK(y.itv() == IntervalVector({{0.,2.},{10.,12.},{-2.,2.}}));
  }

  {
    // [affine-variables-1-beg]
    AffineVariables u(3);                  // 3 unbound components, ]-oo,+oo[ each
    AffineVariables v(IntervalVector({{1,2},{-1,1},{3,4}}));
    AffineVariables w(Vector({1.5,0.,3.5})); // degenerate (point) components
    // [affine-variables-1-end]
    CHECK(u.size() == 3);
    CHECK(u[0].itv() == Interval());
    CHECK(v[0].itv() == Interval(1,2));
    CHECK(w[0].itv() == Interval(1.5,1.5));
  }

  {
    // [affine-variables-2-beg]
    // Each component owns a dedicated, distinct noise symbol, exposed by
    // noise_index(). Re-assigning a component from an Interval or a double
    // replaces its value but keeps that same noise_index() unchanged.
    AffineVariables v(IntervalVector({{1,2},{-1,1}}));
    Index i0 = v[0].noise_index(); // 0
    Index i1 = v[1].noise_index(); // 1
    v[0] = Interval(5,6);          // v[0] now encloses [5,6]...
    // [affine-variables-2-end]
    CHECK(i0 == 0);
    CHECK(i1 == 1);
    CHECK(v[0].noise_index() == i0); // ...but keeps its own symbol
    CHECK(v[0].itv() == Interval(5,6));
  }

  {
    // [affine-variables-assignment-beg]
    AffineVariables v(IntervalVector({{1.,2.},{-1.,1.}}));

    Index i0 = v[0].noise_index();
    v[0] = Interval(5.,6.);
    v[1] = 3.;

    // The numeric values change, but the declared identities remain.
    // [affine-variables-assignment-end]
    CHECK(v[0].noise_index() == i0);
    CHECK(v[1].noise_index() == 1);
    CHECK(v[0].itv() == Interval(5.,6.));
    CHECK(v[1].itv() == Interval(3.));
  }

  {
    // [affine-variables-dependency-beg]
    AffineVariables v(IntervalVector({{1.,2.}}));
    Affine a = v[0];
    Affine b = v[0];
    Affine correlated = a - b;

    // These two plain Affine objects are created independently.
    Affine p(Interval(1.,2.));
    Affine q(Interval(1.,2.));
    Affine independent = p - q;
    // [affine-variables-dependency-end]
    CHECK(correlated.itv() == Interval(0.0));
    CHECK(independent.itv() == Approx<Interval>(Interval(-1.,1.)));
  }

  {
    // [affine-variables-3-beg]
    // resize() drops any existing dependency information and reconstructs
    // every component as a fresh, independent, unbounded variable.
    // conservativeResize() instead keeps each surviving component's
    // interval enclosure (not its correlations), and initializes any newly
    // added component to ]-oo,+oo[.
    AffineVariables v(IntervalVector({{1,2},{-1,1},{3,4}}));
    Index old_i0 = v[0].noise_index();

    v.conservativeResize(4);
    // v[0..2] still enclose their original intervals, v[3] is unbounded.
    // Their previous noise-symbol identities are intentionally not kept.
    Index new_i0 = v[0].noise_index();

    // resize() discards even the surviving interval enclosures and creates
    // a completely fresh set of unbounded variables.
    v.resize(2);
    // [affine-variables-3-end]
    CHECK(v.size() == 2);
    CHECK(v[0].itv() == Interval());
    CHECK(old_i0 == new_i0);
    CHECK(v[0].noise_index() == new_i0);
    CHECK(v[1].itv() == Interval());
  }

  {
    // [affine-variables-4-beg]
    // init(Interval) broadcasts the same interval to every component,
    // while preserving each component's own dedicated noise symbol.
    AffineVariables v(3);
    v.init(Interval(1,2));
    // [affine-variables-4-end]
    for(Index i = 0; i < v.size(); ++i)
    {
      CHECK(v[i].itv() == Interval(1,2));
      CHECK(v[i].noise_index() == i);
    }
  }

  {
    // [affine-variables-5-beg]
    // To use the declared variables in a computation, convert them to a
    // plain AffineVector (or let arithmetic operators do it implicitly):
    // the conversion goes through AffineMain's own copy constructor and
    // therefore preserves every noise symbol exactly.
    AffineVariables v(IntervalVector({{1,2},{-1,1},{3,4}}));
    AffineVector x(v); // safe: element-wise AffineMain copy construction
    // [affine-variables-5-end]
    CHECK(x[0].itv() == v[0].itv());
  }

  {
    // [affine-variables-6-beg]
    // Whole-vector assignment between two AffineVariables (w = v;) falls
    // through to Eigen's generic, component-wise matrix assignment, which
    // calls AffineVarMain::operator=(const AffineVarMain&) per component.
    // This copies each component's noise symbol identity along with its
    // value: after the assignment, w[i] and v[i] refer to the exact same
    // noise symbol, and are therefore seen as fully correlated by any
    // later computation mixing both.
    AffineVariables v(IntervalVector({{1,2},{-1,1}}));
    AffineVariables w(2);
    w = v;
    Affine d = w[0] - v[0];
    // [affine-variables-6-end]
    CHECK(w[0].noise_index() == v[0].noise_index());
    CHECK(d.itv() == Interval(0.0));

    // A copy is therefore not a declaration of independent uncertainty:
    // w[0] and v[0] remain fully correlated.
  }

  {
    // [affine-variables-7-beg]
    // operator=(const IntervalVector&) reassigns every component from the
    // given interval vector, resizing this vector first if the sizes
    // differ. Every component is rebuilt fresh (like resize(), not
    // conservativeResize()): no dependency from before the assignment
    // survives.
    AffineVariables v(2);
    v = IntervalVector({{1,2},{-1,1},{3,4}}); // v is resized to 3
    // [affine-variables-7-end]
    CHECK(v.size() == 3);
    CHECK(v[2].itv() == Interval(3,4));
    CHECK(v[2].noise_index() == 2);
  }

  {
    // [affine-variables-8-beg]
    // Compound assignment, and re-binding a component or the whole vector
    // to a foreign Affine/AffineVector, are deleted on purpose: they would
    // either mutate a declared variable in place, or hand it a value that
    // fabricates a dependency it never had.
    static_assert(!has_plus_eq_double<AffineVariables>);
    static_assert(!has_init_from_affine<AffineVariables>);
    // [affine-variables-8-end]
  }

  {
    // [affine-variables-9-beg]
    // AffineVariables can appear directly as an operand of a matrix
    // product: it is converted internally to an AffineVector, so the
    // result always carries Affine (not AffineVarMain) coefficients.
    Matrix A(2,2);
    A(0,0) = 1.; A(0,1) = 1.;
    A(1,0) = 0.; A(1,1) = 1.;
    AffineVariables v(IntervalVector({{-1,1},{-1,1}}));
    AffineVector y = A * v;
    // [affine-variables-9-end]
    CHECK(y.itv() == IntervalVector({{-2,2},{-1,1}}));
  }

  {
    // [affine-variables-10-beg]
    // Unary minus also converts to a plain AffineVector: negating a
    // declared variable no longer represents "the" variable itself, but a
    // derived expression.
    AffineVariables v(IntervalVector({{1,2},{-1,1}}));
    AffineVector m = -v;
    // [affine-variables-10-end]
    CHECK(m.itv() == IntervalVector({{-2,-1},{-1,1}}));
  }

  {
    // [affine-variables-11-beg]
    AffineVariables v(IntervalVector({{1,2},{-1,1}}));
    std::ostringstream oss;
    oss << v;
    // [affine-variables-11-end]
    CHECK(!oss.str().empty());
  }

  {
    // [affine-variables-12-beg]
    // Linear combination of two DIFFERENT declared variables: both noise
    // symbols show up in the result, each with its own coefficient.
    AffineVariables v(IntervalVector({{1,2},{3,4}}));
    Affine z = v[0] + v[1];
    // Expected `cout << z` output:
    //   [4, 6] : 5 + 0.5 eps_0 + 0.5 eps_1 + 0 [-1,1]
    // [affine-variables-12-end]
    CHECK(z.itv() == Interval(4,6));
    CHECK(z.noise(0) == 0.5);
    CHECK(z.noise(1) == 0.5);
    CHECK(z.err() == 0.0);
  }

  {
    // [affine-variables-13-beg]
    // Linear combination with real coefficients and a constant term: the
    // constant only shifts the midpoint, each variable's own coefficient
    // is scaled independently.
    AffineVariables v(IntervalVector({{0,2},{0,4}}));
    Affine z = 2.*v[0] - 0.5*v[1] + 1.;
    // Expected `cout << z` output:
    //   [-1, 5] : 2 + 2 eps_0 + -1 eps_1 + 0 [-1,1]
    // [affine-variables-13-end]
    CHECK(z.itv() == Interval(-1,5));
    CHECK(z.noise(0) == 2.0);
    CHECK(z.noise(1) == -1.0);
  }

  {
    // [affine-variables-14-beg]
    // Reusing the same declared variable several times still collapses
    // exactly, whatever the linear combination (see also the AffineVector
    // page for the same property at vector level):
    // 4x - 2(x+x) + (x-x) = 4x - 4x + 0 = 0 for any x.
    AffineVariables v(IntervalVector({{-2,3}}));
    Affine x = v[0];
    Affine z = 4.*x - 2.*(x + x) + (x - x);
    // Expected `cout << z` output:
    //   [0, 0] : 0 + 0 eps_0 + 0 [-1,1]
    // [affine-variables-14-end]
    CHECK(z.itv() == Interval(0.0));
    CHECK(z.noise(0) == 0.0);
    CHECK(z.err() == 0.0);
  }

  {
    // [affine-variables-15-beg]
    // Nonlinear function of a single declared variable: sqr(x) = x*x is
    // computed by directly expanding the product rather than through a
    // generic Chebyshev linearization, so it stays exact here (the
    // remainder term only has to account for the eps_0^2 <= 1 cross term,
    // not for a derivative approximation).
    AffineVariables v(IntervalVector({{1,3}})); // x0 = 2, radius = 1
    Affine x = v[0];
    Affine y = sqr(x);
    // Expected `cout << y` output (hand-traced from the source: new
    // midpoint x0^2 + 0.5*coeff^2 = 4 + 0.5, new coefficient 2*x0*coeff = 4,
    // remainder 0.5 plus a floating-point safety margin far below display
    // precision):
    //   [0, 9] : 4.5 + 4 eps_0 + 0.5 [-1,1]
    // [affine-variables-15-end]
    CHECK(y.itv().is_superset(Interval(1,9))); // true range of x^2 on [1,3]
    CHECK(y  == Approx<Affine>(Interval(0,9),1.e-6));
  }

  {
    // [affine-variables-16-beg]
    // Transcendental functions (sin, cos, exp, ...) go through a
    // Chebyshev/MinRange linearization with its own floating-point error
    // tracking; unlike the linear and sqr() cases above, the exact printed
    // coefficients are not reproduced here (they were not hand-verified
    // against the source) -- only the guaranteed containment is checked.
    // Run this snippet to see the actual printed values.
    AffineVariables v(IntervalVector({{0, 1.5707963267948966}})); // [0, pi/2]
    Affine y = sin(v[0]);  //  [-0.0544991, 1.25618] : 0.600837 + 0.52725 eps_0 + 0.128086 [-1,1]
    // The exact coefficients depend on the selected linearization mode.
    // [affine-variables-16-end]
    CHECK(y.itv().is_superset(Interval(0,1)));
  }
}


TEST_CASE("AffineVector class - manual")
{
  {
    // [affine-vector-1-beg]
    AffineVariables v(IntervalVector({{1,2},{-1,1},{3,4}}));
    AffineVector x(v);          // materialized from declared variables
    AffineVector y = x + x;     // ordinary Eigen vector expression
    IntervalVector box = x.itv();
    // [affine-vector-1-end]
    CHECK(box == IntervalVector({{1,2},{-1,1},{3,4}}));
    CHECK(y.itv() == IntervalVector({{2,4},{-2,2},{6,8}}));
  }

  {
    // [affine-vector-2-beg]
    // Dependency is preserved through vector-level linear combinations
    // exactly as it is for a single Affine (see "Why affine arithmetic" in
    // the Affine class page): reusing the same declared variables several
    // times does not widen the result the way plain interval arithmetic
    // would.
    AffineVariables v(IntervalVector({{-2,3},{0,1}}));
    AffineVector x(v);
    AffineVector z = 4.*x - 2.*(x + x) + (x - x);
    // [affine-vector-2-end]
    CHECK(z.itv() == IntervalVector({{0,0},{0,0}}));
  }

  {
    // [affine-vector-3-beg]
    // Unlike AffineVariables, compound assignment operators are available
    // on a plain AffineVector, since its scalar type (Affine, not a bound
    // AffineVarMain) does not delete them.
    AffineVariables v(IntervalVector({{1,2},{3,4}}));
    AffineVector x(v);
    x += 1.;
    x *= 2.;
    // [affine-vector-3-end]
    CHECK(x.itv() == IntervalVector({{4,6},{8,10}}));
  }

  {
    // [affine-vector-4-beg]
    // AffineVector is a plain Eigen column vector, so the usual block
    // accessors are available and return an AffineVector over the
    // requested range.
    AffineVariables v(IntervalVector({{1,2},{-1,1},{3,4},{5,6}}));
    AffineVector x(v);
    AffineVector first_two = x.head(2);
    AffineVector last_two  = x.tail(2);
    AffineVector middle    = x.segment(1,2);
    // [affine-vector-4-end]
    CHECK(first_two.itv() == IntervalVector({{1,2},{-1,1}}));
    CHECK(last_two.itv()  == IntervalVector({{3,4},{5,6}}));
    CHECK(middle.itv()    == IntervalVector({{-1,1},{3,4}}));
  }

  {
    // [affine-vector-5-beg]
    // transpose() turns an AffineVector into an AffineRow, like it would
    // for any Eigen column vector.
    AffineVariables v(IntervalVector({{1,2},{-1,1}}));
    AffineVector x(v);
    AffineRow r = x.transpose();
    // [affine-vector-5-end]
    CHECK(r(0).itv() == Interval(1,2));
    CHECK(r(1).itv() == Interval(-1,1));
  }

  {
    // [affine-vector-6-beg]
    // Applying an elementary function component-wise: there is no
    // vectorized shortcut for nonlinear functions (unlike the linear
    // compound assignments above), so it is done component by component.
    AffineVariables v(IntervalVector({{0,1},{1,2}}));
    AffineVector x(v);
    AffineVector y(x.size());
    for(Index i = 0; i < x.size(); ++i)
      y[i] = sin(x[i]);
    // [affine-vector-6-end]
    CHECK(y[0].itv().is_superset(Interval(0, 0.8)));    // true range: [0, sin(1)]
    CHECK(y[1].itv().is_superset(Interval(0.85, 0.99))); // true range: [sin(1), 1] (pi/2 in [1,2])
  }

  {
    // [affine-vector-7-beg]
    AffineVariables v(IntervalVector({{1,2},{3,4}}));
    AffineVector x(v);
    std::ostringstream oss;
    oss << x;
    // Expected `cout << x` output (each component prints its own full
    // Affine representation, joined by " ; " and wrapped in "[ ... ]" --
    // this is why itv() is usually preferred for a quick summary):
    //   [ [1, 2] : 1.5 + 0.5 eps_0 + 0 eps_1 + 0 [-1,1]  ;
    //     [3, 4] : 3.5 + 0 eps_0 + 0.5 eps_1 + 0 [-1,1]  ]
    // [affine-vector-7-end]
    CHECK(!oss.str().empty());
  }
}

TEST_CASE("AffineMatrix class - manual")
{
  {
    // [affine-matrix-1-beg]
    Matrix A(2,2);
    A(0,0) = 1.; A(0,1) = 1.;
    A(1,0) = 0.; A(1,1) = 1.;

    AffineVariables v(IntervalVector({{-1,1},{-1,1}}));
    AffineVector x(v);

    AffineVector y = A * x;    // real matrix times a vector of affine forms
    IntervalVector box = y.itv();
    // [affine-matrix-1-end]
    CHECK(box == IntervalVector({{-2,2},{-1,1}}));
  }

  {
    // [affine-matrix-2-beg]
    // Dependency preservation extends to matrix products: evaluating the
    // same declared variables through two matrices and recombining them
    // linearly does not lose correlation, which is typically where affine
    // arithmetic pays off compared to interval arithmetic in practice
    // (e.g. repeated evaluation of the same uncertain state vector through
    // several matrices, as in a state-space / control computation).
    Matrix A(2,2);
    A(0,0) = 1.; A(0,1) = 0.;
    A(1,0) = 0.; A(1,1) = 1.;

    AffineVariables v(IntervalVector({{-1,1},{-1,1}}));
    AffineVector x(v);

    AffineVector diff = AffineVector(A * x) - AffineVector(A * x);
    // [affine-matrix-2-end]
    CHECK(diff.itv() == IntervalVector({{0,0},{0,0}}));
  }

  {
    // [affine-matrix-3-beg]
    // AffineMainMatrix<T> (alias AffineMatrix) is a plain Eigen matrix
    // alias, like AffineVector; itv() converts it to an IntervalMatrix.
    AffineMatrix M(2,2);
    AffineVariables var(IntervalVector({{1,2},{0},{0},{1,2}}));
    M(0,0) = var[0];
    M(0,1) = var[1];
    M(1,0) = var[2];
    M(1,1) = var[3];
    IntervalMatrix box = M.itv();
    // [affine-matrix-3-end]
    CHECK(box(0,0) == Interval(1,2));
    CHECK(box(0,1) == Interval(0,0));
  }

  {
    // [affine-matrix-shared-beg]
    AffineVariables v(2);
    v[0] = Interval(-1,1);
    v[1] = Interval(2,3);

    AffineMatrix M(2,2);
    M(0,0) = v[0];
    M(0,1) = v[1];
    M(1,0) = v[0];
    M(1,1) = v[1];
    // [affine-matrix-shared-end]
  }

  {
    // [affine-matrix-4-beg]
    // transpose(), like on any Eigen matrix, swaps rows and columns.
    AffineMatrix M(2,2);
    M(0,0) = Affine(Interval(1,2)); M(0,1) = Affine(3.);
    M(1,0) = Affine(-1.);           M(1,1) = Affine(Interval(0,1));
    AffineMatrix Mt = M.transpose();
    // [affine-matrix-4-end]
    CHECK(Mt(0,1).itv() == M(1,0).itv());
    CHECK(Mt(1,0).itv() == M(0,1).itv());
  }

  {
    // [affine-matrix-5-beg]
    // Row and column extraction return an AffineRow / AffineVector.
    AffineMatrix M(2,2);
    M(0,0) = Affine(Interval(1,2)); M(0,1) = Affine(3.);
    M(1,0) = Affine(-1.);           M(1,1) = Affine(Interval(0,1));
    AffineRow    row0 = M.row(0);
    AffineVector col1 = M.col(1);
    // [affine-matrix-5-end]
    CHECK(row0(0).itv() == Interval(1,2));
    CHECK(row0(1).itv() == Interval(3,3));
    CHECK(col1(0).itv() == Interval(3,3));
    CHECK(col1(1).itv() == Interval(0,1));
  }

  {
    // [affine-matrix-products-beg]
    AffineVariables v(IntervalVector({{-1,1},{2,3}}));

    AffineMatrix M(2,2);
    M(0,0) = v[0];
    M(0,1) = v[1];
    M(1,0) = v[0];
    M(1,1) = v[1];

    AffineVector x(v);
    AffineVector y1 = M * x;

    AffineMatrix N(2,2);
    N.setIdentity();
    AffineVector y2 = N * x;
    AffineMatrix P = M * N;
    // [affine-matrix-products-end]
    CHECK(y1.size() == 2);
    CHECK(y2.size() == 2);
    CHECK(P.rows() == 2);
    CHECK(P.cols() == 2);
  }

  {
    // [affine-matrix-6-beg]
    AffineMatrix M(1,2);
    M(0,0) = Affine(Interval(1,2));
    M(0,1) = Affine(3.);
    std::ostringstream oss;
    oss << M;
    // Expected `cout << M` format (rows separated by newlines, each
    // component printing its own full Affine representation -- verbose,
    // as for AffineVector; itv() is usually preferred for a quick view).
    // Both M(0,0) and M(0,1) were built as freestanding Affine(Interval)/
    // Affine(double) constants (not through AffineVariables), so neither
    // carries an eps_i term: all uncertainty sits in the error bound.
    //   [ [1, 2] : 1.5 + 0.5 [-1,1]  , [3, 3] : 3 + 0 [-1,1]  ]
    // [affine-matrix-6-end]
    CHECK(!oss.str().empty());
  }
}

