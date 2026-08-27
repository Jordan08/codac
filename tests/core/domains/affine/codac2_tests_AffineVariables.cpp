/**
 * \file codac2_tests_AffineVariables.cpp
 * \brief Tests specific to AffineVarMain<T>/AffineVarMainVector<T>, the
 *        "affine variable" type that introduces fresh, uniquely-indexed
 *        noise symbols: the intended assignment asymmetry with
 *        AffineMain<T> (AffineMain <- AffineVarMain is allowed via
 *        slicing, the converse is deleted), distinct noise symbols per
 *        component, copy/assignment across different sizes or indices,
 *        the unsupported scalar compound assignments, and transparent
 *        conversion to AffineMainVector.
 */

#include <catch2/catch_test_macros.hpp>

#include <concepts>
#include <memory>
#include <sstream>
#include <type_traits>

#include "codac2_Affine.h"
#include "codac2_AffineMatrix.h"
#include "codac2_AffineRow.h"
#include "codac2_AffineVector.h"
#include "codac2_IntervalMatrix.h"
#include "codac2_IntervalRow.h"
#include "codac2_IntervalVector.h"
#include "codac2_Approx.h"
#include <iostream>

using namespace codac2;

namespace {

using Model = AF_Default;
using AffineT = AffineMain<Model>;
using AffineTMatrix = AffineMainMatrix<Model>;
using AffineTVector = AffineMainVector<Model>;
using AffineTVarVector = AffineVarMainVector<Model>;
using AffineTRow = AffineMainRow<Model>;

typedef AF_Default AA;

AffineTVarVector make_variable_vector_5()
{
    return AffineTVarVector(IntervalVector({
        {-2.0, -1.0},
        {0.0, 2.0},
        {3.0, 4.0},
        {-1.0, 1.0},
        {5.0, 7.0}
    }));
}

AffineTVarVector make_variable_vector_3()
{
    return AffineTVarVector(
        IntervalVector({{1.0, 2.0}, {-1.0, 1.0}, {3.0, 4.0}}));
}

} // namespace

// AffineTVarMain deliberately deletes scalar compound assignments. Keep
// this assertion here so that this API contract is explicit and cannot
// accidentally become enabled by the Eigen plugin.
template<class T> concept has_scalar_plus_eq  = requires(T& v, double d) { v += d; };
template<class T> concept has_scalar_minus_eq = requires(T& v, double d) { v -= d; };
template<class T> concept has_scalar_times_eq = requires(T& v, double d) { v *= d; };
template<class T> concept has_scalar_div_eq   = requires(T& v, double d) { v /= d; };

template<class T>
concept can_assign_AffineT_from_var =
    requires(AffineMain<T>& dst, const AffineVarMain<T>& src) { dst = src; };

template<class T>
concept can_assign_var_from_AffineT =
    requires(AffineVarMain<T>& dst, const AffineMain<T>& src) { dst = src; };

TEST_CASE(
  "AffineVarMain : affectction entre indices différents"
)
{

// Une affectction entre indices différents conserve l’identité
// de la destination et ne crée pas une fausse dépendance.
AffineTVarVector x(
    IntervalVector({{1.0, 2.0}, {3.0, 4.0}})
);

x[0] = x[1];

CHECK(x[0].itv().is_superset(Interval(3.0, 4.0)));

AffineT result = x[0] - x[1];
CHECK(result == Interval::zero());
}



TEST_CASE(  "AffineVarMain : copie avec des tailles différentes")
{
AffineTVarVector small( IntervalVector({{1.0, 2.0}}));

AffineTVarVector large( IntervalVector({{3.0, 4.0}, {5.0, 6.0}, {7.0, 8.0}}));

small[0] = large[2];

CHECK(small[0].noise_count() == 3);
CHECK(small[0].itv()==large[2].itv());
}



TEST_CASE("Eigen AffineTVarMainVector itv()", "[Eigen][itv][AffineTVarMainVector]")
{
    const AffineTVarVector variables = make_variable_vector_5();
    const auto boxes = variables.itv();

    REQUIRE(boxes.size() == 5);
    CHECK(boxes[0] == Interval(-2.0, -1.0));
    CHECK(boxes[1] == Interval(0.0, 2.0));
    CHECK(boxes[2] == Interval(3.0, 4.0));
    CHECK(boxes[3] == Interval(-1.0, 1.0));
    CHECK(boxes[4] == Interval(5.0, 7.0));
}


TEST_CASE("AffineTVarMain cannot be assigned from AffineTMain", "[AffineT][AffineTVar][regression]")
{
    static_assert(can_assign_AffineT_from_var<Model>,
                  "AffineTMain should remain assignable from AffineTVarMain (slicing)");
    static_assert(!can_assign_var_from_AffineT<Model>,
                  "AffineTVarMain must not be assignable from AffineTMain");
}


TEST_CASE("AffineVarMain empty() factory produces an empty affine variable")
{
    const AffineVarMain<AA> empty_var = AffineVarMain<AA>::empty();

    CHECK(empty_var.is_empty());
    CHECK(empty_var.itv().is_empty());
    CHECK(empty_var.noise_index() == -1);
}


template<class T>
concept can_init_from_interval =
    requires(AffineTVarVector& v, const Interval& x) { v.init(x); };

template<class T>
concept can_init_from_var =
    requires(AffineTVarVector& v, const AffineVarMain<T>& x) { v.init(x); };

template<class T>
concept can_init_from_affine =
    requires(AffineTVarVector& v, const AffineMain<T>& x) { v.init(x); };

TEST_CASE("AffineVarMainVector::init is only enabled for an Interval argument",
          "[AffineTVarVector][init][regression]")
{
    static_assert(can_init_from_interval<Model>,
                  "init(const Interval&) must remain available");
    static_assert(!can_init_from_var<Model>,
                  "init(const AffineVarMain<T>&) must stay deleted: broadcasting it "
                  "would require duplicating its noise symbol across components");
    static_assert(!can_init_from_affine<Model>,
                  "init(const AffineMain<T>&) must stay deleted: it would require an "
                  "implicit conversion to AffineVarMain<T> that is intentionally absent");
}


TEST_CASE("AffineVarMainVector::resize discards previous values, unlike conservativeResize",
          "[AffineTVarVector][resize][regression]")
{
    const IntervalVector box({{1.0, 2.0}, {3.0, 4.0}});

    AffineTVarVector preserved(box);
    preserved.conservativeResize(3);
    CHECK(preserved[0].itv() == Interval(1.0, 2.0));
    CHECK(preserved[1].itv() == Interval(3.0, 4.0));
    CHECK(preserved[2].itv() == Interval());

    AffineTVarVector reset(box);
    reset.resize(3);
    CHECK(reset[0].itv() == Interval());
    CHECK(reset[1].itv() == Interval());
    CHECK(reset[2].itv() == Interval());

    // resize() also reassigns fresh, sequential noise symbols.
    for (Index i = 0; i < reset.size(); ++i)
        CHECK(reset[i].noise_index() == i);
}


TEST_CASE("AffineVarMainVector operator<< streams the interval box")
{
    AffineTVarVector x(IntervalVector({{1.0, 2.0}, {3.0, 4.0}}));
    std::ostringstream stream;
    stream << x;
    CAPTURE(stream.str());
    CHECK_FALSE(stream.str().empty());
    CHECK(stream.str().find("empty") == std::string::npos);

    AffineTVarVector e(1);
    e[0] = Interval::empty();
    std::ostringstream empty_stream;
    empty_stream << e;
    CHECK(empty_stream.str() == "[ empty 1d box ]");
}


TEST_CASE("AffineTVarVector init(Interval) preserves distinct noise symbols",
          "[AffineT][AffineTVarVector][init][regression]")
{
    AffineTVarVector v(3);
    v.init(Interval(1.0, 2.0));

    REQUIRE(v.size() == 3);

    const double rad = Interval(1.0, 2.0).rad(); // 0.5

    for (int i = 0; i < 3; ++i)
    {
        CAPTURE(i, v[i].itv());

        // Every component carries the same interval enclosure...
        CHECK(v[i].itv() == Interval(1.0, 2.0));
        CHECK(v[i].rad() == rad);
        CHECK(v[i].noise_count() == 3);
        // ...but each keeps its own, distinct noise symbol: component i has
        // a coefficient equal to the interval's radius on eps_i, and zero
        // on every other eps_j. A regression that broadcast a single shared
        // symbol (or degraded to a plain interval constant) would show up
        // here as a non-diagonal pattern.
        for (int j = 0; j < 3; ++j)
            CHECK(v[i].noise(j) == (i == j ? rad : 0.0));
    }

    v.conservativeResize(5);

    for (int i = 0; i < 3; ++i)
    {
        CAPTURE(i, v[i].itv());

        // Every component carries the same interval enclosure...
        CHECK(v[i].itv() == Interval(1.0, 2.0));
        CHECK(v[i].rad() == rad);
        CHECK(v[i].noise_count() == 5);

        // ...but each keeps its own, distinct noise symbol: component i has
        // a coefficient equal to the interval's radius on eps_i, and zero
        // on every other eps_j. A regression that broadcast a single shared
        // symbol (or degraded to a plain interval constant) would show up
        // here as a non-diagonal pattern.
        for (int j = 0; j < 5; ++j)
            CHECK(v[i].noise(j) == (i == j ? rad : 0.0));
    }

    for (int i = 3; i < 5; ++i)
    {
        CAPTURE(i, v[i].itv());

        // Every component carries the same interval enclosure...
        CHECK(v[i].itv() == Interval());
        CHECK(v[i].noise_count() == 5);
    }
}



TEST_CASE("Eigen AffineTVarMainVector segment conversion")
{
    const AffineTVarVector variables = make_variable_vector_5();
    const AffineTVector segment = variables.segment(1, 3);

    REQUIRE(segment.size() == 3);
    CHECK(segment[0].itv() == variables[1].itv());
    CHECK(segment[1].itv() == variables[2].itv());
    CHECK(segment[2].itv() == variables[3].itv());
}


TEST_CASE("Eigen scalar compound assignments - AffineTVarMainVector are unsupported",
          "[Eigen][compound-assignment][AffineTVarMainVector]")
{


	static_assert(!has_scalar_plus_eq<AffineTVarVector>);
	static_assert(!has_scalar_minus_eq<AffineTVarVector>);
	static_assert(!has_scalar_times_eq<AffineTVarVector>);
	static_assert(!has_scalar_div_eq<AffineTVarVector>);
}



TEST_CASE("AffineTVarMainVector converts transparently to AffineTMainVector")
{
    const AffineTVarVector variables = make_variable_vector_3();

    SECTION("construction")
    {
        const AffineTVector values(variables);
        REQUIRE(values.size() == variables.size());
        for(Eigen::Index i = 0; i < values.size(); ++i)
            CHECK(values(i).itv() == variables(i).itv());
    }

    SECTION("assignment")
    {
        AffineTVector values;
        values = variables;
        REQUIRE(values.size() == variables.size());
        for(Eigen::Index i = 0; i < values.size(); ++i)
            CHECK(values(i).itv() == variables(i).itv());
    }
}



TEST_CASE("AffineVarMainVector: construction and assignment from IntervalVector")
{
	AffineVarMainVector<AA> x(2);

	const IntervalVector values({
	    {-1.0, 1.0},
	    {2.0, 3.0},
	    {4.0, 5.0}
	});

	x = values;

	CHECK(x.size() == 3);

	for (Index i = 0; i < x.size(); ++i) {
	    CHECK(x[i].noise_count() == 3);
	    CHECK(x[i].itv() == values[i]);
	}
}




TEST_CASE(
    "AffineVarMainVector does not inherit Eigen constructors",
    "[AffineVarMainVector]"
)
{
    AffineVarMainVector<AF_fAF2> x(3);

    CHECK(x[0].noise_index() == 0);
    CHECK(x[1].noise_index() == 1);
    CHECK(x[2].noise_index() == 2);


    AffineVarMainVector<AF_fAF2> y(IntervalVector({{-1., 1.}, {2., 3.}, {4., 5.}}));
    CHECK(y[0].noise_index() == 0);
    CHECK(y[1].noise_index() == 1);
    CHECK(y[2].noise_index() == 2);

    AffineVarMainVector<AF_fAF2> z(Vector({-1., 2., 3.}));
    CHECK(z[0].noise_index() == 0);
    CHECK(z[1].noise_index() == 1);
    CHECK(z[2].noise_index() == 2);

    AffineVarMainVector<AF_fAF2> w(0);
    w.resize(3);
    CHECK(w[0].noise_index() == 0);
    CHECK(w[1].noise_index() == 1);
    CHECK(w[2].noise_index() == 2);
}


