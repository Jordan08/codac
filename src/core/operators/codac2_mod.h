/** 
 *  \file codac2_mod.h
 * ----------------------------------------------------------------------------
 *  \date       2024
 *  \author     Simon Rohou, Damien Massé
 *  \copyright  Copyright 2024 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include "codac2_Interval.h"
#include "codac2_AnalyticType.h"
#include "codac2_AnalyticExprWrapper.h"
#include "codac2_arith_sub.h"

namespace codac2
{
  struct ModOp
  {
    template<typename X1,typename P>
    static std::string str(const X1& x1, const P& p)
    {
      return "mod(" + x1->str() + "," + p->str() + ")";
    }

    template<typename X1,typename P>
    static std::pair<Index,Index> output_shape([[maybe_unused]] const X1& s1, [[maybe_unused]] const P& p)
    {
      return {1,1};
    }

    static Interval fwd(const Interval& x1, const Interval& p);
    static ScalarType fwd_natural(const ScalarType& x1, const ScalarType& p);
    static ScalarType fwd_centered(const ScalarType& x1, const ScalarType& p);
    static void bwd(const Interval& y, Interval& x1, Interval& p);

    static void fwd_bwd(Interval& x1, Interval& x2, double p);
    static void fwd_bwd(Interval& x1, Interval& x2, Interval& p);
  };
  // Analytic operator
  // The following function can be used to build analytic expressions.

  inline ScalarExpr
  mod(const ScalarExpr& x1, const ScalarExpr& p)
  {
    return { std::make_shared<AnalyticOperationExpr<ModOp,ScalarType,ScalarType,ScalarType>>(x1,p) };
  }

  // Inline functions

  inline Interval ModOp::fwd(const Interval& x1, const Interval& p)
  {
    Interval x2, x1_(x1), p_(p);
    ModOp::fwd_bwd(x2,x1_,p_);
    return x2;
  }

  inline ScalarType ModOp::fwd_natural(const ScalarType& x1, const ScalarType& p)
  {
    return {
      fwd(x1.a, p.a),
      x1.def_domain && p.def_domain
    };
  }

  inline ScalarType ModOp::fwd_centered(const ScalarType& x1, const ScalarType& p)
  {
    if(centered_form_not_available_for_args(x1,p))
      return fwd_natural(x1,p);

    return {
      fwd(x1.m,p.m),
      fwd(x1.a,p.a),
      IntervalMatrix(0,0), // not supported yet for auto diff
      x1.def_domain && p.def_domain
    };
  }

  inline void ModOp::bwd(const Interval& y, Interval& x1, Interval& p)
  {
    Interval y_(y);
    ModOp::fwd_bwd(y_,x1,p);
  }

  inline void ModOp::fwd_bwd(Interval& y, Interval& x1, double p) // y = x1 mod(p)
  {
    // The content of this function comes from the IBEX library.
    // See ibex::Interval (IBEX lib, main author: Gilles Chabert)
    //   https://ibex-lib.readthedocs.io

    assert_release(p > 0. && "Modulo needs a strictly positive period p.");

    if(!(x1.diam() > p || y.diam() > p))
    {
      Interval r = (y-x1)/p;
      Interval ir = integer(r);

      if(ir.is_empty()) // additional protection because an empty interval is considered degenerated.
      {
        y.set_empty();
        x1.set_empty();
      }

      else
      {
        if(ir.is_degenerated())
          SubOp::bwd(ir*p,y,x1);

        else if(ir.diam() == 1.)
        {
          Interval y_1 = y; Interval y_2 = y;
          Interval x1_1 = x1; Interval x1_2 = x1;
          SubOp::bwd(Interval(ir.lb()*p),y_1,x1_1);
          SubOp::bwd(Interval(ir.ub()*p),y_2,x1_2);
          y = y_1 | y_2;
          x1 = x1_1 | x1_2;
        }

        else
        {
          assert_release_constexpr(false && "Modulo diameter error.");
        }
      }
    }
  }

  inline void ModOp::fwd_bwd(Interval& y, Interval& x1, Interval& p) // x = y mod(p)
  {
    assert_release(p.is_degenerated() && "ModOp::fwd_bwd(y,x1,p) (with y and x1 intervals) not implemented yet");
    ModOp::fwd_bwd(y, x1, p.mid());
  }
}
