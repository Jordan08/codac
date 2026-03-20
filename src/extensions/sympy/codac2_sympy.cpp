/**
 *  codac2_sympy.cpp
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Simon Rohou
 *  \copyright  Copyright 2026 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include "codac2_sympy.h"
#include "codac2_sympy_bridge.h"

namespace codac2
{
  namespace
  {
    using symbolic::detail::FlatSymbolTable;
    using symbolic::detail::ScalarBridgeContext;
    using symbolic::detail::SympyTransform;

    bool same_flattened_domain(const FunctionArgsList& args_f, const FunctionArgsList& args_g)
    {
      return args_f.total_size() == args_g.total_size();
    }

    Index flat_input_index_of_scalar_var(const FunctionArgsList& args, const ScalarVar& x, const char* ctx)
    {
      Index flat = 0;
      const auto xid = x.unique_id();

      for(const auto& arg : args)
      {
        if(auto s = std::dynamic_pointer_cast<ScalarVar>(arg))
        {
          if(s->unique_id() == xid)
            return flat;
          ++flat;
        }

        else if(auto v = std::dynamic_pointer_cast<VectorVar>(arg))
          flat += v->size();

        else if(auto m = std::dynamic_pointer_cast<MatrixVar>(arg))
          flat += m->rows() * m->cols();

        else
          throw SymbolicDiffError(std::string(ctx) + ": unsupported variable type in argument list");
      }

      throw SymbolicDiffError(std::string(ctx) + ": scalar variable is not part of the function arguments");
    }

    pybind11::object remap_to_reference_symbols(
      pybind11::object expr,
      const FlatSymbolTable& source_symbols,
      const FlatSymbolTable& reference_symbols)
    {
      if(source_symbols.size() != reference_symbols.size())
        throw SymbolicDiffError("sympy_equal: inconsistent flattened symbol tables");

      pybind11::dict subs;
      for(Index k = 0 ; k < source_symbols.size() ; ++k)
        subs[source_symbols.by_flat_index(k)] = reference_symbols.by_flat_index(k);

      return expr.attr("subs")(subs);
    }

    bool sympy_zero(const pybind11::object& expr)
    {
      const pybind11::object& sympy = symbolic::detail::import_sympy();
      pybind11::object simplified = sympy.attr("simplify")(expr);

      try
      {
        simplified = sympy.attr("nsimplify")(simplified, pybind11::arg("rational")=true);
      }
      catch(const pybind11::error_already_set&)
      {
        // Keep the simplified expression if SymPy cannot rationalize it.
      }

      pybind11::object is_zero = simplified.attr("is_zero");
      if(!is_zero.is_none())
        return pybind11::cast<bool>(is_zero);

      pybind11::object equals_zero = simplified.attr("equals")(pybind11::int_(0));
      if(!equals_zero.is_none())
        return pybind11::cast<bool>(equals_zero);

      return false;
    }

    bool sympy_equal_scalar_expr(
      const FunctionArgsList& args_f,
      const ScalarExpr& yf,
      const FunctionArgsList& args_g,
      const ScalarExpr& yg)
    {
      if(!same_flattened_domain(args_f, args_g))
        return false;

      symbolic::detail::ensure_python_runtime();
      pybind11::gil_scoped_acquire gil;
      ScalarBridgeContext reference_ctx(args_f);
      ScalarBridgeContext source_ctx(args_g);

      pybind11::object ef = reference_ctx.export_scalar(yf);
      pybind11::object eg = source_ctx.export_scalar(yg);
      eg = remap_to_reference_symbols(eg, source_ctx.symbols(), reference_ctx.symbols());

      return sympy_zero(ef - eg);
    }

    ScalarExpr sympy_partial_diff_expr(
      const ScalarBridgeContext& ctx,
      const ScalarExpr& y,
      Index flat_input_index)
    {
      if(flat_input_index < 0 || flat_input_index >= ctx.symbols().size())
        throw SymbolicDiffError("sympy_partial_diff_expr: flat_input_index out of bounds");

      return ctx.transform_scalar_expr(
        y,
        [flat_input_index](const pybind11::object& sympy, const pybind11::object& ys, const FlatSymbolTable& symbols)
        {
          return sympy.attr("diff")(ys, symbols.by_flat_index(flat_input_index));
        });
    }

    ScalarExpr sympy_partial_diff_expr(
      const FunctionArgsList& args,
      const ScalarExpr& y,
      Index flat_input_index)
    {
      ScalarBridgeContext ctx(args);
      return sympy_partial_diff_expr(ctx, y, flat_input_index);
    }

    AnalyticFunction<VectorType>
    transform_vector_function(
      const AnalyticFunction<VectorType>& f,
      const SympyTransform& transform,
      bool do_expand = true)
    {
      const auto shape = f.output_shape();
      if(shape.second != 1)
        throw SymbolicDiffError("Only column-vector outputs are supported");

      ScalarBridgeContext ctx(f.args());
      VectorExpr y(f.expr());
      std::vector<ScalarExpr> entries;
      entries.reserve(static_cast<std::size_t>(shape.first));

      for(Index i = 0 ; i < shape.first ; ++i)
        entries.push_back(ctx.transform_scalar_expr(y[i], transform, do_expand));

      return AnalyticFunction<VectorType>(f.args(), vec(entries));
    }

    AnalyticFunction<MatrixType>
    transform_matrix_function(
      const AnalyticFunction<MatrixType>& f,
      const SympyTransform& transform,
      bool do_expand = true)
    {
      const auto shape = f.output_shape();
      ScalarBridgeContext ctx(f.args());
      MatrixExpr y(f.expr());
      std::vector<VectorExpr> cols;
      cols.reserve(static_cast<std::size_t>(shape.second));

      for(Index j = 0 ; j < shape.second ; ++j)
      {
        std::vector<ScalarExpr> col_entries;
        col_entries.reserve(static_cast<std::size_t>(shape.first));
        for(Index i = 0 ; i < shape.first ; ++i)
          col_entries.push_back(ctx.transform_scalar_expr(y(i,j), transform, do_expand));
        cols.push_back(vec(col_entries));
      }

      return AnalyticFunction<MatrixType>(f.args(), mat(cols));
    }
  }

  AnalyticFunction<ScalarType>
  sympy_simplify(const AnalyticFunction<ScalarType>& f)
  {
    return AnalyticFunction<ScalarType>(
      f.args(),
      symbolic::detail::transform_scalar_expr(
        f.args(), ScalarExpr(f.expr()),
        [](const pybind11::object& sympy, const pybind11::object& ys, const FlatSymbolTable&)
        {
          return sympy.attr("simplify")(ys);
        }));
  }

  AnalyticFunction<VectorType>
  sympy_simplify(const AnalyticFunction<VectorType>& f)
  {
    return transform_vector_function(
      f,
      [](const pybind11::object& sympy, const pybind11::object& ys, const FlatSymbolTable&)
      {
        return sympy.attr("simplify")(ys);
      });
  }

  AnalyticFunction<MatrixType>
  sympy_simplify(const AnalyticFunction<MatrixType>& f)
  {
    return transform_matrix_function(
      f,
      [](const pybind11::object& sympy, const pybind11::object& ys, const FlatSymbolTable&)
      {
        return sympy.attr("simplify")(ys);
      });
  }

  AnalyticFunction<ScalarType>
  sympy_horner(const AnalyticFunction<ScalarType>& f)
  {
    return AnalyticFunction<ScalarType>(
      f.args(),
      symbolic::detail::transform_scalar_expr(
        f.args(), ScalarExpr(f.expr()),
        [](const pybind11::object&, const pybind11::object& ys, const FlatSymbolTable&)
        {
          return symbolic::detail::import_polyfuncs().attr("horner")(ys);
        },
        false));
  }

  AnalyticFunction<VectorType>
  sympy_horner(const AnalyticFunction<VectorType>& f)
  {
    return transform_vector_function(
      f,
      [](const pybind11::object&, const pybind11::object& ys, const FlatSymbolTable&)
      {
        return symbolic::detail::import_polyfuncs().attr("horner")(ys);
      },
      false);
  }

  AnalyticFunction<MatrixType>
  sympy_horner(const AnalyticFunction<MatrixType>& f)
  {
    return transform_matrix_function(
      f,
      [](const pybind11::object&, const pybind11::object& ys, const FlatSymbolTable&)
      {
        return symbolic::detail::import_polyfuncs().attr("horner")(ys);
      },
      false);
  }

  AnalyticFunction<ScalarType>
  sympy_partial_diff(const AnalyticFunction<ScalarType>& f, Index flat_input_index)
  {
    if(flat_input_index < 0 || flat_input_index >= f.input_size())
      throw SymbolicDiffError("sympy_partial_diff: flat_input_index out of bounds");

    ScalarBridgeContext ctx(f.args());
    return AnalyticFunction<ScalarType>(
      f.args(),
      sympy_partial_diff_expr(ctx, ScalarExpr(f.expr()), flat_input_index));
  }

  AnalyticFunction<ScalarType>
  sympy_partial_diff(const AnalyticFunction<ScalarType>& f, const ScalarVar& x)
  {
    return sympy_partial_diff(f, flat_input_index_of_scalar_var(f.args(), x, "sympy_partial_diff"));
  }

  AnalyticFunction<ScalarType>
  sympy_diff(const AnalyticFunction<ScalarType>& f)
  {
    if(f.input_size() != 1)
      throw SymbolicDiffError(
        "sympy_diff(f): this overload supports only scalar functions with one flattened input. "
        "Use sympy_partial_diff(f,j), sympy_diff(f,order), or sympy_gradient(f) for multivariate scalar functions.");

    return sympy_partial_diff(f, 0);
  }

  AnalyticFunction<ScalarType>
  sympy_diff(const AnalyticFunction<ScalarType>& f, const ScalarVar& x)
  {
    return sympy_partial_diff(f, x);
  }

  AnalyticFunction<ScalarType>
  sympy_diff(const AnalyticFunction<ScalarType>& f, Index order)
  {
    if(order < 0)
      throw SymbolicDiffError("sympy_diff(f,order): order must be nonnegative");

    if(order == 0)
      return f;

    if(f.input_size() != 1)
      throw SymbolicDiffError("sympy_diff(f,order): this overload supports only scalar functions with one flattened input");

    ScalarBridgeContext ctx(f.args());
    ScalarExpr y(f.expr());
    for(Index k = 0 ; k < order ; ++k)
      y = sympy_partial_diff_expr(ctx, y, 0);

    return AnalyticFunction<ScalarType>(f.args(), y);
  }

  AnalyticFunction<ScalarType>
  sympy_diff(const AnalyticFunction<ScalarType>& f, const ScalarVar& x, Index order)
  {
    if(order < 0)
      throw SymbolicDiffError("sympy_diff(f,x,order): order must be nonnegative");

    if(order == 0)
      return f;

    const Index flat_input_index = flat_input_index_of_scalar_var(f.args(), x, "sympy_diff");

    ScalarBridgeContext ctx(f.args());
    ScalarExpr y(f.expr());
    for(Index k = 0 ; k < order ; ++k)
      y = sympy_partial_diff_expr(ctx, y, flat_input_index);

    return AnalyticFunction<ScalarType>(f.args(), y);
  }

  AnalyticFunction<VectorType>
  sympy_gradient(const AnalyticFunction<ScalarType>& f)
  {
    ScalarBridgeContext ctx(f.args());
    std::vector<ScalarExpr> entries;
    entries.reserve(static_cast<std::size_t>(f.input_size()));
    const ScalarExpr y(f.expr());

    for(Index j = 0 ; j < f.input_size() ; ++j)
      entries.push_back(sympy_partial_diff_expr(ctx, y, j));

    return AnalyticFunction<VectorType>(f.args(), vec(entries));
  }

  AnalyticFunction<MatrixType>
  sympy_hessian(const AnalyticFunction<ScalarType>& f)
  {
    ScalarBridgeContext ctx(f.args());
    std::vector<VectorExpr> cols;
    cols.reserve(static_cast<std::size_t>(f.input_size()));
    const ScalarExpr y(f.expr());

    for(Index j = 0 ; j < f.input_size() ; ++j)
    {
      const ScalarExpr dj = sympy_partial_diff_expr(ctx, y, j);
      std::vector<ScalarExpr> col_entries;
      col_entries.reserve(static_cast<std::size_t>(f.input_size()));
      for(Index i = 0 ; i < f.input_size() ; ++i)
        col_entries.push_back(sympy_partial_diff_expr(ctx, dj, i));
      cols.push_back(vec(col_entries));
    }

    return AnalyticFunction<MatrixType>(f.args(), mat(cols));
  }

  AnalyticFunction<MatrixType>
  sympy_diff(const AnalyticFunction<VectorType>& f)
  {
    const auto shape = f.output_shape();
    if(shape.second != 1)
      throw SymbolicDiffError("sympy_diff(VectorType): only column-vector outputs are supported");

    ScalarBridgeContext ctx(f.args());
    std::vector<VectorExpr> cols;
    cols.reserve(static_cast<std::size_t>(f.input_size()));
    VectorExpr y(f.expr());

    for(Index j = 0 ; j < f.input_size() ; ++j)
    {
      std::vector<ScalarExpr> col_entries;
      col_entries.reserve(static_cast<std::size_t>(shape.first));
      for(Index i = 0 ; i < shape.first ; ++i)
        col_entries.push_back(sympy_partial_diff_expr(ctx, y[i], j));
      cols.push_back(vec(col_entries));
    }

    return AnalyticFunction<MatrixType>(f.args(), mat(cols));
  }

  AnalyticFunction<ScalarType>
  sympy_series(const AnalyticFunction<ScalarType>& f, double center, Index order)
  {
    if(order < 0)
      throw SymbolicDiffError("sympy_series: order must be nonnegative");

    if(f.input_size() != 1)
      throw SymbolicDiffError("sympy_series: this overload supports only scalar functions with one flattened input");

    return AnalyticFunction<ScalarType>(
      f.args(),
      symbolic::detail::transform_scalar_expr(
        f.args(), ScalarExpr(f.expr()),
        [center, order](const pybind11::object& sympy, const pybind11::object& ys, const FlatSymbolTable& symbols)
        {
          pybind11::object x = symbols.by_flat_index(0);
          pybind11::object out = sympy.attr("series")(ys, x, pybind11::float_(center), pybind11::int_(order+1));
          return out.attr("removeO")();
        }));
  }

  AnalyticFunction<ScalarType>
  sympy_series(const AnalyticFunction<ScalarType>& f, const ScalarVar& x, double center, Index order)
  {
    if(order < 0)
      throw SymbolicDiffError("sympy_series: order must be nonnegative");

    const Index flat_input_index = flat_input_index_of_scalar_var(f.args(), x, "sympy_series");

    return AnalyticFunction<ScalarType>(
      f.args(),
      symbolic::detail::transform_scalar_expr(
        f.args(), ScalarExpr(f.expr()),
        [center, order, flat_input_index](const pybind11::object& sympy, const pybind11::object& ys, const FlatSymbolTable& symbols)
        {
          pybind11::object xj = symbols.by_flat_index(flat_input_index);
          pybind11::object s = sympy.attr("series")(ys, xj, pybind11::float_(center), pybind11::int_(order+1));
          pybind11::object poly = s.attr("removeO")();
          return sympy.attr("expand")(poly);
        }));
  }

  bool
  sympy_equal(const AnalyticFunction<ScalarType>& f, const AnalyticFunction<ScalarType>& g)
  {
    return sympy_equal_scalar_expr(f.args(), ScalarExpr(f.expr()), g.args(), ScalarExpr(g.expr()));
  }

  bool
  sympy_equal(const AnalyticFunction<VectorType>& f, const AnalyticFunction<VectorType>& g)
  {
    if(!same_flattened_domain(f.args(), g.args()))
      return false;

    const auto shape_f = f.output_shape();
    const auto shape_g = g.output_shape();
    if(shape_f != shape_g)
      return false;

    VectorExpr yf(f.expr());
    VectorExpr yg(g.expr());

    for(Index i = 0 ; i < shape_f.first ; ++i)
    {
      if(!sympy_equal_scalar_expr(f.args(), yf[i], g.args(), yg[i]))
        return false;
    }

    return true;
  }

  bool
  sympy_equal(const AnalyticFunction<MatrixType>& f, const AnalyticFunction<MatrixType>& g)
  {
    if(!same_flattened_domain(f.args(), g.args()))
      return false;

    const auto shape_f = f.output_shape();
    const auto shape_g = g.output_shape();
    if(shape_f != shape_g)
      return false;

    MatrixExpr yf(f.expr());
    MatrixExpr yg(g.expr());

    for(Index j = 0 ; j < shape_f.second ; ++j)
    {
      for(Index i = 0 ; i < shape_f.first ; ++i)
      {
        if(!sympy_equal_scalar_expr(f.args(), yf(i,j), g.args(), yg(i,j)))
          return false;
      }
    }

    return true;
  }
}
