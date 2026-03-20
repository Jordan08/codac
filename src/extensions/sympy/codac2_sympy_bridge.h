/** 
 *  \file codac2_sympy_bridge.h
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Simon Rohou
 *  \copyright  Copyright 2026 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <functional>
#include <pybind11/embed.h>
#include "codac2_sympy.h"

namespace codac2::symbolic::detail
{
  struct FlatInputBinding
  {
    enum class Kind
    {
      Scalar,
      Vector,
      Matrix
    };

    Kind kind = Kind::Scalar;
    Index offset = 0;
    Index rows = 1;
    Index cols = 1;
  };

  class FlatSymbolTable
  {
  public:
    explicit FlatSymbolTable(const FunctionArgsList& args);

    pybind11::object by_flat_index(Index k) const;
    ScalarExpr codac_expr_by_name(const std::string& name) const;
    Index size() const;
    pybind11::object for_scalar_var(const ScalarVar& x) const;
    pybind11::object for_vector_component(const VectorVar& x, Index i) const;
    pybind11::object for_matrix_component(const MatrixVar& x, Index i, Index j) const;

  private:
    static std::string make_symbol_name(Index flat_index);
    const FlatInputBinding& binding_of(const ExprID& id) const;

  private:
    std::vector<std::string> _names;
    std::unordered_map<std::string, ScalarExpr> _codac_scalars;
    std::unordered_map<Index, FlatInputBinding> _bindings;
  };

  class SympyExporter
  {
  public:
    explicit SympyExporter(const FlatSymbolTable& symbols);
    pybind11::object export_scalar(const ScalarExpr& e) const;

  private:
    pybind11::object export_node(const std::shared_ptr<ExprBase>& e) const;
    pybind11::object export_vector_component(const std::shared_ptr<ExprBase>& e, Index i) const;
    pybind11::object export_matrix_component(const std::shared_ptr<ExprBase>& e, Index i, Index j) const;
    static double scalar_const_value(const ConstValueExpr<ScalarType>& c);
    static std::shared_ptr<ExprBase> child_at(
      const std::vector<std::shared_ptr<ExprBase>>& children, Index i, const char* ctx);

  private:
    const FlatSymbolTable& _symbols;
  };

  class SympyImporter
  {
  public:
    explicit SympyImporter(const FlatSymbolTable& symbols);
    ScalarExpr import_scalar(const pybind11::handle& obj) const;

  private:
    ScalarExpr import_add(const pybind11::handle& obj) const;
    ScalarExpr import_mul(const pybind11::handle& obj) const;
    ScalarExpr import_pow(const pybind11::handle& obj) const;
    ScalarExpr import_function(const pybind11::handle& obj) const;

  private:
    const FlatSymbolTable& _symbols;
  };

  using SympyTransform = std::function<pybind11::object(
    const pybind11::object&, const pybind11::object&, const FlatSymbolTable&)>;

  class ScalarBridgeContext
  {
  public:
    explicit ScalarBridgeContext(const FunctionArgsList& args);

    const pybind11::object& sympy() const;
    const FlatSymbolTable& symbols() const;

    pybind11::object export_scalar(const ScalarExpr& y) const;
    ScalarExpr import_scalar(const pybind11::handle& obj) const;
    ScalarExpr transform_scalar_expr(const ScalarExpr& y, const SympyTransform& transform, bool do_expand = true) const;

  private:
    FlatSymbolTable _symbols;
    SympyExporter _exporter;
    SympyImporter _importer;
  };

  void ensure_python_runtime();
  const pybind11::object& import_sympy();
  const pybind11::object& import_polyfuncs();
  const pybind11::object& import_builtins();
  pybind11::object normalize_sympy_expr(const pybind11::object& sympy, pybind11::object expr, bool do_expand = true);

  ScalarExpr transform_scalar_expr(
    const FunctionArgsList& args,
    const ScalarExpr& y,
    const SympyTransform& transform,
    bool do_expand = true);
}
