/** 
 *  \file codac2_analytic_flat_input_layout.cpp
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Simon Rohou
 *  \copyright  Copyright 2026 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#include "codac2_analytic_flat_input_layout.h"

#include "codac2_assert.h"
#include "codac2_analytic_variables.h"

namespace codac2
{
  bool FlatInputBinding::is_scalar() const
  {
    return type == std::type_index(typeid(ScalarType));
  }

  bool FlatInputBinding::is_vector() const
  {
    return type == std::type_index(typeid(VectorType));
  }

  bool FlatInputBinding::is_matrix() const
  {
    return type == std::type_index(typeid(MatrixType));
  }

  FlatInputLayout::FlatInputLayout(const FunctionArgsList& args)
  {
    Index flat = 0;

    for(const auto& arg : args)
    {
      if(std::dynamic_pointer_cast<ScalarVar>(arg))
      {
        _bindings.emplace(arg->unique_id().id(),
          FlatInputBinding{typeid(ScalarType), flat, 1, 1});
        ++flat;
      }

      else if(auto v = std::dynamic_pointer_cast<VectorVar>(arg))
      {
        _bindings.emplace(arg->unique_id().id(),
          FlatInputBinding{typeid(VectorType), flat, v->size(), 1});
        flat += v->size();
      }

      else if(auto m = std::dynamic_pointer_cast<MatrixVar>(arg))
      {
        _bindings.emplace(arg->unique_id().id(),
          FlatInputBinding{typeid(MatrixType), flat, m->rows(), m->cols()});
        flat += m->rows() * m->cols();
      }

      else
        assert(false && "unsupported variable type in function argument list");
    }

    _size = flat;
  }

  Index FlatInputLayout::size() const
  {
    return _size;
  }

  bool FlatInputLayout::same_domain_as(const FlatInputLayout& other) const
  {
    return size() == other.size();
  }

  bool FlatInputLayout::same_domain_as(const FunctionArgsList& other_args) const
  {
    return same_domain_as(FlatInputLayout(other_args));
  }

  Index FlatInputLayout::flat_index_of(const ScalarVar& x) const
  {
    const auto& b = binding_of(x.unique_id());
    assert(b.is_scalar() && "expected a scalar input variable");
    return b.offset;
  }

  const FlatInputBinding& FlatInputLayout::binding_of(const ExprID& id) const
  {
    auto it = _bindings.find(id.id());
    assert(it != _bindings.end() && "unknown input expression identifier");
    return it->second;
  }
}