/** 
 *  \file codac2_analytic_flat_input_layout.h
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Simon Rohou
 *  \copyright  Copyright 2026 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <typeindex>
#include <unordered_map>

#include "codac2_AnalyticType.h"
#include "codac2_FunctionArgsList.h"
#include "codac2_ExprBase.h"

namespace codac2
{
  class ScalarVar;

  /**
   * \brief Binding information associated with one input argument in a flattened input domain.
   *
   * An analytic function input argument may be scalar, vectorial or matricial.
   * In a flattened input domain, each argument is represented by a contiguous
   * block of scalar inputs.
   *
   * The type field stores the analytic category of the argument (ScalarType, VectorType, MatrixType)
   */
  struct FlatInputBinding
  {
    std::type_index type = typeid(ScalarType); //!< Analytic type of the bound input argument.
    Index offset = 0;                          //!< First scalar index of the argument in the flattened domain.
    Index rows = 1;                            //!< Number of rows of the argument block.
    Index cols = 1;                            //!< Number of columns of the argument block.

    /**
     * \brief Tests whether the binding corresponds to a scalar input argument.
     *
     * \return True if the binding corresponds to a scalar input argument.
     */
    bool is_scalar() const;

    /**
     * \brief Tests whether the binding corresponds to a vector input argument.
     *
     * \return True if the binding corresponds to a vector input argument.
     */
    bool is_vector() const;

    /**
     * \brief Tests whether the binding corresponds to a matrix input argument.
     *
     * \return True if the binding corresponds to a matrix input argument.
     */
    bool is_matrix() const;
  };

  /**
   * \brief Flattened layout associated with an analytic function input domain.
   *
   * This class provides a canonical flattened representation of a function
   * argument list. Each scalar, vector or matrix input argument is assigned
   * a contiguous block of scalar indices.
   */
  class FlatInputLayout
  {
    public:

      /**
       * \brief Builds the flattened layout associated with a function argument list.
       *
       * \param args Function argument list.
       */
      explicit FlatInputLayout(const FunctionArgsList& args);

      /**
       * \brief Returns the total number of scalar inputs in the flattened domain.
       *
       * \return Flattened input domain size.
       */
      Index size() const;

      /**
       * \brief Tests whether two layouts describe the same flattened input domain.
       *
       * Here, only the total number of flattened scalar inputs matters.
       *
       * \param other Other flattened layout.
       * \return True if both layouts describe the same flattened input domain.
       */
      bool same_domain_as(const FlatInputLayout& other) const;

      /**
       * \brief Tests whether this layout matches the flattened domain induced by a function argument list.
       *
       * \param other_args Function argument list.
       * \return True if both domains are identical once flattened.
       */
      bool same_domain_as(const FunctionArgsList& other_args) const;

      /**
       * \brief Returns the flat index associated with a scalar input variable.
       *
       * This method is restricted to scalar variables appearing directly in the
       * function argument list.
       *
       * \param x Scalar input variable.
       * \return Flat index associated with \p x.
       */
      Index flat_index_of(const ScalarVar& x) const;

      /**
       * \brief Returns the binding associated with an input expression identifier.
       *
       * \param id Input expression identifier.
       * \return Binding associated with \p id.
       */
      const FlatInputBinding& binding_of(const ExprID& id) const;

    private:
    
      std::unordered_map<Index,FlatInputBinding> _bindings; //!< Bindings indexed by input expression identifier.
      Index _size = 0;                                      //!< Total number of scalar inputs in the flattened domain.
  };
}