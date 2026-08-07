/**
 *  \file codac2_AffineVector.h
 *
 *  Dynamic-size vector wrapper for affine forms.
 *
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Jordan Ninin
 *  \copyright  Copyright 2026 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <ostream>
#include "codac2_matrices.h"
#include "codac2_Vector.h"
#include "codac2_IntervalVector.h"
#include "codac2_AffineMain.h"

namespace codac2 {

template<class T>  class AffineVarMainVector;
/**
 * \ingroup arithmetic
 *
 * \brief Vector of Affine form
 *
 * By convention an empty vector has a dimension. A vector becomes empty
 * when one of its component becomes empty and all the components
 * are set to the empty Interval.
 */
/**
 * \brief Alias for a dynamic-size column vector of intervals.
 *
 * Represents a column vector with a dynamic number of rows,
 * where each element is an Interval object.
 *
 * This type alias is based on Eigen's matrix template and corresponds to ``Eigen::Matrix<Interval,-1,1>``.
 */



template<class T>
class AffineMainVector : public Eigen::Matrix<AffineMain<T>, -1, 1> {

public:

	using Base = Eigen::Matrix<AffineMain<T>, -1, 1>;

	// Inherit the constructors and the assignment operators of the Eigen base
	// class (dynamic size constructor, expression assignment, etc.).
	// Without these "using" declarations, the derived class only gets the
	// implicitly-generated default/copy/move constructors and assignment
	// operators, and every other Eigen::Matrix constructor/assignment
	// (e.g. AffineMainVector<T>(n) or v = v1+v2) would fail to compile.
	using Base::Base;
	using Base::operator=;

	AffineMainVector(const AffineVarMainVector<T>& x) : Base(x.size(), 1)
	{
		for(Index i = 0; i < x.size(); ++i)
			(*this)[i] = (x[i]);
	}

	AffineMainVector& operator=(const AffineVarMainVector<T>& x)
	{
		if (this != &x) {
			this->resize(x.size());
			for(Index i = 0; i < x.size(); ++i)
				(*this)[i] = (x[i]);
		}
		return *this;
	}
	/**
	 * \brief Returns the interval hull of each affine component.
	 *
	 * Build an \c IntervalVector with the same size as \c *this, where
	 * the i-th entry is \c (*this)[i].itv().
	 * \pre (*this) must be nonempty
	 */
	IntervalVector itv() const;

	void set_empty();

protected:
	friend class Eigen::Matrix< AffineMain<T>,-1,-1>;
	/**
	 * \brief Creates a vector of size \p n.
	 *
	 * Each component is set to the default (unbounded) affine form, i.e.
	 * the equivalent of ``Interval()`` (]-oo,+oo[).
	 *
	 * \note This constructor hides (replaces) the generic
	 * ``Eigen::Matrix(int n)`` constructor inherited above via
	 * ``using Base::Base``. That generic constructor always tries to
	 * zero-initialize non-interval scalar types through ``init(0.)``,
	 * which does not compile for ``AffineMain<T>`` since there is no
	 * implicit conversion from ``double`` to ``AffineMain<T>``.
	 *
	 * \param n vector size
	 */
	explicit AffineMainVector(int n) : Base(n,1)
	{ }
};




template<class T>
IntervalVector AffineMainVector<T>::itv() const {
	IntervalVector tmp(this->size());
	for (int i = 0; i < this->size(); i++) {
		tmp[i] = (*this)[i].itv();
	}
	return tmp;
}

template<class T>
void AffineMainVector<T>::set_empty() {
	for (int i = 0; i < this->size(); i++) {
		((*this)[i]).set_empty();
	}
}

/**
 * \brief Stream output operator for ``IntervalVector`` objects.
 *
 * \param os The output stream to write to.
 * \param x The interval vector whose contents are to be printed.
 * \return A reference to the modified output stream.
 */

template<class T>
inline std::ostream& operator<<(std::ostream& os, const AffineMainVector<T>& x)
{
	if(x.is_empty())
		return os << "[ empty " << x.size() << "d box ]";

	else
	{
		os << x.format(codac_vector_fmt());
		return os;
	}

}

} // namespace codac2
