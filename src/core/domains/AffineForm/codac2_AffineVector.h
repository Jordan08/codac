/**
 *  \file codac2_AffineVector.h
 *
 *  Dynamic-size vector wrapper for affine forms.
 *
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Jordan Ninin
 *  \copyright  Copyright 2020 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <ostream>
#include "codac2_Vector.h"

namespace codac2 {

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


template<class T=AF_Default>
class AffineMainVector : public Eigen::Matrix<AffineMain<T>, -1, 1> {

public:
	/**
	 * \brief Returns the interval hull of each affine component.
	 *
	 * Build an \c IntervalVector with the same size as \c *this, where
	 * the i-th entry is \c (*this)[i].itv().
	 * \pre (*this) must be nonempty
	 */
	IntervalVector itv() const;

};



typedef AffineMainVector<AF_Default> Affine2Vector;
//typedef AffineMainVector<AF_Other> Affine3Vector;




template<class T>
IntervalVector AffineMainVector<T>::itv() const {
	IntervalVector tmp(this->size());
	for (int i = 0; i < this->size(); i++) {
		tmp[i] = (*this)[i].itv();
	}
	return tmp;
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
