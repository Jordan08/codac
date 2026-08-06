/**
 *  \file codac2_Affine.h
 *
 *  Aggregate header for affine-domain classes.
 *
 * ----------------------------------------------------------------------------
 *  \date       2026
 *  \author     Jordan Ninin
 *  \copyright  Copyright 2026 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once


#include "codac2_AffineMain.h"
#include "codac2_AffineVector.h"
#include "codac2_AffineMatrix.h"
#include "codac2_AffineVar.h"
#include "codac2_AffineVarVector.h"

namespace codac2 {

	//typedef AF_fAF1  AF_Default;
	typedef AF_fAF2  AF_Default;
	//typedef AF_fAF2_fma  AF_Default;
	//typedef AF_iAF  AF_Default;
	//typedef AF_sAF  AF_Default;
	//typedef AF_No  AF_Default;
	//typedef AF_fAFFullI AF_Default;

	/** \brief Default affine type based on \c AF_Default. */
	typedef AffineMain<AF_Default> 			Affine2;
	typedef AffineVarMainVector<AF_Default>	Affine2Variables;
	typedef AffineMainVector<AF_Default> 	Affine2Vector;
	typedef AffineMainMatrix<AF_Default> 	Affine2Matrix;

} // namespace codac2

