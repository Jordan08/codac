/**
 * \file codac2_py_Affine.h
 *
 * Python bindings for the public CODAC affine arithmetic API.
 */

#pragma once

#include <pybind11/pybind11.h>

namespace codac2 {

/** \brief Export the public Affine type and affine arithmetic functions. */
void export_Affine(pybind11::module_& m);

/** \brief Export the public AffineVector type. */
void export_AffineVector(pybind11::module_& m);

/** \brief Export the public AffineRow type. */
void export_AffineRow(pybind11::module_& m);

/** \brief Export the public AffineMatrix type. */
void export_AffineMatrix(pybind11::module_& m);

/** \brief Export the public AffineVariables type. */
void export_AffineVariables(pybind11::module_& m);

/**
 * \brief Export all public affine-domain bindings.
 */
inline void export_affine(pybind11::module_& m)
{
  export_Affine(m);
  export_AffineVariables(m);
  export_AffineVector(m);
  export_AffineRow(m);
  export_AffineMatrix(m);
}

} // namespace codac2
