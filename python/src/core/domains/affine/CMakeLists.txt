# Python bindings for CODAC affine arithmetic.
#
# This directory follows the same role as python/src/core/domains/interval.
# The sources are added to the existing _core target.

target_sources(_core PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/codac2_py_Affine.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/codac2_py_AffineVector.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/codac2_py_AffineRow.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/codac2_py_AffineMatrix.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/codac2_py_AffineVariables.cpp
)
