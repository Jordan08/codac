# ==================================================================
#  Codac - cmake configuration file
# ==================================================================

# Finding Catch2

option(ENABLE_FIND_PACKAGE_CATCH2 "ENABLE_FIND_PACKAGE_CATCH2" ON)
if(ENABLE_FIND_PACKAGE_CATCH2)
  find_package(Catch2 3.0.0 QUIET) # Use Catch2_DIR or CMAKE_PREFIX_PATH to specify a custom Catch2 installation path...
endif()
if(Catch2_FOUND)
  message(STATUS "Found Catch2 version ${Catch2_VERSION} in ${Catch2_DIR}")
else()
  message(STATUS "Attempt to install Catch2 automatically...")
  FetchContent_Declare(Catch2 URL https://github.com/catchorg/Catch2/archive/refs/tags/v3.6.0.zip)
  #FetchContent_Declare(Catch2 URL ${CMAKE_CURRENT_SOURCE_DIR}/../3rd/Catch2-3.6.0.zip) # If needed to be self-contained...
  FetchContent_MakeAvailable(Catch2)
endif()
# Adds Catch2::Catch2WithMain


# Test sources
# ==================================================================

list(APPEND SRC_TESTS

  # ----------------------------------------------------------------
  # Core / Affine
  # ----------------------------------------------------------------

  core/domains/affine/codac2_tests_AffineForm_base
  core/domains/affine/codac2_tests_AffineVector
  core/domains/affine/codac2_tests_AffineMatrix_base
  core/domains/affine/codac2_tests_AffineMatrix_product
  core/domains/affine/codac2_tests_AffineForm_Eigen
  core/domains/affine/codac2_tests_AffineForm_coverage

  ../doc/manual/manual/affine/src


  # ----------------------------------------------------------------
  # 3rd
  # ----------------------------------------------------------------

  core/3rd/codac2_tests_eigen


  # ----------------------------------------------------------------
  # Actions
  # ----------------------------------------------------------------

  core/actions/codac2_tests_OctaSym


  # ----------------------------------------------------------------
  # Contractors
  # ----------------------------------------------------------------

  core/contractors/codac2_tests_CtcAction
  core/contractors/codac2_tests_CtcCartProd
  core/contractors/codac2_tests_CtcCtcBoundary
  core/contractors/codac2_tests_CtcDeriv
  core/contractors/codac2_tests_CtcEval
  core/contractors/codac2_tests_CtcFixpoint
  core/contractors/codac2_tests_CtcInter
  core/contractors/codac2_tests_CtcInverse
  core/contractors/codac2_tests_CtcInverseNotIn
  core/contractors/codac2_tests_CtcLazy
  core/contractors/codac2_tests_CtcLohner
  core/contractors/codac2_tests_CtcPolygon
  core/contractors/codac2_tests_CtcSegment
  core/contractors/codac2_tests_CtcVisible
  core/contractors/codac2_tests_linear_ctc

  ../doc/manual/manual/contractors/geometric/src
  ../doc/manual/manual/contractors/analytic/src
  ../doc/manual/manual/contractors/set/src


  # ----------------------------------------------------------------
  # Domains
  # ----------------------------------------------------------------

  core/domains/codac2_tests_BoolInterval

  core/domains/ellipsoid/codac2_tests_Ellipsoid

  core/domains/interval/codac2_tests_Interval
  core/domains/interval/codac2_tests_Interval_operations
  core/domains/interval/codac2_tests_IntervalMatrix
  core/domains/interval/codac2_tests_IntervalVector

  ../doc/manual/manual/intervals/src

  core/domains/zonotope/codac2_tests_Parallelepiped_eval
  core/domains/zonotope/codac2_tests_Parallelepiped
  core/domains/zonotope/codac2_tests_Zonotope

  core/domains/tube/codac2_tests_TDomain
  core/domains/tube/codac2_tests_Slice
  core/domains/tube/codac2_tests_Slice_polygon
  core/domains/tube/codac2_tests_SlicedTube
  core/domains/tube/codac2_tests_SlicedTube_integral

  ../doc/manual/manual/tubes/src


  # ----------------------------------------------------------------
  # Functions
  # ----------------------------------------------------------------

  core/functions/analytic/codac2_tests_AnalyticFunction

  ../doc/manual/manual/functions/analytic/src


  # ----------------------------------------------------------------
  # Geometry
  # ----------------------------------------------------------------

  core/geometry/codac2_tests_ConvexPolygon
  core/geometry/codac2_tests_geometry
  core/geometry/codac2_tests_Polygon
  core/geometry/codac2_tests_Segment

  ../doc/manual/manual/geometry/src


  # ----------------------------------------------------------------
  # Matrices
  # ----------------------------------------------------------------

  core/matrices/codac2_tests_arithmetic_add
  core/matrices/codac2_tests_arithmetic_div
  core/matrices/codac2_tests_arithmetic_mul
  core/matrices/codac2_tests_arithmetic_sub
  core/matrices/codac2_tests_cart_prod
  core/matrices/codac2_tests_Matrix
  core/matrices/codac2_tests_Vector
  core/matrices/codac2_tests_hull
  core/matrices/codac2_tests_inversion
  core/matrices/codac2_tests_IntFullPivLU
  core/matrices/codac2_tests_GaussJordan

  ../doc/manual/manual/linear/src


  # ----------------------------------------------------------------
  # Operators
  # ----------------------------------------------------------------

  core/operators/codac2_tests_operators


  # ----------------------------------------------------------------
  # Peibos
  # ----------------------------------------------------------------

  core/peibos/codac2_tests_peibos


  # ----------------------------------------------------------------
  # Separators
  # ----------------------------------------------------------------

  core/separators/codac2_tests_SepCartProd
  core/separators/codac2_tests_SepCtcBoundary
  core/separators/codac2_tests_SepInverse
  core/separators/codac2_tests_SepPolygon
  core/separators/codac2_tests_SepProj
  core/separators/codac2_tests_SepTransform
  core/separators/codac2_tests_SepVisible


  # ----------------------------------------------------------------
  # Tools
  # ----------------------------------------------------------------

  core/tools/codac2_tests_Approx
  core/tools/codac2_tests_serialization
  core/tools/codac2_tests_transformations
  core/tools/codac2_tests_trunc
  core/tools/ibex/codac2_tests_ibex

  ../doc/manual/manual/tools/src


  # ----------------------------------------------------------------
  # Trajectory
  # ----------------------------------------------------------------

  core/trajectory/codac2_tests_AnalyticTraj
  core/trajectory/codac2_tests_SampledTraj


  # ----------------------------------------------------------------
  # Graphics
  # ----------------------------------------------------------------

  graphics/styles/codac2_tests_Color
)


# Python files required by tests
# ==================================================================

file(COPY
  ${CMAKE_CURRENT_SOURCE_DIR}/core/domains/tube/codac2_tests_predefined_tubes.py
  # Add here other Python files that are not tests but need to be exported to Python package
  # ...
  
  DESTINATION
  ${CMAKE_BINARY_DIR}/python/python_package/codac/tests
)


# Libraries
# ==================================================================

set(CODAC_LIBRARIES ${PROJECT_NAME}-core ${PROJECT_NAME}-graphics)

# Sympy
# ==================================================================

option(BUILD_SYMPY_EMBED_TESTS "Build C++ tests that require embedded Python for Sympy" ON)
if(WITH_PYTHON AND DEFINED PYBIND11_FINDPYTHON AND NOT PYBIND11_FINDPYTHON)
  message(STATUS "Disabling Sympy C++ embed tests because PYBIND11_FINDPYTHON=OFF")
  set(BUILD_SYMPY_EMBED_TESTS OFF)
endif()

# Sympy test
if (WITH_PYTHON AND BUILD_SYMPY_EMBED_TESTS)
  list(APPEND SRC_TESTS
    extensions/sympy/codac2_tests_sympy
    ../doc/manual/manual/extensions/sympy/src
  )
endif()


# CAPD
# ==================================================================
if (WITH_CAPD)
  list(APPEND SRC_TESTS
    extensions/capd/codac2_tests_capd
    extensions/capd/codac2_tests_peibos_capd
    ../doc/manual/manual/extensions/capd/src
  )
  set(CODAC_LIBRARIES ${CODAC_LIBRARIES} ${PROJECT_NAME}-capd capd::capd)
endif()

# IBEX test
#if (WITH_IBEX)
  #find_package(IBEX REQUIRED)
  #ibex_init_common() # IBEX should have installed this function
  #message(STATUS "Found IBEX version ${IBEX_VERSION}")
#endif()



# Common Codac include directories
# ==================================================================
#
# IMPORTANT:
# Use the source tree as the canonical header tree.
#
# The build/include directory contains symbolic links to these
# headers, generated by src/core/CMakeLists.txt.
#
# Therefore there is only one physical copy of each header.

set(CODAC_CORE_SOURCE_INCLUDE_DIRS

  ${CMAKE_CURRENT_SOURCE_DIR}/../src/core/actions
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/core/contractors

  ${CMAKE_CURRENT_SOURCE_DIR}/../src/core/domains
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/core/domains/ellipsoid
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/core/domains/interval
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/core/domains/interval/eigen
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/core/domains/affine
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/core/domains/paving
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/core/domains/zonotope
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/core/domains/tube

  ${CMAKE_CURRENT_SOURCE_DIR}/../src/core/functions
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/core/functions/analytic
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/core/functions/set

  ${CMAKE_CURRENT_SOURCE_DIR}/../src/core/geometry

  ${CMAKE_CURRENT_SOURCE_DIR}/../src/core/matrices
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/core/matrices/eigen/Matrix_addons
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/core/matrices/eigen/MatrixBase_addons

  ${CMAKE_CURRENT_SOURCE_DIR}/../src/core/operators
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/core/paver
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/core/peibos
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/core/proj
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/core/separators
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/core/tools
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/core/trajectory
)


set(
  CODAC_HEADERS_DIR
  ${CMAKE_CURRENT_BINARY_DIR}/../include
)



# ------------------------------------------------------------------
# Sanitizer runtime preload for Python tests (Debug builds only)
# ------------------------------------------------------------------
#
# In Debug mode, _core.so (the pybind11 module) is compiled with
# -fsanitize=address,undefined -shared-libsan (see src/core/CMakeLists.txt),
# which makes it depend on Clang's shared sanitizer runtime instead of
# embedding it statically. Since python3 itself has no sanitizer runtime,
# that shared runtime must be preloaded before Python imports codac,
# otherwise symbols like __ubsan_vptr_type_cache stay unresolved.

set(PYTHON_TEST_ENV_ARGS "PYTHONPATH=${PYTHON_PACKAGE_DIR}")

if(WITH_PYTHON AND CMAKE_BUILD_TYPE STREQUAL "Debug"
   AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")

  execute_process(
    COMMAND ${CMAKE_CXX_COMPILER} --print-runtime-dir
    OUTPUT_VARIABLE CLANG_RT_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  if(EXISTS "${CLANG_RT_DIR}")
    file(GLOB SANITIZER_RT_LIBS
      "${CLANG_RT_DIR}/libclang_rt.asan-x86_64.so"
      "${CLANG_RT_DIR}/libclang_rt.ubsan_standalone-x86_64.so"
    )
  endif()

  if(SANITIZER_RT_LIBS)
    list(JOIN SANITIZER_RT_LIBS ":" SANITIZER_RT_LIBS_JOINED)
    list(APPEND PYTHON_TEST_ENV_ARGS "LD_PRELOAD=${SANITIZER_RT_LIBS_JOINED}")

    # LeakSanitizer (part of ASan) flags many long-lived, process-lifetime
    # allocations made by the Python interpreter and by pybind11 module
    # init (interned strings, type objects, strdup'd docstrings...) as
    # "leaks", since it can't trace Python's own object graph. This is a
    # well-known false-positive source for embedded Python interpreters,
    # not something codac's code controls. UBSan and ASan's other checks
    # (use-after-free, buffer overflow, UB...) remain fully active.
    list(APPEND PYTHON_TEST_ENV_ARGS "ASAN_OPTIONS=detect_leaks=0")
    list(APPEND PYTHON_TEST_ENV_ARGS "UBSAN_OPTIONS=print_stacktrace=1")
    message(STATUS "Python tests: LD_PRELOAD=${SANITIZER_RT_LIBS_JOINED}")
  else()
    message(WARNING
      "Debug build with WITH_PYTHON: no shared ASan/UBSan runtime found "
      "in '${CLANG_RT_DIR}'. Make sure the _core target is linked with "
      "-shared-libsan, otherwise Python tests will fail with "
      "'undefined symbol' errors.")
  endif()

endif()

# Build every test
# ==================================================================

foreach(SRC_TEST ${SRC_TESTS})
  string(REPLACE "/" "_" TEST_NAME ${SRC_TEST})
  string(REPLACE "codac2_tests_" "" TEST_NAME ${TEST_NAME})
  set(TEST_NAME codac2_tests_${TEST_NAME})


  # ---------------------------------------------------------------
  # C++ test
  # ---------------------------------------------------------------

  add_executable(
    ${TEST_NAME}
    ${CMAKE_CURRENT_SOURCE_DIR}/${SRC_TEST}.cpp
  )


  # Use the same canonical source headers as codac-core.
  #
  # CODAC_HEADERS_DIR is kept as well because it contains generated
  # umbrella headers such as codac-core.h.
  # set(CODAC_HEADERS_DIR ${CMAKE_CURRENT_BINARY_DIR}/../include)

  target_include_directories(
    ${TEST_NAME}
    SYSTEM PUBLIC

      ${CODAC_HEADERS_DIR}

      ${CODAC_CORE_SOURCE_INCLUDE_DIRS}
  )


  target_link_libraries(
    ${TEST_NAME}
    PUBLIC

      Ibex::ibex
      ${CODAC_LIBRARIES}

    PRIVATE

      Catch2::Catch2WithMain
  )


  # Sympy embedded tests
  if(
    WITH_PYTHON
    AND SRC_TEST MATCHES "extensions/sympy/"
  )

    target_link_libraries(
      ${TEST_NAME}
      PRIVATE

        ${PROJECT_NAME}-sympy
        pybind11::embed
    )

  endif()


  # C++20
  set_target_properties(
    ${TEST_NAME}
    PROPERTIES

      CXX_STANDARD 20
      CXX_STANDARD_REQUIRED ON
      CXX_EXTENSIONS OFF
  )


  add_dependencies(
    check
    ${TEST_NAME}
  )

  add_test(
    NAME ${TEST_NAME}_cpp
    COMMAND ${TEST_NAME}
  )

  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set_tests_properties(${TEST_NAME}_cpp PROPERTIES
      ENVIRONMENT "UBSAN_OPTIONS=print_stacktrace=1"
    )
  endif()

  # ---------------------------------------------------------------
  # Python test
  # ---------------------------------------------------------------

  if(WITH_PYTHON AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${SRC_TEST}.py")
    add_test(
      NAME ${TEST_NAME}_py
      COMMAND
        ${CMAKE_COMMAND} -E env
        ${PYTHON_TEST_ENV_ARGS}
        ${PYTHON_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/${SRC_TEST}.py
    )
  endif()

endforeach()
