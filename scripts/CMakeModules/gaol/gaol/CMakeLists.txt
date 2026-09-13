# ==================================================================
#  Codac - CMake build of GAOL, for codac_gaol_build()
# ==================================================================
#
# GAOL is the interval arithmetic library written by Frederic Goualard,
# https://github.com/goualard-f/GAOL (GNU LGPL v2). It comes with an autotools
# and a meson build, and with configuration files written by hand for Visual
# C++ and MinGW, but without a CMake build. codac_gaol_build()
# (scripts/CMakeModules/codac_gaol.cmake) downloads it, and
# ../codac_gaol_patch.cmake copies this file, gaol_configuration.h.in and
# codac_msvc/ into it, and makes to GAOL's sources the few changes it lists,
# before configuring it.
#
# Origin
# ------
# This is the CMakeLists.txt with which IBEX builds GAOL
# (https://github.com/ibex-team/ibex-lib, GNU LGPL v3,
# interval_lib_wrapper/gaol/3rd/gaol-4.2.3alpha0/CMakeLists.txt), written by
# Cyril Bouvier and updated by Gilles Chabert, with the list of sources of
# GAOL's meson build (gaol/meson.build). It differs in the following respects:
#
#  - GAOL is configured the way IBEX configures it, with no option to do
#    otherwise: see gaol_configuration.h.in.
#
#  - What IBEX tests on CMAKE_SYSTEM_NAME and CMAKE_SYSTEM_PROCESSOR, or by
#    running a program (the byte order, the SSE2 and SSE3 instructions), is
#    decided by the preprocessor or by compilation checks here, which work the
#    same in a cross-compilation. IBEX stops on a system it does not know, and
#    asks for those three to be given on the command line when
#    cross-compiling.
#
#  - The compilation flags are the ones codac_gaol_build() passes on, which are
#    Codac's, with CMake's Release flags. IBEX adds -O3, -Wall, -Wconversion,
#    -funroll-loops, -fomit-frame-pointer, -fexpensive-optimizations and
#    hidden visibility: CMake's Release flags already optimise, and the
#    visibility of the symbols only matters for a shared libgaol.
#
#  - Visual C++, for which IBEX does not build GAOL, is given what the fork of
#    GAOL by Fabrice Le Bars (https://github.com/lebarsfa/GAOL) gives it as
#    well: a <unistd.h> for the lexer (codac_msvc/unistd.h), and the definition
#    of __GAOL_PUBLIC__ that tells GAOL it is not a DLL (see below).
#
#  - mathlib is not linked to libgaol, which is a static library: Codac links
#    both (Codac::gaol links Codac::ultim).

cmake_minimum_required(VERSION 3.14...3.25)

# The version gaol/meson.build and configure.ac declare, at the commit
# codac_gaol_build() downloads
project(gaol VERSION 4.2.3 LANGUAGES C CXX)

# As gaol/meson.build asks: cpp_std=c++11, c_std=c99
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_C_STANDARD 99)

# GAOL sets the rounding direction with <fenv.h> (gaol/gaol_fpu.h)
include(CheckIncludeFileCXX)
check_include_file_cxx(fenv.h HAVE_FENV_H)
if(NOT HAVE_FENV_H)
  message(FATAL_ERROR "GAOL needs <fenv.h>, which this compiler does not provide")
endif()

# gaol/gaol_profile.cpp measures time with getrusage() where it exists, and with
# clock() otherwise
include(CheckSymbolExists)
check_symbol_exists(getrusage "sys/resource.h" HAVE_GETRUSAGE)
set(GETRUSAGE_IN_HEADER ${HAVE_GETRUSAGE})

# With -frounding-math, GAOL negates a double with the minus sign rather than by
# flipping its sign bit (gaol_opposite() in gaol/gaol_port.h): the compiler then
# keeps the negations by which GAOL rounds downward with the rounding direction
# set upward.
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-frounding-math HAVE_ROUNDING_MATH_OPTION)

# On x86 processors, the intervals of GAOL are computed with SSE2 instructions
# (gaol/gaol_interval_sse.h), except on Windows: IBEX builds GAOL there with
# --disable-simd ("it can generate failures during execution"), and so does the
# fork of Fabrice Le Bars. On a 32-bit processor, the check passes when the flags
# enable SSE2, as Codac's do (codac_gaol_compiler_flags()). GAOL's SSE3
# instructions are left out: they only serve gaol::interval2f, intervals of
# floats, which Codac does not use, and without them the processor needs no
# more than SSE2, which every x86_64 processor has.
include(CheckCXXSourceCompiles)
check_cxx_source_compiles("
  #if !(defined(__x86_64__) || defined(__i386__)) || defined(_WIN32)
  #error GAOL computes its intervals with SSE2 on x86 processors outside Windows only
  #endif
  #include <emmintrin.h>
  int main() { return _mm_movemask_pd(_mm_set1_pd(-1.0)); }"
  USING_SSE2_INSTRUCTIONS)

configure_file(gaol_configuration.h.in "${CMAKE_CURRENT_BINARY_DIR}/gaol/gaol_configuration.h")
# The double operations of the mathematical library GAOL is built with, mathlib
# here, as gaol/meson.build and IBEX name them
configure_file(gaol/gaol_double_op_apmathlib.h "${CMAKE_CURRENT_BINARY_DIR}/gaol/gaol_double_op.h" COPYONLY)

# The sources of libgaol, as listed in gaol/meson.build. gaol_interval.cpp
# includes gaol_interval_sse.cpp or gaol_interval_fpu.cpp itself, according to
# USING_SSE2_INSTRUCTIONS.
#
# gaol_exact.c defines get_inexact() and clear_inexact(), which the lexer uses,
# with <fenv.h> -- or, for Visual C++, with 32-bit x86 assembly
# (gaol/sysdeps/gaol_exact_msvc.h), which Visual C++ compiles for no other
# processor. It is left out for Visual C++, for which codac_gaol_patch.cmake has
# gaol/gaol_exact.h define them with <fenv.h> instead, in the lexer itself.
set(GAOL_SOURCES
  gaol_interval.cpp gaol_profile.cpp gaol_common.cpp gaol_parser.cpp
  gaol_expression.cpp gaol_port.cpp gaol_exceptions.cpp gaol_init_cleanup.cpp
  s_nextafter.c gaol_intervalf.cpp gaol_interval_lexer.cpp
  gaol_interval_parser.cpp
)
if(NOT MSVC)
  list(APPEND GAOL_SOURCES gaol_exact.c)
endif()
list(TRANSFORM GAOL_SOURCES PREPEND "gaol/")

add_library(gaol STATIC ${GAOL_SOURCES})
# The generated headers, then GAOL's own: its sources include "gaol/..."
target_include_directories(gaol PRIVATE "${CMAKE_CURRENT_BINARY_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
if(MSVC)
  target_include_directories(gaol PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/codac_msvc")
  # gaol/gaol_config.h declares GAOL's classes and functions
  # __declspec(dllimport) for Visual C++ unless __GAOL_PUBLIC__ is defined, as
  # for a DLL. Codac defines it in the same way for the code including GAOL's
  # headers (codac_gaol_create_targets() in scripts/CMakeModules/codac_gaol.cmake).
  target_compile_definitions(gaol PRIVATE __GAOL_PUBLIC__=)
endif()

install(TARGETS gaol ARCHIVE DESTINATION lib)
# GAOL's headers, including those gaol/gaol_config.h includes for Visual C++
# and MinGW, which gaol/meson.build does not install
file(GLOB GAOL_HEADERS gaol/*.h)
install(FILES
  gaol/gaol
  ${GAOL_HEADERS}
  "${CMAKE_CURRENT_BINARY_DIR}/gaol/gaol_configuration.h"
  "${CMAKE_CURRENT_BINARY_DIR}/gaol/gaol_double_op.h"
  DESTINATION include/gaol
)
install(FILES gaol/sysdeps/gaol_exact_c99.h gaol/sysdeps/gaol_exact_msvc.h DESTINATION include/gaol/sysdeps)
