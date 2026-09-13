# ==================================================================
#  Codac - prepares the sources of GAOL or mathlib for their CMake build
# ==================================================================
#
# Run by codac_gaol_build() (scripts/CMakeModules/codac_gaol.cmake) as the patch
# step of the download of mathlib or GAOL:
#
#   cmake -DCOMPONENT=mathlib|gaol -DSOURCE_DIR=<sources> -P codac_gaol_patch.cmake
#
# Neither mathlib nor GAOL comes with a CMake build. This script adds the one
# Codac provides for each (the files of mathlib/ or gaol/ next to it, copied to
# the root of the sources), routes Visual C++ and MinGW to the configuration
# header this build generates, and makes to GAOL the changes listed below,
# which Codac depends on or without which GAOL does not compile or compute
# right on some systems. mathlib's sources are not modified.
#
# Each change replaces an exact piece of a file, and fails if that piece is not
# found exactly once (or the expected number of times): mathlib is pinned by
# the checksum of its archive and GAOL by a commit, and should either change,
# the changes have to be looked at again rather than applied blindly. The
# original of each file modified is kept next to it, as <file>.codac-orig, and
# the changes are always made from it, so that the script can run again on
# sources it already modified. That happens when the script, or one of the
# files it copies, changes: codac_gaol_build() passes their checksum on the
# command line, which reruns the patch step of ExternalProject.

cmake_minimum_required(VERSION 3.14)

if(NOT COMPONENT MATCHES "^(mathlib|gaol)$" OR NOT IS_DIRECTORY "${SOURCE_DIR}")
  message(FATAL_ERROR "Usage: cmake -DCOMPONENT=mathlib|gaol -DSOURCE_DIR=<sources> -P ${CMAKE_CURRENT_LIST_FILE}")
endif()


################################################################################
# Helpers
################################################################################

# Reads into <var> the original of <file> (a path relative to SOURCE_DIR), which
# is saved as <file>.codac-orig the first time.
function(codac_gaol_read_original file var)
  set(_path "${SOURCE_DIR}/${file}")
  if(NOT EXISTS "${_path}.codac-orig")
    configure_file("${_path}" "${_path}.codac-orig" COPYONLY)
  endif()
  file(READ "${_path}.codac-orig" _content)
  set(${var} "${_content}" PARENT_SCOPE)
endfunction()

# Writes <content> to <file>, unless the file already holds it: an unchanged
# file keeps its date, and what is compiled from it is not compiled again.
function(codac_gaol_write file content)
  set(_path "${SOURCE_DIR}/${file}")
  if(EXISTS "${_path}")
    file(READ "${_path}" _current)
    if("${_current}" STREQUAL "${content}")
      return()
    endif()
  endif()
  file(WRITE "${_path}" "${content}")
endfunction()

# Replaces in the variable <var> the occurrences of <old> by <new>, which have
# to be <count> in number. CODAC_GAOL_FILE names the file in the error message.
function(codac_gaol_replace var count old new)
  string(LENGTH "${old}" _length_old)
  set(_rest "${${var}}")
  set(_found 0)
  string(FIND "${_rest}" "${old}" _position)
  while(NOT _position EQUAL -1)
    math(EXPR _found "${_found} + 1")
    math(EXPR _position "${_position} + ${_length_old}")
    string(SUBSTRING "${_rest}" ${_position} -1 _rest)
    string(FIND "${_rest}" "${old}" _position)
  endwhile()
  if(NOT _found EQUAL count)
    message(FATAL_ERROR "codac_gaol_patch.cmake: ${_found} occurrence(s) instead of ${count} in ${CODAC_GAOL_FILE} of:\n${old}")
  endif()
  string(REPLACE "${old}" "${new}" _content "${${var}}")
  set(${var} "${_content}" PARENT_SCOPE)
endfunction()

# Replaces <file>, a configuration file written by hand for one compiler, by a
# header including <header>, the configuration CMake generates from <template>
# for every compiler.
function(codac_gaol_redirect file template header)
  codac_gaol_read_original(${file} _original)
  codac_gaol_write(${file} "/* Written by Codac (scripts/CMakeModules/gaol/codac_gaol_patch.cmake) in place
   of the configuration file for this compiler: CMake generates the configuration
   for every compiler, from ${template}. */

#include \"${header}\"
")
endfunction()


################################################################################
# The CMake build
################################################################################

set(_files_dir "${CMAKE_CURRENT_LIST_DIR}/${COMPONENT}")
file(GLOB_RECURSE _files RELATIVE "${_files_dir}" "${_files_dir}/*")
foreach(_file ${_files})
  configure_file("${_files_dir}/${_file}" "${SOURCE_DIR}/${_file}" COPYONLY)
endforeach()


################################################################################
# mathlib
################################################################################

if(COMPONENT STREQUAL "mathlib")

  # src/mathlib_config.h includes these instead of mathlib_configuration.h for
  # Visual C++ and MinGW. The originals only select the implementation of
  # Init_Lib() for 32-bit x86, which mathlib_configuration.h.in does for every
  # target.
  foreach(_file src/mathlib_config_msvc.h src/mathlib_config_mingw.h)
    codac_gaol_redirect(${_file} mathlib_configuration.h.in mathlib_configuration.h)
  endforeach()

endif()


################################################################################
# GAOL
################################################################################

if(COMPONENT STREQUAL "gaol")

  # gaol/gaol_config.h and gaol/gaol_version.h include these instead of
  # gaol_configuration.h for Visual C++ and MinGW. The originals configure GAOL
  # for 32-bit x86 (assembly, SSE3, preserved rounding) and give it the version
  # 4.1.0.
  foreach(_file gaol/gaol_config_msvc.h gaol/gaol_config_mingw.h gaol/gaol_version_msvc.h)
    codac_gaol_redirect(${_file} gaol_configuration.h.in gaol/gaol_configuration.h)
  endforeach()


  # --------------------------------------------------------------------------
  # From the patch IBEX applies to GAOL
  # (https://github.com/ibex-team/ibex-lib,
  # interval_lib_wrapper/gaol/3rd/gaol-4.2.3alpha0.all.all.patch), by Gilles
  # Chabert. Its hunks are reproduced as they are, except the one on GAOL_NAN
  # in gaol/gaol_port.h, which replaces a static variable by a compound literal
  # of C99 that Visual C++ does not accept in C++.
  # --------------------------------------------------------------------------

  # At initialisation, GAOL writes to the control word of the floating-point
  # unit the value GAOL_FPU_MASK, which is meant for the x87 of x86 processors.
  # Elsewhere the same field of fenv_t is another register: on ARM64, the
  # rounding direction stays to nearest, and on an ARM64 Mac, where the value
  # enables floating-point traps, the program crashes. The default environment
  # is restored instead, and the rounding direction set to +oo.
  set(CODAC_GAOL_FILE gaol/gaol_common.cpp)
  codac_gaol_read_original(${CODAC_GAOL_FILE} _content)
  codac_gaol_replace(_content 1
[==[
            reset_fpu_cw(GAOL_FPU_MASK); // 53 bits precision, all exceptions masked, rounding to +oo
]==] [==[
	fesetenv(FE_DFL_ENV);
	round_upward();
	// next instruction crashes on MacOS ARM64 platform
	//reset_fpu_cw(GAOL_FPU_MASK); // 53 bits precision, all exceptions masked, rounding to +oo
]==])
  codac_gaol_write(${CODAC_GAOL_FILE} "${_content}")

  # interval::midpoint() and the output of an interval (operator<<) set the
  # rounding direction, and restore it with GAOL_RND_LEAVE(), which does nothing
  # when the rounding direction is not preserved, as here: after them, GAOL
  # would go on rounding to nearest, or downward. GAOL_RND_RESTORE() sets it back
  # to +oo in that case. codac2::Interval::mid() calls interval::midpoint().
  #
  # "operator string" is also qualified with its namespace.
  set(CODAC_GAOL_FILE gaol/gaol_interval.cpp)
  codac_gaol_read_original(${CODAC_GAOL_FILE} _content)
  codac_gaol_replace(_content 1
[==[
    //    double r = ((I.right()==0.0) ? 0.0 : I.right());  // Avoids printing -0
    GAOL_RND_ENTER();
]==] [==[
    //    double r = ((I.right()==0.0) ? 0.0 : I.right());  // Avoids printing -0
    GAOL_RND_PRESERVE();
	round_upward();
]==])
  codac_gaol_replace(_content 1
[==[
    GAOL_RND_LEAVE();
    return os;
]==] [==[
    GAOL_RND_RESTORE();
    return os;
]==])
  codac_gaol_replace(_content 1
[==[
  interval::operator string() const
]==] [==[
  interval::operator std::string() const
]==])
  codac_gaol_replace(_content 1
[==[
    unsigned short int _save_state=get_fpu_cw(); round_nearest();
]==] [==[
    GAOL_RND_PRESERVE();
    round_nearest();
]==])
  codac_gaol_replace(_content 1
[==[
    GAOL_RND_LEAVE();
    return middle;
]==] [==[
    GAOL_RND_RESTORE();
    return middle;
]==])
  codac_gaol_write(${CODAC_GAOL_FILE} "${_content}")

  # GAOL's macro opposite(x) is defined in an installed header, and would
  # replace any call to a function named opposite in the code including it. It
  # is renamed gaol_opposite(x), in gaol/gaol_port.h, which defines it, and in
  # the three versions of the double operations, which use it.
  set(CODAC_GAOL_FILE gaol/gaol_port.h)
  codac_gaol_read_original(${CODAC_GAOL_FILE} _gaol_port_h)
  codac_gaol_replace(_gaol_port_h 2 "#  define opposite(x)" "#  define gaol_opposite(x)")
  foreach(CODAC_GAOL_FILE gaol/gaol_double_op_apmathlib.h gaol/gaol_double_op_crlibm.h gaol/gaol_double_op_m.h)
    codac_gaol_read_original(${CODAC_GAOL_FILE} _content)
    codac_gaol_replace(_content 10 "opposite(" "gaol_opposite(")
    if(CODAC_GAOL_FILE STREQUAL "gaol/gaol_double_op_apmathlib.h")
      set(_gaol_double_op_apmathlib_h "${_content}")
    else()
      codac_gaol_write(${CODAC_GAOL_FILE} "${_content}")
    endif()
  endforeach()


  # --------------------------------------------------------------------------
  # For Visual C++, MinGW and the systems configure.ac does not know, from the
  # fork of GAOL by Fabrice Le Bars (https://github.com/lebarsfa/GAOL)
  # --------------------------------------------------------------------------

  # get_fpu_cw() and reset_fpu_cw() of gaol/gaol_fpu_fenv.h read and write the
  # control word of the floating-point unit through a field of fenv_t, which
  # GAOL names for Linux on x86 and ARM64 and for macOS only (CTRLWORD): GAOL
  # does not compile elsewhere, although it only uses these functions to save
  # and restore the rounding direction when it preserves it, which it does not
  # do here. On the other systems -- Visual C++, MinGW, Linux on 32-bit ARM,
  # BSD... -- they save and restore the rounding direction itself, with the
  # functions of <fenv.h>, which need no knowledge of fenv_t. The fork of GAOL
  # by Fabrice Le Bars names the field for each system instead ("Improve
  # compatibility with arm64 and MinGW 15", 2026), which the name MinGW-w64
  # gave it changing with its version 13 shows to be fragile.
  set(CODAC_GAOL_FILE gaol/gaol_fpu_fenv.h)
  codac_gaol_read_original(${CODAC_GAOL_FILE} _content)
  codac_gaol_replace(_content 1
[==[
INLINE unsigned short int get_fpu_cw()
{
  fenv_t tmp;
  fegetenv(&tmp);
  return CTRLWORD(tmp);
}

INLINE void reset_fpu_cw(unsigned short int st)
{
  fenv_t tmp;
  fegetenv(&tmp);
  CTRLWORD(tmp) = st;
  fesetenv(&tmp);
}
]==] [==[
#ifdef CTRLWORD
INLINE unsigned short int get_fpu_cw()
{
  fenv_t tmp;
  fegetenv(&tmp);
  return CTRLWORD(tmp);
}

INLINE void reset_fpu_cw(unsigned short int st)
{
  fenv_t tmp;
  fegetenv(&tmp);
  CTRLWORD(tmp) = st;
  fesetenv(&tmp);
}
#else
/* Codac: where the control word is not named, the rounding direction
   (see scripts/CMakeModules/gaol/codac_gaol_patch.cmake) */
INLINE unsigned short int get_fpu_cw()
{
  return (unsigned short int)fegetround();
}

INLINE void reset_fpu_cw(unsigned short int st)
{
  fesetround(st);
}
#endif
]==])
  codac_gaol_write(${CODAC_GAOL_FILE} "${_content}")

  # What MEMALIGN() allocates, GAOL frees with free() (gaol/gaol_allocator.h,
  # gaol/gaol_interval_sse.cpp). For MinGW, it allocates with _mm_malloc(),
  # whose memory only _mm_free() can free, and for Visual C++ and the systems
  # configure.ac does not know, GAOL stops with an error. MinGW gets malloc(),
  # as in the fork ("Corrected potential problems with MEMALIGN()", 2022), and
  # so does Visual C++: on Windows, GAOL is built without its SSE2 intervals
  # (see CMakeLists.txt), which are what needs memory aligned on 16 bytes. The
  # other systems get posix_memalign(), as Linux has.
  set(CODAC_GAOL_FILE gaol/gaol_port.h)
  codac_gaol_replace(_gaol_port_h 1
[==[
#if defined (__MINGW32__)
#  include <stdlib.h>
#  include <malloc.h>
#  define MEMALIGN(buf,boundary,size) (!(buf=_mm_malloc(size,boundary)))
]==] [==[
#if defined (__MINGW32__) || defined (_MSC_VER)
/* Codac: malloc(), whose memory free() frees, from the fork of GAOL by
   Fabrice Le Bars (see scripts/CMakeModules/gaol/codac_gaol_patch.cmake) */
#  include <stdlib.h>
#  define MEMALIGN(buf,boundary,size) (!(buf=malloc(size)))
]==])
  codac_gaol_replace(_gaol_port_h 1
[==[
#else
#  error "Don't know how to allocate memory aligned on a boundary"
#endif
]==] [==[
#else
/* Codac: any other POSIX system, as Linux */
#  include <stdlib.h>
#  define MEMALIGN(buf,boundary,size) posix_memalign(&buf,boundary,size)
#endif
]==])
  codac_gaol_write(${CODAC_GAOL_FILE} "${_gaol_port_h}")

  # For Visual C++, GAOL declares mathlib's functions __declspec(dllimport),
  # as for a DLL, which makes the linker warn wherever the static libultim is
  # linked. The definition is kept for when _MATHLIB_DLL_ is not already
  # defined, as in the fork ("Export functions problems", 2023), and
  # gaol_configuration.h.in defines it for Visual C++.
  set(CODAC_GAOL_FILE gaol/gaol_double_op_apmathlib.h)
  codac_gaol_replace(_gaol_double_op_apmathlib_h 1
[==[
#if defined (_MSC_VER)
# define _MATHLIB_DLL_ extern "C" __declspec(dllimport)
#elif defined (__GNUC__)
# define _MATHLIB_DLL_ extern "C" __attribute__ ((visibility("default")))
#else
# define _MATHLIB_DLL_ extern "C"
#endif
]==] [==[
/* Codac: unless already defined, from the fork of GAOL by Fabrice Le Bars
   (see scripts/CMakeModules/gaol/codac_gaol_patch.cmake) */
#ifndef _MATHLIB_DLL_
#if defined (_MSC_VER)
# define _MATHLIB_DLL_ extern "C" __declspec(dllimport)
#elif defined (__GNUC__)
# define _MATHLIB_DLL_ extern "C" __attribute__ ((visibility("default")))
#else
# define _MATHLIB_DLL_ extern "C"
#endif
#endif // _MATHLIB_DLL_
]==])
  codac_gaol_write(${CODAC_GAOL_FILE} "${_gaol_double_op_apmathlib_h}")

  # get_inexact() and clear_inexact(), which the lexer uses, are defined for
  # Visual C++ in 32-bit x86 assembly (gaol/sysdeps/gaol_exact_msvc.h), which
  # Visual C++ compiles for no other processor, and without it GAOL stops with an
  # error. As in the fork ("Improve compatibility with other build tools",
  # 2022), where <fenv.h> exists, Visual C++ gets their <fenv.h> version, which
  # gaol/gaol_exact.h includes so that the lexer defines them itself (gaol_exact.c,
  # which would include the assembly, is not compiled for Visual C++: see
  # CMakeLists.txt).
  set(CODAC_GAOL_FILE gaol/gaol_exact.h)
  codac_gaol_read_original(${CODAC_GAOL_FILE} _content)
  codac_gaol_replace(_content 1
[==[
#endif /* __gaol_exact_h__ */
]==] [==[
/* Codac: for Visual C++, the <fenv.h> version, from the fork of GAOL by
   Fabrice Le Bars (see scripts/CMakeModules/gaol/codac_gaol_patch.cmake) */
#if defined (_MSC_VER) && HAVE_FENV_H
#  include "gaol/sysdeps/gaol_exact_c99.h"
#endif

#endif /* __gaol_exact_h__ */
]==])
  codac_gaol_write(${CODAC_GAOL_FILE} "${_content}")

endif()
