# ==================================================================
#  Codac - cmake module for GAOL, the interval arithmetic library
# ==================================================================
#
# codac2::Interval derives from gaol::interval (see
# src/core/domains/interval/codac2_Interval.h). GAOL is written by Frederic
# Goualard (https://frederic.goualard.net, https://github.com/goualard-f/GAOL)
# and distributed under the GNU LGPL v2. It computes its elementary functions
# with mathlib, the IBM Accurate Portable Mathematical Library (libultim), which
# Frederic Goualard distributes along with it, under the GNU GPL v2 or later.
#
# Codac builds the GAOL of Jordan Ninin's fork, https://github.com/Jordan08/GAOL,
# which adds to GAOL a CMake build, the fixes Codac depends on (those of the
# patch IBEX applies to GAOL, and those Visual C++, MinGW and ARM processors
# need) and tests of its bounds. Fixes of GAOL go into the fork rather than into
# Codac. The fork's CMake build downloads mathlib from Frederic Goualard's site,
# and refuses the compilers that do not compute GAOL's intervals right (see its
# README.md).
#
# Codac used to reach both through IBEX, which bundles them. This module does
# for Codac what IBEX does for itself -- it sets the compilation flags interval
# arithmetic needs, and builds GAOL when none is installed -- and the top-level
# CMakeLists.txt uses it the way it uses FetchContent for Eigen and Catch2:
# find_package(GAOL) first (scripts/CMakeModules/FindGAOL.cmake), a build of
# its own only when that fails.
#
# Origin
# ------
# What follows comes from IBEX (https://github.com/ibex-team/ibex-lib, GNU
# LGPL v3):
#
#  - codac_gaol_compiler_flags() is the part of ibex_init_common()
#    (cmake.utils/ibex-config-utils.cmake) that sets the floating-point flags,
#    in the version of the IBEX fork maintained by Fabrice Le Bars,
#    https://github.com/lebarsfa/ibex-lib, tag ibex-2.8.9.20260819, which Codac
#    was built against until now. ibex_init_common() was written by Cyril
#    Bouvier for the CMake build of IBEX; the flags it sets -- -frounding-math
#    and the others that IEEE 754 double support depends on, and their Visual
#    Studio counterparts -- were added to it by Fabrice Le Bars.
#
#  - codac_gaol_build() is the "not found, install it" branch of
#    interval_lib_wrapper/gaol/CMakeLists.txt, written by Cyril Bouvier: build
#    mathlib and GAOL, and install the result next to the library it serves, as
#    IBEX does in include/ibex/3rd and lib/ibex/3rd. The CMake build of GAOL and
#    mathlib it runs is the fork's, which says where it comes from: from IBEX
#    (Cyril Bouvier, Gilles Chabert), with the portability fixes for Visual
#    C++, MinGW and ARM of the forks of GAOL and mathlib by Fabrice Le Bars.
#    The library paths are computed as lib_get_abspath_from_name() of
#    cmake.utils/IbexUtils.cmake computes them, and the include directory is
#    created before the build so that the generation step accepts it.
#
#  - codac_gaol_create_targets() and codac_gaol_config_snippet() do, for Codac,
#    what create_target_import_and_export() of cmake.utils/IbexUtils.cmake
#    (Cyril Bouvier) does for IBEX: an imported target for the build tree, and
#    the lines that recreate it for a consumer of the installed library.
#
# What differs from IBEX
# ----------------------
# IBEX extracts GAOL and mathlib from archives kept in its own repository, and
# builds them as a part of itself (add_subdirectory()). Here GAOL is downloaded
# from the fork, at the head of its master branch, and built as a project of its
# own (ExternalProject), for the reasons given at codac_gaol_build(); the fork's
# build downloads mathlib, checked against the SHA256 of its archive. The other
# differences are explained where they occur.


# Where a GAOL built by codac_gaol_build() is installed, under the installation
# prefix. IBEX keeps the libraries it builds for itself in include/ibex/3rd and
# lib/ibex/3rd, out of the way of a GAOL installed separately; include/codac
# cannot play that part here, being the name of Codac's umbrella header.
set(CODAC_INSTALL_INCLUDEDIR_3RD "${CMAKE_INSTALL_INCLUDEDIR}/codac-3rd")
set(CODAC_INSTALL_LIBDIR_3RD "${CMAKE_INSTALL_LIBDIR}/codac-3rd")


################################################################################
# codac_gaol_compiler_flags(<outvar>)
################################################################################
#
# Returns in <outvar> the compilation flags interval arithmetic depends on.
# Without -frounding-math and its companions, the compiler is free to evaluate
# a floating-point expression at compile time in the default rounding mode, or
# to contract it, and a bound computed that way no longer encloses anything.
#
# From ibex_init_common() (see the top of this file), which appends the same
# flags to CMAKE_C_FLAGS and CMAKE_CXX_FLAGS. They are returned as a list here,
# one flag per element, because Codac also hands them over to its consumers,
# through CODAC_CXX_FLAGS and codac.pc; for the same reason, the Visual Studio
# "/D NAME" pairs of IBEX are spelt "/DNAME", which target_compile_options()
# cannot mistake for two duplicated "/D" and merge. The rest of
# ibex_init_common() is left out: the installation directories, the build type,
# the C++ standard and the uninstall target are Codac's own decisions, and the
# Debug-only -Wall -DDEBUG (/D DEBUG) have no use here, Codac choosing its own
# warnings and none of its code reading DEBUG.
#
# The two flags IBEX added for filib, which is not GAOL, are kept all the same,
# so that Codac goes on being compiled with exactly the flags it had.
function(codac_gaol_compiler_flags outvar)

  include(CheckCXXCompilerFlag)
  set(flags "")

  # Flags related to portability and IEEE 754 double support
  if(MSVC)
    list(APPEND flags /D_CRT_SECURE_NO_WARNINGS /D_CRT_NONSTDC_NO_WARNINGS /fp:strict /Zc:__cplusplus /Zc:strictStrings-)
  else()
    # Each flag is kept only where the compiler takes it. The check results are
    # cached under the names ibex_init_common() gives them, e.g.
    # COMPILER_SUPPORTS_FROUNDING_MATH for -frounding-math.
    foreach(flag -frounding-math -ffloat-store -fno-fast-math -ffp-contract=off
                 -ffp-mode=full -fp-model=strict -fp:strict -mpc64)
      string(MAKE_C_IDENTIFIER "${flag}" _flag_id)
      string(TOUPPER "COMPILER_SUPPORTS${_flag_id}" _flag_var)
      check_cxx_compiler_flag("${flag}" ${_flag_var})
      if(${_flag_var})
        list(APPEND flags ${flag})
      endif()
    endforeach()
    # Due to warnings on macOS with filib
    if(APPLE)
      check_cxx_compiler_flag("-Wno-undefined-var-template" COMPILER_SUPPORTS_WNO_UNDEFINED_VAR_TEMPLATE)
      if(COMPILER_SUPPORTS_WNO_UNDEFINED_VAR_TEMPLATE)
        list(APPEND flags -Wno-undefined-var-template)
      endif()
    endif()

    # Not from ibex_init_common(), but Codac's own: on a 32-bit x86 processor,
    # doubles are computed in SSE2 rather than on the x87 FPU. Visual Studio,
    # in the branch above, computes them in SSE2 already.
    #
    # Computed on the x87, GAOL's bounds and mathlib's results are only right
    # while the precision of the x87 is set to 53 bits, and nothing keeps it
    # so: mathlib's Init_Lib() sets it where mathlib has a version for 32-bit x86
    # (see cmake/mathlib/mathlib_configuration.h.in in the fork of GAOL), but
    # GAOL, initialised right after, restores the default floating-point
    # environment (gaol::init()), whose precision is 64 bits on Linux and
    # with MinGW. Built for an i686 computing on the x87 (Clang 21 with
    # -mcpu=i686), GAOL returned [1.99975, 1.99975] for exp([1,1]), and the
    # bounds of exp, sin and cos missed the exact value for 4000, 3302 and 3913
    # of 4000 random arguments; built with the two flags below, the same
    # program gave all of its 48000 bounds bit for bit as on x86_64. In SSE2,
    # the precision is not a setting. SSE2 asks nothing more of the processor
    # than GAOL's own builds do: they compile GAOL for SSE2 on these systems.
    include(CheckCXXSourceCompiles)
    check_cxx_source_compiles("
      #if !defined(__i386__)
      #error not a 32-bit x86 target
      #endif
      int main() { return 0; }"
      CODAC_TARGET_IS_X86_32)
    if(CODAC_TARGET_IS_X86_32)
      check_cxx_compiler_flag("-msse2 -mfpmath=sse" COMPILER_SUPPORTS_MSSE2_MFPMATH_SSE)
      if(COMPILER_SUPPORTS_MSSE2_MFPMATH_SSE)
        list(APPEND flags -msse2 -mfpmath=sse)
      endif()
    endif()

    # Not from ibex_init_common() either: a warning when the compiler takes
    # -frounding-math but says it does not honour the rounding direction on the
    # target, as Clang does for 32-bit ARM processors ("overriding currently
    # unsupported rounding mode on this target"). It then optimises the
    # negations by which GAOL rounds downward with the rounding direction set
    # upward, and no flag or change to GAOL can prevent it. Built by Clang 21 for
    # 32-bit ARM, 4556 of 16000 random products, squares and cubes computed by
    # GAOL did not enclose their exact value; built by GCC 15, none. The build
    # of GAOL by codac_gaol_build() refuses such a compiler (the fork's
    # CMakeLists.txt); the warning remains for a GAOL found on the system.
    if(COMPILER_SUPPORTS_FROUNDING_MATH)
      set(CMAKE_REQUIRED_FLAGS "-frounding-math")
      check_cxx_source_compiles("int main() { return 0; }" CODAC_COMPILER_HONOURS_ROUNDING_MATH
                                FAIL_REGEX "unsupported rounding mode")
      unset(CMAKE_REQUIRED_FLAGS)
      if(NOT CODAC_COMPILER_HONOURS_ROUNDING_MATH)
        message(WARNING "${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} does not honour the rounding "
                        "direction on this target (-frounding-math): the intervals computed by Codac may not "
                        "enclose the values they should. Use a compiler that does, such as GCC.")
      endif()
    endif()
  endif()
  # Claim IEEE 754 double compatibility, for filib
  list(APPEND flags -D__STDC_IEC_559__=1)

  set(${outvar} ${flags} PARENT_SCOPE)
endfunction()


################################################################################
# codac_gaol_build()
################################################################################
#
# Downloads GAOL, builds it and mathlib in Release and installs them in the
# build tree; then sets GAOL_INCDIR, GAOL_LIB, MATHLIB_INCDIR, MATHLIB_LIB and
# GAOL_VERSION in the caller's scope, as find_package(GAOL) would have, and
# CODAC_GAOL_BUILT_HERE to TRUE. Both libraries are also installed with Codac,
# in CODAC_INSTALL_INCLUDEDIR_3RD and CODAC_INSTALL_LIBDIR_3RD, since the Codac
# libraries are of no use without them.
#
# GAOL comes from Jordan Ninin's fork (see the top of this file), whose CMake
# build downloads mathlib 2.1.1 from Frederic Goualard's site, the archive
# GAOL's README points to, and builds it along with GAOL.
#
# An ExternalProject, as in IBEX's build of GAOL with its autotools, rather than
# the FetchContent that Eigen and Catch2 are brought in with, which would build
# GAOL as a part of this project. Kept apart, GAOL is compiled with Codac's
# floating-point flags, which are passed on below, but not with its warnings,
# sanitizers and coverage instrumentation, and it is always built in Release,
# whatever the configuration of Codac, which is what the MSVC runtime choice of
# the top-level CMakeLists.txt counts on.
#
# The download happens at build time, when the target first needs it, and so
# does the update of the sources with the fork's master branch, at each build.
function(codac_gaol_build)

  include(ExternalProject)

  set(_gaol_prefix "${CMAKE_BINARY_DIR}/_deps/gaol")
  set(_gaol_install "${_gaol_prefix}/install")

  # Absolute paths of the libraries, as lib_get_abspath_from_name() of IBEX
  # computes them. Both are built as static libraries.
  set(_gaol_lib "${_gaol_install}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}gaol${CMAKE_STATIC_LIBRARY_SUFFIX}")
  set(_mathlib_lib "${_gaol_install}/lib/${CMAKE_STATIC_LIBRARY_PREFIX}ultim${CMAKE_STATIC_LIBRARY_SUFFIX}")

  # What IBEX passes to the builds of mathlib and GAOL is
  # CFLAGS=${CMAKE_C_FLAGS} CXXFLAGS=${CMAKE_CXX_FLAGS}, the interval arithmetic
  # flags of ibex_init_common() included. CMAKE_CXX_FLAGS already holds them
  # here (see the top-level CMakeLists.txt); they are added to the C flags for
  # mathlib, which is C.
  string(REPLACE ";" " " _interval_flags "${CODAC_INTERVAL_CXX_FLAGS}")
  set(_c_flags "${CMAKE_C_FLAGS} ${_interval_flags}")

  # ExternalProject hands the generator of this build to the sub-builds, with
  # its platform (-A) and toolset (-T), but nothing else: the compilers, the
  # flags and the target have to be passed on explicitly, or the sub-builds
  # would build GAOL for another machine than the one Codac is built for.
  #
  # Every value below has to come out the same at each configuration of Codac
  # that changes nothing: the sub-builds are configured again as soon as one of
  # them differs, and everything in Codac that includes an interval is then
  # recompiled. This is why the top-level CMakeLists.txt enables C in project()
  # rather than leaving it to a dependency, which changed CMAKE_C_COMPILER
  # between the first configuration and the second.
  set(_args
    -DCMAKE_BUILD_TYPE:STRING=Release
    # Codac's Python modules link these archives into shared libraries.
    -DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=ON
    -DCMAKE_INSTALL_PREFIX:PATH=${_gaol_install}
    # The paths above, whatever GNUInstallDirs would choose (lib64 on some
    # systems)
    -DCMAKE_INSTALL_LIBDIR:PATH=lib
    -DCMAKE_INSTALL_INCLUDEDIR:PATH=include
    # The fork's continuous integration runs its tests
    -DGAOL_BUILD_TESTS:BOOL=OFF
    # The mathlib the fork downloads and builds, never one installed on this
    # machine: the paths above, and the installation of Codac, count on it
    -DGAOL_FIND_MATHLIB:BOOL=OFF
    -DCMAKE_CXX_FLAGS:STRING=${CMAKE_CXX_FLAGS}
    -DCMAKE_C_FLAGS:STRING=${_c_flags}
  )
  # The Visual Studio and Xcode generators take their compilers from the
  # toolset and ignore these.
  if(NOT CMAKE_GENERATOR MATCHES "Visual Studio|Xcode")
    list(APPEND _args -DCMAKE_CXX_COMPILER:FILEPATH=${CMAKE_CXX_COMPILER})
    list(APPEND _args -DCMAKE_C_COMPILER:FILEPATH=${CMAKE_C_COMPILER})
  endif()
  foreach(_var CMAKE_TOOLCHAIN_FILE CMAKE_MSVC_RUNTIME_LIBRARY
               CMAKE_OSX_DEPLOYMENT_TARGET CMAKE_OSX_SYSROOT)
    if(${_var})
      list(APPEND _args "-D${_var}:STRING=${${_var}}")
    endif()
  endforeach()
  # A list, whose ";" would otherwise split the argument: LIST_SEPARATOR below
  # turns the "|" back into ";" for the sub-builds.
  if(CMAKE_OSX_ARCHITECTURES)
    string(REPLACE ";" "|" _osx_architectures "${CMAKE_OSX_ARCHITECTURES}")
    list(APPEND _args "-DCMAKE_OSX_ARCHITECTURES:STRING=${_osx_architectures}")
  endif()
  # Only when this build is itself a cross-compilation (as the macOS jobs of
  # .github/workflows are, by setting CMAKE_SYSTEM_NAME explicitly): passing
  # the host's own name would make the sub-builds believe they cross-compile.
  if(CMAKE_CROSSCOMPILING)
    list(APPEND _args -DCMAKE_SYSTEM_NAME:STRING=${CMAKE_SYSTEM_NAME})
    if(CMAKE_SYSTEM_PROCESSOR)
      list(APPEND _args -DCMAKE_SYSTEM_PROCESSOR:STRING=${CMAKE_SYSTEM_PROCESSOR})
    endif()
  endif()

  # Release whatever the configuration built: --config is what a
  # multi-configuration generator reads, and what the others ignore.
  set(_build_command ${CMAKE_COMMAND} --build <BINARY_DIR> --config Release)
  set(_install_command ${CMAKE_COMMAND} --build <BINARY_DIR> --config Release --target install)

  # GAOL, from Jordan Ninin's fork, at the head of its master branch (version
  # 4.2.3 of GAOL), so that the fixes pushed to the fork reach Codac without a
  # change here. Cloned with Git, the sources are brought up to date with the
  # branch at each build of Codac (the update step of ExternalProject, which
  # needs network access), and only what a new commit changes is compiled
  # again. Without Git, the archive GitHub makes of the branch is downloaded
  # instead, whose checksum changes with each commit and is not checked, and
  # only once per build directory: a new commit reaches such a build once
  # _deps/gaol is deleted from it. The fork's CMake build downloads mathlib
  # 2.1.1 from Frederic Goualard's site, checks its checksum, and builds and
  # installs it along with GAOL.
  find_package(Git QUIET)
  if(GIT_FOUND)
    set(_gaol_source
      GIT_REPOSITORY https://github.com/Jordan08/GAOL.git
      GIT_TAG master
      GIT_SHALLOW TRUE
      # The files as they are in the repository, whatever core.autocrlf says
      # (true on the Windows runners of GitHub Actions).
      GIT_CONFIG core.autocrlf=false
      LOG_UPDATE 1
    )
  else()
    set(_gaol_source
      URL https://github.com/Jordan08/GAOL/archive/refs/heads/master.zip
    )
    # Dates the extracted files from their extraction, so that a new download
    # rebuilds them. CMake 3.24 introduced the option, and warns when it is not
    # given.
    if(NOT CMAKE_VERSION VERSION_LESS 3.24)
      list(APPEND _gaol_source DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
    endif()
  endif()

  ExternalProject_Add(codac_gaol
    PREFIX "${_gaol_prefix}/gaol"
    ${_gaol_source}
    #URL ${CMAKE_SOURCE_DIR}/3rd/GAOL-master.zip # If needed to be self-contained...
    LIST_SEPARATOR |
    CMAKE_CACHE_ARGS ${_args}
    BUILD_COMMAND ${_build_command}
    INSTALL_COMMAND ${_install_command}
    # Named, so that Ninja knows which step produces the files Codac links.
    BUILD_BYPRODUCTS "${_gaol_lib}" "${_mathlib_lib}"
    LOG_DOWNLOAD 1
    LOG_CONFIGURE 1
    LOG_BUILD 1
    LOG_INSTALL 1
    LOG_OUTPUT_ON_FAILURE 1
  )

  # This directory needs to exist during generation step to avoid CMake errors
  file(MAKE_DIRECTORY "${_gaol_install}/include")

  # Install headers and libraries
  install(DIRECTORY "${_gaol_install}/include/" DESTINATION ${CODAC_INSTALL_INCLUDEDIR_3RD})
  install(FILES "${_gaol_lib}" "${_mathlib_lib}" DESTINATION ${CODAC_INSTALL_LIBDIR_3RD})

  set(GAOL_INCDIR "${_gaol_install}/include" PARENT_SCOPE)
  set(GAOL_LIB "${_gaol_lib}" PARENT_SCOPE)
  set(MATHLIB_INCDIR "${_gaol_install}/include" PARENT_SCOPE)
  set(MATHLIB_LIB "${_mathlib_lib}" PARENT_SCOPE)
  set(GAOL_VERSION "4.2.3" PARENT_SCOPE)
  set(CODAC_GAOL_BUILT_HERE TRUE PARENT_SCOPE)
endfunction()


################################################################################
# codac_gaol_create_targets()
################################################################################
#
# Makes the imported targets Codac::gaol and Codac::ultim out of the variables
# set by find_package(GAOL) or codac_gaol_build(), so that the rest of the
# build links one name whichever of the two provided GAOL. Codac::gaol links
# Codac::ultim, as Ibex::gaol links Ibex::ultim in IBEX.
#
# UNKNOWN rather than STATIC: a GAOL found on the system may just as well be a
# shared library, and CMake links either from the path alone.
function(codac_gaol_create_targets)

  add_library(Codac::ultim UNKNOWN IMPORTED)
  set_target_properties(Codac::ultim PROPERTIES
                        IMPORTED_LOCATION "${MATHLIB_LIB}"
                        INTERFACE_INCLUDE_DIRECTORIES "${MATHLIB_INCDIR}")

  add_library(Codac::gaol UNKNOWN IMPORTED)
  set_target_properties(Codac::gaol PROPERTIES
                        IMPORTED_LOCATION "${GAOL_LIB}"
                        INTERFACE_INCLUDE_DIRECTORIES "${GAOL_INCDIR}"
                        INTERFACE_LINK_LIBRARIES Codac::ultim)

  # For Visual C++, GAOL declares its classes and functions
  # __declspec(dllimport), as for a DLL, unless __GAOL_PUBLIC__ is defined
  # (gaol/gaol_config.h). None of GAOL's builds -- its autotools and meson
  # builds, and the CMake build of the fork -- defines _COMPILING__GAOL_PUBLIC__,
  # with which GAOL would export a DLL: the GAOL Visual C++ links is a static
  # library, whether codac_gaol_build() built it or it was found on the system,
  # which every file including its headers has to be told. The CMake build of
  # the fork tells it the same way to the targets linking gaol::gaol.
  if(MSVC)
    set_target_properties(Codac::gaol PROPERTIES INTERFACE_COMPILE_DEFINITIONS "__GAOL_PUBLIC__=")
  endif()

  # Nothing builds an imported target, but the dependencies given to one are
  # followed by every target linking it: whatever links Codac::gaol waits for
  # the build of GAOL to be over.
  if(CODAC_GAOL_BUILT_HERE)
    add_dependencies(Codac::ultim codac_gaol)
    add_dependencies(Codac::gaol codac_gaol)
  endif()
endfunction()


################################################################################
# codac_gaol_config_snippet(<outvar>)
################################################################################
#
# Returns the lines of codac-config.cmake that define Codac::gaol and
# Codac::ultim for a consumer of the installed Codac, which the exported Codac
# targets name -- the counterpart of the ibex-config-gaol.cmake and
# ibex-config-ultim.cmake files that create_target_import_and_export() writes
# for IBEX. A GAOL found on this machine is named by the paths it was found at;
# a GAOL built by Codac by paths relative to the configuration file, since it
# is installed under the same prefix and the prefix may be moved.
function(codac_gaol_config_snippet outvar)

  if(CODAC_GAOL_BUILT_HERE)
    file(RELATIVE_PATH _to_prefix "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_CMAKE}" "${CMAKE_INSTALL_PREFIX}")
    get_filename_component(_gaol_name "${GAOL_LIB}" NAME)
    get_filename_component(_ultim_name "${MATHLIB_LIB}" NAME)
    set(_prefix_line "get_filename_component(_codac_prefix \"\${CMAKE_CURRENT_LIST_DIR}/${_to_prefix}\" ABSOLUTE)")
    set(_gaol_loc "\${_codac_prefix}/${CODAC_INSTALL_LIBDIR_3RD}/${_gaol_name}")
    set(_ultim_loc "\${_codac_prefix}/${CODAC_INSTALL_LIBDIR_3RD}/${_ultim_name}")
    set(_gaol_inc "\${_codac_prefix}/${CODAC_INSTALL_INCLUDEDIR_3RD}")
    set(_ultim_inc "\${_codac_prefix}/${CODAC_INSTALL_INCLUDEDIR_3RD}")
  else()
    set(_prefix_line "")
    set(_gaol_loc "${GAOL_LIB}")
    set(_ultim_loc "${MATHLIB_LIB}")
    set(_gaol_inc "${GAOL_INCDIR}")
    set(_ultim_inc "${MATHLIB_INCDIR}")
  endif()

  # As in codac_gaol_create_targets()
  set(_gaol_definitions "")
  if(MSVC)
    set(_gaol_definitions "\n                            INTERFACE_COMPILE_DEFINITIONS \"__GAOL_PUBLIC__=\"")
  endif()

  set(${outvar} "
    # GAOL, the interval arithmetic library Codac is built upon, and mathlib
    # (libultim), which GAOL depends on. The exported Codac targets link
    # Codac::gaol, so both have to be defined before they are included.
    ${_prefix_line}
    if(NOT TARGET Codac::ultim)
      add_library(Codac::ultim UNKNOWN IMPORTED)
      set_target_properties(Codac::ultim PROPERTIES
                            IMPORTED_LOCATION \"${_ultim_loc}\"
                            INTERFACE_INCLUDE_DIRECTORIES \"${_ultim_inc}\")
    endif()
    if(NOT TARGET Codac::gaol)
      add_library(Codac::gaol UNKNOWN IMPORTED)
      set_target_properties(Codac::gaol PROPERTIES
                            IMPORTED_LOCATION \"${_gaol_loc}\"
                            INTERFACE_INCLUDE_DIRECTORIES \"${_gaol_inc}\"
                            INTERFACE_LINK_LIBRARIES Codac::ultim${_gaol_definitions})
    endif()
" PARENT_SCOPE)
endfunction()
