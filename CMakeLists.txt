# CMake version check
cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

SET(MaCh3DUNE_VERSION 1.2.0)

# Create project (todo: set configurable CUDA)
project(MaCh3DUNE VERSION ${MaCh3DUNE_VERSION} LANGUAGES CXX)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Changes default install path to be a subdirectory of the build dir.
# Can set build dir at configure time with -DCMAKE_INSTALL_PREFIX=/install/path
if(CMAKE_INSTALL_PREFIX STREQUAL "" OR CMAKE_INSTALL_PREFIX STREQUAL
  "/usr/local")
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}")
elseif(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}")
endif()

# Use the compilers found in the path
find_program(CMAKE_C_COMPILER NAMES $ENV{CC} gcc PATHS ENV PATH NO_DEFAULT_PATH)
find_program(CMAKE_CXX_COMPILER NAMES $ENV{CXX} g++ PATHS ENV PATH NO_DEFAULT_PATH)

LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)

message("CMAKE SOURCE_DIR IS ${CMAKE_SOURCE_DIR}")

################################################################################
#                            Check Dependencies
################################################################################
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake/Modules)

#=================== CPM

# download CPM.cmake
file(
  DOWNLOAD
  https://github.com/cpm-cmake/CPM.cmake/releases/download/v0.43.1/CPM.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake
)
include(${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake)

#=================== ROOT
find_package(ROOT 6 REQUIRED)

#=================== DUNEAnaObj
if(NOT DEFINED DUNE_ANAOBJ_BRANCH)
  set(DUNE_ANAOBJ_BRANCH v03_14_00)
  set(DUNE_ANAOBJ_BRANCH ${DUNE_ANAOBJ_BRANCH} CACHE STRING "Specify which duneanaobj branch to use")
endif()

option(MaCh3_DUNE_USE_SRProxy "Whether to build proxy classes for Standard Record objects" OFF)

if(MaCh3_DUNE_USE_SRProxy)
  find_package(duneanaobj)

  if(NOT duneanaobj_FOUND)
    CPMFindPackage(
      NAME duneanaobj
      GIT_TAG feature/Clingification_rebase2024
      GITHUB_REPOSITORY luketpickering/duneanaobj
      VERSION 4.0.0
    )
  endif()

  if(NOT TARGET duneanaobj::all)
    message(FATAL_ERROR "MaCh3 DUNE Expected dependency target: duneanaobj::all")
  endif()
else()
  include(fetchduneanaobj)
  fetchduneanaobj(${DUNE_ANAOBJ_BRANCH})
endif()

##################################  MaCh3  ######################################
#If MaCh3 was sourced find it, otherwise use CPM
SET(MaCh3_FOUND FALSE)
# don't find mach3 core from a previous build in this directory or sad times occur
LIST(APPEND CMAKE_IGNORE_PATH ${CMAKE_CURRENT_BINARY_DIR})
find_package(MaCh3 2.7.1 EXACT QUIET)
set(MaCh3_CORE_BRANCH v2.7.1 CACHE STRING "Specify the MaCh3 core branch to use")
message(STATUS "Using MaCh3_CORE_BRANCH: ${MaCh3_CORE_BRANCH}")

if(NOT MaCh3_FOUND)
  message(STATUS "Didn't find MaCh3, attempting to use built in MaCh3")

  if(NOT DEFINED MaCh3_GPU_ENABLED)
    set(MaCh3_GPU_ENABLED ON)
  endif()

  if(NOT DEFINED MaCh3_DEBUG_ENABLED)
    set(MaCh3_DEBUG_ENABLED FALSE)
  endif()

  if(NOT DEFINED MaCh3_MULTITHREAD_ENABLED)
    set(MaCh3_MULTITHREAD_ENABLED TRUE)
  endif()

  # Options list construction
  set(MaCh3_OPTIONS
    "MaCh3_GPU_ENABLED ${MaCh3_GPU_ENABLED}"
    "MaCh3_DEBUG_ENABLED ${MaCh3_DEBUG_ENABLED}"
    "MaCh3_MULTITHREAD_ENABLED ${MaCh3_MULTITHREAD_ENABLED}"
    "MaCh3_CORE_BRANCH ${MaCh3_CORE_BRANCH}"
    "MaCh3_WERROR_ENABLED OFF" #required due to small uninitialized issue in MaCh3Core 2.6.0, revert when fixed
  )

  # Add LOG_LEVEL if defined
  if (LOG_LEVEL)
    list(APPEND MaCh3_OPTIONS "LOG_LEVEL ${LOG_LEVEL}")
  endif()

  list(APPEND MaCh3_OPTIONS "NuFastLinear_ENABLED TRUE")
  list(APPEND MaCh3_OPTIONS "CUDAProb3_ENABLED TRUE")
  list(APPEND MaCh3_OPTIONS "CUDAProb3Linear_ENABLED TRUE")

  CPMAddPackage(
    NAME MaCh3
    GIT_TAG ${MaCh3_CORE_BRANCH}
    GITHUB_REPOSITORY mach3-software/MaCh3
    OPTIONS
    ${MaCh3_OPTIONS}
  )
else()
  ##KS: This ensure that all executables that are in core will be moved
  FILE(GLOB MaCh3Exe $ENV{MaCh3_ROOT}/Diagnostics/*)
  FILE(COPY ${MaCh3Exe} DESTINATION ${CMAKE_BINARY_DIR}/Diagnostics/)

  FILE(GLOB MaCh3Exe $ENV{MaCh3_ROOT}/plotting/*)
  FILE(COPY ${MaCh3Exe} DESTINATION ${CMAKE_BINARY_DIR}/plotting/)
endif()

if(NOT TARGET MaCh3::All)
  message(FATAL_ERROR "MaCh3 DUNE Expected dependency target: MaCh3::All, did you source setup.MaCh3.sh")
endif()

############################  C++ Compiler  ####################################
get_property(MACH3_CXX_STANDARD GLOBAL PROPERTY MACH3_CXX_STANDARD)
if(NOT DEFINED MACH3_CXX_STANDARD)
  SET(CMAKE_CXX_STANDARD 14)
else()
  SET(CMAKE_CXX_STANDARD ${MACH3_CXX_STANDARD})
endif()
message(STATUS "CMAKE CXX Standard: ${CMAKE_CXX_STANDARD}")

############################  FLAGS  ####################################
include(CheckCXXCompilerFlag)

# KS: Treat warnings as errors. Only sensible with GCC: clang/AppleClang
# (and the macOS toolchain in general) diagnose things GCC doesn't, and warnings
# leaking out of our dependencies' headers would then break the build.
if(NOT DEFINED MaCh3DUNE_WERROR_ENABLED)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(MaCh3DUNE_WERROR_ENABLED ON)
  else()
    set(MaCh3DUNE_WERROR_ENABLED OFF)
  endif()
endif()
set(MaCh3DUNE_WERROR_ENABLED ${MaCh3DUNE_WERROR_ENABLED} CACHE BOOL "Treat compiler warnings as errors")
message(STATUS "MaCh3DUNE_WERROR_ENABLED: ${MaCh3DUNE_WERROR_ENABLED} (compiler: ${CMAKE_CXX_COMPILER_ID})")

add_library(DUNEMaCh3Warnings INTERFACE)

# Candidate warning flags. Not every compiler knows all of them (e.g. clang has
# no -Wuseless-cast), and with -Werror an unknown -W... option is itself a hard
# error, so each one is probed before being used.
set(MaCh3DUNE_WARNING_FLAGS
  -Wextra                 # Enable extra warning flags
  -Wall                   # Enable all standard warning flags
  -Wshadow                # Warn when a variable declaration shadows one from an outer scope
  -Wuninitialized         # Warn about uninitialized variables
  -Wnon-virtual-dtor      # Warn when a class with virtual functions has a non-virtual destructor
  -Woverloaded-virtual    # Warn when a function declaration hides a virtual function from a base class
  -Wformat=2              # Warn on security issues around functions that format output (ie printf)
  -Wunused                # Warn on anything being unused
  -Wredundant-decls       # Warn about multiple declarations of the same entity. Useful for code cleanup.
  -Wstrict-aliasing=2     # Helps detect potential aliasing issues that could lead to undefined behavior.
  -Wuseless-cast          # Warn if you perform a cast to the same type (only in GCC >= 4.8)
  -Wnull-dereference      # Warn if a null dereference is detected (only in GCC >= 6.0)
  -Wold-style-cast        # Warn for c-style casts
  -Wconversion            # Warn on type conversions that may lose data
  -Wformat-security       # Warn on functions that are potentially insecure for formatting
  -Walloca                # Warn if `alloca` is used, as it can lead to stack overflows
  -Wswitch-enum           # Warn if a `switch` statement on an enum does not cover all values
  #-Wfloat-equal          # Warn if floating-point values are compared directly
  #-Wpadded               # Warn when padding is added to a structure or class for alignment
)

foreach(flag ${MaCh3DUNE_WARNING_FLAGS})
  # make a variable name that is safe to use as a cache entry
  string(REGEX REPLACE "[^A-Za-z0-9]" "_" flag_var "HAVE_CXX_FLAG${flag}")
  check_cxx_compiler_flag("${flag}" ${flag_var})
  if(${flag_var})
    target_compile_options(DUNEMaCh3Warnings INTERFACE ${flag})
  else()
    message(STATUS "Compiler does not support ${flag}, skipping it")
  endif()
endforeach()

if(MaCh3DUNE_WERROR_ENABLED)
  target_compile_options(DUNEMaCh3Warnings INTERFACE -Werror)
endif()

# KS: Warnings coming out of our dependencies' headers are not ours to fix, so
# consume those headers as system includes. Imported targets already behave
# that way; targets we build ourselves (MaCh3 core and duneanaobj, both pulled
# in at configure time) have to be marked by hand.
function(MaCh3DUNEMarkTargetSystem tgt)
  if(NOT TARGET ${tgt})
    return()
  endif()

  # properties can't be set through an ALIAS, so resolve it first
  get_target_property(tgt_real ${tgt} ALIASED_TARGET)
  if(NOT tgt_real)
    set(tgt_real ${tgt})
  endif()

  get_target_property(tgt_imported ${tgt_real} IMPORTED)
  if(tgt_imported) # imported targets are treated as system already
    return()
  endif()

  if(NOT CMAKE_VERSION VERSION_LESS 3.25)
    set_target_properties(${tgt_real} PROPERTIES SYSTEM TRUE)
  else()
    # older CMake has no SYSTEM target property, so declare the interface
    # include dirs as system ones explicitly
    get_target_property(tgt_includes ${tgt_real} INTERFACE_INCLUDE_DIRECTORIES)
    if(tgt_includes)
      set_property(TARGET ${tgt_real} APPEND PROPERTY
        INTERFACE_SYSTEM_INCLUDE_DIRECTORIES ${tgt_includes})
    endif()
  endif()
endfunction()

function(MaCh3DUNEMarkDirectorySystem dir)
  get_property(dir_targets DIRECTORY ${dir} PROPERTY BUILDSYSTEM_TARGETS)
  foreach(tgt ${dir_targets})
    MaCh3DUNEMarkTargetSystem(${tgt})
  endforeach()

  get_property(subdirs DIRECTORY ${dir} PROPERTY SUBDIRECTORIES)
  foreach(subdir ${subdirs})
    MaCh3DUNEMarkDirectorySystem(${subdir})
  endforeach()
endfunction()

if(MaCh3_SOURCE_DIR AND IS_DIRECTORY "${MaCh3_SOURCE_DIR}")
  message(STATUS "Treating headers under ${MaCh3_SOURCE_DIR} as system headers")
  MaCh3DUNEMarkDirectorySystem(${MaCh3_SOURCE_DIR})
endif()

# targets defined outside their own subdirectory (e.g. by fetchduneanaobj)
foreach(dep duneanaobj::StandardRecord duneanaobj::all)
  MaCh3DUNEMarkTargetSystem(${dep})
endforeach()

################################# Features ##################################

LIST(APPEND ALL_FEATURES
)
message(STATUS "MaCh3DUNE Features: ")
foreach(f ${ALL_FEATURES})
  message(STATUS "     ${f}: ${MaCh3DUNE_${f}_ENABLED}")
endforeach()

#KS: Options to print dependency graph
if(NOT DEFINED MaCh3DUNE_DependancyGraph)
  set(MaCh3DUNE_DependancyGraph FALSE)
endif()

if(MaCh3DUNE_DependancyGraph)
  add_custom_target(graphviz ALL
    COMMAND ${CMAKE_COMMAND} "--graphviz=foo.dot" .
    COMMAND dot -Tpng foo.dot -o foo.png
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
  )
endif()

################################# Build MaCh3 ##################################
add_library(MaCh3DUNECompilerOptions INTERFACE)

if(MaCh3_DUNE_USE_SRProxy)
  target_compile_definitions(MaCh3DUNECompilerOptions INTERFACE MaCh3_DUNE_USE_SRProxy=1)
else()
  target_compile_definitions(MaCh3DUNECompilerOptions INTERFACE MaCh3_DUNE_USE_SRProxy=0)
endif()

set_target_properties(MaCh3DUNECompilerOptions PROPERTIES EXPORT_NAME CompilerOptions)

target_include_directories(MaCh3DUNECompilerOptions
  INTERFACE
  $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/Apps>
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>)

install(TARGETS MaCh3DUNECompilerOptions
  EXPORT mach3dune-targets
  LIBRARY DESTINATION lib/)

add_subdirectory(Splines)
add_subdirectory(Samples)
add_subdirectory(Apps)

add_library(MaCh3DUNEAll INTERFACE)
set_target_properties(MaCh3DUNEAll PROPERTIES EXPORT_NAME All)
target_link_libraries(MaCh3DUNEAll INTERFACE SamplesDUNE SplinesDUNE MaCh3DUNECompilerOptions)
add_library(MaCh3DUNE::All ALIAS MaCh3DUNEAll)

if(NOT DEFINED MaCh3_PREFIX OR MaCh3_PREFIX STREQUAL "")
    # CPMAddPackage sets MaCh3_SOURCE_DIR wherever the source actually came
    # from, so this also covers a local checkout supplied via
    # CPM_MaCh3_SOURCE, which does not land in _deps/mach3-src.
    if(MaCh3_SOURCE_DIR AND IS_DIRECTORY "${MaCh3_SOURCE_DIR}")
        set(MaCh3_PREFIX "${MaCh3_SOURCE_DIR}")
        message(STATUS "Detected MaCh3 source dir from CPM: ${MaCh3_PREFIX}")
    # Try to detect if MaCh3 was fetched by CPM
    elseif(EXISTS "${CMAKE_BINARY_DIR}/_deps/mach3-src")
        set(MaCh3_PREFIX "${CMAKE_BINARY_DIR}/_deps/mach3-src")
        message(STATUS "Detected CPM-fetched MaCh3: ${MaCh3_PREFIX}")
    # Try to get from MaCh3 package config
    elseif(DEFINED MaCh3_DIR)
        get_filename_component(MaCh3_PREFIX "${MaCh3_DIR}" DIRECTORY)
        get_filename_component(MaCh3_PREFIX "${MaCh3_PREFIX}" DIRECTORY)
        message(STATUS "Detected installed MaCh3 from MaCh3_DIR: ${MaCh3_PREFIX}")
    # Try to get from target properties
    elseif(TARGET MaCh3::All)
        get_target_property(MACH3_INCLUDE_DIRS MaCh3::All INTERFACE_INCLUDE_DIRECTORIES)
        if(MACH3_INCLUDE_DIRS)
            list(GET MACH3_INCLUDE_DIRS 0 FIRST_INCLUDE_DIR)
            get_filename_component(MaCh3_PREFIX "${FIRST_INCLUDE_DIR}" DIRECTORY)
            message(STATUS "Detected MaCh3 from target include dirs: ${MaCh3_PREFIX}")
        endif()
    endif()

    # Final validation
    if(NOT DEFINED MaCh3_PREFIX OR MaCh3_PREFIX STREQUAL "")
        message(FATAL_ERROR "Could not determine MaCh3_PREFIX. Please set it manually with -DMaCh3_PREFIX=<path>")
    endif()

    if(NOT EXISTS "${MaCh3_PREFIX}")
        message(FATAL_ERROR "MaCh3_PREFIX does not exist: ${MaCh3_PREFIX}")
    endif()

    # Cache it for use in subdirectories
    set(MaCh3_PREFIX "${MaCh3_PREFIX}" CACHE PATH "MaCh3 installation or source directory" FORCE)
    message(STATUS "Set MaCh3_PREFIX: ${MaCh3_PREFIX}")
endif()


if(PYTHON_ENABLED)
	message(STATUS "PYTHON IS ENABLED!!")
	# $ORIGIN is the ELF spelling and means nothing to dyld, which would look
	# for a directory literally named '$ORIGIN'. macOS spells it @loader_path.
	if(APPLE)
		set(MaCh3DUNE_RPATH_ORIGIN "@loader_path")
	else()
		set(MaCh3DUNE_RPATH_ORIGIN "$ORIGIN")
	endif()
	set(CMAKE_INSTALL_RPATH "${MaCh3DUNE_RPATH_ORIGIN}/../lib")
	message(STATUS "CMAKE_INSTALL_RPATH is ${CMAKE_INSTALL_RPATH}")
	set(PYMACH3_PATH ${CMAKE_BINARY_DIR}/python/pyMaCh3)
	add_subdirectory(python)
endif()

install(TARGETS MaCh3DUNEAll
  EXPORT mach3dune-targets
  LIBRARY DESTINATION lib/)

configure_file(cmake/Templates/setup.MaCh3DUNE.sh.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/setup.MaCh3DUNE.sh" @ONLY)
install(FILES
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/setup.MaCh3DUNE.sh" DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/output)

############################  Install  ####################################
install(DIRECTORY Inputs DESTINATION ${CMAKE_BINARY_DIR})
install(DIRECTORY Configs DESTINATION ${CMAKE_BINARY_DIR})
install(DIRECTORY plotting DESTINATION ${CMAKE_BINARY_DIR})
install(DIRECTORY Utils DESTINATION ${CMAKE_BINARY_DIR})

install(EXPORT mach3dune-targets
  FILE MaCh3DUNETargets.cmake
  NAMESPACE MaCh3DUNE::
  DESTINATION ${CMAKE_INSTALL_PREFIX}/
)

include(CMakePackageConfigHelpers)
configure_package_config_file(
  ${CMAKE_CURRENT_LIST_DIR}/cmake/Templates/MaCh3DUNEConfig.cmake.in ${CMAKE_BINARY_DIR}/MaCh3DUNEConfig.cmake
  INSTALL_DESTINATION
  /this/is/ignored/for/some/reason/thanks/kitware
  NO_SET_AND_CHECK_MACRO
  NO_CHECK_REQUIRED_COMPONENTS_MACRO)


write_basic_package_version_file(${CMAKE_BINARY_DIR}/MaCh3DUNEConfigVersion.cmake
  VERSION ${MaCh3DUNE_VERSION}
  COMPATIBILITY AnyNewerVersion)

install(FILES
  ${CMAKE_BINARY_DIR}/MaCh3DUNEConfig.cmake
  ${CMAKE_BINARY_DIR}/MaCh3DUNEConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/cmake/MaCh3DUNE)

include(mach3dune-config)
