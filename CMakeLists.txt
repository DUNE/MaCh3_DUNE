# CMake version check
cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

SET(MaCh3DUNE_VERSION 1.0.0)

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
  https://github.com/cpm-cmake/CPM.cmake/releases/download/v0.40.2/CPM.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake
)
include(${CMAKE_CURRENT_BINARY_DIR}/cmake/CPM.cmake)

CPMFindPackage(
  NAME CMakeModules
  GIT_TAG stable
  GITHUB_REPOSITORY NuHepMC/CMakeModules
  DOWNLOAD_ONLY
)
include(${CMakeModules_SOURCE_DIR}/NuHepMCModules.cmake)

#=================== ROOT
include(ROOT)

if(NOT TARGET ROOT::ROOT)
  cmessage(FATAL_ERROR "MaCh3 Expected dependency target: ROOT::ROOT")
endif()

if(DEFINED ROOT_CXX_STANDARD AND ROOT_CXX_STANDARD GREATER CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD ${ROOT_CXX_STANDARD})
endif()

#=================== DUNEAnaObj
option(Build_NDGAr "Whether to build ndgar SampleHandler class or all others" OFF)

if(NOT DEFINED DUNE_ANAOBJ_BRANCH)
  if (Build_NDGAr)
    set(DUNE_ANAOBJ_BRANCH tag/legacy_ndgar_cafs_2024)
  else()
    set(DUNE_ANAOBJ_BRANCH v03_06_00)
  endif()

  set(DUNE_ANAOBJ_BRANCH ${DUNE_ANAOBJ_BRANCH} CACHE STRING "Specify which duneanaobj branch to use")
endif()

if (Build_NDGAr)
  add_compile_definitions(BUILD_NDGAR)
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
    cmessage(FATAL_ERROR "MaCh3 DUNE Expected dependency target: duneanaobj::all")
  endif()
else()
  CPMAddPackage(
    NAME duneanaobj
    GIT_TAG ${DUNE_ANAOBJ_BRANCH}
    GITHUB_REPOSITORY DUNE/duneanaobj 
    DOWNLOAD_ONLY YES
  )
  include_directories(${duneanaobj_SOURCE_DIR})
  
  ROOT_GENERATE_DICTIONARY(StandardRecordDict ${duneanaobj_SOURCE_DIR}/duneanaobj/StandardRecord/StandardRecord.h
    LINKDEF ${duneanaobj_SOURCE_DIR}/duneanaobj/StandardRecord/classes_def.xml)
  
  file(GLOB SR_IMPL_FILES "${duneanaobj_SOURCE_DIR}/duneanaobj/StandardRecord/*.cxx")
  LIST(APPEND SR_IMPL_FILES ${CMAKE_CURRENT_BINARY_DIR}/StandardRecordDict.cxx)
  
  add_library(duneanaobj_StandardRecord SHARED ${SR_IMPL_FILES})
  target_link_libraries(duneanaobj_StandardRecord PUBLIC ROOT::MathCore)
  
  target_include_directories(duneanaobj_StandardRecord PUBLIC 
    $<BUILD_INTERFACE:${duneanaobj_SOURCE_DIR}> 
    $<INSTALL_INTERFACE:include>)
  target_include_directories(duneanaobj_StandardRecord PRIVATE 
    $<BUILD_INTERFACE:${duneanaobj_SOURCE_DIR}/duneanaobj/StandardRecord> #root puts this in the dictionary
  )
  
  set_target_properties(duneanaobj_StandardRecord PROPERTIES EXPORT_NAME all)
  install(TARGETS duneanaobj_StandardRecord EXPORT mach3dune-targets DESTINATION lib)
  
  file(GLOB SR_HEADER_FILES "${duneanaobj_SOURCE_DIR}/duneanaobj/StandardRecord/*.h")
  
  install(FILES ${SR_HEADER_FILES} DESTINATION include/duneanaobj/StandardRecord)
  install(FILES 
    ${CMAKE_CURRENT_BINARY_DIR}/libStandardRecordDict_rdict.pcm 
    ${CMAKE_CURRENT_BINARY_DIR}/libStandardRecordDict.rootmap 
    DESTINATION lib)
  
  add_library(duneanaobj::all ALIAS duneanaobj_StandardRecord)
  
endif()

##################################  MaCh3  ######################################
#If MaCh3 was sourced find it, otherwise use CPM
SET(MaCh3_FOUND FALSE)
find_package(MaCh3 2.2.0 EXACT QUIET)
set(MaCh3_CORE_BRANCH v2.2.0 CACHE STRING "Specify the MaCh3 core branch to use")
message(STATUS "Using MaCh3_CORE_BRANCH: ${MaCh3_CORE_BRANCH}")
option(MaCh3DUNE_Atmos "Whether working with atmospherics or beam oscillations" OFF)

if(NOT MaCh3_FOUND)
  cmessage(STATUS "Didn't find MaCh3, attempting to use built in MaCh3")
  
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
  )

  # Add LOG_LEVEL if defined
  if (LOG_LEVEL)
    list(APPEND MaCh3_OPTIONS "LOG_LEVEL ${LOG_LEVEL}")
  endif()

  list(APPEND MaCh3_OPTIONS "NuFastLinear_ENABLED TRUE")
  if (MaCh3DUNE_Atmos)
    list(APPEND MaCh3_OPTIONS "CUDAProb3_ENABLED TRUE")
    list(APPEND MaCh3_OPTIONS "CUDAProb3Linear_ENABLED FALSE")
  else()
    list(APPEND MaCh3_OPTIONS "CUDAProb3Linear_ENABLED TRUE")
  endif()
  
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
  cmessage(FATAL_ERROR "MaCh3 DUNE Expected dependency target: MaCh3::All, did you source setup.MaCh3.sh")
endif()

############################  C++ Compiler  ####################################
get_property(MACH3_CXX_STANDARD GLOBAL PROPERTY MACH3_CXX_STANDARD)
if(NOT DEFINED MACH3_CXX_STANDARD)
  SET(CMAKE_CXX_STANDARD 14)
else()
  SET(CMAKE_CXX_STANDARD ${MACH3_CXX_STANDARD})
endif()
cmessage(STATUS "CMAKE CXX Standard: ${CMAKE_CXX_STANDARD}")

############################  FLAGS  ####################################
add_library(DUNEMaCh3Warnings INTERFACE)

target_compile_options(DUNEMaCh3Warnings INTERFACE
  -Werror                 # Treat Warnings as Errors
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

################################# Features ##################################

LIST(APPEND ALL_FEATURES
)
cmessage(STATUS "MaCh3DUNE Features: ")
foreach(f ${ALL_FEATURES})
  cmessage(STATUS "     ${f}: ${MaCh3DUNE_${f}_ENABLED}")
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
