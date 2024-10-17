set(CMAKE_VERBOSE_MAKEFILE ON)

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
include(CPM)

CPMFindPackage(
      NAME CMakeModules
      GIT_TAG v0.2.3
      GITHUB_REPOSITORY NuHepMC/CMakeModules
      DOWNLOAD_ONLY
  )
include(${CMakeModules_SOURCE_DIR}/NuHepMCModules.cmake)

include(ROOT)

if(NOT TARGET ROOT::ROOT)
  cmessage(FATAL_ERROR "MaCh3 Expected dependency target: ROOT::ROOT")
endif()

if(DEFINED ROOT_CXX_STANDARD AND ROOT_CXX_STANDARD GREATER CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD ${ROOT_CXX_STANDARD})
endif()

############################  DUNEAnaObj  ####################################

find_package(duneanaobj)

if(NOT duneanaobj_FOUND)
  CPMFindPackage(
    NAME duneanaobj
    GIT_TAG "patch-1"
    GITHUB_REPOSITORY dbarrow257/duneanaobj
    OPTIONS
    "STANDALONE_BUILD ON"
  )
endif()

if(NOT TARGET duneanaobj::All)
  cmessage(FATAL_ERROR "MaCh3 DUNE Expected dependency target: duneanaobj::All")
endif()

##################################  MaCh3  ######################################
#If MaCh3 was sourced find it, otherwise use CPM
find_package(MaCh3)

if(NOT MaCh3_FOUND)
  CPMFindPackage(
      NAME MaCh3
      GIT_TAG "develop"
      GITHUB_REPOSITORY mach3-software/MaCh3
  )
else()
  ##KS: This ensure that all executables that are in core will be moved
  FILE(GLOB MaCh3Exe $ENV{MaCh3_ROOT}/Diagnostics/*)
  FILE(COPY ${MaCh3Exe} DESTINATION ${CMAKE_BINARY_DIR}/Diagnostics/)
endif()

if(NOT TARGET MaCh3::All)
  cmessage(FATAL_ERROR "MaCh3 DUNE Expected dependency target: MaCh3::All, did you source setup.MaCh3.sh")
endif()

############################  C++ Compiler  ####################################
if (NOT DEFINED CMAKE_CXX_STANDARD OR "${CMAKE_CXX_STANDARD} " STREQUAL " ")
  SET(CMAKE_CXX_STANDARD 11)
endif()

if(DEFINED ROOT_CXX_STANDARD AND ROOT_CXX_STANDARD GREATER CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD ${ROOT_CXX_STANDARD})
endif()
cmessage(STATUS "CMAKE CXX Standard: ${CMAKE_CXX_STANDARD}")


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
if(PSYCHE_FOUND)
  target_compile_definitions(MaCh3DUNECompilerOptions INTERFACE PSYCHESETUP)
endif()
target_include_directories(MaCh3DUNECompilerOptions
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
    $<INSTALL_INTERFACE:include>)

install(TARGETS MaCh3DUNECompilerOptions
    EXPORT mach3dune-targets
    LIBRARY DESTINATION lib/)

add_subdirectory(splines)
add_subdirectory(samplePDFDUNE)
add_subdirectory(src)

add_library(MaCh3DUNEAll INTERFACE)
set_target_properties(MaCh3DUNEAll PROPERTIES EXPORT_NAME All)
target_link_libraries(MaCh3DUNEAll INTERFACE splinesDUNE SamplePDFDUNE MaCh3DUNECompilerOptions)
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
install(DIRECTORY inputs DESTINATION ${CMAKE_BINARY_DIR})
install(DIRECTORY configs DESTINATION ${CMAKE_BINARY_DIR})
install(DIRECTORY utils DESTINATION ${CMAKE_BINARY_DIR})

set(export_destinations
    ${CMAKE_INSTALL_PREFIX}/lib/cmake/
    ${CMAKE_INSTALL_PREFIX}/
)

foreach(dest ${export_destinations})
  install(EXPORT mach3dune-targets
          FILE MaCh3DUNETargets.cmake
          NAMESPACE MaCh3DUNE::
          DESTINATION ${dest}
  )
endforeach()

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
