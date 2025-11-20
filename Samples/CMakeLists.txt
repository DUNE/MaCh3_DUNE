set(HEADERS
  SampleHandlerAtm.h
  SampleHandlerBeamFD.h
  SampleHandlerBeamND.h
  SampleHandlerBeamNDGAr.h
  MaCh3DUNEFactory.h
)
set(SOURCE_FILES
  SampleHandlerAtm.cpp
  SampleHandlerBeamFD.cpp
  SampleHandlerBeamND.cpp
  SampleHandlerBeamNDGAr.cpp
  MaCh3DUNEFactory.cpp
)

add_library(SamplesDUNE SHARED
  ${SOURCE_FILES}
)

set_target_properties(SamplesDUNE PROPERTIES
  PUBLIC_HEADER "${HEADERS}"
  EXPORT_NAME SamplesDUNE)

if(NOT CPU_ONLY)
  set_target_properties(SamplesDUNE PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
endif()

target_link_libraries(SamplesDUNE PUBLIC SplinesDUNE duneanaobj::all MaCh3::All MaCh3DUNECompilerOptions)
target_link_libraries(SamplesDUNE PRIVATE DUNEMaCh3Warnings)

target_include_directories(SamplesDUNE PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
  $<INSTALL_INTERFACE:include>)

install(TARGETS SamplesDUNE
  EXPORT mach3dune-targets
  LIBRARY DESTINATION lib/
  PUBLIC_HEADER DESTINATION include/Samples)

add_library(MACH3DUNE::SamplesDUNE ALIAS SamplesDUNE)
