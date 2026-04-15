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
set(HEADERS_PDSP
  SampleHandlerPDSP.h
)
set(SOURCE_FILES_PDSP
  SampleHandlerPDSP.cpp
)

add_library(SamplesDUNE SHARED
  ${SOURCE_FILES}
)
add_library(SamplesPDSP SHARED
  ${SOURCE_FILES_PDSP}
)

set_target_properties(SamplesDUNE PROPERTIES
  PUBLIC_HEADER "${HEADERS}"
  EXPORT_NAME SamplesDUNE)
set_target_properties(SamplesPDSP PROPERTIES
  PUBLIC_HEADER "${HEADERS_PDSP}"
  EXPORT_NAME SamplesPDSP)

if(NOT CPU_ONLY)
  set_target_properties(SamplesDUNE PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
  set_target_properties(SamplesPDSP PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
endif()

target_link_libraries(SamplesDUNE PUBLIC SplinesDUNE duneanaobj::all MaCh3::All MaCh3DUNECompilerOptions)
target_link_libraries(SamplesDUNE PRIVATE DUNEMaCh3Warnings)

target_link_libraries(SamplesPDSP PUBLIC SplinesDUNE duneanaobj::all MaCh3::All MaCh3DUNECompilerOptions)
target_link_libraries(SamplesPDSP PRIVATE DUNEMaCh3Warnings)

target_include_directories(SamplesDUNE PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
  $<INSTALL_INTERFACE:include>)
target_include_directories(SamplesPDSP PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
  $<INSTALL_INTERFACE:include>)

install(TARGETS SamplesDUNE
  EXPORT mach3dune-targets
  LIBRARY DESTINATION lib/
  PUBLIC_HEADER DESTINATION include/Samples)
install(TARGETS SamplesPDSP
  EXPORT mach3dune-targets
  LIBRARY DESTINATION lib/
  PUBLIC_HEADER DESTINATION include/Samples)

add_library(MACH3DUNE::SamplesDUNE ALIAS SamplesDUNE)
add_library(MACH3DUNE::SamplesPDSP ALIAS SamplesPDSP)
