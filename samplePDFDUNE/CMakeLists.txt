set(HEADERS
#  samplePDFDUNEAtm.h
#  samplePDFDUNEBeamFD.h
#  samplePDFDUNEBeamND.h
  samplePDFDUNEBeamNDGar.h
  MaCh3DUNEFactory.h 
)

add_library(SamplePDFDUNE SHARED
#  samplePDFDUNEAtm.cpp
#  samplePDFDUNEBeamFD.cpp
#  samplePDFDUNEBeamND.cpp
  samplePDFDUNEBeamNDGar.cpp
  MaCh3DUNEFactory.cpp
)

set_target_properties(SamplePDFDUNE PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME SamplePDFDUNE)

if(NOT CPU_ONLY)
  set_target_properties(SamplePDFDUNE PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
endif()

target_link_libraries(SamplePDFDUNE PUBLIC splinesDUNE duneanaobj::all MaCh3::All MaCh3DUNECompilerOptions)
target_link_libraries(SamplePDFDUNE PRIVATE DUNEMaCh3Warnings)

target_include_directories(SamplePDFDUNE PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
  $<INSTALL_INTERFACE:include>)

install(TARGETS SamplePDFDUNE
        EXPORT mach3dune-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/samplePDF)

add_library(MACH3DUNE::SamplePDFDUNE ALIAS SamplePDFDUNE)
