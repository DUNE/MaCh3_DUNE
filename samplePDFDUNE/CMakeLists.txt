if (Build_NDGAr)
	message(STATUS "Building BeamNDGAr sample only")
	set(HEADERS
		samplePDFDUNEBeamNDGAr.h
  	MaCh3DUNEFactory.h
	)
	set(SOURCE_FILES
		samplePDFDUNEBeamNDGAr.cpp
    MaCh3DUNEFactory.cpp
	)
else()
	message(STATUS "Building Atm, BeamFD and BeamND samples")
	set(HEADERS
		samplePDFDUNEAtm.h
  	samplePDFDUNEBeamFD.h
   	samplePDFDUNEBeamND.h
	samplePDFDUNEBeamFD_actuallyFD.h
		MaCh3DUNEFactory.h
	)
	set(SOURCE_FILES
		samplePDFDUNEAtm.cpp
    samplePDFDUNEBeamFD.cpp
    samplePDFDUNEBeamND.cpp
	samplePDFDUNEBeamFD_actuallyFD.cpp
    MaCh3DUNEFactory.cpp
	)
endif()

add_library(SamplePDFDUNE SHARED
	${SOURCE_FILES}
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
