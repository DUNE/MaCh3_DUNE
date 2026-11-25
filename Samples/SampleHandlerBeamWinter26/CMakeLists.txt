FetchContent_Declare(
  duneanasel
  GIT_REPOSITORY https://github.com/DUNE/duneanasel.git
  GIT_TAG        ana/MaCh3Winter26
)
FetchContent_MakeAvailable(duneanasel)

set(SOURCE_FILES
    ReadEvents.cpp
    SampleHandlerBeamWinter26.cpp
)
set(HEADERS
  EventInfo.h
  Projections.h
  ReadEvents.h
  SampleHandlerBeamWinter26.h
)

add_library(SampleHandlerBeamWinter26 SHARED
  ${SOURCE_FILES}
)

set_target_properties(SampleHandlerBeamWinter26 PROPERTIES
  PUBLIC_HEADER "${HEADERS}"
  EXPORT_NAME SampleHandlerBeamWinter26)

target_include_directories(SampleHandlerBeamWinter26 PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
  $<INSTALL_INTERFACE:include>)

install(TARGETS SampleHandlerBeamWinter26
  EXPORT mach3dune-targets
  LIBRARY DESTINATION lib/
  PUBLIC_HEADER DESTINATION include/Samples/beamStandardRecord)

target_link_libraries(SampleHandlerBeamWinter26 PUBLIC SplinesDUNE MaCh3::All MaCh3DUNECompilerOptions)
target_link_libraries(SampleHandlerBeamWinter26 PRIVATE duneanaobj::StandardRecordProxy duneanasel::all DUNEMaCh3Warnings ROOT::ROOTNTuple)

add_executable(BeamWinter26BinSpliner BeamWinter26BinSpliner.cpp)
target_link_libraries(BeamWinter26BinSpliner PRIVATE SamplesDUNE SampleHandlerBeamWinter26)
