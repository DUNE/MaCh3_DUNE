set(HEADERS splinesDUNE.h
)

add_library(splinesDUNE SHARED
    splinesDUNE.cpp
)

set_target_properties(splinesDUNE PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME splinesDUNE)

target_link_libraries(splinesDUNE MaCh3::All MaCh3DUNECompilerOptions)

target_include_directories(splinesDUNE PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
  $<INSTALL_INTERFACE:include>)

install(TARGETS splinesDUNE
        EXPORT mach3dune-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/splines)

add_library(MACH3DUNE::splinesDUNE ALIAS splinesDUNE)
