set(HEADERS BinnedSplineHandlerDUNE.h
)

add_library(SplinesDUNE SHARED
    BinnedSplineHandlerDUNE.cpp
)

set_target_properties(SplinesDUNE PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME SplinesDUNE)

target_link_libraries(SplinesDUNE MaCh3::All MaCh3DUNECompilerOptions)

target_include_directories(SplinesDUNE PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
  $<INSTALL_INTERFACE:include>)

install(TARGETS SplinesDUNE
        EXPORT mach3dune-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/Splines)

add_library(MACH3DUNE::SplinesDUNE ALIAS SplinesDUNE)
