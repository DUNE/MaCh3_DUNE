add_custom_target(MaCh3DUNETests)

# Wet up utils lib
set(HEADERS
    Comparison.h
)

add_library(
    MaCh3DuneTestUtils SHARED
    Comparison.cpp
)

set_target_properties(MaCh3DuneTestUtils PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
    EXPORT_NAME MaCh3DuneTestUtils
)

target_link_libraries(MaCh3DuneTestUtils PUBLIC MaCh3::All MaCh3DUNECompilerOptions)

target_include_directories(MaCh3DuneTestUtils PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../>
  $<INSTALL_INTERFACE:include>)

  install(TARGETS MaCh3DuneTestUtils
  EXPORT mach3dune-targets
  LIBRARY DESTINATION lib/
  PUBLIC_HEADER DESTINATION include/MaCh3DuneTestUtils)

add_library(MACH3DUNE::MaCh3DuneTestUtils ALIAS MaCh3DuneTestUtils)

# Now we set up test apps
foreach(app
    TestEventRates
)

    add_executable(${app} ${app}.cpp)
    target_link_libraries(${app} MaCh3DUNE::All MaCh3DuneTestUtils Minuit)
    add_dependencies(MaCh3DUNEApps ${app})
    install(TARGETS ${app} DESTINATION tests)

endforeach(app)


# Now copy remaining stuff
install(DIRECTORY test_results DESTINATION tests)
