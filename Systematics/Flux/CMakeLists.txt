add_library(FluxSystematics SHARED OffAxisFluxUncertaintyHelper.cxx)
target_link_libraries(FluxSystematics PUBLIC ROOT::Hist)

set_target_properties(FluxSystematics PROPERTIES
    PUBLIC_HEADER "OffAxisFluxUncertaintyHelper.h"
    EXPORT_NAME OffAxisFluxUncertaintyHelper)

target_include_directories(FluxSystematics PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../../>
  $<INSTALL_INTERFACE:include>)

install(TARGETS FluxSystematics
        EXPORT mach3dune-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/Systematics/Flux)
install(FILES flux_variations_FD_and_PRISM_2023.root DESTINATION Systematics/Flux)
