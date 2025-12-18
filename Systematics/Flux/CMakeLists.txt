add_library(OffAxisFluxUncertaintyHelper SHARED OffAxisFluxUncertaintyHelper.cxx)
target_link_libraries(OffAxisFluxUncertaintyHelper PUBLIC ROOT::Hist)

set_target_properties(OffAxisFluxUncertaintyHelper PROPERTIES
    PUBLIC_HEADER "OffAxisFluxUncertaintyHelper.h"
    EXPORT_NAME OffAxisFluxUncertaintyHelper)

target_include_directories(OffAxisFluxUncertaintyHelper PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../../>
  $<INSTALL_INTERFACE:include>)

install(TARGETS OffAxisFluxUncertaintyHelper
        EXPORT mach3dune-targets
        LIBRARY DESTINATION lib/
        PUBLIC_HEADER DESTINATION include/Systematics/Flux)