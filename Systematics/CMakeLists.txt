add_subdirectory(Flux)

add_library(SystematicsDUNE INTERFACE)
target_link_libraries(SystematicsDUNE INTERFACE FluxSystematics)

install(TARGETS SystematicsDUNE
  EXPORT mach3dune-targets)
