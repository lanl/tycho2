target_sources(tycho2 INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/GaussLegendre.hpp>
  $<INSTALL_INTERFACE:
    $<INSTALL_PREFIX>/include/tycho2/quadrature/GaussLegendre.hpp>)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/GaussLegendre)
