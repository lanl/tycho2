target_sources(tycho2 INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/quadrature.hpp>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include/tycho2/quadrature.hpp>)

add_subdirectory(quadrature)
