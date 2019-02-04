target_sources(tycho2 INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/Tetrahedral.hpp>
  $<INSTALL_INTERFACE:
    $<INSTALL_PREFIX>/include/tycho2/quadrature/Tetrahedral.hpp>)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/Tetrahedral)
