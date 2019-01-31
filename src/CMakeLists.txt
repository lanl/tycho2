target_sources(tycho2 INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/tycho2.hpp>
  $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include/tycho2.hpp>)

add_subdirectory(tycho2)
