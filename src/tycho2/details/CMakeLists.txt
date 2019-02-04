target_sources(tycho2 INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/make_array.hpp>
  $<INSTALL_INTERFACE:
    $<INSTALL_PREFIX>/include/tycho2/details/make_array.hpp>)

