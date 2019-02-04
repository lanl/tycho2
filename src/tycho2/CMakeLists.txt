target_sources(tycho2 INTERFACE
  $<BUILD_INTERFACE:
    ${CMAKE_CURRENT_SOURCE_DIR}/quadrature.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mesh.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/details.hpp>
  $<INSTALL_INTERFACE:
    $<INSTALL_PREFIX>/include/tycho2/quadrature.hpp
    $<INSTALL_PREFIX>/include/tycho2/mesh.hpp
    $<INSTALL_PREFIX>/include/tycho2/details.hpp>)

#add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/details)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/quadrature)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/mesh)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/io)

