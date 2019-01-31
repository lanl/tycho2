cmake_minimum_required(VERSION 3.11.1)

if(DEFINED PROJECT_NAME)
  set(tycho2.tests OFF)
else()
  set(tycho2.tests ON)
endif()

project(tycho2 LANGUAGES CXX VERSION 1.0)

include(CTest)

project(tycho2 VERSION 0.0.0)

add_library(tycho2 INTERFACE)
add_library(tycho2::tycho2 ALIAS tycho2)
target_compile_features(tycho2 INTERFACE cxx_std_14)
target_include_directories(tycho2 INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
  $<INSTALL_INTERFACE:include>)

if(tycho2.tests)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/dependencies/Catch2)
endif()

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src/tycho2)

