cmake_minimum_required(VERSION 3.12.1)

project(Tycho2 VERSION 0.1 LANGUAGES CXX)

project( TYCHO2 )

find_package(MPI REQUIRED)

set_property(GLOBAL APPEND PROPERTY CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-O3 -fopenmp -Wall -Wextra -Wpedantic ")

add_subdirectory(src)

option(BUILD_UTILS "Build utility programs" ON)
if (BUILD_UTILS)
  add_subdirectory(util)
endif()

