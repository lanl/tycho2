cmake_minimum_required(VERSION 3.0.1)

project(Tycho2 VERSION 0.1 LANGUAGES CXX)

option(ASSERT_ON "Turn on asserts" ON)
option(USE_PETSC "Use PETSc" OFF)
option(USE_KOKKOS "Use Kokkos" ON)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCES
  src/Assert.cc
  src/Comm.cc
  src/CommSides.cc
  src/Global.cc
  src/GraphTraverser.cc
  src/KeyValueReader.cc
  src/Main.cc
  src/ParallelMesh.cc
  src/Problem.cc
  src/PsiData.cc
  src/Quadrature.cc
  src/SerialMesh.cc
  src/Sweeper.cc
  src/Transport.cc
  src/TychoMesh.cc
  src/TychoMeshIO.cc
  src/Util.cc
  )

add_executable(sweep.x ${SOURCES})

function(bool_to_int OUT IN)
  if (${IN})
    set(${OUT} "1" PARENT_SCOPE)
  else()
    set(${OUT} "0" PARENT_SCOPE)
  endif()
endfunction(bool_to_int)

bool_to_int(ASSERT_ON_INT ${ASSERT_ON})
bool_to_int(USE_PETSC_INT ${USE_PETSC})

target_compile_definitions(sweep.x
    PUBLIC "-DASSERT_ON=${ASSERT_ON_INT}" "-DUSE_PETSC=${USE_PETSC_INT}")
target_include_directories(sweep.x
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")

if (USE_PETSC)
  set(PETSC_DIR "" CACHE PATH "PETSc install directory")
  target_include_directories(sweep.x PUBLIC "${PETSC_DIR}/include")
  find_library(PETSC_LIBRARY NAMES petsc PATHS "${PETSC_DIR}/lib" NO_DEFAULT_PATH)
  target_link_libraries(sweep.x PUBLIC "${PETSC_LIBRARY}")
endif()

option(BUILD_UTILS "Build utility programs" ON)
if (BUILD_UTILS)
  add_subdirectory(util)
endif()

message("USE_KOKKOS ${USE_KOKKOS}")
if (USE_KOKKOS)
  set(KOKKOS_PREFIX "" CACHE PATH "Kokkos install directory")
  find_package(KokkosCore REQUIRED PATHS ${KOKKOS_PREFIX} NO_DEFAULT_PATH)
  message("KokkosCore_CONFIG ${KokkosCore_CONFIG}")
  message("KokkosCore_INCLUDE_DIRS ${KokkosCore_INCLUDE_DIRS}")
  target_link_libraries(sweep.x PUBLIC kokkoscore)
  target_include_directories(sweep.x PUBLIC "${KokkosCore_INCLUDE_DIRS}")
endif()
