
option(USE_PETSC "Use PETSc" OFF)
option(USE_KOKKOS "Use Kokkos" ON)
option(ASSERT_ON "Turn on asserts" ON)

function(bool_to_int OUT IN)
  if (${IN})
    set(${OUT} "1" PARENT_SCOPE)
  else()
    set(${OUT} "0" PARENT_SCOPE)
  endif()
endfunction(bool_to_int)

add_executable(sweep.x 
Assert.cc
Comm.cc
CommSides.cc
Global.cc
GraphTraverser.cc
KeyValueReader.cc
Main.cc
ParallelMesh.cc
Problem.cc
PsiData.cc
Quadrature.cc
SerialMesh.cc
Sweeper.cc
TychoMesh.cc
TychoMeshIO.cc
Util.cc)

bool_to_int(USE_PETSC_INT ${USE_PETSC})
bool_to_int(ASSERT_ON_INT ${ASSERT_ON})
target_compile_definitions(sweep.x
    PUBLIC "-DASSERT_ON=${ASSERT_ON_INT}" "-DUSE_PETSC=${USE_PETSC_INT}")

target_include_directories(sweep.x PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

target_include_directories(sweep.x PUBLIC ${MPI_INCLUDE_PATH})
target_link_libraries(sweep.x PUBLIC ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(sweep.x PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(sweep.x PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

message("USE_PETSC ${USE_PETSC} (${USE_PETSC_INT})")
if (USE_PETSC)
  set(PETSC_DIR "" CACHE PATH "PETSc install directory")
  target_include_directories(sweep.x PUBLIC "${PETSC_DIR}/include")
  find_library(PETSC_LIBRARY NAMES petsc PATHS "${PETSC_DIR}/lib" NO_DEFAULT_PATH)
  target_link_libraries(sweep.x PUBLIC "${PETSC_LIBRARY}")
endif()

message("USE_KOKKOS ${USE_KOKKOS}")
if (USE_KOKKOS)
   ##list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
   SET(KOKKOS_DIR "$ENV{KOKKOS_ROOT}" CACHE LOCATION "where in the world is Kokkos Miranda?")

   message(" ... looking for ${KOKKOS_DIR}")

   find_path(
   KOKKOS_INCLUDE_PATH
   NAMES "include"
   PATHS "${KOKKOS_DIR}"
   NO_DEFAULT_PATH)

   set(KOKKOS_INCLUDE_DIRS "${KOKKOS_INCLUDE_PATH}/include")

   find_library(
   KOKKOS_LIBRARIES
   NAMES "libkokkos.a" "libkokkos.so"
   PATHS "${KOKKOS_DIR}/lib"
   NO_DEFAULT_PATH )

   target_include_directories(sweep.x PUBLIC "${KOKKOS_INCLUDE_DIRS}")
   target_link_libraries(sweep.x PUBLIC ${KOKKOS_LIBRARIES})
endif()

