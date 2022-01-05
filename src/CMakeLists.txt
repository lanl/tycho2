
option(USE_PETSC "Use PETSc" OFF)
option(USE_NVGRAPH "Use nvGraph" ON)
option(ASSERT_ON "Turn on asserts" ON)

function(bool_to_int OUT IN)
  if (${IN})
    set(${OUT} "1" PARENT_SCOPE)
  else()
    set(${OUT} "0" PARENT_SCOPE)
  endif()
endfunction(bool_to_int)

add_executable(tycho2.x 
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
target_compile_definitions(tycho2.x
    PUBLIC "-DASSERT_ON=${ASSERT_ON_INT}" "-DUSE_PETSC=${USE_PETSC_INT}")

target_include_directories(tycho2.x PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

target_include_directories(tycho2.x PUBLIC ${MPI_INCLUDE_PATH})

target_link_libraries(tycho2.x PUBLIC ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(tycho2.x PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(tycho2.x PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

message("USE_PETSC ${USE_PETSC} (${USE_PETSC_INT})")
if (USE_PETSC)
  set(PETSC_DIR "" CACHE PATH "PETSc install directory")
  target_include_directories(tycho2.x PUBLIC "${PETSC_DIR}/include")
  find_library(PETSC_LIBRARY NAMES petsc PATHS "${PETSC_DIR}/lib" NO_DEFAULT_PATH)
  target_link_libraries(tycho2.x PUBLIC "${PETSC_LIBRARY}")
endif()

message("USE_NVGRAPH ${USE_NVGRAPH}")
if (USE_NVGRAPH)
  string(REPLACE ":" ";" CMAKE_PREFIX_PATH "$ENV{CMAKE_PREFIX_PATH}" )
  list(APPEND CMAKE_PREFIX_PATH "/usr/tce/packages/cuda/cuda-10.1.243/lib64")
  list(APPEND CMAKE_PREFIX_PATH "/usr/tce/packages/cuda/cuda-10.1.243/include")
  find_library(NVGLIB nvgraph)
  target_link_libraries(tycho2.x PUBLIC "${NVGLIB}" )
  find_path(NVGINC nvgraph.h)
  target_include_directories(tycho2.x PUBLIC "${NVGINC}" )
endif()

