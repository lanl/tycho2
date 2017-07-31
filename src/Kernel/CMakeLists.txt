set(KERNEL_SOURCES
  MainKernel.cc
  SweeperKernel.cc
  GraphTraverserKernel.cc
  UtilKernel.cc
  TransportKernel.cc
  ../Global.cc
  ../Assert.cc
  ../Comm.cc
  ../Quadrature.cc
  ../TychoMesh.cc
  ../TychoMeshIO.cc
  ../ParallelMesh.cc
  )

add_executable(kernel.x ${KERNEL_SOURCES})
target_include_directories(kernel.x PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/..")
target_link_libraries(kernel.x PUBLIC "${PETSC_LIBRARY}")
target_include_directories(kernel.x PUBLIC "${KokkosCore_INCLUDE_DIRS}")
target_link_libraries(kernel.x PUBLIC kokkoscore)
