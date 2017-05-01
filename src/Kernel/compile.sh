SRC="MainKernel.cc SweeperKernel.cc GraphTraverserKernel.cc UtilKernel.cc TransportKernel.cc ../Global.cc ../Assert.cc ../Comm.cc ../Quadrature.cc ../TychoMesh.cc ../TychoMeshIO.cc ../ParallelMesh.cc"

mpicxx -Wall -Wextra -g -std=c++11 -O3 -fopenmp -I.. $SRC
