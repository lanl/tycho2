#include <string>
#include <cstdio>
#include <cassert>
#include <mpi.h>
#include "Comm.hh"
#include "ParallelMesh.hh"

using namespace std;


/*
    main
*/
int main(int argc, char* argv[])
{
    ParallelMesh parallelMesh;
    ParallelMesh::PartitionData partData;
    string inputFile;
    bool verbose = false;
    
    
    // Start MPI
    MPI_Init(&argc, &argv);


    // Print utility name
    if (Comm::rank() == 0) {
        printf("--- ParallelMeshInfo Utility ---\n");
    }
    

    // Get input file
    if (argc == 2) {
        verbose = false;
    }
    else if (argc == 3 && strcmp(argv[2], "verbose") == 0) {
        verbose = true;
    }
    else {
        if (Comm::rank() == 0) {
            printf("Incorrect number of arguments\n");
            printf("Usage: ./ParallelMeshInfo.x <inputFile> (verbose)\n");
            printf("\n\n\n");
        }
        return 0;
    }
    inputFile = argv[1];
    
    
    // Read in parallel mesh serially and print info
    if (Comm::rank() == 0) {
        printf("Parallel mesh data read serially.\n");
        parallelMesh.read(inputFile);
        parallelMesh.print(verbose);
    }
    Comm::barrier();


    // Read parallel mesh in parallel and print info
    if (Comm::rank() == 0) {
        printf("\n\n\nParallel mesh data read in parallel.\n");
    }
    ParallelMesh::readInParallel(inputFile, partData);
    for (int part = 0; part < Comm::numRanks(); part++) {
        if (Comm::rank() == part) {
            printf("Partition %d\n", Comm::rank());
            ParallelMesh::printPartitionData(partData, verbose);
        }
        Comm::barrier();
    }

    
    // Cleanup
    Comm::barrier();
    if (Comm::rank() == 0) {
        printf("\n\n\n");
    }
    MPI_Finalize();
    return 0;
}



