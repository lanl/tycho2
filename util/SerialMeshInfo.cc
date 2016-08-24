#include <string>
#include <cstdio>
#include <cassert>
#include <cstring>
#include "SerialMesh.hh"

using namespace std;


/*
    main
*/
int main(int argc, char* argv[])
{
    SerialMesh serialMesh;
    string inputFile;
    bool verbose = false;
    
    
    // Print utility name
    printf("--- SerialMeshInfo Utility ---\n");
    
    
    // Get input file
    if (argc == 2) {
        verbose = false;
    }
    else if (argc == 3 && strcmp(argv[2], "verbose") == 0) {
        verbose = true;
    }
    else {
        printf("Incorrect number of arguments\n");
        printf("Usage: ./SerialMeshInfo.x <inputFile> (verbose)\n");
        printf("\n\n\n");
        return 0;
    }
    inputFile = argv[1];
    
    
    // Read in serial mesh and convert to parallel mesh
    serialMesh.read(inputFile);
    serialMesh.print(verbose);
    
    
    // Cleanup
    printf("\n\n\n");
    return 0;
}



