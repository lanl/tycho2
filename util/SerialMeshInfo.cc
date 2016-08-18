#include <string>
#include <cstdio>
#include <cassert>
#include "SerialMesh.hh"

using namespace std;


/*
    main
*/
int main(int argc, char* argv[])
{
    SerialMesh serialMesh;
    string inputFile;
    
    
    // Print utility name
    printf("--- SerialMeshInfo Utility ---\n");
    
    
    // Get input file
    if (argc != 2) {
        printf("Incorrect number of arguments\n");
        printf("Usage: ./SerialMeshInfo.x <inputFile>\n");
        printf("\n\n\n");
        return 0;
    }
    inputFile = argv[1];
    
    
    // Read in serial mesh and convert to parallel mesh
    serialMesh.read(inputFile);
    serialMesh.printAll();
    
    
    // Cleanup
    printf("\n\n\n");
    return 0;
}



