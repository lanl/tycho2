#include <moab/Core.hpp>
#include <string>
#include <cstdio>
#include <cassert>
#include "SerialMesh.hh"

using namespace moab;
using namespace std;


/*
    moabCreateAllFaces
*/
static
void createMoabMesh(Core &moab, const EntityHandle &meshSet, const SerialMesh &mesh)
{
    vector<EntityHandle> cellHandles;
    ErrorCode moabStatus;
    
    
    // Go through cells in serial mesh and create tets for Moab
    for (uint64_t cell = 0; cell < mesh.c_numCells; cell++) {
        EntityHandle vrtxHandles[4];
        for (int vrtx = 0; vrtx < 4; vrtx++) {
            uint64_t node = mesh.c_cellData[cell].boundingNodes[vrtx];
            moabStatus = moab.create_vertex(mesh.c_nodeData[node].coords, vrtxHandles[vrtx]);
            assert(moabStatus == MB_SUCCESS);
        }
        
        EntityHandle cellHandle;
        moabStatus = moab.create_element(MBTET, vrtxHandles, 4, cellHandle);
        assert(moabStatus == MB_SUCCESS);
        cellHandles.push_back(cellHandle);
    }
    
    
    // Add all the tets to Moab mesh
    moabStatus = moab.add_entities(meshSet, cellHandles.data(), cellHandles.size());
    assert(moabStatus == MB_SUCCESS);
}


/*
    getOutputFile
*/
string getOutputFile(const string &inputFile)
{
    int length = inputFile.length();
    assert(inputFile.substr(length - 6) == ".smesh");
    return inputFile.substr(0, length - 6) + ".vtk";
}


/*
    main
*/
int main(int argc, char* argv[])
{
    Core moab;
    string inputFile;
    string outputFile;
    SerialMesh mesh;
    ErrorCode moabStatus;
    EntityHandle meshSet;

    
    // Print utility name
    printf("--- SerialMeshToMoab Utility ---\n");
    
    
    // Get input/output files
    if (argc != 2) {
        printf("Incorrect number of arguments\n");
        printf("Usage: ./SerialMeshToMoab.x <inputFile>\n");
        printf("\n\n\n");
        return 0;
    }
    inputFile = argv[1];
    outputFile = getOutputFile(inputFile);
    
    
    // Convert mesh to moab
    moabStatus = moab.create_meshset(MESHSET_SET, meshSet);
    assert(moabStatus == MB_SUCCESS);
    
    printf("Reading serial mesh from file\n");
    mesh.read(inputFile);
    
    printf("Creating Moab mesh\n");
    createMoabMesh(moab, meshSet, mesh);
    
    printf("Writing Moab mesh to file\n");
    moabStatus = moab.write_mesh(outputFile.c_str(), &meshSet, 1);
    assert(moabStatus == MB_SUCCESS);
    
    
    // Cleanup
    printf("\n\n\n");
    return 0;
}