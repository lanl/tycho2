#include <moab/Core.hpp>
#include <string>
#include <cstdio>
#include <cassert>
#include "ParallelMesh.hh"
#include <iostream>
#include <fstream>

using namespace moab;
using namespace std;


/*
    moabCreateAllFaces
*/
static
void createMoabMesh(Core &moab, const EntityHandle &meshSet, 
                    const ParallelMesh &mesh, const uint64_t part)
{
    vector<EntityHandle> cellHandles;
    ErrorCode moabStatus;
    const ParallelMesh::PartitionData &partData = mesh.c_partitionData[part];
    
    
    // Go through cells in mesh and create tets for Moab
    for (uint64_t cell = 0; cell < partData.numCells; cell++) {
        EntityHandle vrtxHandles[4];
        for (int vrtx = 0; vrtx < 4; vrtx++) {
            uint64_t node = partData.cellData[cell].boundingNodes[vrtx];
            moabStatus = moab.create_vertex(partData.nodeData[node].coords, 
                                            vrtxHandles[vrtx]);
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
    getFilePrefix
*/
string getFilePrefix(const string &inputFile)
{
    int length = inputFile.length();
    assert(inputFile.substr(length - 6) == ".pmesh");
    return inputFile.substr(0, length - 6);
}


/*
    writeVisitFile
*/
void writeVisitFile(const string &prefix, const vector<string> &filenames)
{
    ofstream file;
    file.open(prefix + ".visit");
    file << "!NBLOCKS " << filenames.size() << endl;
    
    for (size_t i = 0; i < filenames.size(); i++) {
        file << filenames[i] << endl;
    }
    
    file.close();
}


/*
    main
*/
int main(int argc, char* argv[])
{
    Core moab;
    string inputFile;
    string outputFilePrefix;
    ParallelMesh mesh;
    ErrorCode moabStatus;
    vector<string> filenames;
    
    
    // Print utility name
    printf("--- ParallelMeshToMoab Utility ---\n");
    
    
    // Get input/output files
    if (argc != 2) {
        printf("Incorrect number of arguments\n");
        printf("Usage: ./ParallelMeshToMoab.x <inputFile>\n");
        printf("\n\n\n");
        return 0;
    }
    inputFile = argv[1];
    outputFilePrefix = getFilePrefix(inputFile);
    
    
    // Convert mesh to moab
    printf("Reading parallel mesh from file\n");
    mesh.read(inputFile);
    
    for (uint64_t part = 0; part < mesh.c_numPartitions; part++) {
        EntityHandle meshSet;
        
        moabStatus = moab.create_meshset(MESHSET_SET, meshSet);
        assert(moabStatus == MB_SUCCESS);
        
        printf("Creating Moab mesh for partition %" PRIu64 "\n", part);
        createMoabMesh(moab, meshSet, mesh, part);
        
        printf("   Writing Moab mesh to file\n");
        string outputFile = outputFilePrefix;
        outputFile += "-";
        outputFile += to_string(part);
        outputFile += "p.vtk";
        moabStatus = moab.write_mesh(outputFile.c_str(), &meshSet, 1);
        assert(moabStatus == MB_SUCCESS);
        
        filenames.push_back(outputFile);
    }
    
    
    // Write main file
    printf("Writing .visit file\n");
    writeVisitFile(outputFilePrefix, filenames);
    
    
    // Cleanup
    printf("\n\n\n");
    return 0;
}