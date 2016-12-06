/*
Copyright (c) 2016, Los Alamos National Security, LLC
All rights reserved.

Copyright 2016. Los Alamos National Security, LLC. This software was produced 
under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National 
Laboratory (LANL), which is operated by Los Alamos National Security, LLC for 
the U.S. Department of Energy. The U.S. Government has rights to use, 
reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS 
ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR 
ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified 
to produce derivative works, such modified software should be clearly marked, 
so as not to confuse it with the version available from LANL.

Additionally, redistribution and use in source and binary forms, with or 
without modification, are permitted provided that the following conditions 
are met:
1.      Redistributions of source code must retain the above copyright notice, 
        this list of conditions and the following disclaimer.
2.      Redistributions in binary form must reproduce the above copyright 
        notice, this list of conditions and the following disclaimer in the 
        documentation and/or other materials provided with the distribution.
3.      Neither the name of Los Alamos National Security, LLC, Los Alamos 
        National Laboratory, LANL, the U.S. Government, nor the names of its 
        contributors may be used to endorse or promote products derived from 
        this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND 
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT 
NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL 
SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED 
OF THE POSSIBILITY OF SUCH DAMAGE.
*/

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
