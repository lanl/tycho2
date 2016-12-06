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
