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

#include "SerialMesh.hh"
#include <cstdio>
#include <cassert>
#include <cinttypes>
#include <cstring>

using namespace std;

static const char MESH_FORMAT_NAME[SerialMesh::MESH_FORMAT_NAME_LEN] = 
    {'T', 'y', 'c', 'h', 'o', ' ', '2', ' ', 'S', 'e', 'r', 'i', 'a', 'l', 
     ' ', 'M', 'e', 's', 'h', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};



/*
    printAll
    
    Prints the entire mesh data.
*/
void SerialMesh::print(bool printVerbose)
{
    // Header Data
    printf("Mesh Format Name: %s\n", c_meshFormatName);
    printf("Mesh Version: %" PRIu64 "\n", c_version);
    printf("Num Cells: %" PRIu64 "\n", c_numCells);
    printf("Num Faces: %" PRIu64 "\n", c_numFaces);
    printf("Num Nodes: %" PRIu64 "\n", c_numNodes);
    
    
    // Exit if not verbose print
    if (!printVerbose)
        return;


    // Cell Data
    for (uint64_t cell = 0; cell < c_numCells; cell++) {
        printf("Cell %" PRIu64 "\n", cell);
        printf("   Faces %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 "\n", 
               c_cellData[cell].boundingFaces[0],
               c_cellData[cell].boundingFaces[1],
               c_cellData[cell].boundingFaces[2],
               c_cellData[cell].boundingFaces[3]);
        printf("   Nodes %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 "\n", 
               c_cellData[cell].boundingNodes[0],
               c_cellData[cell].boundingNodes[1],
               c_cellData[cell].boundingNodes[2],
               c_cellData[cell].boundingNodes[3]);
        printf("   Material Index %" PRIu64 "\n", 
               c_cellData[cell].materialIndex);
    }
    
    
    // Face Data
    for (uint64_t face = 0; face < c_numFaces; face++) {
        printf("Face %" PRIu64 "\n", face);
        printf("   Cells %" PRIu64 " %" PRIu64 "\n", 
               c_faceData[face].boundingCells[0],
               c_faceData[face].boundingCells[1]);
        printf("   Nodes %" PRIu64 " %" PRIu64 " %" PRIu64 "\n", 
               c_faceData[face].boundingNodes[0],
               c_faceData[face].boundingNodes[1],
               c_faceData[face].boundingNodes[2]);
    }
    
    
    // Node Data
    for (uint64_t node = 0; node < c_numNodes; node++) {
        printf("Node %" PRIu64 ": (%f, %f, %f)\n", node, 
               c_nodeData[node].coords[0],
               c_nodeData[node].coords[1],
               c_nodeData[node].coords[2]);
    }
}


/*
    write
    
    Write serial mesh to file.
*/
void SerialMesh::write(const std::string &filename)
{
    FILE *file;
    size_t numWritten;
    vector<uint64_t> bufferUint;
    vector<double> bufferDouble;
    
    
    // Open file
    file = fopen(filename.c_str(), "wb");
    assert(file != NULL);


    // Make sure format name and version are filled in
    memcpy(c_meshFormatName, MESH_FORMAT_NAME, 
           SerialMesh::MESH_FORMAT_NAME_LEN * sizeof(char));
    c_version = SerialMesh::VERSION;
    
    
    // Write mesh format name
    numWritten = fwrite(c_meshFormatName, sizeof(char), 
                        SerialMesh::MESH_FORMAT_NAME_LEN, file);
    assert(numWritten == SerialMesh::MESH_FORMAT_NAME_LEN);
    

    // Buffer header data
    bufferUint.push_back(c_version);
    bufferUint.push_back(c_numCells);
    bufferUint.push_back(c_numFaces);
    bufferUint.push_back(c_numNodes);

    
    // Buffer cell data
    for (uint64_t cell = 0; cell < c_numCells; cell++) {
        bufferUint.push_back(c_cellData[cell].boundingFaces[0]);
        bufferUint.push_back(c_cellData[cell].boundingFaces[1]);
        bufferUint.push_back(c_cellData[cell].boundingFaces[2]);
        bufferUint.push_back(c_cellData[cell].boundingFaces[3]);
        
        bufferUint.push_back(c_cellData[cell].boundingNodes[0]);
        bufferUint.push_back(c_cellData[cell].boundingNodes[1]);
        bufferUint.push_back(c_cellData[cell].boundingNodes[2]);
        bufferUint.push_back(c_cellData[cell].boundingNodes[3]);

        bufferUint.push_back(c_cellData[cell].materialIndex);
    }
    
    
    // Buffer face data
    for (uint64_t face = 0; face < c_numFaces; face++) {
        bufferUint.push_back(c_faceData[face].boundingCells[0]);
        bufferUint.push_back(c_faceData[face].boundingCells[1]);
        
        bufferUint.push_back(c_faceData[face].boundingNodes[0]);
        bufferUint.push_back(c_faceData[face].boundingNodes[1]);
        bufferUint.push_back(c_faceData[face].boundingNodes[2]);
    }


    // Buffer node data
    for (uint64_t node = 0; node < c_numNodes; node++) {
        bufferDouble.push_back(c_nodeData[node].coords[0]);
        bufferDouble.push_back(c_nodeData[node].coords[1]);
        bufferDouble.push_back(c_nodeData[node].coords[2]);
    }
    
    
    // Write uint64_t buffered data
    numWritten = fwrite(bufferUint.data(), sizeof(uint64_t), 
                        bufferUint.size(), file);
    assert(numWritten == bufferUint.size());

    
    // Write double buffered data
    numWritten = fwrite(bufferDouble.data(), sizeof(double), 
                        bufferDouble.size(), file);
    assert(numWritten == bufferDouble.size());

    
    // Close file
    fclose(file);
}


/*
    read
    
    Read serial mesh from file.
*/
void SerialMesh::read(const std::string &filename)
{
    FILE *file;
    size_t numRead;
    vector<uint64_t> bufferUint;
    vector<double> bufferDouble;
    
    
    // Open file
    file = fopen(filename.c_str(), "rb");
    assert(file != NULL);
    
    
    // Read mesh format name
    numRead = fread(c_meshFormatName, sizeof(char), 
                    SerialMesh::MESH_FORMAT_NAME_LEN, file);
    assert(numRead == SerialMesh::MESH_FORMAT_NAME_LEN);
    assert(memcmp(c_meshFormatName, MESH_FORMAT_NAME, 
                  SerialMesh::MESH_FORMAT_NAME_LEN * sizeof(char)) == 0);

    
    // Read rest of header data
    bufferUint.resize(4);
    numRead = fread(bufferUint.data(), sizeof(uint64_t), 4, file);
    assert(numRead == 4);
    c_version  = bufferUint[0];
    c_numCells = bufferUint[1];
    c_numFaces = bufferUint[2];
    c_numNodes = bufferUint[3];
    assert(c_version == SerialMesh::VERSION);

    
    // Cell Data
    bufferUint.resize(c_numCells * 9);
    numRead = fread(bufferUint.data(), sizeof(uint64_t), c_numCells * 9, file);
    assert(numRead == c_numCells * 9);

    c_cellData.resize(c_numCells);
    for (uint64_t cell = 0; cell < c_numCells; cell++) {
        c_cellData[cell].boundingFaces[0] = bufferUint[cell * 9 + 0];
        c_cellData[cell].boundingFaces[1] = bufferUint[cell * 9 + 1];
        c_cellData[cell].boundingFaces[2] = bufferUint[cell * 9 + 2];
        c_cellData[cell].boundingFaces[3] = bufferUint[cell * 9 + 3];

        c_cellData[cell].boundingNodes[0] = bufferUint[cell * 9 + 4];
        c_cellData[cell].boundingNodes[1] = bufferUint[cell * 9 + 5];
        c_cellData[cell].boundingNodes[2] = bufferUint[cell * 9 + 6];
        c_cellData[cell].boundingNodes[3] = bufferUint[cell * 9 + 7];

        c_cellData[cell].materialIndex    = bufferUint[cell * 9 + 8];
    }
    
    
    // Face Data
    bufferUint.resize(c_numFaces * 5);
    numRead = fread(bufferUint.data(), sizeof(uint64_t), c_numFaces * 5, file);
    assert(numRead == c_numFaces * 5);

    c_faceData.resize(c_numFaces);
    for (uint64_t face = 0; face < c_numFaces; face++) {
        c_faceData[face].boundingCells[0] = bufferUint[face * 5 + 0];
        c_faceData[face].boundingCells[1] = bufferUint[face * 5 + 1];

        c_faceData[face].boundingNodes[0] = bufferUint[face * 5 + 2];
        c_faceData[face].boundingNodes[1] = bufferUint[face * 5 + 3];
        c_faceData[face].boundingNodes[2] = bufferUint[face * 5 + 4];
    }
    
    
    // Node Data
    bufferDouble.resize(c_numNodes * 3);
    numRead = fread(bufferDouble.data(), sizeof(double), c_numNodes * 3, file);
    assert(numRead == c_numNodes * 3);

    c_nodeData.resize(c_numNodes);
    for (uint64_t node = 0; node < c_numNodes; node++) {
        c_nodeData[node].coords[0] = bufferDouble[node * 3 + 0];
        c_nodeData[node].coords[1] = bufferDouble[node * 3 + 1];
        c_nodeData[node].coords[2] = bufferDouble[node * 3 + 2];
    }
    
    
    // Close file
    fclose(file);
}




