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

#include "ParallelMesh.hh"
#include <cstdio>
#include <cassert>
#include <cinttypes>
#include <cstring>
#include <set>
#include <map>

using namespace std;

static const char MESH_FORMAT_NAME[ParallelMesh::MESH_FORMAT_NAME_LEN] = 
    {'T', 'y', 'c', 'h', 'o', ' ', '2', ' ', 'P', 'a', 'r', 'a', 'l', 'l', 'e', 'l', 
     ' ', 'M', 'e', 's', 'h', 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};


/*
    print
    
    Prints the entire mesh.
*/
void ParallelMesh::print(bool printVerbose)
{
    // Header Data
    printf("Mesh Format Name: %s\n", c_meshFormatName);
    printf("Mesh Version: %" PRIu64 "\n", c_version);
    printf("Num partitions %" PRIu64 "\n", c_numPartitions);
    printf("Bytes Offset\n");
    for (uint64_t part = 0; part < c_numPartitions + 1; part++) {
        printf("   %" PRIu64 "\n", c_bytesOffset[part]);
    }
    
    
    // Partition Data
    for (uint64_t i = 0; i < c_numPartitions; i++) {
        printf("Partition %" PRIu64 "\n", i);
        printPartitionData(c_partitionData[i], printVerbose);
    }
}


/*
    printPartitionData
    
    Prints the entire mesh.
*/
void ParallelMesh::printPartitionData(const PartitionData &partData, 
                                      bool printVerbose)
{
    const char *boundaryNames[NumBoundaryTypes] = {
        "MeshBoundary", "PartitionBoundary", "NotBoundary"
    };
    
    
    printf("   Num Cells: %" PRIu64 "\n", partData.numCells);
    printf("   Num Faces: %" PRIu64 "\n", partData.numFaces);
    printf("   Num Nodes: %" PRIu64 "\n", partData.numNodes);
    
    
    // Skip the rest if not verbose print
    if (!printVerbose)
        return;

    
    // Cell Data
    for (uint64_t cell = 0; cell < partData.numCells; cell++) {
        printf("   Cell %" PRIu64 "\n", cell);
        printf("      Faces %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 "\n", 
               partData.cellData[cell].boundingFaces[0],
               partData.cellData[cell].boundingFaces[1],
               partData.cellData[cell].boundingFaces[2],
               partData.cellData[cell].boundingFaces[3]);
        printf("      Nodes %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 "\n", 
               partData.cellData[cell].boundingNodes[0],
               partData.cellData[cell].boundingNodes[1],
               partData.cellData[cell].boundingNodes[2],
               partData.cellData[cell].boundingNodes[3]);
        printf("      globalID: %" PRIu64 "\n", 
               partData.cellData[cell].globalID);
        printf("      materialIndex: %" PRIu64 "\n", 
               partData.cellData[cell].materialIndex);
    }
    
    
    // Face Data
    for (uint64_t face = 0; face < partData.numFaces; face++) {
        printf("   Face %" PRIu64 "\n", face);
        printf("      Cells %" PRIu64 " %" PRIu64 "\n", 
               partData.faceData[face].boundingCells[0],
               partData.faceData[face].boundingCells[1]);
        printf("      Nodes %" PRIu64 " %" PRIu64 " %" PRIu64 "\n", 
               partData.faceData[face].boundingNodes[0],
               partData.faceData[face].boundingNodes[1],
               partData.faceData[face].boundingNodes[2]);
        printf("      globalID: %" PRIu64 "\n", 
               partData.faceData[face].globalID);
        printf("      partition: %" PRIu64 " %" PRIu64 "\n", 
               partData.faceData[face].partition[0],
               partData.faceData[face].partition[1]);
        assert(partData.faceData[face].boundaryType < NumBoundaryTypes);
        printf("      boundaryType: %s\n", 
               boundaryNames[partData.faceData[face].boundaryType]);
    }
    
    
    // Node Data
    for (uint64_t node = 0; node < partData.numNodes; node++) {
        printf("   Node %" PRIu64 "\n", node);
        printf("      coords: (%f, %f, %f)\n", 
               partData.nodeData[node].coords[0],
               partData.nodeData[node].coords[1],
               partData.nodeData[node].coords[2]);
        printf("      globalID: %" PRIu64 "\n", 
               partData.nodeData[node].globalID);
    }
}


/*
    write
    
    Write entire mesh to file.
*/
void ParallelMesh::write(const std::string &filename)
{
    FILE *file;
    size_t numWritten;
    vector<uint64_t> bufferUint;
    
    
    // Open file
    file = fopen(filename.c_str(), "wb");
    assert(file != NULL);


    // Make sure format name and version are filled in
    memcpy(c_meshFormatName, MESH_FORMAT_NAME, 
           ParallelMesh::MESH_FORMAT_NAME_LEN * sizeof(char));
    c_version = ParallelMesh::VERSION;


    // Write mesh format name
    numWritten = fwrite(c_meshFormatName, sizeof(char),
                        ParallelMesh::MESH_FORMAT_NAME_LEN, file);
    assert(numWritten == ParallelMesh::MESH_FORMAT_NAME_LEN);


    // Buffer header data and write
    bufferUint.push_back(c_version);
    bufferUint.push_back(c_numPartitions);
    for (uint64_t part = 0; part < c_numPartitions + 1; part++) {
        bufferUint.push_back(c_bytesOffset[part]);
    }

    
    // Partition Data
    for (uint64_t part = 0; part < c_numPartitions; part++) {
        
        // Buffer header
        PartitionData &partData = c_partitionData[part];
        bufferUint.push_back(partData.numCells);
        bufferUint.push_back(partData.numFaces);
        bufferUint.push_back(partData.numNodes);
        
        
        // Buffer Cell Data
        for (uint64_t cell = 0; cell < partData.numCells; cell++) {
            bufferUint.push_back(partData.cellData[cell].boundingFaces[0]);
            bufferUint.push_back(partData.cellData[cell].boundingFaces[1]);
            bufferUint.push_back(partData.cellData[cell].boundingFaces[2]);
            bufferUint.push_back(partData.cellData[cell].boundingFaces[3]);

            bufferUint.push_back(partData.cellData[cell].boundingNodes[0]);
            bufferUint.push_back(partData.cellData[cell].boundingNodes[1]);
            bufferUint.push_back(partData.cellData[cell].boundingNodes[2]);
            bufferUint.push_back(partData.cellData[cell].boundingNodes[3]);

            bufferUint.push_back(partData.cellData[cell].globalID);
            bufferUint.push_back(partData.cellData[cell].materialIndex);
        }
        
        
        // Buffer Face Data
        for (uint64_t face = 0; face < partData.numFaces; face++) {
            bufferUint.push_back(partData.faceData[face].boundingCells[0]);
            bufferUint.push_back(partData.faceData[face].boundingCells[1]);

            bufferUint.push_back(partData.faceData[face].boundingNodes[0]);
            bufferUint.push_back(partData.faceData[face].boundingNodes[1]);
            bufferUint.push_back(partData.faceData[face].boundingNodes[2]);

            bufferUint.push_back(partData.faceData[face].globalID);

            bufferUint.push_back(partData.faceData[face].partition[0]);
            bufferUint.push_back(partData.faceData[face].partition[1]);

            bufferUint.push_back(partData.faceData[face].boundaryType);
        }
        
        
        // Buffer Node Data
        for (uint64_t node = 0; node < partData.numNodes; node++) {
            uint64_t coord0, coord1, coord2;
            
            memcpy(&coord0, &partData.nodeData[node].coords[0], sizeof(double));
            memcpy(&coord1, &partData.nodeData[node].coords[1], sizeof(double));
            memcpy(&coord2, &partData.nodeData[node].coords[2], sizeof(double));
            
            bufferUint.push_back(coord0);
            bufferUint.push_back(coord1);
            bufferUint.push_back(coord2);
            
            bufferUint.push_back(partData.nodeData[node].globalID);
        }
    }


    // Write data
    numWritten = fwrite(bufferUint.data(), sizeof(uint64_t),
                        bufferUint.size(), file);
    assert(numWritten == bufferUint.size());
    
    
    // Close file
    fclose(file);
}


/*
    read
    
    Read entire mesh from file.
*/
void ParallelMesh::read(const std::string &filename)
{
    FILE *file;
    size_t numRead;
    vector<uint64_t> bufferUint;
    
    
    // Open file
    file = fopen(filename.c_str(), "rb");
    assert(file != NULL);


    // Read mesh format name
    numRead = fread(c_meshFormatName, sizeof(char), 
                    ParallelMesh::MESH_FORMAT_NAME_LEN, file);
    assert(numRead == ParallelMesh::MESH_FORMAT_NAME_LEN);
    assert(memcmp(c_meshFormatName, MESH_FORMAT_NAME, 
                  ParallelMesh::MESH_FORMAT_NAME_LEN * sizeof(char)) == 0);
    
    
    // Header Data
    bufferUint.resize(2);
    numRead = fread(bufferUint.data(), sizeof(uint64_t), 2, file);
    assert(numRead == 2);
    c_version = bufferUint[0];
    c_numPartitions = bufferUint[1];
    assert(c_version == ParallelMesh::VERSION);


    // Partition offsets
    c_bytesOffset.resize(c_numPartitions + 1);
    numRead = fread(c_bytesOffset.data(), sizeof(uint64_t), c_numPartitions + 1, 
                    file);
    assert(numRead == c_numPartitions + 1);
    
    
    // Partition Data
    c_partitionData.resize(c_numPartitions);
    for (uint64_t i = 0; i < c_numPartitions; i++) {
        
        // Get header data
        PartitionData &partData = c_partitionData[i];
        bufferUint.resize(3);
        numRead = fread(bufferUint.data(), sizeof(uint64_t), 3, file);
        assert(numRead == 3);
        
        partData.numCells = bufferUint[0];
        partData.numFaces = bufferUint[1];
        partData.numNodes = bufferUint[2];
        
        
        // Cell Data
        bufferUint.resize(partData.numCells * 10);
        numRead = fread(bufferUint.data(), sizeof(uint64_t), 
                        partData.numCells * 10, file);
        assert(numRead == partData.numCells * 10);

        partData.cellData.resize(partData.numCells);
        for (uint64_t cell = 0; cell < partData.numCells; cell++) {
            partData.cellData[cell].boundingFaces[0] = bufferUint[cell * 10 + 0];
            partData.cellData[cell].boundingFaces[1] = bufferUint[cell * 10 + 1];
            partData.cellData[cell].boundingFaces[2] = bufferUint[cell * 10 + 2];
            partData.cellData[cell].boundingFaces[3] = bufferUint[cell * 10 + 3];

            partData.cellData[cell].boundingNodes[0] = bufferUint[cell * 10 + 4];
            partData.cellData[cell].boundingNodes[1] = bufferUint[cell * 10 + 5];
            partData.cellData[cell].boundingNodes[2] = bufferUint[cell * 10 + 6];
            partData.cellData[cell].boundingNodes[3] = bufferUint[cell * 10 + 7];
            
            partData.cellData[cell].globalID         = bufferUint[cell * 10 + 8];
            partData.cellData[cell].materialIndex    = bufferUint[cell * 10 + 9];
        }
        
        
        // Face Data
        bufferUint.resize(partData.numFaces * 9);
        numRead = fread(bufferUint.data(), sizeof(uint64_t), 
                        partData.numFaces * 9, file);
        assert(numRead == partData.numFaces * 9);

        partData.faceData.resize(partData.numFaces);
        for (uint64_t face = 0; face < partData.numFaces; face++) {
            partData.faceData[face].boundingCells[0] = bufferUint[face * 9 + 0];
            partData.faceData[face].boundingCells[1] = bufferUint[face * 9 + 1];

            partData.faceData[face].boundingNodes[0] = bufferUint[face * 9 + 2];
            partData.faceData[face].boundingNodes[1] = bufferUint[face * 9 + 3];
            partData.faceData[face].boundingNodes[2] = bufferUint[face * 9 + 4];
            
            partData.faceData[face].globalID         = bufferUint[face * 9 + 5];
            
            partData.faceData[face].partition[0]     = bufferUint[face * 9 + 6];
            partData.faceData[face].partition[1]     = bufferUint[face * 9 + 7];
            
            partData.faceData[face].boundaryType     = bufferUint[face * 9 + 8];
        }
        
        
        // Node Data
        bufferUint.resize(partData.numNodes * 4);
        numRead = fread(bufferUint.data(), sizeof(uint64_t), 
                        partData.numNodes * 4, file);
        assert(numRead == partData.numNodes * 4);

        partData.nodeData.resize(partData.numNodes);
        for (uint64_t node = 0; node < partData.numNodes; node++) {
            partData.nodeData[node].coords[0] = *((double*)(&bufferUint[node * 4 + 0]));
            partData.nodeData[node].coords[1] = *((double*)(&bufferUint[node * 4 + 1]));
            partData.nodeData[node].coords[2] = *((double*)(&bufferUint[node * 4 + 2]));
            
            partData.nodeData[node].globalID = bufferUint[node * 4 + 3];
        }
    }
    
    
    // Close file
    fclose(file);
}


/*
    createFromSerialMesh

    Creates a ParallelMesh given a SerialMesh and partition info.
*/
void ParallelMesh::createFromSerialMesh(const SerialMesh &serialMesh, 
                                        const vector<uint64_t> &partitionVector,
                                        const int numPartitions)
{
    assert(partitionVector.size() == serialMesh.c_numCells);
    
    
    // Num partitions
    c_numPartitions = numPartitions;
    c_bytesOffset.resize(numPartitions + 1);
    c_partitionData.resize(numPartitions);
    
    
    // Set of cells, faces, nodes for each partition
    vector<set<uint64_t>> cellsInPartition(numPartitions);
    vector<set<uint64_t>> facesInPartition(numPartitions);
    vector<set<uint64_t>> nodesInPartition(numPartitions);
    for (uint64_t cell = 0; cell < serialMesh.c_numCells; cell++) {
        uint64_t partition = partitionVector[cell];
        
        // cells
        cellsInPartition[partition].insert(cell);
        
        // faces
        for (int lface = 0; lface < 4; lface++) {
            int face = serialMesh.c_cellData[cell].boundingFaces[lface];
            facesInPartition[partition].insert(face);
        }
        
        // nodes
        for (int lnode = 0; lnode < 4; lnode++) {
            int node = serialMesh.c_cellData[cell].boundingNodes[lnode];
            nodesInPartition[partition].insert(node);
        }
    }
    
    
    // Global to partition indices
    vector<map<uint64_t,uint64_t>> gpCells(numPartitions);
    vector<map<uint64_t,uint64_t>> gpFaces(numPartitions);
    vector<map<uint64_t,uint64_t>> gpNodes(numPartitions);
    for (uint64_t part = 0; part < c_numPartitions; part++) {
        
        // cells
        uint64_t pcell = 0;
        for (uint64_t gcell : cellsInPartition[part]) {
            gpCells[part].insert(pair<uint64_t,uint64_t>(gcell, pcell));
            pcell++;
        }
        
        // faces
        uint64_t pface = 0;
        for (uint64_t gface : facesInPartition[part]) {
            gpFaces[part].insert(pair<uint64_t,uint64_t>(gface, pface));
            pface++;
        }
        
        // nodes
        uint64_t pnode = 0;
        for (uint64_t gnode : nodesInPartition[part]) {
            gpNodes[part].insert(pair<uint64_t,uint64_t>(gnode, pnode));
            pnode++;
        }
    }
    
    
    // Create partition data
    for (uint64_t part = 0; part < c_numPartitions; part++) {
        
        PartitionData &partData = c_partitionData[part];
        
        
        // Number of cells, faces, nodes
        partData.numCells = cellsInPartition[part].size();
        partData.numFaces = facesInPartition[part].size();
        partData.numNodes = nodesInPartition[part].size();
        
        partData.cellData.resize(partData.numCells);
        partData.faceData.resize(partData.numFaces);
        partData.nodeData.resize(partData.numNodes);
        
        
        // CellData
        for (uint64_t gcell : cellsInPartition[part]) {
            uint64_t pcell = gpCells[part].find(gcell)->second;
            CellData &cellData = partData.cellData[pcell];
            
            cellData.globalID = gcell;
            cellData.materialIndex = serialMesh.c_cellData[gcell].materialIndex;
            for (int lface = 0; lface < 4; lface++) {
                uint64_t gface = serialMesh.c_cellData[gcell].boundingFaces[lface];
                cellData.boundingFaces[lface] = gpFaces[part].find(gface)->second;
            }
            for (int lnode = 0; lnode < 4; lnode++) {
                uint64_t gnode = serialMesh.c_cellData[gcell].boundingNodes[lnode];
                cellData.boundingNodes[lnode] = gpNodes[part].find(gnode)->second;
            }
        }
        
        
        // FaceData
        for (uint64_t gface : facesInPartition[part]) {
            uint64_t pface = gpFaces[part].find(gface)->second;
            FaceData &faceData = partData.faceData[pface];
            
            faceData.globalID = gface;
            faceData.boundaryType = ParallelMesh::NotBoundary;
            // boundary type can change below
            for (int lcell = 0; lcell < 2; lcell++) {
                uint64_t gcell = serialMesh.c_faceData[gface].boundingCells[lcell];
                if (gcell == numeric_limits<uint64_t>::max()) {
                    faceData.boundingCells[lcell] = numeric_limits<uint64_t>::max();
                    faceData.boundaryType = ParallelMesh::MeshBoundary;
                    faceData.partition[lcell] = numeric_limits<uint64_t>::max();
                }
                else if (partitionVector[gcell] != part) {
                    int part2 = partitionVector[gcell];
                    faceData.boundingCells[lcell] = gpCells[part2].find(gcell)->second;
                    faceData.boundaryType = ParallelMesh::PartitionBoundary;
                    faceData.partition[lcell] = part2;
                }
                else {
                    faceData.boundingCells[lcell] = gpCells[part].find(gcell)->second;
                    faceData.partition[lcell] = part;
                }
            }
            for (int lnode = 0; lnode < 3; lnode++) {
                uint64_t gnode = serialMesh.c_faceData[gface].boundingNodes[lnode];
                faceData.boundingNodes[lnode] = gpNodes[part].find(gnode)->second;
            }
        }
        
        
        // NodeData
        for (uint64_t gnode : nodesInPartition[part]) {
            uint64_t pnode = gpNodes[part].find(gnode)->second;
            NodeData &nodeData = partData.nodeData[pnode];
            
            nodeData.globalID = gnode;
            for (int i = 0; i < 3; i++) {
                nodeData.coords[i] = serialMesh.c_nodeData[gnode].coords[i];
            }
        }
    }
    
    
    // Bytes Offset
    c_bytesOffset[0] = 32 * sizeof(char) + (3 + numPartitions) * sizeof(uint64_t);
    for (uint64_t part = 0; part < c_numPartitions; part++) {
        
        // Structure sizes in bytes
        const uint64_t cellDataSize = 10 * sizeof(uint64_t);
        const uint64_t faceDataSize = 9 * sizeof(uint64_t);
        const uint64_t nodeDataSize = 3 * sizeof(double) + 1 * sizeof(uint64_t);
        
        // Size of this partition in bytes
        PartitionData &partData = c_partitionData[part];
        uint64_t sizeOfPartition = 3 * sizeof(uint64_t) + 
                                 partData.numCells * cellDataSize + 
                                 partData.numFaces * faceDataSize + 
                                 partData.numNodes * nodeDataSize;
        
        // Offset in bytes
        c_bytesOffset[part+1] = c_bytesOffset[part] + sizeOfPartition;
    }
}





/*
    Allow only serial read/write by defining PARALLEL_MESH_READ_SERIAL_ONLY.
*/
#ifndef PARALLEL_MESH_READ_SERIAL_ONLY
#include <mpi.h>
#include <Comm.hh>


/*
    readInParallel
    
    Read only 1 partition of mesh per MPI process from file.
*/
void ParallelMesh::readInParallel(const std::string &filename, 
                                  PartitionData &partData)
{
    MPI_File file;
    uint64_t offset = 0;
    vector<uint64_t> bufferUint;
    
    
    // Open file
    Comm::openFileForRead(filename, file);


    // Check header
    if (Comm::rank() == 0) {
        char meshFormatName[32];
        uint64_t version;
        uint64_t numPartitions;
        Comm::seek(file, 0);
        Comm::readChars(file, meshFormatName, 32);
        assert(memcmp(meshFormatName, MESH_FORMAT_NAME, 
                      ParallelMesh::MESH_FORMAT_NAME_LEN * sizeof(char)) == 0);
        
        bufferUint.resize(2);
        Comm::readUint64(file, bufferUint.data(), 2);
        version = bufferUint[0];
        numPartitions = bufferUint[1];
        assert(version == ParallelMesh::VERSION);
        assert(numPartitions == (uint64_t)Comm::numRanks());
    }
    Comm::barrier();
    
    
    // Read in bytes offset
    Comm::seek(file, 32 * sizeof(char) + (2 + Comm::rank()) * sizeof(uint64_t));
    Comm::readUint64(file, offset);
    Comm::seek(file, offset);
    
    
    // Num Cells, Face, Nodes
    bufferUint.resize(3);
    Comm::readUint64(file, bufferUint.data(), 3);
    partData.numCells = bufferUint[0];
    partData.numFaces = bufferUint[1];
    partData.numNodes = bufferUint[2];
    
    
    // Cell Data
    bufferUint.resize(partData.numCells * 10);
    Comm::readUint64(file, bufferUint.data(), partData.numCells * 10);
    partData.cellData.resize(partData.numCells);
    
    for (uint64_t cell = 0; cell < partData.numCells; cell++) {
        partData.cellData[cell].boundingFaces[0] = bufferUint[cell * 10 + 0];
        partData.cellData[cell].boundingFaces[1] = bufferUint[cell * 10 + 1];
        partData.cellData[cell].boundingFaces[2] = bufferUint[cell * 10 + 2];
        partData.cellData[cell].boundingFaces[3] = bufferUint[cell * 10 + 3];

        partData.cellData[cell].boundingNodes[0] = bufferUint[cell * 10 + 4];
        partData.cellData[cell].boundingNodes[1] = bufferUint[cell * 10 + 5];
        partData.cellData[cell].boundingNodes[2] = bufferUint[cell * 10 + 6];
        partData.cellData[cell].boundingNodes[3] = bufferUint[cell * 10 + 7];

        partData.cellData[cell].globalID         = bufferUint[cell * 10 + 8];
        partData.cellData[cell].materialIndex    = bufferUint[cell * 10 + 9];
    }
    
    
    // Face Data
    bufferUint.resize(partData.numFaces * 9);
    Comm::readUint64(file, bufferUint.data(), partData.numFaces * 9);
    partData.faceData.resize(partData.numFaces);
    
    for (uint64_t face = 0; face < partData.numFaces; face++) {
        partData.faceData[face].boundingCells[0] = bufferUint[face * 9 + 0];
        partData.faceData[face].boundingCells[1] = bufferUint[face * 9 + 1];

        partData.faceData[face].boundingNodes[0] = bufferUint[face * 9 + 2];
        partData.faceData[face].boundingNodes[1] = bufferUint[face * 9 + 3];
        partData.faceData[face].boundingNodes[2] = bufferUint[face * 9 + 4];

        partData.faceData[face].globalID         = bufferUint[face * 9 + 5];

        partData.faceData[face].partition[0]     = bufferUint[face * 9 + 6];
        partData.faceData[face].partition[1]     = bufferUint[face * 9 + 7];

        partData.faceData[face].boundaryType     = bufferUint[face * 9 + 8];
        assert(partData.faceData[face].boundaryType < NumBoundaryTypes);
    }
    
    
    // Node Data
    bufferUint.resize(partData.numNodes * 4);
    Comm::readUint64(file, bufferUint.data(), partData.numNodes * 4);
    partData.nodeData.resize(partData.numNodes);
    
    for (uint64_t node = 0; node < partData.numNodes; node++) {
        partData.nodeData[node].coords[0] = *((double*)(&bufferUint[node * 4 + 0]));
        partData.nodeData[node].coords[1] = *((double*)(&bufferUint[node * 4 + 1]));
        partData.nodeData[node].coords[2] = *((double*)(&bufferUint[node * 4 + 2]));

        partData.nodeData[node].globalID  = bufferUint[node * 4 + 3];
    }
    
    
    // Close file
    Comm::closeFile(file);
}

#endif

