/*
    ParallelMesh.cc
    
    Implements a tetrahedral mesh that has been partitioned.
    Look at ParallelMesh.hh for the file format.
    This class can read in, store, and write the entire partitioned mesh.
    Alternatively, it can read in just one partition per MPI process in parallel.
*/

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

using namespace std;


/*
    printAll
    
    Prints the entire mesh.
*/
void ParallelMesh::printAll()
{
    const char *boundaryNames[NumBoundaryTypes] = {
        "MeshBoundary", "PartitionBoundary", "NotBoundary"
    };
    
    
    // Header Data
    printf("Num partitions %" PRIu64 "\n", c_numPartitions);
    printf("Bytes Offset\n");
    for (uint64_t part = 0; part < c_numPartitions + 1; part++) {
        printf("   %" PRIu64 "\n", c_bytesOffset[part]);
    }
    
    
    // Partition Data
    for (uint64_t i = 0; i < c_numPartitions; i++) {
        
        PartitionData &partData = c_partitionData[i];
        printf("Partition %" PRIu64 "\n", i);
        printf("   Num Cells: %" PRIu64 "\n", partData.numCells);
        printf("   Num Faces: %" PRIu64 "\n", partData.numFaces);
        printf("   Num Nodes: %" PRIu64 "\n", partData.numNodes);
        
        
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
}


/*
    printSummary
    
    Prints number of cells, faces, and nodes for each partition.
    Also prints number of faces for each boundary type.
*/
void ParallelMesh::printSummary()
{
    for (uint64_t part = 0; part < c_numPartitions; part++) {
        
        PartitionData &partData = c_partitionData[part];
        uint64_t meshFace = 0;
        uint64_t partFace = 0;
        uint64_t intFace = 0;
        
        for (uint64_t face = 0; face < partData.numFaces; face++) {
            uint64_t boundaryType = partData.faceData[face].boundaryType;
            if (boundaryType == MeshBoundary)
                meshFace++;
            else if (boundaryType == PartitionBoundary)
                partFace++;
            else if (boundaryType == NotBoundary)
                intFace++;
            else
                assert(false);
        }
        
        printf("Partition: %" PRIu64 "\n", part);
        printf("   Cells: %" PRIu64 ", Nodes: %" PRIu64 ", Faces: %" PRIu64 "\n", 
               partData.numCells, partData.numNodes, partData.numFaces);
        printf("   Mesh Boundary Faces:      %" PRIu64 "\n", meshFace);
        printf("   Partition Boundary Faces: %" PRIu64 "\n", partFace);
        printf("   Interior Faces:           %" PRIu64 "\n", intFace);
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
    
    
    // Open file
    file = fopen(filename.c_str(), "wb");
    assert(file != NULL);
    
    
    // Header Data
    numWritten = fwrite(&c_numPartitions, sizeof(uint64_t), 1, file);
    assert(numWritten == 1);
    numWritten = fwrite(&c_bytesOffset[0], sizeof(uint64_t), 
                        c_numPartitions + 1, file);
    assert(numWritten == c_numPartitions + 1);
    
    
    // Partition Data
    for (uint64_t part = 0; part < c_numPartitions; part++) {
        
        PartitionData &partData = c_partitionData[part];
        numWritten = fwrite(&partData.numCells, sizeof(uint64_t), 1, file);
        assert(numWritten == 1);
        numWritten = fwrite(&partData.numFaces, sizeof(uint64_t), 1, file);
        assert(numWritten == 1);
        numWritten = fwrite(&partData.numNodes, sizeof(uint64_t), 1, file);
        assert(numWritten == 1);
        
        
        // Cell Data
        for (uint64_t cell = 0; cell < partData.numCells; cell++) {
            numWritten = fwrite(partData.cellData[cell].boundingFaces, 
                                sizeof(uint64_t), 4, file);
            assert(numWritten == 4);
            numWritten = fwrite(partData.cellData[cell].boundingNodes, 
                                sizeof(uint64_t), 4, file);
            assert(numWritten == 4);
            numWritten = fwrite(&partData.cellData[cell].globalID, 
                                sizeof(uint64_t), 1, file);
            assert(numWritten == 1);
        }
        
        
        // Face Data
        for (uint64_t face = 0; face < partData.numFaces; face++) {
            numWritten = fwrite(partData.faceData[face].boundingCells, 
                                sizeof(uint64_t), 2, file);
            assert(numWritten == 2);
            numWritten = fwrite(partData.faceData[face].boundingNodes, 
                                sizeof(uint64_t), 3, file);
            assert(numWritten == 3);
            numWritten = fwrite(&partData.faceData[face].globalID, 
                                sizeof(uint64_t), 1, file);
            assert(numWritten == 1);
            numWritten = fwrite(partData.faceData[face].partition, 
                                sizeof(uint64_t), 2, file);
            assert(numWritten == 2);
            numWritten = fwrite(&partData.faceData[face].boundaryType, 
                                sizeof(uint64_t), 1, file);
            assert(numWritten == 1);
            assert(partData.faceData[face].boundaryType < NumBoundaryTypes);
        }
        
        
        // Node Data
        for (uint64_t node = 0; node < partData.numNodes; node++) {
            numWritten = fwrite(partData.nodeData[node].coords, 
                                sizeof(double), 3, file);
            assert(numWritten == 3);
            numWritten = fwrite(&partData.nodeData[node].globalID, 
                                sizeof(uint64_t), 1, file);
            assert(numWritten == 1);
        }
    }
    
    
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
    
    
    // Open file
    file = fopen(filename.c_str(), "rb");
    assert(file != NULL);
    
    
    // Header Data
    numRead = fread(&c_numPartitions, sizeof(uint64_t), 1, file);
    assert(numRead == 1);
    c_bytesOffset.resize(c_numPartitions + 1);
    numRead = fread(&c_bytesOffset[0], sizeof(uint64_t), 
                    c_numPartitions + 1, file);
    assert(numRead == c_numPartitions + 1);
    
    
    // Partition Data
    c_partitionData.resize(c_numPartitions);
    for (uint64_t i = 0; i < c_numPartitions; i++) {
        
        PartitionData &partData = c_partitionData[i];
        numRead = fread(&partData.numCells, sizeof(uint64_t), 1, file);
        assert(numRead == 1);
        numRead = fread(&partData.numFaces, sizeof(uint64_t), 1, file);
        assert(numRead == 1);
        numRead = fread(&partData.numNodes, sizeof(uint64_t), 1, file);
        assert(numRead == 1);
        
        
        // Cell Data
        partData.cellData.resize(partData.numCells);
        for (uint64_t cell = 0; cell < partData.numCells; cell++) {
            numRead = fread(partData.cellData[cell].boundingFaces, 
                            sizeof(uint64_t), 4, file);
            assert(numRead == 4);
            numRead = fread(partData.cellData[cell].boundingNodes, 
                            sizeof(uint64_t), 4, file);
            assert(numRead == 4);
            numRead = fread(&partData.cellData[cell].globalID, 
                            sizeof(uint64_t), 1, file);
            assert(numRead == 1);
        }
        
        
        // Face Data
        partData.faceData.resize(partData.numFaces);
        for (uint64_t face = 0; face < partData.numFaces; face++) {
            numRead = fread(partData.faceData[face].boundingCells, 
                            sizeof(uint64_t), 2, file);
            assert(numRead == 2);
            numRead = fread(partData.faceData[face].boundingNodes, 
                            sizeof(uint64_t), 3, file);
            assert(numRead == 3);
            numRead = fread(&partData.faceData[face].globalID, 
                            sizeof(uint64_t), 1, file);
            assert(numRead == 1);
            numRead = fread(partData.faceData[face].partition, 
                            sizeof(uint64_t), 2, file);
            assert(numRead == 2);
            numRead = fread(&partData.faceData[face].boundaryType, 
                            sizeof(uint64_t), 1, file);
            assert(numRead == 1);
            assert(partData.faceData[face].boundaryType < NumBoundaryTypes);
        }
        
        
        // Node Data
        partData.nodeData.resize(partData.numNodes);
        for (uint64_t node = 0; node < partData.numNodes; node++) {
            numRead = fread(&partData.nodeData[node].coords, 
                            sizeof(double), 3, file);
            assert(numRead == 3);
            numRead = fread(&partData.nodeData[node].globalID, 
                            sizeof(uint64_t), 1, file);
            assert(numRead == 1);
        }
    }
    
    
    // Close file
    fclose(file);
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
    vector<uint64_t> buffer;
    
    
    // Open file
    Comm::openFileForRead(filename, file);
    
    
    // Read num partitions (This is just a sanity check)
    if (Comm::rank() == 0) {
        uint64_t numPartitions;
        Comm::seek(file, 0);
        Comm::readUint64(file, numPartitions);
        assert(numPartitions == (uint64_t)Comm::numRanks());
    }
    Comm::barrier();
    
    
    // Read in bytes offset
    Comm::seek(file, (1 + Comm::rank()) * sizeof(uint64_t));
    Comm::readUint64(file, offset);
    Comm::seek(file, offset);
    
    
    // Num Cells, Face, Nodes
    Comm::readUint64(file, partData.numCells);
    Comm::readUint64(file, partData.numFaces);
    Comm::readUint64(file, partData.numNodes);
    
    
    // Cell Data
    buffer.resize(partData.numCells * 9);
    Comm::readUint64(file, buffer.data(), partData.numCells * 9);
    partData.cellData.resize(partData.numCells);
    
    for (uint64_t cell = 0; cell < partData.numCells; cell++) {
        partData.cellData[cell].boundingFaces[0] = buffer[cell * 9 + 0];
        partData.cellData[cell].boundingFaces[1] = buffer[cell * 9 + 1];
        partData.cellData[cell].boundingFaces[2] = buffer[cell * 9 + 2];
        partData.cellData[cell].boundingFaces[3] = buffer[cell * 9 + 3];
        partData.cellData[cell].boundingNodes[0] = buffer[cell * 9 + 4];
        partData.cellData[cell].boundingNodes[1] = buffer[cell * 9 + 5];
        partData.cellData[cell].boundingNodes[2] = buffer[cell * 9 + 6];
        partData.cellData[cell].boundingNodes[3] = buffer[cell * 9 + 7];
        partData.cellData[cell].globalID         = buffer[cell * 9 + 8];
    }
    
    
    // Face Data
    buffer.resize(partData.numFaces * 9);
    Comm::readUint64(file, buffer.data(), partData.numFaces * 9);
    partData.faceData.resize(partData.numFaces);
    
    for (uint64_t face = 0; face < partData.numFaces; face++) {
        partData.faceData[face].boundingCells[0] = buffer[face * 9 + 0];
        partData.faceData[face].boundingCells[1] = buffer[face * 9 + 1];
        partData.faceData[face].boundingNodes[0] = buffer[face * 9 + 2];
        partData.faceData[face].boundingNodes[1] = buffer[face * 9 + 3];
        partData.faceData[face].boundingNodes[2] = buffer[face * 9 + 4];
        partData.faceData[face].globalID         = buffer[face * 9 + 5];
        partData.faceData[face].partition[0]     = buffer[face * 9 + 6];
        partData.faceData[face].partition[1]     = buffer[face * 9 + 7];
        partData.faceData[face].boundaryType     = buffer[face * 9 + 8];
        assert(partData.faceData[face].boundaryType < NumBoundaryTypes);
    }
    
    
    // Node Data
    buffer.resize(partData.numNodes * 4);
    Comm::readUint64(file, buffer.data(), partData.numNodes * 4);
    partData.nodeData.resize(partData.numNodes);
    
    for (uint64_t node = 0; node < partData.numNodes; node++) {
        partData.nodeData[node].coords[0] = *((double*)(&buffer[node * 4 + 0]));
        partData.nodeData[node].coords[1] = *((double*)(&buffer[node * 4 + 1]));
        partData.nodeData[node].coords[2] = *((double*)(&buffer[node * 4 + 2]));
        partData.nodeData[node].globalID  = buffer[node * 4 + 3];
    }
    
    
    // Close file
    Comm::closeFile(file);
}

#endif

