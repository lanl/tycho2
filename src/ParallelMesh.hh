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

#ifndef __PARALLEL_MESH_HH__
#define __PARALLEL_MESH_HH__

/*
    File format on disk.
    
    (Array<char>(32)) "Tycho 2 Parallel Mesh" (Rest of characters set to 0)
    (uint64_t) version of parallel mesh format (version 2 currently)
    (uint64_t) number of partitions
    (Array<uint64_t>(num partitions + 1)) offsets in bytes for the partition data
    (Array<PartitionData>(num partitions)) partition data
    
    Partition Data Format
    (uint64_t) number of cells
    (uint64_t) number of faces
    (uint64_t) number of nodes
    (Array<uint64_t>(num cells)) cell data
    (Array<uint64_t>(num faces)) face data
    (Array<uint64_t>(num nodes)) node data
    
    Cell Data Format:
    (uint64_t[4]) bounding faces
    (uint64_t[4]) bounding nodes
    (uint64_t)    global ID
    (uint64_t)    material index
    
    Face Data Format:
    (uint64_t[2]) bounding cells
    (uint64_t[3]) bounding nodes
    (uint64_t) global ID
    (uint64_t[2]) partition
    (uint64_t) boundary type
    (Note: if only one bounding cell because of a mesh boundary,  
    bounding cell 2 has index INVALID_INDEX)
    (Note: if bounding cell is in another partition, the local cell id for that
    partition is inserted)
    (Note: if on boundary of mesh, partition 2 = INVALID_INDEX)
    
    Node Data Format:
    (double[3]) coordinates
    (uint64_t) global ID
*/

#include "SerialMesh.hh"
#include <string>
#include <vector>
#include <limits>


class ParallelMesh
{
public:
    
    static const uint64_t INVALID_INDEX = UINT64_MAX; // Indicates boundary
    static const uint64_t VERSION = 2;
    static const uint64_t MESH_FORMAT_NAME_LEN = 32;
    
    
    // Sub-structures
    enum BoundaryType
    {
        MeshBoundary, PartitionBoundary, NotBoundary, NumBoundaryTypes
    };
    struct CellData
    {
        uint64_t boundingFaces[4];
        uint64_t boundingNodes[4];
        uint64_t globalID;
        uint64_t materialIndex;
    };
    struct FaceData
    {
        uint64_t boundingCells[2];
        uint64_t boundingNodes[3];
        uint64_t globalID;
        uint64_t partition[2];
        uint64_t boundaryType;
    };
    struct NodeData
    {
        double coords[3];
        uint64_t globalID;
    };
    struct PartitionData
    {
        uint64_t numCells;
        uint64_t numFaces;
        uint64_t numNodes;
        std::vector<CellData> cellData;
        std::vector<FaceData> faceData;
        std::vector<NodeData> nodeData;
    };
    
    
    // Data
    char c_meshFormatName[MESH_FORMAT_NAME_LEN];
    uint64_t c_version;
    uint64_t c_numPartitions;
    std::vector<uint64_t> c_bytesOffset;
    std::vector<PartitionData> c_partitionData;
    
    
    // Read and write functions
    void write(const std::string &filename);
    void read(const std::string &filename);
    void print(bool printVerbose);
    void createFromSerialMesh(const SerialMesh &serialMesh, 
                              const std::vector<uint64_t> &partitionVector,
                              const int numPartitions);
    
    static
    void printPartitionData(const PartitionData &partData, bool printVerbose);

    #ifndef PARALLEL_MESH_READ_SERIAL_ONLY
    static
    void readInParallel(const std::string &filename, PartitionData &partData);
    #endif
};


#endif
