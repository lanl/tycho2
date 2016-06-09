/*
    ParallelMesh.hh
*/

#ifndef __PARALLEL_MESH_HH__
#define __PARALLEL_MESH_HH__

/*
    File format on disk.
    
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
    (uint64_t) global ID
    
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

#include <string>
#include <vector>
#include <limits>


class ParallelMesh
{
public:
    
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
    uint64_t c_numPartitions;
    std::vector<uint64_t> c_bytesOffset;
    std::vector<PartitionData> c_partitionData;
    
    
    // Read and write functions
    void write(const std::string &filename);
    void read(const std::string &filename);
    void printAll();
    void printSummary();
    
    #ifndef PARALLEL_MESH_READ_SERIAL_ONLY
    static
    void readInParallel(const std::string &filename, PartitionData &partData);
    #endif
    
    
    // Used mainly to indicate a boundary
    static const uint64_t INVALID_INDEX = UINT64_MAX;
};


#endif
