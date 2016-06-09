/*
    SerialMesh.hh
*/

#ifndef __SERIAL_MESH_HH__
#define __SERIAL_MESH_HH__

/*
    File format on disk.
    
    (uint64_t) number of cells
    (uint64_t) number of faces
    (uint64_t) number of nodes
    (Array<uint64_t>(num cells)) cell data
    (Array<uint64_t>(num cells)) face data
    (Array<uint64_t>(num cells)) node data
    
    Cell Data Format:
    (uint64_t[4]) bounding faces
    (uint64_t[4]) bounding nodes
    
    Face Data Format:
    (uint64_t[2]) bounding cells
    (uint64_t[3]) bounding nodes
    (Note: if only one bounding cell because of a boundary, the second bounding 
    cell has value INVALID_INDEX)
    
    Node Data Format:
    (double[3]) coordinates
*/

#include <string>
#include <vector>
#include <limits>


class SerialMesh
{
public:
    
    // Sub-structures
    struct CellData
    {
        uint64_t boundingFaces[4];
        uint64_t boundingNodes[4];
    };
    struct FaceData
    {
        uint64_t boundingCells[2];
        uint64_t boundingNodes[3];
    };
    struct NodeData
    {
        double coords[3];
    };
    
    
    // Data
    uint64_t c_numCells;
    uint64_t c_numFaces;
    uint64_t c_numNodes;
    std::vector<CellData> c_cellData;
    std::vector<FaceData> c_faceData;
    std::vector<NodeData> c_nodeData;
    
    
    // Read and write functions
    void write(const std::string &filename);
    void read(const std::string &filename);
    void printAll();
    void printSummary();
    
    
    // Used mainly to indicate a boundary
    static const uint64_t INVALID_INDEX = UINT64_MAX;
};


#endif
