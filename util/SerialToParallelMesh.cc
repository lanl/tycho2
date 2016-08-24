#include <string>
#include <cstdio>
#include <cassert>
#include <set>
#include <map>
#include <metis.h>
#include "SerialMesh.hh"
#include "ParallelMesh.hh"

using namespace std;


/*
    partitionMesh
*/
static
void partitionMesh(const int numPartitions, const SerialMesh &serialMesh, 
                   vector<uint64_t> &partitionVector)
{
    assert(partitionVector.size() == serialMesh.c_numCells);
    
    
    // Check for only 1 partition
    if (numPartitions == 1) {
        for (uint64_t cell = 0; cell < partitionVector.size(); cell++) {
            partitionVector[cell] = 0;
        }
        return;
    }
    
    
    // Parameters for metis
    idx_t nvtxs = serialMesh.c_numCells;
    idx_t ncon = 1;
    vector<idx_t> xadj(nvtxs+1);
    vector<idx_t> adjncy(2 * serialMesh.c_numFaces);  // Allocating more than needed.
    idx_t *vwgt = NULL;
    idx_t *vsize = NULL;
    idx_t *adjwgt = NULL;
    idx_t nparts = numPartitions;
    real_t *tpwgts = NULL;
    real_t *ubvec = NULL;
    idx_t *options = NULL;
    idx_t objval;
    vector<idx_t> part(nvtxs);
    
    
    // Setup adjacency graph
    xadj[0] = 0;
    for (uint64_t cell = 0; cell < serialMesh.c_numCells; cell++) {
        uint64_t numConnections = 0;
        for (uint64_t lface = 0; lface < 4; lface++) {
            uint64_t face = serialMesh.c_cellData[cell].boundingFaces[lface];
            uint64_t cell1 = serialMesh.c_faceData[face].boundingCells[0];
            uint64_t cell2 = serialMesh.c_faceData[face].boundingCells[1];
            uint64_t adjCell = (cell1 == cell) ? cell2 : cell1;
            if (adjCell != numeric_limits<uint64_t>::max()) {
                uint64_t index = xadj[cell] + numConnections;
                adjncy[index] = adjCell;
                numConnections++;
            }
        }
        xadj[cell+1] = xadj[cell] + numConnections;
    }
    
    
    // Partition graph
    int retValue = METIS_PartGraphKway(
                            &nvtxs, &ncon, xadj.data(), adjncy.data(), vwgt, 
                            vsize, adjwgt, &nparts, tpwgts, ubvec, options, 
                            &objval, part.data());
    
    assert(retValue == METIS_OK);
    
    
    // Put partition data in my format
    for (uint64_t cell = 0; cell < partitionVector.size(); cell++) {
        partitionVector[cell] = part[cell];
    }
}


/*
    createParallelMesh
*/
static
void createParallelMesh(const SerialMesh &serialMesh, 
                        const vector<uint64_t> &partitionVector,
                        const int numPartitions, 
                        ParallelMesh &parallelMesh)
{
    assert(partitionVector.size() == serialMesh.c_numCells);
    
    
    // Num partitions
    parallelMesh.c_numPartitions = numPartitions;
    parallelMesh.c_bytesOffset.resize(numPartitions + 1);
    parallelMesh.c_partitionData.resize(numPartitions);
    
    
    // Set of cells, faces, nodes for each partition
    vector<set<uint64_t>> cellsInPartition(numPartitions);
    vector<set<uint64_t>> facesInPartition(numPartitions);
    vector<set<uint64_t>> nodesInPartition(numPartitions);
    for (uint64_t cell = 0; cell < serialMesh.c_numCells; cell++) {
        int partition = partitionVector[cell];
        
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
    for (uint64_t part = 0; part < parallelMesh.c_numPartitions; part++) {
        
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
    for (uint64_t part = 0; part < parallelMesh.c_numPartitions; part++) {
        
        ParallelMesh::PartitionData &partData = parallelMesh.c_partitionData[part];
        
        
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
            ParallelMesh::CellData &cellData = partData.cellData[pcell];
            
            cellData.globalID = gcell;
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
            ParallelMesh::FaceData &faceData = partData.faceData[pface];
            
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
            ParallelMesh::NodeData &nodeData = partData.nodeData[pnode];
            
            nodeData.globalID = gnode;
            for (int i = 0; i < 3; i++) {
                nodeData.coords[i] = serialMesh.c_nodeData[gnode].coords[i];
            }
        }
    }
    
    
    // Bytes Offset
    parallelMesh.c_bytesOffset[0] = 
        32 * sizeof(char) + (3 + numPartitions) * sizeof(uint64_t);
    for (uint64_t part = 0; part < parallelMesh.c_numPartitions; part++) {
        
        // Structure sizes in bytes
        const uint64_t cellDataSize = 9 * sizeof(uint64_t);
        const uint64_t faceDataSize = 9 * sizeof(uint64_t);
        const uint64_t nodeDataSize = 3 * sizeof(double) + 1 * sizeof(uint64_t);
        
        // Size of this partition in bytes
        ParallelMesh::PartitionData &partData = parallelMesh.c_partitionData[part];
        uint64_t sizeOfPartition = 3 * sizeof(uint64_t) + 
                                 partData.numCells * cellDataSize + 
                                 partData.numFaces * faceDataSize + 
                                 partData.numNodes * nodeDataSize;
        
        // Offset in bytes
        parallelMesh.c_bytesOffset[part+1] = 
            parallelMesh.c_bytesOffset[part] + sizeOfPartition;
    }
}


/*
    main
*/
int main(int argc, char* argv[])
{
    SerialMesh serialMesh;
    ParallelMesh parallelMesh;
    string inputFile;
    string outputFile;
    vector<uint64_t> partitionVector;
    int numPartitions = 0;
    
    
    // Print utility name
    printf("--- SerialToParallelMesh Utility ---\n");
    
    
    // Get input/output files
    if (argc != 4) {
        printf("Incorrect number of arguments\n");
        printf("Usage: ./SerialToParallelMesh.x <# partitions> <inputFile> <outputFile>\n");
        printf("\n\n\n");
        return 0;
    }
    numPartitions = atoi(argv[1]);
    inputFile = argv[2];
    outputFile = argv[3];
    assert(numPartitions > 0);
    printf("Partition %s into %d partitions.\n", inputFile.c_str(), numPartitions);
    printf("Write to %s\n", outputFile.c_str());
    
    
    // Read in serial mesh and convert to parallel mesh
    serialMesh.read(inputFile);
    partitionVector.resize(serialMesh.c_numCells);
    partitionMesh(numPartitions, serialMesh, partitionVector);
    createParallelMesh(serialMesh, partitionVector, numPartitions, parallelMesh);
    parallelMesh.write(outputFile);
    
    
    // Cleanup
    printf("\n\n\n");
    return 0;
}



