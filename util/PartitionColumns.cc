#include <string>
#include <cstdio>
#include <cassert>
#include <cinttypes>
#include <set>
#include <map>
#include <algorithm>
#include "SerialMesh.hh"
#include "ParallelMesh.hh"

using namespace std;


/*
    Centroid struct
    
    Holds centroid for a given cell.
*/
struct Centroid
{
    uint64_t cell;
    double coord[3];
};


/*
    splitIntoChunks
    
    Gives indices of chunks up to size where chunks are appoximately same size.
    e.g. size = 10, numChunks = 3
         chunks are {0,1,2,3} {4,5,6} {7,8,9}
         chunkIndices = {0,4,7,10}
*/
static
vector<uint64_t> splitIntoChunks(const uint64_t size, const uint64_t numChunks)
{
    vector<uint64_t> chunkIndices(numChunks + 1);
    chunkIndices[0] = 0;
    
    for (uint64_t i = 0; i < numChunks; i++) {
        chunkIndices[i+1] = chunkIndices[i] + (size / numChunks);
        if (i < (size % numChunks))
            chunkIndices[i+1]++;
    }
    
    return chunkIndices;
}


/*
    compareX
*/
bool compareX(Centroid c1, Centroid c2) 
{
    return c1.coord[0] < c2.coord[0];
}


/*
    compareY
*/
bool compareY(Centroid c1, Centroid c2) 
{
    return c1.coord[1] < c2.coord[1];
}


/*
    partitionMesh
*/
static
void partitionMesh(const int numPartitionsX, const int numPartitionsY, 
                   const SerialMesh &serialMesh, 
                   vector<uint64_t> &partitionVector)
{
    vector<Centroid> centroids(serialMesh.c_numCells);
    vector<uint64_t> chunksX;
    vector<uint64_t> chunksY[numPartitionsX];
    
    
    // Put centroids into vector
    for (uint64_t cell = 0; cell < serialMesh.c_numCells; cell++) {
        Centroid centroid;
        centroid.cell = cell;
        centroid.coord[0] = centroid.coord[1] = centroid.coord[2] = 0.0;
        
        for (int vrtx = 0; vrtx < 4; vrtx++) {
            uint64_t node = serialMesh.c_cellData[cell].boundingNodes[vrtx];
            centroid.coord[0] += serialMesh.c_nodeData[node].coords[0];
            centroid.coord[1] += serialMesh.c_nodeData[node].coords[1];
            centroid.coord[2] += serialMesh.c_nodeData[node].coords[2];
        }
        
        centroid.coord[0] = centroid.coord[0] / 4.0;
        centroid.coord[1] = centroid.coord[1] / 4.0;
        centroid.coord[2] = centroid.coord[2] / 4.0;
        
        centroids[cell] = centroid;
    }
    
    
    // Sort across x-axis
    sort(centroids.begin(), centroids.end(), compareX);
    
    
    // Sort across y-axis in chunks
    chunksX = splitIntoChunks(centroids.size(), numPartitionsX);
    for (int i = 0; i < numPartitionsX; i++) {
        auto indexBegin = centroids.begin() + chunksX[i];
        auto indexEnd = centroids.begin() + chunksX[i+1];
        sort(indexBegin, indexEnd, compareY);
        chunksY[i] = splitIntoChunks(chunksX[i+1] - chunksX[i], numPartitionsY);
    }
    
    
    // Set partitions
    for (int i = 0; i < numPartitionsX; i++) {
    for (int j = 0; j < numPartitionsY; j++) {
        
        uint64_t part = i * numPartitionsY + j;
        uint64_t numInPartition = chunksY[i][j+1] - chunksY[i][j];
        uint64_t indexBegin = chunksX[i] + chunksY[i][j];
        
        if (part < 100) {
            printf("Cells in partition (%d,%d) = %" PRIu64 "\n", 
                   i, j, numInPartition);
        }
        
        for (uint64_t k = 0; k < numInPartition; k++) {
            uint64_t cell = centroids[indexBegin + k].cell;
            partitionVector[cell] = part;
        }
    }}
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
    int numPartitionsX = 0;
    int numPartitionsY = 0;
    
    
    // Print utility name
    printf("--- PartitionColumns Utility ---\n");
    
    
    // Get input/output files
    if (argc != 5) {
        printf("Incorrect number of arguments\n");
        printf("Usage: ./PartitionColumns.x <# partitions x> <# partitions y> "
               "<inputFile> <outputFile>\n");
        printf("\n\n\n");
        return 0;
    }
    numPartitionsX = atoi(argv[1]);
    numPartitionsY = atoi(argv[2]);
    inputFile = argv[3];
    outputFile = argv[4];
    assert(numPartitionsX > 0);
    assert(numPartitionsY > 0);
    printf("Partition %s into %d x %d partitions.\n", 
           inputFile.c_str(), numPartitionsX, numPartitionsY);
    printf("Write to %s\n", outputFile.c_str());
    
    
    // Read in serial mesh and convert to parallel mesh
    serialMesh.read(inputFile);
    partitionVector.resize(serialMesh.c_numCells);
    partitionMesh(numPartitionsX, numPartitionsY, serialMesh, partitionVector);
    createParallelMesh(serialMesh, partitionVector, 
                       numPartitionsX * numPartitionsY, parallelMesh);
    parallelMesh.write(outputFile);
    
    
    // Cleanup
    printf("\n\n\n");
    return 0;
}



