#include <string>
#include <cstdio>
#include <cassert>
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
    parallelMesh.createFromSerialMesh(serialMesh, partitionVector, numPartitions);
    parallelMesh.write(outputFile);
    
    
    // Cleanup
    printf("\n\n\n");
    return 0;
}



