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
    printf("--- PartitionMetis Utility ---\n");
    
    
    // Get input/output files
    if (argc != 4) {
        printf("Incorrect number of arguments\n");
        printf("Usage: ./PartitionMetis.x <# partitions> <inputFile> <outputFile>\n");
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



