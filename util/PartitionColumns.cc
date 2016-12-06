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
#include <cinttypes>
#include <set>
#include <map>
#include <algorithm>
#include <stdint.h>
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
    parallelMesh.createFromSerialMesh(serialMesh, partitionVector, 
                                      numPartitionsX * numPartitionsY);
    parallelMesh.write(outputFile);
    
    
    // Cleanup
    printf("\n\n\n");
    return 0;
}



