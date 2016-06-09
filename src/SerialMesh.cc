/*
    SerialMesh.cc
    
    Implements a tetrahedral mesh.  Look at SerialMesh.hh for file format.
*/

#include "SerialMesh.hh"
#include <cstdio>
#include <cassert>
#include <cinttypes>

using namespace std;


/*
    printAll
    
    Prints the entire mesh data.
*/
void SerialMesh::printAll()
{
    // Header Data
    printf("Num Cells: %" PRIu64 "\n", c_numCells);
    printf("Num Faces: %" PRIu64 "\n", c_numFaces);
    printf("Num Nodes: %" PRIu64 "\n", c_numNodes);
    
    
    // Cell Data
    for (uint64_t cell = 0; cell < c_numCells; cell++) {
        printf("Cell %" PRIu64 "\n", cell);
        printf("   Faces %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 "\n", 
               c_cellData[cell].boundingFaces[0],
               c_cellData[cell].boundingFaces[1],
               c_cellData[cell].boundingFaces[2],
               c_cellData[cell].boundingFaces[3]);
        printf("   Nodes %" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 "\n", 
               c_cellData[cell].boundingNodes[0],
               c_cellData[cell].boundingNodes[1],
               c_cellData[cell].boundingNodes[2],
               c_cellData[cell].boundingNodes[3]);
    }
    
    
    // Face Data
    for (uint64_t face = 0; face < c_numFaces; face++) {
        printf("Face %" PRIu64 "\n", face);
        printf("   Cells %" PRIu64 " %" PRIu64 "\n", c_faceData[face].boundingCells[0],
                                       c_faceData[face].boundingCells[1]);
        printf("   Nodes %" PRIu64 " %" PRIu64 " %" PRIu64 "\n", c_faceData[face].boundingNodes[0],
                                            c_faceData[face].boundingNodes[1],
                                            c_faceData[face].boundingNodes[2]);
    }
    
    
    // Node Data
    for (uint64_t node = 0; node < c_numNodes; node++) {
        printf("Node %" PRIu64 ": (%f, %f, %f)\n", node, c_nodeData[node].coords[0],
                                                  c_nodeData[node].coords[1],
                                                  c_nodeData[node].coords[2]);
    }
}


/*
    printSummary
    
    Prints number of cells, faces, and nodes.
    Also prints number of boundary faces vs interior faces.
*/
void SerialMesh::printSummary()
{
    uint64_t bdryFace = 0;
    uint64_t intFace = 0;
    for (uint64_t face = 0; face < c_numFaces; face++) {
        if (c_faceData[face].boundingCells[1] == INVALID_INDEX)
            bdryFace++;
        else
            intFace++;
    }
    
    printf("Cells: %" PRIu64 ", Nodes: %" PRIu64 ", Faces: %" PRIu64 "\n", 
           c_numCells, c_numNodes, c_numFaces);
    printf("Boundary Faces: %" PRIu64 "   Interior Faces: %" PRIu64 "\n", bdryFace, intFace);
}


/*
    write
    
    Write serial mesh to file.
*/
void SerialMesh::write(const std::string &filename)
{
    FILE *file;
    size_t numWritten;
    
    
    // Open file
    file = fopen(filename.c_str(), "wb");
    assert(file != NULL);
    
    
    // Header Data
    numWritten = fwrite(&c_numCells, sizeof(uint64_t), 1, file);
    assert(numWritten == 1);
    numWritten = fwrite(&c_numFaces, sizeof(uint64_t), 1, file);
    assert(numWritten == 1);
    numWritten = fwrite(&c_numNodes, sizeof(uint64_t), 1, file);
    assert(numWritten == 1);
    
    
    // Cell Data
    for (uint64_t cell = 0; cell < c_numCells; cell++) {
        numWritten = fwrite(c_cellData[cell].boundingFaces, sizeof(uint64_t), 4, file);
        assert(numWritten == 4);
        numWritten = fwrite(c_cellData[cell].boundingNodes, sizeof(uint64_t), 4, file);
        assert(numWritten == 4);
    }
    
    
    // Face Data
    for (uint64_t face = 0; face < c_numFaces; face++) {
        numWritten = fwrite(c_faceData[face].boundingCells, sizeof(uint64_t), 2, file);
        assert(numWritten == 2);
        numWritten = fwrite(c_faceData[face].boundingNodes, sizeof(uint64_t), 3, file);
        assert(numWritten == 3);
    }
    
    
    // Node Data
    for (uint64_t node = 0; node < c_numNodes; node++) {
        numWritten = fwrite(c_nodeData[node].coords, sizeof(double), 3, file);
        assert(numWritten == 3);
    }
    
    
    // Close file
    fclose(file);
}


/*
    read
    
    Read serial mesh from file.
*/
void SerialMesh::read(const std::string &filename)
{
    FILE *file;
    size_t numRead;
    
    
    // Open file
    file = fopen(filename.c_str(), "rb");
    assert(file != NULL);
    
    
    // Header Data
    numRead = fread(&c_numCells, sizeof(uint64_t), 1, file);
    assert(numRead == 1);
    numRead = fread(&c_numFaces, sizeof(uint64_t), 1, file);
    assert(numRead == 1);
    numRead = fread(&c_numNodes, sizeof(uint64_t), 1, file);
    assert(numRead == 1);
    
    
    // Cell Data
    c_cellData.resize(c_numCells);
    for (uint64_t cell = 0; cell < c_numCells; cell++) {
        numRead = fread(c_cellData[cell].boundingFaces, sizeof(uint64_t), 4, file);
        assert(numRead == 4);
        numRead = fread(c_cellData[cell].boundingNodes, sizeof(uint64_t), 4, file);
        assert(numRead == 4);
    }
    
    
    // Face Data
    c_faceData.resize(c_numFaces);
    for (uint64_t face = 0; face < c_numFaces; face++) {
        numRead = fread(c_faceData[face].boundingCells, sizeof(uint64_t), 2, file);
        assert(numRead == 2);
        numRead = fread(c_faceData[face].boundingNodes, sizeof(uint64_t), 3, file);
        assert(numRead == 3);
    }
    
    
    // Node Data
    c_nodeData.resize(c_numNodes);
    for (uint64_t node = 0; node < c_numNodes; node++) {
        numRead = fread(c_nodeData[node].coords, sizeof(double), 3, file);
        assert(numRead == 3);
    }
    
    
    // Close file
    fclose(file);
}




