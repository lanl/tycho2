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

/*
    RefineSerialMesh.cc
    
    Refines a serial mesh by placing a node at the centroid of every cell and
    drawing edges from this new node to each of the 4 enclosing nodes of the 
    cell.  This makes 4 times the number of tetrahedral cells.
*/

#include <string>
#include <cstdio>
#include <cassert>
#include <cfloat>
#include <cmath>
#include <cinttypes>
#include <map>
#include "../src/SerialMesh.hh"

using namespace std;


/*
    NodeCoords structure
*/
struct NodeCoords
{
    double coords[3];
    static double graphDiameter;
    
    // Implement n1 < n2 needed for map structure.
    bool operator()(const NodeCoords &n1, const NodeCoords &n2) const
    {
        double TOL = 1e-14 * graphDiameter;
        if (n1.coords[0] < n2.coords[0] - TOL)
            return true;
        if (n1.coords[0] > n2.coords[0] + TOL)
            return false;
        if (n1.coords[1] < n2.coords[1] - TOL)
            return true;
        if (n1.coords[1] > n2.coords[1] + TOL)
            return false;
        if (n1.coords[2] < n2.coords[2] - TOL)
            return true;
        if (n1.coords[2] > n2.coords[2] + TOL)
            return false;
        return false;
    }
};
double NodeCoords::graphDiameter;


/*
    FaceNodes structure
*/
struct FaceNodes
{
    uint64_t nodes[3];
    
    // Implement fUnordered1 < fUnordered2 needed for map structure.
    bool operator()(const FaceNodes &fUnordered1, const FaceNodes &fUnordered2) const
    {
        // Order nodes on face
        FaceNodes f1 = orderFaceNodes(fUnordered1);
        FaceNodes f2 = orderFaceNodes(fUnordered2);
        
        if (f1.nodes[0] < f2.nodes[0])
            return true;
        if (f1.nodes[0] > f2.nodes[0])
            return false;
        if (f1.nodes[1] < f2.nodes[1])
            return true;
        if (f1.nodes[1] > f2.nodes[1])
            return false;
        if (f1.nodes[2] < f2.nodes[2])
            return true;
        if (f1.nodes[2] > f2.nodes[2])
            return false;
        return false;
    }
    
    // Order nodes
    static
    FaceNodes orderFaceNodes(const FaceNodes &fIn)
    {
        FaceNodes fOut;
    
        fOut.nodes[0] = min(fIn.nodes[0], min(fIn.nodes[1], fIn.nodes[2]));
        fOut.nodes[2] = max(fIn.nodes[0], max(fIn.nodes[1], fIn.nodes[2]));
    
        if (fIn.nodes[0] != fOut.nodes[0] && fIn.nodes[0] != fOut.nodes[2])
            fOut.nodes[1] = fIn.nodes[0];
        else if (fIn.nodes[1] != fOut.nodes[0] && fIn.nodes[1] != fOut.nodes[2])
            fOut.nodes[1] = fIn.nodes[1];
        else
            fOut.nodes[1] = fIn.nodes[2];
    
        return fOut;
    }
};


/*
    Tet structure
*/
struct Tet
{
    NodeCoords nodes[4];
};


/*
    index2To1
    
    Edge midpoints are indexed via a tuple.  Turn the tuple into a single index.
*/
static
int index2To1(int i, int j)
{
    if (i == 0 && j == 1) return 0;
    if (i == 0 && j == 2) return 1;
    if (i == 0 && j == 3) return 2;
    if (i == 1 && j == 2) return 3;
    if (i == 1 && j == 3) return 4;
    if (i == 2 && j == 3) return 5;
    
    assert(false);
    return -1;
}


/*
    index1To2
    
    Edge midpoints are indexed via a tuple.  Turn the single index into a tuple.
*/
static
void index1To2(int index, int &i, int &j)
{
    if (index == 0) {i = 0; j = 1; return; }
    if (index == 1) {i = 0; j = 2; return; }
    if (index == 2) {i = 0; j = 3; return; }
    if (index == 3) {i = 1; j = 2; return; }
    if (index == 4) {i = 1; j = 3; return; }
    if (index == 5) {i = 2; j = 3; return; }
    
    assert(false);
}


/*
    distanceBetweenNodes
*/
static
double distanceBetweenNodes(const NodeCoords &node1, const NodeCoords &node2)
{
    double dist = 0.0;
    
    for (int i = 0; i < 3; i++) {
        double diff = node1.coords[i] - node2.coords[i];
        dist += diff * diff;
    }
    
    return sqrt(dist);
}


/*
    refineTet
    
    Refines the tetrahedron into 8 tetrahedra.
*/
static
void refineTet(const Tet &inTet, Tet outTets[8])
{
    NodeCoords origNodes[4];
    NodeCoords midNodes[6];
    NodeCoords oppositeMidNodes[3][2];
    
    
    // Copy inTet nodes to origNodes
    for (int index = 0; index < 4; index++) {
    for (int coord = 0; coord < 3; coord++) {
        origNodes[index].coords[coord] = inTet.nodes[index].coords[coord];
    }}
    
    
    // Calculate midpoint nodes on each edge
    for (int index = 0; index < 6; index++) {
    for (int coord = 0; coord < 3; coord++) {
        int i, j;
        index1To2(index, i, j);
        midNodes[index].coords[coord] = 
            0.5 * (origNodes[i].coords[coord] + origNodes[j].coords[coord]);
    }}
    
    
    // Mid nodes on opposite sides of middle octahedron
    // There are 3 pairs: (0,1) across from (2,3)
    //                    (0,2) across from (1,3)
    //                    (0,3) across from (1,2)
    oppositeMidNodes[0][0] = midNodes[index2To1(0,1)];
    oppositeMidNodes[0][1] = midNodes[index2To1(2,3)];
    
    oppositeMidNodes[1][0] = midNodes[index2To1(0,2)];
    oppositeMidNodes[1][1] = midNodes[index2To1(1,3)];
    
    oppositeMidNodes[2][0] = midNodes[index2To1(0,3)];
    oppositeMidNodes[2][1] = midNodes[index2To1(1,2)];
    
    
    // Find min distance between opposite midNodes
    double minDist = numeric_limits<double>::max();
    int minIndex = 0;
    for (int i = 0; i < 3; i++) {
        double dist = distanceBetweenNodes(oppositeMidNodes[i][0], oppositeMidNodes[i][1]);
        if (dist < minDist) {
            minIndex = i;
            minDist = dist;
        }
    }
    
    
    // 4 Tets from corners of original Tet
    outTets[0].nodes[0] = origNodes[0];
    outTets[0].nodes[1] = midNodes[index2To1(0,1)];
    outTets[0].nodes[2] = midNodes[index2To1(0,2)];
    outTets[0].nodes[3] = midNodes[index2To1(0,3)];
    
    outTets[1].nodes[0] = origNodes[1];
    outTets[1].nodes[1] = midNodes[index2To1(0,1)];
    outTets[1].nodes[2] = midNodes[index2To1(1,2)];
    outTets[1].nodes[3] = midNodes[index2To1(1,3)];
    
    outTets[2].nodes[0] = origNodes[2];
    outTets[2].nodes[1] = midNodes[index2To1(0,2)];
    outTets[2].nodes[2] = midNodes[index2To1(1,2)];
    outTets[2].nodes[3] = midNodes[index2To1(2,3)];
    
    outTets[3].nodes[0] = origNodes[3];
    outTets[3].nodes[1] = midNodes[index2To1(0,3)];
    outTets[3].nodes[2] = midNodes[index2To1(1,3)];
    outTets[3].nodes[3] = midNodes[index2To1(2,3)];
    
    
    // 4 Tets from inside octahedron
    if (minIndex == 0) {
        outTets[4].nodes[0] = midNodes[index2To1(2,3)];
        outTets[5].nodes[0] = midNodes[index2To1(2,3)];
        outTets[6].nodes[0] = midNodes[index2To1(0,1)];
        outTets[7].nodes[0] = midNodes[index2To1(0,1)];
    }
    else if (minIndex == 1) {
        outTets[4].nodes[0] = midNodes[index2To1(1,3)];
        outTets[5].nodes[0] = midNodes[index2To1(0,2)];
        outTets[6].nodes[0] = midNodes[index2To1(1,3)];
        outTets[7].nodes[0] = midNodes[index2To1(0,2)];
    }
    else if (minIndex == 2) {
        outTets[4].nodes[0] = midNodes[index2To1(1,2)];
        outTets[5].nodes[0] = midNodes[index2To1(0,3)];
        outTets[6].nodes[0] = midNodes[index2To1(0,3)];
        outTets[7].nodes[0] = midNodes[index2To1(1,2)];
    }
    
    outTets[4].nodes[1] = midNodes[index2To1(0,1)];
    outTets[4].nodes[2] = midNodes[index2To1(0,2)];
    outTets[4].nodes[3] = midNodes[index2To1(0,3)];
    
    outTets[5].nodes[1] = midNodes[index2To1(0,1)];
    outTets[5].nodes[2] = midNodes[index2To1(1,2)];
    outTets[5].nodes[3] = midNodes[index2To1(1,3)];
    
    outTets[6].nodes[1] = midNodes[index2To1(0,2)];
    outTets[6].nodes[2] = midNodes[index2To1(1,2)];
    outTets[6].nodes[3] = midNodes[index2To1(2,3)];
    
    outTets[7].nodes[1] = midNodes[index2To1(0,3)];
    outTets[7].nodes[2] = midNodes[index2To1(1,3)];
    outTets[7].nodes[3] = midNodes[index2To1(2,3)];
}


/*
    findFace
*/
uint64_t findFace(const FaceNodes &faceNodes, 
                const std::map<FaceNodes, uint64_t, FaceNodes> &faceMap)
{
    auto it = faceMap.find(faceNodes);
    if (it == faceMap.end())
        return numeric_limits<uint64_t>::max();
    return it->second;
}


/*
    findNode
*/
uint64_t findNode(const NodeCoords &nodeCoords, 
                const std::map<NodeCoords, uint64_t, NodeCoords> &nodeMap)
{
    auto it = nodeMap.find(nodeCoords);
    if (it == nodeMap.end())
        return numeric_limits<uint64_t>::max();
    return it->second;
}


/*
    refineMesh
*/
static
void refineMesh(const SerialMesh &coarseMesh, SerialMesh &refinedMesh)
{
    // Maps to quickly find already created nodes/faces
    std::map<NodeCoords, uint64_t, NodeCoords> nodeMap;
    std::map<FaceNodes, uint64_t, FaceNodes> faceMap;
    
    
    // Get graph diameter per coordinate
    // Used for a relative tolerance in implementing < for NodeCoords
    const double DOUBLE_MAX = numeric_limits<double>::max();
    double minCoords[3] = {DOUBLE_MAX, DOUBLE_MAX, DOUBLE_MAX};
    double maxCoords[3] = {-DOUBLE_MAX, -DOUBLE_MAX, -DOUBLE_MAX};
    for (uint64_t coarseCell = 0; coarseCell < coarseMesh.c_numCells; coarseCell++) {
        
        for (int lnode = 0; lnode < 4; lnode++) {
            uint64_t node = coarseMesh.c_cellData[coarseCell].boundingNodes[lnode];
            for (int i = 0; i < 3; i++) {
                double coord = coarseMesh.c_nodeData[node].coords[i];
                minCoords[i] = (coord < minCoords[i]) ? coord : minCoords[i];
                maxCoords[i] = (coord > maxCoords[i]) ? coord : maxCoords[i];
            }
        }
    }
    NodeCoords::graphDiameter = max(maxCoords[0] - minCoords[0], 
        max(maxCoords[1] - minCoords[1], maxCoords[2] - minCoords[2]));
        
    
    // Create the refined mesh.
    for (uint64_t coarseCell = 0; coarseCell < coarseMesh.c_numCells; coarseCell++) {
        
        if (coarseCell % (coarseMesh.c_numCells / 1000) == 0)
            printf("coarseCell: %" PRIu64 "\n", coarseCell);
        
        
        // Get coarse tet for this coarse cell.
        Tet coarseTet;
        for (int lnode = 0; lnode < 4; lnode++) {
            uint64_t node = coarseMesh.c_cellData[coarseCell].boundingNodes[lnode];
            for (int i = 0; i < 3; i++) {
                coarseTet.nodes[lnode].coords[i] = 
                    coarseMesh.c_nodeData[node].coords[i];
            }
        }
        
        
        // Refine the coarse tet.
        Tet refinedTets[8];
        refineTet(coarseTet, refinedTets);
        
        
        // Add cell, face, and node data by going through the refined tets.
        for (int tet = 0; tet < 8; tet++) {
            
            // Fill in cell data below
            SerialMesh::CellData cellData;
            
            
            // Bounding Nodes
            for (int lnode = 0; lnode < 4; lnode++) {
                
                SerialMesh::NodeData nodeData;
                NodeCoords nodeCoords;
                
                // Get the node coords
                for (int i = 0; i < 3; i++) {
                    double coord = refinedTets[tet].nodes[lnode].coords[i];
                    nodeData.coords[i] = coord;
                    nodeCoords.coords[i] = coord;
                }
                
                // Find the node if it was created earlier
                uint64_t node = findNode(nodeCoords, nodeMap);
                
                // If we didn't find the node
                if (node == numeric_limits<uint64_t>::max()) {
                    node = refinedMesh.c_nodeData.size();
                    pair<NodeCoords,uint64_t> nodePair(nodeCoords, node);
                    nodeMap.insert(nodePair);
                    cellData.boundingNodes[lnode] = node;
                    refinedMesh.c_nodeData.push_back(nodeData);
                }
                // If we found the node
                else {
                    cellData.boundingNodes[lnode] = node;
                }
            }
            
            
            // Bounding Faces
            for (int lface = 0; lface < 4; lface++) {
                
                SerialMesh::FaceData faceData;
                FaceNodes faceNodes;
                uint64_t cell = refinedMesh.c_cellData.size();
                
                // Get the bounding face nodes.
                // lnode: index for face,   lnode2: index for cell
                // Node i is across from face i in local indices
                for (int lnode = 0; lnode < 3; lnode++) {
                    int lnode2 = (lnode < lface) ? lnode : lnode + 1;
                    uint64_t node = cellData.boundingNodes[lnode2];
                    faceData.boundingNodes[lnode] = node;
                    faceNodes.nodes[lnode] = node;
                }
                
                // Find the face if it was created earlier
                uint64_t face = findFace(faceNodes, faceMap);
                
                // If we didn't find the face
                if (face == numeric_limits<uint64_t>::max()) {
                    face = refinedMesh.c_faceData.size();
                    pair<FaceNodes,uint64_t> facePair(faceNodes, face);
                    faceMap.insert(facePair);
                    cellData.boundingFaces[lface] = face;
                    faceData.boundingCells[0] = cell;
                    faceData.boundingCells[1] = numeric_limits<uint64_t>::max();
                    refinedMesh.c_faceData.push_back(faceData);
                }
                // If we found the face
                else {
                    cellData.boundingFaces[lface] = face;
                    refinedMesh.c_faceData[face].boundingCells[1] = cell;
                }
            }


            // Material Index
            cellData.materialIndex = 
                coarseMesh.c_cellData[coarseCell].materialIndex;
            
            
            // Add new cell
            refinedMesh.c_cellData.push_back(cellData);
        }
    }
    
    
    // Set number of cells, faces, nodes
    refinedMesh.c_numCells = refinedMesh.c_cellData.size();
    refinedMesh.c_numFaces = refinedMesh.c_faceData.size();
    refinedMesh.c_numNodes = refinedMesh.c_nodeData.size();
}


/*
    printMeshData
*/
static
void printMeshData(const SerialMesh &mesh)
{
    // Print number of interior/exterior faces
    uint64_t extFace = 0;
    uint64_t intFace = 0;
    for (uint64_t face = 0; face < mesh.c_numFaces; face++) {
        if (mesh.c_faceData[face].boundingCells[1] == numeric_limits<uint64_t>::max())
            extFace++;
        else
            intFace++;
    }
    
    printf("Cells: %" PRIu64 ", Nodes: %" PRIu64 ", Faces: %" PRIu64 "\n", mesh.c_numCells, 
           mesh.c_numNodes, mesh.c_numFaces);
    printf("Exterior Faces: %" PRIu64 "   Interior Faces: %" PRIu64 "\n", extFace, intFace);
}


/*
    main
*/
int main(int argc, char* argv[])
{
    SerialMesh coarseMesh;
    SerialMesh refinedMesh;
    string inputFile;
    string outputFile;
    
    
    // Print utility name
    printf("--- RefineSerialMesh Utility ---\n");
    
    
    // Get input/output files
    if (argc != 3) {
        printf("Incorrect number of arguments\n");
        printf("Usage: ./RefineSerialMesh.x <inputFile> <outputFile>\n");
        printf("\n\n\n");
        return 0;
    }
    inputFile = argv[1];
    outputFile = argv[2];
    printf("Input:  %s\n", inputFile.c_str());
    printf("Output: %s\n", outputFile.c_str());
    
    
    // Read in serial mesh, refine, and write out refined mesh.
    coarseMesh.read(inputFile);
    refineMesh(coarseMesh, refinedMesh);
    refinedMesh.write(outputFile);
    
    
    // Print coarse and fine mesh data
    printf("\nCoarse Mesh\n");
    printMeshData(coarseMesh);
    printf("\nRefined Mesh\n");
    printMeshData(refinedMesh);
    
    
    // Cleanup
    printf("\n\n\n");
    return 0;
}



