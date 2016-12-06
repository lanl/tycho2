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

#include <moab/Core.hpp>
#include <string>
#include <cstdio>
#include <cassert>
#include "SerialMesh.hh"

using namespace moab;
using namespace std;

static const EntityHandle ENTIRE_MESH = 0;
static const int ONE_ENTITY = 1;
static const int CELL_DIMENSION = 3;
static const int FACE_DIMENSION = 2;
static const int NODE_DIMENSION = 0;
static const bool CREATE_MISSING = true;


enum ProblemType
{
    PROBLEM_ONE_MATERIAL, PROBLEM_TWO_NESTED_MATERIALS, PROBLEM_NONE
};


/*
    getCoords
*/
static
void getCoords(Core &moab, EntityHandle nodeHandle, double coords[3])
{
    ErrorCode moabStatus;
    moabStatus = moab.get_coords(&nodeHandle, ONE_ENTITY, coords);
    assert(moabStatus == MB_SUCCESS);
}


/*
    createSerialMesh
*/
static
void createSerialMesh(Core &moab, SerialMesh &mesh, ProblemType problemType)
{
    ErrorCode moabStatus;
    Range cellRange, faceRange, nodeRange;
    map<EntityHandle, uint64_t> cellMap, faceMap, nodeMap;
    
    
    // Get ranges.
    moabStatus = moab.get_entities_by_dimension(ENTIRE_MESH, CELL_DIMENSION, cellRange);
    assert(moabStatus == MB_SUCCESS);
    moabStatus = moab.get_entities_by_dimension(ENTIRE_MESH, FACE_DIMENSION, faceRange);
    assert(moabStatus == MB_SUCCESS);
    moabStatus = moab.get_entities_by_dimension(ENTIRE_MESH, NODE_DIMENSION, nodeRange);
    assert(moabStatus == MB_SUCCESS);
    
    
    // Number of cells, faces, nodes
    mesh.c_numCells = cellRange.size();
    mesh.c_numFaces = faceRange.size();
    mesh.c_numNodes = nodeRange.size();
    mesh.c_cellData.resize(mesh.c_numCells);
    mesh.c_faceData.resize(mesh.c_numFaces);
    mesh.c_nodeData.resize(mesh.c_numNodes);
    printf("Num Cells: %llu\n", mesh.c_numCells);
    printf("Num Faces: %llu\n", mesh.c_numFaces);
    printf("Num Nodes: %llu\n", mesh.c_numNodes);
    
    
    // Make map of EntityHandle to cell, face, node
    for (uint64_t cell = 0; cell < mesh.c_numCells; cell++) {
        cellMap[cellRange[cell]] = cell;
    }
    for (uint64_t face = 0; face < mesh.c_numFaces; face++) {
        faceMap[faceRange[face]] = face;
    }
    for (uint64_t node = 0; node < mesh.c_numNodes; node++) {
        nodeMap[nodeRange[node]] = node;
    }
    
    
    // Put in cell data
    for (uint64_t cell = 0; cell < mesh.c_numCells; cell++) {
        EntityHandle cellHandle = cellRange[cell];
        Range localFaceRange, localNodeRange;
        SerialMesh::CellData &cellData = mesh.c_cellData[cell];
        
        // Faces
        moabStatus = moab.get_adjacencies(&cellHandle, ONE_ENTITY, 
            FACE_DIMENSION, CREATE_MISSING, localFaceRange);
        assert(moabStatus == MB_SUCCESS);
        assert(localFaceRange.size() == 4);
        
        for (uint64_t lface = 0; lface < localFaceRange.size(); lface++) {
            cellData.boundingFaces[lface] = 
                faceMap.find(localFaceRange[lface])->second;
        }
        
        // Nodes
        moabStatus = moab.get_adjacencies(&cellHandle, ONE_ENTITY, 
            NODE_DIMENSION, CREATE_MISSING, localNodeRange);
        assert(moabStatus == MB_SUCCESS);
        assert(localNodeRange.size() == 4);
        
        for (uint64_t lnode = 0; lnode < localNodeRange.size(); lnode++) {
            cellData.boundingNodes[lnode] = 
                nodeMap.find(localNodeRange[lnode])->second;
        }

        // Material Index (default to 1)
        cellData.materialIndex = 1;
        if (problemType == PROBLEM_TWO_NESTED_MATERIALS) {
            double centroid[3] = {0.0, 0.0, 0.0};
            double coords[3]   = {0.0, 0.0, 0.0};

            for (uint64_t lnode = 0; lnode < localNodeRange.size(); lnode++) {
                getCoords(moab, localNodeRange[lnode], coords);
                centroid[0] += coords[0];
                centroid[1] += coords[1];
                centroid[2] += coords[2];
            }

            centroid[0] = centroid[0] / localNodeRange.size();
            centroid[1] = centroid[1] / localNodeRange.size();
            centroid[2] = centroid[2] / localNodeRange.size();

            if (centroid[0] > 25.0 && centroid[0] < 75.0 && 
                centroid[1] > 25.0 && centroid[1] < 75.0 && 
                centroid[2] > 25.0 && centroid[2] < 75.0)
            {
                cellData.materialIndex = 2;
            }
        }
    }
    
    
    // Put in face data
    for (uint64_t face = 0; face < mesh.c_numFaces; face++) {
        EntityHandle faceHandle = faceRange[face];
        Range localCellRange, localNodeRange;
        SerialMesh::FaceData &faceData = mesh.c_faceData[face];
        
        // Cells
        moabStatus = moab.get_adjacencies(&faceHandle, ONE_ENTITY, 
            CELL_DIMENSION, CREATE_MISSING, localCellRange);
        assert(moabStatus == MB_SUCCESS);
        assert(localCellRange.size() == 1 || localCellRange.size() == 2);
        
        for (uint64_t lcell = 0; lcell < localCellRange.size(); lcell++) {
            faceData.boundingCells[lcell] = 
                cellMap.find(localCellRange[lcell])->second;
        }
        
        if (localCellRange.size() == 1) {
            faceData.boundingCells[1] = numeric_limits<uint64_t>::max();
        }
        
        // Nodes
        moabStatus = moab.get_adjacencies(&faceHandle, ONE_ENTITY, 
            NODE_DIMENSION, CREATE_MISSING, localNodeRange);
        assert(moabStatus == MB_SUCCESS);
        assert(localNodeRange.size() == 3);
        
        for (uint64_t lnode = 0; lnode < localNodeRange.size(); lnode++) {
            faceData.boundingNodes[lnode] = 
                nodeMap.find(localNodeRange[lnode])->second;
        }
    }
    
    
    // Put in node data
    for (uint64_t node = 0; node < mesh.c_numNodes; node++) {
        double coords[3];
        //EntityHandle nodeHandle = nodeRange[node];
        SerialMesh::NodeData &nodeData = mesh.c_nodeData[node];
        
        //moabStatus = moab.get_coords(&nodeHandle, ONE_ENTITY, coords);
        //assert(moabStatus == MB_SUCCESS);
        getCoords(moab, nodeRange[node], coords);
        
        nodeData.coords[0] = coords[0];
        nodeData.coords[1] = coords[1];
        nodeData.coords[2] = coords[2];
    }
}


/*
    moabCreateAllFaces
*/
static
void moabCreateAllFaces(Core &moab)
{
    // Get range for cells.
    Range cellRange;
    ErrorCode moabStatus;
    moabStatus = moab.get_entities_by_dimension(ENTIRE_MESH, CELL_DIMENSION, cellRange);
    assert(moabStatus == MB_SUCCESS);
    
    // Create faces by getting adjacent faces to each cell.
    for (uint64_t cell = 0; cell < cellRange.size(); cell++) {
        Range faceRange;
        EntityHandle cellHandle = cellRange[cell];
        moabStatus = moab.get_adjacencies(&cellHandle, ONE_ENTITY, 
            FACE_DIMENSION, CREATE_MISSING, faceRange);
        assert(moabStatus == MB_SUCCESS);
    }
}


/*
    main
*/
int main(int argc, char* argv[])
{
    Core *moab = NULL;
    string inputFile;
    string outputFile;
    ErrorCode moabStatus;
    SerialMesh mesh;
    ProblemType problemType;
    int menuChoice;
    
    
    // Print utility name
    printf("--- MoabToSerialMesh Utility ---\n");
    
    
    // Init MOAB
    moab = new Core();
    assert(moab != NULL);
    
    
    // Get input/output files
    if (argc != 3) {
        printf("Incorrect number of arguments\n");
        printf("Usage: ./MoabToSerialMesh.x <inputFile> <outputFile>\n");
        printf("\n\n\n");
        return 0;
    }
    inputFile = argv[1];
    outputFile = argv[2];
    
    
    // Ask for material problem type
    printf("1 - Uniform material\n");
    printf("2 - Two nested materials\n"); 
    scanf("%d", &menuChoice);
    switch (menuChoice) {
        case 1:
            problemType = PROBLEM_ONE_MATERIAL;
            break;
        case 2:
            problemType = PROBLEM_TWO_NESTED_MATERIALS;
            break;
        default:
            printf("Menu option not recognized\n");
            printf("Exiting without processing mesh\n");
            problemType = PROBLEM_NONE;
            break;
    }
    

    // Do mesh conversion
    if (problemType != PROBLEM_NONE) {
        
        // Load file
        printf("Loading file %s ...\n", inputFile.c_str());
        moabStatus = moab->load_file(inputFile.c_str());
        assert(moabStatus == MB_SUCCESS);

        // Create mesh and write it to file
        moabCreateAllFaces(*moab);
        createSerialMesh(*moab, mesh, problemType);
        mesh.write(outputFile);
    }

    
    // Clean up
    printf("\n\n\n");
    return 0;
}



