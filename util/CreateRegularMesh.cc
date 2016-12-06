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




// Each cube has 8 nodes, 6 cells, 18 faces when triangularized with tets
struct LocalToGlobalIndices
{
    uint64_t lgNodes[8];
    uint64_t lgFaces[18];
    uint64_t lgCells[6];
}


struct OwningHexIndices
{
    uint64_t nodes[8][3];
    uint64_t faces[18][3];
    uint64_t cells[6][3];
}



void createMesh(uint64_t nx, uint64_t ny, uint64_t nz, SerialMesh &mesh)
{
    const double dx = 100.0 / nx;
    const double dy = 100.0 / ny;
    const double dx = 100.0 / nz;
    
    LocalToGlobalIndices lgIndices[nx][ny][nz];
    OwningHexIndices owningHexIndices[nx][ny][nz];
    
    SerialMesh::FaceData localFaceData[18];
    SerialMesh::CellData localCellData[6];
    uint64_t boundingCellOffset[18][3];
    uint64_t otherHexIndex[18];
    
    
    // Local Cell Data
    localCellData[0].boundingNodes = {2,3,4,6};
    localCellData[1].boundingNodes = {3,4,6,7};
    localCellData[2].boundingNodes = {3,4,5,7};
    localCellData[3].boundingNodes = {0,2,3,4};
    localCellData[4].boundingNodes = {0,1,3,4};
    localCellData[5].boundingNodes = {1,3,4,5};
    
    localCellData[0].boundingFaces = {0,3,12,17};
    localCellData[1].boundingFaces = {1,4,12,13};
    localCellData[2].boundingFaces = {5,6,13,14};
    localCellData[3].boundingFaces = {2,9,16,17};
    localCellData[4].boundingFaces = {8,10,15,16};
    localCellData[5].boundingFaces = {7,11,14,15};
    
    
    // Local Face Data
    localFaceData[0].boundingNodes  = {2,3,6};
    localFaceData[1].boundingNodes  = {3,6,7};
    localFaceData[2].boundingNodes  = {0,2,4};
    localFaceData[3].boundingNodes  = {2,4,6};
    localFaceData[4].boundingNodes  = {4,6,7};
    localFaceData[5].boundingNodes  = {4,5,7};
    localFaceData[6].boundingNodes  = {3,5,7};
    localFaceData[7].boundingNodes  = {1,3,5};
    localFaceData[8].boundingNodes  = {0,1,3};
    localFaceData[9].boundingNodes  = {0,2,3};
    localFaceData[10].boundingNodes = {0,1,4};
    localFaceData[11].boundingNodes = {1,4,5};
    localFaceData[12].boundingNodes = {3,4,6};
    localFaceData[13].boundingNodes = {3,4,7};
    localFaceData[14].boundingNodes = {3,4,5};
    localFaceData[15].boundingNodes = {1,3,4};
    localFaceData[16].boundingNodes = {0,3,4};
    localFaceData[17].boundingNodes = {2,3,4};
    
    localFaceData[0].boundingCells  = {0,4};
    localFaceData[1].boundingCells  = {1,5};
    localFaceData[2].boundingCells  = {3,5};
    localFaceData[3].boundingCells  = {0,2};
    localFaceData[4].boundingCells  = {1,3};
    localFaceData[5].boundingCells  = {2,4};
    localFaceData[6].boundingCells  = {2,0};
    localFaceData[7].boundingCells  = {5,3};
    localFaceData[8].boundingCells  = {4,2};
    localFaceData[9].boundingCells  = {3,1};
    localFaceData[10].boundingCells = {4,0};
    localFaceData[11].boundingCells = {5,1};
    localFaceData[12].boundingCells = {0,1};
    localFaceData[13].boundingCells = {1,2};
    localFaceData[14].boundingCells = {2,5};
    localFaceData[15].boundingCells = {4,5};
    localFaceData[16].boundingCells = {3,4};
    localFaceData[17].boundingCells = {0,3};
    
    otherHexIndex[0]  = 10;
    otherHexIndex[1]  = 11;
    otherHexIndex[2]  = 7;
    otherHexIndex[3]  = 6;
    otherHexIndex[4]  = 9;
    otherHexIndex[5]  = 8;
    otherHexIndex[6]  = 3;
    otherHexIndex[7]  = 2;
    otherHexIndex[8]  = 5;
    otherHexIndex[9]  = 4;
    otherHexIndex[10] = 0;
    otherHexIndex[11] = 1;
    otherHexIndex[12] = SerialMesh::INVALID_INDEX;
    otherHexIndex[13] = SerialMesh::INVALID_INDEX;
    otherHexIndex[14] = SerialMesh::INVALID_INDEX;
    otherHexIndex[15] = SerialMesh::INVALID_INDEX;
    otherHexIndex[16] = SerialMesh::INVALID_INDEX;
    otherHexIndex[17] = SerialMesh::INVALID_INDEX;
    
    boundingCellOffset[0]  = { 0, +1,  0};
    boundingCellOffset[1]  = { 0, +1,  0};
    boundingCellOffset[2]  = {-1,  0,  0};
    boundingCellOffset[3]  = {-1,  0,  0};
    boundingCellOffset[4]  = { 0,  0, +1};
    boundingCellOffset[5]  = { 0,  0, +1};
    boundingCellOffset[6]  = {+1,  0,  0};
    boundingCellOffset[7]  = {+1,  0,  0};
    boundingCellOffset[8]  = { 0,  0, -1};
    boundingCellOffset[9]  = { 0,  0, -1};
    boundingCellOffset[10] = { 0, -1,  0};
    boundingCellOffset[11] = { 0, -1,  0};
    boundingCellOffset[12] = { 0,  0,  0};
    boundingCellOffset[13] = { 0,  0,  0};
    boundingCellOffset[14] = { 0,  0,  0};
    boundingCellOffset[15] = { 0,  0,  0};
    boundingCellOffset[16] = { 0,  0,  0};
    boundingCellOffset[17] = { 0,  0,  0};
    
    
    // lgIndices
    for (uint64_t i = 0; i < nx; i++) {
    for (uint64_t j = 0; j < ny; j++) {
    for (uint64_t k = 0; k < nz; k++) {
        
        // Cells
        for (int lcell = 0; lcell < 6; lcell++) {
            SerialMesh::CellData cellData;
            lgIndices[i][j][k].lgCells[lcell] = mesh.c_cellData.size();
            mesh.c_cellData.push_back(cellData);
            owningHexIndices[i][j][k].cells[lcell] = {i, j, k};
        }
        
        
        // Faces
        for (int lface = 0; lface < 18; lface++) {
            
            // Owning hex index
            uint64_t i1 = i;
            uint64_t j1 = j;
            uint64_t k1 = k;
            if ((lface ==  2 || lface ==  3) && i1 > 0)  i1--;
            if ((lface == 10 || lface == 11) && j1 > 0)  j1--;
            if ((lface ==  8 || lface ==  9) && k1 > 0)  k1--;
            owningHexIndices[i][j][k].faces[lface] = {i1, j1, k1};
            
            // If this hex own face, add the face
            if (i1 == i && j1 == j && k1 == k) {
                SerialMesh::FaceData faceData;
                lgIndices[i][j][k].lgFaces[lface] = mesh.c_faceData.size();
                mesh.c_faceData.push_back(faceData);
            }
            
            // Otherwise just update lgIndices
            else {
                localIndex = otherHexIndex[lface];
                lgIndices[i][j][k].lgFaces[lface] = 
                    lgIndices[i1][j1][k1].lgFaces[localIndex];
            }
        }
        
        
        // Nodes
        for (int lnode = 0; lnode < 8; lnode++) {
            
            // Owning hex index
            uint64_t i1 = i;
            uint64_t j1 = j;
            uint64_t k1 = k;
            if ((lnode % 2) < 1 && i1 > 0)  i1--;
            if ((lnode % 4) < 2 && j1 > 0)  j1--;
            if ((lnode % 8) < 4 && k1 > 0)  k1--;
            owningHexIndices[i][j][k].nodes[lnode] = {i1, j1, k1};
            
            // If this hex owns node, add the node
            if (i1 == i && j1 == j && k1 == k) {
                SerialMesh::NodeData nodeData;
                lgIndices[i][j][k].lgNodes[lnode] = mesh.c_nodeData.size();
                mesh.c_nodeData.push_back(nodeData);
            }
            // Otherwise just update lgIndices
            else {
                localIndex = lnode;
                if (i1 != i)  localIndex += 1;
                if (j1 != j)  localIndex += 2;
                if (k1 != k)  localIndex += 4;
                lgIndices[i][j][k].lgFaces[lface] = 
                    lgIndices[i1][j1][k1].lgFaces[localIndex];
            }
        }
    }}}
    
    
    // Node Data
    for (uint64_t i = 0; i < nx; i++) {
    for (uint64_t j = 0; j < ny; j++) {
    for (uint64_t k = 0; k < nz; k++) {
        
        int ox[8] = {-1,  0, -1,  0, -1,  0, -1,  0};
        int oy[8] = {-1, -1,  0,  0, -1, -1,  0,  0};
        int oz[8] = {-1, -1, -1, -1,  0,  0,  0,  0};
        
        
        // Node Data
        for (int lnode = 0; lnode < 8; lnode++) {
            
            uint64_t nodeIndex = lgIndices[i][j][k].lgCells[lnode]
            SerialMesh::NodeData &nodeData = mesh.c_nodeData[nodeIndex];
            
            // Owning hex of node
            int i1 = i + ox[lnode];
            int j1 = j + oy[lnode];
            int k1 = k + oz[lnode];
            
            
            // Coords of node
            nodeData.coords = {(i1 + 1) * dx, (j1 + 1) * dy, (k1 + 1) * dz};
        }
    }}}
    
        
    // Face Data
    for (uint64_t i = 0; i < nx; i++) {
    for (uint64_t j = 0; j < ny; j++) {
    for (uint64_t k = 0; k < nz; k++) {
        
        SerialMesh::FaceData faceData[18];
        int ox[18] = { 0,  0, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0};
        int oy[18] = { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1, -1,  0,  0,  0,  0,  0,  0};
        int oz[18] = { 0,  0,  0,  0,  0,  0,  0,  0, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0};
        
        
        // Face Data
        for (int lface = 0; lface < 18; lface++) {
            
            // Owning hex of face
            int i1 = i + ox[lface];
            int j1 = j + oy[lface];
            int k1 = k + oz[lface];
            
            if (i1 < 0)  i1 = 0;
            if (j1 < 0)  j1 = 0;
            if (k1 < 0)  k1 = 0;
            
            
            // Bounding nodes
            for (int lnode = 0; lnode < 3; lnode++) {
                localNodeIndex = localFaceData[lface].boundingNodes[lnode];
                faceData[lface].boundingNodes[lnode] = 
                    lgIndices[i][j][k].lgNodes[localNodeIndex];
            }
            
            
            // Bounding cells
            for (int lcell = 0; lcell < 2; lcell++) {
                int ci = i;
                int cj = j;
                int ck = k;
                
                if (lcell == 1) {
                    ci += boundingCellOffset[lcell][0];
                    cj += boundingCellOffset[lcell][1];
                    ck += boundingCellOffset[lcell][2];
                }
                
                localCellIndex = localFaceData.boundingCells[lcell];
                if (ci < 0 || cj < 0 || ck < 0) {
                    faceData[lface].boundingCells[lcell] = 
                        SerialMesh::INVALID_INDEX;
                }
                else {
                    faceData[lface].boundingCells[lcell] = 
                        lgIndices[ci][cj][ck].lgCells[localCellIndex];
                }
            }
            
            
            // Update mesh's FaceData and lgIndices for faces
            if (i1 < i || j1 < j || k1 < k) {
                lgIndices[i][j][k].lgFaces[lface] = 
                    lgIndices[i1][j1][k1].lgFaces[lface];
            }
            else {
                lgIndices[i][j][k].lgFaces[lface] = mesh.c_faceData.size();
                mesh.c_faceData.push_back(faceData[lface]);
            }
        }
    }}}
    
    
    // Cell Data
    for (uint64_t i = 0; i < nx; i++) {
    for (uint64_t j = 0; j < ny; j++) {
    for (uint64_t k = 0; k < nz; k++) {
        
        SerialMesh::FaceData faceData[18];
        int ox[18] = { 0,  0, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0};
        int oy[18] = { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1, -1,  0,  0,  0,  0,  0,  0};
        int oz[18] = { 0,  0,  0,  0,  0,  0,  0,  0, -1, -1,  0,  0,  0,  0,  0,  0,  0,  0};
        
        
        // Face Data
        for (int lcell = 0; lcell < 6; lcell++) {
            
            uint64_t cellIndex = lgIndices[i][j][k].lgCells[lcell]
            SerialMesh::CellData &cellData = mesh.c_cellData[cellIndex];
            
            
            // Bounding nodes
            for (int lnode = 0; lnode < 4; lnode++) {
                localNodeIndex = localCellData[lcell].boundingNodes[lnode];
                cellData.boundingNodes[lnode] = 
                    lgIndices[i][j][k].lgNodes[localNodeIndex];
            }
            
            
            // Bounding faces
            for (int lface = 0; lface < 4; lface++) {
                localFaceIndex = localCellData[lcell].boundingFaces[lface];
                cellData.boundingFaces[lface] = 
                    lgIndices[i][j][k].lgFaces[localFaceIndex];
            }
        }
    }}}
    
    
    // Number of mesh elements
    mesh.c_numCells = mesh.c_cellData.size();
    mesh.c_numFaces = mesh.c_faceData.size();
    mesh.c_numNodes = mesh.c_nodeData.size();
}


