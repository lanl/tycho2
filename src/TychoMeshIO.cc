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

#include "Mat.hh"
#include "TychoMesh.hh"
#include "Global.hh"
#include "Assert.hh"
#include "Comm.hh"
#include "ParallelMesh.hh"
#include <memory>
#include <stddef.h>
#include <utility>

using namespace std;


/*
    isV1LessThanV2

    An ordering on vertices so that vertices from the same face
    but different cells contain the same ordering.
    This is especially important for neighboring cells in different
    partitions of the mesh.
*/
static
bool isV1LessThanV2(const double v1[3], const double v2[3])
{
    for(UINT i = 0; i < 3; i++) {
        if(v1[i] < v2[i]) return true;
        if(v1[i] > v2[i]) return false;
    }
    
    // Shouldn't get here if the points are different.
    Assert(true);
    return true;
}


/*
    order3Vertices

    Puts the vertices of a face in order.
*/
static 
void order3Vertices(const double coords0[3], const double coords1[3], 
                    const double coords2[3], UINT order[3])
{
    if (isV1LessThanV2(coords0, coords1) && isV1LessThanV2(coords0, coords2)) {
        order[0] = 0;
        if (isV1LessThanV2(coords1, coords2)) {
            order[1] = 1;
            order[2] = 2;
        }
        else {
            order[1] = 2;
            order[2] = 1;
        }
        return;
    }
    if (isV1LessThanV2(coords1, coords0) && isV1LessThanV2(coords1, coords2)) {
        order[0] = 1;
        if (isV1LessThanV2(coords0, coords2)) {
            order[1] = 0;
            order[2] = 2;
        }
        else {
            order[1] = 2;
            order[2] = 0;
        }
        return;
    }
    if (isV1LessThanV2(coords2, coords0) && isV1LessThanV2(coords2, coords1)) {
        order[0] = 2;
        if (isV1LessThanV2(coords0, coords1)) {
            order[1] = 0;
            order[2] = 1;
        }
        else {
            order[1] = 1;
            order[2] = 0;
        }
        return;
    }
}


/*
    getLFaceIndex

    Gets the face index as the index of the vertex not in the face.
*/
static
UINT getLFaceIndex(uint64_t cellBoundingNodes[4], uint64_t faceBoundingNodes[3])
{
    // Find which vrtx is not in the face.
    for(UINT cnode = 0; cnode < g_nVrtxPerCell; cnode++) {
        if(cellBoundingNodes[cnode] != faceBoundingNodes[0] && 
           cellBoundingNodes[cnode] != faceBoundingNodes[1] && 
           cellBoundingNodes[cnode] != faceBoundingNodes[2])
        {
            return cnode;
        }
    }
    
    
    // Should never reach this point.
    Assert(true);
    return SIZE_MAX;
}


/*
    readTychoMesh

    Reads in the TychoMesh
*/
void TychoMesh::readTychoMesh(const std::string &filename)
{
    UINT side;
    Mat2<UINT> cellFaceHandles;
    ParallelMesh::PartitionData partData;
    
    
    // Read mesh
    ParallelMesh::readInParallel(filename, partData);
    
    
    // g_nCells, c_nNodes
    g_nCells = partData.numCells;
    c_nNodes = partData.numNodes;
    
    
    // c_nodeCoords
    c_nodeCoords.resize(c_nNodes, 3);
    for(UINT node = 0; node < c_nNodes; node++) {
        c_nodeCoords(node, 0) = partData.nodeData[node].coords[0];
        c_nodeCoords(node, 1) = partData.nodeData[node].coords[1];
        c_nodeCoords(node, 2) = partData.nodeData[node].coords[2];
    }


    // c_lGCells, c_cellMaterial
    c_lGCells.resize(g_nCells);
    c_cellMaterial.resize(g_nCells);
    for (UINT cell = 0; cell < g_nCells; cell++) {
        c_lGCells(cell) = partData.cellData[cell].globalID;
        c_cellMaterial(cell) = partData.cellData[cell].materialIndex;
    }
    
    
    // c_cellNodes
    c_cellNodes.resize(g_nCells, g_nVrtxPerCell);
    for(UINT cell = 0; cell < g_nCells; cell++) {
    for(UINT vrtx = 0; vrtx < g_nVrtxPerCell; vrtx++) {
        c_cellNodes(cell, vrtx) = partData.cellData[cell].boundingNodes[vrtx];
    }}
    
    
    // c_adjCell, c_adjFace, c_nSides
    // This is tricky because the face indexing from ParallelMesh is not
    // the same as the face indexing for a TychoMesh.
    c_nSides = 0;
    c_adjCell.resize(g_nCells, g_nFacePerCell);
    c_adjFace.resize(g_nCells, g_nFacePerCell);
    cellFaceHandles.resize(g_nCells, g_nFacePerCell);
    for(UINT cell = 0; cell < g_nCells; cell++) {
    for(UINT faceIndex = 0; faceIndex < g_nFacePerCell; faceIndex++) {
        UINT face = partData.cellData[cell].boundingFaces[faceIndex];
        UINT lface = getLFaceIndex(partData.cellData[cell].boundingNodes, 
                                   partData.faceData[face].boundingNodes);
        
        cellFaceHandles(cell, lface) = faceIndex;
        
        if(partData.faceData[face].boundaryType != ParallelMesh::NotBoundary) {
            c_adjCell(cell, lface) = BOUNDARY_FACE;
            c_adjFace(cell, lface) = BOUNDARY_FACE;
            c_nSides++;
        }
        else {
            UINT cell1 = partData.faceData[face].boundingCells[0];
            UINT cell2 = partData.faceData[face].boundingCells[1];
            UINT adjCell = (cell1 == cell) ? cell2 : cell1;
            
            c_adjCell(cell, lface) = adjCell;
            c_adjFace(cell, lface) = 
                getLFaceIndex(partData.cellData[adjCell].boundingNodes, 
                              partData.faceData[face].boundingNodes);
        }
    }}
    
    
    // c_sideCell, c_side, c_lGSides, c_gLSides
    side = 0;
    c_sideCell.resize(c_nSides);
    c_side.resize(g_nCells, g_nFacePerCell);
    c_lGSides.resize(c_nSides);
    for(UINT cell = 0; cell < g_nCells; cell++) {
    for(UINT lface = 0; lface < g_nFacePerCell; lface++) {
        if(c_adjCell(cell, lface) == BOUNDARY_FACE) {
            UINT cellFaceHandle = cellFaceHandles(cell, lface);
            UINT face = partData.cellData[cell].boundingFaces[cellFaceHandle];
            UINT gside = partData.faceData[face].globalID;
            
            c_sideCell(side) = cell;
            c_side(cell, lface) = side;
            c_lGSides(side) = gside;
            c_gLSides.insert(make_pair(gside, side));
            side++;
        }
        else {
            c_side(cell, lface) = NOT_BOUNDARY_FACE;
        }
    }}
    
    
    // c_adjProc
    c_adjProc.resize(g_nCells, g_nFacePerCell);
    for(UINT cell = 0; cell < g_nCells; cell++) {
    for(UINT lface = 0; lface < g_nFacePerCell; lface++) {
        UINT cellFaceHandle = cellFaceHandles(cell, lface);
        UINT face = partData.cellData[cell].boundingFaces[cellFaceHandle];
        UINT proc1 = partData.faceData[face].partition[0];
        UINT proc2 = partData.faceData[face].partition[1];
        UINT thisProc = Comm::rank();
        
        Assert(proc1 == thisProc || proc2 == thisProc);
        c_adjProc(cell, lface) = (proc1 == thisProc) ? proc2 : proc1;
        
        if (c_adjProc(cell, lface) == ParallelMesh::INVALID_INDEX) 
            c_adjProc(cell, lface) = BAD_RANK;
    }}
    
    
    // c_faceToCellVrtx, c_cellToFaceVrtx
    c_faceToCellVrtx.resize(g_nCells, g_nFacePerCell, g_nVrtxPerFace);
    c_cellToFaceVrtx.resize(g_nCells, g_nFacePerCell, g_nVrtxPerCell);
    for(UINT cell = 0; cell < g_nCells; cell++) {
    for(UINT lface = 0; lface < g_nFacePerCell; lface++) {
        UINT order[3];
        UINT cellFaceHandle = cellFaceHandles(cell, lface);
        UINT face = partData.cellData[cell].boundingFaces[cellFaceHandle];
        UINT node0 = partData.faceData[face].boundingNodes[0];
        UINT node1 = partData.faceData[face].boundingNodes[1];
        UINT node2 = partData.faceData[face].boundingNodes[2];
        
        order3Vertices(partData.nodeData[node0].coords,
                       partData.nodeData[node1].coords,
                       partData.nodeData[node2].coords,
                       order);
        for(UINT fvrtx = 0; fvrtx < g_nVrtxPerFace; fvrtx++) {
        for(UINT cvrtx = 0; cvrtx < g_nVrtxPerCell; cvrtx++) {
            if(partData.faceData[face].boundingNodes[order[fvrtx]] == 
               partData.cellData[cell].boundingNodes[cvrtx])
            {
                c_faceToCellVrtx(cell, lface, order[fvrtx]) = cvrtx;
                c_cellToFaceVrtx(cell, lface, cvrtx) = order[fvrtx];
            }
        }}
    }}
    
    
    // c_adjCellFromSide, c_adjFaceFromSide
    vector<MPI_Request> mpiRequests;
    for(UINT  cell = 0; cell < g_nCells; cell++) {
    for(UINT face = 0; face < g_nFacePerCell; face++) {
        UINT adjProc = c_adjProc(cell, face);
        UINT adjCell = c_adjCell(cell, face);
        
        if(adjCell == BOUNDARY_FACE && adjProc != BAD_RANK) {
            MPI_Request request;
            UINT cellFace[2];
            cellFace[0] = cell;
            cellFace[1] = face;
            MPI_Isend(cellFace, 2 * sizeof(UINT), MPI_CHAR, adjProc, 0, 
                      MPI_COMM_WORLD, &request);
            mpiRequests.push_back(request);
        }
    }}
    
    c_adjCellFromSide.resize(c_nSides);
    c_adjFaceFromSide.resize(c_nSides);
    for(UINT cell = 0; cell < g_nCells; cell++) {
    for(UINT face = 0; face < g_nFacePerCell; face++) {
        UINT adjProc = c_adjProc(cell, face);
        UINT adjCell = c_adjCell(cell, face);
        
        if(adjCell == BOUNDARY_FACE && adjProc != BAD_RANK) {
            UINT cellFace[2];
            MPI_Recv(cellFace, 2 * sizeof(UINT), MPI_CHAR, adjProc, 0, 
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            c_adjCellFromSide(c_side(cell, face)) = cellFace[0];
            c_adjFaceFromSide(c_side(cell, face)) = cellFace[1];
        }
    }}
    
    MPI_Waitall(mpiRequests.size(), &mpiRequests[0], MPI_STATUSES_IGNORE);
}


