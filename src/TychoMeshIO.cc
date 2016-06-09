//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh/TychoMeshIO.cc
 * \author Kris Garrett
 * \date   Mon Nov 9 2015
 * \brief  \link rtt_mesh::TychoMeshIO TychoMeshIO \endlink IO for TychoMesh
 */
//---------------------------------------------------------------------------//
// $Id: TychoMeshIO.cc,v 1.6 2000/04/06 21:45:16 pautz Exp $
//---------------------------------------------------------------------------//

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
    order 3 vertices
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
    Read in the TychoMesh
*/
void TychoMesh::readTychoMesh(const std::string &filename)
{
    UINT side;
    Mat2<UINT> cellFaceHandles;
    ParallelMesh::PartitionData partData;
    
    
    // Read mesh
    ParallelMesh::readInParallel(filename, partData);
    
    
    // c_nCells, c_nNodes
    c_nCells = partData.numCells;
    c_nNodes = partData.numNodes;
    
    
    // c_nodeCoords
    c_nodeCoords = Mat2<double>(c_nNodes, 3);
    for(UINT node = 0; node < c_nNodes; node++) {
        c_nodeCoords(node, 0) = partData.nodeData[node].coords[0];
        c_nodeCoords(node, 1) = partData.nodeData[node].coords[1];
        c_nodeCoords(node, 2) = partData.nodeData[node].coords[2];
    }
    
    
    // c_cellNodes
    c_cellNodes = Mat2<UINT>(c_nCells, g_nVrtxPerCell);
    for(UINT cell = 0; cell < c_nCells; cell++) {
    for(UINT vrtx = 0; vrtx < g_nVrtxPerCell; vrtx++) {
        c_cellNodes(cell, vrtx) = partData.cellData[cell].boundingNodes[vrtx];
    }}
    
    
    // c_adjCell, c_adjFace, c_nSides
    c_nSides = 0;
    c_adjCell = Mat2<UINT>(c_nCells, g_nFacePerCell);
    c_adjFace = Mat2<UINT>(c_nCells, g_nFacePerCell);
    cellFaceHandles = Mat2<UINT>(c_nCells, g_nFacePerCell);
    for(UINT cell = 0; cell < c_nCells; cell++) {
    for(UINT faceRangeIndex = 0; faceRangeIndex < g_nFacePerCell; faceRangeIndex++) {
        UINT face = partData.cellData[cell].boundingFaces[faceRangeIndex];
        UINT lface = getLFaceIndex(partData.cellData[cell].boundingNodes, 
                                  partData.faceData[face].boundingNodes);
        
        cellFaceHandles(cell, lface) = faceRangeIndex;
        
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
            c_adjFace(cell, lface) = getLFaceIndex(partData.cellData[adjCell].boundingNodes, 
                                                   partData.faceData[face].boundingNodes);
        }
    }}
    
    
    // c_sideCell, c_side, c_lGSides, c_gLSides
    side = 0;
    c_sideCell = Mat1<UINT>(c_nSides);
    c_side = Mat2<UINT>(c_nCells, g_nFacePerCell);
    c_lGSides = Mat1<UINT>(c_nSides);
    for(UINT cell = 0; cell < c_nCells; cell++) {
    for(UINT lface = 0; lface < g_nFacePerCell; lface++) {
        if(c_adjCell(cell, lface) == BOUNDARY_FACE) {
            UINT face = partData.cellData[cell].boundingFaces[cellFaceHandles(cell,lface)];
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
    c_adjProc = Mat2<UINT>(c_nCells, g_nFacePerCell);
    for(UINT cell = 0; cell < c_nCells; cell++) {
    for(UINT lface = 0; lface < g_nFacePerCell; lface++) {
        UINT face = partData.cellData[cell].boundingFaces[cellFaceHandles(cell,lface)];
        UINT proc1 = partData.faceData[face].partition[0];
        UINT proc2 = partData.faceData[face].partition[1];
        UINT thisProc = Comm::rank();
        
        Assert(proc1 == thisProc || proc2 == thisProc);
        c_adjProc(cell, lface) = (proc1 == thisProc) ? proc2 : proc1;
        
        if (c_adjProc(cell, lface) == ParallelMesh::INVALID_INDEX) 
            c_adjProc(cell, lface) = BAD_RANK;
    }}
    
    
    // c_faceToCellVrtx, c_cellToFaceVrtx
    c_faceToCellVrtx = Mat3<UINT>(c_nCells, g_nFacePerCell, g_nVrtxPerFace);
    c_cellToFaceVrtx = Mat3<UINT>(c_nCells, g_nFacePerCell, g_nVrtxPerCell);
    for(UINT cell = 0; cell < c_nCells; cell++) {
    for(UINT lface = 0; lface < g_nFacePerCell; lface++) {
        UINT order[3];
        UINT face = partData.cellData[cell].boundingFaces[cellFaceHandles(cell,lface)];
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
    for(UINT  cell = 0; cell < c_nCells; cell++) {
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
    
    c_adjCellFromSide = Mat1<UINT>(c_nSides);
    c_adjFaceFromSide = Mat1<UINT>(c_nSides);
    for(UINT cell = 0; cell < c_nCells; cell++) {
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


