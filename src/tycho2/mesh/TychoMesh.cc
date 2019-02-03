/*
Copyright (c) 2016, Los Alamos National Security, LLC
All rights reserved.
*/

#include "TychoMesh.hh"
#include "Global.hh"
#include "Quadrature.hh"
#include <stdio.h>
#include <math.h>

using namespace std;



/*
    Constructor
*/
TychoMesh::TychoMesh(const std::string &filename)
{
    // Read in the TychoMesh
    readTychoMesh(filename);
    
    
    // Cell volumes and face areas
    c_cellVolume.resize(g_nCells);
    c_faceArea.resize(g_nCells, g_nFacePerCell);
    for(UINT cell = 0; cell < g_nCells; cell++) {
        c_cellVolume(cell) = getTetVolume(getCellVrtxCoords(cell));
        for(UINT face = 0; face < g_nFacePerCell; face++) {
            c_faceArea(cell, face) = 
                getTriangleArea(getFaceVrtxCoords(cell, face));
        }
    }
    
    
    // Calculates Omega dot N for every face
    c_omegaDotN.resize(g_nAngles, g_nCells, g_nFacePerCell);
    for (UINT angle = 0; angle < g_nAngles; ++angle) {
        const vector<double> omega = g_quadrature->getOmega(angle);
        for (UINT cell = 0; cell < g_nCells; ++cell) {
        for (UINT face = 0; face < g_nFacePerCell; ++face) {
            vector<double> normal = getNormal(getFaceVrtxCoords(cell, face), 
                                              getCellVrtxCoords(cell));
            c_omegaDotN(angle, cell, face) = omega[0] * normal[0] + 
                                             omega[1] * normal[1] + 
                                             omega[2] * normal[2];
        }}
    }
    
    
    // CHECK getCellToFaceVrtx and getFaceToCellVrtx
    for(UINT cell = 0; cell < g_nCells; cell++) {
    for(UINT face = 0; face < g_nFacePerCell; face++) {
        for(UINT cvrtx = 0; cvrtx < g_nVrtxPerCell; cvrtx++) {
            if(cvrtx != face) {
                UINT fvrtx = getCellToFaceVrtx(cell, face, cvrtx);
                UINT cvrtx1 = getFaceToCellVrtx(cell, face, fvrtx);
                Insist(cvrtx == cvrtx1, "Error 1");
            }
        }
        for(UINT fvrtx = 0; fvrtx < g_nVrtxPerFace; fvrtx++) {
            UINT cvrtx = getFaceToCellVrtx(cell, face, fvrtx);
            UINT fvrtx1 = getCellToFaceVrtx(cell, face, cvrtx);
            Insist(fvrtx == fvrtx1, "Error 2");
        }
    }}
    

    // Neighbor vertices
    c_neighborVrtx.resize(g_nCells, g_nFacePerCell, g_nVrtxPerFace);
    for(UINT cell = 0; cell < g_nCells; cell++) {
    for(UINT face = 0; face < g_nFacePerCell; face++) {
    for(UINT fvrtx = 0; fvrtx < g_nVrtxPerFace; fvrtx++) {
        UINT neighborCell = getAdjCell(cell, face);
        
        if(neighborCell == BOUNDARY_FACE) {
            c_neighborVrtx(cell, face, fvrtx) = UINT64_MAX;
            continue;
        }
        
        UINT vrtx = getFaceToCellVrtx(cell, face, fvrtx);
        UINT node = getCellNode(cell, vrtx);
        c_neighborVrtx(cell, face, fvrtx) = getCellVrtx(neighborCell, node);
    }}}
}


/*
    Cell vertex coords
*/
TychoMesh::CellCoords TychoMesh::getCellVrtxCoords(UINT cell) const
{
    CellCoords cellVrtxCoords;
    
    for (UINT vrtx = 0; vrtx < g_nVrtxPerCell; ++vrtx) {
    for (UINT dim = 0; dim < g_ndim; ++dim) {
        cellVrtxCoords.c[vrtx][dim] =
            c_nodeCoords(c_cellNodes(cell, vrtx), dim);
    }}
    
    return cellVrtxCoords;
}


/*
    Face vertex coords
*/
TychoMesh::FaceCoords TychoMesh::getFaceVrtxCoords(UINT cell, UINT face) const 
{
    FaceCoords faceVrtxCoords;

    for (UINT fvrtx = 0; fvrtx < g_nVrtxPerFace; ++fvrtx) {
        UINT cvrtx = getFaceToCellVrtx(cell, face, fvrtx);
        UINT node = getCellNode(cell, cvrtx);
        for (UINT dim = 0; dim < g_ndim; ++dim)
            faceVrtxCoords.c[fvrtx][dim] = c_nodeCoords(node, dim);
    }

    return faceVrtxCoords;
}


/*
    getCellVrtx
*/
UINT TychoMesh::getCellVrtx(const UINT cell, const UINT node) const
{
    for (UINT vrtx = 0; vrtx < g_nVrtxPerCell; ++vrtx) {
        if(getCellNode(cell, vrtx) == node)
            return vrtx;
    }
    
    // Should never reach this point.
    Assert(true);
    return 0;
}







