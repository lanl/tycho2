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

#include "TychoMesh.hh"
#include "Global.hh"
#include "Quadrature.hh"
#include <stdio.h>
#include <math.h>

using namespace std;


/*
    Returns c = a x b
*/
static 
void getCrossProduct(double a[3], double b[3], double c[3])
{
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}


/*
    Area of a triangle is 1/2 | A x B | with 
    A = points1 - points0 and B = points2 - points0.
    
    points is indexed by points(index, dim)
*/
static 
double getTriangleArea(TychoMesh::FaceCoords points)
{
    double a[3], b[3];
    double crossProduct[3];
    
    a[0] = points.c[1][0] - points.c[0][0];
    a[1] = points.c[1][1] - points.c[0][1];
    a[2] = points.c[1][2] - points.c[0][2];
    
    b[0] = points.c[2][0] - points.c[0][0];
    b[1] = points.c[2][1] - points.c[0][1];
    b[2] = points.c[2][2] - points.c[0][2];
    
    getCrossProduct(a, b, crossProduct);
    
    return 0.5 * sqrt(crossProduct[0] * crossProduct[0] + 
                      crossProduct[1] * crossProduct[1] + 
                      crossProduct[2] * crossProduct[2]);
}


/*
    Volume of a tetrahedron is 1/6 |C . (A x B)| with
    A = points1-points0, B = points2-points0, C=points3-points0
    
    points is indexed by points(index, dim)
*/
static 
double getTetVolume(TychoMesh::CellCoords points)
{
    double a[3], b[3], c[3];
    double crossProduct[3];
    
    a[0] = points.c[1][0] - points.c[0][0];
    a[1] = points.c[1][1] - points.c[0][1];
    a[2] = points.c[1][2] - points.c[0][2];
    
    b[0] = points.c[2][0] - points.c[0][0];
    b[1] = points.c[2][1] - points.c[0][1];
    b[2] = points.c[2][2] - points.c[0][2];
    
    c[0] = points.c[3][0] - points.c[0][0];
    c[1] = points.c[3][1] - points.c[0][1];
    c[2] = points.c[3][2] - points.c[0][2];
    
    getCrossProduct(a, b, crossProduct);
    
    return 1.0 / 6.0 * fabs(crossProduct[0] * c[0] + crossProduct[1] * c[1] + 
                            crossProduct[2] * c[2]);
}


/*
    Get the outward normal for a face given face points and cell points.
    The extra cell point gives the outward normal direction.
*/
static
vector<double> getNormal(TychoMesh::FaceCoords facePoints, 
                         TychoMesh::CellCoords cellPoints)
{
    double a[3], b[3];
    double crossProduct[3];
    double norm;
    UINT vertIndex;
    vector<double> normal(3);
    double dotProduct;
    
    
    // Get two vectors to cross
    a[0] = facePoints.c[1][0] - facePoints.c[0][0];
    a[1] = facePoints.c[1][1] - facePoints.c[0][1];
    a[2] = facePoints.c[1][2] - facePoints.c[0][2];
    
    b[0] = facePoints.c[2][0] - facePoints.c[0][0];
    b[1] = facePoints.c[2][1] - facePoints.c[0][1];
    b[2] = facePoints.c[2][2] - facePoints.c[0][2];
    
    
    // Get a normal vector
    getCrossProduct(a, b, crossProduct);
    
    norm = sqrt(crossProduct[0] * crossProduct[0] + 
                crossProduct[1] * crossProduct[1] + 
                crossProduct[2] * crossProduct[2]);
    
    normal[0] = crossProduct[0] / norm;
    normal[1] = crossProduct[1] / norm;
    normal[2] = crossProduct[2] / norm;
    
    
    // Find index of point in cellPoints not in facePoints
    for(vertIndex = 0; vertIndex < 4; vertIndex++) {
        bool pointEqual = false;
        for(UINT vert = 0; vert < 3; vert++) {
            if(cellPoints.c[vertIndex][0] == facePoints.c[vert][0] && 
               cellPoints.c[vertIndex][1] == facePoints.c[vert][1] && 
               cellPoints.c[vertIndex][2] == facePoints.c[vert][2])
            {
                pointEqual = true;
            }
        }
        
        if(pointEqual == false)
            break;
    }
    
    
    // Make sure normal vector is outward normal
    dotProduct = normal[0] * (cellPoints.c[vertIndex][0] - facePoints.c[0][0]) + 
                 normal[1] * (cellPoints.c[vertIndex][1] - facePoints.c[0][1]) + 
                 normal[2] * (cellPoints.c[vertIndex][2] - facePoints.c[0][2]);
    if(dotProduct > 0) {
        normal[0] = -normal[0];
        normal[1] = -normal[1];
        normal[2] = -normal[2];
    }
    
    
    return normal;
}


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







