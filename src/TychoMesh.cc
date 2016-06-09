//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mesh/TychoMesh.cc
 * \author Shawn Pautz
 * \date   Fri Jan 14 16:21:46 2000
 * \brief  \link rtt_mesh::TychoMesh TychoMesh \endlink unstructured mesh
 *         class implementation file.
 */
//---------------------------------------------------------------------------//
// $Id: TychoMesh.cc,v 1.6 2000/04/06 21:45:16 pautz Exp $
//---------------------------------------------------------------------------//

#include "TychoMesh.hh"
#include "Global.hh"
#include "Quadrature.hh"

#include <map>
#include <set>
#include <utility>
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
double getTriangleArea(Mat2<double> points)
{
    double a[3], b[3];
    double crossProduct[3];
    
    a[0] = points(1,0) - points(0,0);
    a[1] = points(1,1) - points(0,1);
    a[2] = points(1,2) - points(0,2);
    
    b[0] = points(2,0) - points(0,0);
    b[1] = points(2,1) - points(0,1);
    b[2] = points(2,2) - points(0,2);
    
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
double getTetVolume(Mat2<double> points)
{
    double a[3], b[3], c[3];
    double crossProduct[3];
    
    a[0] = points(1,0) - points(0,0);
    a[1] = points(1,1) - points(0,1);
    a[2] = points(1,2) - points(0,2);
    
    b[0] = points(2,0) - points(0,0);
    b[1] = points(2,1) - points(0,1);
    b[2] = points(2,2) - points(0,2);
    
    c[0] = points(3,0) - points(0,0);
    c[1] = points(3,1) - points(0,1);
    c[2] = points(3,2) - points(0,2);
    
    getCrossProduct(a, b, crossProduct);
    
    return 1.0 / 6.0 * fabs(crossProduct[0] * c[0] + crossProduct[1] * c[1] + 
                            crossProduct[2] * c[2]);
}


/*
    Get the outward normal for a face given face points and cell points.
    The extra cell point gives the outward normal direction.
*/
static
vector<double> getNormal(Mat2<double> facePoints, Mat2<double> cellPoints)
{
    double a[3], b[3];
    double crossProduct[3];
    double norm;
    UINT vertIndex;
    vector<double> normal(3);
    double dotProduct;
    
    
    // Get two vectors to cross
    a[0] = facePoints(1,0) - facePoints(0,0);
    a[1] = facePoints(1,1) - facePoints(0,1);
    a[2] = facePoints(1,2) - facePoints(0,2);
    
    b[0] = facePoints(2,0) - facePoints(0,0);
    b[1] = facePoints(2,1) - facePoints(0,1);
    b[2] = facePoints(2,2) - facePoints(0,2);
    
    
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
            if(cellPoints(vertIndex,0) == facePoints(vert,0) && 
               cellPoints(vertIndex,1) == facePoints(vert,1) && 
               cellPoints(vertIndex,2) == facePoints(vert,2))
            {
                pointEqual = true;
            }
        }
        
        if(pointEqual == false)
            break;
    }
    
    
    // Make sure normal vector is outward normal
    dotProduct = normal[0] * (cellPoints(vertIndex, 0) - facePoints(0,0)) + 
                 normal[1] * (cellPoints(vertIndex, 1) - facePoints(0,1)) + 
                 normal[2] * (cellPoints(vertIndex, 2) - facePoints(0,2));
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
    c_cellVolume = Mat1<double>(c_nCells);
    c_faceArea = Mat2<double>(c_nCells, g_nFacePerCell);
    for(size_t cell = 0; cell < c_nCells; cell++) {
        c_cellVolume(cell) = getTetVolume(getCellVrtxCoords(cell));
        for(UINT face = 0; face < g_nFacePerCell; face++) {
            c_faceArea(cell, face) = getTriangleArea(getFaceVrtxCoords(cell, face));
            Mat2<double> fc = getFaceVrtxCoords(cell, face);
        }
    }
    
    
    // Calculates Omega dot N for every face
    c_omegaDotN = Mat3<double>(g_quadrature->getNumAngles(), c_nCells, g_nFacePerCell);
    for (UINT angle = 0; angle < g_quadrature->getNumAngles(); ++angle) {
        const vector<double> omega = g_quadrature->getOmega(angle);
        for (size_t cell = 0; cell < c_nCells; ++cell) {
        for (UINT face = 0; face < g_nFacePerCell; ++face) {
            vector<double> normal = getNormal(getFaceVrtxCoords(cell, face), 
                                              getCellVrtxCoords(cell));
            for (UINT dim = 0; dim < g_ndim; ++dim)
                c_omegaDotN(angle, cell, face) += omega[dim] * normal[dim];
        }}
    }
    
    
    // CHECK !!!!!!!!
    for(size_t cell = 0; cell < c_nCells; cell++) {
    for(UINT face = 0; face < g_nFacePerCell; face++) {
        for(UINT cvrtx = 0; cvrtx < g_nVrtxPerCell; cvrtx++) {
            if(cvrtx != face) {
                UINT fvrtx = getCellToFaceVrtx(cell, face, cvrtx);
                UINT cvrtx1 = getFaceToCellVrtx(cell, face, fvrtx);
                if(cvrtx != cvrtx1)
                    printf("Error1\n");
            }
        }
        for(UINT fvrtx = 0; fvrtx < g_nVrtxPerFace; fvrtx++) {
            UINT cvrtx = getFaceToCellVrtx(cell, face, fvrtx);
            UINT fvrtx1 = getCellToFaceVrtx(cell, face, cvrtx);
            if(fvrtx != fvrtx1)
                printf("Error2\n");
        }
        }
    }
    
    // Neighbor vertices
    c_neighborVrtx = Mat3<UINT>(c_nCells, g_nFacePerCell, g_nVrtxPerFace);
    for(size_t cell = 0; cell < c_nCells; cell++) {
    for(UINT face = 0; face < g_nFacePerCell; face++) {
    for(UINT fvrtx = 0; fvrtx < g_nVrtxPerFace; fvrtx++) {
        size_t neighborCell = getAdjCell(cell, face);
        
        if(neighborCell == BOUNDARY_FACE) {
            c_neighborVrtx(cell, face, fvrtx) = -1;
            continue;
        }
        
        UINT vrtx = getFaceToCellVrtx(cell, face, fvrtx);
        size_t node = getCellNode(cell, vrtx);
        c_neighborVrtx(cell, face, fvrtx) = getCellVrtx(neighborCell, node);
    }}}
}


/*
    Cell vertex coords
*/
Mat2<double> TychoMesh::getCellVrtxCoords(UINT cell) const
{
    Mat2<double> cellVrtxCoords(g_nVrtxPerCell, g_ndim);
    
    for (UINT vrtx = 0; vrtx < g_nVrtxPerCell; ++vrtx) {
    for (UINT dim = 0; dim < g_ndim; ++dim) {
        cellVrtxCoords(vrtx, dim) =
            c_nodeCoords(c_cellNodes(cell, vrtx), dim);
    }}
    
    return cellVrtxCoords;
}


/*
    Face vertex coords
*/
Mat2<double> TychoMesh::getFaceVrtxCoords(UINT cell, UINT face) const 
{
    Mat2<double> faceVrtxCoords(g_nVrtxPerFace, g_ndim);

    for (UINT fvrtx = 0; fvrtx < g_nVrtxPerFace; ++fvrtx) {
        UINT cvrtx = getFaceToCellVrtx(cell, face, fvrtx);
        size_t node = getCellNode(cell, cvrtx);
        for (UINT dim = 0; dim < g_ndim; ++dim)
            faceVrtxCoords(fvrtx, dim) = c_nodeCoords(node, dim);
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







