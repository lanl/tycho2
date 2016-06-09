//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   transport/Transport.cc
 * \author Shawn Pautz
 * \date   Fri Jan 14 16:30:59 2000
 * \brief \link rtt_transport::Transport Transport \endlink discretization
 *        method class implementation file.
 */
//---------------------------------------------------------------------------//
// $Id: Transport.cc,v 1.13 2000/10/17 16:36:45 pautz Exp $
//---------------------------------------------------------------------------//

#include "Transport.hh"
#include "Assert.hh"
#include "Global.hh"
#include "TychoMesh.hh"
#include "PsiData.hh"
#include <iostream>
#include <string>
#include <cmath>
#include <stdio.h>

// Namespaces
using namespace std;


/*
    calcLocalSource
*/
static
void calcLocalSource(const UINT cell,
                     const Mat2<double> &localSource,
                     double cellSource[g_nVrtxPerCell],
                     const UINT group) 
{
    double volume = g_spTychoMesh->getCellVolume(cell);
    
    double q0 = localSource(0, group);
    double q1 = localSource(1, group);
    double q2 = localSource(2, group);
    double q3 = localSource(3, group);
    
    cellSource[0] = volume / 20.0 * (2.0 * q0 + q1 + q2 + q3);
    cellSource[1] = volume / 20.0 * (q0 + 2.0 * q1 + q2 + q3);
    cellSource[2] = volume / 20.0 * (q0 + q1 + 2.0 * q2 + q3);
    cellSource[3] = volume / 20.0 * (q0 + q1 + q2 + 2.0 * q3);
}


/*
    streamPlusColl
*/
static
void streamPlusColl(const UINT cell, const UINT angle,
                    const double sigmaTotal,
                    double matrix[g_nVrtxPerCell][g_nVrtxPerCell]) 
{
    double volume, area[4];

    // Get cell volume and face areas
    volume = g_spTychoMesh->getCellVolume(cell);
    for(UINT face = 0; face < 4; face++) {
        area[face] = g_spTychoMesh->getFaceArea(cell, face) * 
                     g_spTychoMesh->getOmegaDotN(angle, cell, face);
    }
    
    // Setup LHS volume integrals
    for(UINT row = 0; row < 4; row++) {
    for(UINT col = 0; col < 4; col++) {
        double factor = (row == col) ? 2.0 : 1.0;
        matrix[row][col] = area[row] / 12.0 + sigmaTotal * volume / 20.0 * factor;
    }}
}


/*
    faceDependence
    Takes a lot of computational time, so it has been optimized.
    Beware of changing order of operations.
*/
static
void faceDependence(const UINT cell, const UINT angle,
                    double matrix[g_nVrtxPerCell][g_nVrtxPerCell],
                    const Mat3<double> &localPsiBound, 
                    double cellSource[g_nVrtxPerCell], const UINT group) 
{
    double area[4];
    UINT cellToFaceVrtx[4][4];
    UINT indices[4][3];
    
    
    // Get face areas
    for(UINT face = 0; face < 4; face++) {
        area[face] = g_spTychoMesh->getFaceArea(cell, face) * 
                     g_spTychoMesh->getOmegaDotN(angle, cell, face);
    }
    
    
    // Populate cellToFaceVrtx
    for(UINT face = 0; face < 4; face++) {
    for(UINT col = 0; col < 4; col++) {
        if(face != col)
            cellToFaceVrtx[face][col] = g_spTychoMesh->getCellToFaceVrtx(cell, face, col);
    }}
    
    
    // Populate indices
    indices[0][0] = 1; indices[0][1] = 2; indices[0][2] = 3;
    indices[1][0] = 0; indices[1][1] = 2; indices[1][2] = 3;
    indices[2][0] = 0; indices[2][1] = 1; indices[2][2] = 3;
    indices[3][0] = 0; indices[3][1] = 1; indices[3][2] = 2;
    
    
    // Update for fluxes.
    for(UINT face = 0; face < 4; face++) {
        
        // Outgoing flux
        if(area[face] > 0) {
            for(UINT rowIndex = 0; rowIndex < 3; rowIndex++) {
            for(UINT colIndex = 0; colIndex < 3; colIndex++) {
                UINT row = indices[face][rowIndex];
                UINT col = indices[face][colIndex];
                double factor = (row == col) ? 2.0 : 1.0;
                matrix[row][col] += area[face] / 12.0 * factor;
            }}
        }
        
        // Incoming flux
        else {
            for(UINT rowIndex = 0; rowIndex < 3; rowIndex++) {
            for(UINT colIndex = 0; colIndex < 3; colIndex++) {
                UINT row = indices[face][rowIndex];
                UINT col = indices[face][colIndex];
                double factor = (row == col) ? 2.0 : 1.0;
                UINT faceVertex = cellToFaceVrtx[face][col];
                double psiNeighbor = localPsiBound(faceVertex, face, group);
                cellSource[row] -= area[face] / 12.0 * psiNeighbor * factor;
            }}
        }
    }
}


/*
    gaussElim4
*/
static 
void gaussElim4(double A[4][4], double b[4])
{
    const int n = 4;
    
    for (int column = 0; column < n-1; ++column) {
        int rowmax = column;
        double colmax = fabs(A[column][column]);
        for (int row = column+1; row < n; ++row) {
            double temp = fabs(A[row][column]);
            if (temp > colmax)  {
                rowmax = row;
                colmax = temp;
            }
        }

        if (rowmax != column) {
            for (int column2 = 0; column2 < n; ++column2) {
                double temp = A[rowmax][column2];
                A[rowmax][column2] = A[column][column2];
                A[column][column2] = temp;
            }
            double temp = b[rowmax];
            b[rowmax] = b[column];
            b[column] = temp;
        }

        Assert(A[column][column] != 0.);
        A[column][column] = 1./A[column][column];

        for (int row = column+1; row < n; ++row)
            A[row][column] *= A[column][column];

        for (int column2 = 0; column2 <= column; ++column2) {
        for (int row = column2+1; row < n; ++row) {
            A[row][column+1] -= A[row][column2]*A[column2][column+1];
        }}
    }

    Assert(A[n-1][n-1] != 0.);
    A[n-1][n-1] = 1./A[n-1][n-1];

    for (int column = 0; column < n-1; ++column) {
    for (int row = column+1; row < n; ++row) {
        b[row] -= A[row][column]*b[column];
    }}

    for (int column = n-1; column >= 0; --column) {
        b[column] *= A[column][column];
        for (int row = column-1; row >= 0; --row)
            b[row] -= A[row][column]*b[column];
    }
}


// Global functions
namespace Transport
{

/*
    solve
*/
void solve(const UINT cell, const UINT angle, const double sigmaTotal,
           const Mat3<double> &localPsiBound, const Mat2<double> &localSource,
           Mat2<double> &localPsi)
{
    for(UINT group = 0; group < g_nGroups; group++) {
        double matrix[g_nVrtxPerCell][g_nVrtxPerCell] = {0.0};
        double cellSource[g_nVrtxPerCell] = {0.0};
        double solution[g_nVrtxPerCell];
        
        // form local source term
        calcLocalSource(cell, localSource, cellSource, group);
        
        // form streaming-plus-collision portion of matrix
        streamPlusColl(cell, angle, sigmaTotal, matrix);
        
        // form dependencies on incoming (outgoing) faces
        faceDependence(cell, angle, matrix, localPsiBound, cellSource, group);
        
        // solve matrix
        for (UINT vertex = 0; vertex < g_nVrtxPerCell; ++vertex)
            solution[vertex] = cellSource[vertex];
        gaussElim4(matrix, solution);
        
        // put local solution onto global solution
        for (UINT vertex = 0; vertex < g_nVrtxPerCell; ++vertex)
            localPsi(vertex, group) = solution[vertex];
    }
}


} // End namespace Transport


