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
    calcSource
*/
static
KOKKOS_INLINE_FUNCTION
void calcSourceKokkos(
    const double volume, 
    const double localSource[g_nVrtxPerCell][g_nMaxGroups],
    double cellSource[g_nVrtxPerCell],
    const UINT group) 
{
    double q0 = localSource[0][group];
    double q1 = localSource[1][group];
    double q2 = localSource[2][group];
    double q3 = localSource[3][group];
    
    cellSource[0] = volume / 20.0 * (2.0 * q0 + q1 + q2 + q3);
    cellSource[1] = volume / 20.0 * (q0 + 2.0 * q1 + q2 + q3);
    cellSource[2] = volume / 20.0 * (q0 + q1 + 2.0 * q2 + q3);
    cellSource[3] = volume / 20.0 * (q0 + q1 + q2 + 2.0 * q3);
}


/*
    calcVolumeIntegrals
*/
static
KOKKOS_INLINE_FUNCTION
void calcVolumeIntegralsKokkos(
    const double volume, 
    const double area[g_nFacePerCell],
    const double sigmaTotal,
    double matrix[g_nVrtxPerCell][g_nVrtxPerCell]) 
{
    matrix[0][0] = area[0] / 12.0 + 2.0 * sigmaTotal * volume / 20.0;
    matrix[0][1] = area[0] / 12.0 + 1.0 * sigmaTotal * volume / 20.0;
    matrix[0][2] = area[0] / 12.0 + 1.0 * sigmaTotal * volume / 20.0;
    matrix[0][3] = area[0] / 12.0 + 1.0 * sigmaTotal * volume / 20.0;
    
    matrix[1][0] = area[1] / 12.0 + 1.0 * sigmaTotal * volume / 20.0;
    matrix[1][1] = area[1] / 12.0 + 2.0 * sigmaTotal * volume / 20.0;
    matrix[1][2] = area[1] / 12.0 + 1.0 * sigmaTotal * volume / 20.0;
    matrix[1][3] = area[1] / 12.0 + 1.0 * sigmaTotal * volume / 20.0;
    
    matrix[2][0] = area[2] / 12.0 + 1.0 * sigmaTotal * volume / 20.0;
    matrix[2][1] = area[2] / 12.0 + 1.0 * sigmaTotal * volume / 20.0;
    matrix[2][2] = area[2] / 12.0 + 2.0 * sigmaTotal * volume / 20.0;
    matrix[2][3] = area[2] / 12.0 + 1.0 * sigmaTotal * volume / 20.0;
    
    matrix[3][0] = area[3] / 12.0 + 1.0 * sigmaTotal * volume / 20.0;
    matrix[3][1] = area[3] / 12.0 + 1.0 * sigmaTotal * volume / 20.0;
    matrix[3][2] = area[3] / 12.0 + 1.0 * sigmaTotal * volume / 20.0;
    matrix[3][3] = area[3] / 12.0 + 2.0 * sigmaTotal * volume / 20.0;
}


/*
    calcOutgoingFlux
*/
static
KOKKOS_INLINE_FUNCTION
void calcOutgoingFluxKokkos(
    const double area[g_nFacePerCell],
    double matrix[g_nVrtxPerCell][g_nVrtxPerCell])
{
    if (area[0] > 0) {
        matrix[1][1] += 2.0 * area[0] / 12.0;
        matrix[1][2] += 1.0 * area[0] / 12.0;
        matrix[1][3] += 1.0 * area[0] / 12.0;

        matrix[2][1] += 1.0 * area[0] / 12.0;
        matrix[2][2] += 2.0 * area[0] / 12.0;
        matrix[2][3] += 1.0 * area[0] / 12.0;

        matrix[3][1] += 1.0 * area[0] / 12.0;
        matrix[3][2] += 1.0 * area[0] / 12.0;
        matrix[3][3] += 2.0 * area[0] / 12.0;
    }
    
    if (area[1] > 0) {
        matrix[0][0] += 2.0 * area[1] / 12.0;
        matrix[0][2] += 1.0 * area[1] / 12.0;
        matrix[0][3] += 1.0 * area[1] / 12.0;

        matrix[2][0] += 1.0 * area[1] / 12.0;
        matrix[2][2] += 2.0 * area[1] / 12.0;
        matrix[2][3] += 1.0 * area[1] / 12.0;

        matrix[3][0] += 1.0 * area[1] / 12.0;
        matrix[3][2] += 1.0 * area[1] / 12.0;
        matrix[3][3] += 2.0 * area[1] / 12.0;
    }
    
    if (area[2] > 0) {
        matrix[0][0] += 2.0 * area[2] / 12.0;
        matrix[0][1] += 1.0 * area[2] / 12.0;
        matrix[0][3] += 1.0 * area[2] / 12.0;

        matrix[1][0] += 1.0 * area[2] / 12.0;
        matrix[1][1] += 2.0 * area[2] / 12.0;
        matrix[1][3] += 1.0 * area[2] / 12.0;

        matrix[3][0] += 1.0 * area[2] / 12.0;
        matrix[3][1] += 1.0 * area[2] / 12.0;
        matrix[3][3] += 2.0 * area[2] / 12.0;
    }
    
    if (area[3] > 0) {
        matrix[0][0] += 2.0 * area[3] / 12.0;
        matrix[0][1] += 1.0 * area[3] / 12.0;
        matrix[0][2] += 1.0 * area[3] / 12.0;

        matrix[1][0] += 1.0 * area[3] / 12.0;
        matrix[1][1] += 2.0 * area[3] / 12.0;
        matrix[1][2] += 1.0 * area[3] / 12.0;

        matrix[2][0] += 1.0 * area[3] / 12.0;
        matrix[2][1] += 1.0 * area[3] / 12.0;
        matrix[2][2] += 2.0 * area[3] / 12.0;
    }
}


/*
    calcIncomingFlux
*/
static
KOKKOS_INLINE_FUNCTION
void calcIncomingFluxKokkos(
    const UINT cell, 
    const double area[g_nFacePerCell],
    const double localPsiBound[g_nVrtxPerFace][g_nFacePerCell][g_nMaxGroups],
    double cellSource[g_nVrtxPerCell],
    const UINT group,
    device_mat3_t<UINT> const &cell_to_face_vertex)
{
    UINT faceVertex0, faceVertex1, faceVertex2, faceVertex3;
    double psiNeighbor0, psiNeighbor1, psiNeighbor2, psiNeighbor3;

    if (area[0] < 0) {
        faceVertex1 = cell_to_face_vertex(cell, 0, 1);
        faceVertex2 = cell_to_face_vertex(cell, 0, 2);
        faceVertex3 = cell_to_face_vertex(cell, 0, 3);
        
        psiNeighbor1 = localPsiBound[faceVertex1][0][group];
        psiNeighbor2 = localPsiBound[faceVertex2][0][group];
        psiNeighbor3 = localPsiBound[faceVertex3][0][group];
        
        cellSource[1] -= 2.0 * area[0] / 12.0 * psiNeighbor1;
        cellSource[1] -= 1.0 * area[0] / 12.0 * psiNeighbor2;
        cellSource[1] -= 1.0 * area[0] / 12.0 * psiNeighbor3;
        
        cellSource[2] -= 1.0 * area[0] / 12.0 * psiNeighbor1;
        cellSource[2] -= 2.0 * area[0] / 12.0 * psiNeighbor2;
        cellSource[2] -= 1.0 * area[0] / 12.0 * psiNeighbor3;
        
        cellSource[3] -= 1.0 * area[0] / 12.0 * psiNeighbor1;
        cellSource[3] -= 1.0 * area[0] / 12.0 * psiNeighbor2;
        cellSource[3] -= 2.0 * area[0] / 12.0 * psiNeighbor3;
    }

    if (area[1] < 0) {
        faceVertex0 = cell_to_face_vertex(cell, 1, 0);
        faceVertex2 = cell_to_face_vertex(cell, 1, 2);
        faceVertex3 = cell_to_face_vertex(cell, 1, 3);
        
        psiNeighbor0 = localPsiBound[faceVertex0][1][group];
        psiNeighbor2 = localPsiBound[faceVertex2][1][group];
        psiNeighbor3 = localPsiBound[faceVertex3][1][group];
        
        cellSource[0] -= 2.0 * area[1] / 12.0 * psiNeighbor0;
        cellSource[0] -= 1.0 * area[1] / 12.0 * psiNeighbor2;
        cellSource[0] -= 1.0 * area[1] / 12.0 * psiNeighbor3;
        
        cellSource[2] -= 1.0 * area[1] / 12.0 * psiNeighbor0;
        cellSource[2] -= 2.0 * area[1] / 12.0 * psiNeighbor2;
        cellSource[2] -= 1.0 * area[1] / 12.0 * psiNeighbor3;
        
        cellSource[3] -= 1.0 * area[1] / 12.0 * psiNeighbor0;
        cellSource[3] -= 1.0 * area[1] / 12.0 * psiNeighbor2;
        cellSource[3] -= 2.0 * area[1] / 12.0 * psiNeighbor3;
    }

    if (area[2] < 0) {
        faceVertex0 = cell_to_face_vertex(cell, 2, 0);
        faceVertex1 = cell_to_face_vertex(cell, 2, 1);
        faceVertex3 = cell_to_face_vertex(cell, 2, 3);
        
        psiNeighbor0 = localPsiBound[faceVertex0][2][group];
        psiNeighbor1 = localPsiBound[faceVertex1][2][group];
        psiNeighbor3 = localPsiBound[faceVertex3][2][group];
        
        cellSource[0] -= 2.0 * area[2] / 12.0 * psiNeighbor0;
        cellSource[0] -= 1.0 * area[2] / 12.0 * psiNeighbor1;
        cellSource[0] -= 1.0 * area[2] / 12.0 * psiNeighbor3;
        
        cellSource[1] -= 1.0 * area[2] / 12.0 * psiNeighbor0;
        cellSource[1] -= 2.0 * area[2] / 12.0 * psiNeighbor1;
        cellSource[1] -= 1.0 * area[2] / 12.0 * psiNeighbor3;
        
        cellSource[3] -= 1.0 * area[2] / 12.0 * psiNeighbor0;
        cellSource[3] -= 1.0 * area[2] / 12.0 * psiNeighbor1;
        cellSource[3] -= 2.0 * area[2] / 12.0 * psiNeighbor3;
    }

    if (area[3] < 0) {
        faceVertex0 = cell_to_face_vertex(cell, 3, 0);
        faceVertex1 = cell_to_face_vertex(cell, 3, 1);
        faceVertex2 = cell_to_face_vertex(cell, 3, 2);
        
        psiNeighbor0 = localPsiBound[faceVertex0][3][group];
        psiNeighbor1 = localPsiBound[faceVertex1][3][group];
        psiNeighbor2 = localPsiBound[faceVertex2][3][group];
        
        cellSource[0] -= 2.0 * area[3] / 12.0 * psiNeighbor0;
        cellSource[0] -= 1.0 * area[3] / 12.0 * psiNeighbor1;
        cellSource[0] -= 1.0 * area[3] / 12.0 * psiNeighbor2;
        
        cellSource[1] -= 1.0 * area[3] / 12.0 * psiNeighbor0;
        cellSource[1] -= 2.0 * area[3] / 12.0 * psiNeighbor1;
        cellSource[1] -= 1.0 * area[3] / 12.0 * psiNeighbor2;
        
        cellSource[2] -= 1.0 * area[3] / 12.0 * psiNeighbor0;
        cellSource[2] -= 1.0 * area[3] / 12.0 * psiNeighbor1;
        cellSource[2] -= 2.0 * area[3] / 12.0 * psiNeighbor2;
    }
}


/*
    gaussElim4
*/
static 
KOKKOS_INLINE_FUNCTION
void gaussElim4Kokkos(double A[4][4], double b[4])
{
    double tmp;
    
    // Normalize first row
    tmp = 1.0/A[0][0];
    A[0][0] = 1.0;
    A[0][1] = A[0][1] * tmp;
    A[0][2] = A[0][2] * tmp;
    A[0][3] = A[0][3] * tmp;
    b[0] = b[0] * tmp;

    // Set column zero to 0.0
    tmp = A[1][0];
    A[1][0] = 0.0;
    A[1][1] = A[1][1] - A[0][1] * tmp;
    A[1][2] = A[1][2] - A[0][2] * tmp;
    A[1][3] = A[1][3] - A[0][3] * tmp;
    b[1] = b[1] - b[0] * tmp;   

    tmp = A[2][0];
    A[2][0] = 0.0;
    A[2][1] = A[2][1] - A[0][1] * tmp;
    A[2][2] = A[2][2] - A[0][2] * tmp;
    A[2][3] = A[2][3] - A[0][3] * tmp;
    b[2] = b[2] - b[0] * tmp;
    
    tmp = A[3][0];
    A[3][0] = 0.0;
    A[3][1] = A[3][1] - A[0][1] * tmp;
    A[3][2] = A[3][2] - A[0][2] * tmp;
    A[3][3] = A[3][3] - A[0][3] * tmp;
    b[3] = b[3] - b[0] * tmp;
    
    // Normalize second row
    tmp = 1.0/A[1][1];
    A[1][1] = 1.0;
    A[1][2] = A[1][2] * tmp;
    A[1][3] = A[1][3] * tmp;
    b[1] = b[1] * tmp;
    
    // Set column one to 0.0
    tmp = A[2][1];
    A[2][1] = 0.0;
    A[2][2] = A[2][2] - A[1][2] * tmp;
    A[2][3] = A[2][3] - A[1][3] * tmp;
    b[2] = b[2] - b[1] * tmp;
    
    tmp = A[3][1];
    A[3][1] = 0.0;
    A[3][2] = A[3][2] - A[1][2] * tmp;
    A[3][3] = A[3][3] - A[1][3] * tmp;
    b[3] = b[3] - b[1] * tmp;
    
    // Normalize third row
    tmp = 1.0/A[2][2];
    A[2][2] = 1.0;
    A[2][3] = A[2][3] * tmp;
    b[2] = b[2] * tmp;
    
    // Set column two to 0.0
    tmp = A[3][2];
    A[3][2] = 0.0;
    A[3][3] = A[3][3] - A[2][3] * tmp;
    b[3] = b[3] - b[2] * tmp;

    // Backward Solve
    b[3] = b[3]/A[3][3];    
    b[2] = b[2] - A[2][3]*b[3];
    b[1] = b[1] - A[1][3]*b[3] - A[1][2]*b[2];
    b[0] = b[0] - A[0][3]*b[3] - A[0][2]*b[2] - A[0][1]*b[1]; 
} 


/*
    solve
*/
static 
KOKKOS_INLINE_FUNCTION
void solveKokkos(
    const UINT cell, 
    const UINT angle, 
    const double sigmaTotal,
    const double localPsiBound[g_nVrtxPerFace][g_nFacePerCell][g_nMaxGroups],
    const double localSource[g_nVrtxPerCell][g_nMaxGroups],
    double localPsi[g_nVrtxPerCell][g_nMaxGroups],
    device_mat1_t<double> const &cell_volume,
    device_mat2_t<double> const &face_area,
    device_mat3_t<double> const &omega_dot_n,
    device_mat3_t<UINT> const &cell_to_face_vertex,
    UINT nGroups)
{
    double volume, area[g_nFacePerCell];

    
    // Get cell volume and face areas
    volume = cell_volume(cell);
    
    area[0] = face_area(cell, 0) * 
              omega_dot_n(angle, cell, 0);
    area[1] = face_area(cell, 1) * 
              omega_dot_n(angle, cell, 1);
    area[2] = face_area(cell, 2) * 
              omega_dot_n(angle, cell, 2);
    area[3] = face_area(cell, 3) * 
              omega_dot_n(angle, cell, 3);
    
    
    // Solve local transport problem for each group
    for (UINT group = 0; group < nGroups; group++) {
        
        double cellSource[g_nVrtxPerCell] = {0.0};
        double matrix[g_nVrtxPerCell][g_nVrtxPerCell] = {0.0};
        double solution[g_nVrtxPerCell];
    
        // form local source term
        calcSourceKokkos(volume, localSource, cellSource, group);
        
        // form streaming-plus-collision portion of matrix
        calcVolumeIntegralsKokkos(volume, area, sigmaTotal, matrix);
        
        // form dependencies on incoming (outgoing) faces
        calcOutgoingFluxKokkos(area, matrix);
        calcIncomingFluxKokkos(cell, area, localPsiBound, cellSource, group,
                               cell_to_face_vertex);
        
        // solve matrix
        for (UINT vertex = 0; vertex < g_nVrtxPerCell; ++vertex)
            solution[vertex] = cellSource[vertex];
        gaussElim4Kokkos(matrix, solution);
        
        // put local solution onto global solution
        for (UINT vertex = 0; vertex < g_nVrtxPerCell; ++vertex)
            localPsi[vertex][group] = solution[vertex];
    }
}

/*
    populateLocalPsiBound
    
    Put data from neighboring cells into localPsiBound(fvrtx, face, group).
*/
static 
KOKKOS_INLINE_FUNCTION
void populateLocalPsiBoundKokkos(
    const UINT angle, 
    const UINT cell, 
    device_psi_data_t const &psi, 
    device_psi_data_t const &psiBound, 
    double localPsiBound[g_nVrtxPerFace][g_nFacePerCell][g_nMaxGroups],
    device_mat3_t<double> const& omega_dot_n,
    device_mat2_t<UINT> const& adj_cell,
    device_mat3_t<UINT> const& neighbor_vertex,
    device_mat2_t<UINT> const& adj_proc,
    device_mat2_t<UINT> const& sides,
    UINT nGroups)
{
    for (UINT i = 0; i < g_nVrtxPerFace; i++) {
    for (UINT j = 0; j < g_nFacePerCell; j++) {
    for (UINT k = 0; k < nGroups; k++) {
        localPsiBound[i][j][k] = 0.0;
    }}}
    
    // Populate if incoming flux
    #pragma omp simd
    for (UINT group = 0; group < nGroups; group++) {
    for (UINT face = 0; face < g_nFacePerCell; face++) {
        if (omega_dot_n(angle, cell, face) <= 0) {
            UINT neighborCell = adj_cell(cell, face);
            
            // In local mesh
            if (neighborCell != TychoMesh::BOUNDARY_FACE) {
                for (UINT fvrtx = 0; fvrtx < g_nVrtxPerFace; fvrtx++) {
                    UINT neighborVrtx = neighbor_vertex(cell, face, fvrtx);
                    localPsiBound[fvrtx][face][group] = 
                        psi(group, neighborVrtx, angle, neighborCell);
                }
            }
            
            // Not in local mesh
            else if (adj_proc(cell, face) != TychoMesh::BAD_RANK) {
                for (UINT fvrtx = 0; fvrtx < g_nVrtxPerFace; fvrtx++) {
                    UINT side = sides(cell, face);
                    localPsiBound[fvrtx][face][group] = 
                        psiBound(group, fvrtx, angle, side);
                }
            }
        }
    }}
}


// Global functions
namespace Transport
{

/*
    update
    
    Does a transport update for the given cell/angle pair.
*/
KOKKOS_INLINE_FUNCTION
void updateKokkos(
    const UINT cell, 
    const UINT angle,
    device_psi_data_t const &device_source, 
    device_psi_data_t const &device_psi_bound, 
    device_psi_data_t device_psi,       // WHY NO REFERENCE HERE???
    UINT nGroups,
    device_mat3_t<double> const& device_omega_dot_n,
    device_mat2_t<UINT> const& device_adj_cell,
    device_mat3_t<UINT> const& device_neighbor_vertex,
    device_mat2_t<UINT> const& device_adj_proc,
    device_mat2_t<UINT> const& device_side,
    device_mat1_t<double> const &device_cell_volume,
    device_mat1_t<double> const &device_sigma_t,
    device_mat2_t<double> const &device_face_area,
    device_mat3_t<UINT> const &device_cell_to_face_vertex)
{
    double localSource[g_nVrtxPerCell][g_nMaxGroups];
    double localPsi[g_nVrtxPerCell][g_nMaxGroups];
    double localPsiBound[g_nVrtxPerFace][g_nFacePerCell][g_nMaxGroups];

    
    // Populate localSource
    for (UINT group = 0; group < nGroups; group++) {
    for (UINT vrtx = 0; vrtx < g_nVrtxPerCell; vrtx++) {
        localSource[vrtx][group] = device_source(group, vrtx, angle, cell);
    }}
    
    
    // Populate localPsiBound
    populateLocalPsiBoundKokkos(angle, cell, device_psi, device_psi_bound, 
                                localPsiBound, device_omega_dot_n,
                                device_adj_cell, device_neighbor_vertex, 
                                device_adj_proc, device_side, nGroups);
    
    
    // Transport solve
    solveKokkos(cell, angle, device_sigma_t(cell), localPsiBound, localSource, 
          localPsi, device_cell_volume, device_face_area, device_omega_dot_n,
          device_cell_to_face_vertex, nGroups);
    
    
    // localPsi -> psi
    for (UINT group = 0; group < nGroups; group++) {
    for (UINT vrtx = 0; vrtx < g_nVrtxPerCell; vrtx++) {
        device_psi(group, vrtx, angle, cell) = localPsi[vrtx][group];
    }}
}



} // End namespace Transport


