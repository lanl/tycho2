#if USE_PETSC

//static char help[] = "Solves using Schur Complement.\n\n";

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
#include "SweeperSchur.hh"
#include "SourceIteration.hh"
#include "Global.hh"
#include "TraverseGraph.hh"
#include "Priorities.hh"
#include "Transport.hh"
#include "PsiData.hh"
#include "Comm.hh"
#include "SweepData.hh"
#include "CommSides.hh"
#include <omp.h>
#include <vector>
#include <math.h>
#include <limits>
#include <string.h>
#include <petscmat.h>
#include <petscvec.h>
#include <petscksp.h>


// Initial definition of static variables
static CommSides *s_commSides;
static SweepData *s_sweepData;
static const bool s_doComm = false;
static const UINT s_maxComputePerStep = std::numeric_limits<uint64_t>::max();

static SweeperSchurOuter *s_sweeperSchurOuter;
static PsiData *s_psi;
static PsiData *s_source;


/*
    psiBoundToArray 
*/
static
void psiBoundToArray(PetscScalar Array[], const PsiBoundData &psiBound)
{
    // Start an array index
    int ArrayIndex = 0;
    
    // Incoming flux
    for (UINT angle = 0; angle < g_nAngles; ++angle) {
    for (UINT cell = 0; cell < g_nCells; ++cell) {
    for (UINT group = 0; group < g_nGroups; group++) {
    for (UINT face = 0; face < g_nFacePerCell; face++) {
        if (g_tychoMesh->isIncoming(angle, cell, face)) {
            UINT adjCell = g_tychoMesh->getAdjCell(cell, face);
            UINT adjRank = g_tychoMesh->getAdjRank(cell, face);
            
            // On internal boundary
            if (adjCell == TychoMesh::BOUNDARY_FACE && 
                adjRank != TychoMesh::BAD_RANK)
            {
                UINT side = g_tychoMesh->getSide(cell, face);
                for (UINT fvrtx = 0; fvrtx < g_nVrtxPerFace; fvrtx++) {
                    Array[ArrayIndex] = psiBound(group, fvrtx, angle, side);
                    ArrayIndex++;
                }
            }
        }
    }}}}
}


/*
    arrayToPsiBound
*/
static
void arrayToPsiBound(const PetscScalar Array[], PsiBoundData &psiBound)
{
    // Start an array index
    int ArrayIndex = 0;
    
    // Incoming flux
    for (UINT angle = 0; angle < g_nAngles; ++angle) {
    for (UINT cell = 0; cell < g_nCells; ++cell) {
    for (UINT group = 0; group < g_nGroups; group++) {
    for (UINT face = 0; face < g_nFacePerCell; face++) {
        if (g_tychoMesh->isIncoming(angle, cell, face)) {
            UINT adjCell = g_tychoMesh->getAdjCell(cell, face);
            UINT adjRank = g_tychoMesh->getAdjRank(cell, face);
            
            // On internal boundary
            if (adjCell == TychoMesh::BOUNDARY_FACE && 
                adjRank != TychoMesh::BAD_RANK)
            {
                UINT side = g_tychoMesh->getSide(cell, face);
                for (UINT fvrtx = 0; fvrtx < g_nVrtxPerFace; fvrtx++) {
                    psiBound(group, fvrtx, angle, side) = Array[ArrayIndex];
                    ArrayIndex++;
                }
            }
        }
    }}}}
}


/*
    getVecSize
*/
static
int getVecSize()
{
    // Start an array index
    int size = 0;
    

    // Incoming flux
    for (UINT angle = 0; angle < g_nAngles; ++angle) {
    for (UINT cell = 0; cell < g_nCells; ++cell) {
    for (UINT group = 0; group < g_nGroups; group++) {
    for (UINT face = 0; face < g_nFacePerCell; face++) {
        
        if (g_tychoMesh->isIncoming(angle, cell, face)) {
            
            UINT neighborCell = g_tychoMesh->getAdjCell(cell, face);
            UINT adjRank = g_tychoMesh->getAdjRank(cell, face);
            
            // In local mesh
            if (neighborCell == TychoMesh::BOUNDARY_FACE &&  
                adjRank != TychoMesh::BAD_RANK)
            {
                for (UINT fvrtx = 0; fvrtx < g_nVrtxPerFace; fvrtx++) {
                    size++;
                }
            }
        }
    }}}}

    return size;
}


/*
    Schur 

    This performs the sweep and returns the boundary
*/
static
PetscErrorCode Schur(Mat mat, Vec x, Vec b)
{
    UNUSED_VARIABLE(mat);


    // Vector -> array
    const PetscScalar *In;
    VecGetArrayRead(x, &In);
    arrayToPsiBound(In, s_sweepData->getSideData());
    VecRestoreArrayRead(x, &In);


    // Traverse the graph to get the values on the outward facing boundary
    // call commSides to transfer boundary data
    traverseGraph(s_maxComputePerStep, *s_sweepData, s_doComm, MPI_COMM_WORLD, 
                  Direction_Forward);
    s_commSides->commSides(*s_sweepData);

    
    // Take the values of s_psi from the sweep and put them in b
    PetscScalar *Out;
    VecGetArray(b, &Out);
    psiBoundToArray(Out, s_sweepData->getSideData());
    VecRestoreArray(b, &Out);


    // b = x - b (aka Out = In - Out)
    VecAXPBY(b, 1, -1, x);

    
    return(0);
}


/*
    SchurOuter 

    This performs the sweep and returns the boundary
*/
static
PetscErrorCode SchurOuter(Mat mat, Vec x, Vec b)
{
    UNUSED_VARIABLE(mat);


    // Vector -> array
    const PetscScalar *In;
    VecGetArrayRead(x, &In);
    arrayToPsiBound(In, s_sweepData->getSideData());
    VecRestoreArrayRead(x, &In);


    // Traverse the graph to get the values on the outward facing boundary
    // call commSides to transfer boundary data
    SourceIteration::solve(s_sweeperSchurOuter, *s_psi, *s_source, true);
    s_commSides->commSides(*s_sweepData);

    
    // Take the values of s_psi from the sweep and put them in b
    PetscScalar *Out;
    VecGetArray(b, &Out);
    psiBoundToArray(Out, s_sweepData->getSideData());
    VecRestoreArray(b, &Out);


    // b = x - b (aka Out = In - Out)
    VecAXPBY(b, 1, -1, x);

    
    return(0);
}


/*
    petscInit
*/
void petscInit(Mat &A, Vec &x, Vec &b, KSP &ksp, void (*func)(void))
{
    int argc = 0;
    char **args = NULL;
    PetscInt vecSize;
    PetscInt totalVecSize;


    // Start up petsc
    PetscInitialize(&argc, &args, (char*)0, NULL);

    
    // Local and global vector sizes
    vecSize = getVecSize();
    UINT totalSize = vecSize;
    Comm::gsum(totalSize);
    totalVecSize = totalSize;


    // Create vectors
    VecCreate(MPI_COMM_WORLD, &x);
    VecSetSizes(x, vecSize, totalVecSize);
    VecSetType(x, VECMPI);
    VecDuplicate(x, &b);
    

    // Create matrix shell and define it as the operator
    MatCreateShell(MPI_COMM_WORLD, vecSize, vecSize, totalVecSize, totalVecSize,
                   (void*)(NULL), &A);
    MatShellSetOperation(A, MATOP_MULT, func);
    

    // Set solver context
    KSPCreate(MPI_COMM_WORLD, &ksp);
    

    // Set operator to KSP context. No preconditioning will actually
    // be used due to the PCNONE option.
    KSPSetOperators(ksp, A, A);
    KSPSetTolerances(ksp, g_ddErrMax, PETSC_DEFAULT, PETSC_DEFAULT, 
                     g_ddIterMax);
}


/*
    petscEnd
*/
static
void petscEnd(Mat &A, Vec &x, Vec &b, KSP &ksp)
{
    // Destroy vectors, ksp, and matrices to free work space
    VecDestroy(&x);
    VecDestroy(&b);
    MatDestroy(&A);
    KSPDestroy(&ksp);
    
    
    // Destroy petsc
    PetscFinalize();
}


/*
    solve
*/
void SweeperSchur::solve()
{
    // Init petsc
    Mat A;
    petscInit(A, c_x, c_b, c_ksp, (void(*)(void))Schur);

    
    // Solve
    PsiData psi;
    PsiData totalSource;
    SourceIteration::solve(this, psi, totalSource);

    
    // End petsc
    petscEnd(A, c_x, c_b, c_ksp);
}


/*
    sweep
    
    Run the Krylov solver
*/
void SweeperSchur::sweep(PsiData &psi, const PsiData &source)
{
    // Initialize variables
    PsiData zeroSource;
    zeroSource.setToValue(0.0);
    
    Mat2<UINT> priorities(g_nCells, g_nAngles);
    SweepData zeroSideSweepData(psi, source, g_sigmaTotal, priorities);
    SweepData zeroSourceSweepData(psi, zeroSource, g_sigmaTotal, priorities);
    SweepData sweepData(psi, source, g_sigmaTotal, priorities);
    
    PetscScalar *BIn;
    PetscScalar *XIn;
    PetscScalar *XOut;
    PetscReal rnorm;
    PetscInt its;

    
    // Set static variables
    s_sweepData = &zeroSourceSweepData;
    s_commSides = &c_commSides;
    
    
    // Initial guess (currently zero)
    zeroSideSweepData.zeroSideData();
    VecGetArray(c_x, &XIn);  
    psiBoundToArray(XIn, zeroSideSweepData.getSideData());
    VecRestoreArray(c_x, &XIn);
    
    
    // Do a sweep on the source 
    if (Comm::rank() == 0) {
        printf("    Sweeping Source\n");
    }
    traverseGraph(s_maxComputePerStep, zeroSideSweepData, s_doComm, 
                  MPI_COMM_WORLD, Direction_Forward);
    if (Comm::rank() == 0) {
        printf("    Source Swept\n");
    }


    // Input source into BIn
    c_commSides.commSides(zeroSideSweepData);
    VecGetArray(c_b, &BIn);
    psiBoundToArray(BIn, zeroSideSweepData.getSideData());
    VecRestoreArray(c_b, &BIn);


    // Solve the system (x is the solution, b is the RHS)
    if (Comm::rank() == 0) {
        printf("    Starting Krylov Solve on Boundary\n");
    }
    KSPSolve(c_ksp, c_b, c_x);
    

    // Print some stats
    KSPGetIterationNumber(c_ksp, &its);
    KSPGetResidualNorm(c_ksp, &rnorm);
    if (Comm::rank() == 0) {
        printf("    Krylov iterations: %u with Rnorm: %e\n", its, rnorm);
    }


    // Put x in XOut and output the answer from XOut to psi
    VecGetArray(c_x, &XOut);
    arrayToPsiBound(XOut, sweepData.getSideData());
    VecRestoreArray(c_x, &XOut);


    // Sweep to solve for the non-boundary values
    if (Comm::rank() == 0) {
        printf("    Sweeping to solve non-boundary values\n");
    }
    traverseGraph(s_maxComputePerStep, sweepData, s_doComm, MPI_COMM_WORLD, 
                  Direction_Forward);
    if (Comm::rank() == 0) {
        printf("    Non-boundary values swept\n");
    }
}

void hatSource(const double sigmaT, const double sigmaS, PsiData &source);

/*
    solve
*/
void SweeperSchurOuter::solve()
{
    // Init petsc
    Mat A;
    petscInit(A, c_x, c_b, c_ksp, (void(*)(void))SchurOuter);

    
    // Initialize variables
    Mat2<UINT> priorities(g_nCells, g_nAngles);
    SweepData zeroSideSweepData(c_psi, c_source, g_sigmaTotal, priorities);
    SweepData sweepData(c_psi, c_source, g_sigmaTotal, priorities);
    
    PetscScalar *BIn;
    PetscScalar *XIn;
    PetscScalar *XOut;
    PetscReal rnorm;
    PetscInt its;

    
    // Set static variables
    s_sweepData = &c_sweepData;
    s_commSides = &c_commSides;
    
    s_sweeperSchurOuter = this;
    s_psi = &c_psi;
    s_source = &c_source;
    hatSource(g_sigmaTotal, g_sigmaScat, c_source);
    
    
    // Initial guess (currently zero)
    zeroSideSweepData.zeroSideData();
    VecGetArray(c_x, &XIn);
    psiBoundToArray(XIn, zeroSideSweepData.getSideData());
    VecRestoreArray(c_x, &XIn);
    
    
    // Do a sweep on the source 
    if (Comm::rank() == 0) {
        printf("    Sweeping Source\n");
    }
    traverseGraph(s_maxComputePerStep, zeroSideSweepData, s_doComm, 
                  MPI_COMM_WORLD, Direction_Forward);
    if (Comm::rank() == 0) {
        printf("    Source Swept\n");
    }


    // Input source into BIn
    c_commSides.commSides(zeroSideSweepData);
    VecGetArray(c_b, &BIn);
    psiBoundToArray(BIn, zeroSideSweepData.getSideData());
    VecRestoreArray(c_b, &BIn);


    // Solve the system (x is the solution, b is the RHS)
    if (Comm::rank() == 0) {
        printf("    Starting Krylov Solve on Boundary\n");
    }
    KSPSolve(c_ksp, c_b, c_x);
    

    // Print some stats
    KSPGetIterationNumber(c_ksp, &its);
    KSPGetResidualNorm(c_ksp, &rnorm);
    if (Comm::rank() == 0) {
        printf("    Krylov iterations: %u with Rnorm: %e\n", its, rnorm);
    }


    // Put x in XOut and output the answer from XOut to psi
    VecGetArray(c_x, &XOut);
    arrayToPsiBound(XOut, sweepData.getSideData());
    VecRestoreArray(c_x, &XOut);


    // Sweep to solve for the non-boundary values
    hatSource(g_sigmaTotal, g_sigmaScat, c_source);
    if (Comm::rank() == 0) {
        printf("    Sweeping to solve non-boundary values\n");
    }
    traverseGraph(s_maxComputePerStep, sweepData, s_doComm, MPI_COMM_WORLD, 
                  Direction_Forward);
    if (Comm::rank() == 0) {
        printf("    Non-boundary values swept\n");
    }

    
    // End petsc
    petscEnd(A, c_x, c_b, c_ksp);
}


/*
    sweep
    
    Run the Krylov solver
*/
void SweeperSchurOuter::sweep(PsiData &psi, const PsiData &source)
{
    UNUSED_VARIABLE(psi);
    UNUSED_VARIABLE(source);
    
    // Do 1 graph traversal
    traverseGraph(s_maxComputePerStep, c_sweepData, s_doComm, MPI_COMM_WORLD,
                  Direction_Forward);
}


#endif
