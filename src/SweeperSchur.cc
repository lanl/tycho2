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


//Initial definition of static variables
static CommSides *s_commSides;
static SweepData *s_sweepData;
static const bool doComm = false;
static const UINT maxComputePerStep = std::numeric_limits<uint64_t>::max();


/*
    psiBoundToArray 
*/
void psiBoundToArray(PetscScalar Array[], const PsiBoundData &psiBound)
{
    //Start an array index
    int ArrayIndex = 0;
     
    // Incoming flux
    for (UINT angle = 0; angle < g_nAngles; ++angle) {
    for (UINT cell = 0; cell < g_nCells; ++cell) {
    for (UINT group = 0; group < g_nGroups; group++) {
    for (UINT face = 0; face < g_nFacePerCell; face++) {
        if (g_tychoMesh->isIncoming(angle, cell, face)) {
            UINT adjCell = g_tychoMesh->getAdjCell(cell, face);
            UINT adjRank = g_tychoMesh->getAdjRank(cell,face);
            
            // In local mesh
            if (adjCell == TychoMesh::BOUNDARY_FACE && adjRank != TychoMesh::BAD_RANK) {
                for (UINT fvrtx = 0; fvrtx < g_nVrtxPerFace; fvrtx++) {
                    Array[ArrayIndex] = 
                        psiBound(group, fvrtx, angle, g_tychoMesh->getSide(cell,face));
                    ArrayIndex++;
                }
            }
        }
    }}}}
}


/*
    arrayToPsiBound
*/
void arrayToPsiBound(const PetscScalar Array[], PsiBoundData &psiBound)
{
    //Start an array index
    int ArrayIndex = 0;
      
    //Incoming flux
    for (UINT angle = 0; angle < g_nAngles; ++angle) {
    for (UINT cell = 0; cell < g_nCells; ++cell) {
    for (UINT group = 0; group < g_nGroups; group++) {
    for (UINT face = 0; face < g_nFacePerCell; face++) {
        if (g_tychoMesh->isIncoming(angle, cell, face)) {
            UINT adjCell = g_tychoMesh->getAdjCell(cell, face);
            UINT adjRank = g_tychoMesh->getAdjRank(cell,face);
            
            // In local mesh
            if (adjCell == TychoMesh::BOUNDARY_FACE && adjRank != TychoMesh::BAD_RANK) {
                for (UINT fvrtx = 0; fvrtx < g_nVrtxPerFace; fvrtx++) {
                    psiBound(group, fvrtx, angle, g_tychoMesh->getSide(cell,face)) = Array[ArrayIndex];
                    ArrayIndex++;
                }
            }
        }
    }}}}
}


/*
    getVecSize
*/
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
    traverseGraph(maxComputePerStep, *s_sweepData, doComm, MPI_COMM_WORLD, 
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
    solve
*/
void SweeperSchur::solve()
{
    PsiData psi;
    PsiData totalSource;
    SourceIteration::solve(this, psi, totalSource);
}


/*
    sweep
    
    Run the Krylov solver
*/
void SweeperSchur::sweep(PsiData &psi, const PsiData &source)
{
    // Initalize variables
    Mat A;
    Vec x, b;
    KSP ksp;
    int argc = 0;
    char **args = NULL;
    PetscInt vecSize;
    PetscInt totalVecSize;
    
    PsiData zeroSource;
    Mat2<UINT> priorities(g_nCells, g_nAngles);

    
    // Make a zero source for use in the inner loop
    zeroSource.setToValue(0.0);


    // The sweep data for Schur iterations
    s_sweepData = new SweepData(psi, zeroSource, g_sigmaTotal, priorities);


    // Set variables
    s_commSides = &c_commSides;
    
    
    /* Petsc setup:
      
      Set up petsc for the Krylov solve.
      In this section, the vectors are x and b are intialized so that 
      they can be used in the Krylov solve, and a matrix shell A is 
      set up so that it runs a matrix free method from the Operator.cc file. 
    */

    
    // Local and global vector sizes
    vecSize = getVecSize();
    UINT totalSize = vecSize;
    Comm::gsum(totalSize);
    totalVecSize = totalSize;


    // Start up petsc
    PetscInitialize(&argc, &args, (char*)0, NULL);


    // Create vectors
    VecCreate(MPI_COMM_WORLD, &x);
    VecSetSizes(x, vecSize, totalVecSize);
    VecSetType(x, VECMPI);
    VecDuplicate(x, &b);
    

    // Create matrix shell and define it as the operator
    MatCreateShell(MPI_COMM_WORLD, vecSize, vecSize, totalVecSize, totalVecSize,
                   (void*)(NULL), &A);
    MatShellSetOperation(A, MATOP_MULT, (void(*)(void))Schur);
    

    // Set solver context
    KSPCreate(MPI_COMM_WORLD, &ksp);
    
 
    // Set operator to KSP context. No preconditioning will actually
    // be used due to the PCNONE option.
    KSPSetOperators(ksp, A, A);
    KSPSetTolerances(ksp, g_ddErrMax, PETSC_DEFAULT, PETSC_DEFAULT, g_ddIterMax);
                

    // Input psi into XIn
    PetscScalar *XIn;
    VecGetArray(x, &XIn);  
    psiBoundToArray(XIn, s_sweepData->getSideData());
    VecRestoreArray(x, &XIn);
    
    
    // Do a sweep on the source 
    SweepData sweepSourceData(psi, source, g_sigmaTotal, priorities);
    sweepSourceData.zeroSideData();
    if (Comm::rank() == 0) {
        printf("    Sweeping Source\n");
    }
    traverseGraph(maxComputePerStep, sweepSourceData, doComm, MPI_COMM_WORLD, 
                  Direction_Forward);
    if (Comm::rank() == 0) {
        printf("    Source Swept\n");
    }


    // Input source into BIn
    c_commSides.commSides(sweepSourceData);
    PetscScalar *BIn;
    VecGetArray(b, &BIn);
    psiBoundToArray(BIn, sweepSourceData.getSideData());
    VecRestoreArray(b, &BIn);


    // Solve the system (x is the solution, b is the RHS)
    if (Comm::rank() == 0) {
        printf("    Starting Krylov Solve on Boundary\n");
    }
    KSPSolve(ksp, b, x);
    

    // Print some stats
    PetscReal rnorm;
    PetscInt its;
    KSPGetIterationNumber(ksp, &its);
    KSPGetResidualNorm(ksp, &rnorm);
    if (Comm::rank() == 0) {
        printf("    Krylov iterations: %u with Rnorm: %lf\n", its, rnorm);
    }


    // Put x in XOut and output the answer from XOut to psi
    PetscScalar *XOut;
    VecGetArray(x, &XOut);
    arrayToPsiBound(XOut, s_sweepData->getSideData());
    VecRestoreArray(x, &XOut);


    // Sweep to solve for the non-boundary values
    if (Comm::rank() == 0) {
        printf("    Sweeping to solve non-boundary values\n");
    }
    traverseGraph(maxComputePerStep, *s_sweepData, doComm, MPI_COMM_WORLD, 
                  Direction_Forward);
    if (Comm::rank() == 0) {
        printf("    Non-boundary values swept\n");
    }
    

    // Destroy vectors, ksp, and matrices to free work space
    VecDestroy(&x);
    VecDestroy(&b);
    MatDestroy(&A);
    KSPDestroy(&ksp);

    
    // Destroy Petsc
    PetscFinalize();
}

#endif
