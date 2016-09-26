#if USE_PETSC

static char help[] = "Solves using Schur Complement.\n\n";

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
#include <omp.h>
#include <vector>
#include <math.h>
#include <string.h>
#include <petscmat.h>
#include <petscvec.h>
#include <petscksp.h>
#include "SweepData.hh"
#include "CommSides.hh"


//Initial definition of static variables
static std::vector<UINT> s_adjRanks;
static std::vector<std::vector<CommSides::MetaData>> s_sendMetaData;
static std::vector<UINT> s_numSendPackets;
static std::vector<UINT> s_numRecvPackets; 
static PsiData s_psi;
static PsiBoundData s_psiBound;
static double s_sigmaTotal;
static const bool doComm = false;
static const UINT maxComputePerStep = std::numeric_limits<uint64_t>::max();
static int argc = 0;
static char **args = NULL;
static PsiData ZeroSource;
static PetscInt VecSize = 0;


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
    //Start an array index
    int ArrayIndex = 0;
      
    // Incoming flux
    for (UINT angle = 0; angle < g_nAngles; ++angle) {
    for (UINT cell = 0; cell < g_nCells; ++cell) {
    for (UINT group = 0; group < g_nGroups; group++) {
    for (UINT face = 0; face < g_nFacePerCell; face++) {
        if (g_tychoMesh->isIncoming(angle, cell, face)) {
            UINT neighborCell = g_tychoMesh->getAdjCell(cell, face);
            UINT adjRank = g_tychoMesh->getAdjRank(cell,face);
            // In local mesh
            if (neighborCell == TychoMesh::BOUNDARY_FACE &&  adjRank != TychoMesh::BAD_RANK) {
                for (UINT fvrtx = 0; fvrtx < g_nVrtxPerFace; fvrtx++) {
                    ArrayIndex++;
                }
            }
        }
    }}}}

    return(ArrayIndex);
}


/*
  Schur 
 
  This performs the sweep and returns the boundary
 */
PetscErrorCode Schur(Mat mat, Vec x, Vec b)
{
   
    //Vector -> array
    const PetscScalar *In;
    VecGetArrayRead(x, &In);
    arrayToPsiBound(In, s_psiBound);
    VecRestoreArrayRead(x,&In);
    
 
    //Create a sweepData
    SweepData s_sweepData(s_psi, ZeroSource, s_psiBound, s_sigmaTotal);
   

    //Traverse the graph to get the values on the outward facing boundary, call commSides to transfer boundary data
    traverseGraph(maxComputePerStep, s_sweepData, doComm, MPI_COMM_WORLD, Direction_Forward);     
    CommSides::commSides(s_adjRanks, s_sendMetaData, s_numSendPackets, 
                         s_numRecvPackets, s_sweepData);	

     
    //Take the values of s_psi from the sweep and put them in b
    PetscScalar *Out;
    VecGetArray(b, &Out);    
    psiBoundToArray(Out, s_sweepData.getSideData());
    VecRestoreArray(b, &Out);


    //b = x - b (aka Out = In - Out)
    VecAXPBY(b,1,-1,x);
  
    
    //Free allocated memory
    free(Out);
   

    return(0);
}



/*
 
 Run the Krylov solver
 
*/
void SweeperSchur::sweep(PsiData &psi, PsiData &source)
{

    //Initalize variables
    Mat A;     //This will be used as a shell for the Operator
    Vec x, b;  //Guess, RHS    
    KSP ksp;    
    PC pc; //preconditioner (will be set to "None" in this case)
    PetscMPIInt size;

    //Make a ZeroSource for use in the inner loop
    ZeroSource = source;
    ZeroSource.setToValue(0.0);  


    //Get rank numbers
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
       
    //Set initial guess for psiBound
    PsiBoundData psiBound;  
    psiBound.setToValue(0.0);
 
    //Put old psi values into PsiBound as a guess
  //  for (UINT angle = 0; angle < g_nAngles; ++angle) {
 //   for (UINT cell = 0; cell < g_nCells; ++cell) {
  //  for (UINT group = 0; group < g_nGroups; group++) {
   // for (UINT face = 0; face < g_nFacePerCell; face++) {
    //    if (g_tychoMesh->isIncoming(angle, cell, face)) {
    //        UINT adjCell = g_tychoMesh->getAdjCell(cell, face);
     //       UINT adjRank = g_tychoMesh->getAdjRank(cell,face);
            
            // In local mesh
     //       if (adjCell == TychoMesh::BOUNDARY_FACE && adjRank != TychoMesh::BAD_RANK) {
      //          for (UINT fvrtx = 0; fvrtx < g_nVrtxPerFace; fvrtx++) {
      //              psiBound(group, fvrtx, angle, g_tychoMesh->getSide(cell,face)) = psi(group, g_tychoMesh->getFaceToCellVrtx(cell,face,fvrtx), angle, cell);
        //        }
        //    }
       // }
   // }}}}
     

    PsiBoundData zeroPsiBound;   


    //Get adjacent ranks
    std::vector<UINT> adjRanks;
    for (UINT cell = 0; cell < g_nCells; cell++) {
    for (UINT face = 0; face < g_nFacePerCell; face++) {
        
        UINT adjRank = g_tychoMesh->getAdjRank(cell, face);
        UINT adjCell = g_tychoMesh->getAdjCell(cell, face);
        
        if (adjCell == TychoMesh::BOUNDARY_FACE && 
            adjRank != TychoMesh::BAD_RANK &&
            std::count(adjRanks.begin(), adjRanks.end(), adjRank) == 0)
        {
            adjRanks.push_back(adjRank);
        }
    }}
    
    
    // Populate sendMetaData, numSendPackets, and numRecvPackets
    std::vector<std::vector<CommSides::MetaData>> sendMetaData(adjRanks.size());
    std::vector<UINT> numSendPackets(adjRanks.size());
    std::vector<UINT> numRecvPackets(adjRanks.size());
    
    for (UINT rankIndex = 0; rankIndex < adjRanks.size(); rankIndex++) {
        
        numSendPackets[rankIndex] = 0;
        numRecvPackets[rankIndex] = 0;
        for (UINT cell = 0; cell < g_nCells; cell++) {
        for (UINT face = 0; face < g_nFacePerCell; face++) {
        
            UINT adjRank = g_tychoMesh->getAdjRank(cell, face);        
            if (adjRank == adjRanks[rankIndex]) {
                for (UINT angle = 0; angle < g_nAngles; angle++) {
                    if (g_tychoMesh->isOutgoing(angle, cell, face)) {
                        CommSides::MetaData md;
                        UINT side = g_tychoMesh->getSide(cell, face);
                        md.gSide = g_tychoMesh->getLGSide(side);
                        md.angle = angle;
                        md.cell  = cell;
                        md.face  = face;
                        sendMetaData[rankIndex].push_back(md);
                        
                        numSendPackets[rankIndex]++;
                    }
                    else {
                        numRecvPackets[rankIndex]++;
                    }
                }
            }
        }}
    }

   
    //Set variables
    s_adjRanks = adjRanks;
    s_sendMetaData = sendMetaData;
    s_numSendPackets = numSendPackets;
    s_numRecvPackets = numRecvPackets; 
    s_psi = psi;
    s_psiBound = psiBound;
    s_sigmaTotal = g_sigmaTotal;

    /* Petsc setup:
      
      Set up petsc for the Krylov solve. In this section, the vectors are x and b are intialized so that 
      they can be used in the Krylov solve, and a matrix shell A is set up so that it runs a matrix free method from the Operator.cc file.     */

  
    VecSize = getVecSize();   
    PetscInt TotalVecSize = VecSize;
    Comm::gsum(TotalVecSize);


    //Start up petsc
    PetscInitialize(&argc,&args,(char*)0, help);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

  
    //Create vectors
    VecCreate(MPI_COMM_WORLD,&x);
    VecSetSizes(x,VecSize,TotalVecSize); 
    VecSetType(x, VECMPI);
    VecDuplicate(x,&b);
        

    //Create matrix shell and define it as the operator
    MatCreateShell(MPI_COMM_WORLD,VecSize,VecSize,TotalVecSize,TotalVecSize,(void*)(NULL),&A);
    MatShellSetOperation(A, MATOP_MULT, (void(*)(void))Schur);
    

    //Set solver context
    KSPCreate(MPI_COMM_WORLD,&ksp);
    
 
    //Set operator to KSP context. No preconditioning will actually
    //be used due to the PCNONE option.
    KSPSetOperators(ksp,A,A);
    KSPSetTolerances(ksp,1.e-5,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
                

    //Input psi into XIn
    PetscScalar *XIn;
    VecGetArray(x,&XIn);  
    psiBoundToArray(XIn, psiBound);
    VecRestoreArray(x, &XIn);
    
    
    //Set zeroPsiBound to 0
    zeroPsiBound.setToValue(0.0); 


    //Do a sweep on the source 
    PsiData sourceCopy = source;//I think this is unnecessary, can just use psi
    SweepData sourceData;
    if (rank==0){
        printf("    Sweeping Source\n");
    }
    traverseGraph(maxComputePerStep, sourceData, doComm, MPI_COMM_WORLD, Direction_Forward);
    if (rank == 0){
        printf("    Source Swept\n");  
    }


    //Input source into BIn
    CommSides::commSides(adjRanks, sendMetaData, numSendPackets, 
                         numRecvPackets, sourceData);    
    PetscScalar *BIn; 
    VecGetArray(b,&BIn);
    psiBoundToArray(BIn, sourceData.getSideData()); 
    VecRestoreArray(b,&BIn);
  
 
    //Solve the system (x is the solution, b is the RHS)
    if (rank==0){
        printf("    Starting Krylov Solve on Boundary\n");
    }    
    KSPSolve(ksp,b,x);
    

    //Print some stats
    PetscReal rnorm;    
    PetscInt its;
    KSPGetIterationNumber(ksp, &its); 
    KSPGetResidualNorm(ksp,&rnorm);
    if (rank==0){
        printf("    Krylov iterations: %u with Rnorm: %lf\n", its, rnorm);
    }

    //Put x in XOut and output the answer from XOut to psi  
    PetscScalar *XOut ;
    VecGetArray(x,&XOut);    
    arrayToPsiBound(XOut, psiBound); 
    VecRestoreArray(x,&XOut);    


    //Sweep to solve for the non-boundary values
    if (rank==0){
        printf("    Sweeping to solve non-boundary values\n");
    }
    SweepData sweepData(psi, source, psiBound, g_sigmaTotal);
    traverseGraph(maxComputePerStep, sweepData, doComm, MPI_COMM_WORLD, Direction_Forward);      	   
    if (rank==0){
        printf("    Non-boundary values swept\n");
    }
    

    //Destroy vectors, ksp, and matrices to free work space
    VecDestroy(&x);
    VecDestroy(&b);
    MatDestroy(&A);
    KSPDestroy(&ksp);

    
    //Destroy Petsc
    PetscFinalize();
    

    
}

#endif
