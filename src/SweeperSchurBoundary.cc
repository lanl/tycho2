static char help[] = "Solves using Schur Complement.\n\n";

/*     SweeperSchurBoundary.cc
    
    Implements Schur Complement solver.
*/

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
#include "SweeperSchurBoundary.hh"
#include "Global.hh"
#include "TraverseGraph.hh"
#include "Priorities.hh"
#include "Transport.hh"
#include "PsiData.hh"
#include "Comm.hh"
#include "Typedef.hh"
#include <omp.h>
#include <vector>
#include <math.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <petscmat.h>
#include <petscvec.h>
#include <petscksp.h>
#include "SweepData2.hh"
#include "CommSides.hh"


//Initial definition of static variables
static std::vector<UINT> s_adjRanks;
static std::vector<std::vector<MetaData>> s_sendMetaData;
static std::vector<UINT> s_numSendPackets;
static std::vector<UINT> s_numRecvPackets; 
static PsiData s_psi(g_nVrtxPerCell, 0, 0,g_nGroups);
static PsiData s_psiBound(g_nVrtxPerCell, 0,0,g_nGroups);
static double s_sigmaTotal;
static const bool doComm = false;
static const UINT maxComputePerStep = std::numeric_limits<uint64_t>::max();
static int argc = 0;
static char **args = NULL;
static PsiData ZeroSource(0,0,0,0);
static PetscInt VecSize = 0;

/*
    Constructor
*/
SweeperSchurBoundary::SweeperSchurBoundary(const double sigmaTotal)
{
    c_sigmaTotal = sigmaTotal;
}

//Returns an array of the local boundary values 
void psiBoundToArray(PetscScalar Array[], const PsiData &psiBound)
{
    //Start an array index
    int ArrayIndex = 0;
      
    // Incoming flux
    for (UINT angle = 0; angle < g_quadrature->getNumAngles(); ++angle) {
    for (UINT cell = 0; cell < g_spTychoMesh->getNCells(); ++cell) {
    for (UINT group = 0; group < g_nGroups; group++) {
    for (UINT face = 0; face < g_nFacePerCell; face++) {
        if (g_spTychoMesh->isIncoming(angle, cell, face)) {
            UINT adjCell = g_spTychoMesh->getAdjCell(cell, face);
            UINT adjRank = g_spTychoMesh->getAdjRank(cell,face);
            
            // In local mesh
            if (adjCell == TychoMesh::BOUNDARY_FACE && adjRank != TychoMesh::BAD_RANK) {
                for (UINT fvrtx = 0; fvrtx < g_nVrtxPerFace; fvrtx++) {
                    Array[ArrayIndex] = 
                        psiBound(g_spTychoMesh->getSide(cell,face), angle, fvrtx, group);
                    ArrayIndex++;
                }
            }
        }
    }}}}
}

//Sets local boundary values to the values in array
void arrayToPsiBound(const PetscScalar Array[], PsiData &psiBound)
{
    //Start an array index
    int ArrayIndex = 0;
      
    //Incoming flux
    for (UINT angle = 0; angle < g_quadrature->getNumAngles(); ++angle) {
    for (UINT cell = 0; cell < g_spTychoMesh->getNCells(); ++cell) {
    for (UINT group = 0; group < g_nGroups; group++) {
    for (UINT face = 0; face < g_nFacePerCell; face++) {
        if (g_spTychoMesh->isIncoming(angle, cell, face)) {
            UINT adjCell = g_spTychoMesh->getAdjCell(cell, face);
            UINT adjRank = g_spTychoMesh->getAdjRank(cell,face);
            
            // In local mesh
            if (adjCell == TychoMesh::BOUNDARY_FACE && adjRank != TychoMesh::BAD_RANK) {
                for (UINT fvrtx = 0; fvrtx < g_nVrtxPerFace; fvrtx++) {
                    psiBound(g_spTychoMesh->getSide(cell,face), angle, fvrtx, group) = Array[ArrayIndex];
                    ArrayIndex++;
                }
            }
        }
    }}}}
}

//Returns the Vector Size for a Psi and psiData
int GetVecSize()
{
    //Start an array index
    int ArrayIndex = 0;
      
    // Incoming flux
    for (UINT angle = 0; angle < g_quadrature->getNumAngles(); ++angle) {
    for (UINT cell = 0; cell < g_spTychoMesh->getNCells(); ++cell) {
    for (UINT group = 0; group < g_nGroups; group++) {
    for (UINT face = 0; face < g_nFacePerCell; face++) {
        if (g_spTychoMesh->isIncoming(angle, cell, face)) {
            UINT neighborCell = g_spTychoMesh->getAdjCell(cell, face);
            UINT adjRank = g_spTychoMesh->getAdjRank(cell,face);
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


/*Schur complement operator. This performs the sweep and returns the boundary */
PetscErrorCode Schur(Mat mat, Vec x, Vec b){

   
    //Vector -> array
    const PetscScalar *In;
    VecGetArrayRead(x, &In);
    arrayToPsiBound(In, s_psiBound);
    VecRestoreArrayRead(x,&In);
    
 
    //Create a sweepData
    SweepData2 s_sweepData(s_psi, ZeroSource, s_psiBound, s_sigmaTotal);
   
    //Traverse the graph to get the values on the outward facing boundary, call commSides to transfer boundary data
    traverseGraph(maxComputePerStep, s_sweepData, doComm, MPI_COMM_WORLD, Direction_Forward);     
    commSides(s_adjRanks, s_sendMetaData, s_numSendPackets, s_numRecvPackets, s_sweepData);	


     
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
void SweeperSchurBoundary::sweep(PsiData &psi, PsiData &source)
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
    
       
 /* 
       Set up the sweep data to place into x and b  */
    //Set initial guess for psiBound
    PsiData psiBound(g_spTychoMesh->getNSides(), g_quadrature->getNumAngles(), 
                     g_nVrtxPerFace, g_nGroups);  
        
    psiBound.setToValue(0.0);
 
    //Put old psi values into PsiBound as a guess
  //  for (UINT angle = 0; angle < g_quadrature->getNumAngles(); ++angle) {
 //   for (UINT cell = 0; cell < g_spTychoMesh->getNCells(); ++cell) {
  //  for (UINT group = 0; group < g_nGroups; group++) {
   // for (UINT face = 0; face < g_nFacePerCell; face++) {
    //    if (g_spTychoMesh->isIncoming(angle, cell, face)) {
    //        UINT adjCell = g_spTychoMesh->getAdjCell(cell, face);
     //       UINT adjRank = g_spTychoMesh->getAdjRank(cell,face);
            
            // In local mesh
     //       if (adjCell == TychoMesh::BOUNDARY_FACE && adjRank != TychoMesh::BAD_RANK) {
      //          for (UINT fvrtx = 0; fvrtx < g_nVrtxPerFace; fvrtx++) {
      //              psiBound(g_spTychoMesh->getSide(cell,face), angle, fvrtx, group) = psi(g_spTychoMesh->getFaceToCellVrtx(cell,face,fvrtx), angle, cell, group);
        //        }
        //    }
       // }
   // }}}}
     

    PsiData zeroPsiBound(g_spTychoMesh->getNSides(), g_quadrature->getNumAngles(), 
                     g_nVrtxPerFace, g_nGroups);   


    //Get adjacent ranks
    std::vector<UINT> adjRanks;
    for (UINT cell = 0; cell < g_spTychoMesh->getNCells(); cell++) {
    for (UINT face = 0; face < g_nFacePerCell; face++) {
        
        UINT adjRank = g_spTychoMesh->getAdjRank(cell, face);
        UINT adjCell = g_spTychoMesh->getAdjCell(cell, face);
        
        if (adjCell == TychoMesh::BOUNDARY_FACE && 
            adjRank != TychoMesh::BAD_RANK &&
            std::count(adjRanks.begin(), adjRanks.end(), adjRank) == 0)
        {
            adjRanks.push_back(adjRank);
        }
    }}
    
    
    // Populate sendMetaData, numSendPackets, and numRecvPackets
    std::vector<std::vector<MetaData>> sendMetaData(adjRanks.size());
    std::vector<UINT> numSendPackets(adjRanks.size());
    std::vector<UINT> numRecvPackets(adjRanks.size());
    
    for (UINT rankIndex = 0; rankIndex < adjRanks.size(); rankIndex++) {
        
        numSendPackets[rankIndex] = 0;
        numRecvPackets[rankIndex] = 0;
        for (UINT cell = 0; cell < g_spTychoMesh->getNCells(); cell++) {
        for (UINT face = 0; face < g_nFacePerCell; face++) {
        
            UINT adjRank = g_spTychoMesh->getAdjRank(cell, face);        
            if (adjRank == adjRanks[rankIndex]) {
                for (UINT angle = 0; angle < g_quadrature->getNumAngles(); angle++) {
                    if (g_spTychoMesh->isOutgoing(angle, cell, face)) {
                        MetaData md;
                        UINT side = g_spTychoMesh->getSide(cell, face);
                        md.gSide = g_spTychoMesh->getLGSide(side);
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
    s_sigmaTotal = c_sigmaTotal;

    /* Petsc setup:
      
      Set up petsc for the Krylov solve. In this section, the vectors are x and b are intialized so that 
      they can be used in the Krylov solve, and a matrix shell A is set up so that it runs a matrix free method from the Operator.cc file.     */

  
    VecSize = GetVecSize();   
    PetscInt g_VecSize = VecSize;
    Comm::gsum(g_VecSize);


    //Start up petsc
    PetscInitialize(&argc,&args,(char*)0, help);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

  
    //Create vectors
    VecCreate(MPI_COMM_WORLD,&x);
    VecSetSizes(x,VecSize,g_VecSize); 
    VecSetType(x, VECMPI);
    VecDuplicate(x,&b);
        

    //Create matrix shell and define it as the operator
    MatCreateShell(MPI_COMM_WORLD,VecSize,VecSize,g_VecSize,g_VecSize,(void*)(NULL),&A);
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
    SweepData2 sourceData(sourceCopy, source, zeroPsiBound, c_sigmaTotal);
    if (rank==0){
    printf("    Sweeping Source\n");
    }
    traverseGraph(maxComputePerStep, sourceData, doComm, MPI_COMM_WORLD, Direction_Forward);
    if (rank == 0){
    printf("    Source Swept\n");  
    }


    //Input source into BIn
    commSides(adjRanks, sendMetaData, numSendPackets, numRecvPackets, sourceData);    
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
    printf("    Krylov iterations: %u for Rank: %u with Rnorm: %lf\n", its, rank, rnorm);


    //Put x in XOut and output the answer from XOut to psi  
    PetscScalar *XOut ;
    VecGetArray(x,&XOut);    
    arrayToPsiBound(XOut, psiBound); 
    VecRestoreArray(x,&XOut);    


    //Sweep to solve for the non-boundary values
    if (rank==0){
    printf("    Sweeping to solve non-boundary values\n");
    }
    SweepData2 sweepData(psi, source, psiBound, c_sigmaTotal);
    traverseGraph(maxComputePerStep, sweepData, doComm, MPI_COMM_WORLD, Direction_Forward);      	   

    //Destroy vectors, ksp, and matrices to free work space
    VecDestroy(&x);
    VecDestroy(&b);
    MatDestroy(&A);
    KSPDestroy(&ksp);

    
    //Destroy Petsc
    PetscFinalize();
    

    
}


void SweeperSchurBoundary::write(PsiData &psi, const PsiData &source)
{
    std::ofstream outputfile("tests/testSchurKrylov.txt");
    
    for (UINT group = 0; group < g_nGroups; ++group) {
    for (UINT cell = 0; cell < g_spTychoMesh->getNCells(); ++cell) {
    for (UINT angle = 0; angle < g_quadrature->getNumAngles(); ++angle) {
    for (UINT vertex = 0; vertex < g_nVrtxPerCell; ++vertex) {
       outputfile << psi(vertex, angle, cell, group) << '\n' ;

	
    }}}}

    outputfile.close();
}



