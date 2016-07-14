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
#include "Operator.hh"
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
#include <functional>



PetscErrorCode PassFunc(Mat A, Vec x, Vec b, MPI_Comm mpiComm, const std::vector<UINT> &adjRanks, const std::vector<std::vector<MetaData>> &sendMetaData, const std::vector<UINT> &numSendPackets,const std::vector<UINT> &numRecvPackets, SweepDataSchur &sweepData, PsiData &psi)
{Operator OpGlobal(mpiComm, adjRanks, sendMetaData, numSendPackets, numRecvPackets, sweepData, psi);
return OpGlobal.Schur(A,x,b);
};


/*
    Constructor
*/
SweeperSchurBoundary::SweeperSchurBoundary(const double sigmaTotal)
{
    c_sigmaTotal = sigmaTotal;
}


/*
 
 Run the Krylov solver
 
*/
void SweeperSchurBoundary::sweep(PsiData &psi, const PsiData &source)
{
    
    //Initalize variables
    Mat A;     //This will be used as a shell for the Operator
    Vec x, b;  //Guess, RHS    
    KSP ksp;    
    PC pc; //preconditioner (will be set to "None" in this case)
    PetscMPIInt size;
    PetscScalar  one = 1.0;
    PetscBool  nonzeroguess = PETSC_FALSE;
    PetscInt n = g_nGroups*(g_spTychoMesh->getNCells())*(g_quadrature->getNumAngles())*g_nVrtxPerCell;
    int argc = 0;
    char **args = NULL;
    //const UINT maxIter = 100;
    //const double tolerance = 1e-5; //Equivalent below
     
      
 /* 
       Set up the sweep data to place into x and b  */

    //Set initial guess for psiBound
    PsiData psiBound(g_spTychoMesh->getNSides(), g_quadrature->getNumAngles(), 
                     g_nVrtxPerFace, g_nGroups);
    psiBound.setToValue(0.0);  // Change to something more reasonable.
    psi.setToValue(0.0);                               
                                       
                                                                                                                                      
    // Create SweepData for traversal
    SweepDataSchur sweepData(psi, source, psiBound, c_sigmaTotal);                    
    
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

    

    /* Petsc setup:
      
      Set up petsc for the Krylov solve. In this section, the vectors are x and b are intialized so that 
      they can be used in the Krylov solve, and a matrix shell A is set up so that it runs a matrix free method from the Operator.cc file.     */

    //Start up petsc
    PetscInitialize(&argc,&args,(char*)0, help); //argc and args were removed!!!!!!
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    //if (size != 1) SETERRQ(MPI_COMM_WORLD,1,"This is a uniprocessor example only!"); //Remove this line? MPI should be using more than one processor. Don't want mixed up
    PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL); //WHAT IS N=20 ???!!!!!
    PetscOptionsGetBool(NULL,NULL,"-nonzero_guess",&nonzeroguess,NULL);


    //Create vectors
    VecCreate(PETSC_COMM_WORLD,&x);
    PetscObjectSetName((PetscObject) x, "Solution"); //Change name? Necessary?
    VecSetSizes(x,PETSC_DECIDE,n); 
    VecSetFromOptions(x);
    VecDuplicate(x,&b);
        

    //Create an operator from Operator.cc for the matrix-free method
    PetscErrorCode (Operator::*ptrtry)(Mat, Vec, Vec) = NULL;
    ptrtry = &Operator::Schur;
    Operator Op(MPI_COMM_WORLD, adjRanks, sendMetaData, numSendPackets, numRecvPackets, sweepData, psi);

    //Create matrix shell and define it as the operator
    MatCreateShell(PETSC_COMM_WORLD,n,n,n,n,(void*)(NULL),&A);
    //MatSetType(A, MATMPIAIJ);
    MatShellSetOperation(A, MATOP_MULT, (void(*)(void))(Op.*ptrtry));
    //MatShellSetOperation(A, MATOP_MULT, (void(*)(void))PassFunc(A, x, b, MPI_COMM_WORLD,adjRanks, sendMetaData, numSendPackets, numRecvPackets, sweepData, psi));
   

    //Set solver context
    KSPCreate(PETSC_COMM_WORLD,&ksp);
    
 
    //Set operator to KSP context. No preconditioning will actually
    //be used due to the PCNONE option.
    KSPSetOperators(ksp,A,A);
    KSPGetPC(ksp,&pc);
    PCSetType(pc,PCNONE);
    KSPSetTolerances(ksp,1.e-5,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);


    //
    //KSPSetFromOptions(ksp); //CHECK ON THIS!!!!!
    //if (nonzeroguess) {
        // PetscScalar p = .5;
        // VecSet(x,p);
         //KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
        // }
    
    //Input psi into temp
    PetscInt count = 0;
    UINT arraySize = g_nGroups*(g_spTychoMesh->getNCells())*(g_quadrature->getNumAngles())*g_nVrtxPerCell;
    PetscScalar *temp, p;
    VecGetArray(x,&temp);    //Gets temp to correct size
    for (UINT group = 0; group < g_nGroups; ++group) {
    for (UINT cell = 0; cell < g_spTychoMesh->getNCells(); ++cell) {
    for (UINT angle = 0; angle < g_quadrature->getNumAngles(); ++angle) {
    for (UINT vertex = 0; vertex < g_nVrtxPerCell; ++vertex) {
        p = psi(vertex, angle, cell, group);
	temp[count] = p;
        count += 1;
    }}}}
    
    //Create an index for x
    //PetscInt ind[count];
    //for (UINT i=0;i<count;i++){
    //	ind[i] = count;
    //}

    //Input temp into x
    //VecSetValues(x, count, ind, temp, INSERT_VALUES);
    VecRestoreArray(x,&temp);

                                               
    //Assemble vector x
    VecAssemblyBegin(x);
    VecAssemblyEnd(x);
    
    //Input source into tempb !!Check below
    count = 0;
    PetscScalar *tempb, pb;
    VecGetArray(b,&tempb);    //Gets tempb to correct size
    for (UINT group = 0; group < g_nGroups; ++group) {
    for (UINT cell = 0; cell < g_spTychoMesh->getNCells(); ++cell) {
    for (UINT angle = 0; angle < g_quadrature->getNumAngles(); ++angle) {
    for (UINT vertex = 0; vertex < g_nVrtxPerCell; ++vertex) {
       pb = source(vertex, angle, cell, group);
       tempb[count] = pb;
       count += 1;
    }}}}//!!!!!!!!!!!!!!!!!!!!!!THIS ONLY NEEDS TO HAPPEN ON THE FIRST RUN THROUGH -> WOULD SAVE TIME TO ELIMINATE
    

    //Input tempb into b
    //VecSetValues(x, count, ind, temp, INSERT_VALUES);
    VecRestoreArray(b,&tempb);
    

    //Assemble vectors
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);
    //MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY); !!!!
    //MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY); !!!!
    
   
    //Solve the system (x is the solution, b is the RHS)
    KSPSolve(ksp,b,x);
    

    //Put x in temp
    VecGetArray(x,&temp);
        

    //Output the answer from temp to psi    
    for (UINT group = 0; group < g_nGroups; ++group) {
    for (UINT cell = 0; cell < g_spTychoMesh->getNCells(); ++cell) {
    for (UINT angle = 0; angle < g_quadrature->getNumAngles(); ++angle) {
    for (UINT vertex = 0; vertex < g_nVrtxPerCell; ++vertex) {
         psi(vertex, angle, cell, group) = temp[count];
         count += 1;
    }}}}
    

    //View solver info; we could instead use the option -ksp_view to
    // print this info to the screen at the conclusion of KSPSolve().
    //KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);


    //Destroy vectors, ksp, and matrices to free work space
    VecDestroy(&x);
    VecDestroy(&b);
    MatDestroy(&A);
    KSPDestroy(&ksp);

    
    //Destroy Petsc
    PetscFinalize();
    

}


