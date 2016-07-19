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
#include "SweepDataSchur.hh"
#include "Operator.hh"


//Initial definition of static variables
std::vector<UINT> Operator::c_adjRanks = {0};
std::vector<std::vector<MetaData>> Operator::c_sendMetaData = {{{0,0,0,0}}};
std::vector<UINT> Operator::c_numSendPackets = {0};
std::vector<UINT> Operator::c_numRecvPackets = {0};  
PsiData Operator::c_psi(0,0,0,0);
PsiData Operator::c_psiBound(0,0,0,0);
PsiData Operator::c_psiSource(0,0,0,0);
double Operator::c_sigmaTotal = 0;
MPI_Comm Operator::c_mpiComm = NULL;


/*
    Constructor
*/
SweeperSchurBoundary::SweeperSchurBoundary(const double sigmaTotal)
{
    c_sigmaTotal = sigmaTotal;
}


/*
    commSides
*/
void commSides(const std::vector<UINT> &adjRanks,
               const std::vector<std::vector<MetaData>> &sendMetaData,
               const std::vector<UINT> &numSendPackets,
               const std::vector<UINT> &numRecvPackets,
               SweepDataSchur &sweepData)
{
    int mpiError;
    UINT numToRecv;
    UINT numAdjRanks = adjRanks.size();
    UINT packetSize = 2 * sizeof(UINT) + sweepData.getDataSize();
    std::vector<MPI_Request> mpiRecvRequests(numAdjRanks);
    std::vector<MPI_Request> mpiSendRequests(numAdjRanks);
    std::vector<std::vector<char>> dataToSend(numAdjRanks);
    std::vector<std::vector<char>> dataToRecv(numAdjRanks);
    
    
    // Data structures to send/recv packets
    for (UINT rankIndex = 0; rankIndex < numAdjRanks; rankIndex++) {
        dataToSend[rankIndex].resize(packetSize * numSendPackets[rankIndex]);
        dataToRecv[rankIndex].resize(packetSize * numRecvPackets[rankIndex]);
    }
    
    
    // Irecv data
    numToRecv = 0;
    for (UINT rankIndex = 0; rankIndex < numAdjRanks; rankIndex++) {
        
        if (dataToRecv[rankIndex].size() > 0) {
            int tag = 0;
            int adjRank = adjRanks[rankIndex];
            MPI_Request request;
            mpiError = MPI_Irecv(dataToRecv[rankIndex].data(), 
                                 dataToRecv[rankIndex].size(), 
                                 MPI_BYTE, adjRank, tag, MPI_COMM_WORLD,
                                 &request);
            Insist(mpiError == MPI_SUCCESS, "");
            mpiRecvRequests[rankIndex] = request;
            numToRecv++;
        }
        
        else {
            mpiRecvRequests[rankIndex] = MPI_REQUEST_NULL;
        }
    }
    
    
    // Update data to send and Isend it
    for (UINT rankIndex = 0; rankIndex < numAdjRanks; rankIndex++) {
        
        if (dataToSend[rankIndex].size() > 0) {
            for (UINT metaDataIndex = 0; 
                 metaDataIndex < sendMetaData[rankIndex].size(); 
                 metaDataIndex++)
            {
                UINT gSide = sendMetaData[rankIndex][metaDataIndex].gSide;
                UINT angle = sendMetaData[rankIndex][metaDataIndex].angle;
                UINT cell  = sendMetaData[rankIndex][metaDataIndex].cell;
                UINT face  = sendMetaData[rankIndex][metaDataIndex].face;
                const char *data = sweepData.getData(cell, face, angle);
                
                char *ptr = &dataToSend[rankIndex][metaDataIndex * packetSize];
                memcpy(ptr, &gSide, sizeof(UINT));
                ptr += sizeof(UINT);
                memcpy(ptr, &angle, sizeof(UINT));
                ptr += sizeof(UINT);
                memcpy(ptr, data, sweepData.getDataSize());
            }
            
            int tag = 0;
	    int adjRank = adjRanks[rankIndex];
            MPI_Request request;
            mpiError = MPI_Isend(dataToSend[rankIndex].data(), 
                                 dataToSend[rankIndex].size(), 
                                 MPI_BYTE, adjRank, tag, MPI_COMM_WORLD, 
                                 &request);
            Insist(mpiError == MPI_SUCCESS, "");
            mpiSendRequests[rankIndex] = request;
        }
        
        else {
            mpiSendRequests[rankIndex] = MPI_REQUEST_NULL;
        }
    }
    
    
    // Get data from Irecv
    for (UINT numWaits = 0; numWaits < numToRecv; numWaits++) {
        
        // Wait for a data packet to arrive
        int rankIndex;
        mpiError = MPI_Waitany(mpiRecvRequests.size(), mpiRecvRequests.data(), 
                               &rankIndex, MPI_STATUS_IGNORE);
        Insist(mpiError == MPI_SUCCESS, "");
        
        
        // Process Data
        UINT numPackets = dataToRecv[rankIndex].size() / packetSize;
        for (UINT packetIndex = 0; packetIndex < numPackets; packetIndex++) {
            char *ptr = &dataToRecv[rankIndex][packetIndex * packetSize];
            UINT gSide = 0;
            UINT angle = 0;
            memcpy(&gSide, ptr, sizeof(UINT));
            ptr += sizeof(UINT);
            memcpy(&angle, ptr, sizeof(UINT));
            ptr += sizeof(UINT);
            UINT side = g_spTychoMesh->getGLSide(gSide);
            sweepData.setSideData(side, angle, ptr);
        }
    }
    
    
    // Wait on send to complete
    if (mpiSendRequests.size() > 0) {
        mpiError = MPI_Waitall(mpiSendRequests.size(), mpiSendRequests.data(), 
                               MPI_STATUSES_IGNORE);
        Insist(mpiError == MPI_SUCCESS, "");
    }
}




/*Schur complement operator. This performs the sweep and returns the outgoing boundary */
PetscErrorCode Schur(Mat mat, Vec x, Vec b){


    //Retrieve variables
    Operator Op;
    MPI_Comm s_mpiComm = Op.getmpiComm();
    const std::vector<UINT> s_adjRanks = Op.getadjRanks();
    const std::vector<std::vector<MetaData>> s_sendMetaData = Op.getsendMetaData();
    const std::vector<UINT> s_numSendPackets = Op.getnumSendPackets();
    const std::vector<UINT> s_numRecvPackets = Op.getnumRecvPackets();  
    PsiData s_psi = Op.getpsi();
    PsiData s_psiBound = Op.getpsiBound();
    PsiData s_psiSource = Op.getpsiSource();
    double s_sigmaTotal = Op.getsigmaTotal();
    SweepDataSchur s_sweepData(s_psi, s_psiSource, s_psiBound, s_sigmaTotal);
    
    PsiData psiOld = Op.getpsi();//!!!!!!!
    

    //Declare variables
    PetscErrorCode ierr;
    const bool doComm = false;
    const UINT maxComputePerStep = std::numeric_limits<uint64_t>::max(); ;
    PetscInt size;
    

    //Vector -> array
    VecGetSize(x, &size);
    const PetscScalar *tempIn;
    VecGetArrayRead(x, &tempIn);

    
    //Set data back
    VecRestoreArrayRead(x,&tempIn);
    
      
    //Traverse the graph to get the values on the outward facing boundary, call commSides to transfer boundary data
    traverseGraph(maxComputePerStep, s_sweepData, doComm, s_mpiComm, Direction_Forward);      	
    commSides(s_adjRanks, s_sendMetaData, s_numSendPackets, s_numRecvPackets, s_sweepData);	



//    count = 0;
//    PetscScalar p1;
//    for (UINT group = 0; group < g_nGroups; ++group) {
//    for (UINT cell = 0; cell < g_spTychoMesh->getNCells(); ++cell) {
//    for (UINT face = 0; face < g_nFacePerCell; ++face) {
//	     UINT adjCell = g_spTychoMesh->getAdjCell(cell, face);
//	     UINT adjRank = g_spTychoMesh->getAdjRank(cell, face);
//             for (UINT angle = 0; angle < g_quadrature->getNumAngles(); ++angle) {
//		     if ((adjCell == TychoMesh::BOUNDARY_FACE && g_spTychoMesh->isOutgoing(angle, cell, face))||(adjCell == TychoMesh::BOUNDARY_FACE && adjRank == TychoMesh::BAD_RANK)){
//		             for (UINT vertex = 0; vertex < g_nVrtxPerCell; ++vertex) {
//		   		    p1 = s_psi(vertex, angle, cell, group);
//  				    //ierr = VecSetValues(b,1,&count,&p1,INSERT_VALUES); //n=1??!!!
 //   				    tempOut[count] = p1;
   //                                 count += 1;
          
//	}}}}}}
     

    //Initialize an array to put into b
    PetscScalar Out[size];


    //Take the values of s_psi from the sweep and put them in b
    PetscInt count = 0;    
    for (UINT group = 0; group < g_nGroups; ++group) {
    for (UINT cell = 0; cell < g_spTychoMesh->getNCells(); ++cell) {
    for (UINT face = 0; face < g_nFacePerCell; face++) {
         UINT adjCell = g_spTychoMesh->getAdjCell(cell, face);
         for (UINT angle = 0; angle < g_quadrature->getNumAngles(); ++angle) {
         if (adjCell == TychoMesh::BOUNDARY_FACE && g_spTychoMesh->isIncoming(angle, cell, face)){
             for (UINT vertex = 0; vertex < g_nVrtxPerCell; ++vertex){
                  Out[count] = s_psi(vertex, angle, cell, group);
	          count += 1;
    }}}}}}
     
  
    //Produce b
    PetscInt ind;
    PetscScalar tempOut;
    for (ind=0; ind<count; ind++){
        tempOut = Out[ind];
        VecSetValue(b, ind, tempOut, INSERT_VALUES);
    }


    //Set values back into Op
    Op.setmpiComm(s_mpiComm);
    Op.setadjRanks(s_adjRanks);
    Op.setsendMetaData(s_sendMetaData);
    Op.setnumSendPackets(s_numSendPackets);
    Op.setnumRecvPackets(s_numRecvPackets);  
    Op.setpsi(s_psi);
    Op.setpsiBound(s_psiBound);
    Op.setpsiSource(s_psiSource);
    Op.setsigmaTotal(s_sigmaTotal);


        // Check tolerance and set psi0 = psi1
        double errL1 = 0.0;
        double normL1 = 0.0;
        for (UINT group = 0; group < g_nGroups; ++group) {
        for (UINT cell = 0; cell < g_spTychoMesh->getNCells(); ++cell) {
        for (UINT angle = 0; angle < g_quadrature->getNumAngles(); ++angle) {
        for (UINT vertex = 0; vertex < g_nVrtxPerCell; ++vertex) {
            double p0 = s_sweepData.getData2(vertex, angle, cell, group);
            double p1 = s_psi(vertex, angle, cell, group);
             
            errL1  += fabs(p0 - p1);
            normL1 += fabs(p1);
        }}}}
    
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if(rank==0){
	printf("	Total errL1: %f \n", errL1);
	printf("	Total normL1: %f \n", normL1);
	printf("	errL1/normL1: %f \n", errL1/normL1);
        }
        
   
    return(0);
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
    PetscBool  nonzeroguess = PETSC_FALSE;
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
    Operator Op;
    Op.setmpiComm(MPI_COMM_WORLD);
    Op.setadjRanks(adjRanks);
    Op.setsendMetaData(sendMetaData);
    Op.setnumSendPackets(numSendPackets);
    Op.setnumRecvPackets(numRecvPackets); 
    Op.setpsi(psi);
    Op.setpsiBound(psiBound);
    Op.setpsiSource(source);
    Op.setsigmaTotal(c_sigmaTotal);
    

    /* Petsc setup:
      
      Set up petsc for the Krylov solve. In this section, the vectors are x and b are intialized so that 
      they can be used in the Krylov solve, and a matrix shell A is set up so that it runs a matrix free method from the Operator.cc file.     */


    //Get vec size (would save time to eliminate this step somehow, or could do it once for each rank)
    PetscInt VecSize = 0;
    for (UINT group = 0; group < g_nGroups; ++group) {
    for (UINT cell = 0; cell < g_spTychoMesh->getNCells(); ++cell) {
    for (UINT face = 0; face < g_nFacePerCell; face++) {
        UINT adjCell = g_spTychoMesh->getAdjCell(cell, face);
        for (UINT angle = 0; angle < g_quadrature->getNumAngles(); ++angle) {
        if (adjCell == TychoMesh::BOUNDARY_FACE && g_spTychoMesh->isIncoming(angle, cell, face)){
        for (UINT vertex = 0; vertex < g_nVrtxPerCell; ++vertex) {
            VecSize += 1;
    }}}}}}
    

    //Start up petsc
    PetscInitialize(&argc,&args,(char*)0, help); //argc and args were removed!!!!!!
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    PetscOptionsGetInt(NULL,NULL,"-n",&VecSize,NULL); 
    PetscOptionsGetBool(NULL,NULL,"-nonzero_guess",&nonzeroguess,NULL);

  
    //Create vectors
    VecCreate(PETSC_COMM_SELF,&x);
    VecSetSizes(x,PETSC_DECIDE,VecSize); 
    VecSetFromOptions(x);
    VecDuplicate(x,&b);
        

    //Create matrix shell and define it as the operator
    MatCreateShell(PETSC_COMM_SELF,VecSize,VecSize,VecSize,VecSize,(void*)(NULL),&A);
    //MatSetType(A, MATMPIAIJ);
    MatShellSetOperation(A, MATOP_MULT, (void(*)(void))Schur);
    

    //Set solver context
    KSPCreate(PETSC_COMM_SELF,&ksp);
    
 
    //Set operator to KSP context. No preconditioning will actually
    //be used due to the PCNONE option.
    KSPSetOperators(ksp,A,A);
    KSPGetPC(ksp,&pc);
    PCSetType(pc,PCNONE);
    KSPSetTolerances(ksp,1.e-5,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);


    //
    //KSPSetFromOptions(ksp); //CHECK ON THIS!!!!!
    //if (nonzeroguess) {
    //    PetscScalar p = .5;
    //    VecSet(x,p);
     //    KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
     //    }
    

    //Input psi into temp
    PetscInt count = 0;
    PetscScalar *temp, p;
    VecGetArray(x,&temp);  
    for (UINT group = 0; group < g_nGroups; ++group) {
    for (UINT cell = 0; cell < g_spTychoMesh->getNCells(); ++cell) {
    for (UINT face = 0; face < g_nFacePerCell; face++) {
        UINT adjCell = g_spTychoMesh->getAdjCell(cell, face);
        for (UINT angle = 0; angle < g_quadrature->getNumAngles(); ++angle) {
        if (adjCell == TychoMesh::BOUNDARY_FACE && g_spTychoMesh->isIncoming(angle, cell, face)){
        for (UINT vertex = 0; vertex < g_nVrtxPerCell; ++vertex) {
            p = psi(vertex, angle, cell, group);
	    temp[count] = p;
            count += 1;
    }}}}}}
    

    //Input temp into x
    //VecSetValues(x, count, ind, temp, INSERT_VALUES);
    VecRestoreArray(x,&temp);

                                               
    //Assemble vector x
    VecAssemblyBegin(x);
    VecAssemblyEnd(x);
    


    //Take the values of s_psi from the sweep and put them in b
    //PetscInt count = 0;    
    //for (UINT group = 0; group < g_nGroups; ++group) {
    //for (UINT cell = 0; cell < g_spTychoMesh->getNCells(); ++cell) {
    //for (UINT face = 0; face < g_nFacePerCell; face++) {
     //    UINT adjCell = g_spTychoMesh->getAdjCell(cell, face);
     //    for (UINT angle = 0; angle < g_quadrature->getNumAngles(); ++angle) {
     //    if (adjCell == TychoMesh::BOUNDARY_FACE && g_spTychoMesh->isIncoming(angle, cell, face)){
     //        for (UINT vertex = 0; vertex < g_nVrtxPerCell; ++vertex){
     //             Out[count] = s_psi(vertex, angle, cell, group);
//	          count += 1;
  //  }}}}}}



    //Input source into tempb
    count = 0;
    PetscScalar *tempb, pb;
    VecGetArray(b,&tempb);    
    for (UINT group = 0; group < g_nGroups; ++group) {
    for (UINT cell = 0; cell < g_spTychoMesh->getNCells(); ++cell) {
    for (UINT face = 0; face < g_nFacePerCell; face++) {
         UINT adjCell = g_spTychoMesh->getAdjCell(cell, face);
         for (UINT angle = 0; angle < g_quadrature->getNumAngles(); ++angle) {
         if (adjCell == TychoMesh::BOUNDARY_FACE && g_spTychoMesh->isIncoming(angle, cell, face)){
         for (UINT vertex = 0; vertex < g_nVrtxPerCell; ++vertex) {
             pb = source(vertex, angle, cell, group);
             tempb[count] = pb;
             count += 1;
    }}}}}}//!!!!!!!!!!!!!!!!!!!!!!THIS ONLY NEEDS TO HAPPEN ON THE FIRST RUN THROUGH -> WOULD SAVE TIME TO ELIMINATE
    

    //Input tempb into b
    //VecSetValues(x, count, ind, temp, INSERT_VALUES);
    VecRestoreArray(b,&tempb);
    

    //Assemble vectors
    VecAssemblyBegin(b);
    VecAssemblyEnd(b);
    
   
    //Solve the system (x is the solution, b is the RHS)
    KSPSolve(ksp,b,x);
    
    
    //Get number of iterations
    PetscInt its;
    KSPGetIterationNumber(ksp, &its);
    printf("\n\n%u\n\n", its);


    //Put x in temp
    VecGetArray(x,&temp);
        

    //Output the answer from temp to psi   
    count = 0; 
    for (UINT group = 0; group < g_nGroups; ++group) {
    for (UINT cell = 0; cell < g_spTychoMesh->getNCells(); ++cell) {
    for (UINT angle = 0; angle < g_quadrature->getNumAngles(); ++angle) {
    for (UINT vertex = 0; vertex < g_nVrtxPerCell; ++vertex) {
         psi(vertex, angle, cell, group) = temp[count];
         count += 1;
    }}}}
    

    //View solver info; we could instead use the option -ksp_view to
    // print this info to the screen at the conclusion of KSPSolve().
    KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);


    //Destroy vectors, ksp, and matrices to free work space
    VecDestroy(&x);
    VecDestroy(&b);
    MatDestroy(&A);
    KSPDestroy(&ksp);

    
    //Destroy Petsc
    PetscFinalize();


    //Retrieve variables
    adjRanks = Op.getadjRanks();
    sendMetaData = Op.getsendMetaData();
    numSendPackets = Op.getnumSendPackets();
    numRecvPackets = Op.getnumRecvPackets();  
    PsiData psiNew = Op.getpsi();
    psiBound = Op.getpsiBound();
    c_sigmaTotal = Op.getsigmaTotal();
    
}

