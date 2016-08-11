/*
    SweeperPBJ.cc
    
    Implements parallel block Jacobi solver.
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

#include "SweeperPBJ.hh"
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
#include <fstream>
#include <iostream>
#include "SweepData2.hh"
#include "CommSides.hh"

using namespace std;


/*
    Constructor
*/
SweeperPBJ::SweeperPBJ(const double sigmaTotal)
{
    c_sigmaTotal = sigmaTotal;
}


/*
    sweep
*/
void SweeperPBJ::sweep(PsiData &psi, const PsiData &source)
{
    const bool doComm = false;
    const UINT maxComputePerStep = std::numeric_limits<uint64_t>::max(); ;
    const UINT maxIter = 100;
    const double tolerance = 1e-5;
    
    
    // Set initial guess for psiBound
    PsiData psiBound(g_spTychoMesh->getNSides(), g_quadrature->getNumAngles(), 
                     g_nVrtxPerFace, g_nGroups);
    psiBound.setToValue(0.0);  // Change to something more reasonable.
    
    
    // Set psi0
    PsiData psi0(g_nVrtxPerCell, g_quadrature->getNumAngles(), 
                 g_spTychoMesh->getNCells(), g_nGroups);

    //Set zeroed versions of source and psiBound for computing the residual
    PsiData zeroPsiBound=psiBound; 
    PsiData zeroSource = source;
    zeroSource.setToValue(0.0);
    
    
    // Create SweepData for traversal
    SweepData2 sweepData(psi, source, psiBound, c_sigmaTotal);
    
    
    // Get adjacent ranks
    vector<UINT> adjRanks;
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
    vector<vector<MetaData>> sendMetaData(adjRanks.size());
    vector<UINT> numSendPackets(adjRanks.size());
    vector<UINT> numRecvPackets(adjRanks.size());
    
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
    
    //Sweep to get the swept source for the error convergence (could be moved out a loop also, to before the source iteration)
    PsiData sourceCopy = source;
    SweepData2 sourceData(sourceCopy, source, zeroPsiBound, c_sigmaTotal);
    traverseGraph(maxComputePerStep, sourceData, doComm, MPI_COMM_WORLD, Direction_Forward);
    commSides(adjRanks, sendMetaData, numSendPackets, numRecvPackets, sourceData);
    printf("Source Swept for Error\n");
    

    // Sweep till converged
    UINT iter = 0;
    while (iter < maxIter) {

      psi0 = psi;
        
      // Sweep
      traverseGraph(maxComputePerStep, sweepData, doComm, MPI_COMM_WORLD, 
                      Direction_Forward);      

      //Sweep psi again with source of zero for the residual
      PsiData psiNew = psi;
      SweepData2 resSweepData(psiNew, zeroSource, psiBound, c_sigmaTotal);
      traverseGraph(maxComputePerStep, resSweepData, doComm, MPI_COMM_WORLD, 
                    Direction_Forward);        
      commSides(adjRanks, sendMetaData, numSendPackets, numRecvPackets, 
                  resSweepData);

     // Check tolerance
        double bnorm = 0.0;
        double rnorm = 0.0;
        for (UINT group = 0; group < g_nGroups; ++group) {
        for (UINT cell = 0; cell < g_spTychoMesh->getNCells(); ++cell) {
        for (UINT face=0; face<g_nFacePerCell; face++){
        for (UINT angle = 0; angle < g_quadrature->getNumAngles(); ++angle) {
        if (g_spTychoMesh->isOutgoing(angle, cell, face)) {
            UINT adjCell = g_spTychoMesh->getAdjCell(cell, face);
            UINT adjRank = g_spTychoMesh->getAdjRank(cell,face);
        if (adjCell == TychoMesh::BOUNDARY_FACE && adjRank != TychoMesh::BAD_RANK) {
        for (UINT vertex = 0; vertex < g_nVrtxPerCell; ++vertex) {
            double p0 = psi0(vertex,angle,cell,group);
            double bSwept = sourceCopy(vertex, angle, cell, group);
            double p1 = psiNew(vertex, angle, cell, group);
            
            bnorm +=  (bSwept)*(bSwept);
            rnorm +=  (p0 - p1 - bSwept)*(p0 - p1 - bSwept);       
            
        }}}}}}}

            Comm::gsum(bnorm);
            Comm::gsum(rnorm);
                
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        
        // Communicate
        commSides(adjRanks, sendMetaData, numSendPackets, numRecvPackets, 
                  sweepData);
        
        //First argument is relative tolerance, second is absolute size of the residual norm
        if (rnorm < tolerance*tolerance*bnorm || tolerance*tolerance > rnorm){
        	break;}

        
        // Increment iter
        iter++;
    }
    
    
    // Print statistics
    if (Comm::rank() == 0) {
        printf("      PBJ Iters: %" PRIu64 "\n", iter);
    }
}

/*
    SweeperPBJ::write

    writes psi to a file
 
 */

void SweeperPBJ::write(PsiData &psi, const PsiData &source)
{
    ofstream outputfile ("tests/testPBJ.txt");
    
    for (UINT group = 0; group < g_nGroups; ++group) {
    for (UINT cell = 0; cell < g_spTychoMesh->getNCells(); ++cell) {
    for (UINT angle = 0; angle < g_quadrature->getNumAngles(); ++angle) {
    for (UINT vertex = 0; vertex < g_nVrtxPerCell; ++vertex) {
       outputfile << psi(vertex, angle, cell, group) << '\n' ;

	
    }}}}

    outputfile.close();
}

