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

static PsiData resOld(0,0,0,0);
using namespace std;


/*
    SweepData
    
    Holds psi and other data for the sweep.
*/
class SweepDataPBJ : public TraverseData
{
public:
    
    /*
        Constructor
    */
    SweepDataPBJ(PsiData &psi, const PsiData &source, PsiData &psiBound, 
              const double sigmaTotal)
    : c_psi(psi), c_psiBound(psiBound), 
      c_source(source), c_sigmaTotal(sigmaTotal),
      c_localFaceData(g_nThreads)
    {
        for (UINT angleGroup = 0; angleGroup < g_nThreads; angleGroup++) {
            c_localFaceData[angleGroup] = Mat2<double>(g_nVrtxPerFace, g_nGroups);
        }
    }
    
    
    /*
        data
        
        Return face data for (cell, angle) pair.
    */
    virtual const char* getData(UINT cell, UINT face, UINT angle)
    {
        Mat2<double> &localFaceData = c_localFaceData[omp_get_thread_num()];
        
        for (UINT group = 0; group < g_nGroups; group++) {
        for (UINT fvrtx = 0; fvrtx < g_nVrtxPerFace; fvrtx++) {
            UINT vrtx = g_spTychoMesh->getFaceToCellVrtx(cell, face, fvrtx);
            localFaceData(fvrtx, group) = c_psi(vrtx, angle, cell, group);
        }}
        
        return (char*) (&localFaceData[0]);
    }
    
    
    /*
        getDataSize
    */
    virtual size_t getDataSize()
    {
        return g_nGroups * g_nVrtxPerFace * sizeof(double);
    }
    
    
    /*
        sideData
        
        Set face data for (side, angle) pair.
    */
    virtual void setSideData(UINT side, UINT angle, const char *data)
    {
        Mat2<double> localFaceData(g_nVrtxPerFace, g_nGroups, (double*)data);
        
        for (UINT group = 0; group < g_nGroups; group++) {
        for (UINT fvrtx = 0; fvrtx < g_nVrtxPerFace; fvrtx++) {
            c_psiBound(side, angle, fvrtx, group) = localFaceData(fvrtx, group);
        }}
    }
    
    
    /*
        getPriority
        
        Return a priority for the cell/angle pair.
        Not needed for this class, so it is just set to a constant.
    */
    virtual UINT getPriority(UINT cell, UINT angle)
    {
        UNUSED_VARIABLE(cell);
        UNUSED_VARIABLE(angle);
        return 1;
    }
    
    
    /*
        update
        
        Updates psi for a given (cell, angle) pair.
    */
    virtual void update(UINT cell, UINT angle, 
                        UINT adjCellsSides[g_nFacePerCell], 
                        BoundaryType bdryType[g_nFacePerCell])
    {
        UNUSED_VARIABLE(adjCellsSides);
        UNUSED_VARIABLE(bdryType);
        
        Mat2<double> localSource(g_nVrtxPerCell, g_nGroups);
        Mat2<double> localPsi(g_nVrtxPerCell, g_nGroups);
        Mat3<double> localPsiBound(g_nVrtxPerFace, g_nFacePerCell, g_nGroups);
        
        
        // Populate localSource
        for (UINT group = 0; group < g_nGroups; group++) {
        for (UINT vrtx = 0; vrtx < g_nVrtxPerCell; vrtx++) {
            localSource(vrtx, group) = c_source(vrtx, angle, cell, group);
        }}
        
        
        // Populate localPsiBound
        Transport::populateLocalPsiBound(angle, cell, c_psi, c_psiBound, 
                                         localPsiBound);
        
        
        // Transport solve
        Transport::solve(cell, angle, c_sigmaTotal, 
                         localPsiBound, localSource, localPsi);
        
        
        // localPsi -> psi
        for (UINT group = 0; group < g_nGroups; group++) {
        for (UINT vrtx = 0; vrtx < g_nVrtxPerCell; vrtx++) {
            c_psi(vrtx, angle, cell, group) = localPsi(vrtx, group);
        }}
    }
    
private:
    PsiData &c_psi;
    PsiData &c_psiBound;
    const PsiData &c_source;
    const double c_sigmaTotal;
    vector<Mat2<double>> c_localFaceData;
};


/*
    MetaData struct
*/
struct MetaData
{
    UINT gSide;
    UINT angle;
    UINT cell;
    UINT face;
};

/*
    commSides
*/
void commSides(const vector<UINT> &adjRanks,
               const vector<vector<MetaData>> &sendMetaData,
               const vector<UINT> &numSendPackets,
               const vector<UINT> &numRecvPackets,
               SweepDataPBJ &sweepData)
{
    int mpiError;
    UINT numToRecv;
    UINT numAdjRanks = adjRanks.size();
    UINT packetSize = 2 * sizeof(UINT) + sweepData.getDataSize();
    vector<MPI_Request> mpiRecvRequests(numAdjRanks);
    vector<MPI_Request> mpiSendRequests(numAdjRanks);
    vector<vector<char>> dataToSend(numAdjRanks);
    vector<vector<char>> dataToRecv(numAdjRanks);
    
    
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

    //resOld = psi0;
    PsiData res(g_nVrtxPerCell, g_quadrature->getNumAngles(), 
                g_spTychoMesh->getNCells(), g_nGroups);
    
    PsiData zeroPsiBound=psiBound; 
    PsiData zeroSource = source;
    zeroSource.setToValue(0.0);

    for (UINT group = 0; group < g_nGroups; ++group) {
    for (UINT cell = 0; cell < g_spTychoMesh->getNCells(); ++cell) {
    for (UINT angle = 0; angle < g_quadrature->getNumAngles(); ++angle) {
    for (UINT vertex = 0; vertex < g_nVrtxPerCell; ++vertex) {
        psi0(vertex, angle, cell, group) = psi(vertex, angle, cell, group);
        //resOld(vertex, angle, cell, group) = 0;
        //res(vertex, angle, cell, group) = 0;
    }}}}
    
    
    // Create SweepData for traversal
    SweepDataPBJ sweepData(psi, source, psiBound, c_sigmaTotal);
    
    
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
    SweepDataPBJ sourceData(sourceCopy, source, zeroPsiBound, c_sigmaTotal);
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
      SweepDataPBJ resSweepData(psiNew, zeroSource, psiBound, c_sigmaTotal);
      traverseGraph(maxComputePerStep, resSweepData, doComm, MPI_COMM_WORLD, 
                    Direction_Forward);        
      commSides(adjRanks, sendMetaData, numSendPackets, numRecvPackets, 
                  resSweepData);

     // Check tolerance and set psi0 = psi1
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
            double b = source(vertex, angle, cell, group);
            double p1 = psiNew(vertex, angle, cell, group);
            
            bnorm +=  (bSwept)*(bSwept);
            rnorm +=  (p0 - p1 - bSwept)*(p0 - p1 - bSwept);       
            
            //resOld(vertex, angle, cell, group) = res(vertex,angle,cell,group);
        }}}}}}}

            Comm::gsum(bnorm);
            Comm::gsum(rnorm);
            //double sqrtbnorm = std::pow(bnorm, 0.5);
            //double sqrtrnorm = std::pow(rnorm, 0.5);
                
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        
        // Communicate
        commSides(adjRanks, sendMetaData, numSendPackets, numRecvPackets, 
                  sweepData);
        
        //First argument is relative tolerance, second is absolute size of the residual norm
        if (rnorm < tolerance*tolerance*bnorm || 1e-5*1e-5 > rnorm){
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

