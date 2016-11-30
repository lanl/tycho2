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

#include "Sweeper.hh"
#include "Problem.hh"
#include "SourceIteration.hh"
#include "Quadrature.hh"
#include "Global.hh"
#include "SweepSchedule.hh"
#include "Comm.hh"
#include "Transport.hh"
#include "PsiData.hh"
#include "Timer.hh"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <omp.h>
#include <mpi.h>

using namespace std;


/*
    VG
    
    Converts (vrtx,group) index to single index.
*/
static inline
UINT VG(UINT vrtx, UINT group)
{
    return vrtx + g_nVrtxPerFace * group;
}


/*
    getBoundData
    
    Puts psiBound data into communication data structure.
*/
static
void getBoundData(const PsiBoundData &psiBound, const UINT side,
                  const UINT angle, vector<double> &psiSide)
{
    for (UINT vertex = 0; vertex < g_nVrtxPerFace; ++vertex) {
    for (UINT group = 0; group < g_nGroups; ++group) {
        psiSide[VG(vertex, group)] = psiBound(group, vertex, angle, side);
    }}
}


/*
    setBoundData
    
    Sets psiBound from the communication data structures.
*/
static
void setBoundData(PsiBoundData &psiBound, const vector<UINT> &commSides,
                  const vector<UINT> &commAngles, const vector<double> &commPsi) 
{
    for(unsigned i = 0; i < commSides.size(); i++) {
        UINT side = commSides[i];
        UINT angle = commAngles[i];
        UINT commPsiIndex = i * g_nVrtxPerFace * g_nGroups;
        for (UINT vrtx = 0; vrtx < g_nVrtxPerFace; ++vrtx) {
        for (UINT group = 0; group < g_nGroups; ++group) {
            psiBound(group, vrtx, angle, side) = commPsi[commPsiIndex + VG(vrtx, group)];
        }}
    }
}


/*
    send
    
    Send side data.
*/
static
void send(const UINT step, const UINT angleGroup,
          Mat2<vector<UINT>> &commSidesAngles, Mat2<vector<double>> &commPsi,
          vector<MPI_Request> &mpiRequests)
{
    if (step < g_sweepSchedule[angleGroup]->nSteps()) {
        for (UINT proc : g_sweepSchedule[angleGroup]->getSendProcs(step)) {
            
            MPI_Request mpiRequest;
            UINT nSides = commSidesAngles(angleGroup, proc).size() / 2;
            UINT nData = commPsi(angleGroup, proc).size();
            vector<UINT> nSidesData = {nSides, nData};
            int tag0 = angleGroup * 3 + 0;
            int tag1 = angleGroup * 3 + 1;
            int tag2 = angleGroup * 3 + 2;
    
            Comm::iSendUIntVector(nSidesData, proc, tag0, mpiRequest);
            mpiRequests.push_back(mpiRequest);
            if (nSides > 0) {
                Comm::iSendUIntVector(commSidesAngles(angleGroup, proc), proc, 
                                     tag1, mpiRequest);
                mpiRequests.push_back(mpiRequest);
                Comm::iSendDoubleVector(commPsi(angleGroup, proc), proc, 
                                        tag2, mpiRequest);
                mpiRequests.push_back(mpiRequest);
            }
            
            commSidesAngles(angleGroup, proc).clear();
            commPsi(angleGroup, proc).clear();
        }
    }
}


/*
    recv
    
    Receive side data.
*/
static
void recv(const UINT step, const UINT angleGroup, PsiBoundData &psiBound)
{
    if (step < g_sweepSchedule[angleGroup]->nSteps()) {
        for (UINT proc : g_sweepSchedule[angleGroup]->getRecvProcs(step)) {
            
            // Get num sides and data size
            vector<UINT> nSidesData(2);
            int tag0 = angleGroup * 3 + 0;
            int tag1 = angleGroup * 3 + 1;
            int tag2 = angleGroup * 3 + 2;
    
            Comm::recvUIntVector(nSidesData, proc, tag0);
            UINT nSides = nSidesData[0];
            UINT nData = nSidesData[1];
    
            if (nSides > 0) {
    
                // Receive data
                vector<UINT> commSidesAngles(2*nSides);
                vector<double> commPsi(nData);
                Comm::recvUIntVector(commSidesAngles, proc, tag1);
                Comm::recvDoubleVector(commPsi, proc, tag2);
    
                // Pull apart side and angle data
                // Convert side data to local side indexing
                vector<UINT> commSides(nSides);
                vector<UINT> commAngles(nSides);
                for (UINT side = 0; side < nSides; ++side) {
                    commSides[side] = 
                        g_tychoMesh->getGLSide(commSidesAngles[2*side]);
                    commAngles[side] = commSidesAngles[2*side+1];
                }

                // Set the boundary data
                setBoundData(psiBound, commSides, commAngles, commPsi);
            }
        }
    }
}


/*
    updateComm
    
    Updates data structures for communicating data between meshes.
*/
static
void updateComm(const UINT cell, const UINT angle,
                const PsiBoundData &psiBound,
                Mat2<vector<UINT>> &commSidesAngles,
                Mat2<vector<double>> &commPsi)
{
    UINT angleGroup = omp_get_thread_num();
    vector<double> psiSide(g_nVrtxPerFace * g_nGroups);
    for (UINT face = 0; face < g_nFacePerCell; ++face) {
        size_t neighborCell = g_tychoMesh->getAdjCell(cell, face);
        UINT proc = g_tychoMesh->getAdjRank(cell, face);
        
        if (neighborCell == g_tychoMesh->BOUNDARY_FACE && 
            g_tychoMesh->isOutgoing(angle, cell, face) &&
            proc != TychoMesh::BAD_RANK)
        {
            UINT side = g_tychoMesh->getSide(cell, face);
            UINT globalSide = g_tychoMesh->getLGSide(side);
            getBoundData(psiBound, side, angle, psiSide);
            commPsi(angleGroup, proc).insert(
                commPsi(angleGroup, proc).end(), psiSide.begin(), psiSide.end());
            commSidesAngles(angleGroup, proc).push_back(globalSide);
            commSidesAngles(angleGroup, proc).push_back(angle);
        }
    }
}


/*
    updateBoundData
    
    Updates psiBound after doing transport on a set of work.
*/
static
void updateBoundData(const UINT cell, const UINT angle, PsiBoundData &psiBound,
                     const Mat2<double> &localPsi)
{
    for (UINT group = 0; group < g_nGroups; group++) {
    for (UINT face = 0; face < g_nFacePerCell; ++face) {
        
        size_t neighborCell = g_tychoMesh->getAdjCell(cell, face);
        if (neighborCell == g_tychoMesh->BOUNDARY_FACE) {
            
            // if outgoing face
            if (g_tychoMesh->isOutgoing(angle, cell, face)) {
                UINT side = g_tychoMesh->getSide(cell, face);
                for (UINT vertex = 0; vertex < g_nVrtxPerFace; ++vertex) {
                    UINT cellVrtx = 
                        g_tychoMesh->getFaceToCellVrtx(cell, face, vertex);
                    psiBound(group, vertex, angle, side) = 
                        localPsi(cellVrtx, group);
                }
            }
        }
    }}
}


/*
    doComputation
    
    Overall computation part of the sweeper.
    Split out so multiple sweeper schedules can be tested.
*/
static
void doComputation(const UINT step,
                   const UINT angleGroup, 
                   const PsiData &source, 
                   PsiData &psi, 
                   Mat2<vector<UINT>> &commSidesAngles,
                   Mat2<vector<double>> &commPsi,
                   PsiBoundData &psiBound)
{
    Mat2<double> localSource(g_nVrtxPerCell, g_nGroups);
    Mat2<double> localPsi(g_nVrtxPerCell, g_nGroups);
    Mat3<double> localPsiBound(g_nVrtxPerFace, g_nFacePerCell, g_nGroups);            
    
    // Do work
    if (step < g_sweepSchedule[angleGroup]->nSteps()) {
        // get work to solve in this step
        for (const SweepSchedule::Work &work : 
             g_sweepSchedule[angleGroup]->getWork(step)) {
            
            UINT cell = work.getCell();
            UINT angle = work.getAngle();
            
            // Populate localSource
            for (UINT group = 0; group < g_nGroups; group++) {
            for (UINT vrtx = 0; vrtx < g_nVrtxPerCell; vrtx++) {
                localSource(vrtx, group) = source(group, vrtx, angle, cell);
            }}
            
            // Populate localPsiBound
            Transport::populateLocalPsiBound(angle, cell, psi, psiBound, 
                                             localPsiBound);
            
            // Transport solve
            Transport::solve(cell, angle, g_sigmaT[cell], 
                             localPsiBound, localSource, localPsi);
            
            // localPsi -> psi
            for (UINT group = 0; group < g_nGroups; group++) {
            for (UINT vrtx = 0; vrtx < g_nVrtxPerCell; vrtx++) {
                psi(group, vrtx, angle, cell) = localPsi(vrtx, group);
            }}
            
            // Update psiBound and comm variables
            updateBoundData(cell, angle, psiBound, localPsi);
            updateComm(cell, angle, psiBound, commSidesAngles, commPsi);
        }
    }
}


/*
    Constructor

    Splits all the angles into sets of angle groups.
    One angle group per OMP thread.
    Creates an array of SweepSchedules, one entry for each angle group.
*/
Sweeper::Sweeper()
{
    // SweepSchedule for each angle group
    g_sweepSchedule = new SweepSchedule*[g_nAngleGroups];
    

    // Get the angle indices for each angle group
    vector<UINT> angleBdryIndices(g_nAngleGroups + 1);
    angleBdryIndices[0] = 0;
    for (UINT angleGroup = 0; angleGroup < g_nAngleGroups; angleGroup++) {
        UINT numAngles = g_nAngles / g_nAngleGroups;
        if (angleGroup < g_nAngles % g_nAngleGroups)
            numAngles++;
        angleBdryIndices[angleGroup+1] = 
            angleBdryIndices[angleGroup] + numAngles;
    }
    

    // Create a SweepSchedule for each angle group
    for (UINT angleGroup = 0; angleGroup < g_nAngleGroups; angleGroup++) {
        UINT numAngles = angleBdryIndices[angleGroup+1] - 
                            angleBdryIndices[angleGroup];
        vector<UINT> angles(numAngles);
        for (UINT angle = 0; angle < numAngles; angle++) {
            angles[angle] = angle + angleBdryIndices[angleGroup];
        }
        g_sweepSchedule[angleGroup] = 
            new SweepSchedule(angles, g_maxCellsPerStep, g_intraAngleP, 
                              g_interAngleP);
    }
}


/*
    solve
*/
void Sweeper::solve()
{
    Problem::getSource(c_source);
    c_psi.setToValue(0.0);

    if (g_useSourceIteration)
        SourceIteration::fixedPoint(*this, c_psi, c_source);
    else
        SourceIteration::krylov(*this, c_psi, c_source);
}


/*
    sweep
    
    Does an Sn transport sweep.
*/
void Sweeper::sweep(PsiData &psi, const PsiData &source, bool zeroPsiBound)
{
    UNUSED_VARIABLE(zeroPsiBound);


    // Time the sweep
    Timer totalTimer;
    totalTimer.start();
    
    
    // Get max steps for all OpenMP threads
    UINT maxNSteps = g_sweepSchedule[0]->nSteps();
    for(UINT angleGroup = 1; angleGroup < g_nAngleGroups; angleGroup++) {
        maxNSteps = max(maxNSteps, g_sweepSchedule[angleGroup]->nSteps());
    }
    
    
    // Communication variables
    Mat2<vector<UINT>> commSidesAngles(g_nAngleGroups, Comm::numRanks());
    Mat2<vector<double>> commPsi(g_nAngleGroups, Comm::numRanks());
    PsiBoundData psiBound;
    
    
    // Time computation for each thread
    vector<double> computationTimes(g_nAngleGroups);
    computationTimes.assign(g_nAngleGroups, 0.0);
    
    
    // Sweep Type 0
    // Splits computation and communication
    // Contains an expensive omp barrier
    if (g_sweepType == SweepType_OriginalTycho1) {
        
        // Do the sweep
        for (UINT step = 0; step < maxNSteps; ++step) {
            
            // Computation
            #pragma omp parallel
            {
                Timer timer1;
                timer1.start();
                UINT angleGroup = omp_get_thread_num();
                doComputation(step, angleGroup, source, psi, 
                              commSidesAngles, commPsi, psiBound);
                timer1.stop();
                computationTimes[angleGroup] += timer1.wall_clock();
            }
            
            
            // Communication (Non blocking send followed by blocking recv)
            vector<MPI_Request> mpiRequests;
            
            for (UINT angleGroup = 0; angleGroup < g_nAngleGroups; angleGroup++) {    
                send(step, angleGroup, commSidesAngles, commPsi, mpiRequests);
            }
            
            for (UINT angleGroup = 0; angleGroup < g_nAngleGroups; angleGroup++) {
                recv(step, angleGroup, psiBound);
            }
            
            MPI_Waitall(mpiRequests.size(), &mpiRequests[0], MPI_STATUSES_IGNORE);
            //Comm::barrier();
        }
    }
    
    
    // Sweep Type 1
    // Overlaps computation and communication between threads
    // No omp barrier
    if (g_sweepType == SweepType_OriginalTycho2) {
        
        // Thread ID allowed to communicate
        UINT commNumber = 0;
        
        // Do the sweep
        #pragma omp parallel
        {
            UINT angleGroup = omp_get_thread_num();
            
            for (UINT step = 0; step < maxNSteps; ++step) {
                
                // Computation
                Timer timer1;
                timer1.start();
                doComputation(step, angleGroup, source, psi, 
                              commSidesAngles, commPsi, psiBound);
                timer1.stop();
                computationTimes[angleGroup] += timer1.wall_clock();
                
                
                // Require communication in thread id order
                while(true) {
                    double temp;
                    #pragma omp atomic read
                    temp = commNumber;
                    if (temp == angleGroup)
                        break;
                }
                
                
                // Nonblocking send and blocking recv
                vector<MPI_Request> mpiRequests;
                send(step, angleGroup, commSidesAngles, commPsi, mpiRequests);
                recv(step, angleGroup, psiBound);
                MPI_Waitall(mpiRequests.size(), &mpiRequests[0], 
                            MPI_STATUSES_IGNORE);
                
                
                // Barrier all MPI ranks for thread 'commNumber'
                Comm::barrier();
                
                
                // Allow next thread to communicate
                #pragma omp atomic write
                commNumber = (angleGroup + 1) % g_nAngleGroups;
                
            }
        }
    }
    
    
    // Stop timing the sweep
    totalTimer.stop();
    
    
    // Timing output
    double totalTime = totalTimer.wall_clock();
    Comm::gmax(totalTime);
    
    for (UINT i = 0; i < g_nAngleGroups; i++) {
        Comm::gmax(computationTimes[i]);
    }
    
    if (Comm::rank() == 0) {
        double computationTime = 0.0;
        for (UINT i = 0; i < g_nAngleGroups; i++) {
            computationTime = max(computationTime, computationTimes[i]);
        }
        printf("     Sweeper Timer (computation): %fs\n", computationTime);
        printf("     Sweeper Timer (other): %fs\n", totalTime - computationTime);
        printf("     Sweeper Timer (total time): %fs\n", totalTime);
    }
}

