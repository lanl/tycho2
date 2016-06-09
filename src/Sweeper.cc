/*
    Inner.cc
    
    Implements source iteration to solve the transport problem.
*/

#include "Sweeper.hh"
/*
    Sweeper.cc
    
    Implements the Sn transport sweep using the scheduling from
    SweepScheduler.
*/

#include "Quadrature.hh"
#include "Global.hh"
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
void getBoundData(const PsiData &psiBound, const UINT side,
                  const UINT angle, vector<double> &psiSide)
{
    for (UINT vertex = 0; vertex < g_nVrtxPerFace; ++vertex) {
    for (UINT group = 0; group < g_nGroups; ++group) {
        psiSide[VG(vertex, group)] = psiBound(side, angle, vertex, group);
    }}
}


/*
    setBoundData
    
    Sets psiBound from the communication data structures.
*/
static
void setBoundData(PsiData &psiBound, const vector<UINT> &commSides,
                  const vector<UINT> &commAngles, const vector<double> &commPsi) 
{
    for(unsigned i = 0; i < commSides.size(); i++) {
        UINT side = commSides[i];
        UINT angle = commAngles[i];
        UINT commPsiIndex = i * g_nVrtxPerFace * g_nGroups;
        for (UINT vrtx = 0; vrtx < g_nVrtxPerFace; ++vrtx) {
        for (UINT group = 0; group < g_nGroups; ++group) {
            psiBound(side, angle, vrtx, group) = commPsi[commPsiIndex + VG(vrtx, group)];
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
    if (step < g_spSweepSchedule[angleGroup]->nSteps()) {
        for (UINT proc : g_spSweepSchedule[angleGroup]->getSendProcs(step)) {
            
            MPI_Request mpiRequest;
            UINT nSides = commSidesAngles(angleGroup, proc).size()/2;
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
void recv(const UINT step, const UINT angleGroup, PsiData &psiBound)
{
    if (step < g_spSweepSchedule[angleGroup]->nSteps()) {
        for (UINT proc : g_spSweepSchedule[angleGroup]->getRecvProcs(step)) {
            
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
                    commSides[side] = g_spTychoMesh->getGLSide(commSidesAngles[2*side]);
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
                const PsiData &psiBound,
                Mat2<vector<UINT>> &commSidesAngles,
                Mat2<vector<double>> &commPsi)
{
    UINT angleGroup = omp_get_thread_num();
    vector<double> psiSide(g_nVrtxPerFace * g_nGroups);
    for (UINT face = 0; face < g_nFacePerCell; ++face) {
        size_t neighborCell = g_spTychoMesh->getAdjCell(cell, face);
        UINT proc = g_spTychoMesh->getAdjRank(cell, face);
        
        if (neighborCell == g_spTychoMesh->BOUNDARY_FACE && 
            g_spTychoMesh->isOutgoing(angle, cell, face) &&
            proc != TychoMesh::BAD_RANK)
        {
            UINT side = g_spTychoMesh->getSide(cell, face);
            UINT globalSide = g_spTychoMesh->getLGSide(side);
            getBoundData(psiBound, side, angle, psiSide);
            commPsi(angleGroup, proc).insert(commPsi(angleGroup, proc).end(), psiSide.begin(), psiSide.end());
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
void updateBoundData(const UINT cell, const UINT angle, PsiData &psiBound,
                     const Mat2<double> &localPsi)
{
    for (UINT group = 0; group < g_nGroups; group++) {
    for (UINT face = 0; face < g_nFacePerCell; ++face) {
        
        size_t neighborCell = g_spTychoMesh->getAdjCell(cell, face);
        if (neighborCell == g_spTychoMesh->BOUNDARY_FACE) {
            
            // if outgoing face
            if (g_spTychoMesh->isOutgoing(angle, cell, face)) {
                UINT side = g_spTychoMesh->getSide(cell, face);
                for (UINT vertex = 0; vertex < g_nVrtxPerFace; ++vertex) {
                    psiBound(side, angle, vertex, group) =
                        localPsi(g_spTychoMesh->getFaceToCellVrtx(cell, face, vertex), group);
                }
            }
        }
    }}
}


/*
    populateLocalPsiBound
    
    Put data from neighboring cells into localPsiBound(fvrtx, face, group).
*/
static
void populateLocalPsiBound(const UINT angle, const UINT cell, 
                           const PsiData &psi, const PsiData &psiBound,
                           Mat3<double> &localPsiBound)
{
    // Default to 0.0
    localPsiBound = 0.0;
    
    // Populate if incoming flux
    for (UINT group = 0; group < g_nGroups; group++) {
    for (UINT face = 0; face < g_nFacePerCell; face++) {
        if (g_spTychoMesh->isIncoming(angle, cell, face)) {
            size_t neighborCell = g_spTychoMesh->getAdjCell(cell, face);
            
            // In local mesh
            if (neighborCell != TychoMesh::BOUNDARY_FACE) {
                for (UINT fvrtx = 0; fvrtx < g_nVrtxPerFace; fvrtx++) {
                    UINT neighborVrtx = 
                        g_spTychoMesh->getNeighborVrtx(cell, face, fvrtx);
                    localPsiBound(fvrtx, face, group) = 
                        psi(neighborVrtx, angle, neighborCell, group);
                }
            }
            
            // Not in local mesh
            else if (g_spTychoMesh->getAdjRank(cell, face) != TychoMesh::BAD_RANK) {
                for (UINT fvrtx = 0; fvrtx < g_nVrtxPerFace; fvrtx++) {
                    UINT side = g_spTychoMesh->getSide(cell, face);
                    localPsiBound(fvrtx, face, group) = 
                        psiBound(side, angle, fvrtx, group);
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
                   const double sigmaTotal,
                   const PsiData &source, 
                   PsiData &psi, 
                   Mat2<vector<UINT>> &commSidesAngles,
                   Mat2<vector<double>> &commPsi,
                   PsiData &psiBound)
{
    Mat2<double> localSource(g_nVrtxPerCell, g_nGroups);
    Mat2<double> localPsi(g_nVrtxPerCell, g_nGroups);
    Mat3<double> localPsiBound(g_nVrtxPerFace, g_nFacePerCell, g_nGroups);            
    
    // Do work
    if (step < g_spSweepSchedule[angleGroup]->nSteps()) {
        // get work to solve in this step
        for (const SweepSchedule::Work &work : 
             g_spSweepSchedule[angleGroup]->getWork(step)) {
            
            UINT cell = work.getCell();
            UINT angle = work.getAngle();
            
            // Populate localSource
            for (UINT group = 0; group < g_nGroups; group++) {
            for (UINT vrtx = 0; vrtx < g_nVrtxPerCell; vrtx++) {
                localSource(vrtx, group) = source(vrtx, angle, cell, group);
            }}
            
            // Populate localPsiBound
            populateLocalPsiBound(angle, cell, psi, psiBound, localPsiBound);
            
            // Transport solve
            Transport::solve(cell, angle, sigmaTotal, 
                             localPsiBound, localSource, localPsi);
            
            // localPsi -> psi
            for (UINT group = 0; group < g_nGroups; group++) {
            for (UINT vrtx = 0; vrtx < g_nVrtxPerCell; vrtx++) {
                psi(vrtx, angle, cell, group) = localPsi(vrtx, group);
            }}
            
            // Update psiBound and comm variables
            updateBoundData(cell, angle, psiBound, localPsi);
            updateComm(cell, angle, psiBound, commSidesAngles, commPsi);
        }
    }
}


namespace Sweeper {

/*
    sweep
    
    Does an Sn transport sweep.
*/
void sweep(PsiData &psi, const PsiData &source, const double sigmaTotal)
{
    // Time the sweep
    Timer totalTimer;
    totalTimer.start();
    
    
    // Get max steps for all OpenMP threads
    UINT maxNSteps = g_spSweepSchedule[0]->nSteps();
    for(UINT angleGroup = 1; angleGroup < g_nAngleGroups; angleGroup++) {
        maxNSteps = max(maxNSteps, g_spSweepSchedule[angleGroup]->nSteps());
    }
    
    
    // Communication variables
    Mat2<vector<UINT>> commSidesAngles(g_nAngleGroups, Comm::numRanks());
    Mat2<vector<double>> commPsi(g_nAngleGroups, Comm::numRanks());
    PsiData psiBound(g_spTychoMesh->getNSides(), g_quadrature->getNumAngles(), 
                         g_nVrtxPerFace, g_nGroups);
    
    
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
                doComputation(step, angleGroup, sigmaTotal, source, psi, 
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
                doComputation(step, angleGroup, sigmaTotal, source, psi, 
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

} // End namespace

