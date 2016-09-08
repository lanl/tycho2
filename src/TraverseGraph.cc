/*
    TraverseGraph.cc
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

#include "TraverseGraph.hh"
#include "Mat.hh"
#include "Global.hh"
#include "Comm.hh"
#include "Timer.hh"
#include <vector>
#include <set>
#include <queue>
#include <utility>
#include <omp.h>
#include <limits.h>
#include <string.h>

using namespace std;

    Timer timer;

/*
    Tuple class
*/
class Tuple
{
private:
    UINT c_cell;
    UINT c_angle;
    UINT c_priority;
    
public:
    Tuple(UINT cell, UINT angle, UINT priority)
        : c_cell(cell), c_angle(angle), c_priority(priority) {}
    
    UINT getCell() const { return c_cell; }
    UINT getAngle() const { return c_angle; }
    
    // Comparison operator to determine relative priorities
    // Needed for priority_queue
    bool operator<(const Tuple &rhs) const
    {
        return c_priority < rhs.c_priority;
    }
};


/*
    splitPacket
    
    Packet is (global side, angle, data)
*/
static
void splitPacket(char **packet, UINT &globalSide, UINT &angle)
{
    memcpy(&globalSide, *packet, sizeof(UINT));
    *packet += sizeof(UINT);
    memcpy(&angle, *packet, sizeof(UINT));
    *packet += sizeof(UINT);
}


/*
    createPacket
    
    Packet is (global side, angle, data)
*/
static
void createPacket(vector<char> &packet, UINT globalSide, UINT angle, 
                  UINT dataSize, const char *data)
{
    packet.resize(2 * sizeof(UINT) + dataSize);
    char *p = packet.data();
    
    memcpy(p, &globalSide, sizeof(UINT));
    p += sizeof(UINT);
    memcpy(p, &angle, sizeof(UINT));
    p += sizeof(UINT);
    memcpy(p, data, dataSize);
}


/*
    isIncoming
    
    Determines whether data is incoming to the cell 
    depending on sweep direction.
*/
static
bool isIncoming(UINT angle, UINT cell, UINT face, Direction direction)
{
    if (direction == Direction_Forward)
        return g_tychoMesh->isIncoming(angle, cell, face);
    else if (direction == Direction_Backward)
        return g_tychoMesh->isOutgoing(angle, cell, face);
    
    // Should never get here
    Assert(false);
    return false;
}


/*
    angleGroupIndex
    
    Gets angle groups for angle index.
    e.g. 20 angles numbered 0...19 with 3 threads.
    Split into 3 angle chunks of size 7,7,6:  0...6  7...13  14...19
    If angle in 0...6,   return 0
    If angle in 7...13,  return 1
    If angle in 14...19, return 2
*/
UINT angleGroupIndex(UINT angle)
{
    UINT numAngles = g_quadrature->getNumAngles();
    UINT chunkSize = numAngles / g_nThreads;
    UINT numChunksBigger = numAngles % g_nThreads;
    UINT lowIndex = 0;
    
    
    // Find angleGroup
    for (UINT angleGroup = 0; angleGroup < g_nThreads; angleGroup++) {
        
        UINT nextLowIndex = lowIndex + chunkSize;
        if (angleGroup < numChunksBigger)
            nextLowIndex++;
        
        if (angle < nextLowIndex) 
            return angleGroup;
        
        lowIndex = nextLowIndex;
    }
    
    
    // Should never get here
    Assert(false);
    return 0;
}


/*
    sendAndRecvData()
    
    The algorithm is
    - Irecv data size for all adjacent ranks
    - ISend data size and then data (if any)
    - Wait on recv of data size, then blocking recv for data if there is any.
    
    Data is sent in two steps to each adjacent rank.
    First is the number of bytes of data that will be sent.
    Second is the raw data in bytes.
    The tag for the first send is 0.
    The tag for the second send is 1.
    
    The raw data is made of data packets containing:
    globalSide, angle, and data to send.
    The data can have different meanings depending on the TraverseData 
    subclass.
    
    When done traversing local graph, you want to stop all communication.
    This is done by setting killComm to true.
    To mark killing communication, sendSize is set to UINT64_MAX.
    In this event, commDark[rank] is set to true on the receiving rank 
    so we no longer look for communication from this rank.
*/
static
void sendAndRecvData(const vector<vector<char>> &sendBuffers, 
                     const vector<UINT> &adjRankIndexToRank, 
                     TraverseData &traverseData, 
                     set<pair<UINT,UINT>> &sideRecv,
                     vector<bool> &commDark, const bool killComm,
                     MPI_Comm mpiComm)
{
    // Check input
    Assert(adjRankIndexToRank.size() == sendBuffers.size());
    Assert(adjRankIndexToRank.size() == commDark.size());
    
    
    // Variables
    UINT numAdjRanks = adjRankIndexToRank.size();
    int mpiError;
    
    vector<UINT> recvSizes(numAdjRanks);
    vector<UINT> sendSizes(numAdjRanks);
    vector<MPI_Request> mpiRecvRequests(numAdjRanks);
    vector<MPI_Request> mpiSendRequests;
    UINT numRecv = numAdjRanks;
    
    
    // Setup recv of data size
    for (UINT index = 0; index < numAdjRanks; index++) {
        
        // No recv if adjRank is no longer communicating
        if (commDark[index]) {
            mpiRecvRequests[index] = MPI_REQUEST_NULL;
            numRecv--;
            continue;
        }
        
        // Irecv data size
        int numDataToRecv = 1;
        int adjRank = adjRankIndexToRank[index];
        int tag0 = 0;
        
        mpiError = MPI_Irecv(&recvSizes[index], numDataToRecv, MPI_UINT64_T, 
                             adjRank, tag0, mpiComm, &mpiRecvRequests[index]);
        Insist(mpiError == MPI_SUCCESS, "");
    }
    
    
    // Send data size and data
    for (UINT index = 0; index < numAdjRanks; index++) {
        
        // Don't send if adjRank is no longer communicating
        if (commDark[index])
            continue;
        
        const vector<char> &sendBuffer = sendBuffers[index];
        int numDataToSend = 1;
        int adjRank = adjRankIndexToRank[index];
        int tag0 = 0;
        int tag1 = 1;
        
        
        // Send data size
        MPI_Request request;
        sendSizes[index] = sendBuffer.size();
        if (killComm)
            sendSizes[index] = UINT64_MAX;
        
        mpiError = MPI_Isend(&sendSizes[index], numDataToSend, MPI_UINT64_T, 
                             adjRank, tag0, mpiComm, &request);
        Insist(mpiError == MPI_SUCCESS, "");
        mpiSendRequests.push_back(request);
        
        
        // Send data
        if (sendSizes[index] > 0 && sendSizes[index] != UINT64_MAX) {
            MPI_Request request;
            Assert(sendBuffer.size() < INT_MAX);
            
            mpiError = MPI_Isend(const_cast<char*>(sendBuffer.data()), 
                                 sendBuffer.size(), 
                                 MPI_BYTE, adjRank, tag1, 
                                 mpiComm, &request);
            Insist(mpiError == MPI_SUCCESS, "");
            mpiSendRequests.push_back(request);
        }
    }
    
    
    // Recv data size and data
    for (UINT numWaits = 0; numWaits < numRecv; numWaits++) {
        
        // Wait for a data size to arrive
        int index;

		timer.start();
        mpiError = MPI_Waitany(mpiRecvRequests.size(), mpiRecvRequests.data(), 
                               &index, MPI_STATUS_IGNORE);
        Insist(mpiError == MPI_SUCCESS, "");
        timer.stop();
        
        // Recv data
        if (recvSizes[index] > 0 && recvSizes[index] != UINT64_MAX) {
            
            int adjRank = adjRankIndexToRank[index];
            int tag1 = 1;
            vector<char> dataPackets(recvSizes[index]);
            
            mpiError = MPI_Recv(dataPackets.data(), recvSizes[index], 
                                MPI_BYTE, adjRank, tag1, mpiComm, 
                                MPI_STATUS_IGNORE);
            Insist(mpiError == MPI_SUCCESS, "");
            
            UINT packetSize = 2 * sizeof(UINT) + traverseData.getDataSize();
            UINT numPackets = recvSizes[index] / packetSize;
            Assert(recvSizes[index] % packetSize == 0);
            
            for (UINT i = 0; i < numPackets; i++) {
                char *packet = &dataPackets[i * packetSize];
                UINT globalSide;
                UINT angle;
                splitPacket(&packet, globalSide, angle);
                
                UINT localSide = g_tychoMesh->getGLSide(globalSide);
                traverseData.setSideData(localSide, angle, packet);
                sideRecv.insert(make_pair(localSide,angle));
            }
        }
        
        
        // Stop communication with this rank
        if (recvSizes[index] == UINT64_MAX) {
            commDark[index] = true;
        }
    }
    
    
    // Make sure all sends are done
    if (mpiSendRequests.size() > 0) {
        mpiError = MPI_Waitall(mpiSendRequests.size(), mpiSendRequests.data(), 
                               MPI_STATUSES_IGNORE);
        Insist(mpiError == MPI_SUCCESS, "");
    }

}


/*
    traverseGraph
    
    Traverses g_tychoMesh.
    If doComm is true, this is done globally.
    If doComm is false, each mesh partition is traversed locally with no
    consideration for boundaries between partitions.
*/
void traverseGraph(const UINT maxComputePerStep,
                   TraverseData &traverseData, bool doComm,
                   MPI_Comm mpiComm, Direction direction)
{
    UINT numCells = g_tychoMesh->getNCells();
    UINT numAngles = g_quadrature->getNumAngles();
    vector<priority_queue<Tuple>> canCompute(g_nThreads);
    Mat2<UINT> numDependencies(numCells, numAngles);
    UINT numCellAnglePairsToCalculate = numAngles * numCells;
    vector<UINT> adjRankIndexToRank;
    map<UINT,UINT> adjRankToRankIndex;
    set<pair<UINT,UINT>> sideRecv;
    Mat2<vector<char>> sendBuffers;
    vector<vector<char>> sendBuffers1;
    vector<bool> commDark;
    Timer totalTimer;
    Timer commTimer;
    Timer updateTimer;
    
    // Start total timer
    totalTimer.start();
    
    
    // Calc num dependencies for each (cell, angle) pair
    for (UINT angle = 0; angle < numAngles; angle++) {
    for (UINT cell = 0; cell < g_tychoMesh->getNCells(); cell++) {
        
        numDependencies(cell, angle) = 0;
        for (UINT face = 0; face < g_nFacePerCell; face++) {
            
            bool incoming = isIncoming(angle, cell, face, direction);
            UINT adjRank = g_tychoMesh->getAdjRank(cell, face);
            UINT adjCell = g_tychoMesh->getAdjCell(cell, face);
            
            if (doComm && incoming && adjRank != TychoMesh::BAD_RANK) {
                numDependencies(cell, angle)++;
            }
            else if (!doComm && incoming && adjCell != TychoMesh::BOUNDARY_FACE) {
                numDependencies(cell, angle)++;
            }
        }
    }}
    
    
    // Get adjacent ranks
    for (UINT cell = 0; cell < g_tychoMesh->getNCells(); cell++) {
    for (UINT face = 0; face < g_nFacePerCell; face++) {
        
        UINT adjRank = g_tychoMesh->getAdjRank(cell, face);
        UINT adjCell = g_tychoMesh->getAdjCell(cell, face);
        
        if (adjCell == TychoMesh::BOUNDARY_FACE && 
            adjRank != TychoMesh::BAD_RANK &&
            adjRankToRankIndex.count(adjRank) == 0)
        {
            UINT rankIndex = adjRankIndexToRank.size();
            adjRankToRankIndex.insert(make_pair(adjRank, rankIndex));
            adjRankIndexToRank.push_back(adjRank);
        }
    }}
    
    
    // Set size of sendBuffers and commDark
    UINT numAdjRanks = adjRankIndexToRank.size();
    sendBuffers = Mat2<vector<char>>(g_nThreads, numAdjRanks);
    sendBuffers1.resize(numAdjRanks);
    commDark.resize(numAdjRanks, false);
    
    
    // Initialize canCompute queue
    for (UINT angle = 0; angle < numAngles; angle++) {
    for (UINT cell = 0; cell < g_tychoMesh->getNCells(); cell++) {
        if (numDependencies(cell, angle) == 0) {
            UINT priority = traverseData.getPriority(cell, angle);
            canCompute[angleGroupIndex(angle)].push(Tuple(cell, angle, priority));
        }
    }}
    
    
    // Traverse the graph
    while (numCellAnglePairsToCalculate > 0) {
        
        // Do local traversal
        #pragma omp parallel
        {
            UINT stepsTaken = 0;
            UINT angleGroup = omp_get_thread_num();
            while (canCompute[angleGroup].size() > 0 && 
                   stepsTaken < maxComputePerStep)
            {
                // Get cell/angle pair to compute
                Tuple cellAnglePair = canCompute[angleGroup].top();
                canCompute[angleGroup].pop();
                UINT cell = cellAnglePair.getCell();
                UINT angle = cellAnglePair.getAngle();
                #pragma omp atomic
                numCellAnglePairsToCalculate--;
                stepsTaken++;
                
                
                // Get boundary type and adjacent cell/side data for each face
                BoundaryType bdryType[g_nFacePerCell];
                UINT adjCellsSides[g_nFacePerCell];
                for (UINT face = 0; face < g_nFacePerCell; face++) {
                    
                    UINT adjCell = g_tychoMesh->getAdjCell(cell, face);
                    UINT adjRank = g_tychoMesh->getAdjRank(cell, face);
                    adjCellsSides[face] = adjCell;
                    
                    if (g_tychoMesh->isOutgoing(angle, cell, face)) {
                        
                        if (adjCell == TychoMesh::BOUNDARY_FACE && 
                            adjRank != TychoMesh::BAD_RANK)
                        {
                            bdryType[face] = BoundaryType_OutIntBdry;
                            adjCellsSides[face] = g_tychoMesh->getSide(cell, face);
                        }
                        
                        else if (adjCell == TychoMesh::BOUNDARY_FACE && 
                                 adjRank == TychoMesh::BAD_RANK)
                        {
                            bdryType[face] = BoundaryType_OutExtBdry;
                        }
                        
                        else {
                            bdryType[face] = BoundaryType_OutInt;
                        }
                    }
                    else {
                        
                        if (adjCell == TychoMesh::BOUNDARY_FACE && 
                            adjRank != TychoMesh::BAD_RANK)
                        {
                            bdryType[face] = BoundaryType_InIntBdry;
                            adjCellsSides[face] = g_tychoMesh->getSide(cell, face);
                        }
                        
                        else if (adjCell == TychoMesh::BOUNDARY_FACE && 
                                 adjRank == TychoMesh::BAD_RANK)
                        {
                            bdryType[face] = BoundaryType_InExtBdry;
                        }
                        
                        else {
                            bdryType[face] = BoundaryType_InInt;
                        }
                    }
                }
                
                
                // Update data for this cell-angle pair
                #pragma omp master
				updateTimer.start();
				traverseData.update(cell, angle, adjCellsSides, bdryType);
                #pragma omp master
				updateTimer.stop();
                
                // Update dependency for children
                for (UINT face = 0; face < g_nFacePerCell; face++) {
                    
                    if (!isIncoming(angle, cell, face, direction)) {
                        
                        UINT adjCell = g_tychoMesh->getAdjCell(cell, face);
                        UINT adjRank = g_tychoMesh->getAdjRank(cell, face);
                        
                        if (adjCell != TychoMesh::BOUNDARY_FACE) {
                            numDependencies(adjCell, angle)--;
                            if (numDependencies(adjCell, angle) == 0) {
                                UINT priority = traverseData.getPriority(adjCell, angle);
                                canCompute[angleGroup].push(Tuple(adjCell, angle, priority));
                            }
                        }
                        
                        else if (doComm && adjRank != TychoMesh::BAD_RANK) {
                            UINT rankIndex = adjRankToRankIndex.at(adjRank);
                            UINT side = g_tychoMesh->getSide(cell, face);
                            UINT globalSide = g_tychoMesh->getLGSide(side);
                            
                            vector<char> packet;
                            createPacket(packet, globalSide, angle, 
                                         traverseData.getDataSize(), 
                                         traverseData.getData(cell, face, angle));
                            
                            sendBuffers(angleGroup, rankIndex).insert(
                                sendBuffers(angleGroup, rankIndex).end(), 
                                packet.begin(), packet.end());
                        }
                    }
                }
            }
        }
        
        
        // Put together sendBuffers from different angleGroups
        vector<vector<char>> sendBuffers1(numAdjRanks);
        for (UINT angleGroup = 0; angleGroup < g_nThreads; angleGroup++) {
        for (UINT rankIndex = 0; rankIndex < numAdjRanks; rankIndex++) {
            sendBuffers1[rankIndex].insert(
                sendBuffers1[rankIndex].end(), 
                sendBuffers(angleGroup, rankIndex).begin(), 
                sendBuffers(angleGroup, rankIndex).end());
        }}
               
        
        // Do communication
        commTimer.start();
        if (doComm) {
            
            // Send/Recv
            sideRecv.clear();
            
            const bool killComm = false;
            sendAndRecvData(sendBuffers1, adjRankIndexToRank, traverseData, 
                            sideRecv, commDark, killComm, mpiComm);
            
            
            // Clear send buffers for next iteration
            for (UINT angleGroup = 0; angleGroup < g_nThreads; angleGroup++) {
            for (UINT rankIndex = 0; rankIndex < numAdjRanks; rankIndex++) {
                sendBuffers(angleGroup, rankIndex).clear();
            }}
            
            for (UINT rankIndex = 0; rankIndex < numAdjRanks; rankIndex++) {
                sendBuffers1[rankIndex].clear();
            }
            
            
            // Update dependency for parents using received side data
            for (auto sideAngle : sideRecv) {
                UINT side = sideAngle.first;
                UINT angle = sideAngle.second;
                UINT cell = g_tychoMesh->getSideCell(side);
                numDependencies(cell, angle)--;
                if (numDependencies(cell, angle) == 0) {
                    UINT priority = traverseData.getPriority(cell, angle);
                    canCompute[angleGroupIndex(angle)].push(Tuple(cell, angle, priority));
                }
            }
        }
        commTimer.stop();
    }
    
    
    // Send kill comm signal to adjacent ranks
    commTimer.start();
    if (doComm) {
        const bool killComm = true;
        sendAndRecvData(sendBuffers1, adjRankIndexToRank, traverseData, 
                        sideRecv, commDark, killComm, mpiComm);
    }
    commTimer.stop();
    
    
    // Print times
    totalTimer.stop();
    double totalTime = totalTimer.wall_clock();
    Comm::gmax(totalTime, mpiComm);
    double commTime = commTimer.sum_wall_clock();
    Comm::gmax(commTime, mpiComm);
    
	double time = timer.sum_wall_clock();
	Comm::gmax(time, mpiComm);

    double updateTime = updateTimer.sum_wall_clock();
    Comm::gmax(updateTime, mpiComm);
    if (Comm::rank(mpiComm) == 0) {
        printf("     Traverse Timer (communication): %fs\n", commTime);
        printf("     Traverse Timer (total time): %fs\n", totalTime);
    }

	if (Comm::rank(mpiComm) == 0) {
		printf("	commdiff %fs\n", time);
	}
	timer.reset();

	if (Comm::rank(mpiComm) == 0) {
		printf("	Traverse Timer update  %fs\n", updateTime);
	}
	
}



