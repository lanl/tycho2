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

#include "GraphTraverser.hh"
#include "Mat.hh"
#include "Global.hh"
#include "TychoMesh.hh"
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


/*
    Tuple class
*/
namespace {

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
};}


/*
    splitPacket
    
    Packet is (global side, angle, data)
*/
static
void splitPacket(char *packet, UINT &globalSide, UINT &angle, char **data)
{
    memcpy(&globalSide, packet, sizeof(UINT));
    packet += sizeof(UINT);
    memcpy(&angle, packet, sizeof(UINT));
    packet += sizeof(UINT);
    *data = packet;
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
    UINT numAngles = g_nAngles;
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
    OneSidedImpl implementation
*/
GraphTraverser::OneSidedImpl::
OneSidedImpl(UINT numAdjRanks,
             std::vector<UINT> offRankHeaderOffsets,
             std::vector<UINT> offRankDataOffsets,
             std::vector<UINT> onRankHeaderOffsets,
             std::vector<UINT> onRankDataOffsets,
             UINT dataSizePerChunk,
             MPI_Win mpiWin) : 
        c_offRankHeaderOffsets(offRankHeaderOffsets),
        c_offRankDataOffsets(offRankDataOffsets),
        c_onRankHeaderOffsets(onRankHeaderOffsets),
        c_onRankDataOffsets(onRankDataOffsets),
        c_dataSizePerChunk(dataSizePerChunk),
        c_mpiWin(mpiWin)
{
    c_send_numWrittenVector[0].resize(numAdjRanks);
    c_send_numWrittenVector[1].resize(numAdjRanks);
    c_send_currentDataChunkVector.resize(numAdjRanks);

    c_recv_numReadVector[0].resize(numAdjRanks);
    c_recv_numReadVector[1].resize(numAdjRanks);
    c_recv_headerDataVector.resize(4 * numAdjRanks);
    c_recv_currentDataChunkVector.resize(numAdjRanks);
}


/*
    Send specific functions
*/
uint32_t GraphTraverser::OneSidedImpl::
send_getNumWritten(int adjRankIndex)
{
    uint32_t dataChunk = c_send_currentDataChunkVector[adjRankIndex];
    return c_send_numWrittenVector[dataChunk][adjRankIndex];
}

void GraphTraverser::OneSidedImpl::
send_setNumWritten(int adjRankIndex, uint32_t numBytesWritten)
{
    uint32_t dataChunk = c_send_currentDataChunkVector[adjRankIndex];
    c_send_numWrittenVector[dataChunk][adjRankIndex] = numBytesWritten;
}

void GraphTraverser::OneSidedImpl::
send_switchDataChunk(int adjRankIndex)
{
    uint32_t dataChunk = c_send_currentDataChunkVector[adjRankIndex];
    dataChunk = (dataChunk + 1) % 2;
    c_send_currentDataChunkVector[adjRankIndex] = dataChunk;
}

UINT GraphTraverser::OneSidedImpl::
send_getLockOffset(int adjRankIndex)
{
    uint32_t dataChunk = c_send_currentDataChunkVector[adjRankIndex];
    return c_offRankHeaderOffsets[adjRankIndex] + 
           sizeof(uint32_t) * dataChunk;
}

UINT GraphTraverser::OneSidedImpl::
send_getNumWrittenOffset(int adjRankIndex)
{
    uint32_t dataChunk = c_send_currentDataChunkVector[adjRankIndex];
    return c_offRankHeaderOffsets[adjRankIndex] + 
           sizeof(uint32_t) * (2 + dataChunk);
}

UINT GraphTraverser::OneSidedImpl::
send_getDataOffset(int adjRankIndex, uint32_t bytesOffset)
{
    uint32_t dataChunk = c_send_currentDataChunkVector[adjRankIndex];
    return c_offRankDataOffsets[adjRankIndex] + 
           c_dataSizePerChunk * dataChunk + bytesOffset;
}


/*
    Recv specific functions
*/
uint32_t GraphTraverser::OneSidedImpl::
recv_getNumRead(int adjRankIndex)
{
    uint32_t dataChunk = c_recv_currentDataChunkVector[adjRankIndex];
    return c_recv_numReadVector[dataChunk][adjRankIndex];
}

void GraphTraverser::OneSidedImpl::
recv_setNumRead(int adjRankIndex, uint32_t numBytesRead)
{
    uint32_t dataChunk = c_recv_currentDataChunkVector[adjRankIndex];
    c_recv_numReadVector[dataChunk][adjRankIndex] = numBytesRead;
}

void GraphTraverser::OneSidedImpl::
recv_switchDataChunk(int adjRankIndex)
{
    uint32_t dataChunk = c_recv_currentDataChunkVector[adjRankIndex];
    dataChunk = (dataChunk + 1) % 2;
    c_recv_currentDataChunkVector[adjRankIndex] = dataChunk;
}

UINT GraphTraverser::OneSidedImpl::
recv_getBaseHeaderOffset()
{
    return c_onRankHeaderOffsets[0];
}

int GraphTraverser::OneSidedImpl::
recv_getHeaderCount()
{
    return c_recv_headerDataVector.size();
}

uint32_t* GraphTraverser::OneSidedImpl::
recv_getHeaderPointer()
{
    return c_recv_headerDataVector.data();
}

UINT GraphTraverser::OneSidedImpl::
recv_getLockOffset(int adjRankIndex)
{
    uint32_t dataChunk = c_recv_currentDataChunkVector[adjRankIndex];
    return c_onRankHeaderOffsets[adjRankIndex] + 
           sizeof(uint32_t) * dataChunk;
}

UINT GraphTraverser::OneSidedImpl::
recv_getNumWrittenOffset(int adjRankIndex)
{
    uint32_t dataChunk = c_recv_currentDataChunkVector[adjRankIndex];
    return c_onRankHeaderOffsets[adjRankIndex] + 
           sizeof(uint32_t) * (2 +  dataChunk);
}

UINT GraphTraverser::OneSidedImpl::
recv_getDataOffset(int adjRankIndex, uint32_t bytesOffset)
{
    uint32_t dataChunk = c_recv_currentDataChunkVector[adjRankIndex];
    return c_onRankDataOffsets[adjRankIndex] + 
           c_dataSizePerChunk * dataChunk + bytesOffset;
}

uint32_t GraphTraverser::OneSidedImpl::
recv_getLock(int adjRankIndex)
{
    uint32_t dataChunk = c_recv_currentDataChunkVector[adjRankIndex];
    return c_recv_headerDataVector[4 * adjRankIndex + dataChunk];
}

uint32_t GraphTraverser::OneSidedImpl::
recv_getNumWritten(int adjRankIndex)
{
    uint32_t dataChunk = c_recv_currentDataChunkVector[adjRankIndex];
    return c_recv_headerDataVector[4 * adjRankIndex + 2 + dataChunk];
}


/*
    sendData1Sided

    Implements one-sided MPI for sending data.
*/
void GraphTraverser::
sendData1Sided(const vector<vector<char>> &sendBuffers) const
{
    // Send data to each adjacent rank
    int numAdjRanks = c_adjRankIndexToRank.size();
    for (int index = 0; index < numAdjRanks; index++) {
        
        // Make sure there is data to send
        if (sendBuffers[index].size() == 0)
            continue;
        

        // Useful values
        int mpiError;
        UINT offset;
        MPI_Win mpiWin = c_oneSidedImpl->getWin();
        int adjRank = c_adjRankIndexToRank[index];
        const vector<char> &sendBuffer = sendBuffers[index];
        uint32_t numBytesWritten = 
            c_oneSidedImpl->send_getNumWritten(index);
        
        
        // Check to see if there is NOT room to write data in current data chunk
        UINT dataSizePerChunk = c_oneSidedImpl->getDataSizePerChunk();
        if (numBytesWritten + sendBuffer.size() >= dataSizePerChunk) {
            
            // Flag indicating we can't write to this data chunk anymore
            // Must be removed by recvData
            uint32_t one = 1;
            offset = c_oneSidedImpl->send_getLockOffset(index);
            mpiError = MPI_Accumulate(&one, 1, MPI_UINT32_T, adjRank, offset, 
                                      1, MPI_UINT32_T, MPI_REPLACE, mpiWin);
            Insist(mpiError == MPI_SUCCESS, "");
            
            
            // change data chunk we're working on
            numBytesWritten = 0;
            c_oneSidedImpl->send_setNumWritten(index, numBytesWritten);
            c_oneSidedImpl->send_switchDataChunk(index);


            // Spin wait for availability to write to new data chunk
            while (true) {
                uint32_t locked = 1;
                offset = c_oneSidedImpl->send_getLockOffset(index);
                mpiError = MPI_Fetch_and_op(NULL, &locked, MPI_UINT32_T,
                                            adjRank, offset, MPI_NO_OP, mpiWin);
                Insist(mpiError == MPI_SUCCESS, "");
                
                mpiError = MPI_Win_flush_local(adjRank, mpiWin);
                Insist(mpiError == MPI_SUCCESS, "");

                if (locked == 0)
                    break;
            }
        }


        // Write the data
        offset = c_oneSidedImpl->send_getDataOffset(index, numBytesWritten);
        mpiError = MPI_Put(sendBuffer.data(), sendBuffer.size(), 
                           MPI_BYTE, adjRank, offset, 
                           sendBuffer.size(), MPI_BYTE, mpiWin);
        Insist(mpiError == MPI_SUCCESS, "");
        

        // Flush window to know the put is done on both ends
        mpiError = MPI_Win_flush(adjRank, mpiWin);
        Insist(mpiError == MPI_SUCCESS, "");

        
        // Update number of bytes written
        offset = c_oneSidedImpl->send_getNumWrittenOffset(index);
        numBytesWritten += sendBuffer.size();
        c_oneSidedImpl->send_setNumWritten(index, numBytesWritten);
        
        mpiError = MPI_Accumulate(&numBytesWritten, 1, MPI_UINT32_T, adjRank, 
                                  offset, 1, MPI_UINT32_T, MPI_REPLACE, mpiWin);
        Insist(mpiError == MPI_SUCCESS, "");

        
        // Flush window for local variables
        mpiError = MPI_Win_flush_local(adjRank, mpiWin);
        Insist(mpiError == MPI_SUCCESS, "");
    }
}


/*
    recvData1Sided

    Implements one-sided MPI for receiving data.
*/
void GraphTraverser::recvData1Sided(vector<char> &dataPackets) const
{
    // Useful values
    MPI_Win mpiWin = c_oneSidedImpl->getWin();
    int numAdjRanks = c_adjRankIndexToRank.size();
    int myRank = Comm::rank();
    int mpiError;
    UINT offset;
    int count;


    // Get header data for all adjacent ranks
    offset = c_oneSidedImpl->recv_getBaseHeaderOffset();
    count = c_oneSidedImpl->recv_getHeaderCount();
    mpiError = MPI_Get_accumulate(NULL, count, MPI_UINT32_T, 
                                  c_oneSidedImpl->recv_getHeaderPointer(), 
                                  count, MPI_UINT32_T, myRank, offset, count, 
                                  MPI_UINT32_T, MPI_NO_OP, mpiWin);
    Insist(mpiError == MPI_SUCCESS, "");

    mpiError = MPI_Win_flush_local(myRank, mpiWin);
    Insist(mpiError == MPI_SUCCESS, "");
    
    
    // Get the data
    for (int index = 0; index < numAdjRanks; index++) {
    
        // Useful values
        uint32_t numBytesRead = c_oneSidedImpl->recv_getNumRead(index);
        uint32_t numBytesWritten = c_oneSidedImpl->recv_getNumWritten(index);

        
        // Read data
        uint32_t numBytesToRecv = numBytesWritten - numBytesRead;
        if (numBytesToRecv > 0) {
            
            // Get data
            offset = c_oneSidedImpl->recv_getDataOffset(index, numBytesRead);
            UINT originalSize = dataPackets.size();

            dataPackets.resize(originalSize + numBytesToRecv);
            mpiError = MPI_Get(&dataPackets[originalSize], numBytesToRecv, 
                               MPI_BYTE, myRank, offset, numBytesToRecv, 
                               MPI_BYTE, mpiWin);
            Insist(mpiError == MPI_SUCCESS, "");

            mpiError = MPI_Win_flush_local(myRank, mpiWin);
            Insist(mpiError == MPI_SUCCESS, "");


            // Update numRead
            c_oneSidedImpl->recv_setNumRead(index, numBytesWritten);
        }
    
    
        // Check to see if we need to move to other data chunk
        if (c_oneSidedImpl->recv_getLock(index) == 1) {
            
            // Put zeros in header
            uint32_t zero1 = 0;
            offset = c_oneSidedImpl->recv_getNumWrittenOffset(index);
            mpiError = MPI_Accumulate(&zero1, 1, MPI_UINT32_T, myRank, offset, 
                                      1, MPI_UINT32_T, MPI_REPLACE, mpiWin);
            Insist(mpiError == MPI_SUCCESS, "");

            mpiError = MPI_Win_flush(myRank, mpiWin);
            Insist(mpiError == MPI_SUCCESS, "");
            

            // Set locked to zero
            uint32_t zero2 = 0;
            offset = c_oneSidedImpl->recv_getLockOffset(index);
            mpiError = MPI_Accumulate(&zero2, 1, MPI_UINT32_T, myRank, offset, 
                                      1, MPI_UINT32_T, MPI_REPLACE, mpiWin);
            Insist(mpiError == MPI_SUCCESS, "");
            
            mpiError = MPI_Win_flush_local(myRank, mpiWin);
            Insist(mpiError == MPI_SUCCESS, "");
            

            // Update internal state
            c_oneSidedImpl->recv_setNumRead(index, 0);
            c_oneSidedImpl->recv_switchDataChunk(index);
        }
    }
}


/*
    sendData2Sided

    Implements two-sided MPI for sending data.
*/
void GraphTraverser::sendData2Sided(
    const vector<vector<char>> &sendBuffers) const
{
    vector<MPI_Request> mpiSendRequests;
    UINT numAdjRanks = c_adjRankIndexToRank.size();
    int mpiError;


    // Send the data
    for (UINT index = 0; index < numAdjRanks; index++) {
        
        const vector<char> &sendBuffer = sendBuffers[index];
        if (sendBuffer.size() > 0) {
            MPI_Request request;
            const int adjRank = c_adjRankIndexToRank[index];
            const int tag = 0;
            
            mpiError = MPI_Isend(const_cast<char*>(sendBuffer.data()), 
                                 sendBuffer.size(), 
                                 MPI_BYTE, adjRank, tag, 
                                 MPI_COMM_WORLD, &request);
            Insist(mpiError == MPI_SUCCESS, "");
            mpiSendRequests.push_back(request);
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
    recvData2Sided

    Implements two-sided MPI for receiving data.
*/
void GraphTraverser::recvData2Sided(vector<char> &dataPackets) const
{
    UINT numAdjRanks = c_adjRankIndexToRank.size();
    int mpiError;

    for (UINT index = 0; index < numAdjRanks; index++) {
        
        const int adjRank = c_adjRankIndexToRank[index];
        const int tag = 0;
        int flag = 0;
        MPI_Status mpiStatus;
        int recvCount;

        // Probe for new message
        mpiError = MPI_Iprobe(adjRank, tag, MPI_COMM_WORLD, &flag, &mpiStatus);
        Insist(mpiError == MPI_SUCCESS, "");
        mpiError = MPI_Get_count(&mpiStatus, MPI_BYTE, &recvCount);
        Insist(mpiError == MPI_SUCCESS, "");


        // Recv message if there is one
        if (flag) {
            UINT originalSize = dataPackets.size();
            dataPackets.resize(originalSize + recvCount);
            
            mpiError = MPI_Recv(&dataPackets[originalSize], recvCount, 
                                MPI_BYTE, adjRank, tag, MPI_COMM_WORLD, 
                                MPI_STATUS_IGNORE);
            Insist(mpiError == MPI_SUCCESS, "");
        }
    }
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
void GraphTraverser::sendAndRecvData(const vector<vector<char>> &sendBuffers, 
                                     vector<char> &dataPackets,
                                     vector<bool> &commDark, 
                                     const bool killComm) const
{
    // Check input
    Assert(c_adjRankIndexToRank.size() == sendBuffers.size());
    Assert(c_adjRankIndexToRank.size() == commDark.size());
    
    
    // Variables
    UINT numAdjRanks = c_adjRankIndexToRank.size();
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
        int adjRank = c_adjRankIndexToRank[index];
        int tag0 = 0;
        
        mpiError = MPI_Irecv(&recvSizes[index], numDataToRecv, MPI_UINT64_T, 
                             adjRank, tag0, MPI_COMM_WORLD, 
                             &mpiRecvRequests[index]);
        Insist(mpiError == MPI_SUCCESS, "");
    }
    
    
    // Send data size and data
    for (UINT index = 0; index < numAdjRanks; index++) {
        
        // Don't send if adjRank is no longer communicating
        if (commDark[index])
            continue;
        
        const vector<char> &sendBuffer = sendBuffers[index];
        int numDataToSend = 1;
        int adjRank = c_adjRankIndexToRank[index];
        int tag0 = 0;
        int tag1 = 1;
        
        
        // Send data size
        MPI_Request request;
        sendSizes[index] = sendBuffer.size();
        if (killComm)
            sendSizes[index] = UINT64_MAX;
        
        mpiError = MPI_Isend(&sendSizes[index], numDataToSend, MPI_UINT64_T, 
                             adjRank, tag0, MPI_COMM_WORLD, &request);
        Insist(mpiError == MPI_SUCCESS, "");
        mpiSendRequests.push_back(request);
        
        
        // Send data
        if (sendSizes[index] > 0 && sendSizes[index] != UINT64_MAX) {
            MPI_Request request;
            Assert(sendBuffer.size() < INT_MAX);
            
            mpiError = MPI_Isend(const_cast<char*>(sendBuffer.data()), 
                                 sendBuffer.size(), 
                                 MPI_BYTE, adjRank, tag1, 
                                 MPI_COMM_WORLD, &request);
            Insist(mpiError == MPI_SUCCESS, "");
            mpiSendRequests.push_back(request);
        }
    }
    
    
    // Recv data size and data
    for (UINT numWaits = 0; numWaits < numRecv; numWaits++) {
        
        // Wait for a data size to arrive
        int index;
        mpiError = MPI_Waitany(mpiRecvRequests.size(), mpiRecvRequests.data(), 
                               &index, MPI_STATUS_IGNORE);
        Insist(mpiError == MPI_SUCCESS, "");
        

        // Recv data
        if (recvSizes[index] > 0 && recvSizes[index] != UINT64_MAX) {
            
            int adjRank = c_adjRankIndexToRank[index];
            int tag1 = 1;
            UINT originalSize = dataPackets.size();
            dataPackets.resize(originalSize + recvSizes[index]);
            
            mpiError = MPI_Recv(&dataPackets[originalSize], recvSizes[index], 
                                MPI_BYTE, adjRank, tag1, MPI_COMM_WORLD, 
                                MPI_STATUS_IGNORE);
            Insist(mpiError == MPI_SUCCESS, "");
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
    GraphTraverser
    
    If doComm is true, graph traversal is global.
    If doComm is false, each mesh partition is traversed locally with no
    consideration for boundaries between partitions.
*/
GraphTraverser::GraphTraverser(Direction direction, bool doComm, 
                               UINT dataSizeInBytes)
    : c_direction(direction), c_doComm(doComm), 
      c_dataSizeInBytes(dataSizeInBytes)
{
    // Get adjacent ranks
    for (UINT cell = 0; cell < g_nCells; cell++) {
    for (UINT face = 0; face < g_nFacePerCell; face++) {
        
        UINT adjRank = g_tychoMesh->getAdjRank(cell, face);
        UINT adjCell = g_tychoMesh->getAdjCell(cell, face);
        
        if (adjCell == TychoMesh::BOUNDARY_FACE && 
            adjRank != TychoMesh::BAD_RANK &&
            c_adjRankToRankIndex.count(adjRank) == 0)
        {
            UINT rankIndex = c_adjRankIndexToRank.size();
            c_adjRankToRankIndex.insert(make_pair(adjRank, rankIndex));
            c_adjRankIndexToRank.push_back(adjRank);
        }
    }}
    
    
    // Calc num dependencies for each (cell, angle) pair
    c_initNumDependencies.resize(g_nAngles, g_nCells);
    for (UINT cell = 0; cell < g_nCells; cell++) {
    for (UINT angle = 0; angle < g_nAngles; angle++) {
        
        c_initNumDependencies(angle, cell) = 0;
        for (UINT face = 0; face < g_nFacePerCell; face++) {
            
            bool incoming = isIncoming(angle, cell, face, c_direction);
            UINT adjRank = g_tychoMesh->getAdjRank(cell, face);
            UINT adjCell = g_tychoMesh->getAdjCell(cell, face);
            
            if (c_doComm && incoming && adjRank != TychoMesh::BAD_RANK) {
                c_initNumDependencies(angle, cell)++;
            }
            else if (!c_doComm && incoming && adjCell != TychoMesh::BOUNDARY_FACE) {
                c_initNumDependencies(angle, cell)++;
            }
        }
    }}


    // Setup one-sided MPI
    if (g_mpiType == MPIType_OneSided) {
        setupOneSidedMPI();
    }
}


/*
    setupOneSidedMPI

    Memory is laid out as:
    (lock1,lock2) * numAdjRanks
    (numWritten1,numWritten2,...,numWrittenN,writable1,writable2,...,writableN) * numAdjRanks
    (data) * numAdjRanks
*/
void GraphTraverser::setupOneSidedMPI()
{
    int mpiError;
    UINT maxMessages = 10;
    UINT maxPackets = maxMessages * g_maxCellsPerStep;
    
    
    // Allocate MPI_Win
    char *mpiWinMemory;
    MPI_Win mpiWin;
    UINT numAdjRanks = c_adjRankIndexToRank.size();
    UINT packetSize = 2 * sizeof(UINT) + c_dataSizeInBytes;
    UINT headerSize = sizeof(uint32_t) * 4 * numAdjRanks;
    UINT dataSize = 2 * maxPackets * packetSize * numAdjRanks;
    UINT windowSizeInBytes = headerSize + dataSize;

    MPI_Info mpiInfo;
    MPI_Info_create(&mpiInfo);
    MPI_Info_set(mpiInfo, "accumulate_ops", "same_op_no_op");
    //MPI_Info_set(mpiInfo, "accumulate_ordering", "none");
    MPI_Win_allocate(windowSizeInBytes, 1, mpiInfo, 
                     MPI_COMM_WORLD, &mpiWinMemory, &mpiWin);
    MPI_Info_free(&mpiInfo);


    // Setup onRankOffsets
    vector<UINT> onRankHeaderOffsets(numAdjRanks);
    vector<UINT> onRankDataOffsets(numAdjRanks);
    for (UINT i = 0; i < numAdjRanks; i++) {
        onRankHeaderOffsets[i] = i * 4 * sizeof(uint32_t);
        onRankDataOffsets[i] = headerSize + i * 2 * maxPackets * packetSize;
    }


    // Setup offRankOffsets
    vector<UINT> offRankHeaderOffsets(numAdjRanks);
    vector<UINT> offRankDataOffsets(numAdjRanks);
    vector<MPI_Request> mpiRequests(4 * numAdjRanks);
    for (UINT i = 0; i < numAdjRanks; i++) {
        UINT adjRank = c_adjRankIndexToRank[i];
        int tagHeader = 1;
        int tagData = 2;
        
        // Send header offset
        mpiError = MPI_Isend(&onRankHeaderOffsets[i], 1, MPI_UINT64_T, adjRank, 
                             tagHeader, MPI_COMM_WORLD, &mpiRequests[4*i+0]);
        Insist(mpiError == MPI_SUCCESS, "");

        // Send data offset
        mpiError = MPI_Isend(&onRankDataOffsets[i], 1, MPI_UINT64_T, adjRank, 
                             tagData, MPI_COMM_WORLD, &mpiRequests[4*i+1]);
        Insist(mpiError == MPI_SUCCESS, "");

        // Recv header offset
        mpiError = MPI_Irecv(&offRankHeaderOffsets[i], 1, MPI_UINT64_T, adjRank, 
                             tagHeader, MPI_COMM_WORLD, &mpiRequests[4*i+2]);
        Insist(mpiError == MPI_SUCCESS, "");
        
        // Recv data offset
        mpiError = MPI_Irecv(&offRankDataOffsets[i], 1, MPI_UINT64_T, adjRank, 
                             tagData, MPI_COMM_WORLD, &mpiRequests[4*i+3]);
        Insist(mpiError == MPI_SUCCESS, "");
    }
    

    // Wait for messages to send/recv
    if (numAdjRanks > 0) {
        mpiError = MPI_Waitall(mpiRequests.size(), mpiRequests.data(), 
                               MPI_STATUSES_IGNORE);
        Insist(mpiError == MPI_SUCCESS, "");
    }


    // Lock the window to start RMA operations and set memory to zero
    MPI_Win_lock_all(MPI_MODE_NOCHECK, mpiWin);
    memset(mpiWinMemory, 0, windowSizeInBytes);
    MPI_Win_sync(mpiWin);
    Comm::barrier();

    c_oneSidedImpl = new OneSidedImpl(numAdjRanks, 
                                      offRankHeaderOffsets, 
                                      offRankDataOffsets,
                                      onRankHeaderOffsets, 
                                      onRankDataOffsets,
                                      dataSize / 2 / numAdjRanks,
                                      mpiWin);
}


/*
    ~GraphTraverser
*/
GraphTraverser::~GraphTraverser()
{
    if (g_mpiType == MPIType_OneSided) {
        MPI_Win mpiWin = c_oneSidedImpl->getWin();
        MPI_Win_unlock_all(mpiWin);
        MPI_Win_free(&mpiWin);
    }
}


/*
    traverse
    
    Traverses g_tychoMesh.
*/
void GraphTraverser::traverse(const UINT maxComputePerStep,
                              TraverseData &traverseData)
{
    vector<priority_queue<Tuple>> canCompute(g_nThreads);
    Mat2<UINT> numDependencies(g_nAngles, g_nCells);
    UINT numCellAnglePairsToCalculate = g_nAngles * g_nCells;
    Mat2<vector<char>> sendBuffers;
    vector<vector<char>> sendBuffers1;
    vector<bool> commDark;
    Timer totalTimer;
    Timer setupTimer;
    Timer commTimer;
    Timer sendTimer;
    Timer recvTimer;
    

    // Start total timer
    totalTimer.start();
    setupTimer.start();
    
    
    // Calc num dependencies for each (cell, angle) pair
    for (UINT cell = 0; cell < g_nCells; cell++) {
    for (UINT angle = 0; angle < g_nAngles; angle++) {
        numDependencies(angle, cell) = c_initNumDependencies(angle, cell);
    }}
    
    
    // Set size of sendBuffers and commDark
    UINT numAdjRanks = c_adjRankIndexToRank.size();
    sendBuffers.resize(g_nThreads, numAdjRanks);
    sendBuffers1.resize(numAdjRanks);
    commDark.resize(numAdjRanks, false);
    
    
    // Initialize canCompute queue
    for (UINT cell = 0; cell < g_nCells; cell++) {
    for (UINT angle = 0; angle < g_nAngles; angle++) {
        if (numDependencies(angle, cell) == 0) {
            UINT priority = traverseData.getPriority(cell, angle);
            UINT angleGroup = angleGroupIndex(angle);
            canCompute[angleGroup].push(Tuple(cell, angle, priority));
        }
    }}


    // End setup timer
    setupTimer.stop();
    
    
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
                stepsTaken++;
                
                #pragma omp atomic
                numCellAnglePairsToCalculate--;
                
                
                // Get boundary type and adjacent cell/side data for each face
                BoundaryType bdryType[g_nFacePerCell];
                UINT adjCellsSides[g_nFacePerCell];
                bool isOutgoingWrtDirection[g_nFacePerCell];
                for (UINT face = 0; face < g_nFacePerCell; face++) {
                    
                    UINT adjCell = g_tychoMesh->getAdjCell(cell, face);
                    UINT adjRank = g_tychoMesh->getAdjRank(cell, face);
                    adjCellsSides[face] = adjCell;
                    
                    if (g_tychoMesh->isOutgoing(angle, cell, face)) {
                        
                        if (adjCell == TychoMesh::BOUNDARY_FACE && 
                            adjRank != TychoMesh::BAD_RANK)
                        {
                            bdryType[face] = BoundaryType_OutIntBdry;
                            adjCellsSides[face] = 
                                g_tychoMesh->getSide(cell, face);
                        }
                        
                        else if (adjCell == TychoMesh::BOUNDARY_FACE && 
                                 adjRank == TychoMesh::BAD_RANK)
                        {
                            bdryType[face] = BoundaryType_OutExtBdry;
                        }
                        
                        else {
                            bdryType[face] = BoundaryType_OutInt;
                        }
                        
                        if (c_direction == Direction_Forward) {
                            isOutgoingWrtDirection[face] = true;
                        }
                        else {
                            isOutgoingWrtDirection[face] = false;
                        }
                    }
                    else {
                        
                        if (adjCell == TychoMesh::BOUNDARY_FACE && 
                            adjRank != TychoMesh::BAD_RANK)
                        {
                            bdryType[face] = BoundaryType_InIntBdry;
                            adjCellsSides[face] = 
                                g_tychoMesh->getSide(cell, face);
                        }
                        
                        else if (adjCell == TychoMesh::BOUNDARY_FACE && 
                                 adjRank == TychoMesh::BAD_RANK)
                        {
                            bdryType[face] = BoundaryType_InExtBdry;
                        }
                        
                        else {
                            bdryType[face] = BoundaryType_InInt;
                        }
                        
                        if (c_direction == Direction_Forward) {
                            isOutgoingWrtDirection[face] = false;
                        }
                        else {
                            isOutgoingWrtDirection[face] = true;
                        }
                    }
                }
                
                
                // Update data for this cell-angle pair
                traverseData.update(cell, angle, adjCellsSides, bdryType);
                
                
                // Update dependency for children
                for (UINT face = 0; face < g_nFacePerCell; face++) {
                    
                    if (isOutgoingWrtDirection[face]) {

                        UINT adjCell = g_tychoMesh->getAdjCell(cell, face);
                        UINT adjRank = g_tychoMesh->getAdjRank(cell, face);
                        
                        if (adjCell != TychoMesh::BOUNDARY_FACE) {
                            numDependencies(angle, adjCell)--;
                            if (numDependencies(angle, adjCell) == 0) {
                                UINT priority = 
                                    traverseData.getPriority(adjCell, angle);
                                Tuple tuple(adjCell, angle, priority);
                                canCompute[angleGroup].push(tuple);
                            }
                        }
                        
                        else if (c_doComm && adjRank != TychoMesh::BAD_RANK) {
                            UINT rankIndex = c_adjRankToRankIndex.at(adjRank);
                            UINT side = g_tychoMesh->getSide(cell, face);
                            UINT globalSide = g_tychoMesh->getLGSide(side);
                            
                            vector<char> packet;
                            createPacket(packet, globalSide, angle, 
                                         c_dataSizeInBytes, 
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
        for (UINT angleGroup = 0; angleGroup < g_nThreads; angleGroup++) {
        for (UINT rankIndex = 0; rankIndex < numAdjRanks; rankIndex++) {
            sendBuffers1[rankIndex].insert(
                sendBuffers1[rankIndex].end(), 
                sendBuffers(angleGroup, rankIndex).begin(), 
                sendBuffers(angleGroup, rankIndex).end());
        }}
               
        
        // Do communication
        commTimer.start();
        if (c_doComm) {
            
            // Send/Recv
            UINT packetSizeInBytes = 2 * sizeof(UINT) + c_dataSizeInBytes;
            vector<char> dataPackets;
            
            if (g_mpiType == MPIType_TychoTwoSided) {
                const bool killComm = false;
                sendAndRecvData(sendBuffers1, dataPackets, commDark, killComm);
            }
            else if (g_mpiType == MPIType_CapsaicinTwoSided) {
                sendTimer.start();
                sendData2Sided(sendBuffers1);
                sendTimer.stop();

                recvTimer.start();
                recvData2Sided(dataPackets);
                recvTimer.stop();
            }
            else if (g_mpiType == MPIType_OneSided) {
                recvTimer.start();
                recvData1Sided(dataPackets);
                recvTimer.stop();

                sendTimer.start();
                sendData1Sided(sendBuffers1);
                sendTimer.stop();
            }
            else {
                Insist(false, "MPI type not recognized.");
            }

            
            // Clear send buffers for next iteration
            for (UINT angleGroup = 0; angleGroup < g_nThreads; angleGroup++) {
            for (UINT rankIndex = 0; rankIndex < numAdjRanks; rankIndex++) {
                sendBuffers(angleGroup, rankIndex).clear();
            }}
            
            for (UINT rankIndex = 0; rankIndex < numAdjRanks; rankIndex++) {
                sendBuffers1[rankIndex].clear();
            }
            
            
            // Unpack packets
            UINT numPackets = dataPackets.size() / packetSizeInBytes;
            Assert(dataPackets.size() % packetSizeInBytes == 0);

            for (UINT i = 0; i < numPackets; i++) {
                char *packet = &dataPackets[i * packetSizeInBytes];
                UINT globalSide;
                UINT angle;
                char *packetData;
                splitPacket(packet, globalSide, angle, &packetData);
                
                UINT localSide = g_tychoMesh->getGLSide(globalSide);
                traverseData.setSideData(localSide, angle, packetData);

                UINT cell = g_tychoMesh->getSideCell(localSide);
                numDependencies(angle, cell)--;
                if (numDependencies(angle, cell) == 0) {
                    UINT priority = traverseData.getPriority(cell, angle);
                    Tuple tuple(cell, angle, priority);
                    canCompute[angleGroupIndex(angle)].push(tuple);
                }
            }
        }
        commTimer.stop();
    }
    
    
    // Send kill comm signal to adjacent ranks
    if (g_mpiType == MPIType_TychoTwoSided) {
        commTimer.start();
        if (c_doComm) {
            vector<char> dataPackets;
            const bool killComm = true;
            sendAndRecvData(sendBuffers1, dataPackets, commDark, killComm);
        }
        commTimer.stop();
    }

    
    // Print times
    totalTimer.stop();

    double totalTime = totalTimer.wall_clock();
    Comm::gmax(totalTime);

    double setupTime = setupTimer.wall_clock();
    Comm::gmax(setupTime);

    double commTime = commTimer.sum_wall_clock();
    Comm::gmax(commTime);
    
    double sendTime = sendTimer.sum_wall_clock();
    Comm::gmax(sendTime);
    
    double recvTime = recvTimer.sum_wall_clock();
    Comm::gmax(recvTime);
    
    if (Comm::rank() == 0) {
        printf("      Traverse Timer (comm):    %fs\n", commTime);
        printf("      Traverse Timer (send):    %fs\n", sendTime);
        printf("      Traverse Timer (recv):    %fs\n", recvTime);
        printf("      Traverse Timer (setup):   %fs\n", setupTime);
        printf("      Traverse Timer (total):   %fs\n", totalTime);
    }
}



