/*
 CommSides.cc

Communicates boundary values between the sides
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
#include "CommSides.hh"


/*
    commSides
*/
void commSides(const std::vector<UINT> &adjRanks,
               const std::vector<std::vector<MetaData>> &sendMetaData,
               const std::vector<UINT> &numSendPackets,
               const std::vector<UINT> &numRecvPackets,
               SweepData2 &sweepData)
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


