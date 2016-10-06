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


#include "CommSides.hh"
#include "Global.hh"
#include "Comm.hh"
#include <vector>
#include <algorithm>
#include <string.h>


/*
    Constructor
*/
CommSides::CommSides()
{
    // Get adjacent ranks
    for (UINT cell = 0; cell < g_nCells; cell++) {
    for (UINT face = 0; face < g_nFacePerCell; face++) {
        
        UINT adjRank = g_tychoMesh->getAdjRank(cell, face);
        UINT adjCell = g_tychoMesh->getAdjCell(cell, face);
        
        if (adjCell == TychoMesh::BOUNDARY_FACE && 
            adjRank != TychoMesh::BAD_RANK &&
            std::count(c_adjRanks.begin(), c_adjRanks.end(), adjRank) == 0)
        {
            c_adjRanks.push_back(adjRank);
        }
    }}
    
    
    // Populate sendMetaData, numSendPackets, and numRecvPackets
    c_sendMetaData.resize(c_adjRanks.size());
    c_numSendPackets.resize(c_adjRanks.size());
    c_numRecvPackets.resize(c_adjRanks.size());
    
    for (UINT rankIndex = 0; rankIndex < c_adjRanks.size(); rankIndex++) {
        
        c_numSendPackets[rankIndex] = 0;
        c_numRecvPackets[rankIndex] = 0;
        
        for (UINT cell = 0; cell < g_nCells; cell++) {
        for (UINT face = 0; face < g_nFacePerCell; face++) {
        
            UINT adjRank = g_tychoMesh->getAdjRank(cell, face);        
            if (adjRank == c_adjRanks[rankIndex]) {
                for (UINT angle = 0; angle < g_nAngles; angle++) {
                    if (g_tychoMesh->isOutgoing(angle, cell, face)) {
                        CommSides::MetaData md;
                        UINT side = g_tychoMesh->getSide(cell, face);
                        md.gSide = g_tychoMesh->getLGSide(side);
                        md.angle = angle;
                        md.cell  = cell;
                        md.face  = face;
                        c_sendMetaData[rankIndex].push_back(md);
                        
                        c_numSendPackets[rankIndex]++;
                    }
                    else {
                        c_numRecvPackets[rankIndex]++;
                    }
                }
            }
        }}
    }
}


/*
    getDataSize

    Returns data size of 1 cell/face of data to send.
*/
static UINT getDataSize()
{
    return g_nGroups * g_nVrtxPerFace * sizeof(double);
}


/*
    commSides
*/
void CommSides::commSides(PsiData &psi, PsiBoundData &psiBound)
{
    int mpiError;
    UINT numToRecv;
    UINT numAdjRanks = c_adjRanks.size();
    UINT packetSize = 2 * sizeof(UINT) + getDataSize();
    std::vector<MPI_Request> mpiRecvRequests(numAdjRanks);
    std::vector<MPI_Request> mpiSendRequests(numAdjRanks);
    std::vector<std::vector<char>> dataToSend(numAdjRanks);
    std::vector<std::vector<char>> dataToRecv(numAdjRanks);
    Mat2<double> localFaceData(g_nVrtxPerFace, g_nGroups);
    
    
    // Data structures to send/recv packets
    for (UINT rankIndex = 0; rankIndex < numAdjRanks; rankIndex++) {
        dataToSend[rankIndex].resize(packetSize * c_numSendPackets[rankIndex]);
        dataToRecv[rankIndex].resize(packetSize * c_numRecvPackets[rankIndex]);
    }
    
    
    // Irecv data
    numToRecv = 0;
    for (UINT rankIndex = 0; rankIndex < numAdjRanks; rankIndex++) {
        
        if (dataToRecv[rankIndex].size() > 0) {
            int tag = 0;
            int adjRank = c_adjRanks[rankIndex];
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
                 metaDataIndex < c_sendMetaData[rankIndex].size(); 
                 metaDataIndex++)
            {
                UINT gSide = c_sendMetaData[rankIndex][metaDataIndex].gSide;
                UINT angle = c_sendMetaData[rankIndex][metaDataIndex].angle;
                UINT cell  = c_sendMetaData[rankIndex][metaDataIndex].cell;
                UINT face  = c_sendMetaData[rankIndex][metaDataIndex].face;
                
                for (UINT group = 0; group < g_nGroups; group++) {
                for (UINT fvrtx = 0; fvrtx < g_nVrtxPerFace; fvrtx++) {
                    UINT vrtx = g_tychoMesh->getFaceToCellVrtx(cell, face, fvrtx);
                    localFaceData(fvrtx, group) = psi(group, vrtx, angle, cell);
                }}

                const char *data = (char*) (&localFaceData[0]);
                
                char *ptr = &dataToSend[rankIndex][metaDataIndex * packetSize];
                memcpy(ptr, &gSide, sizeof(UINT));
                ptr += sizeof(UINT);
                memcpy(ptr, &angle, sizeof(UINT));
                ptr += sizeof(UINT);
                memcpy(ptr, data, getDataSize());
            }
            
            int tag = 0;
            int adjRank = c_adjRanks[rankIndex];
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
            UINT side = g_tychoMesh->getGLSide(gSide);

            localFaceData.setData((double*)ptr);
            for (UINT fvrtx = 0; fvrtx < g_nVrtxPerFace; fvrtx++) {
            for (UINT group = 0; group < g_nGroups; group++) {
                psiBound(group, fvrtx, angle, side) = localFaceData(fvrtx, group);
            }}
        }
    }
    
    
    // Wait on send to complete
    if (mpiSendRequests.size() > 0) {
        mpiError = MPI_Waitall(mpiSendRequests.size(), mpiSendRequests.data(), 
                               MPI_STATUSES_IGNORE);
        Insist(mpiError == MPI_SUCCESS, "");
    }
}

