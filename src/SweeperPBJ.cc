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
#include "SweepData2.hh"
#include "CommSides.hh"
#include <omp.h>
#include <vector>
#include <math.h>
#include <string.h>

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
    SweepDataPBJ(PsiData &psi, const PsiData &source, PsiBoundData &psiBound, 
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
            localFaceData(fvrtx, group) = c_psi(group, vrtx, angle, cell);
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
            c_psiBound(group, fvrtx, angle, side) = localFaceData(fvrtx, group);
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
            localSource(vrtx, group) = c_source(group, vrtx, angle, cell);
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
            c_psi(group, vrtx, angle, cell) = localPsi(vrtx, group);
        }}
    }
    
private:
    PsiData &c_psi;
    PsiBoundData &c_psiBound;
    const PsiData &c_source;
    const double c_sigmaTotal;
    vector<Mat2<double>> c_localFaceData;
};



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
    PsiBoundData psiBound;
    psiBound.setToValue(0.0);  // Change to something more reasonable.
    
    
    // Set psi0
    PsiData psi0;
    for (UINT i = 0; i < psi.size(); i++) {
        psi0[i] = psi[i];
    }
    
    
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
    
    
    // Sweep till converged
    UINT iter = 1;
    while (iter < maxIter) {
        
        // Sweep
        traverseGraph(maxComputePerStep, sweepData, doComm, MPI_COMM_WORLD, 
                      Direction_Forward);
        
        
        // Check tolerance and set psi0 = psi1
        double errL1 = 0.0;
        double normL1 = 0.0;
        for (UINT i = 0; i < psi.size(); i++) {
            errL1  += fabs(psi0[i] - psi[i]);
            normL1 += fabs(psi[i]);
            
            psi0[i] = psi[i];
        }
        
        Comm::gsum(errL1);
        Comm::gsum(normL1);

        if (errL1 / normL1 < tolerance)
            break;
        
        
        // Communicate
        commSides(adjRanks, sendMetaData, numSendPackets, numRecvPackets, 
                  sweepData);
        
        
        // Increment iter
        iter++;
    }
    
    
    // Print statistics
    if (Comm::rank() == 0) {
        printf("      PBJ Iters: %" PRIu64 "\n", iter);
    }
}


