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

#ifndef __SWEEP_DATA_HH__
#define __SWEEP_DATA_HH__


#include "Assert.hh"
#include "GraphTraverser.hh"
#include "Transport.hh"
#include "Global.hh"
#include <stddef.h>
#include <omp.h>

/*
    SweepData
    
    Holds psi and other data for the sweep.
*/
class SweepData : public TraverseData
{
public:
    
    SweepData(PsiData &psi, const PsiData &source, PsiBoundData &psiBound,  
               const Mat2<UINT> &priorities)
    : c_psi(psi), c_psiBound(psiBound), c_source(source), 
      c_priorities(priorities), c_localFaceData(g_nThreads),
      c_localSource(g_nThreads), c_localPsi(g_nThreads),
      c_localPsiBound(g_nThreads)
    {
        for (UINT angleGroup = 0; angleGroup < g_nThreads; angleGroup++) {
            c_localFaceData[angleGroup].resize(g_nVrtxPerFace, g_nGroups);
            c_localSource[angleGroup].resize(g_nVrtxPerCell, g_nGroups);
            c_localPsi[angleGroup].resize(g_nVrtxPerCell, g_nGroups);
            c_localPsiBound[angleGroup].resize(g_nVrtxPerFace, g_nFacePerCell, 
                                               g_nGroups);
        }
    }
    

    /*
        getDataSizeInBytes
    */
    static
    size_t getDataSizeInBytes()
    {
        return g_nGroups * g_nVrtxPerFace * sizeof(double);
    }
    
    
    /*
        data
        
        Return psi for vertices and groups at the given (cell,face,angle) tuple
    */
    virtual const char* getData(UINT cell, UINT face, UINT angle)
    {
        Mat2<double> &localFaceData = c_localFaceData[omp_get_thread_num()];
        
        for (UINT group = 0; group < g_nGroups; group++) {
        for (UINT fvrtx = 0; fvrtx < g_nVrtxPerFace; fvrtx++) {
            UINT vrtx = g_tychoMesh->getFaceToCellVrtx(cell, face, fvrtx);
            localFaceData(fvrtx, group) = c_psi(group, vrtx, angle, cell);
        }}
        
        return (char*) (&localFaceData[0]);
    }
       
        
    /*
        sideData
        
        Set psiBound for the (side, angle) pair.
    */
    virtual void setSideData(UINT side, UINT angle, const char *data)
    {
        Mat2<double> localFaceData(g_nVrtxPerFace, g_nGroups);
        localFaceData.setData((double*)data);
        
        for (UINT fvrtx = 0; fvrtx < g_nVrtxPerFace; fvrtx++) {
        for (UINT group = 0; group < g_nGroups; group++) {
            c_psiBound(group, fvrtx, angle, side) = localFaceData(fvrtx, group);
        }}
    }


    /*
        getPriority
        
        Return a priority for the cell/angle pair.
    */
    virtual UINT getPriority(UINT cell, UINT angle)
    {
        return c_priorities(cell, angle);
    }
    
    
    /*
        update
        
        Does a transport update for the given cell/angle pair.
    */
    virtual void update(UINT cell, UINT angle, 
                        UINT adjCellsSides[g_nFacePerCell], 
                        BoundaryType bdryType[g_nFacePerCell])
    {
        UNUSED_VARIABLE(adjCellsSides);
        UNUSED_VARIABLE(bdryType);
        
        Mat2<double> &localSource = c_localSource[omp_get_thread_num()];
        Mat2<double> &localPsi = c_localPsi[omp_get_thread_num()];
        Mat3<double> &localPsiBound = c_localPsiBound[omp_get_thread_num()];

        
        // Populate localSource
        #pragma omp simd
        for (UINT group = 0; group < g_nGroups; group++) {
        for (UINT vrtx = 0; vrtx < g_nVrtxPerCell; vrtx++) {
            localSource(vrtx, group) = c_source(group, vrtx, angle, cell);
        }}
        
        
        // Populate localPsiBound
        Transport::populateLocalPsiBound(angle, cell, c_psi, c_psiBound, 
                                         localPsiBound);
        
        
        // Transport solve
        Transport::solve(cell, angle, g_sigmaT[cell],
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
    const Mat2<UINT> &c_priorities;
    std::vector<Mat2<double>> c_localFaceData;
    std::vector<Mat2<double>> c_localSource;
    std::vector<Mat2<double>> c_localPsi;
    std::vector<Mat3<double>> c_localPsiBound;
};

#endif

