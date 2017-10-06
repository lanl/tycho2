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

#ifndef __GRAPH_TRAVERSER_HH__
#define __GRAPH_TRAVERSER_HH__

#include "Global.hh"
#include "Mat.hh"
#include "PsiData.hh"
#include "Transport.hh"
#include <mpi.h>
#include <omp.h>
#include <vector>
#include <map>


/*
    SweepData
    
    Holds psi and other data for the sweep.
*/
class SweepData
{
public:
    
    SweepData(PsiData &psi, const PsiData &source, PsiBoundData &psiBound)
    : c_psi(psi), c_psiBound(psiBound), c_source(source), 
      c_localSource(g_nThreads), c_localPsi(g_nThreads),
      c_localPsiBound(g_nThreads)
    {
        for (UINT angleGroup = 0; angleGroup < g_nThreads; angleGroup++) {
            c_localSource[angleGroup].resize(g_nVrtxPerCell, g_nGroups);
            c_localPsi[angleGroup].resize(g_nVrtxPerCell, g_nGroups);
            c_localPsiBound[angleGroup].resize(g_nVrtxPerFace, g_nFacePerCell, 
                                               g_nGroups);
        }
    }
    

    /*
        update
        
        Does a transport update for the given cell/angle pair.
    */
    void update(UINT cell, UINT angle) 
    {
        
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
    std::vector<Mat2<double>> c_localSource;
    std::vector<Mat2<double>> c_localPsi;
    std::vector<Mat3<double>> c_localPsiBound;
};




class GraphTraverser
{
public:
    GraphTraverser();
    void traverse(SweepData &traverseData);

private:

    Mat2<UINT> c_initNumDependencies;
};

#endif
