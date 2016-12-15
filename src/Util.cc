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

#include "Util.hh"
#include "Global.hh"
#include "Comm.hh"
#include "SweepData.hh"
#include "CommSides.hh"
#include <math.h>
#include <limits>


namespace Util
{

/*
    diffBetweenGroups
*/
double diffBetweenGroups(const PsiData &psi)
{
    double maxDiff = 0.0;
    double maxEntry = 0.0;
    
    for(UINT cell = 0; cell < g_nCells; cell++) {
    for(UINT angle = 0; angle < g_nAngles; angle++) {
    for(UINT vrtx = 0; vrtx < g_nVrtxPerCell; vrtx++) {
        
        double psi0 = psi(0, vrtx, angle, cell);
        
        if(fabs(psi0) > maxEntry)
            maxEntry = fabs(psi0);
        
        for(UINT group = 1; group < g_nGroups; group++) {
            double psi1 = psi(group, vrtx, angle, cell);
            if (fabs(psi0 - psi1) > maxDiff)
                maxDiff = fabs(psi0 - psi1);
        }
    }}}
    
    maxDiff = maxDiff / maxEntry;
    Comm::gmax(maxDiff);
    return maxDiff;
}


/*
    psiToPhi
*/
void psiToPhi(PhiData &phi, const PsiData &psi) 
{
    phi.setToValue(0.0);
    
    #pragma omp parallel for
    for (UINT cell = 0; cell < g_nCells; ++cell) {
    for (UINT angle = 0; angle < g_nAngles; ++angle) {
    for (UINT vertex = 0; vertex < g_nVrtxPerCell; ++vertex) {
    for (UINT group = 0; group < g_nGroups; ++group) {
        phi(group, vertex, cell) +=
            psi(group, vertex, angle, cell) * g_quadrature->getWt(angle);
    }}}}
}


/*
    phiToPsi
*/
void phiToPsi(const PhiData &phi, PsiData &psi) 
{
    #pragma omp parallel for
    for (UINT cell = 0; cell < g_nCells; ++cell) {
    for (UINT angle = 0; angle < g_nAngles; ++angle) {
    for (UINT vertex = 0; vertex < g_nVrtxPerCell; ++vertex) {
    for (UINT group = 0; group < g_nGroups; ++group) {
        psi(group, vertex, angle, cell) = phi(group, vertex, cell);
    }}}}
}


/*
    calcTotalSource
*/
void calcTotalSource(const PsiData &source, const PhiData &phi, 
                     PsiData &totalSource)
{
    #pragma omp parallel for
    for (UINT cell = 0; cell < g_nCells; ++cell) {
    for (UINT angle = 0; angle < g_nAngles; ++angle) {
    for (UINT vertex = 0; vertex < g_nVrtxPerCell; ++vertex) {
    for (UINT group = 0; group < g_nGroups; ++group) {
        totalSource(group, vertex, angle, cell) = 
            source(group, vertex, angle, cell) + 
            g_sigmaS[cell] / (4.0 * M_PI) *  phi(group, vertex, cell);
    }}}}
}


/*
    sweepLocal

    Solves L_I Psi = L_B Psi_B + Q
*/
void sweepLocal(PsiData &psi, const PsiData &source, PsiBoundData &psiBound)
{
    Mat2<UINT> priorities(g_nCells, g_nAngles);
    const UINT maxComputePerStep = std::numeric_limits<uint64_t>::max();
    SweepData sweepData(psi, source, psiBound, priorities);
    
    g_graphTraverserForward->traverse(maxComputePerStep, sweepData);
}


/*
    operatorS
*/
void operatorS(const PhiData &phi1, PhiData &phi2)
{
    for (UINT cell = 0; cell < g_nCells; ++cell) {
    for (UINT vertex = 0; vertex < g_nVrtxPerCell; ++vertex) {
    for (UINT group = 0; group < g_nGroups; ++group) {
        phi2(group,vertex,cell) = g_sigmaS[cell] / (4.0 * M_PI) * 
                                  phi1(group,vertex,cell);
    }}}
}


} // End namespace Util
