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

#include "Solver.hh"
#include "Global.hh"
#include "PsiData.hh"
#include "Sweeper.hh"
#include "Sweeper2.hh"
#include "Quadrature.hh"
#include "Comm.hh"
#include "Timer.hh"
#include <math.h>
#include <vector>


using namespace std;
static const double cubeSize = 100.0;


/*
    source_exponential
*/
// static
// void exponentialSource(const double sigmaT, const double sigmaS, PsiData &source)
// {
//     double L = 32.0 / (cubeSize / 2.0 * cubeSize / 2.0);
//     
//     for(UINT vrtx = 0; vrtx < g_nVrtxPerCell; vrtx++) {
//     for(UINT angle = 0; angle < g_quadrature->getNumAngles(); angle++) {
//     for(UINT cell = 0; cell < g_spTychoMesh->getNCells(); cell++) {
//     for(UINT group = 0; group < g_nGroups; group++) {
//         UINT node = g_spTychoMesh->getCellNode(cell, vrtx);
//         double x = g_spTychoMesh->getNodeCoord(node, 0) - cubeSize / 2.0;
//         double y = g_spTychoMesh->getNodeCoord(node, 1) - cubeSize / 2.0;
//         double z = g_spTychoMesh->getNodeCoord(node, 2) - cubeSize / 2.0;
//         
//         double xi  = g_quadrature->getXi(angle);
//         double eta = g_quadrature->getEta(angle);
//         double mu  = g_quadrature->getMu(angle);
//         
//         double e = exp(-(x*x + y*y + z*z) / (2.0 * L));
//         
//         source(vrtx, angle, cell, group) = 
//             - (x * xi * e / L) - (y * eta * e / L) - (z * mu * e / L)
//             + (sigmaT - sigmaS) * e;
//     }}}}
// }


/*
    hatSource
*/
static
void hatSource(const double sigmaT, const double sigmaS, PsiData &source)
{
    for(UINT vrtx = 0; vrtx < g_nVrtxPerCell; vrtx++) {
    for(UINT angle = 0; angle < g_quadrature->getNumAngles(); angle++) {
    for(UINT cell = 0; cell < g_spTychoMesh->getNCells(); cell++) {
        
        UINT node = g_spTychoMesh->getCellNode(cell, vrtx);
        double x = g_spTychoMesh->getNodeCoord(node, 0) - cubeSize / 2.0;
        double y = g_spTychoMesh->getNodeCoord(node, 1) - cubeSize / 2.0;
        double z = g_spTychoMesh->getNodeCoord(node, 2) - cubeSize / 2.0;
        
        double xi  = g_quadrature->getXi(angle);
        double eta = g_quadrature->getEta(angle);
        double mu  = g_quadrature->getMu(angle);
        
        double c = sqrt(x*x + y*y + z*z);
        
        for(UINT group = 0; group < g_nGroups; group++) {
            if(c <= 30.0) {
                source(vrtx, angle, cell, group) = 
                    - x / (30.0*c) * xi - y / (30.0*c) * eta - z / (30.0*c) * mu
                    + (sigmaT - sigmaS) * (1.0 - c / 30.0);
            }
            else {
                source(vrtx, angle, cell, group) = 0.0;
            }
        }
    }}}
}


/*
    exponentialL2Error
*/
// static
// double exponentialL2Error(const PsiData &psi)
// {
//     double L = 32.0 / (cubeSize / 2.0 * cubeSize / 2.0);
//     double diff = 0.0;
//     double norm = 0.0;
//     
//     for(UINT angle = 0; angle < g_quadrature->getNumAngles(); angle++) {
//     for(UINT cell = 0; cell < g_spTychoMesh->getNCells(); cell++) {
//     for(UINT group = 0; group < g_nGroups; group++) {
//         double x = 0.0;
//         double y = 0.0;
//         double z = 0.0;
//         double psiVal = 0.0;
//         
//         for(UINT vrtx = 0; vrtx < g_nVrtxPerCell; vrtx++) {
//             UINT node = g_spTychoMesh->getCellNode(cell, vrtx);
//             x += g_spTychoMesh->getNodeCoord(node, 0) - cubeSize / 2.0;
//             y += g_spTychoMesh->getNodeCoord(node, 1) - cubeSize / 2.0;
//             z += g_spTychoMesh->getNodeCoord(node, 2) - cubeSize / 2.0;
//             psiVal += psi(vrtx, angle, cell, group);
//         }
//         
//         x = x / g_nVrtxPerCell;
//         y = y / g_nVrtxPerCell;
//         z = z / g_nVrtxPerCell;
//         psiVal = psiVal / g_nVrtxPerCell;
//         double psiAct = exp(-(x*x + y*y + z*z) / (2.0 * L));
//         
//         double localDiff = psiVal - psiAct;
//         norm += psiVal * psiVal;
//         diff += localDiff * localDiff;
//     }}}
//     
//     return sqrt(diff / norm);
// }


/*
    hatL2Error
*/
static
double hatL2Error(const PsiData &psi)
{
    double diff = 0.0;
    double norm = 0.0;
    
    
    for(UINT angle = 0; angle < g_quadrature->getNumAngles(); angle++) {
    for(UINT cell = 0; cell < g_spTychoMesh->getNCells(); cell++) {
    for(UINT group = 0; group < g_nGroups; group++) {
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        double psiVal = 0.0;
        
        for(UINT vrtx = 0; vrtx < g_nVrtxPerCell; vrtx++) {
            UINT node = g_spTychoMesh->getCellNode(cell, vrtx);
            x += g_spTychoMesh->getNodeCoord(node, 0) - cubeSize / 2.0;
            y += g_spTychoMesh->getNodeCoord(node, 1) - cubeSize / 2.0;
            z += g_spTychoMesh->getNodeCoord(node, 2) - cubeSize / 2.0;
            psiVal += psi(vrtx, angle, cell, group);
        }
        
        x = x / g_nVrtxPerCell;
        y = y / g_nVrtxPerCell;
        z = z / g_nVrtxPerCell;
        psiVal = psiVal / g_nVrtxPerCell;
        
        double c = sqrt(x*x + y*y + z*z);
        double psiAct = 0.0;
        if(c <= 30.0) {
            psiAct = 1.0 - c / 30.0;
        }
        
        double localDiff = psiVal - psiAct;
        norm += psiVal * psiVal;
        diff += localDiff * localDiff;
    }}}
    
    
    Comm::gsum(norm);
    Comm::gsum(diff);
    
    return sqrt(diff / norm);
}


/*
    diffBetweenGroups
*/
static
double diffBetweenGroups(const PsiData &psi)
{
    double maxDiff = 0.0;
    double maxEntry = 0.0;
    
    for(UINT vrtx = 0; vrtx < g_nVrtxPerCell; vrtx++) {
    for(UINT angle = 0; angle < g_quadrature->getNumAngles(); angle++) {
    for(UINT cell = 0; cell < g_spTychoMesh->getNCells(); cell++) {
        
        double psi0 = psi(vrtx, angle, cell, 0);
        
        if(fabs(psi0) > maxEntry)
            maxEntry = psi0;
        
        for(UINT group = 1; group < g_nGroups; group++) {
            double psi1 = psi(vrtx, angle, cell, group);
            if (fabs(psi0 - psi1) > maxDiff)
                maxDiff = fabs(psi0 - psi1);
        }
    }}}
    
    maxDiff = maxDiff / maxEntry;
    Comm::gsum(maxDiff);
    return maxDiff;
}


/*
    calcMass
*/
static
double calcMass(const PsiData &psi)
{
    double mass = 0.0;
    
    for(UINT cell = 0; cell < g_spTychoMesh->getNCells(); cell++) {
        double sumVrtxMass = 0.0;
        for(UINT vrtx = 0; vrtx < g_nVrtxPerCell; vrtx++) {
            double localMass = 0.0;
            for(UINT angle = 0; angle < g_quadrature->getNumAngles(); angle++) {
                localMass += psi(vrtx, angle, cell, 0) * g_quadrature->getWt(angle);
            }
            sumVrtxMass += localMass;
        }
        mass += sumVrtxMass / g_nVrtxPerCell * g_spTychoMesh->getCellVolume(cell);
    }
    
    Comm::gsum(mass);
    return mass;
}


/*
    psiToPhi
*/
static
void psiToPhi(PhiData &phi, const PsiData &psi) 
{
    phi.setToValue(0.0);
    
    //for (UINT angle = 0; angle < g_quadrature->getNumAngles(); ++angle) {
    #pragma omp parallel for
    for (UINT cell = 0; cell < g_spTychoMesh->getNCells(); ++cell) {
    for (UINT angle = 0; angle < g_quadrature->getNumAngles(); ++angle) {
    for (UINT vertex = 0; vertex < g_nVrtxPerCell; ++vertex) {
    for (UINT group = 0; group < g_nGroups; ++group) {
        phi(vertex, cell, group) +=
            psi(vertex, angle, cell, group) * g_quadrature->getWt(angle);
    }}}}
}


/*
    calcScatterSource
*/
static
void calcScatterSource(PsiData &scatSource, const double sigmaScat,
                       const PhiData &phiOld) 
{
    //for (UINT angle = 0; angle < g_quadrature->getNumAngles(); ++angle) {
    #pragma omp parallel for
    for (UINT cell = 0; cell < g_spTychoMesh->getNCells(); ++cell) {
    for (UINT angle = 0; angle < g_quadrature->getNumAngles(); ++angle) {
    for (UINT vertex = 0; vertex < g_nVrtxPerCell; ++vertex) {
    for (UINT group = 0; group < g_nGroups; ++group) {
        // need to loop over moments in the future
        scatSource(vertex, angle, cell, group) =
            sigmaScat / (4.0 * M_PI) *  phiOld(vertex, cell, group);
    }}}}
}



static
void splitAnglesAcrossThreads(vector<vector<UINT>> &anglesVector)
{
    // Get the angle indices for each angle group
    vector<UINT> angleBdryIndices(g_nAngleGroups + 1);
    angleBdryIndices[0] = 0;
    for (UINT angleGroup = 0; angleGroup < g_nAngleGroups; angleGroup++) {
        UINT numAngles = g_quadrature->getNumAngles() / g_nAngleGroups;
        if (angleGroup < g_quadrature->getNumAngles() % g_nAngleGroups)
            numAngles++;
        angleBdryIndices[angleGroup+1] = angleBdryIndices[angleGroup] + numAngles;
    }
    
    // Create a SweepSchedule for each angle group
    for (UINT angleGroup = 0; angleGroup < g_nAngleGroups; angleGroup++) {
        UINT numAngles = angleBdryIndices[angleGroup+1] - angleBdryIndices[angleGroup];
        anglesVector[angleGroup].resize(numAngles);
        vector<UINT> &angles = anglesVector[angleGroup];
        for (UINT angle = 0; angle < numAngles; angle++) {
            angles[angle] = angle + angleBdryIndices[angleGroup];
        }
    }
}


namespace Solver
{

/*
    Solve problem
*/
void solve(const double sigmaTotal, const double sigmaScat, 
           const UINT iterMax, const double errMax)
{
    // Data for problem
    PsiData fixedSource(g_nVrtxPerCell, g_quadrature->getNumAngles(), 
                        g_spTychoMesh->getNCells(), g_nGroups);
    PsiData scatSource(g_nVrtxPerCell,  g_quadrature->getNumAngles(), 
                       g_spTychoMesh->getNCells(), g_nGroups);
    PsiData totalSource(g_nVrtxPerCell, g_quadrature->getNumAngles(), 
                        g_spTychoMesh->getNCells(), g_nGroups);
    
    PsiData psi(g_nVrtxPerCell, g_quadrature->getNumAngles(), 
                g_spTychoMesh->getNCells(), g_nGroups);
    
    PhiData phiNew(g_nVrtxPerCell, g_spTychoMesh->getNCells(), g_nGroups);
    PhiData phiOld(g_nVrtxPerCell, g_spTychoMesh->getNCells(), g_nGroups);
    
    hatSource(sigmaTotal, sigmaScat, fixedSource);
    double mass = calcMass(fixedSource);
    if(Comm::rank() == 0) {
        printf("Mass of Q: %e\n", mass);
    }
    
    
    // Volume of cube
    double volume = 0.0;
    for(UINT cell = 0; cell < g_spTychoMesh->getNCells(); cell++) {
        volume += g_spTychoMesh->getCellVolume(cell);
    }
    Comm::gsum(volume);
    if(Comm::rank() == 0) {
        printf("Volume of cube: %e\n", volume);
    }
    
    
    // Solver 2 class
    Sweeper2 *sweeper2 = NULL;
    if (g_sweepType == SweepType_TraverseGraph) {
        vector<vector<UINT>> anglesVector(g_nAngleGroups);
        splitAnglesAcrossThreads(anglesVector);
        sweeper2 = new Sweeper2(g_maxCellsPerStep,
                                g_intraAngleP, g_interAngleP, sigmaTotal);
    }
    
    
    // Source iteration
    UINT iter = 0;
    double error = 1.0;
    Timer innerTimer;
    innerTimer.start();
    while (iter < iterMax && error > errMax)
    {
        Timer timer, timer2, timer3, timer4;
        double wallClockTime = 0.0;
        double clockTime2 = 0.0;
        double clockTime3 = 0.0;
        //double clockTime4 = 0.0;
        timer.start();
        
        for(UINT i = 0; i < phiOld.size(); i++) {
            phiOld[i] = phiNew[i];
        }
        
        // Need to update fixed source with scattering source here
        timer2.start();
        calcScatterSource(scatSource, sigmaScat, phiOld);
        #pragma omp parallel for
        for(UINT i = 0; i < totalSource.size(); i++) {
            totalSource[i] = fixedSource[i] + scatSource[i];
        }
        timer2.stop();
        clockTime2 = timer2.wall_clock();
        Comm::gmax(clockTime2);
        if(Comm::rank() == 0) {
            printf("\n   Source time: %f\n", clockTime2);
        }
        
        // Sweep
        switch (g_sweepType) {
            case SweepType_OriginalTycho1:
            case SweepType_OriginalTycho2:
                Sweeper::sweep(psi, totalSource, sigmaTotal);
                break;
            case SweepType_TraverseGraph:
                sweeper2->sweep(psi, totalSource);
                break;
            default:
                Insist(false, "Sweep type not recognized.");
                break;
        }
        
        // Need to incorporate convergence tests here to iterate as necessary
        timer3.start();
        psiToPhi(phiNew, psi);
        error = 0.0;
        for (UINT element = 0; element < phiNew.size(); ++element) {
            error = max(error, 
                fabs(phiNew[element] - phiOld[element]) / phiNew[element]);
        }
        Comm::gmax(error);
        timer3.stop();
        clockTime3 = timer3.wall_clock();
        Comm::gmax(clockTime3);
        if(Comm::rank() == 0) {
            printf("   psiToPhi Time: %f\n", clockTime3);
        }
        
        timer.stop();
        wallClockTime = timer.wall_clock();
        Comm::gmax(wallClockTime);
        if(Comm::rank() == 0) {
            printf("   iteration: %" PRIu64 "   error: %f   time: %f\n", 
                   iter, error, wallClockTime);
        }
        
        ++iter;
    }
    innerTimer.stop();
    double clockTime = innerTimer.wall_clock();
    Comm::gmax(clockTime);
    if(Comm::rank() == 0) {
        printf("Total solve time: %.2f\n\n",
               clockTime);
        printf("Solve time per iteration: %.2f\n\n",
               clockTime / iter);
    }
    
    
    // Print tests
    mass = calcMass(psi);
    double psiError = hatL2Error(psi);
    double diffGroups = diffBetweenGroups(psi);
    if(Comm::rank() == 0) {
        printf("Mass of psi: %e\n", mass);
        printf("L2 Relative Error: %e\n", psiError);
        printf("Diff between groups: %e\n", diffGroups);
    }
}

}
