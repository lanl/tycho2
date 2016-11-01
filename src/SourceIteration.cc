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

#include "SourceIteration.hh"
#include "Global.hh"
#include "PsiData.hh"
#include "Quadrature.hh"
#include "Comm.hh"
#include "Timer.hh"
#include <math.h>
#include <vector>
#include <petscmat.h>
#include <petscvec.h>
#include <petscksp.h>


using namespace std;


// cubeSize assumes using cube meshes in util folder
static const double cubeSize = 100.0;


/*
    hatSource
*/
static
void hatSource(const double sigmaT, const double sigmaS, PsiData &source)
{
    for(UINT cell = 0; cell < g_nCells; cell++) {
    for(UINT angle = 0; angle < g_nAngles; angle++) {
    for(UINT vrtx = 0; vrtx < g_nVrtxPerCell; vrtx++) {
        
        UINT node = g_tychoMesh->getCellNode(cell, vrtx);
        double x = g_tychoMesh->getNodeCoord(node, 0) - cubeSize / 2.0;
        double y = g_tychoMesh->getNodeCoord(node, 1) - cubeSize / 2.0;
        double z = g_tychoMesh->getNodeCoord(node, 2) - cubeSize / 2.0;
        
        double xi  = g_quadrature->getXi(angle);
        double eta = g_quadrature->getEta(angle);
        double mu  = g_quadrature->getMu(angle);
        
        double c = sqrt(x*x + y*y + z*z);
        
        for(UINT group = 0; group < g_nGroups; group++) {
            if(c <= 30.0) {
                source(group, vrtx, angle, cell) = 
                    - x / (30.0*c) * xi - y / (30.0*c) * eta - z / (30.0*c) * mu
                    + (sigmaT - sigmaS) * (1.0 - c / 30.0);
            }
            else {
                source(group, vrtx, angle, cell) = 0.0;
            }
        }
    }}}
}


/*
    hatL2Error
*/
static
double hatL2Error(const PsiData &psi)
{
    double diff = 0.0;
    double norm = 0.0;
    
    
    for(UINT cell = 0; cell < g_nCells; cell++) {
    for(UINT angle = 0; angle < g_nAngles; angle++) {
    for(UINT group = 0; group < g_nGroups; group++) {
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        double psiVal = 0.0;
        
        for(UINT vrtx = 0; vrtx < g_nVrtxPerCell; vrtx++) {
            UINT node = g_tychoMesh->getCellNode(cell, vrtx);
            x += g_tychoMesh->getNodeCoord(node, 0) - cubeSize / 2.0;
            y += g_tychoMesh->getNodeCoord(node, 1) - cubeSize / 2.0;
            z += g_tychoMesh->getNodeCoord(node, 2) - cubeSize / 2.0;
            psiVal += psi(group, vrtx, angle, cell);
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
    
    for(UINT cell = 0; cell < g_nCells; cell++) {
    for(UINT angle = 0; angle < g_nAngles; angle++) {
    for(UINT vrtx = 0; vrtx < g_nVrtxPerCell; vrtx++) {
        
        double psi0 = psi(0, vrtx, angle, cell);
        
        if(fabs(psi0) > maxEntry)
            maxEntry = psi0;
        
        for(UINT group = 1; group < g_nGroups; group++) {
            double psi1 = psi(group, vrtx, angle, cell);
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
    
    for(UINT cell = 0; cell < g_nCells; cell++) {
        double sumVrtxMass = 0.0;
        for(UINT vrtx = 0; vrtx < g_nVrtxPerCell; vrtx++) {
            double localMass = 0.0;
            for(UINT angle = 0; angle < g_nAngles; angle++) {
                localMass += psi(0, vrtx, angle, cell) * g_quadrature->getWt(angle);
            }
            sumVrtxMass += localMass;
        }
        mass += sumVrtxMass / g_nVrtxPerCell * g_tychoMesh->getCellVolume(cell);
    }
    
    Comm::gsum(mass);
    return mass;
}


/*
    calcMass
*/
/*static
double calcMass(const PhiData &phi)
{
    double mass = 0.0;
    
    for(UINT cell = 0; cell < g_nCells; cell++) {
        double sumVrtxMass = 0.0;
        for(UINT vrtx = 0; vrtx < g_nVrtxPerCell; vrtx++) {
            sumVrtxMass += phi(0, vrtx, cell);
        }
        mass += sumVrtxMass / g_nVrtxPerCell * g_tychoMesh->getCellVolume(cell);
    }
    
    Comm::gsum(mass);
    return mass;
}*/


/*
    psiToPhi
*/
static
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
static
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
static
void calcTotalSource(const PsiData &fixedSource, const PhiData &phiOld, 
                     PsiData &totalSource, bool onlyScatSource) 
{
    if (onlyScatSource) {
        #pragma omp parallel for
        for (UINT cell = 0; cell < g_nCells; ++cell) {
        for (UINT angle = 0; angle < g_nAngles; ++angle) {
        for (UINT vertex = 0; vertex < g_nVrtxPerCell; ++vertex) {
        for (UINT group = 0; group < g_nGroups; ++group) {
            totalSource(group, vertex, angle, cell) = 
                g_sigmaScat / (4.0 * M_PI) *  phiOld(group, vertex, cell);
        }}}}
    }
    else {
        #pragma omp parallel for
        for (UINT cell = 0; cell < g_nCells; ++cell) {
        for (UINT angle = 0; angle < g_nAngles; ++angle) {
        for (UINT vertex = 0; vertex < g_nVrtxPerCell; ++vertex) {
        for (UINT group = 0; group < g_nGroups; ++group) {
            totalSource(group, vertex, angle, cell) = 
                fixedSource(group, vertex, angle, cell) + 
                g_sigmaScat / (4.0 * M_PI) *  phiOld(group, vertex, cell);
        }}}}
    }
}


bool g_useZeroPsiBound = false;
// Global functions
namespace SourceIteration
{


/*
    getProblemSource
*/
void getProblemSource(PsiData &source)
{
    hatSource(g_sigmaTotal, g_sigmaScat, source);
}


/*
    Fixed point iteration (typical source iteration)
*/
UINT fixedPoint(SweeperAbstract *sweeper, PsiData &psi, const PsiData &source, 
                bool onlyScatSource)
{
    // Data for problem
    PsiData totalSource;
    PhiData phiNew;
    PhiData phiOld;
    
    
    // Get phi
    psiToPhi(phiNew, psi);
    
    
    // Volume of mesh
    double volume = 0.0;
    for(UINT cell = 0; cell < g_nCells; cell++) {
        volume += g_tychoMesh->getCellVolume(cell);
    }
    Comm::gsum(volume);
    if(Comm::rank() == 0) {
        printf("Volume of mesh: %e\n", volume);
    }


    // Source iteration
    UINT iter = 0;
    double error = 1.0;
    Timer innerTimer;
    innerTimer.start();
    while (iter < g_iterMax && error > g_errMax)
    {
        Timer timer, timer2, timer3, timer4;
        double wallClockTime = 0.0;
        double clockTime2 = 0.0;
        double clockTime3 = 0.0;
        timer.start();
        

        // phiOld = phiNew
        for(UINT i = 0; i < phiOld.size(); i++) {
            phiOld[i] = phiNew[i];
        }

        
        // totalSource = fixedSource + phiOld
        timer2.start();
        calcTotalSource(source, phiOld, totalSource, onlyScatSource);
        timer2.stop();
        
        clockTime2 = timer2.wall_clock();
        Comm::gmax(clockTime2);
        if(Comm::rank() == 0) {
            printf("\n   Source time: %f\n", clockTime2);
        }

        
        // Sweep
        sweeper->sweep(psi, totalSource);
        
        
        // Calculate L_inf relative error for phi
        timer3.start();
        psiToPhi(phiNew, psi);
        error = 0.0;
        for (UINT i = 0; i < phiNew.size(); i++) {
            error = max(error, fabs(phiNew[i] - phiOld[i]) / phiNew[i]);
        }
        Comm::gmax(error);
        timer3.stop();
        clockTime3 = timer3.wall_clock();
        Comm::gmax(clockTime3);
        if(Comm::rank() == 0) {
            printf("   psiToPhi Time: %f\n", clockTime3);
        }
        

        // Print iteration stats
        timer.stop();
        wallClockTime = timer.wall_clock();
        Comm::gmax(wallClockTime);
        if(Comm::rank() == 0) {
            printf("   iteration: %" PRIu64 "   error: %f   time: %f\n", 
                   iter, error, wallClockTime);
        }
        

        // Increment iteration
        ++iter;
    }
    
    
    // Print tests 
    double mass = calcMass(psi);
    double mass2 = calcMass(totalSource);
    double psiError = hatL2Error(psi);
    double diffGroups = diffBetweenGroups(psi);
    if(Comm::rank() == 0) {
        printf("Mass of psi: %e %e\n", mass, mass2);
        printf("L2 Relative Error: %e\n", psiError);
        printf("Diff between groups: %e\n", diffGroups);
    }


    // Time total solve
    innerTimer.stop();
    double clockTime = innerTimer.wall_clock();
    Comm::gmax(clockTime);
    if(Comm::rank() == 0) {
        printf("\nTotal source iteration time: %.2f\n",
               clockTime);
        printf("Solve time per iteration: %.2f\n\n",
               clockTime / iter);
    }


    // Return number of iterations
    return iter;
}


/*
    LHSData
    
    Data needed by the GMRES solver
*/
class LHSData
{
public:
    PsiData &c_psi;
    PsiData &c_source;
    SweeperAbstract &c_sweeper;

    LHSData(PsiData &psi, PsiData &source, SweeperAbstract &sweeper) :
    c_psi(psi), c_source(source), c_sweeper(sweeper)    
    {}
};


/*
    lhsOperator

    This performs the sweep and returns the boundary
*/
static
PetscErrorCode lhsOperator(Mat mat, Vec x, Vec b)
{
    // Get data for the solve
    void *dataVoid;
    MatShellGetContext(mat, &dataVoid);
    LHSData *data = (LHSData*) dataVoid;


    // Copy b to x
    VecCopy(x, b);

    // Get raw array from Petsc
    double *bArray;
    VecGetArray(b, &bArray);
    PhiData phi(bArray);
    
    // S operator
    for (UINT i = 0; i < phi.size(); i++) {
        phi[i] = g_sigmaScat / (4.0 * M_PI) * phi[i];
    }

    // M operator
    phiToPsi(phi, data->c_source);

    // L^-1 operator
    g_useZeroPsiBound = true;
    data->c_sweeper.sweep(data->c_psi, data->c_source);

    // D operator
    psiToPhi(phi, data->c_psi);

    // Give array back to Petsc
    VecRestoreArray(b, &bArray);

    // b = x - b
    VecAXPBY(b, 1, -1, x);

    
    return 0;
}


/*
    Krylov solver
*/
UINT krylov(SweeperAbstract *sweeper, PsiData &psi, const PsiData &source, 
            bool onlyScatSource)
{
    UNUSED_VARIABLE(onlyScatSource);
    
    KSP ksp;
    Mat mat;
    Vec x, b;
    PetscInt vecSize;
    PetscInt totalVecSize;
    int its;
    double rnorm;
    double *bArray;
    double *xArray;
    PsiData tempSource;


    // Local and global vector sizes
    vecSize = g_nCells * g_nVrtxPerCell * g_nGroups;
    UINT totalSize = vecSize;
    Comm::gsum(totalSize);
    totalVecSize = totalSize;


    // Create vectors
    VecCreate(MPI_COMM_WORLD, &x);
    VecSetSizes(x, vecSize, totalVecSize);
    VecSetType(x, VECMPI);
    VecDuplicate(x, &b);
    

    // Create matrix shell and define it as the operator
    MatCreateShell(MPI_COMM_WORLD, vecSize, vecSize, totalVecSize, totalVecSize,
                   (void*)(NULL), &mat);
    MatShellSetOperation(mat, MATOP_MULT, (void(*)(void))lhsOperator);
    

    // Set solver context
    KSPCreate(MPI_COMM_WORLD, &ksp);
    

    // Set operator to KSP context. No preconditioning will actually
    // be used due to the PCNONE option.
    KSPSetOperators(ksp, mat, mat);
    KSPSetTolerances(ksp, g_errMax, PETSC_DEFAULT, PETSC_DEFAULT, g_iterMax);


    // Setup RHS (b = D L^{-1} Q)
    if (Comm::rank() == 0)
        printf("Krylov source\n");
    sweeper->sweep(psi, source);

    VecGetArray(b, &bArray);
    PhiData phiB(bArray);
    psiToPhi(phiB, psi);
    VecRestoreArray(b, &bArray);


    // Set LHSData
    LHSData data(psi, tempSource, *sweeper);
    MatShellSetContext(mat, &data);


    // Solve
    if (Comm::rank() == 0)
        printf("Begin Solve\n");
    KSPSolve(ksp, b, x);
    if (Comm::rank() == 0)
        printf("End Solve\n");


    // Get Psi from Phi
    // Psi = L^{-1} (MS Phi + Q)
    VecGetArray(x, &xArray);
    PhiData phiX(xArray);
    for (UINT i = 0; i < phiX.size(); i++) {
        phiX[i] = g_sigmaScat / (4.0 * M_PI) * phiX[i];
    }
    phiToPsi(phiX, psi);
    VecRestoreArray(x, &xArray);
    for (UINT i = 0; i < psi.size(); i++) {
        tempSource[i] = source[i] + psi[i];
    }

    sweeper->sweep(psi, tempSource);
    
    
    // Print some stats
    KSPGetIterationNumber(ksp, &its);
    KSPGetResidualNorm(ksp, &rnorm);
    if (Comm::rank() == 0) {
        printf("   Krylov iterations: %u with Rnorm: %e\n", its, rnorm);
    }


    // Destroy vectors, ksp, and matrices to free work space
    VecDestroy(&x);
    VecDestroy(&b);
    MatDestroy(&mat);
    KSPDestroy(&ksp);
    
    
    // Return number of iterations
    return its;
}
}//End namespace SourceIteration


