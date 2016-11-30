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
#include "Comm.hh"
#include "Timer.hh"
#include "Util.hh"
#include "KrylovSolver.hh"
#include <math.h>


namespace
{

/*
    LHSData
    
    Data needed by the Krylov solver
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

    Performs b = (I - D L^{-1} M S) x
*/
void lhsOperator(const double *x, double *b, void *voidData)
{
    // Get data for the solve
    LHSData *data = (LHSData*) voidData;
    UINT vecSize = g_nCells * g_nVrtxPerCell * g_nGroups;
    const bool zeroPsiBound = true;


    // b = x
    for (UINT i = 0; i < vecSize; i++) {
        b[i] = x[i];
    }


    // Get raw array from Petsc
    PhiData phi(b);
    

    // S operator
    Util::operatorS(phi, phi);


    // M operator
    Util::phiToPsi(phi, data->c_source);


    // L^-1 operator
    data->c_sweeper.sweep(data->c_psi, data->c_source, zeroPsiBound);


    // D operator
    Util::psiToPhi(phi, data->c_psi);


    // b = x - b
    for (UINT i = 0; i < vecSize; i++) {
        b[i] = x[i] - b[i];
    }
}

} // End anonymous namespace


// Global functions
namespace SourceIteration
{

/*
    Fixed point iteration (typically called source iteration)
    L Psi^{n+1} = MS \Phi^n + Q
*/
UINT fixedPoint(SweeperAbstract &sweeper, PsiData &psi, const PsiData &source)
{
    // Data for problem
    PsiData totalSource;
    PhiData phiNew;
    PhiData phiOld;
    
    
    // Get phi
    Util::psiToPhi(phiNew, psi);
    
    
    // Source iteration
    UINT iter = 0;
    double error = 1.0;
    Timer totalTimer;
    totalTimer.start();
    while (iter < g_iterMax && error > g_errMax)
    {
        Timer timer;
        double wallClockTime = 0.0;
        double norm = 0.0;
        timer.start();
        

        // phiOld = phiNew
        for(UINT i = 0; i < phiOld.size(); i++) {
            phiOld[i] = phiNew[i];
        }

        
        // totalSource = fixedSource + phiOld
        Util::calcTotalSource(source, phiOld, totalSource);

        
        // Sweep
        sweeper.sweep(psi, totalSource);
        
        
        // Calculate L_1 relative error for phi
        Util::psiToPhi(phiNew, psi);
        error = 0.0;
        for (UINT i = 0; i < phiNew.size(); i++) {
            error += fabs(phiNew[i] - phiOld[i]);
            norm += fabs(phiNew[i]);
        }
        Comm::gsum(error);
        Comm::gsum(norm);
        error = error / norm;
        

        // Print iteration stats
        timer.stop();
        wallClockTime = timer.wall_clock();
        Comm::gmax(wallClockTime);
        if(Comm::rank() == 0) {
            printf("   iteration: %" PRIu64 "   error: %e   time: %f\n", 
                   iter, error, wallClockTime);
        }
        

        // Increment iteration
        ++iter;
    }
    
    
    // Time total solve
    totalTimer.stop();
    double clockTime = totalTimer.wall_clock();
    Comm::gmax(clockTime);
    if(Comm::rank() == 0) {
        printf("\nTotal source iteration time: %.2f\n",
               clockTime);
        printf("Average source iteration time: %.2f\n\n",
               clockTime / iter);
    }


    // Return number of iterations
    return iter;
}


/*
    Krylov solver

    Solves (I - DL^{-1}MS) \Phi = DL^{-1} Q.
    L could be the full sweep or a local sweep L_I.
*/
UINT krylov(SweeperAbstract &sweeper, PsiData &psi, const PsiData &source)
{
    UINT vecSize;
    int its;
    double rnorm;
    double *bArray;
    double *xArray;
    PsiData tempSource;
    Timer totalTimer;
    totalTimer.start();


    // Create the Krylov solver
    vecSize = g_nCells * g_nVrtxPerCell * g_nGroups;
    LHSData data(psi, tempSource, sweeper);
    KrylovSolver krylovSolver(vecSize, g_errMax, g_iterMax, lhsOperator);
    krylovSolver.setData(&data);


    // Setup RHS (b = D L^{-1} Q)
    if (Comm::rank() == 0)
        printf("Krylov source\n");
    sweeper.sweep(psi, source);

    bArray = krylovSolver.getB();
    PhiData phiB(bArray);
    Util::psiToPhi(phiB, psi);
    krylovSolver.releaseB();


    // Solve
    if (Comm::rank() == 0)
        printf("Krylov Begin Solve\n");
    krylovSolver.solve();
    if (Comm::rank() == 0)
        printf("Krylov End Solve\n");


    // Get Psi from Phi
    // Psi = L^{-1} (MS Phi + Q)
    xArray = krylovSolver.getX();
    PhiData phiX(xArray);
    Util::calcTotalSource(source, phiX, tempSource);
    krylovSolver.releaseX();
    sweeper.sweep(psi, tempSource);
    
    
    // Print some stats
    its = krylovSolver.getNumIterations();
    rnorm = krylovSolver.getResidualNorm();
    if (Comm::rank() == 0) {
        printf("Krylov iterations: %u with Rnorm: %e\n", its, rnorm);
    }


    // Time total solve
    totalTimer.stop();
    double clockTime = totalTimer.wall_clock();
    Comm::gmax(clockTime);
    if(Comm::rank() == 0) {
        printf("\nTotal Krylov time: %.2f\n",
               clockTime);
        printf("Average Krylov time: %.2f\n\n",
               clockTime / its);
    }


    // Return number of iterations
    return its;
}
}//End namespace SourceIteration


