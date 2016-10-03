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
#include "PsiData.hh"
#include "Comm.hh"
#include "SweepData.hh"
#include "CommSides.hh"
#include "SourceIteration.hh"
#include <math.h>
#include <limits>

using namespace std;



////////////////////////////////////////////////////////////////////////////////
//            SweeperPBJOuter functions
////////////////////////////////////////////////////////////////////////////////

/*
    solve
*/
void SweeperPBJOuter::solve()
{
    PsiData psi0;
    vector<UINT> sourceIts;
    
    
    // Source iterate till converged
    UINT iter = 1;
    while (iter < g_ddIterMax) {
        
        // Sweep
        UINT its = SourceIteration::solve(this, c_psi, c_source);
        sourceIts.push_back(its);
        

        // Check tolerance and set psi0 = psi1
        double errL1 = 0.0;
        double normL1 = 0.0;
        for (UINT i = 0; i < c_psi.size(); i++) {
            errL1  += fabs(psi0[i] - c_psi[i]);
            normL1 += fabs(c_psi[i]);
            
            psi0[i] = c_psi[i];
        }
        
        Comm::gsum(errL1);
        Comm::gsum(normL1);

        if (errL1 / normL1 < g_ddErrMax)
            break;
        
        
        // Communicate
        c_commSides.commSides(c_sweepData);
        
        
        // Increment iter
        iter++;
    }
    
    
    // Print statistics
    if (Comm::rank() == 0) {
        printf("PBJ Iters: %" PRIu64 "\n", iter);
        printf("Num source iterations:");
        for (UINT i = 0; i < sourceIts.size(); i++)
            printf(" %" PRIu64, sourceIts[i]);
        printf("\n");
    }
}


/*
    sweep
*/
void SweeperPBJOuter::sweep(PsiData &psi, const PsiData &source)
{
    UNUSED_VARIABLE(psi);
    UNUSED_VARIABLE(source);
    
    const bool doComm = false;
    const UINT maxComputePerStep = std::numeric_limits<uint64_t>::max();

    
    // Do 1 graph traversal
    traverseGraph(maxComputePerStep, c_sweepData, doComm, MPI_COMM_WORLD,
                  Direction_Forward);
}


////////////////////////////////////////////////////////////////////////////////
//            SweeperPBJ functions
////////////////////////////////////////////////////////////////////////////////

/*
    solve
*/
void SweeperPBJ::solve()
{
    SourceIteration::solve(this, c_psi, c_source);
}


/*
    sweep
*/
void SweeperPBJ::sweep(PsiData &psi, const PsiData &source)
{
    const bool doComm = false;
    const UINT maxComputePerStep = std::numeric_limits<uint64_t>::max();
    
    
    // Set psi0
    PsiData psi0;
    for (UINT i = 0; i < psi.size(); i++) {
        psi0[i] = psi[i];
    }
    
    
    // Create SweepData for traversal
    // Use a dummy set of priorities
    Mat2<UINT> priorities(g_nCells, g_nAngles);
    SweepData sweepData(psi, source, c_psiBoundPrev, g_sigmaTotal, priorities);
    //SweepData sweepData(psi, source, g_sigmaTotal, priorities);
    
    
    // Sweep till converged
    UINT iter = 1;
    while (iter < g_ddIterMax) {
        
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

        if (errL1 / normL1 < g_ddErrMax)
            break;
        
        
        // Communicate
        c_commSides.commSides(sweepData);
        
        
        // Increment iter
        iter++;
    }


    PsiBoundData &psiBound = sweepData.getSideData();
    for (UINT i = 0; i < psiBound.size(); i++) {
        c_psiBoundPrev[i] = psiBound[i];
    }
    
    
    // Print statistics
    if (Comm::rank() == 0) {
        printf("      PBJ Iters: %" PRIu64 "\n", iter);
    }
}


