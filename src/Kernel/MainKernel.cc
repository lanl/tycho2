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

#include "TychoMesh.hh"
#include "Quadrature.hh"
#include "Comm.hh"
#include "Global.hh"
#include "Assert.hh"
#include "Timer.hh"
#include "SweeperKernel.hh"
#include "GraphTraverserKernel.hh"
#include <omp.h>


using namespace std;


/*
    setInput
*/
static
void setInput(double &sigmaT1, double &sigmaS1,
              double &sigmaT2, double &sigmaS2)
{
    g_snOrder = 8;
    g_iterMax = 10;
    g_errMax = 1e-4;
    g_maxCellsPerStep = 100;
    g_intraAngleP = 3;
    g_interAngleP = 1;
    g_nGroups = 1;
    sigmaT1 = 10.0;
    sigmaS1 = 5.0;
    sigmaT2 = 10.0;
    sigmaS2 = 5.0;
    g_outputFile = false;
    g_outputFilename = "nothing.psi";
    g_ddIterMax = 10;
    g_ddErrMax = 1e-5;
    g_useSourceIteration = true;
    g_mpiType = MPIType_CapsaicinTwoSided;
    g_sweepType = SweepType_PBJSI;
    g_gaussElim = GaussElim_NoPivot;
}


/*
    main
    
    Start of the program.
*/
int main(int argc, char *argv[])
{
    double sigmaT1, sigmaS1, sigmaT2, sigmaS2;

    
    // Init MPI
    int required = MPI_THREAD_SINGLE;
    int provided = MPI_THREAD_SINGLE;
    int mpiResult = MPI_Init_thread(&argc, &argv, required, &provided);
    Insist (mpiResult == MPI_SUCCESS, "MPI_Init failed.");
    Insist (required <= provided, "");


    // Input data.
    if (argc < 2) {
        if (Comm::rank() == 0) {
            printf("Usage: ./sweep.x <parallel mesh>\n");
        }
        MPI_Finalize();
        return 0;
    }
    setInput(sigmaT1, sigmaS1, sigmaT2, sigmaS2);
    

    // Print initial stuff
    if(Comm::rank() == 0) {
        printf("\n\n--- Initiating test of parallel sweeps. ---\n");
        printf("sigmaT1: %lf   sigmaS1: %lf\n", sigmaT1, sigmaS1);
        printf("sigmaT2: %lf   sigmaS2: %lf\n", sigmaT2, sigmaS2);
        printf("ASSERT_ON: %d\n", ASSERT_ON);
    }
    
    
    // Get number of angle groups
    // This is the same as the number of OpenMP threads
    g_nAngleGroups = 1;
    #pragma omp parallel
    {
        if(omp_get_thread_num() == 0)
            g_nAngleGroups = omp_get_num_threads();
    }
    g_nThreads = g_nAngleGroups;
    if (Comm::rank() == 0)
        printf("Num angle groups: %" PRIu64 "\n", g_nAngleGroups);
            
    
    // Create quadrature
    g_quadrature = new Quadrature(g_snOrder);
    
    
    // Create Tycho mesh
    Comm::barrier();
    Timer meshTimer;
    if(Comm::rank() == 0) {
        printf("Create Tycho Mesh\n");
        meshTimer.start();
    }
    g_tychoMesh = new TychoMesh(argv[1]);
    Comm::barrier();
    if(Comm::rank() == 0) {
        meshTimer.stop();
        printf("Create Tycho Mesh Done: %fs\n", meshTimer.wall_clock());
    }


    // Create cross sections for each cell
    g_sigmaT.resize(g_nCells, sigmaT1);
    g_sigmaS.resize(g_nCells, sigmaS1);
    
    
    // Setup sweeper
    SweeperKernel sweeper;
    g_graphTraverserKernel = new GraphTraverserKernel();

    
    // Solve
    Timer timer;
    timer.start();
    sweeper.solve();
    timer.stop();
    
    
    // Time total solve
    double clockTime = timer.wall_clock();
    Comm::gmax(clockTime);
    if(Comm::rank() == 0) {
        printf("Total time: %.2f\n", clockTime);
    }
    
    
    MPI_Finalize();
    return 0;
}

