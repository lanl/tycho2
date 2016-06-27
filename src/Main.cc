/*
    Main.cc
    
    Beginning of the program.  Sets up everything and starts the solver.
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

#include "TychoMesh.hh"
#include "Solver.hh"
#include "Quadrature.hh"
#include "Comm.hh"
#include "KeyValueReader.hh"
#include "Global.hh"
#include "Assert.hh"
#include "Timer.hh"
#include <signal.h>
#include <execinfo.h>
#include <omp.h>
#include <unistd.h>
#include <vector>


using namespace std;


/*
    signalHandler
    
    Prints a backtrace.  Used for debugging.
*/
static 
void signalHandler(int sig)
{
    void *array[20];
    size_t size;
    
    // get void*'s for all entries on the stack
    size = backtrace(array, 20);
    
    // print out all the frames to stderr
    fprintf(stderr, "Error: signal %d:\n", sig);
    backtrace_symbols_fd(array, size, STDERR_FILENO);
    exit(1);
}


/*
    createSweepSchedule
    
    Splits all the angles into sets of angle groups.
    One angle group per OMP thread.
    Creates an array of SweepSchedules, one entry for each angle group.
*/
static
void createSweepSchedule()
{
    // SweepSchedule for each angle group
    g_spSweepSchedule = new SweepSchedule*[g_nAngleGroups];
    
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
        vector<UINT> angles(numAngles);
        for (UINT angle = 0; angle < numAngles; angle++) {
            angles[angle] = angle + angleBdryIndices[angleGroup];
        }
        g_spSweepSchedule[angleGroup] = 
            new SweepSchedule(angles, g_maxCellsPerStep, g_intraAngleP, 
                              g_interAngleP);
    }
}


/*
    readInput
*/
static
void readInput(const string &inputFileName)
{
    // Read data
    CKG_Utils::KeyValueReader kvr;
    kvr.readFile(inputFileName);
    
    
    // Reader only reads int and not UINT type
    int snOrder, iterMax, maxCellsPerStep, intraAngleP, interAngleP, nGroups;
    
    
    // Get data
    kvr.getInt("snOrder", snOrder);
    kvr.getInt("iterMax", iterMax);
    kvr.getDouble("errMax", g_errMax);
    kvr.getInt("maxCellsPerStep", maxCellsPerStep);
    kvr.getInt("intraAngleP", intraAngleP);
    kvr.getInt("interAngleP", interAngleP);
    kvr.getInt("nGroups", nGroups);
    
    g_snOrder = snOrder;
    g_iterMax = iterMax;
    g_maxCellsPerStep = maxCellsPerStep;
    g_intraAngleP = intraAngleP;
    g_interAngleP = interAngleP;
    g_nGroups = nGroups;
    
    
    string sweepType;
    kvr.getString("SweepType", sweepType);
    if (sweepType == "OriginalTycho1")
        g_sweepType = SweepType_OriginalTycho1;
    else if (sweepType == "OriginalTycho2")
        g_sweepType = SweepType_OriginalTycho2;
    else if (sweepType == "TraverseGraph")
        g_sweepType = SweepType_TraverseGraph;
    else if (sweepType == "PBJ")
        g_sweepType = SweepType_PBJ;
    else
        Insist(false, "Sweep type not recognized.");
}


/*
    main
    
    Start of the program.
*/
int main( int argc, char *argv[] )
{
    static const double sigmaTotal = 10.0;
    static const double sigmaScat = 5.0;
    
    
    // For Debugging (prints a backtrace)
    signal(SIGSEGV, signalHandler);
    signal(SIGABRT, signalHandler);
    
    
    // Init MPI
    int required = MPI_THREAD_SINGLE;
    int provided = MPI_THREAD_SINGLE;
    int mpiResult = MPI_Init_thread(&argc, &argv, required, &provided);
    Insist (mpiResult == MPI_SUCCESS, "MPI_Init failed.");
    Insist (required <= provided, "");
    
    
    // Check inputs
    if (Comm::rank() == 0) {
        if (argc != 3) {
            printf("Incorrect number of arguments\n");
            printf("Usage: ./sweep.x <.pmesh file> <input.deck file>\n\n\n");
            MPI_Abort(MPI_COMM_WORLD, 6);
        }
    }
    
    
    // Print initial stuff
    if(Comm::rank() == 0) {
        printf("\n\n--- Initiating test of parallel sweeps. ---\n");
        printf("sigmaTotal: %f   sigmaScat: %f\n", sigmaTotal, sigmaScat);
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
    
    
    // Input data.
    readInput(argv[2]);
        
    
    // Create quadrature
    g_quadrature = new Quadrature(g_snOrder);
    
    
    // Create Tycho mesh
    Comm::barrier();
    Timer meshTimer;
    if(Comm::rank() == 0) {
        printf("Create Tycho Mesh\n");
        meshTimer.start();
    }
    g_spTychoMesh = new TychoMesh(argv[1]);
    Comm::barrier();
    if(Comm::rank() == 0) {
        meshTimer.stop();
        printf("Create Tycho Mesh Done: %fs\n", meshTimer.wall_clock());
    }
    
    
    // Make Sweep schedule object.
    if (g_sweepType == SweepType_OriginalTycho1 || 
        g_sweepType == SweepType_OriginalTycho2)
    {
        Comm::barrier();
        if(Comm::rank() == 0)
            printf("Create sweep schedule.\n");
        createSweepSchedule();
        Comm::barrier();
        if(Comm::rank() == 0)
            printf("Create sweep schedule done.\n");
    }
    
    
    // Do source iterations.
    Solver::solve(sigmaTotal, sigmaScat, g_iterMax, g_errMax);

    
    // End program
    MPI_Finalize();
    return 0;
}

