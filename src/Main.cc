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
#include "Problem.hh"
#include "SourceIteration.hh"
#include "Util.hh"
#include "Quadrature.hh"
#include "Comm.hh"
#include "KeyValueReader.hh"
#include "GraphTraverser.hh"
#include "Global.hh"
#include "Assert.hh"
#include "Timer.hh"
#include "SweepData.hh"
#include "SweeperAbstract.hh"
#include "Sweeper.hh"
#include "SweeperTraverse.hh"
#include "SweeperPBJ.hh"
#include "SweeperSchur.hh"
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
    readInput
*/
static
void readInput(const string &inputFileName, 
               double &sigmaT1, double &sigmaS1,
               double &sigmaT2, double &sigmaS2)
{
    // Read data
    CKG_Utils::KeyValueReader kvr;
    kvr.readFile(inputFileName);
    
    
    // Reader only reads int and not UINT type
    int snOrder, iterMax, maxCellsPerStep, intraAngleP, interAngleP, nGroups;
    int ddIterMax;
    
    
    // Get data
    kvr.getInt("snOrder", snOrder);
    kvr.getInt("iterMax", iterMax);
    kvr.getDouble("errMax", g_errMax);
    kvr.getInt("maxCellsPerStep", maxCellsPerStep);
    kvr.getInt("intraAngleP", intraAngleP);
    kvr.getInt("interAngleP", interAngleP);
    kvr.getInt("nGroups", nGroups);
    kvr.getDouble("sigmaT1", sigmaT1);
    kvr.getDouble("sigmaS1", sigmaS1);
    kvr.getDouble("sigmaT2", sigmaT2);
    kvr.getDouble("sigmaS2", sigmaS2);
    kvr.getBool("OutputFile", g_outputFile);
    kvr.getString("OutputFilename", g_outputFilename);
    kvr.getInt("DD_IterMax", ddIterMax);
    kvr.getDouble("DD_ErrMax", g_ddErrMax);
    kvr.getBool("SourceIteration", g_useSourceIteration);
    kvr.getBool("OneSidedMPI", g_useOneSidedMPI);
       
    g_snOrder = snOrder;
    g_iterMax = iterMax;
    g_maxCellsPerStep = maxCellsPerStep;
    g_intraAngleP = intraAngleP;
    g_interAngleP = interAngleP;
    g_nGroups = nGroups;
    g_ddIterMax = ddIterMax;
    
    
    
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
    else if (sweepType == "PBJOuter")
        g_sweepType = SweepType_PBJOuter;
    else if (sweepType == "Schur")
        g_sweepType = SweepType_Schur;
    else if (sweepType == "SchurOuter")
        g_sweepType = SweepType_SchurOuter;
    else if (sweepType == "PBJSI")
        g_sweepType = SweepType_PBJSI;
    else if (sweepType == "SchurKrylov")
        g_sweepType = SweepType_SchurKrylov;
    else
        Insist(false, "Sweep type not recognized.");


    string gaussElimMethod;
    kvr.getString("GaussElim", gaussElimMethod);
    if(gaussElimMethod == "Original")
        g_gaussElim = GaussElim_Original;
    else if (gaussElimMethod == "NoPivot")
        g_gaussElim = GaussElim_NoPivot;
    else if (gaussElimMethod == "CramerGlu")
        g_gaussElim = GaussElim_CramerGlu;
    else if (gaussElimMethod == "CramerIntel")
        g_gaussElim = GaussElim_CramerIntel;

}


/*
    main
    
    Start of the program.
*/
int main(int argc, char *argv[])
{
    double sigmaT1, sigmaS1, sigmaT2, sigmaS2;

    
    // For Debugging (prints a backtrace)
    signal(SIGSEGV, signalHandler);
    signal(SIGABRT, signalHandler);
    
    
    // Init MPI
    int required = MPI_THREAD_SINGLE;
    int provided = MPI_THREAD_SINGLE;
    int mpiResult = MPI_Init_thread(&argc, &argv, required, &provided);
    Insist (mpiResult == MPI_SUCCESS, "MPI_Init failed.");
    Insist (required <= provided, "");


    // Startup Petsc
    #if USE_PETSC
    PetscInitialize(&argc, &argv, (char*)0, NULL);
    #endif


    // Input data.
    if (argc < 3) {
        if (Comm::rank() == 0) {
            printf("Usage: ./sweep.x <parallel mesh> <input deck>\n");
        }
        MPI_Finalize();
        return 0;
    }
    readInput(argv[2], sigmaT1, sigmaS1, sigmaT2, sigmaS2);
    

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
    Problem::createCrossSections(g_sigmaT, g_sigmaS, sigmaT1, sigmaS1, 
                                 sigmaT2, sigmaS2);
    
    
    // Setup sweeper
    SweeperAbstract *sweeper = NULL;
    switch (g_sweepType) {
        case SweepType_OriginalTycho1:
        case SweepType_OriginalTycho2:
            g_graphTraverserForward = NULL;
            sweeper = new Sweeper();
            break;
        case SweepType_TraverseGraph:
            g_graphTraverserForward = 
                new GraphTraverser(Direction_Forward, true, 
                                   SweepData::getDataSizeInBytes());
            sweeper = new SweeperTraverse();
            break;
        case SweepType_Schur:
            g_graphTraverserForward = 
                new GraphTraverser(Direction_Forward, false, 
                                   SweepData::getDataSizeInBytes());
            sweeper = new SweeperSchur();
            break;
        case SweepType_SchurOuter:
            g_graphTraverserForward = 
                new GraphTraverser(Direction_Forward, false, 
                                   SweepData::getDataSizeInBytes());
            sweeper = new SweeperSchurOuter();
            break;
        case SweepType_SchurKrylov:
            g_graphTraverserForward = 
                new GraphTraverser(Direction_Forward, false, 
                                   SweepData::getDataSizeInBytes());
            sweeper = new SweeperSchurKrylov();
            break;
        case SweepType_PBJ:
            g_graphTraverserForward = 
                new GraphTraverser(Direction_Forward, false, 
                                   SweepData::getDataSizeInBytes());
            sweeper = new SweeperPBJ();
            break;
        case SweepType_PBJOuter:
            g_graphTraverserForward = 
                new GraphTraverser(Direction_Forward, false, 
                                   SweepData::getDataSizeInBytes());
            sweeper = new SweeperPBJOuter();
            break;
        case SweepType_PBJSI:
            g_graphTraverserForward = 
                new GraphTraverser(Direction_Forward, false, 
                                   SweepData::getDataSizeInBytes());
            sweeper = new SweeperPBJSI();
            break;
        default:
            Insist(false, "Sweep type not recognized.");
            break;
    }

    
    // Solve
    Timer timer;
    timer.start();
    sweeper->solve();
    timer.stop();
    
    
    // Time total solve
    double clockTime = timer.wall_clock();
    Comm::gmax(clockTime);
    if(Comm::rank() == 0) {
        printf("Total time: %.2f\n", clockTime);
    }
    
    
    // Print tests 
    double psiError = Problem::hatL2Error(sweeper->getPsi());
    double diffGroups = Util::diffBetweenGroups(sweeper->getPsi());
    if(Comm::rank() == 0) {
        printf("L2 Relative Error: %e\n", psiError);
        printf("Diff between groups: %e\n", diffGroups);
    }


    // Output psi to file
    if(g_outputFile)
        sweeper->writePsiToFile(g_outputFilename);

    
    // End program
    #if USE_PETSC
    PetscFinalize();
    #endif

    MPI_Finalize();
    return 0;
}

