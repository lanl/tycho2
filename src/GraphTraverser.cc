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

#include "GraphTraverser.hh"
#include "Mat.hh"
#include "Global.hh"
#include "TychoMesh.hh"
#include "Comm.hh"
#include "Timer.hh"
#include "Transport.hh"
#include <vector>
#include <queue>
#include <omp.h>
#include <Kokkos_Core.hpp>

using namespace std;


/*
    Tuple class
*/
namespace {
struct CellAnglePair
{
    UINT cell;
    UINT angle;
};
}


/*
    angleGroupIndex
    
    Gets angle groups for angle index.
    e.g. 20 angles numbered 0...19 with 3 threads.
    Split into 3 angle chunks of size 7,7,6:  0...6  7...13  14...19
    If angle in 0...6,   return 0
    If angle in 7...13,  return 1
    If angle in 14...19, return 2
*/
UINT angleGroupIndex(UINT angle)
{
    UINT numAngles = g_nAngles;
    UINT chunkSize = numAngles / g_nThreads;
    UINT numChunksBigger = numAngles % g_nThreads;
    UINT lowIndex = 0;
    
    
    // Find angleGroup
    for (UINT angleGroup = 0; angleGroup < g_nThreads; angleGroup++) {
        
        UINT nextLowIndex = lowIndex + chunkSize;
        if (angleGroup < numChunksBigger)
            nextLowIndex++;
        
        if (angle < nextLowIndex) 
            return angleGroup;
        
        lowIndex = nextLowIndex;
    }
    
    
    // Should never get here
    Assert(false);
    return 0;
}


/*
    GraphTraverser
*/
GraphTraverser::GraphTraverser(
        PsiData &psi, 
        const PsiData &source, 
        PsiBoundData &psiBound)
: c_initNumDependencies(g_nAngles, g_nCells),
  c_psi(psi), 
  c_psiBound(psiBound), 
  c_source(source)
{
    // Calc num dependencies for each (cell, angle) pair
    for (UINT cell = 0; cell < g_nCells; cell++) {
    for (UINT angle = 0; angle < g_nAngles; angle++) {
        
        c_initNumDependencies(angle, cell) = 0;
        for (UINT face = 0; face < g_nFacePerCell; face++) {
            
            bool incoming = g_tychoMesh->isIncoming(angle, cell, face);
            UINT adjCell = g_tychoMesh->getAdjCell(cell, face);
            
            if (incoming && adjCell != TychoMesh::BOUNDARY_FACE) {
                c_initNumDependencies(angle, cell)++;
            }
        }
    }}
}


/*
    update
    
    Does a transport update for the given cell/angle pair.
*/
void GraphTraverser::update(UINT cell, UINT angle) 
{
    double localSource[g_nVrtxPerCell][g_nMaxGroups];
    double localPsi[g_nVrtxPerCell][g_nMaxGroups];
    double localPsiBound[g_nVrtxPerFace][g_nFacePerCell][g_nMaxGroups];

    
    // Populate localSource
    #pragma omp simd
    for (UINT group = 0; group < g_nGroups; group++) {
    for (UINT vrtx = 0; vrtx < g_nVrtxPerCell; vrtx++) {
        localSource[vrtx][group] = c_source(group, vrtx, angle, cell);
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
        c_psi(group, vrtx, angle, cell) = localPsi[vrtx][group];
    }}
}


/*
    traverse
    
    Traverses g_tychoMesh.
*/
void GraphTraverser::traverse()
{
    Timer totalTimer;
    Timer setupTimer;
    

    // Start total timer
    totalTimer.start();
    

    // Kokkos Version
    if (g_useKokkos) {
        
        // Start setup timer
        setupTimer.start();


        // Get dependencies
        using space = Kokkos::DefaultExecutionSpace;
        auto nitems = int(g_nCells * g_nAngles);
        Kokkos::View<int*> counts("counts", nitems);
        Kokkos::parallel_for(nitems, KOKKOS_LAMBDA(int item) {
            int cell = item % g_nCells;
            int angle = item / g_nCells;
            for (UINT face = 0; face < g_nFacePerCell; ++face) {
                bool is_out = g_tychoMesh->isOutgoing(angle, cell, face);
                if (!is_out) continue;
                auto adjCell = g_tychoMesh->getAdjCell(cell, face);
                if (adjCell == TychoMesh::BOUNDARY_FACE) continue;
                counts(item) += 1;
            }
        });

        
        // Get the row_map
        Kokkos::View<int*> row_map;
        Kokkos::get_crs_row_map_from_counts(row_map, counts);
        auto nedges = row_map(row_map.size() - 1);
        Kokkos::View<int*> entries("entries", nedges);
        Kokkos::parallel_for(nitems, KOKKOS_LAMBDA(int item) {
            int cell = item % g_nCells;
            int angle = item / g_nCells;
            int j = 0;
            for (UINT face = 0; face < g_nFacePerCell; ++face) {
                bool is_out = g_tychoMesh->isOutgoing(angle, cell, face);
                if (!is_out) continue;
                auto adjCell = g_tychoMesh->getAdjCell(cell, face);
                if (adjCell == TychoMesh::BOUNDARY_FACE) continue;
                entries(row_map(item) + j) = adjCell + g_nCells * angle;
                ++j;
            }
            Assert(j + row_map(item) == row_map(item + 1));
        });


        // Get the policy
        auto graph = Kokkos::Crs<int,space,void,int>(row_map, entries);
        auto policy = Kokkos::WorkGraphPolicy<space,int>(graph);

        
        // End setup timer
        setupTimer.stop();


        // Traverse DAG
        auto lambda = KOKKOS_LAMBDA(int item) {
            int cell = item % g_nCells;
            int angle = item / g_nCells;
          
            // Update data for this cell-angle pair
            update(cell, angle);
        };
        Kokkos::parallel_for(policy, lambda);
    }

    // Non Kokkos version
    else {
        
        // Variables for non-kokkos version
        vector<queue<CellAnglePair>> canCompute(g_nThreads);
        Mat2<UINT> numDependencies(g_nAngles, g_nCells);

        // Start setup timer
        setupTimer.start();


        // Calc num dependencies for each (cell, angle) pair
        for (UINT cell = 0; cell < g_nCells; cell++) {
        for (UINT angle = 0; angle < g_nAngles; angle++) {
            numDependencies(angle, cell) = c_initNumDependencies(angle, cell);
        }}
        
        
        // Initialize canCompute queue
        for (UINT cell = 0; cell < g_nCells; cell++) {
        for (UINT angle = 0; angle < g_nAngles; angle++) {
            if (numDependencies(angle, cell) == 0) {
                UINT angleGroup = angleGroupIndex(angle);
                CellAnglePair cellAnglePair{cell, angle};
                canCompute[angleGroup].push(cellAnglePair);
            }
        }}


        // End setup timer
        setupTimer.stop();
        
        
        // Do local traversal
        #pragma omp parallel
        {
            UINT angleGroup = omp_get_thread_num();
            while (canCompute[angleGroup].size() > 0)
            {
                // Get cell/angle pair to compute
                CellAnglePair cellAnglePair = canCompute[angleGroup].front();
                canCompute[angleGroup].pop();
                UINT cell = cellAnglePair.cell;
                UINT angle = cellAnglePair.angle;
                
                
                // Update data for this cell-angle pair
                update(cell, angle);
                
                
                // Update dependency for children
                for (UINT face = 0; face < g_nFacePerCell; face++) {
                    
                    if (g_tychoMesh->isOutgoing(angle, cell, face)) {

                        UINT adjCell = g_tychoMesh->getAdjCell(cell, face);
                        
                        if (adjCell != TychoMesh::BOUNDARY_FACE) {
                            numDependencies(angle, adjCell)--;
                            if (numDependencies(angle, adjCell) == 0) {
                                CellAnglePair cellAnglePair{adjCell, angle};
                                canCompute[angleGroup].push(cellAnglePair);
                            }
                        }
                    }
                }
            }
        }
    }
    
    
    // Print times
    totalTimer.stop();

    double totalTime = totalTimer.wall_clock();
    Comm::gmax(totalTime);

    double setupTime = setupTimer.wall_clock();
    Comm::gmax(setupTime);

    if (Comm::rank() == 0) {
        printf("      Traverse Timer (setup):   %fs\n", setupTime);
        printf("      Traverse Timer (total):   %fs\n", totalTime);
    }
}

