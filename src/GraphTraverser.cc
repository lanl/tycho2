
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
#include "Transport.cc.hh"
#include <nvgraph.h>
#include <algorithm>

/*
    GraphTraverser
*/
GraphTraverser::GraphTraverser()
{
}

/*
    traverseNV
    
    Traverses g_tychoMesh using nvgraph.
*/

void check_status(nvgraphStatus_t status){
  if((int) status !=0) {
    printf("nvGraph ERROR: %d\n", status);
    exit(0);
  }
}

void GraphTraverser::traverseNV(
    const PsiData &source,
    const PsiBoundData &psiBound,
    PsiData &psi)
{
    Timer totalTimer;
    Timer setupTimer;

    // Start total timer
    totalTimer.start();

    // Start setup timer
    setupTimer.start();

    // Get dependencies
    int nitems = int(g_nCells * g_nAngles);

    int counts[nitems];
    for (int i = 0; i < nitems; ++i) {
      counts[i] = 0;
    }

    // Get row_map from counts (number of outgoing cell faces (Omega . n > 0) 
    //   along a given angle in mesh)
    int sum = 0;
    for(int item = 0; item  < nitems; ++item) {
        int cell = item % g_nCells;
        int angle = item / g_nCells;
        for (UINT face = 0; face < g_nFacePerCell; ++face) {
            bool is_out = g_tychoMesh->isOutgoing(angle, cell, face);
            if (!is_out) continue;
            auto adjCell = g_tychoMesh->getAdjCell(cell, face);
            if (adjCell == TychoMesh::BOUNDARY_FACE) continue;
            counts[item] += 1;
        }
        sum += counts[item];
    }

    int nrows = nitems + 1;
    int row_map[nrows];

    row_map[0] = 0;

    for (int j=1; j < nrows; ++j){
       row_map[j] = counts[j-1] + row_map[j-1];        
    }        

    // Get entries
    int nedges = row_map[nitems];
#if 0
    if (Comm::rank() == 0) {

      printf(" g_nCells: %lu, g_nAngles: %lu \n", g_nCells, g_nAngles);
      printf(" nitems: %d, nrows: %d, nedges: %d \n", nitems, nrows, nedges);
      printf(" sum: %d \n", sum);
      printf(" counts[%d]: %d \n", nitems-1, counts[nitems-1]);
    }
#endif

    int entries[nedges];

    for(int item = 0; item < nitems; ++item){
      int cell = item % g_nCells;
      int angle = item / g_nCells;
      int j = 0;
      for(UINT face = 0; face < g_nFacePerCell; ++face){
        bool is_out = g_tychoMesh->isOutgoing(angle, cell, face);
        if(!is_out)
          continue;
        auto adjCell = g_tychoMesh->getAdjCell(cell, face);
        if(adjCell == TychoMesh::BOUNDARY_FACE)
          continue;
        entries[row_map[item] + j] = adjCell + g_nCells * angle;
        ++j;
      }
      Assert(j + row_map[item]  == row_map[item + 1]);
    }
          
#if 0
    // Printing count, row_map and entries info
    if (Comm::rank() == 0) {

      for(int i = 0; i < nitems+1; i++){
        printf("  row_map[%d]: %d \n", i, row_map[i]);
      }
      printf("\n");

      for(int i = 0; i < nedges; i++){
        printf("   entries[%d]: %d \n", i, entries[i]);
      }
      printf("\n");

    }
#endif

    // End setup timer
    setupTimer.stop();

   // Begin nvgraph calculations

   // Store distances from source and where to store
   // predecessors in search tree

    int bfs_distances_h[nitems];
    //int bfs_predecessors_h[nitems];

    //const size_t vertex_numsets = 2;
    const size_t vertex_numsets = 1;
   
    nvgraphHandle_t handle;
    nvgraphGraphDescr_t graph;
    nvgraphCSRTopology32I_t CSR_input;
    cudaDataType_t* vertex_dimT;
    size_t distances_index = 0;
    //size_t predecessors_index = 1;
    vertex_dimT = (cudaDataType_t*)malloc(vertex_numsets*sizeof(cudaDataType_t));
    vertex_dimT[distances_index] = CUDA_R_32I;
    //vertex_dimT[predecessors_index] = CUDA_R_32I;

    //Creating nvgraph objects
    check_status(nvgraphCreate (&handle));
    check_status(nvgraphCreateGraphDescr (handle, &graph));

    // Set graph connectivity and properties (tranfers)
    CSR_input = (nvgraphCSRTopology32I_t) malloc(sizeof(struct nvgraphCSCTopology32I_st));
    CSR_input->nvertices = nitems;
    CSR_input->nedges = nedges;
    CSR_input->source_offsets = row_map;
    CSR_input->destination_indices = entries;

    check_status(nvgraphSetGraphStructure(handle, graph, (void*)CSR_input, NVGRAPH_CSR_32));
    check_status(nvgraphAllocateVertexData(handle, graph, vertex_numsets, vertex_dimT));

    nvgraphTraversalParameter_t traversal_param;
    nvgraphTraversalParameterInit(&traversal_param);
    nvgraphTraversalSetDistancesIndex(&traversal_param, distances_index);
    //nvgraphTraversalSetPredecessorsIndex(&traversal_param, predecessors_index);
    nvgraphTraversalSetUndirectedFlag(&traversal_param, false);
    
    // Setup queue for solution process which has a maximum possible number of bLevels equal
    // to maxBLevels, the largest possible distance to any leaf in the entire graph.
    // nvGraph stores these in the bfs_distances_h vector and are converted to Q values based 
    // on the bLevel (priority) and item number; i.e. Q = Q[maxBLevels][nitems], where 
    // nitems = g_nAngles * g_nCells. maxBLevels is estimated below, but may need increased based on
    // the number of cells in the problem.
    // 

    int maxBLevels = 100;
    int **Q = new int*[maxBLevels];

    // Initialize Q. Initially, assume all cells for each angle all
    // have priority (bLevel) 0, so preload Q[0][item] = cell.
    Q[0] = new int[nitems];

    for(int p = 1; p < maxBLevels; ++p){
      Q[p] = new int [nitems];
      for(UINT ang=0; ang < g_nAngles; ++ang){
         for(UINT cell = 0; cell < g_nCells; ++cell){
           int item = cell + ang*g_nCells;
              Q[p][item] = -1;
          }
       }
    }

    for(int item = 0; item < nitems; item++)
       Q[0][item] = item;
   
    // Do the graph traverse and do the updates 
    // 
    int actMaxBL = 0;

    for (int item = 0; item < nitems; item++){   // nvgraph traverses

      int source_vert = item;

      //printf(" source_vert: %d \n", source_vert);
        
      // Computing traversal using BFS algorithm
      check_status(nvgraphTraversal(handle, graph, NVGRAPH_TRAVERSAL_BFS, &source_vert, traversal_param));

      // Get result
     check_status(nvgraphGetVertexData(handle, graph, (void*)bfs_distances_h, distances_index));
     //check_status(nvgraphGetVertexData(handle, graph, (void*)bfs_predecessors_h, predecessors_index));

      //int cell = item % g_nCells;
      int angle = item / g_nCells;
#if 0
      UINT NN = 0;
      for(UINT i = 0; i < nitems; ++i){
        if(i / g_nCells == NN){
          printf(" angle: %lu \n", NN);
          NN++;
        }
        if(bfs_distances_h[i] != pow(2,31)-1)
          printf("Distance to vertex %d: %i\n",i, bfs_distances_h[i]);
      }
       printf("\n");
#endif
      // Fill the Q array with cells for a given bValue 
      int MAXBL = pow(2,31)-1;
      int bItem = angle*g_nCells;
      int eItem = bItem + g_nCells;
      //printf("bItem: %d, eItem: %d \n", bItem, eItem);

      for(int it = bItem; it < eItem; ++it){
         int idx = bfs_distances_h[it]; 
         if(idx == MAXBL)
           continue;
         else{
           Q[idx][it] = it;
           actMaxBL = std::max(actMaxBL, idx);
         }
      }
#if 0
      NN = 0;
      for(int i = 0; i < nitems; ++i){
        if(i / g_nCells == NN){
          printf(" angle: %d \n", NN);
          NN++;
        }
        //if(bfs_predecessors_h[i] > -1)
          printf("Predecessor of vertex %d: %i\n",i, bfs_predecessors_h[i]);
      }
      printf("\n");

#endif
     }  // End of nvgraph traverses

     //printf("actMaxBL: %d \n", actMaxBL);

      // Select highest bLevel in each queue
      for(int it = 0; it < nitems; it++){
         for(int p = actMaxBL+1; p-- > 0; ){
            //printf("Q[%d][%d]: %d \n", p, it, Q[p][it]);
            if(Q[p][it] != -1){
               for(int pp = p; pp-- > 0; ){
                 Q[pp][it] = -1;
               }
            }
          }
      }
#if 0
  // Print  Q values
  printf("\n \n");
  for(int p = 0; p < actMaxBL+1; ++p)
     for(int item = 0; item < nitems; ++item)
         printf("Q[%d][%d]: %d \n", p, item, Q[p][item]);
  printf("\n \n");
#endif  

    // Setup pointers to solution data
    PsiData* d_psi = &psi;
    const PsiData* d_source = &source;
    const PsiBoundData* d_psiBound = &psiBound;

    // Setup update method for queue
    auto lambda = [=](int item) {
        int cell = item % g_nCells;
        int angle = item / g_nCells;

        Transport::updateNVgraph(
            cell,
            angle,
            d_source,
            d_psiBound,
            d_psi,
            g_nGroups,
            g_tychoMesh->c_omegaDotN,
            g_tychoMesh->c_adjCell,
            g_tychoMesh->c_neighborVrtx,
            g_tychoMesh->c_adjProc,
            g_tychoMesh->c_side,
            g_tychoMesh->c_cellVolume,
            g_sigmaT,
            g_tychoMesh->c_faceArea,
            g_tychoMesh->c_cellToFaceVrtx);
    };
    
    // Perform sweep

    // These values were taken from the Kokkos version of the solution for a 12 Cell, 8 angle test
    //  mesh found in /util.  They were used as a comparison for the Q values used below.
    //int queue[nitems] = {1,12, 34,41,50,51,64,77,89,2,5,6,20,21,22,24,29,33,27,40,46,49,59,52,63,65,68,73,76,82,85,88,94,3,11,7,9,16,19,18,17,32,28,25,30,38,42,36,45,48,55,56,53,62,64,60,67,75,80,72,81,86,90,92,84,93,4,8,0,10,23,13,31,27,26,39,47,43,44, 54,58,61,69,70,76,72,83,79,78,95,14,15,35,57,91,4,8,0,10,23,13,31,27,26,39,47,43,44, 54,58,61,69,70,76,72,83,79,78,95,14,15,35,57,91,4,8,0,10,23,13,31,27,26,39,47,43,44, 54,58,61,69,70,76,72,83,79,78,95,14,15,35,57,91,4,8,0,10,23,13,31,27,26,39,47,43,44, 54,58,61,69,70,76,72,83,79,78,95,14,15,35,57,91};

    // Do sweep 
    for(int p = 0; p < actMaxBL+1; p++){
       for(UINT ang = 0; ang < g_nAngles; ang++){
          int bItem = ang*g_nCells;
          int eItem = bItem + g_nCells;
          for(int item = bItem; item < eItem; item++){
            if(Q[p][item] == -1)
              continue;
            int Item = Q[p][item];
            //printf("angle: %lu, Item: %d \n", ang, Item);
            lambda(Item);
          }
       }
   }

   for(int p = 0; p < maxBLevels; p++)
         delete [] Q[p]; 
   delete [] Q;
 
   free(vertex_dimT);
   free(CSR_input);

   check_status(nvgraphDestroyGraphDescr (handle, graph));
   check_status(nvgraphDestroy (handle));

   // Print times
   totalTimer.stop();
   double totalTime = totalTimer.wall_clock();
   Comm::gmax(totalTime);
   //
   double setupTime = setupTimer.wall_clock();
     Comm::gmax(setupTime);
     if (Comm::rank() == 0) {
        printf("      Traverse Timer (setup):   %fs\n", setupTime);
        printf("      Traverse Timer (total):   %fs\n", totalTime);
    }
}
