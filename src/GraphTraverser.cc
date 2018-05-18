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
#include <Kokkos_Core.hpp>


/*
    GraphTraverser
*/
GraphTraverser::GraphTraverser()
{
  Kokkos::Profiling::pushRegion("GraphTraverser::GraphTraverser");
    // Put constant data on GPU
    auto host_omega_dot_n = host_mat3_t<double>(
        g_tychoMesh->c_omegaDotN.data(),
        g_nAngles,
        g_nCells,
        g_nFacePerCell);
    auto host_adj_cell = host_mat2_t<UINT>(
        g_tychoMesh->c_adjCell.data(),
        g_nCells,
        g_nFacePerCell);
    auto host_neighbor_vertex = host_mat3_t<UINT>(
        g_tychoMesh->c_neighborVrtx.data(),
        g_nCells,
        g_nFacePerCell,
        g_nVrtxPerFace);
    auto host_adj_proc = host_mat2_t<UINT>(
        g_tychoMesh->c_adjProc.data(),
        g_nCells,
        g_nFacePerCell);
    auto host_side = host_mat2_t<UINT>(
        g_tychoMesh->c_side.data(),
        g_nCells,
        g_nFacePerCell);
    auto host_sigma_t = host_mat1_t<double>(
        g_sigmaT.data(),
        g_sigmaT.size());
    auto host_cell_volume = host_mat1_t<double>(
        g_tychoMesh->c_cellVolume.data(),
        g_nCells);
    auto host_face_area = host_mat2_t<double>(
        g_tychoMesh->c_faceArea.data(),
        g_nCells,
        g_nFacePerCell);
    auto host_cell_to_face_vertex = host_mat3_t<UINT>(
        g_tychoMesh->c_cellToFaceVrtx.data(),
        g_nCells,
        g_nFacePerCell,
        g_nVrtxPerCell);

    c_device_omega_dot_n = Kokkos::create_mirror_view_and_copy(
        device(), host_omega_dot_n);
    c_device_adj_cell = Kokkos::create_mirror_view_and_copy(
        device(), host_adj_cell);
    c_device_neighbor_vertex = Kokkos::create_mirror_view_and_copy(
        device(), host_neighbor_vertex);
    c_device_adj_proc = Kokkos::create_mirror_view_and_copy(
        device(), host_adj_proc);
    c_device_side = Kokkos::create_mirror_view_and_copy(
        device(), host_side);
    c_device_sigma_t = Kokkos::create_mirror_view_and_copy(
        device(), host_sigma_t);
    c_device_cell_volume = Kokkos::create_mirror_view_and_copy(
        device(), host_cell_volume);
    c_device_face_area = Kokkos::create_mirror_view_and_copy(
        device(), host_face_area);
    c_device_cell_to_face_vertex = Kokkos::create_mirror_view_and_copy(
        device(), host_cell_to_face_vertex);

    // Allocate space on GPU
    // NOT SAFE?
    c_device_source = device_psi_data_t(
        "source",
        g_nGroups,
        g_nVrtxPerCell,
        g_nAngles,
        g_nCells);
    c_device_psi_bound = device_psi_data_t(
        "psi_bound",
        g_nGroups,
        g_nVrtxPerFace,
        g_nAngles,
        g_tychoMesh->getNSides());
    c_device_psi = device_psi_data_t(
        "psi",
        g_nGroups,
        g_nVrtxPerCell,
        g_nAngles,
        g_nCells);

    // Get dependencies
    auto nitems = int(g_nCells * g_nAngles);
    Kokkos::View<int*, host> counts("counts", nitems);
    Kokkos::parallel_for(Kokkos::RangePolicy<int, host>(0,nitems), [=](int item) {
        int cell = item % g_nCells;
        int angle = item / g_nCells;
        for (UINT face = 0; face < g_nFacePerCell; ++face) {
            bool is_out = g_tychoMesh->isOutgoing(angle, cell, face);
            if (!is_out) continue;
            auto adjCell = g_tychoMesh->getAdjCell(cell, face);
            if (adjCell == TychoMesh::BOUNDARY_FACE) continue;
            counts(item) += 1;
        }
    }, "count-dependencies");

    
    // Get the row_map
    Kokkos::View<int*, host> row_map;
    Kokkos::get_crs_row_map_from_counts(row_map, counts);
    auto nedges = row_map(row_map.size() - 1);
    Kokkos::View<int*, host> entries("entries", nedges);
    Kokkos::parallel_for(Kokkos::RangePolicy<int, host>(0,nitems), [=](int item) {
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
    }, "fill-dependencies");

    // Get the policy
    auto device_row_map = Kokkos::create_mirror_view_and_copy(device(), row_map);
    auto device_entries = Kokkos::create_mirror_view_and_copy(device(), entries);
    c_device_graph = Kokkos::Crs<int,device,void,int>(device_row_map, device_entries);
    Kokkos::Profiling::popRegion();
}


/*
    traverse
    
    Traverses g_tychoMesh.
*/
void GraphTraverser::traverse(
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

    auto policy = Kokkos::WorkGraphPolicy<device,int>(c_device_graph);

    // End setup timer
    setupTimer.stop();

    // Copy data to device views
    auto host_source = host_psi_data_t(
        const_cast<PsiData&>(source).data(),
        g_nGroups,
        g_nVrtxPerCell,
        g_nAngles,
        g_nCells);
    auto host_psi_bound = host_psi_data_t(
        const_cast<PsiBoundData&>(psiBound).data(),
        g_nGroups,
        g_nVrtxPerFace,
        g_nAngles,
        g_tychoMesh->getNSides());
    auto host_psi = host_psi_data_t(
        psi.data(),
        g_nGroups,
        g_nVrtxPerCell,
        g_nAngles,
        g_nCells);

    Kokkos::deep_copy(c_device_source, host_source);
    Kokkos::deep_copy(c_device_psi_bound, host_psi_bound);
    Kokkos::deep_copy(c_device_psi, host_psi);


    // Actually do graph traversal
    // THIS IS WEIRD!!!
    auto nCells = g_nCells;
    auto nGroups = g_nGroups;
    auto device_omega_dot_n = c_device_omega_dot_n;
    auto device_adj_cell = c_device_adj_cell;
    auto device_neighbor_vertex = c_device_neighbor_vertex;
    auto device_adj_proc = c_device_adj_proc;
    auto device_side = c_device_side;
    auto device_sigma_t = c_device_sigma_t;
    auto device_cell_volume = c_device_cell_volume;
    auto device_face_area = c_device_face_area;
    auto device_cell_to_face_vertex = c_device_cell_to_face_vertex;
    auto device_source = c_device_source;
    auto device_psi_bound = c_device_psi_bound;
    auto device_psi = c_device_psi;

    auto lambda = KOKKOS_LAMBDA(int item) {
        int cell = item % nCells;
        int angle = item / nCells;

        Transport::updateKokkos(
            cell,
            angle,
            device_source,
            device_psi_bound,
            device_psi,
            nGroups,
            device_omega_dot_n,
            device_adj_cell,
            device_neighbor_vertex,
            device_adj_proc,
            device_side,
            device_cell_volume,
            device_sigma_t,
            device_face_area,
            device_cell_to_face_vertex);
    };
    

    // Perform sweep and copy Psi back to host
    Kokkos::parallel_for(policy, lambda, "traverse-dag");
    Kokkos::deep_copy(host_psi, device_psi);

    
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

