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


#include "SweeperSchur.hh"
#include "Util.hh"
#include "Problem.hh"
#include "SourceIteration.hh"
#include "Global.hh"
#include "PsiData.hh"
#include "Comm.hh"
#include "CommSides.hh"
#include <vector>
#include <math.h>
#include <string.h>


namespace
{

/*
    Data needed for Schur, SchurOuter, and SchurKrylov functions
*/
struct SchurData
{
    UINT psiBoundSize;
    CommSides *commSides;
    PsiData *psi;
    PsiBoundData *psiBound;
    PsiData *source;
    
    // Only needed for SchurKrylov
    PhiData *phi;

    // Only needed for SchurOuter
    std::vector<UINT> *sourceIts;
    SweeperSchurOuter *sweeperSchurOuter;
};


/*
    psiBoundToVec 
*/
void psiBoundToVec(double *x, const PsiBoundData &psiBound)
{
    int xArrayIndex = 0;

    // Set array
    for (UINT angle = 0; angle < g_nAngles; ++angle) {
    for (UINT cell = 0; cell < g_nCells; ++cell) {
    for (UINT group = 0; group < g_nGroups; group++) {
    for (UINT face = 0; face < g_nFacePerCell; face++) {
        if (g_tychoMesh->isIncoming(angle, cell, face)) {
            UINT adjCell = g_tychoMesh->getAdjCell(cell, face);
            UINT adjRank = g_tychoMesh->getAdjRank(cell, face);
            
            // On internal boundary
            if (adjCell == TychoMesh::BOUNDARY_FACE && 
                adjRank != TychoMesh::BAD_RANK)
            {
                UINT side = g_tychoMesh->getSide(cell, face);
                for (UINT fvrtx = 0; fvrtx < g_nVrtxPerFace; fvrtx++) {
                    x[xArrayIndex] = psiBound(group, fvrtx, angle, side);
                    xArrayIndex++;
                }
            }
        }
    }}}}
}


/*
    vecToPsiBound
*/
void vecToPsiBound(const double *x, PsiBoundData &psiBound)
{
    int xArrayIndex = 0;

    // Set psiBound
    for (UINT angle = 0; angle < g_nAngles; ++angle) {
    for (UINT cell = 0; cell < g_nCells; ++cell) {
    for (UINT group = 0; group < g_nGroups; group++) {
    for (UINT face = 0; face < g_nFacePerCell; face++) {
        if (g_tychoMesh->isIncoming(angle, cell, face)) {
            UINT adjCell = g_tychoMesh->getAdjCell(cell, face);
            UINT adjRank = g_tychoMesh->getAdjRank(cell, face);
            
            // On internal boundary
            if (adjCell == TychoMesh::BOUNDARY_FACE && 
                adjRank != TychoMesh::BAD_RANK)
            {
                UINT side = g_tychoMesh->getSide(cell, face);
                for (UINT fvrtx = 0; fvrtx < g_nVrtxPerFace; fvrtx++) {
                    psiBound(group, fvrtx, angle, side) = x[xArrayIndex];
                    xArrayIndex++;
                }
            }
        }
    }}}}
}


/*
    getPsiBoundSize
*/
UINT getPsiBoundSize()
{
    // Start an array index
    UINT size = 0;

    // Incoming flux
    for (UINT angle = 0; angle < g_nAngles; ++angle) {
    for (UINT cell = 0; cell < g_nCells; ++cell) {
    for (UINT group = 0; group < g_nGroups; group++) {
    for (UINT face = 0; face < g_nFacePerCell; face++) {
        
        if (g_tychoMesh->isIncoming(angle, cell, face)) {
            
            UINT neighborCell = g_tychoMesh->getAdjCell(cell, face);
            UINT adjRank = g_tychoMesh->getAdjRank(cell, face);
            
            // In local mesh
            if (neighborCell == TychoMesh::BOUNDARY_FACE &&  
                adjRank != TychoMesh::BAD_RANK)
            {
                for (UINT fvrtx = 0; fvrtx < g_nVrtxPerFace; fvrtx++) {
                    size++;
                }
            }
        }
    }}}}

    return size;
}


} // End anonymous namespace




////////////////////////////////////////////////////////////////////////////////
//              SweeperSchur functions
////////////////////////////////////////////////////////////////////////////////


/*
    Schur 

    This performs the sweep and returns the boundary
    Performs b = (I - W L_I^{-1} L_B) x
*/
static
void Schur(const double *x, double *b, void *voidData)
{
    // Unpack data
    SchurData *data = (SchurData*) voidData;


    // x -> psiBound data type
    vecToPsiBound(x, *data->psiBound);


    // Perform W L_I^{-1} L_B
    Util::sweepLocal(*data->psi, *data->source, *data->psiBound);
    data->commSides->commSides(*data->psi, *data->psiBound);

    
    // psiBound -> b
    psiBoundToVec(b, *data->psiBound);


    // b = x - b
    for (UINT i = 0; i < data->psiBoundSize; i++) {
        b[i] = x[i] - b[i];
    }
}


/*
    solve
*/
void SweeperSchur::solve()
{
    // Setup problem
    KrylovSolver ks(getPsiBoundSize(), g_ddErrMax, g_ddIterMax, Schur);
    c_krylovSolver = &ks;
    c_iters = 0;
    Problem::getSource(c_source);
    c_psi.setToValue(0.0);


    // Solve
    if (g_useSourceIteration)
        SourceIteration::fixedPoint(*this, c_psi, c_source);
    else
        SourceIteration::krylov(*this, c_psi, c_source);

    
    // Print data
    if (Comm::rank() == 0) {
        printf("Schur: Num local sweeps: %" PRIu64 "\n", c_iters);
    }
}


/*
    sweep
    
    Run the Krylov solver
*/
void SweeperSchur::sweep(PsiData &psi, const PsiData &source, bool zeroPsiBound)
{
    UNUSED_VARIABLE(zeroPsiBound);
    
    // Initialize variables
    PsiData zeroSource;
    zeroSource.setToValue(0.0);
    PsiBoundData psiBound;
    
    double rnorm;
    UINT its;
    double *x;
    double *b;

    
    // Set static variables
    SchurData data;
    data.commSides = &c_commSides;
    data.psi = &psi;
    data.psiBound = &psiBound;
    data.source = &zeroSource;
    data.psiBoundSize = getPsiBoundSize();
    c_krylovSolver->setData(&data);
    
    
    // Initial guess
    x = c_krylovSolver->getX();
    psiBoundToVec(x, c_psiBoundPrev);
    c_krylovSolver->releaseX();
    c_krylovSolver->setInitialGuessNonzero();
    
    
    // Set RHS
    if (Comm::rank() == 0) {
        printf("      Schur: Set RHS\n");
    }
    psiBound.setToValue(0.0);
    Util::sweepLocal(psi, source, psiBound);

    c_commSides.commSides(psi, psiBound);
    b = c_krylovSolver->getB();
    psiBoundToVec(b, psiBound);
    c_krylovSolver->releaseB();


    // Solve the system
    if (Comm::rank() == 0) {
        printf("      Schur: Krylov solve\n");
    }
    c_krylovSolver->solve();
    

    // Got Psi_B.  Now do a local sweep for Psi.
    if (Comm::rank() == 0) {
        printf("      Schur: Post solve processing\n");
    }
    x = c_krylovSolver->getX();
    vecToPsiBound(x, psiBound);
    vecToPsiBound(x, c_psiBoundPrev);
    c_krylovSolver->releaseX();
    Util::sweepLocal(psi, source, psiBound);
    
    
    // Print some stats
    its = c_krylovSolver->getNumIterations();
    rnorm = c_krylovSolver->getResidualNorm();
    if (Comm::rank() == 0) {
        printf("      Schur: Krylov iterations: %" PRIu64 "\n", its);
        printf("      Schur: Residual norm: %e\n", rnorm);
    }

    
    // Increment number of sweep calls
    c_iters += 2 + its;
}


////////////////////////////////////////////////////////////////////////////////
//              SweeperSchurOuter functions
////////////////////////////////////////////////////////////////////////////////


/*
    SchurOuter 

    Performs b = (I - W (L_i - MSD)^{-1} L_B) x
*/
static
void SchurOuter(const double *x, double *b, void *voidData)
{
    // Get data for the solve
    SchurData *data = (SchurData*) voidData;


    // x -> psiBound
    vecToPsiBound(x, *data->psiBound);


    // Perform W (L_i - MSD)^{-1} L_B
    UINT its;
    data->source->setToValue(0.0);
    if (g_useSourceIteration)
        its = SourceIteration::fixedPoint(*data->sweeperSchurOuter, *data->psi,
                                          *data->source);
    else
        its = SourceIteration::krylov(*data->sweeperSchurOuter, *data->psi,
                                      *data->source);
    data->sourceIts->push_back(its);
    data->commSides->commSides(*data->psi, *data->psiBound);

    
    // psiBound -> b
    psiBoundToVec(b, *data->psiBound);


    // b = x - b
    for (UINT i = 0; i < data->psiBoundSize; i++) {
        b[i] = x[i] - b[i];
    }
}


/*
    solve
*/
void SweeperSchurOuter::solve()
{
    // Variables
    double rnorm;
    UINT its;
    UINT sourceIts1, sourceIts3;
    std::vector<UINT> sourceItsVec;
    double *x;
    double *b;


    // Initialize class variables
    Problem::getSource(c_source);
    c_psi.setToValue(0.0);
    c_psiBound.setToValue(0.0);
    KrylovSolver ks(getPsiBoundSize(), g_ddErrMax, g_ddIterMax, SchurOuter);
    c_krylovSolver = &ks;
    
    
    // Set data for Krylov solver
    SchurData data;
    data.commSides = &c_commSides;
    data.psi = &c_psi;
    data.psiBound = &c_psiBound;
    data.source = &c_source;
    data.sourceIts = &sourceItsVec;
    data.sweeperSchurOuter = this;
    data.psiBoundSize = getPsiBoundSize();
    c_krylovSolver->setData(&data);


    // RHS
    if (Comm::rank() == 0) {
        printf("SchurOuter: Calculate RHS\n");
    }
    
    if (g_useSourceIteration)
        sourceIts1 = SourceIteration::fixedPoint(*this, c_psi, c_source);
    else
        sourceIts1 = SourceIteration::krylov(*this, c_psi, c_source);
    
    c_commSides.commSides(c_psi, c_psiBound);
    b = c_krylovSolver->getB();
    psiBoundToVec(b, c_psiBound);
    c_krylovSolver->releaseB();


    // Initial values
    c_psiBound.setToValue(0.0);
    c_psi.setToValue(0.0);
    c_source.setToValue(0.0);
    
    
    // Solve the system
    if (Comm::rank() == 0) {
        printf("SchurOuter: Krylov solve\n");
    }
    c_krylovSolver->solve();
    

    // Got Psi_B.  Now calculate Psi.
    x = c_krylovSolver->getX();
    vecToPsiBound(x, c_psiBound);
    c_krylovSolver->releaseX();
    
    Problem::getSource(c_source);
    if (g_useSourceIteration)
        sourceIts3 = SourceIteration::fixedPoint(*this, c_psi, c_source);
    else
        sourceIts3 = SourceIteration::krylov(*this, c_psi, c_source);

    
    // Print some stats
    its = c_krylovSolver->getNumIterations();
    rnorm = c_krylovSolver->getResidualNorm();
    
    if (Comm::rank() == 0) {
        printf("SchurOuter: Krylov iterations: %" PRIu64 "\n", its);
        printf("SchurOuter: Residual Norm: %e\n", rnorm);
        printf("SchurOuter: Num sweeps Q: %" PRIu64 "\n", sourceIts1);
        printf("SchurOuter: Num sweeps KSP:");
        for (UINT i = 0; i < sourceItsVec.size(); i++)
            printf(" %" PRIu64, sourceItsVec[i]);
        printf("\n");
        printf("SchurOuter: Num sweeps END: %" PRIu64 "\n", sourceIts3);
    }
}


/*
    sweep
    
    Run the Krylov solver
*/
void SweeperSchurOuter::sweep(PsiData &psi, const PsiData &source, 
                              bool zeroPsiBound)
{
    if (zeroPsiBound) {
        Util::sweepLocal(psi, source, c_zeroPsiBound);
    }
    else {
        Util::sweepLocal(psi, source, c_psiBound);
    }

}


////////////////////////////////////////////////////////////////////////////////
//              SweeperSchurKrylov functions
////////////////////////////////////////////////////////////////////////////////


/*
    SchurKrylov

    Performs (I - W L_I^{-1} L_B) Psi_B -   (W L_I^{-1} M S)   Phi
              -(D L_I^{-1} L_B)   Psi_B + (I - D L_I^{-1} M S) Phi
*/
static
void SchurKrylov(const double *x, double *b, void *voidData)
{
    // Unpack data for the solve
    SchurData *data = (SchurData*) voidData;
    UINT psiBoundSize = data->psiBoundSize;
    CommSides &commSides = *(data->commSides);
    PsiData &psi = *(data->psi);
    PsiBoundData &psiBound = *(data->psiBound);
    PsiData &source = *(data->source);
    PhiData &phi = *(data->phi);


    // x -> (Psi_B, Phi)
    vecToPsiBound(x, psiBound);
    const PhiData phi1(const_cast<double*>(&x[psiBoundSize]));
    Util::operatorS(phi1, phi);
    Util::phiToPsi(phi, source);


    // Perform most of the operator
    Util::sweepLocal(psi, source, psiBound);
    commSides.commSides(psi, psiBound);
    Util::psiToPhi(phi, psi);

    
    // (Psi_B, Phi) -> b
    psiBoundToVec(b, psiBound);
    for (UINT i = 0; i < phi.size(); i++) {
        b[i+psiBoundSize] = phi[i];
    }


    // b = x - b
    for (UINT i = 0; i < psiBoundSize + phi.size(); i++) {
        b[i] = x[i] - b[i];
    }
}


/*
    solve
*/
void SweeperSchurKrylov::solve()
{
    
    double *x;
    double *b;
    double rnorm;
    UINT its;
    PhiData phi;
    SchurData data;
    
    UINT psiBoundSize = getPsiBoundSize();
    UINT vecSize = psiBoundSize + phi.size();
    const bool zeroPsiBound = false;

    KrylovSolver ks(vecSize, g_ddErrMax, g_ddIterMax, SchurKrylov);
    c_krylovSolver = &ks;

    
    // Set data for Krylov solver
    data.psiBoundSize = psiBoundSize;
    data.commSides = &c_commSides;
    data.psi = &c_psi;
    data.psiBound = &c_psiBound;
    data.source = &c_source;
    data.phi = &phi;
    c_krylovSolver->setData(&data);


    // Calculate RHS
    Problem::getSource(c_source);
    c_psi.setToValue(0.0);
    c_psiBound.setToValue(0.0);
    
    if (Comm::rank() == 0) {
        printf("SchurKrylov: Calculate RHS\n");
    }
    sweep(c_psi, c_source, zeroPsiBound);
    
    c_commSides.commSides(c_psi, c_psiBound);
    Util::psiToPhi(phi, c_psi);
    
    b = c_krylovSolver->getB();
    psiBoundToVec(b, c_psiBound);
    for (UINT i = 0; i < phi.size(); i++) {
        b[i+psiBoundSize] = phi[i];
    }
    c_krylovSolver->releaseB();


    // Initial data
    c_psiBound.setToValue(0.0);
    c_psi.setToValue(0.0);
    c_source.setToValue(0.0);
    
    
    // Solve the system
    if (Comm::rank() == 0) {
        printf("SchurKrylov: Solve system\n");
    }
    c_krylovSolver->solve();
    

    // Got Psi_B and Phi.  Solve for Psi.
    if (Comm::rank() == 0) {
        printf("SchurKrylov: Postprocess\n");
    }
    x = c_krylovSolver->getX();
    vecToPsiBound(x, c_psiBound);
    for (UINT i = 0; i < phi.size(); i++) {
        phi[i] = x[i+psiBoundSize];
    }
    Problem::getSource(c_source);
    Util::calcTotalSource(c_source, phi, c_source);
    c_krylovSolver->releaseX();
    
    sweep(c_psi, c_source, zeroPsiBound);
    
    
    // Print some stats
    its = c_krylovSolver->getNumIterations();
    rnorm = c_krylovSolver->getResidualNorm();

    if (Comm::rank() == 0) {
        printf("SchurKrylov: Krylov iterations: %" PRIu64 "\n", its);
        printf("SchurKrylov: Residual Norm: %e\n", rnorm);
    }
}


/*
    sweep
*/
void SweeperSchurKrylov::sweep(PsiData &psi, const PsiData &source, 
                               bool zeroPsiBound)
{
    UNUSED_VARIABLE(zeroPsiBound);
    Util::sweepLocal(psi, source, c_psiBound);
}

