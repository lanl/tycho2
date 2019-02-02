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

#if USE_PETSC

#include "KrylovSolver.hh"
#include "Comm.hh"


static
PetscErrorCode lhsPetsc(Mat mat, Vec x, Vec b)
{
    void *voidData;
    void *userData;
    KrylovSolver::Data *krylovData;
    KrylovSolver::Function lhsOperator;
    const double *xArray;
    double *bArray;
    
    
    // Get data for the solve
    MatShellGetContext(mat, &voidData);
    krylovData = (KrylovSolver::Data*) voidData;
    userData = krylovData->data;
    lhsOperator = krylovData->lhsOperator;
    

    // Get vector arrays
    VecGetArrayRead(x, &xArray);
    VecGetArray(b, &bArray);
    

    // Call user function
    lhsOperator(xArray, bArray, userData);


    // Restore vectors
    VecRestoreArrayRead(x, NULL);
    VecRestoreArray(b, NULL);

    return 0;
}


KrylovSolver::KrylovSolver(UINT localVecSize, double rtol, UINT iterMax, 
                           Function lhsOperator)
{
    // Create the vectors
    UINT globalVecSize = localVecSize;
    Comm::gsum(globalVecSize);
    VecCreate(MPI_COMM_WORLD, &c_x);
    VecSetSizes(c_x, localVecSize, globalVecSize);
    VecSetType(c_x, VECMPI);
    VecDuplicate(c_x, &c_b);
    

    // Create matrix shell and define it as the operator
    MatCreateShell(MPI_COMM_WORLD, localVecSize, localVecSize, 
                   globalVecSize, globalVecSize, (void*)(NULL), &c_mat);
    MatShellSetOperation(c_mat, MATOP_MULT, (void(*)(void))lhsPetsc);
    

    // Create ksp and set tolerances
    KSPCreate(MPI_COMM_WORLD, &c_ksp);
    KSPSetOperators(c_ksp, c_mat, c_mat);
    KSPSetTolerances(c_ksp, rtol, PETSC_DEFAULT, PETSC_DEFAULT, iterMax);
    KSPSetType(c_ksp, KSPGMRES);


    // Set data
    c_krylovData.lhsOperator = lhsOperator;
}


void KrylovSolver::setData(void *data)
{
    c_krylovData.data = data;
    MatShellSetContext(c_mat, &c_krylovData);
}


void KrylovSolver::setInitialGuessNonzero()
{
    KSPSetInitialGuessNonzero(c_ksp, PETSC_TRUE);
}


double* KrylovSolver::getB()
{
    double *bArray;
    VecGetArray(c_b, &bArray);
    return bArray;
}


void KrylovSolver::releaseB()
{
    VecRestoreArray(c_b, NULL);
}


double* KrylovSolver::getX()
{
    double *xArray;
    VecGetArray(c_x, &xArray);
    return xArray;
}


void KrylovSolver::releaseX()
{
    VecRestoreArray(c_x, NULL);
}


UINT KrylovSolver::getNumIterations()
{
    int iters;
    KSPGetIterationNumber(c_ksp, &iters);
    return iters;
}


double KrylovSolver::getResidualNorm()
{
    double norm;
    KSPGetResidualNorm(c_ksp, &norm);
    return norm;
}


void KrylovSolver::solve()
{
    KSPSolve(c_ksp, c_b, c_x);
}


KrylovSolver::~KrylovSolver()
{
    // Destroy vectors, ksp, and matrices to free work space
    VecDestroy(&c_x);
    VecDestroy(&c_b);
    MatDestroy(&c_mat);
    KSPDestroy(&c_ksp);
}

#endif
