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


