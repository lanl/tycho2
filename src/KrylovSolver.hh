#ifndef __KRYLOV_SOLVER_HH__
#define __KRYLOV_SOLVER_HH__


#include "Global.hh"
#include <petscmat.h>
#include <petscvec.h>
#include <petscksp.h>

class KrylovSolver
{
public:
    typedef void (*Function)(const double*, double*, void*);


    KrylovSolver(UINT localVecSize, double rtol, UINT iterMax,
                 Function lhsOperator);
    ~KrylovSolver();


    void solve();
    double* getB();
    void releaseB();
    double* getX();
    void releaseX();
    UINT getNumIterations();
    double getResidualNorm();
    void setData(void *data);
    void setInitialGuessNonzero();


    struct Data
    {
        void *data;
        Function lhsOperator;
    };

private:
    KSP c_ksp;
    Mat c_mat;
    Vec c_x, c_b;
    Data c_krylovData;
};

#endif

