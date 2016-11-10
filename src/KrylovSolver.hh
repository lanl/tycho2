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

#ifndef __KRYLOV_SOLVER_HH__
#define __KRYLOV_SOLVER_HH__


/*
    If using PETSc.
*/
#if USE_PETSC
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







/*
    If not using PETSc.
*/
#else
#include "Global.hh"

class KrylovSolver
{
public:
    typedef void (*Function)(const double*, double*, void*);


    KrylovSolver(UINT localVecSize, double rtol, UINT iterMax,
                 Function lhsOperator)
    {
        UNUSED_VARIABLE(localVecSize);
        UNUSED_VARIABLE(rtol);
        UNUSED_VARIABLE(iterMax);
        UNUSED_VARIABLE(lhsOperator);
        Insist(false, "KrylovSolver: PETSc not in use.");
    }
    ~KrylovSolver() {}


    void solve()                    { }
    double* getB()                  { return NULL; }
    void releaseB()                 {}
    double* getX()                  { return NULL; }
    void releaseX()                 {}
    UINT getNumIterations()         { return 0; }
    double getResidualNorm()        { return 0.0; }
    void setData(void *data)        { UNUSED_VARIABLE(data); }
    void setInitialGuessNonzero()   {}


    struct Data
    {
        void *data;
        Function lhsOperator;
    };

private:
    
};



#endif




#endif

