#ifndef __SOLVER_HH__
#define __SOLVER_HH__

#include "Typedef.hh"

namespace Solver 
{
    void solve(const double sigmaTotal, const double sigmaScat, 
               const UINT iterMax, const double errMax);
};


#endif