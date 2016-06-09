#ifndef __SWEEPER2_HH__
#define __SWEEPER2_HH__

#include "Mat.hh"
#include "PsiData.hh"
#include "Typedef.hh"
#include <vector>


class Sweeper2
{
public:
    Sweeper2(const UINT maxComputePerStep,
             const UINT intraAngleP, 
             const UINT interAngleP, 
             const double sigmaTotal);
    
    void sweep(PsiData &psi, const PsiData &source);

private:
    Mat2<UINT> c_priorities;
    UINT c_maxComputePerStep;
    double c_sigmaTotal;
};

#endif