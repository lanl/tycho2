#ifndef __SWEEPER_PBJ_HH__
#define __SWEEPER_PBJ_HH__

#include "Mat.hh"
#include "PsiData.hh"
#include "Typedef.hh"
#include <vector>


class SweeperPBJ
{
public:
    SweeperPBJ(const double sigmaTotal);
    
    void sweep(PsiData &psi, const PsiData &source);

private:
    double c_sigmaTotal;
};

#endif