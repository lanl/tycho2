/*
    Sweeper.hh
*/

#ifndef __SWEEPER_HH__
#define __SWEEPER_HH__

#include "PsiData.hh"

namespace Sweeper 
{
    void sweep(PsiData &psi, const PsiData &source, const double sigmaTotal);
};


#endif
