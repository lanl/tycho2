#ifndef __PRIORITIES_HH__
#define __PRIORITIES_HH__

#include "Mat.hh"
#include "Typedef.hh"
#include <vector>

namespace Priorities
{

void calcPriorities(const UINT maxComputePerStep,
                    const UINT intraAngleP, const UINT interAngleP, 
                    Mat2<UINT> &priorities);

}

#endif