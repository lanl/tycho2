#ifndef __USER_MULT_HH__
#define __USER_MULT_HH__

#include <petscmat.h>
#include <petscvec.h>

class UserMult
{
public:
    UserMult();

    static void mult(Mat mat, Vec x, Vec y);
    
};

#endif
