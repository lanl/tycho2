/*
    Assert.cc
    
    Implements assertions for code that can be turned on and off.
    These should only be accessed by the macros defined in Assert.hh
*/

#include "Assert.hh"
#include <stdio.h>
#include <stdlib.h>
//#include <mpi.h>


namespace Assert
{

/*
    my_assert
*/
void my_assert(const char *cond, const char *file, const int line)
{
    printf("Assertion: %s, failed in %s, line %d.\n\n\n", cond, file, line);
    abort();
}


/*
    my_insist
*/
void my_insist(const char *cond, const char *msg, const char *file, const int line)
{
    printf("Insist: %s, failed in %s, line %d.\n", cond, file, line);
    printf("The following message was provided: %s\n\n\n", msg);
    abort();
}

} // end of Assert namespace

