/*
    Assert.cc
    
    Implements assertions for code that can be turned on and off.
    These should only be accessed by the macros defined in Assert.hh
    
    Authors:
        Kris Garrett - ckgarrett@lanl.gov
*/

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
