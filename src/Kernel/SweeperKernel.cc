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

#include "SweeperKernel.hh"
#include "UtilKernel.hh"
#include "Global.hh"
#include "PsiData.hh"
#include "Comm.hh"
//#include "CommSides.hh"
#include "SourceIteration.hh"
#include <math.h>

using namespace std;


////////////////////////////////////////////////////////////////////////////////
//            SweeperPBJSI functions
////////////////////////////////////////////////////////////////////////////////

/*
    solve
*/
void SweeperKernel::solve()
{
    PhiData phi0;
    PhiData phi1;
    PsiData totalSource;
    PsiBoundData psiBound0;
    vector<UINT> sourceIts;
    
    
    // Initialize source and psi
    c_source.setToValue(1.0);
    c_psi.setToValue(0.0);
    c_psiBound.setToValue(0.0);
    phi0.setToValue(0.0);

    
    // Source iterate till converged
    UINT iter = 1;
    while (iter < g_ddIterMax) {
        
        UtilKernel::calcTotalSource(c_source, phi0, totalSource);
        UtilKernel::sweepLocal(c_psi, totalSource, c_psiBound);
        UtilKernel::psiToPhi(phi1, c_psi);
        //c_commSides.commSides(c_psi, c_psiBound);
        

        // Check tolerance and set phi0 = phi1
        double errL1 = 0.0;
        double normL1 = 0.0;
        for (UINT i = 0; i < phi1.size(); i++) {
            errL1  += fabs(phi1[i] - phi0[i]);
            normL1 += fabs(phi1[i]);
            phi0[i] = phi1[i];
        }
        
        Comm::gsum(errL1);
        Comm::gsum(normL1);

        if (Comm::rank() == 0) {
            printf("Relative error: %e\n", errL1 / normL1);
        }

        //if (errL1 / normL1 < g_ddErrMax)
        //    break;
        
        
        // Increment iter
        iter++;
    }
    
    
    // Print statistics
    if (Comm::rank() == 0) {
        printf("PBJ Iters: %" PRIu64 "\n", iter);
        printf("Num source iterations:");
        for (UINT i = 0; i < sourceIts.size(); i++)
            printf(" %" PRIu64, sourceIts[i]);
        printf("\n");
    }
}





