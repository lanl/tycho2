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


#ifndef __OPERATOR_HH__
#define __OPERATOR_HH__

#include <petscmat.h>
#include <petscvec.h>
#include "Mat.hh"
#include "PsiData.hh"
#include "Typedef.hh"
#include <vector>
#include "Global.hh"
#include "SweeperSchurBoundary.hh"
#include "SweepDataSchur.hh"


class Operator {

public:

    Operator();//MPI_Comm mpiComm, const std::vector<UINT> &adjRanks, const std::vector<std::vector<MetaData>> &sendMetaData, const std::vector<UINT> &numSendPackets, const std::vector<UINT> &numRecvPackets, SweepDataSchur &sweepData, PsiData &psi); 
    
    //PetscErrorCode Schur(Mat mat, Vec x, Vec y);

    //void commSides(const std::vector<UINT> &adjRanks, const std::vector<std::vector<MetaData>> &sendMetaData,const std::vector<UINT> &numSendPackets,const std::vector<UINT> &numRecvPackets,SweepDataSchur &sweepData);

    MPI_Comm getmpiComm();
    const std::vector<UINT> getadjRanks();
    const std::vector<std::vector<MetaData>> getsendMetaData();
    const std::vector<UINT> getnumSendPackets();
    const std::vector<UINT> getnumRecvPackets();  
    PsiData getpsi();
    PsiData getpsiBound();
    PsiData getpsiSource();
    double getsigmaTotal();
        
     
    void setmpiComm(MPI_Comm MPI_IN);
    void setadjRanks(std::vector<UINT> RANKS_IN);
    void setsendMetaData(std::vector<std::vector<MetaData>> META_IN);
    void setnumSendPackets(std::vector<UINT> SEND_IN);
    void setnumRecvPackets(std::vector<UINT> RECV_IN); 
    void setpsi(PsiData PSI_IN);
    void setpsiBound(PsiData BOUND_IN);
    void setpsiSource(PsiData SOURCE_IN);
    void setsigmaTotal(double SIGMA_IN);
    
    

    static MPI_Comm c_mpiComm;
    static std::vector<UINT> c_adjRanks;
    static std::vector<std::vector<MetaData>> c_sendMetaData;
    static std::vector<UINT> c_numSendPackets;
    static std::vector<UINT> c_numRecvPackets;  
    static PsiData c_psi;
    static PsiData c_psiBound;
    static PsiData c_psiSource;
    static double c_sigmaTotal;
    
};

#endif



