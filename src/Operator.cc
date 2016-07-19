/* Operators: right now only one is Schur, which takes in x, a guess for psi on the incoming boundary, and returns b, the corresponding psi values for the outgoing boundary */

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

#include "Operator.hh"
#include "Global.hh"
#include <vector>
#include "TraverseGraph.hh"
#include "Typedef.hh"
#include "Comm.hh"
#include <petscmat.h>
#include <petscvec.h>
#include <petscksp.h>
#include <omp.h>
#include "SweepDataSchur.hh"

//using namespace std;




/*
    Constructor
*/
Operator::Operator()//MPI_Comm mpiComm, const std::vector<UINT> &adjRanks, const std::vector<std::vector<MetaData>> &sendMetaData, const std::vector<UINT> &numSendPackets, const std::vector<UINT> &numRecvPackets, SweepDataSchur &sweepData, PsiData &psi)
//: c_mpiComm(mpiComm), c_adjRanks(adjRanks), c_sendMetaData(sendMetaData), c_numSendPackets(numSendPackets), c_numRecvPackets(numRecvPackets) ,c_sweepData(sweepData), c_psi(psi)
{

}


/*
   getData 
*/
MPI_Comm Operator::getmpiComm()
{return Operator::c_mpiComm;}

const std::vector<UINT> Operator::getadjRanks()
{return Operator::c_adjRanks;}

const std::vector<std::vector<MetaData>> Operator::getsendMetaData()
{return Operator::c_sendMetaData;}

const std::vector<UINT> Operator::getnumSendPackets()
{return Operator::c_numSendPackets;}

const std::vector<UINT> Operator::getnumRecvPackets()
{return Operator::c_numRecvPackets;}
    
PsiData Operator::getpsi()
{return Operator::c_psi;}

PsiData Operator::getpsiBound()
{return Operator::c_psiBound;}

PsiData Operator::getpsiSource()
{return Operator::c_psiSource;}

double Operator::getsigmaTotal()
{return Operator::c_sigmaTotal;}
/*
   setData 
*/
void Operator::setmpiComm(MPI_Comm MPI_IN)
{Operator::c_mpiComm = MPI_IN;}

void Operator::setadjRanks(const std::vector<UINT> RANKS_IN)
{Operator::c_adjRanks = RANKS_IN;}

void Operator::setsendMetaData(const std::vector<std::vector<MetaData>> META_IN)
{Operator::c_sendMetaData = META_IN;}

void Operator::setnumSendPackets(const std::vector<UINT> SEND_IN)
{Operator::c_numSendPackets = SEND_IN;}

void Operator::setnumRecvPackets(const std::vector<UINT> RECV_IN)
{Operator::c_numRecvPackets = RECV_IN;}
   
void Operator::setpsi(PsiData PSI_IN)
{Operator::c_psi = PSI_IN;}

void Operator::setpsiBound(PsiData BOUND_IN)
{Operator::c_psiBound = BOUND_IN;}

void Operator::setpsiSource(PsiData SOURCE_IN)
{Operator::c_psiSource = SOURCE_IN;}

void Operator::setsigmaTotal(double SIGMA_IN)
{Operator::c_sigmaTotal = SIGMA_IN;}


