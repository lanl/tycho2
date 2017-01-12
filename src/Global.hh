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

#ifndef __GLOBAL_HH__
#define __GLOBAL_HH__

#include <string>
#include <cinttypes>
#include <vector>


// Hack so I don't have to redefine extern variables in Global.cc
#ifdef NO_EXTERN
#define EXTERN
#else
#define EXTERN extern
#endif


// Forward declaration of classes needed for global pointers below
class Quadrature;
class TychoMesh;
class SweepSchedule;
class GraphTraverser;


// Macro to get around some warnings
#define UNUSED_VARIABLE(x) (void)(x)


// Shorter version of uint64_t
// Also allows changing the UINT type
typedef uint64_t UINT;


// Global constants
static const UINT g_ndim = 3;
static const UINT g_nVrtxPerCell = 4;
static const UINT g_nVrtxPerFace = 3;
static const UINT g_nFacePerCell = 4;


// Enum types
enum SweepType
{
    SweepType_OriginalTycho1,
    SweepType_OriginalTycho2,
    SweepType_TraverseGraph,
    SweepType_PBJ,
    SweepType_PBJOuter,
    SweepType_Schur,
    SweepType_SchurOuter,
    SweepType_PBJSI,
    SweepType_SchurKrylov
};

enum GaussElim
{
    GaussElim_Original,
    GaussElim_NoPivot,
    GaussElim_CramerGlu,
    GaussElim_CramerIntel
};


// Global variables
EXTERN UINT g_nAngleGroups;
EXTERN UINT g_nThreads;
EXTERN UINT g_nGroups;
EXTERN UINT g_snOrder;
EXTERN UINT g_iterMax;
EXTERN double g_errMax;
EXTERN std::vector<double> g_sigmaT;
EXTERN std::vector<double> g_sigmaS;
EXTERN UINT g_maxCellsPerStep;
EXTERN UINT g_intraAngleP;
EXTERN UINT g_interAngleP;
EXTERN SweepType g_sweepType;
EXTERN TychoMesh *g_tychoMesh;
EXTERN SweepSchedule **g_sweepSchedule;
EXTERN Quadrature *g_quadrature;
EXTERN GraphTraverser *g_graphTraverserForward;
EXTERN GaussElim g_gaussElim;
EXTERN bool g_outputFile;
EXTERN std::string g_outputFilename;
EXTERN UINT g_nAngles;
EXTERN UINT g_nCells;
EXTERN double g_ddErrMax;
EXTERN UINT g_ddIterMax;
EXTERN bool g_useSourceIteration;
EXTERN bool g_useOneSidedMPI;

#endif

