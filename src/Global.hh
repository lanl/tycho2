/*
    Global.hh
*/

#ifndef __GLOBAL_HH__
#define __GLOBAL_HH__

#include "TychoMesh.hh"
#include "SweepSchedule.hh"
#include "Quadrature.hh"
#include "Typedef.hh"

#ifdef NO_EXTERN
#define EXTERN
#else
#define EXTERN extern
#endif

static const UINT g_ndim = 3;
static const UINT g_nVrtxPerCell = 4;
static const UINT g_nVrtxPerFace = 3;
static const UINT g_nFacePerCell = 4;

enum SweepType
{
    SweepType_OriginalTycho1,
    SweepType_OriginalTycho2,
    SweepType_TraverseGraph,
};

EXTERN UINT g_nAngleGroups;
EXTERN UINT g_nThreads;
EXTERN UINT g_nGroups;
EXTERN UINT g_snOrder;
EXTERN UINT g_iterMax;
EXTERN double g_errMax;
EXTERN UINT g_maxCellsPerStep;
EXTERN UINT g_intraAngleP;
EXTERN UINT g_interAngleP;
EXTERN SweepType g_sweepType;
EXTERN TychoMesh *g_spTychoMesh;
EXTERN SweepSchedule **g_spSweepSchedule;
EXTERN Quadrature *g_quadrature;

#endif

