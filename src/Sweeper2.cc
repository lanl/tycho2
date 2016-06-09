/*
    Sweeper2.cc
*/

#include "Sweeper2.hh"
#include "Global.hh"
#include "TraverseGraph.hh"
#include "Priorities.hh"
#include "Transport.hh"
#include "PsiData.hh"
#include <omp.h>
#include <vector>

using namespace std;


/*
    populateLocalPsiBound
    
    Put data from neighboring cells into localPsiBound(fvrtx, face, group).
*/
static
void populateLocalPsiBound(const UINT angle, const UINT cell, 
                           const PsiData &psi, const PsiData &psiBound,
                           Mat3<double> &localPsiBound)
{
    // Default to 0.0
    localPsiBound = 0.0;
    
    // Populate if incoming flux
    for (UINT group = 0; group < g_nGroups; group++) {
    for (UINT face = 0; face < g_nFacePerCell; face++) {
        if (g_spTychoMesh->isIncoming(angle, cell, face)) {
            UINT neighborCell = g_spTychoMesh->getAdjCell(cell, face);
            
            // In local mesh
            if (neighborCell != TychoMesh::BOUNDARY_FACE) {
                for (UINT fvrtx = 0; fvrtx < g_nVrtxPerFace; fvrtx++) {
                    UINT neighborVrtx = 
                        g_spTychoMesh->getNeighborVrtx(cell, face, fvrtx);
                    localPsiBound(fvrtx, face, group) = 
                        psi(neighborVrtx, angle, neighborCell, group);
                }
            }
            
            // Not in local mesh
            else if (g_spTychoMesh->getAdjRank(cell, face) != TychoMesh::BAD_RANK) {
                for (UINT fvrtx = 0; fvrtx < g_nVrtxPerFace; fvrtx++) {
                    UINT side = g_spTychoMesh->getSide(cell, face);
                    localPsiBound(fvrtx, face, group) = 
                        psiBound(side, angle, fvrtx, group);
                }
            }
        }
    }}
}



/*
    SweepData
    
    Holds psi and other data for the sweep.
    Note: OpenMP assumes threading across angle.
          Without this assumption, there could be a race condition.
*/
class SweepData : public TraverseData
{
public:
    
    /*
        Constructor
    */
    SweepData(PsiData &psi, const PsiData &source,
              const double sigmaTotal,
              const Mat2<UINT> &priorities)
    : c_psi(psi), c_psiBound(g_spTychoMesh->getNSides(), 
                             g_quadrature->getNumAngles(), 
                             g_nVrtxPerFace, g_nGroups), 
      c_source(source), c_sigmaTotal(sigmaTotal), c_priorities(priorities),
      c_localFaceData(g_nThreads)
    {
        for (UINT angleGroup = 0; angleGroup < g_nThreads; angleGroup++) {
            c_localFaceData[angleGroup] = Mat2<double>(g_nVrtxPerFace, g_nGroups);
        }
    }
    
    
    /*
        data
        
        Return psi for vertices and groups at the given (cell,face,angle) tuple
    */
    virtual const char* getData(UINT cell, UINT face, UINT angle)
    {
        Mat2<double> &localFaceData = c_localFaceData[omp_get_thread_num()];
        
        for (UINT group = 0; group < g_nGroups; group++) {
        for (UINT fvrtx = 0; fvrtx < g_nVrtxPerFace; fvrtx++) {
            UINT vrtx = g_spTychoMesh->getFaceToCellVrtx(cell, face, fvrtx);
            localFaceData(fvrtx, group) = c_psi(vrtx, angle, cell, group);
        }}
        
        return (char*) (&localFaceData[0]);
    }
    
    
    /*
        getDataSize
    */
    virtual size_t getDataSize()
    {
        return g_nGroups * g_nVrtxPerFace * sizeof(double);
    }
    
    
    /*
        sideData
        
        Set psiBound for the (side, angle) pair.
    */
    virtual void setSideData(UINT side, UINT angle, const char *data)
    {
        Mat2<double> localFaceData(g_nVrtxPerFace, g_nGroups, (double*)data);
        
        for (UINT group = 0; group < g_nGroups; group++) {
        for (UINT fvrtx = 0; fvrtx < g_nVrtxPerFace; fvrtx++) {
            c_psiBound(side, angle, fvrtx, group) = localFaceData(fvrtx, group);
        }}
    }
    
    
    /*
        getPriority
        
        Return priority for the cell/angle pair.
    */
    virtual UINT getPriority(UINT cell, UINT angle)
    {
        return c_priorities(cell, angle);
    }
    
    
    /*
        update
        
        Does a transport update for the given cell/angle pair.
    */
    virtual void update(UINT cell, UINT angle, 
                        UINT adjCellsSides[g_nFacePerCell], 
                        BoundaryType bdryType[g_nFacePerCell])
    {
        UNUSED_VARIABLE(adjCellsSides);
        UNUSED_VARIABLE(bdryType);
        
        Mat2<double> localSource(g_nVrtxPerCell, g_nGroups);
        Mat2<double> localPsi(g_nVrtxPerCell, g_nGroups);
        Mat3<double> localPsiBound(g_nVrtxPerFace, g_nFacePerCell, g_nGroups);
        
        
        // Populate localSource
        for (UINT group = 0; group < g_nGroups; group++) {
        for (UINT vrtx = 0; vrtx < g_nVrtxPerCell; vrtx++) {
            localSource(vrtx, group) = c_source(vrtx, angle, cell, group);
        }}
        
        
        // Populate localPsiBound
        populateLocalPsiBound(angle, cell, c_psi, c_psiBound, localPsiBound);
        
        
        // Transport solve
        Transport::solve(cell, angle, c_sigmaTotal, 
                         localPsiBound, localSource, localPsi);
        
        
        // localPsi -> psi
        for (UINT group = 0; group < g_nGroups; group++) {
        for (UINT vrtx = 0; vrtx < g_nVrtxPerCell; vrtx++) {
            c_psi(vrtx, angle, cell, group) = localPsi(vrtx, group);
        }}
    }
    
    
private:
    PsiData &c_psi;
    PsiData c_psiBound;
    const PsiData &c_source;
    const double c_sigmaTotal;
    const Mat2<UINT> &c_priorities;
    vector<Mat2<double>> c_localFaceData;
};


/*
    Sweeper2 constructor
*/
Sweeper2::Sweeper2(const UINT maxComputePerStep,
                   const UINT intraAngleP, 
                   const UINT interAngleP, 
                   const double sigmaTotal)
{
    c_maxComputePerStep = maxComputePerStep;
    c_sigmaTotal = sigmaTotal;
    c_priorities = 
        Mat2<UINT>(g_spTychoMesh->getNCells(), g_quadrature->getNumAngles());
    Priorities::calcPriorities(c_maxComputePerStep, intraAngleP, interAngleP, 
                               c_priorities);
}


/*
    Sweeper2::sweep
    
    Sweep by traversing graph.
*/
void Sweeper2::sweep(PsiData &psi, const PsiData &source)
{
    const bool doComm = true;
    SweepData sweepData(psi, source, c_sigmaTotal, c_priorities);
    traverseGraph(c_maxComputePerStep, sweepData, doComm, MPI_COMM_WORLD, 
                  Direction_Forward);
}

