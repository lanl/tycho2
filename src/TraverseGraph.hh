#ifndef __TRAVERSE_GRAPH_HH__
#define __TRAVERSE_GRAPH_HH__

#include "Global.hh"
#include <mpi.h>
#include <vector>

/*
    Boundary Type for faces of a cell.
    They are split into incoming and outgoing wrt sweep direction.
    Interior means adjacent cell is on same proc.
    Interior boundary means adj cell is on different proc.
    Exterior boundary means it is a boundary for the whole mesh.
*/
enum BoundaryType
{
    // Outgoing
    BoundaryType_OutIntBdry,    // Interior boundary
    BoundaryType_OutExtBdry,    // Exterior boundary
    BoundaryType_OutInt,        // Interior
    
    // Incoming
    BoundaryType_InIntBdry,     // Interior boundary
    BoundaryType_InExtBdry,     // Exterior boundary
    BoundaryType_InInt          // Interior
};


enum Direction
{
    Direction_Forward,
    Direction_Backward
};


/*
    TraverseData class
    
    Abstract class defining the methods needed to traverse a graph.
*/
class TraverseData
{
public:
    virtual const char* getData(UINT cell, UINT face, UINT angle) = 0;
    virtual size_t getDataSize() = 0;
    virtual void setSideData(UINT side, UINT angle, const char *data) = 0;
    virtual UINT getPriority(UINT cell, UINT angle) = 0;
    virtual void update(UINT cell, UINT angle, 
                        UINT adjCellsSides[g_nFacePerCell], 
                        BoundaryType bdryType[g_nFacePerCell]) = 0;

protected:
    // Don't allow construction of this base class.
    TraverseData() { }
};


void traverseGraph(const UINT maxComputePerStep,
                   TraverseData &traverseData, bool doComm,
                   MPI_Comm mpiComm, Direction direction);

#endif