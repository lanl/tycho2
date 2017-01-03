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

#include "Priorities.hh"
#include "GraphTraverser.hh"
#include "Mat.hh"
#include "Global.hh"
#include "TychoMesh.hh"
#include "Comm.hh"
#include <vector>
#include <algorithm>

using namespace std;


/*
    BLevelData
    
    Calculates b-levels when traversing a graph.
    Note: OpenMP assumes threading across angle.
          Without this assumption, there could be a race condition.
*/
class BLevelData : public TraverseData
{
public:
    
    /*
        Constructor
    */
    BLevelData(Mat2<UINT> &bLevels, Mat2<UINT> &sideBLevels)
    : c_bLevels(bLevels), c_sideBLevels(sideBLevels)
    {
        // Initialize all to 0 b-level
        c_maxBLevel = 0;
        bLevels.setAll(0);
        sideBLevels.setAll(0);
        c_bLevels.setAll(0);
        c_sideBLevels.setAll(0);
    }
    
    
    /*
        getData
        
        Return b-level data given (cell, angle) pair.
    */
    virtual const char* getData(UINT cell, UINT face, UINT angle)
    {
        UNUSED_VARIABLE(face);
        return (char*) (&c_bLevels(cell, angle));
    }
    
    
    /*
        setSideData
        
        Return b-level data given (side, angle) pair.
    */
    virtual void setSideData(UINT side, UINT angle, const char *data)
    {
        c_sideBLevels(side, angle) = *((UINT*)(data));
    }
    
    
    /*
        getPriority
        
        Return a priority for the cell/angle pair.
        Not needed for this class, so it is just set to a constant.
    */
    virtual UINT getPriority(UINT cell, UINT angle)
    {
        UNUSED_VARIABLE(cell);
        UNUSED_VARIABLE(angle);
        return 1;
    }
    
    
    /*
        update
        
        Updates b-level information for a given (cell, angle) pair.
    */
    virtual void update(UINT cell, UINT angle, 
                        UINT adjCellsSides[g_nFacePerCell], 
                        BoundaryType bdryType[g_nFacePerCell])
    {
        c_bLevels(cell, angle) = 0;
        for (UINT face = 0; face < g_nFacePerCell; face++) {
            
            if (bdryType[face] == BoundaryType_OutIntBdry) {
                UINT adjSide = adjCellsSides[face];
                UINT adjBLevel = c_sideBLevels(adjSide, angle);
                c_bLevels(cell, angle) = max(c_bLevels(cell, angle), adjBLevel + 1);
            }
            
            else if (bdryType[face] == BoundaryType_OutInt) {
                UINT adjCell = adjCellsSides[face];
                UINT adjBLevel = c_bLevels(adjCell, angle);
                c_bLevels(cell, angle) = max(c_bLevels(cell, angle), adjBLevel + 1);
            }
        }
        
        #pragma omp critical
        {
            c_maxBLevel = max(c_maxBLevel, c_bLevels(cell, angle));
        }
    }
    
    
    /*
        getMaxBLevel
    */
    virtual UINT getMaxBLevel()
    {
        UINT maxBLevel;
        
        #pragma omp atomic read
        maxBLevel = c_maxBLevel;
        
        return maxBLevel;
    }
    
private:
    Mat2<UINT> &c_bLevels;
    Mat2<UINT> &c_sideBLevels;
    UINT c_maxBLevel;
};


/*
    NeighborPriorityData
    
    Calculates priorities when traversing a graph.
    Can calculate BFDS, DFDS, or DFHDS depending on 
    boundScale, boundShift, parentShift.
    
            boundScale   boundShift   parentShift
    BFDS  =     1             0             0
    DFDS  =     1         maxBLevel        -1
    DFHDS =  maxBLevel    maxBLevel        -1
    
    Note: OpenMP assumes threading across angle.
          Without this assumption, there could be a race condition.
*/
class NeighborPriorityData : public TraverseData
{
public:
    
    /*
        Constructor
    */
    NeighborPriorityData(Mat2<UINT> &priorities, const Mat2<UINT> &sideBLevels, 
                         const UINT boundScale, const UINT boundShift, 
                         const int parentShift)
    : c_priorities(priorities), c_sideBLevels(sideBLevels), 
      c_boundScale(boundScale), c_boundShift(boundShift), c_parentShift(parentShift)
    {
        c_priorities.setAll(0);
    }
    
    
    /*
        data
        
        Return priority data given (cell, angle) pair.
    */
    virtual const char* getData(UINT cell, UINT face, UINT angle)
    {
        UNUSED_VARIABLE(face);
        return (char*) (&c_priorities(cell, angle));
    }
    
    
    /*
        getPriority
        
        Return a priority for the cell/angle pair.
        Not needed for this class, so it is just set to a constant.
    */
    virtual UINT getPriority(UINT cell, UINT angle)
    {
        UNUSED_VARIABLE(cell);
        UNUSED_VARIABLE(angle);
        return 1;
    }
    
    
    /*
        update
        
        Updates priority information for a given (cell, angle) pair.
    */
    virtual void update(UINT cell, UINT angle, 
                        UINT adjCellsSides[g_nFacePerCell], 
                        BoundaryType bdryType[g_nFacePerCell])
    {
        c_priorities(cell, angle) = 0;
        
        for (UINT face = 0; face < g_nFacePerCell; face++) {
            
            UINT priority = 0;
            
            if (bdryType[face] == BoundaryType_OutIntBdry) {
                UINT adjSide = adjCellsSides[face];
                priority = c_sideBLevels(adjSide, angle) * c_boundScale + 
                           c_boundShift;
            }
            
            else if (bdryType[face] == BoundaryType_OutInt) {
                UINT adjCell = adjCellsSides[face];
                priority = c_priorities(adjCell, angle) + c_parentShift;
            }
            
            c_priorities(cell, angle) = max(c_priorities(cell, angle), priority);
        }
    }
    
    
    /*
        These should never be called.
        They are only used when communication is involved in traversing the
        graph.
    */
    virtual void setSideData(UINT side, UINT angle, const char *data)
    {
        UNUSED_VARIABLE(side);
        UNUSED_VARIABLE(angle);
        UNUSED_VARIABLE(data);
        Assert(false);
    }
    
private:
    Mat2<UINT> &c_priorities;
    const Mat2<UINT> &c_sideBLevels;
    const UINT c_boundScale;
    const UINT c_boundShift;
    const int c_parentShift;
};


/*
    calcBLevels
*/
static
UINT calcBLevels(Mat2<UINT> &bLevels, Mat2<UINT> &sideBLevels,
                 GraphTraverser *graphTraverser)
{
    BLevelData bLevelData(bLevels, sideBLevels);
    graphTraverser->traverse(g_maxCellsPerStep, bLevelData);
    
    UINT maxBLevel = bLevelData.getMaxBLevel();
    Comm::gmax(maxBLevel);
    return maxBLevel;
}


/*
    randomPriorities
*/
static 
void randomPriorities(const UINT numAngles, Mat2<UINT> &priorities)
{
    srand(0);
    for (UINT angle = 0; angle < numAngles; ++angle) {
    for (UINT cell = 0; cell < g_nCells; ++cell) {
        priorities(cell, angle) = rand();
    }}
}


/*
    levelPriorities
*/
static 
void levelPriorities(const Mat2<UINT> &bLevels, const UINT numAngles, 
                     Mat2<UINT> &priorities)
{
    for (UINT angle = 0; angle < numAngles; ++angle) {
    for (UINT cell = 0; cell < g_nCells; ++cell) {
        priorities(cell, angle) = bLevels(cell, angle);
    }}
}


/*
    neighborPriorities
*/
static 
void neighborPriorities(const Mat2<UINT> &sideBLevels, 
                        const UINT boundScale, 
                        const UINT boundShift, 
                        const int parentShift,  // Can be -1 or 0
                        Mat2<UINT> &priorities,
                        GraphTraverser *graphTraverser)
{
    NeighborPriorityData priorityData(priorities, sideBLevels, 
                                      boundScale, boundShift, parentShift);
    graphTraverser->traverse(g_maxCellsPerStep, priorityData);
}


/*
    anglePriorities
*/
static
void anglePriorities(const UINT numAngles, const UINT interAngleP, 
                     const UINT nlevels, Mat2<UINT> &priorities)
{
    vector<UINT> highPriorities(numAngles, 0);
    for (UINT angle = 0; angle < numAngles; ++angle) {
    for (UINT cell = 0; cell < g_nCells; ++cell) {
        highPriorities[angle] =
            max(highPriorities[angle], priorities(cell, angle));
    }}

    UINT ncells = g_nCells;
    Comm::gsum(ncells);

    switch (interAngleP)
    {
        case 0:  // interleaved
            break;
        case 1:  // globally prioritized angles
            for (UINT angle = 0; angle < numAngles; ++angle) {
            for (UINT cell = 0; cell < g_nCells; ++cell) {
                priorities(cell, angle) += angle * nlevels * ncells;
            }}
            break;
        case 2:  // locally prioritized angles
        {
            vector<pair<UINT, UINT> > orderedAngles;
            for (UINT angle = 0; angle < numAngles; ++angle) {
                orderedAngles.push_back(make_pair(highPriorities[angle], angle));
            }
            std::sort(orderedAngles.begin(), orderedAngles.end());
            UINT order = 0;
            for (auto iter = orderedAngles.rbegin();
                 iter != orderedAngles.rend(); ++iter)
            {
                ++order;
                UINT angle = iter->second;
                for (UINT cell = 0; cell < g_nCells; ++cell) {
                    priorities(cell, angle) +=
                        (numAngles-order) * nlevels * ncells;
                }
            }
            break;
        }
    }
}


namespace Priorities
{

/*
    calcPriorities
*/
void calcPriorities(Mat2<UINT> &priorities)
{
    const bool doComm = false;
    GraphTraverser graphTraverser(Direction_Backward, doComm, sizeof(UINT));
    
    UINT numAngles = g_nAngles;
    Mat2<UINT> bLevels(g_nCells, numAngles);
    Mat2<UINT> sideBLevels(g_tychoMesh->getNSides(), numAngles);
    
    UINT maxBLevel = calcBLevels(bLevels, sideBLevels, 
                                 &graphTraverser);
    
    
    // Calculate intra-angle priorities
    switch (g_intraAngleP)
    {
      case 0:  // random
        randomPriorities(numAngles, priorities);
        break;
      case 1:  // b-levels
        levelPriorities(bLevels, numAngles, priorities);
        break;
      case 2:  // breadth-first dependent seeking
        neighborPriorities(sideBLevels, 1, 0, 0, 
                           priorities, &graphTraverser);
        break;
      case 3:  // depth-first dependent seeking
        neighborPriorities(sideBLevels, 1, maxBLevel, -1, 
                           priorities, &graphTraverser);
        break;
      case 4:  // strict depth-first dependent seeking
        neighborPriorities(sideBLevels, maxBLevel, maxBLevel, -1, 
                           priorities, &graphTraverser);
        break;
    }
    
    
    // Calculate inter-angle priorities
    anglePriorities(numAngles, g_interAngleP, maxBLevel, priorities);
}

} // End namespace


