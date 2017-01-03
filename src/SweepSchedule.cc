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

#include "SweepSchedule.hh"
#include "Global.hh"
#include "Comm.hh"
#include "TychoMesh.hh"
#include "Mat.hh"
#include "Priorities.hh"

#include <cmath>
#include <algorithm>
#include <algorithm>
#include <utility>
#include <iostream>
#include <set>
#include <queue>


using namespace std;
const static double TOL = 1.e-10;


/*
    PriorityWork class
    
    Specialization of Work that contains a priority of solution.
    Used to determine a good sweep ordering.
*/
class PriorityWork
{
private:
    UINT c_priority;
    UINT c_cell;
    UINT c_angle;
    
public:
    PriorityWork(UINT cell, UINT angle, double priority)
        : c_priority(priority), c_cell(cell), c_angle(angle) {}
    
    UINT getCell() const { return c_cell; }
    UINT getAngle() const { return c_angle; }
    
    // Comparison operator to determine relative priorities
    // Needed for priority_queue
    bool operator<(const PriorityWork &rhs) const
        {return c_priority < rhs.c_priority;}
};


/*
    calcNeighborProcs
*/
static
void calcNeighborProcs(set<UINT> &neighborProcs)
{
    for (UINT cell = 0; cell < g_nCells; cell++) {
    for (UINT face = 0; face < g_nFacePerCell; face++) {
        UINT proc = g_tychoMesh->getAdjRank(cell, face);
        UINT adjCell = g_tychoMesh->getAdjCell(cell, face);
        if (adjCell == TychoMesh::BOUNDARY_FACE && proc != TychoMesh::BAD_RANK) {
            neighborProcs.insert(proc);
        }
    }}
}


/*
    calcDependencies
*/
static
void calcDependencies(const vector<UINT> &angles, 
                      Mat3<bool> &hasChild, Mat3<bool> &hasParent)
{
    for (unsigned angle = 0; angle < angles.size(); ++angle) {
    for (UINT cell = 0; cell < g_nCells; ++cell) {
    for (UINT face = 0; face < g_nFacePerCell; ++face) {
        hasChild(angle, cell, face) = false;
        
        if (g_tychoMesh->getOmegaDotN(angles[angle], cell, face) > TOL) {
            if (g_tychoMesh->getAdjRank(cell, face) != TychoMesh::BAD_RANK)
                hasChild(angle, cell, face) = true;
        }
        else if (g_tychoMesh->getOmegaDotN(angles[angle], cell, face) < -TOL) {
            if (g_tychoMesh->getAdjRank(cell, face) != TychoMesh::BAD_RANK)
                hasParent(angle, cell, face) = true;
        }
    }}}
}


/*
    calcNumDependents
*/
static
void calcNumDependents(UINT numAngles, const Mat3<bool> &hasChild, 
                       Mat2<UINT> &nDependents)
{
    for (UINT angle = 0; angle < numAngles; ++angle) {
    for (UINT cell = 0; cell < g_nCells; ++cell) {
    for (UINT face = 0; face < g_nFacePerCell; ++face) {
        if (hasChild(angle, cell, face))
            ++nDependents(angle, cell);
    }}}
}


/*
    initializeWorkQ
*/
static
void initializeWorkQ(const UINT numAngles, const Mat2<UINT> &nNeeded,
                     const Mat2<UINT> &priorities,
                     priority_queue<PriorityWork> &initialWork) 
{
    for (UINT angle = 0; angle < numAngles; ++angle) {
    for (UINT cell = 0; cell < g_nCells; ++cell) {
        if (nNeeded(angle, cell) == 0) {
            initialWork.push(PriorityWork(cell, angle, priorities(angle, cell)));
        }
    }}
}


/*
    allDone
*/
static
bool allDone(const priority_queue<PriorityWork> &availableWork)
{
    bool done;
    UINT doneInt;
    
    if (availableWork.empty())
        doneInt = 0;
    else
        doneInt = 1;

    Comm::gmax(doneInt);
    if (doneInt == 1)
        done = false;
    else
        done = true;

    return done;
}


/*
    partialTopoSort
*/
static
void partialTopoSort(priority_queue<PriorityWork> &availableWork, 
                     Mat2<UINT> &nNeeded,
                     const Mat3<bool> &hasChild, 
                     const Mat2<UINT> &priorities, 
                     const unsigned maxCellsPerStep,
                     const vector<UINT> &angles,
                     vector<SweepSchedule::Work> &workDone)
{
    while (!availableWork.empty() && workDone.size() < maxCellsPerStep) {
        PriorityWork cellWork = availableWork.top();
        availableWork.pop();
        UINT cell = cellWork.getCell();
        UINT angle = cellWork.getAngle();
        workDone.push_back(SweepSchedule::Work(cell, angles[angle]));

        for (UINT face = 0; face < g_nFacePerCell; ++face) {
            if (hasChild(angle, cell, face)) {
                UINT childCell = g_tychoMesh->getAdjCell(cell, face);
                
                // on processor
                if (childCell != TychoMesh::BOUNDARY_FACE) {
                    --nNeeded(angle, childCell);
                    if (nNeeded(angle, childCell) == 0) {
                        availableWork.push(PriorityWork(childCell, angle,
                                           priorities(angle, childCell)));
                    }
                }
            }
        }
    }
}


/*
    glAngle
*/
static
UINT glAngle(const vector<UINT> &angles, const UINT gAngle)
{
    unsigned angle = 0;
    for (angle = 0; angle < angles.size(); ++angle) {
        if(angles[angle] == gAngle)
            break;
    }
    
    return angle;
}


/*
    updateLevels
*/
static
void updateLevels(const vector<SweepSchedule::Work> &workDone,
                  const vector<UINT> &angles, 
                  Mat2<UINT> &cellLevels,
                  const Mat3<bool> &hasParent)
{
    for (const SweepSchedule::Work &work : workDone) {
        UINT cell = work.getCell();
        UINT angle = glAngle(angles, work.getAngle());
        for (UINT face = 0; face < g_nFacePerCell; ++face) {
            UINT adjCell = g_tychoMesh->getAdjCell(cell, face);
            if (hasParent(angle, cell, face) &&
                adjCell != TychoMesh::BOUNDARY_FACE) // on processor
            {
                cellLevels(angle, adjCell) = 
                    max(cellLevels(angle, adjCell), cellLevels(angle, cell)+1);
            }
        }
    }
}


/*
    sendLevels
*/
static
void sendLevels(const vector<SweepSchedule::Work> &workDone,
                const vector<UINT> &angles, 
                const Mat2<UINT> &cellLevels,
                const Mat3<bool> &hasParent,
                const set<UINT> &neighborProcs)
{
    vector<vector<UINT> > commSidesAnglesLevels(Comm::numRanks());
    
    for (const SweepSchedule::Work &work : workDone) {
        UINT cell = work.getCell();
        UINT angle = glAngle(angles, work.getAngle());
        for (UINT face = 0; face < g_nFacePerCell; ++face) {
            UINT adjCell = g_tychoMesh->getAdjCell(cell, face);
            if (hasParent(angle, cell, face) &&
                adjCell == TychoMesh::BOUNDARY_FACE) // off processor
            {
                UINT side = g_tychoMesh->getSide(cell, face);
                UINT gSide = g_tychoMesh->getLGSide(side);
                UINT neighborProc = g_tychoMesh->getAdjRank(cell, face);
                commSidesAnglesLevels[neighborProc].push_back(gSide);
                commSidesAnglesLevels[neighborProc].push_back(angle);
                commSidesAnglesLevels[neighborProc].push_back
                    (cellLevels(angle, cell));
            }
        }
    }
    
    for(UINT proc : neighborProcs) {
        UINT nSides = commSidesAnglesLevels[proc].size()/3;
        Comm::sendUInt(nSides, proc);
        if (nSides > 0) {
            Comm::sendUIntVector(commSidesAnglesLevels[proc], proc);
        }
    }
}


/*
    recvLevels
*/
static
void recvLevels(Mat2<UINT> &cellLevels, Mat2<UINT> &sideLevels,
                   priority_queue<PriorityWork> &availableWork,
                   const Mat2<UINT> &priorities,
                   Mat2<UINT> &nNeeded,
                   const set<UINT> &neighborProcs)
{
    for (UINT proc : neighborProcs) {
        UINT nSides;
        Comm::recvUInt(nSides, proc);
        if (nSides > 0) {
            vector<UINT> commSidesAnglesLevels(3*nSides);
            Comm::recvUIntVector(commSidesAnglesLevels, proc);

            for (UINT child = 0; child < nSides; ++child) {
                const UINT globalSide = commSidesAnglesLevels[child * 3];
                const UINT angle = commSidesAnglesLevels[child * 3 + 1];
                const UINT level = commSidesAnglesLevels[child * 3 + 2];
                const UINT localSide = g_tychoMesh->getGLSide(globalSide);
                const UINT cell = g_tychoMesh->getSideCell(localSide);
                
                --nNeeded(angle, cell);  // dependencies must be in sync
                if (nNeeded(angle, cell) == 0) {
                    availableWork.push(PriorityWork(cell, angle, 
                            priorities(angle, cell)));
                }
                cellLevels(angle, cell) = max(cellLevels(angle, cell), level+1);
                sideLevels(angle, localSide) = level;
            }
        }
    }
}


/*
    sendOrders
*/
static
void sendOrders(const vector<SweepSchedule::Work> &workDone,
                const vector<UINT> &angles, 
                const Mat3<bool> &hasChild,
                const set<UINT> &neighborProcs,
                vector<vector<UINT> > &sendProcs)
{
    vector<vector<UINT> > commSidesAngles(Comm::numRanks());
    
    for (const SweepSchedule::Work &work : workDone) {
        UINT cell = work.getCell();
        UINT angle = glAngle(angles, work.getAngle());
        for (UINT face = 0; face < g_nFacePerCell; ++face) {
            UINT adjCell = g_tychoMesh->getAdjCell(cell, face);
            if (hasChild(angle, cell, face) &&
                adjCell == TychoMesh::BOUNDARY_FACE) // off processor
            {
                UINT side = g_tychoMesh->getSide(cell, face);
                UINT gSide = g_tychoMesh->getLGSide(side);
                UINT neighborProc = g_tychoMesh->getAdjRank(cell, face);
                commSidesAngles[neighborProc].push_back(gSide);
                commSidesAngles[neighborProc].push_back(angle);
            }
        }
    }

    vector<UINT> stepSendProcs;
    for (UINT proc : neighborProcs) {
        UINT nSides = commSidesAngles[proc].size()/2;
        Comm::sendUInt(nSides, proc);
        if (nSides > 0) {
            Comm::sendUIntVector(commSidesAngles[proc], proc);
            stepSendProcs.push_back(proc);
        }
    }
    sendProcs.push_back(stepSendProcs);
}


/*
    recvOrders
*/
static
void recvOrders(priority_queue<PriorityWork> &availableWork,
                const Mat2<UINT> &priorities,
                Mat2<UINT> &nNeeded,
                const set<UINT> &neighborProcs,
                vector<vector<UINT> > &recvProcs)
{
    vector<UINT> stepRecvProcs;
    
    for (UINT proc : neighborProcs) {
        UINT nSides;
        Comm::recvUInt(nSides, proc);
        if (nSides > 0) {
            vector<UINT> commSidesAngles(2*nSides);
            Comm::recvUIntVector(commSidesAngles, proc);
            stepRecvProcs.push_back(proc);
            
            for (UINT parent = 0; parent < nSides; ++parent) {
                const UINT globalSide = commSidesAngles[2 * parent];
                const UINT localSide = g_tychoMesh->getGLSide(globalSide);
                const UINT cell = g_tychoMesh->getSideCell(localSide);
                const UINT angle = commSidesAngles[2 * parent + 1];
                --nNeeded(angle, cell);  // dependencies must be in sync
                if (nNeeded(angle, cell) == 0) {
                    availableWork.push(PriorityWork(cell, 
                            angle, priorities(angle, cell)));
                }
            }
        }
    }

    recvProcs.push_back(stepRecvProcs);
}


/*
    calcOrdering
*/
static
void calcOrdering(const Mat3<bool> &hasChild,
                  const Mat3<bool> &hasParent,
                  const Mat2<UINT> &priorities,
                  const set<UINT> &neighborProcs,
                  const UINT numAngles,
                  const UINT maxCellsPerStep, 
                  const vector<UINT> &angles,
                  vector<vector<SweepSchedule::Work> > &workOrders,
                  vector<vector<UINT> > &sendProcs,
                  vector<vector<UINT> > &recvProcs)
{
    Mat2<UINT> nParentsNeeded(numAngles, g_nCells);
    calcNumDependents(numAngles, hasParent, nParentsNeeded);
    priority_queue<PriorityWork> availableWork;
    initializeWorkQ(numAngles, nParentsNeeded, priorities, availableWork);

    UINT step = 0;
    UINT elapsed = 0;
    while (!allDone(availableWork)) {
        ++step;
        vector<SweepSchedule::Work> workDone;
        partialTopoSort(availableWork, nParentsNeeded, hasChild, priorities, 
                        maxCellsPerStep, angles, workDone);
        workOrders.push_back(workDone);
        sendOrders(workDone, angles, hasChild, neighborProcs, sendProcs);
        recvOrders(availableWork, priorities, nParentsNeeded, neighborProcs, 
                   recvProcs);

        UINT stepSize = workDone.size();
        Comm::gmax(stepSize);
        elapsed += stepSize;
    }
    
    UINT totalNCells = g_nCells;
    Comm::gsum(totalNCells);
    if (Comm::rank() == 0) {
        printf("   Number of sweep steps: %" PRIu64 "\n", step);
        double stepEff = static_cast<double>
            (numAngles * totalNCells) / (Comm::numRanks() * elapsed);
        printf("   Step efficiency: %f\n", stepEff);
    }
}


/*
    calcLevels
*/
static
UINT calcLevels(Mat2<UINT> &cellLevels, 
               Mat2<UINT> &sideLevels,
               const Mat3<bool> &hasChild, 
               const Mat3<bool> &hasParent,
               const set<UINT> &neighborProcs, 
               const vector<UINT> angles, 
               const UINT numAngles, 
               const UINT maxCellsPerStep)
{
    Mat2<UINT> nChildrenNeeded(numAngles, g_nCells);
    calcNumDependents(numAngles, hasChild, nChildrenNeeded);
    Mat2<UINT> priorities(numAngles, g_nCells);
    priorities.setAll(0.0);
    
    priority_queue<PriorityWork> availableWork;
    initializeWorkQ(numAngles, nChildrenNeeded, priorities, availableWork);

    // This may fail if dependencies are not in sync between processors
    while (!allDone(availableWork)) {
        vector<SweepSchedule::Work> workDone;
        partialTopoSort(availableWork, nChildrenNeeded, hasParent, priorities, 
                        maxCellsPerStep, angles, workDone);
        updateLevels(workDone, angles, cellLevels, hasParent);
        sendLevels(workDone, angles, cellLevels, hasParent, neighborProcs);
        recvLevels(cellLevels, sideLevels, availableWork, priorities,
                   nChildrenNeeded, neighborProcs);
    }
    
    UINT nlevels = cellLevels[0];
    for(unsigned i = 1; i < cellLevels.size(); i++) {
        nlevels = std::max(nlevels, cellLevels[i]);
    }
    Comm::gmax(nlevels);
    
    return nlevels;
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
        priorities(angle, cell) = rand();
    }}
}


/*
    levelPriorities
*/
static
void levelPriorities(const Mat2<UINT> &cellLevels, const UINT numAngles,
                     Mat2<UINT> &priorities)
{
    for (UINT angle = 0; angle < numAngles; ++angle) {
    for (UINT cell = 0; cell < g_nCells; ++cell) {
        priorities(angle, cell) = cellLevels(angle, cell);
    }}
}


/*
    updatePriorities
*/
static
void updatePriorities(Mat2<UINT> &priorities,
                      const vector<SweepSchedule::Work> &workDone,
                      const vector<UINT> &angles, 
                      const Mat3<bool> &hasParent,
                      const UINT parentShift)
{
    for (const SweepSchedule::Work &work : workDone) {
        UINT cell = work.getCell();
        UINT angle = glAngle(angles, work.getAngle());
        for (UINT face = 0; face < g_nFacePerCell; ++face) {
            UINT adjCell = g_tychoMesh->getAdjCell(cell, face);
            if (hasParent(angle, cell, face) &&
                adjCell != TychoMesh::BOUNDARY_FACE) // on processor
            {
                priorities(angle, adjCell) =
                    max(priorities(angle, adjCell),
                        priorities(angle, cell) + parentShift);
            }
        }
    }
}


/*
    neighborPriorities
*/
static
void neighborPriorities(const Mat2<UINT> &sideLevels, 
                        const Mat3<bool> &hasChild,
                        const Mat3<bool> &hasParent, 
                        const vector<UINT> &angles, 
                        const UINT boundScale, 
                        const UINT boundShift, 
                        const int parentShift,  // Can be -1 or 0
                        const UINT numAngles, 
                        const UINT maxCellsPerStep,
                        Mat2<UINT> &priorities)
{
    Mat2<UINT> nChildrenNeeded(numAngles, g_nCells);
    calcNumDependents(numAngles, hasChild, nChildrenNeeded);

    for (UINT angle = 0; angle < numAngles; ++angle) {
    for (UINT cell = 0; cell < g_nCells; ++cell) {
    for (UINT face = 0; face < g_nFacePerCell; ++face) {
        if (hasChild(angle, cell, face) &&
            g_tychoMesh->getAdjCell(cell, face) == TychoMesh::BOUNDARY_FACE)
        {
            UINT side = g_tychoMesh->getSide(cell, face);

            // Internal boundary
            if (g_tychoMesh->getAdjRank(cell, face) != TychoMesh::BAD_RANK)
            {
                priorities(angle, cell) =
                    max(priorities(angle, cell), 
                    (sideLevels(angle, side)*boundScale + boundShift));
                --nChildrenNeeded(angle, cell);
            }
        }
    }}}

    Mat2<UINT> dummyPriorities(numAngles, g_nCells);
    dummyPriorities.setAll(0.0);
    priority_queue<PriorityWork> availableWork;
    initializeWorkQ(numAngles, nChildrenNeeded, dummyPriorities, availableWork);
    UINT nSolved = 0;
    while (nSolved != numAngles * g_nCells)
    {
        vector<SweepSchedule::Work> workDone;
        partialTopoSort(availableWork, nChildrenNeeded, hasParent, 
                        dummyPriorities, maxCellsPerStep, angles, workDone);
        nSolved += workDone.size();
        updatePriorities(priorities, workDone, angles, hasParent, parentShift);
    }
}


/*
    anglePriorities
*/
static
void anglePriorities(Mat2<UINT> &priorities, const UINT numAngles, 
                     const UINT interAngleP, const UINT nlevels)
{
    vector<UINT> highPriorities(numAngles, 0);
    for (UINT angle = 0; angle < numAngles; ++angle) {
    for (UINT cell = 0; cell < g_nCells; ++cell) {
        highPriorities[angle] =
            max(highPriorities[angle], priorities(angle, cell));
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
                priorities(angle, cell) += angle * nlevels * ncells;
            }}
            break;
        case 2:  // locally prioritized angles
        {
            vector<pair<double, UINT> > orderedAngles;
            for (UINT angle = 0; angle < numAngles; ++angle) {
                orderedAngles.push_back(make_pair(highPriorities[angle], angle));
            }
            std::sort(orderedAngles.begin(), orderedAngles.end());
            UINT order = 0;
            for (vector<pair<double, UINT> >::reverse_iterator iter =
                 orderedAngles.rbegin();
                 iter != orderedAngles.rend(); ++iter)
            {
                ++order;
                UINT angle = iter->second;
                for (UINT cell = 0; cell < g_nCells; ++cell) {
                    priorities(angle, cell) +=
                        (numAngles-order) * nlevels * ncells;
                }
            }
            break;
        }
    }
}


/*
    Constructor
    Note: function to break cyclic dependencies is not implemented
*/
SweepSchedule::SweepSchedule(const std::vector<UINT> &angles,
                             const UINT maxCellsPerStep,
                             const UINT intraAngleP, 
                             const UINT interAngleP)
{
    // Get processors neighboring this processor
    set<UINT> neighborProcs;
    calcNeighborProcs(neighborProcs);
    
    
    // Get which faces are parents/children for specific angles
    Mat3<bool> hasChild(angles.size(), g_nCells,
                        g_nFacePerCell);
    Mat3<bool> hasParent(angles.size(), g_nCells,
                         g_nFacePerCell);
    calcDependencies(angles, hasChild, hasParent);
    
    
    // Calculate B-Levels
    Mat2<UINT> cellLevels(angles.size(), g_nCells);
    Mat2<UINT> sideLevels(angles.size(), g_tychoMesh->getNSides());
    UINT nlevels = calcLevels(cellLevels, sideLevels, hasChild, hasParent, 
                             neighborProcs, angles, angles.size(), maxCellsPerStep);
    
    
    // Calculate intra-angle priorities
    Mat2<UINT> priorities(angles.size(), g_nCells);
    switch (intraAngleP)
    {
      case 0:  // random
        randomPriorities(angles.size(), priorities);
        break;
      case 1:  // directed graph levels
        levelPriorities(cellLevels, angles.size(), priorities);
        break;
      case 2:  // breadth-first dependent seeking
        neighborPriorities(sideLevels, hasChild, hasParent, angles, 
                           1, 0, 0, 
                           angles.size(), maxCellsPerStep, priorities);
        break;
      case 3:  // depth-first dependent seeking
        neighborPriorities(sideLevels, hasChild, hasParent, angles, 
                           1, nlevels, -1, 
                           angles.size(), maxCellsPerStep, priorities);
        break;
      case 4:  // strict depth-first dependent seeking
        neighborPriorities(sideLevels, hasChild, hasParent, angles, 
                           nlevels, nlevels, -1, 
                           angles.size(), maxCellsPerStep, priorities);
        break;
    }
    
    
    // Calculate inter-angle priorities
    anglePriorities(priorities, angles.size(), interAngleP, nlevels);
    
    
    // Calculate the Ordering
    calcOrdering(hasChild, hasParent, priorities, neighborProcs, angles.size(), 
                 maxCellsPerStep, angles, c_workOrders, c_sendProcs, c_recvProcs);
}

