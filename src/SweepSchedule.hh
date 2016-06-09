//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   topologicalSort/SweepSchedule.hh
 * \author Shawn Pautz
 * \date   Mon Apr  3 12:49:24 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id: SweepSchedule.hh,v 1.4 2002/04/12 16:15:39 pautz Exp $
//---------------------------------------------------------------------------//

#ifndef __topologicalSort_SweepSchedule_hh__
#define __topologicalSort_SweepSchedule_hh__

#include "Typedef.hh"
#include <vector>


class SweepSchedule 
{
  public:
    
    // Class holding a unit of work for the sweeper
    class Work
    {
      private:
        UINT c_cell;
        UINT c_angle;

      public:
        Work(UINT cell, UINT angle) : c_cell(cell), c_angle(angle) {}
        UINT getCell() const { return c_cell; }
        UINT getAngle() const { return c_angle; }
    };
    
    
    // SweepSchedule implementation
    SweepSchedule(const std::vector<UINT> &angles,
                  const UINT maxCellsPerStep,
                  const UINT intraAngleP,
                  const UINT interAngleP);

    UINT nSteps() const 
        { return c_workOrders.size(); }
    const std::vector<Work>& getWork(const UINT step) const
        { return c_workOrders[step]; }
    const std::vector<UINT>& getSendProcs(const UINT step) const
        { return c_sendProcs[step]; }
    const std::vector<UINT>& getRecvProcs(const UINT step) const
        { return c_recvProcs[step]; }
    
    
  private:
    // work to be performed in each step
    std::vector<std::vector<Work> > c_workOrders;
    // processors to send data to after each step
    std::vector<std::vector<UINT> > c_sendProcs;
    // processors to receive data from after each step
    std::vector<std::vector<UINT> > c_recvProcs;
};


#endif                          // __topologicalSort_SweepSchedule_hh__

//---------------------------------------------------------------------------//
//                              end of topologicalSort/SweepSchedule.hh
//---------------------------------------------------------------------------//
