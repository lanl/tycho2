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

#ifndef __SWEEP_SCHEDULE_HH__
#define __SWEEP_SCHEDULE_HH__

#include "Global.hh"
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


#endif

