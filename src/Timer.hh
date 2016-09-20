//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/Timer.hh
 * \author Thomas M. Evans
 * \date   Mon Mar 25 17:35:07 2002
 * \brief  Timer class.
 */
//---------------------------------------------------------------------------//
// $Id: Timer.hh,v 1.5 2003/07/23 20:44:12 kellyt Exp $
//---------------------------------------------------------------------------//

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

#ifndef __c4_Timer_hh__
#define __c4_Timer_hh__

#include "Global.hh"
#include "Assert.hh"
#include <chrono>



class Timer 
{
private:
    std::chrono::high_resolution_clock::time_point beginHighRes;
    std::chrono::high_resolution_clock::time_point endHighRes;
    bool timer_on;
    double sum_wall;
    UINT num_intervals;


public:
    //! Constructor.
    Timer()
    {
        timer_on = false;
        reset();
    }

    //! Set the beginning of the time interval.
    inline
    void start()
    {
        Assert(! timer_on);
        timer_on = true;
        ++num_intervals;
        beginHighRes = std::chrono::high_resolution_clock::now();
    }

    //! Set the end of the time interval.
    inline
    void stop()
    {
        Assert(timer_on);
        endHighRes = std::chrono::high_resolution_clock::now();
        timer_on = false;
        sum_wall += wall_clock();
    }

    //! Return the wall clock time in seconds, for the last interval.
    inline
    double wall_clock() const
    {
        Assert(! timer_on);
        std::chrono::duration<double> timeSpan = 
            std::chrono::duration_cast<std::chrono::duration<double>>(endHighRes - beginHighRes);
        return timeSpan.count();
    }

    //! Return the wall clock time in seconds, summed over all intervals.
    inline
    double sum_wall_clock() const 
    {
        Assert(! timer_on);
        return sum_wall;
    }

    //! Return the number of time intervals used in the sums.
    inline
    int intervals() const
    {
        Assert(! timer_on);
        return num_intervals;
    }

    //! Reset the interval sums.
    inline
    void reset()
    {
        Assert(! timer_on);
        num_intervals = 0;
        sum_wall      = 0.0;
    }
};



#endif                          // __c4_Timer_hh__

//---------------------------------------------------------------------------//
//                              end of c4/Timer.hh
//---------------------------------------------------------------------------//
