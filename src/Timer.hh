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

#ifndef __TIMER_HH__
#define __TIMER_HH__

#include "Global.hh"
#include "Assert.hh"
#include <chrono>


/*
    Timer

    Implements a high resolution timer.
*/
class Timer 
{
private:
    std::chrono::high_resolution_clock::time_point begin;
    std::chrono::high_resolution_clock::time_point end;
    bool timer_on;
    double sum_wall;
    UINT num_intervals;


public:
    // Constructor
    Timer()
    {
        timer_on = false;
        reset();
    }

    // Start the timer
    inline
    void start()
    {
        Assert(!timer_on);
        timer_on = true;
        ++num_intervals;
        begin = std::chrono::high_resolution_clock::now();
    }

    // Stop the timer (does not reset the timer)
    inline
    void stop()
    {
        Assert(timer_on);
        end = std::chrono::high_resolution_clock::now();
        timer_on = false;
        sum_wall += wall_clock();
    }

    // Get last timer stop() - start()
    inline
    double wall_clock() const
    {
        Assert(!timer_on);
        std::chrono::duration<double> timeSpan = end - begin;
        return timeSpan.count();
    }

    // Return sum of stop() - start() times
    inline
    double sum_wall_clock() const 
    {
        Assert(!timer_on);
        return sum_wall;
    }

    // Return number of times timer has been used
    inline
    int intervals() const
    {
        Assert(!timer_on);
        return num_intervals;
    }

    // Reset timer
    inline
    void reset()
    {
        Assert(!timer_on);
        num_intervals = 0;
        sum_wall      = 0.0;
    }
};



#endif

