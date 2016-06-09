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

#ifndef __c4_Timer_hh__
#define __c4_Timer_hh__

#include "Typedef.hh"
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
