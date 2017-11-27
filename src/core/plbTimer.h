/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2017 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at 
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/** \file
 * A timer class for benchmarking program parts -- header file.
 */
#ifndef PLB_TIMER_H
#define PLB_TIMER_H

#include <ctime>

#ifdef PLB_USE_POSIX
#include <unistd.h>
#endif

namespace plb {

namespace global {

/// A cumulative timer for benchmarking program parts and summing up the times.
/** In serial, the C function clock() is used. In parallel, the MPI
 *  function is used.
 */
class PlbTimer {
public:
    PlbTimer();
    /// Proceed with time measurement.
    void start();
    /// Reset timer to zero and start time measurement.
    void restart();
    /// Interrupt time measurement ( you can still proceed with start() ).
    /* \return Current cumulative time.
     */
    double stop();
    /// Reset timer to zero.
    void reset();
    /// Get current cumulative time.
    double getTime() const;
private:
    double cumulativeTime;
    bool   isOn;
#ifdef PLB_MPI_PARALLEL
    double startTime;
#else
#if defined PLB_USE_POSIX && defined _POSIX_TIMERS && (_POSIX_TIMERS > 0) && !defined(PLB_NGETTIME)
    double startTime;
#else
    clock_t startClock;
#endif
#endif
friend PlbTimer& timer(std::string nameOfTimer);
friend PlbTimer& plbTimer(std::string nameOfTimer);
};

// Global instance of timer objects, for public use.
PlbTimer& timer(std::string nameOfTimer);

// Global instance of timer objects, for internal use.
PlbTimer& plbTimer(std::string nameOfTimer);

/// A cumulative counter for benchmarking program parts and summing up occurrences
///   of events.
class PlbCounter {
public:
    PlbCounter();
    /// Increment the counter.
    plint increment(plint value=1);
    /// Reset counter to zero.
    void reset();
    /// Return current count.
    plint getCount() const;
private:
    plint count;
friend PlbCounter& counter(std::string nameOfCounter);
friend PlbCounter& plbCounter(std::string nameOfCounter);
};

// Global instance of counter objects, for public use.
PlbCounter& counter(std::string nameOfCounter);

// Global instance of counter objects, for internal use.
PlbCounter& plbCounter(std::string nameOfCounter);

}  // namespace global

}  // namespace plb

#endif  // PLB_TIMER_H
