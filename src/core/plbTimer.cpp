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

#include "parallelism/mpiManager.h"
#include "core/plbTimer.h"
#include <map>

#include <ctime>

#ifdef PLB_USE_POSIX
#include <unistd.h>
#endif

namespace plb {

namespace global {

/* ************** Timer ***************************************** */

PlbTimer::PlbTimer()
    : cumulativeTime(0.),
      isOn(false)
{ }

void PlbTimer::start() {
#ifdef PLB_MPI_PARALLEL
    startTime = mpi().getTime();
#else
#if defined PLB_USE_POSIX && defined _POSIX_TIMERS && (_POSIX_TIMERS > 0) && !defined(PLB_NGETTIME)
    timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    startTime = (double) ts.tv_nsec * (double) 1.0e-9;
#else
    startClock = clock();
#endif
#endif
    isOn = true;
}

void PlbTimer::restart() {
    reset();
    start();
}

double PlbTimer::stop() {
    cumulativeTime = getTime();
    isOn = false;
    return cumulativeTime;
}

void PlbTimer::reset() {
    cumulativeTime = 0.;
}

double PlbTimer::getTime() const {
    if (isOn) {
#ifdef PLB_MPI_PARALLEL
        return cumulativeTime + mpi().getTime()-startTime;
#else
#if defined PLB_USE_POSIX && defined _POSIX_TIMERS && (_POSIX_TIMERS > 0) && !defined(PLB_NGETTIME)
        timespec ts;
        clock_gettime(CLOCK_REALTIME, &ts);
        double endTime = (double) ts.tv_nsec * (double) 1.0e-9;
        return cumulativeTime + endTime-startTime;
#else
        return cumulativeTime + (double)(clock()-startClock)
                              / (double)CLOCKS_PER_SEC;
#endif
#endif
    }
    else {
        return cumulativeTime;
    }
}

PlbTimer& timer(std::string nameOfTimer) {
    static std::map<std::string, PlbTimer> timerCollection;
    return timerCollection[nameOfTimer];
}

PlbTimer& plbTimer(std::string nameOfTimer) {
    static std::map<std::string, PlbTimer> timerCollection;
    PlbTimer& answer=timerCollection[nameOfTimer];
    return answer;
}

/* ************** Counter ***************************************** */

PlbCounter::PlbCounter()
    : count(0)
{ }

plint PlbCounter::increment(plint value) {
    count += value;
    return count;
}

void PlbCounter::reset() {
    count = 0;
}

plint PlbCounter::getCount() const {
    return count;
}

PlbCounter& counter(std::string nameOfCounter) {
    static std::map<std::string, PlbCounter> counterCollection;
    return counterCollection[nameOfCounter];
}

PlbCounter& plbCounter(std::string nameOfCounter) {
    static std::map<std::string, PlbCounter> counterCollection;
    return counterCollection[nameOfCounter];
}

}  // namespace global

}  // namespace plb
