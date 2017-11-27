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

#ifndef PLB_PROFILER_H
#define PLB_PROFILER_H

#include "core/globalDefs.h"
#include "core/plbTimer.h"
#include "io/plbFiles.h"
#include "libraryInterfaces/TINYXML_xmlIO.h"
#include <string>
#include <set>

namespace plb {

namespace global {

/**
 * Counters:
 * =========
 * "collStreamCells":                Number of coll-stream cells.
 * "iterations":                     Number of iterations.
 * "mpiSendChar":                    Number of bytes sent by MPI.
 * "mpiReceiveChar":                 Number of bytes received by MPI.
 *
 * Timers:
 * =======
 * "collStream":                     Time for raw execution of collision-streaming.
 * "cycle":                          Time for full cycles, including communication + data processors.
 * "dataProcessor":                  Time for data processor calls.
 * "envelope-update":                Time for update of envelopes, including MPI communication.
 * "mpiCommunication":               Total Time for MPI communication.
 * "io":                             Time spent for I/O operations.
 * "totalTime":                      Total time.
**/
class Profiler {
public:
    void turnOn();
    void turnOff();
    void automaticCycling();
    void manualCycling();
    void cycle();
    bool cyclingIsAutomatic() const {
        return !manualCycleFlag;
    }
    bool doProfiling() const {
        return profilingFlag;
    }
    void start(char const* timer) {
        if (doProfiling()) {
            verifyTimer(timer);
            plbTimer(timer).start();
        }
    }
    void stop(char const* timer) {
        if (doProfiling()) {
            verifyTimer(timer);
            plbTimer(timer).stop();
        }
    }
    void increment(char const* counter) {
        if (doProfiling()) {
            verifyCounter(counter);
            plbCounter(counter).increment();
        }
    }
    void increment(char const* counter, plint value) {
        if (doProfiling()) {
            verifyCounter(counter);
            plbCounter(counter).increment(value);
        }
    }
    plint getCounter(char const* counter) {
        verifyCounter(counter);
        return plbCounter(counter).getCount();
    }
    double getTimer(char const* timer) {
        verifyTimer(timer);
        return plbTimer(timer).getTime();
    }
    void setReportFile(FileName const& reportFile_);
    void writeReport();
private:
    void verifyTimer(std::string const& timer);
    void verifyCounter(std::string const& counter);
    void addStatisticalValue(XMLwriter& writer, std::string name, double value);
    void addMainProcValue(XMLwriter& writer, std::string name, plint value);

    Profiler();
private:
    bool profilingFlag;
    bool manualCycleFlag;
    FileName reportFile;
    std::set<std::string> validTimers;
    std::set<std::string> validCounters;
friend Profiler& profiler();
};

inline Profiler& profiler() {
    static Profiler instance;
    return instance;
}

}  // namespace global

}  // namespace plb

#endif  // PLB_PROFILER_H
