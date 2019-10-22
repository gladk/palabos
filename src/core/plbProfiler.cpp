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

#include "core/plbProfiler.h"
#include "parallelism/mpiManager.h"
#include "core/runTimeDiagnostics.h"
#include "algorithm/statistics.h"
#include "libraryInterfaces/TINYXML_xmlIO.hh"

namespace plb {

namespace global {

Profiler::Profiler() {
    turnOff();
    automaticCycling();
    setReportFile("plbProfile");

    validCounters.insert("collStreamCells");
    validCounters.insert("iterations");
    validCounters.insert("mpiSendChar");
    validCounters.insert("mpiReceiveChar");
    
    validTimers.insert("collStream");
    validTimers.insert("cycle");
    validTimers.insert("dataProcessor");
    validTimers.insert("mpiCommunication");
    validTimers.insert("io");
    validTimers.insert("totalTime");
}

void Profiler::turnOn() {
    profilingFlag = true;
    start("totalTime");
}

void Profiler::turnOff() {
    profilingFlag = false;
}

void Profiler::automaticCycling() {
    manualCycleFlag = false;
}

void Profiler::manualCycling() {
    manualCycleFlag = true;
}

void Profiler::cycle() {
    increment("iterations");
}


void Profiler::writeReport() {
    plint collStreamCells = getCounter("collStreamCells");
    plint iterations = getCounter("iterations");
    //plint mpiSendChar = getCounter("mpiSendChar");
    //plint mpiReceiveChar = getCounter("mpiReceiveChar");

    double t_collStream = getTimer("collStream");
    double t_cycle = getTimer("cycle");
    //double t_dataProcessor = getTimer("dataProcessor");
    double t_mpiCommunication = getTimer("mpiCommunication");
    double t_io = getTimer("io");
    //double t_totalTime = getTimer("totalTime");

    XMLwriter writer;
    XMLwriter& globalSection(writer["Global"]);

    addMainProcValue(globalSection, "NumProcesses", global::mpi().getSize());
    addMainProcValue(globalSection, "NumIterations", iterations);
    addStatisticalValue(globalSection, "TotalCells", collStreamCells);
    addStatisticalValue(globalSection, "SuS_raw", (double)collStreamCells / t_collStream);
    addStatisticalValue(globalSection, "SuS_overall", (double)collStreamCells / t_cycle);
    addStatisticalValue(globalSection, "Time_per_cycle", t_cycle / (double)iterations);
    addStatisticalValue(globalSection, "Relative_communication_time", t_mpiCommunication / t_cycle);
    addStatisticalValue(globalSection, "Total_io_time", t_io);
    addStatisticalValue(globalSection, "Relative_io_time", t_io / (t_cycle+t_io));

    writer.print(reportFile);
}

void Profiler::addStatisticalValue(XMLwriter& writer, std::string name, double value) {
    std::vector<double> allValues(global::mpi().getSize());
    std::fill(allValues.begin(), allValues.end(), 0.);
    allValues[global::mpi().getRank()] = value;
#ifdef PLB_MPI_PARALLEL
    global::mpi().allReduceVect<double>(allValues, MPI_SUM);
#endif
    util::Stats stats(allValues);
    writer[name]["Mean"].set(stats.getMean());
    writer[name]["StdDev"].set(stats.getStddev());
    writer[name]["Min"].set(stats.getMin());
    writer[name]["Max"].set(stats.getMax());
    writer[name]["Values"].set(allValues);
}

void Profiler::addMainProcValue(XMLwriter& writer, std::string name, plint value) {
    writer[name].set(value);
}


void Profiler::verifyTimer(std::string const& timer) {
    if (validTimers.find(timer)==validTimers.end()) {
        plbLogicError("Invalid timer for profiling: "+timer);
    }
}

void Profiler::verifyCounter(std::string const& counter) {
    if (validCounters.find(counter)==validCounters.end()) {
        plbLogicError("Invalid counter for profiling: "+counter);
    }
}

void Profiler::setReportFile(FileName const& reportFile_) {
    reportFile = reportFile_;
    reportFile.defaultPath(directories().getOutputDir());
    reportFile.defaultExt("xml");
}

}  // namespace global

}  // namespace plb
