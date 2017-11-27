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
 * The BlockStatistics class -- header file.
 */
#ifndef BLOCK_STATISTICS_H
#define BLOCK_STATISTICS_H

#include "core/globalDefs.h"
#include <vector>
#include <algorithm>

namespace plb {

// Forward declaration

class BlockStatistics;

/// A counter for keeping track of discrete time evolution
class TimeCounter {
public:
    TimeCounter() : latticeTime(0) { }
    TimeCounter(plint iniTime) : latticeTime(iniTime) { }
    /// Increment the value of time step in this lattice
    void incrementTime() { ++latticeTime; };
    /// Reset the value of time step in this lattice
    void resetTime(plint value=0) { latticeTime=value; }
    /// Get the value of time step in this lattice
    plint getTime() const { return latticeTime; }
private:
    plint latticeTime;
};

/// A polymorphic type to handle subscription to a BlockStatistics class
class StatSubscriber {
public:
    virtual ~StatSubscriber() { }
    /// Subscribe a new observable for which the average value is computed.
    virtual plint subscribeAverage() =0;
    /// Subscribe a new observable for which the sum is computed.
    virtual plint subscribeSum() =0;
    /// Subscribe a new observable for which the maximum is computed.
    virtual plint subscribeMax() =0;
    /// Subscribe a new integer observable for which the sum is computed.
    virtual plint subscribeIntSum() =0;
};

/// Store instances of observables, and compute their statistics.
/** This class is not intended to be inherited from. It is not a template, because 
 *  statistics should always be in double-precision, even if the simulation
 *  is in single-precision. Remember that the statistics makes reductions over
 *  large data sets, in which case single-precision arithmetics can rapidly be
 *  insufficient.
 */
class BlockStatistics {
public:
    BlockStatistics();
    BlockStatistics(BlockStatistics const& rhs);
    
    void swap(BlockStatistics& rhs);
    /// evaluate() must be called after each lattice iteration.
    void evaluate();
    /// Attribute a value to the public statistics, and reset running statistics to default.
    void evaluate(std::vector<double> const& average, std::vector<double> const& sum,
                  std::vector<double> const& max, std::vector<plint> const& intSum, pluint numCells_);
    /// Contribute the values of the current cell to the statistics of an "average observable"
    void gatherAverage(plint whichAverage, double value);
    /// Contribute the values of the current cell to the statistics of a "sum observable"
    void gatherSum(plint whichSum, double value);
    /// Contribute the values of the current cell to the statistics of a "max observable"
    void gatherMax(plint whichMax, double value);
    /// Contribute the values of the current cell to the statistics of an integer "sum observable"
    void gatherIntSum(plint whichSum, plint value);
    /// Call this function once all statistics for a cell have been added
    void incrementStats();
    /// Return number of cells for which statistics have been added so far
    pluint const& getNumCells() const { return numCells; }

    /// Get the public value for any "average observable"
    double getAverage(plint whichAverage) const;
    /// Get the public value for any "sum observable"
    double getSum(plint whichSum) const;
    /// Get the public value for any "max observable"
    double getMax(plint whichMax) const;
    /// Get the public value for any integer "sum observable"
    plint getIntSum(plint whichSum) const;

    /// Get a handle to the vector with all "average observables"
    std::vector<double>& getAverageVect() { return averageVect; }
    /// Get a handle to the vector with all "sum observables"
    std::vector<double>& getSumVect() { return sumVect; }
    /// Get a handle to the vector with all "max observables"
    std::vector<double>& getMaxVect() { return maxVect; }
    /// Get a handle to the vector with all integer "sum observables"
    std::vector<plint>& getIntSumVect() { return intSumVect; }
    /// Get all real-valued statistics, in their order of subscription.
    std::vector<double> getDoubleVect() const;

    /// Subscribe a new observable for which the average value is computed.
    plint subscribeAverage();
    /// Subscribe a new observable for which the sum is computed.
    plint subscribeSum();
    /// Subscribe a new observable for which the maximum is computed.
    plint subscribeMax();
    /// Subscribe a new integer observable for which the sum is computed.
    plint subscribeIntSum();
    /// Rescale the floating-point values computed so far.
    void rescale(std::vector<double> const& scales);
private:
    /// Give a name to floating point reductions to keep track of their order of subscription.
    enum DoubleReductions {averageRed, sumRed, maxRed};
    /// Variables to store running statistics of type double.
    std::vector<double> tmpAv, tmpSum, tmpMax;
    /// Variables to store summed integer observables
    std::vector<plint> tmpIntSum;
    /// Running value for number of cells over which statistics has been computed
    pluint tmpNumCells;
    /// Keep track of the order of subscriptions for floating point reductions.
    std::vector<DoubleReductions> doubleReductions;
    /// Variables containing the public result of type double
    std::vector<double> averageVect, sumVect, maxVect;
    /// Variables containing the public result for the summed integer observables
    std::vector<plint> intSumVect;
    /// Public result for number of cells over which statistics has been computed
    pluint numCells;
};

void combine(std::vector<BlockStatistics*>& components, BlockStatistics& result);

}  // namespace plb

#endif

