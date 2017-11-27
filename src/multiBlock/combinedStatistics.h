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
 * The CombinedStatistics class -- header file.
 */
#ifndef COMBINED_STATISTICS_H
#define COMBINED_STATISTICS_H

#include "core/globalDefs.h"
#include "core/blockStatistics.h"
#include <vector>

namespace plb {

class CombinedStatistics {
public:
    virtual ~CombinedStatistics();
    virtual CombinedStatistics* clone() const =0;
    void combine (
            std::vector<BlockStatistics const*>& individualStatistics,
            BlockStatistics& result ) const;
protected:
    virtual void reduceStatistics (
            std::vector<double>& averageObservables,
            std::vector<double>& sumWeights,
            std::vector<double>& sumObservables,
            std::vector<double>& maxObservables,
            std::vector<plint>& intSumObservables ) const =0;
private:
    void computeLocalAverage (
            std::vector<BlockStatistics const*> const& individualStatistics,
            std::vector<double>& averageObservables,
            std::vector<double>& sumWeights ) const;
    void computeLocalSum (
            std::vector<BlockStatistics const*> const& individualStatistics,
            std::vector<double>& sumObservables ) const;
    void computeLocalMax (
            std::vector<BlockStatistics const*> const& individualStatistics,
            std::vector<double>& maxObservables ) const;
    void computeLocalIntSum (
            std::vector<BlockStatistics const*> const& individualStatistics,
            std::vector<plint>& intSumObservables ) const;
};

class SerialCombinedStatistics : public CombinedStatistics {
public:
    virtual SerialCombinedStatistics* clone() const;
protected:
    virtual void reduceStatistics (
            std::vector<double>& averageObservables,
            std::vector<double>& sumWeights,
            std::vector<double>& sumObservables,
            std::vector<double>& maxObservables,
            std::vector<plint>& intSumObservables ) const;
};

}  // namespace plb

#endif
