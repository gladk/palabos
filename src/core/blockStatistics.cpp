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
 * The BlockStatistics class -- generic implementation.
 */
#include "core/blockStatistics.h"
#include "core/util.h"
#include "core/plbDebug.h"
#include <cmath>
#include <numeric>
#include <limits>
#include <iostream>

namespace plb {

////////////////////// Class BlockStatistics /////////////////

BlockStatistics::BlockStatistics()
  : tmpNumCells(0)
{ }

BlockStatistics::BlockStatistics(BlockStatistics const& rhs)
  : tmpAv(rhs.tmpAv),
    tmpSum(rhs.tmpSum),
    tmpMax(rhs.tmpMax),
    tmpIntSum(rhs.tmpIntSum),
    tmpNumCells(rhs.tmpNumCells),
    doubleReductions(rhs.doubleReductions),
    averageVect(rhs.averageVect),
    sumVect(rhs.sumVect),
    maxVect(rhs.maxVect),
    intSumVect(rhs.intSumVect),
    numCells(rhs.numCells)
{}

void BlockStatistics::swap(BlockStatistics& rhs)
{
    tmpAv.swap    (rhs.tmpAv);
    tmpSum.swap   (rhs.tmpSum);
    tmpMax.swap   (rhs.tmpMax);
    tmpIntSum.swap(rhs.tmpIntSum);
    std::swap(tmpNumCells, rhs.tmpNumCells);

    doubleReductions.swap(rhs.doubleReductions);
    averageVect.swap(rhs.averageVect);
    sumVect.swap    (rhs.sumVect);
    maxVect.swap    (rhs.maxVect);
    intSumVect.swap (rhs.intSumVect);
    std::swap(numCells, rhs.numCells);
}

/** In the reset function, running statistics are copied to public statistics,
 *  and running statistics are reset to zero.
 */
void BlockStatistics::evaluate() {
    // First step: copy running statistics to public statistics

    // Avoid division by zero while evaluating average: if no cell has
    //   been accounted for so far, simply reset to zero.
    if (tmpNumCells == 0) {   
        for (pluint iVect=0; iVect<averageVect.size(); ++iVect) {
            averageVect[iVect] = 0.;
        }
    }
    else {
        for (pluint iVect=0; iVect<averageVect.size(); ++iVect) {
            averageVect[iVect] = tmpAv[iVect] / (double)tmpNumCells;
        }
    }
    for (pluint iVect=0; iVect<sumVect.size(); ++iVect) {
        sumVect[iVect]     = tmpSum[iVect];
    }
    for (pluint iVect=0; iVect<maxVect.size(); ++iVect) {
        maxVect[iVect]     = tmpMax[iVect];
    }
    for (pluint iVect=0; iVect<intSumVect.size(); ++iVect) {
        intSumVect[iVect]  = tmpIntSum[iVect];
    }
    numCells = tmpNumCells;

    // Second step: reset the running statistics, in order to be ready
    //   for next lattice iteration
    for (pluint iVect=0; iVect<averageVect.size(); ++iVect) {
        tmpAv[iVect]     = 0.;
    }
    for (pluint iVect=0; iVect<sumVect.size(); ++iVect) {
        tmpSum[iVect]    = 0.;
    }
    for (pluint iVect=0; iVect<maxVect.size(); ++iVect) {
        // Use -max() instead of min(), because min<float> yields a positive value close to zero.
        tmpMax[iVect]    = -std::numeric_limits<double>::max();
    }
    for (pluint iVect=0; iVect<intSumVect.size(); ++iVect) {
        tmpIntSum[iVect] = 0;
    }

    tmpNumCells = 0;
}

void BlockStatistics::evaluate (
        std::vector<double> const& average, std::vector<double> const& sum,
        std::vector<double> const& max, std::vector<plint> const& intSum, pluint numCells_ )
{
    PLB_PRECONDITION ( averageVect.size() == average.size() );
    PLB_PRECONDITION ( sumVect.size()     == sum.size() );
    PLB_PRECONDITION ( maxVect.size()     == max.size() );
    PLB_PRECONDITION ( intSumVect.size()  == intSum.size() );

    for (pluint iAverage=0; iAverage<averageVect.size(); ++iAverage) {
        averageVect[iAverage] = average[iAverage];
        tmpAv[iAverage] = 0.;
    }
    for (pluint iSum=0; iSum<sumVect.size(); ++iSum) {
        sumVect[iSum] = sum[iSum];
        tmpSum[iSum]  = 0.;
    }
    for (pluint iMax=0; iMax<maxVect.size(); ++iMax) {
        maxVect[iMax] = max[iMax];
        // Use -max() instead of min(), because min<float> yields a positive value close to zero.
        tmpMax[iMax]  = -std::numeric_limits<double>::max();
    }
    for (pluint iSum=0; iSum<intSumVect.size(); ++iSum) {
        intSumVect[iSum] = intSum[iSum];
        tmpIntSum[iSum]  = 0;
    }
    numCells    = numCells_;
    tmpNumCells = 0;
}

/** \return Identifier for this observable, to be used for gatherAverage() and getAverage().
 */
plint BlockStatistics::subscribeAverage() {
    doubleReductions.push_back(averageRed);
    plint newSize = tmpAv.size()+1;
    tmpAv.resize(newSize);
    averageVect.resize(newSize);
    plint newIndex = newSize-1;
    tmpAv[newIndex] = 0.;
    averageVect[newIndex] = 0.;
    return newIndex;
}

/** \return Identifier for this observable, to be used for gatherSum() and getSum().
 */
plint BlockStatistics::subscribeSum() {
    doubleReductions.push_back(sumRed);
    plint newSize = tmpSum.size()+1;
    tmpSum.resize(newSize);
    sumVect.resize(newSize);
    plint newIndex = newSize-1;
    tmpSum[newIndex] = 0.;
    sumVect[newIndex] = 0.;
    return newIndex;
}

/** \return Identifier for this observable, to be used for gatherMax() and getMax().
 */
plint BlockStatistics::subscribeMax() {
    doubleReductions.push_back(maxRed);
    plint newSize = tmpMax.size()+1;
    tmpMax.resize(newSize);
    maxVect.resize(newSize);
    plint newIndex = newSize-1;
    // Use -max() instead of min(), because min<float> yields a positive value close to zero.
    tmpMax[newIndex] = -std::numeric_limits<double>::max();
    // Use -max() instead of min(), because min<float> yields a positive value close to zero.
    maxVect[newIndex] = -std::numeric_limits<double>::max();
    return newIndex;
}

/** \return Identifier for this observable, to be used for gatherIntSum() and getIntSum().
 */
plint BlockStatistics::subscribeIntSum() {
    plint newSize = tmpIntSum.size()+1;
    tmpIntSum.resize(newSize);
    intSumVect.resize(newSize);
    plint newIndex = newSize-1;
    tmpIntSum[newIndex] = 0;
    intSumVect[newIndex] = 0;
    return newIndex;
}

void BlockStatistics::gatherAverage(plint whichAverage, double value) {
    PLB_PRECONDITION( whichAverage < (plint) tmpAv.size() );
    tmpAv[whichAverage] += value;
}

void BlockStatistics::gatherSum(plint whichSum, double value) {
    PLB_PRECONDITION( whichSum < (plint) tmpSum.size() );
    tmpSum[whichSum] += value;
}

void BlockStatistics::gatherMax(plint whichMax, double value) {
    PLB_PRECONDITION( whichMax < (plint) tmpMax.size() );
    if (value > tmpMax[whichMax]) {
        tmpMax[whichMax] = value;
    }
}

void BlockStatistics::gatherIntSum(plint whichSum, plint value) {
    PLB_PRECONDITION( whichSum < (plint) tmpIntSum.size() );
    tmpIntSum[whichSum] += value;
}

void BlockStatistics::incrementStats() {
    ++tmpNumCells;
}

double BlockStatistics::getAverage(plint whichAverage) const {
    PLB_PRECONDITION( whichAverage < (plint) tmpAv.size() );
    return averageVect[whichAverage];
}

double BlockStatistics::getSum(plint whichSum) const {
    PLB_PRECONDITION( whichSum < (plint) tmpSum.size() );
    return sumVect[whichSum];
}

double BlockStatistics::getMax(plint whichMax) const {
    PLB_PRECONDITION( whichMax < (plint) tmpMax.size() );
    return maxVect[whichMax];
}

plint BlockStatistics::getIntSum(plint whichSum) const {
    PLB_PRECONDITION( whichSum < (plint) tmpIntSum.size() );
    return intSumVect[whichSum];
}


void BlockStatistics::rescale(std::vector<double> const& scales) {
    PLB_PRECONDITION( scales.size() == doubleReductions.size() );
    pluint iAverage=0;
    pluint iSum=0;
    pluint iMax=0;
    for(pluint iRed=0; iRed<doubleReductions.size(); ++iRed) {
        switch(doubleReductions[iRed]) {
            case averageRed: averageVect[iAverage++]*=scales[iRed]; break;
            case sumRed:     sumVect[iSum++]*=scales[iRed]; break;
            case maxRed:     maxVect[iMax++]*=scales[iRed]; break;
        }
    }
}

void combine(std::vector<BlockStatistics*>& components, BlockStatistics& result)
{
    PLB_PRECONDITION( components.size()>0 );
    std::vector<double> averageVect(components[0]->getAverageVect());
    std::vector<double> sumVect(components[0]->getSumVect());
    std::vector<double> maxVect(components[0]->getMaxVect());
    std::vector<plint>  intSumVect(components[0]->getIntSumVect());

    for (pluint iComp=1; iComp<components.size(); ++iComp) {
        // TODO: Here, averages are treated like sums. This needs to be
        //   corrected.
        for (pluint iAve=0; iAve<averageVect.size(); ++iAve) {
            averageVect[iAve] += components[iComp]->getAverageVect()[iAve];
        }
        for (pluint iSum=0; iSum<sumVect.size(); ++iSum) {
            sumVect[iSum] += components[iComp]->getSumVect()[iSum];
        }
        for (pluint iMax=0; iMax<maxVect.size(); ++iMax) {
            maxVect[iMax] = std::max(maxVect[iMax],components[iComp]->getMaxVect()[iMax]);
        }
        for (pluint iIntSum=0; iIntSum<intSumVect.size(); ++iIntSum) {
            intSumVect[iIntSum] += components[iComp]->getIntSumVect()[iIntSum];
        }
    }


    result.evaluate(averageVect, sumVect, maxVect, intSumVect, 0);
}

}  // namespace plb
