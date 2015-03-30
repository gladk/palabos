/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2015 FlowKit Sarl
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

/* Main author: Daniel Lagrava
 */

#include "multiGrid/parallelizer3D.h"
#include "parallelism/mpiManager.h"
#include "multiGrid/multiScale.h"
#include "io/parallelIO.h"

namespace plb {
    
void Parallelizer3D::parallelizeLevel(plint whichLevel,
                                      std::vector<std::vector<Box3D> > const& originalBlocks,
                                      std::vector<Box3D> const& parallelRegions,
                                      std::vector<plint> const& regionIDs )
{
    PLB_PRECONDITION( parallelRegions.size() == regionIDs.size() );
    PLB_PRECONDITION( whichLevel < (plint)originalBlocks.size() );
    std::vector<Box3D> newBlocks;
    // IDs are going to be reattributed at the level whichLevel.
    if (finalMpiDistribution.size() <= (pluint)whichLevel) {
        finalMpiDistribution.resize(whichLevel+1);
    }
    for (pluint iRegion=0; iRegion<parallelRegions.size(); ++iRegion) {
        plint currentId = regionIDs[iRegion];
        for (pluint iBlock=0; iBlock<originalBlocks[whichLevel].size(); ++iBlock) {
            Box3D intersection;
            if ( intersect( originalBlocks[whichLevel][iBlock],
                            parallelRegions[iRegion], intersection ) )
            {
                newBlocks.push_back(intersection);
                finalMpiDistribution[whichLevel].push_back(currentId);
            }
        }
    }

    recomputedBlocks[whichLevel].insert( recomputedBlocks[whichLevel].end(),newBlocks.begin(),newBlocks.end() );
}

plint Parallelizer3D::computeCost(std::vector<std::vector<Box3D> > const& originalBlocks, Box3D box){
    plint totalCost = 0;
    plint numLevels = originalBlocks.size();
    
    for (plint iLevel=(plint)originalBlocks.size()-1; iLevel>=0; --iLevel){
        // convert the box to the current level
        Box3D levelBox = global::getDefaultMultiScaleManager().scaleBox(box,iLevel-(numLevels-1));
        for (pluint iComp=0; iComp<originalBlocks[iLevel].size(); ++iComp){
            Box3D currentBox;
            if (intersect(originalBlocks[iLevel][iComp], levelBox, currentBox)){
                plint volume = currentBox.getNx()*currentBox.getNy()*currentBox.getNz();
                totalCost += (plint) util::twoToThePower(iLevel) * volume;
            }
        }
    }
    
    return totalCost;
}

/* ************* ParallellizeByCubes3D **************** */
ParallellizeByCubes3D::ParallellizeByCubes3D(std::vector<std::vector<Box3D> > const& originalBlocks_,
                                     Box3D finestBoundingBox_, plint xTiles_, plint yTiles_, plint zTiles_ )
    : originalBlocks(originalBlocks_), finestBoundingBox(finestBoundingBox_),
      processorNumber(global::mpi().getSize()),xTiles(xTiles_),yTiles(yTiles_),zTiles(zTiles_) 
{
    // divide the finest bounding box in xTiles by yTiles by zTiles cubes
    computeFinestDivision();
}

void ParallellizeByCubes3D::computeFinestDivision(){
    std::vector<std::pair<plint,plint> > rangesX;
    std::vector<std::pair<plint,plint> > rangesY;
    std::vector<std::pair<plint,plint> > rangesZ;
    
    plint nx = finestBoundingBox.x1;
    plint ny = finestBoundingBox.y1;
    plint nz = finestBoundingBox.z1;

    util::linearRepartition(0, nx, xTiles, rangesX);
    util::linearRepartition(0, ny, yTiles, rangesY);
    util::linearRepartition(0, nz, zTiles, rangesZ);

    finestDivision.resize(0);
    mpiDistribution.resize(xTiles*yTiles*zTiles);
    
    for (plint iX=0; iX<(plint)rangesX.size(); ++iX){
        for (plint iY=0; iY<(plint)rangesY.size(); ++iY){
            for (plint iZ=0; iZ<(plint)rangesZ.size(); ++iZ){
                // create a Box3D with coordinates for the new sector
                Box3D box(rangesX[iX].first,rangesX[iX].second,rangesY[iY].first,rangesY[iY].second,
                      rangesZ[iZ].first,rangesZ[iZ].second);
                // put it in the finestDivision
                finestDivision.push_back(box);
            }
        }
    }
}

void ParallellizeByCubes3D::parallelize(){
    // we must have the same number of blocks as processors
    PLB_PRECONDITION(xTiles*yTiles*zTiles == processorNumber);
    
    plint totalCost = computeCost(originalBlocks, finestBoundingBox);
    plint idealCostPerProcessor = totalCost/processorNumber;
    
    pcout << "Total cost of computations = " << totalCost << std::endl;
    pcout << "We are using " << processorNumber << " processors...\n";
    pcout << "Ideal cost per processor = " << idealCostPerProcessor << std::endl;
    
    std::vector<plint> totalCosts(processorNumber);
    
    // greedy load balancing part
//     plint currentProcessor = 0;
//     bool allAssigned = false;
//     plint iBlock = 0;
//     plint maxBlockCost = 0;
//     while ( (currentProcessor<processorNumber) && !allAssigned) {
//         plint blockCost = computeCost(originalBlocks, finestDivision[iBlock]);
//         if (blockCost > maxBlockCost) maxBlockCost = blockCost;
//         if ( totalCosts[currentProcessor] < idealCostPerProcessor ){
//             totalCosts[currentProcessor] += blockCost;
//             mpiDistribution[iBlock] = currentProcessor;
//             ++iBlock;
//         }
//         else {
//             currentProcessor++;
//         }
//         if (iBlock>=(plint)finestDivision.size()){
//             allAssigned = true;
//         }
//     } 
//     
//     if (maxBlockCost > idealCostPerProcessor){
//         pcout << "There is a problem : maxBlockCost=" << maxBlockCost << " and ideal cost=" << idealCostPerProcessor
//               << std::endl;
//     }
    
    plint total = 0;
    for (plint iProc=0; iProc<processorNumber; ++iProc){
        plint blockCost = computeCost(originalBlocks,finestDivision[iProc]);
        totalCosts[iProc] += blockCost;
        mpiDistribution[iProc] = iProc;
        total += blockCost;
    }
    
    pcout << "---- Costs Per Processor ----\n";
    
    for (pluint i=0; i<totalCosts.size(); ++i){
        pcout << i << " : " << totalCosts[i] << std::endl;
        // check if everyone is doing something
        if (totalCosts[i] == 0){
            pcout << "\t >> processor " << i << " does not have work to do. Exiting.....\n";
            std::exit(1);
        }
    }
    
    pcout << "*******************************\n";
    pcout << "Sum of all costs = " << total << std::endl;
    pcout << "*******************************\n";
    
    // convert the original blocks to the new blocks
    recomputedBlocks.resize(originalBlocks.size());
    finalMpiDistribution.resize(originalBlocks.size());
    
    plint finestLevel= (plint)originalBlocks.size()-1;
    for (plint iLevel=finestLevel; iLevel>=0; --iLevel) {
        parallelizeLevel(iLevel, originalBlocks,finestDivision, mpiDistribution);
        // Adapt the regions to the next-coarser level.
        for (pluint iRegion=0; iRegion<finestDivision.size(); ++iRegion) {
            finestDivision[iRegion] = finestDivision[iRegion].divideAndFitSmaller(2);
        }
    }
}

} // namespace plb

