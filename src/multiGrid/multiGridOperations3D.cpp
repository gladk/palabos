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

/* Main author: Daniel Lagrava, adapted by Helen Morrison
 **/

/** \file
 * Operations on the 3D multigrid -- implementation
 */

#include "core/globalDefs.h"
#include "multiGrid/multiGridOperations3D.h"
#include "multiBlock/multiBlockOperations3D.h"
#include "atomicBlock/dataProcessor3D.h"
#include "multiGrid/multiGridUtil.h"
#include "multiGrid/multiScale.h"

#include "multiBlock/multiBlockInfo3D.h"
#include "io/parallelIO.h"

#include <vector>
#include <memory>

namespace plb {

/// Execute a data processor over several MultiGrid3D
void executeDataProcessor( DataProcessorGenerator3D const& generator,
                           std::vector<MultiGrid3D*> multiGrids,
                           plint referenceLevel )
{
    if(multiGrids.empty()) return;
    
    for (plint iLevel=0; iLevel<(plint)multiGrids[0]->getNumLevels(); ++iLevel) {
        std::auto_ptr<DataProcessorGenerator3D> localGenerator(generator.clone());
        plint dxScale = referenceLevel - iLevel;
        plint dtScale = dxScale;  // TODO: here, we assume convective scaling; general case should be considered.
        
        int boxRescaleFactor = util::roundToInt(util::twoToThePower(std::abs(referenceLevel-iLevel)));
        if (dxScale > 0) // if we go to a coarser grid          // I (Helen, 2014) changed the sign here... this seems to work for my particular case, but I'm not entirely sure about this!!
            localGenerator->divide(boxRescaleFactor);  
        else  // otherwise we go to a finer grid
            localGenerator->multiply(boxRescaleFactor);
        
        localGenerator->setscale(dxScale,dtScale);
        
        std::vector<MultiBlock3D*> localBlocks(multiGrids.size());
        for (plint iBlock=0; iBlock<(plint)localBlocks.size(); ++iBlock) {
            localBlocks[iBlock] =  &multiGrids[iBlock]->getComponent(iLevel);
        }
        executeDataProcessor(*localGenerator, localBlocks);
    }
}


/// Execute a data processor over a single multiGrid3D
void executeDataProcessor( DataProcessorGenerator3D const& generator,
                           MultiGrid3D& object, plint referenceLevel )
{
    std::vector<MultiGrid3D*> blocks;
    blocks.push_back(&object);
    executeDataProcessor(generator,blocks,referenceLevel);
}

/// Execute a data processor over two multiGrid3D
void executeDataProcessor( DataProcessorGenerator3D const& generator,
                           MultiGrid3D& object1, MultiGrid3D& object2,
                           plint referenceLevel )
{
    std::vector<MultiGrid3D*> blocks;
    blocks.push_back(&object1);
    blocks.push_back(&object2);
    executeDataProcessor(generator,blocks,referenceLevel);
}


void executeDataProcessor( ReductiveDataProcessorGenerator3D& generator,
                           std::vector<MultiGrid3D*> multiGrids,
                           plint referenceLevel )
{
    if(multiGrids.empty()) return;
    plint numLevels = (plint)multiGrids[0]->getNumLevels();
    std::vector<int> dimensionsX, dimensionsT;
    generator.getDimensionsX(dimensionsX);
    generator.getDimensionsT(dimensionsT);
    std::vector<ReductiveDataProcessorGenerator3D*> localGenerators(numLevels);
    std::vector<BlockStatistics*> localStatistics(numLevels);

    for (plint iLevel=0; iLevel<numLevels; ++iLevel) {
        int dxScale = referenceLevel - iLevel;
        int dtScale = dxScale;  // TODO: here, we assume convective scaling; general case could be considered.
        localGenerators[iLevel] = generator.clone();
        
        int boxRescaleFactor = util::roundToInt(util::twoToThePower(std::abs(referenceLevel-iLevel)));
        if (dxScale > 0)
            generator.divide(boxRescaleFactor);
        else
            generator.multiply(boxRescaleFactor);

        localGenerators[iLevel]->setscale(dxScale,dtScale);

        std::vector<MultiBlock3D*> localBlocks(multiGrids.size());
        for (plint iBlock=0; iBlock<(plint)localBlocks.size(); ++iBlock) {
            localBlocks[iBlock] = &multiGrids[iBlock]->getComponent(iLevel);
        }
        executeDataProcessor(*localGenerators[iLevel], localBlocks);
//        std::vector<double> scales(dimensionsX.size());
//        for (pluint iScale=0; iScale<scales.size(); ++iScale) {
//            scales[iScale] = scaleToReference(dxScale, dimensionsX[iScale], dtScale, dimensionsT[iScale]);
//        }
//        localGenerators[iLevel]->getStatistics().rescale(scales);
        localStatistics[iLevel] = &(localGenerators[iLevel]->getStatistics());
    }

    combine(localStatistics, generator.getStatistics());
}

void executeDataProcessor( ReductiveDataProcessorGenerator3D& generator,
                           MultiGrid3D& object,
                           plint referenceLevel )
{
    std::vector<MultiGrid3D*> blocks;
    blocks.push_back(&object);
    executeDataProcessor(generator,blocks,referenceLevel);
}

void executeDataProcessor( ReductiveDataProcessorGenerator3D& generator,
                           MultiGrid3D& object1, MultiGrid3D& object2,
                           plint referenceLevel )
{
    std::vector<MultiGrid3D*> blocks;
    blocks.push_back(&object1);
    blocks.push_back(&object2);
    executeDataProcessor(generator,blocks,referenceLevel);
}


void addInternalProcessor( DataProcessorGenerator3D const& generator,
                           std::vector<MultiGrid3D*> multiGrids,
                           plint referenceLevel,
                           plint processorLevel)
{
    if(multiGrids.empty()) return;
    for (plint iLevel=0; iLevel<(plint)multiGrids[0]->getNumLevels(); ++iLevel) {
        std::auto_ptr<DataProcessorGenerator3D> localGenerator(generator.clone());
        plint dxScale = referenceLevel - iLevel;
        plint dtScale = dxScale;  // TODO: here, we assume convective scaling; general case should be considered.
        
        localGenerator -> setscale(dxScale,dtScale);
        int boxRescaleFactor = util::roundToInt(util::twoToThePower(std::abs(referenceLevel-iLevel)));
        // if dxScale < 0 this means we go to a coarser grid
        if (dxScale < 0)
            localGenerator->divide(boxRescaleFactor);  
        else
            localGenerator->multiply(boxRescaleFactor);
        
        std::vector<MultiBlock3D*> localBlocks(multiGrids.size());
        for (plint iBlock=0; iBlock<(plint)localBlocks.size(); ++iBlock) {
            localBlocks[iBlock] =  &multiGrids[iBlock]->getComponent(iLevel);
        }
        addInternalProcessor(*localGenerator, localBlocks, processorLevel);
    }
}

void addInternalProcessor( DataProcessorGenerator3D const& generator,
                           MultiGrid3D& object,
                           plint referenceLevel,
                           plint processorLevel)
{
    std::vector<MultiGrid3D*> blocks;
    blocks.push_back(&object);
    addInternalProcessor(generator,blocks,referenceLevel, processorLevel);
}

void addInternalProcessor( DataProcessorGenerator3D const& generator,
                           MultiGrid3D& object1, MultiGrid3D& object2,
                           plint referenceLevel,
                           plint processorLevel)
{
    std::vector<MultiGrid3D*> blocks;
    blocks.push_back(&object1);
    blocks.push_back(&object2);
    addInternalProcessor(generator,blocks,referenceLevel,processorLevel);
}

} // namespace plb
