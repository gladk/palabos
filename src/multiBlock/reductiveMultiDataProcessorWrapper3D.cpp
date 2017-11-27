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

#include "multiBlock/reductiveMultiDataProcessorWrapper3D.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiDataField3D.h"
#include "multiBlock/multiBlockOperations3D.h"
#include "atomicBlock/dataProcessor3D.h"
#include "core/plbDebug.h"

namespace plb {

/* *************** BoxProcessing3D, general case *************************** */

void applyProcessingFunctional(ReductiveBoxProcessingFunctional3D& functional,
                               Box3D domain, std::vector<MultiBlock3D*> multiBlocks)
{
    // Let the generator get off with a clone of the functional (because the generator
    //   owns its functional, whereas we still need our copy afterwards to get at the
    //   reducted values).
    ReductiveBoxProcessorGenerator3D generator(functional.clone(), domain);
    executeDataProcessor(generator, multiBlocks );
    // Recover reducted values from the generator's functional.
    functional.getStatistics() = generator.getFunctional().getStatistics();
}


/* *************** DotProcessing, general case ***************************** */

void applyProcessingFunctional(ReductiveDotProcessingFunctional3D& functional,
                               DotList3D const& dotList,
                               std::vector<MultiBlock3D*> multiBlocks)
{
    // Let the generator get off with a clone of the functional (because the generator
    //   owns its functional, whereas we still need our copy afterwards to get at the
    //   reducted values).
    ReductiveDotProcessorGenerator3D generator(functional.clone(), dotList);
    executeDataProcessor(generator, multiBlocks);
    // Recover reducted values from the generator's functional.
    functional.getStatistics() = generator.getFunctional().getStatistics();
}

/* *************** BoundedReductiveBoxProcessing3D, general case *********** */

void applyProcessingFunctional(BoundedReductiveBoxProcessingFunctional3D& functional,
                               Box3D domain, std::vector<MultiBlock3D*> multiBlocks,
                               plint boundaryWidth )
{
    std::vector<ReductiveBoxProcessorGenerator3D*> generators;
    functional.getGenerators(domain, boundaryWidth, generators);
    std::vector<BlockStatistics const*> individualStatistics(generators.size());
    for (pluint iGen=0; iGen<generators.size(); ++iGen) {
        executeDataProcessor( *generators[iGen], multiBlocks );
        individualStatistics[iGen] = &(generators[iGen]->getStatistics());
    }
    SerialCombinedStatistics().combine(individualStatistics, functional.getStatistics());
    for (pluint iGen=0; iGen<generators.size(); ++iGen) {
        delete generators[iGen];
    }
}

}  // namespace plb
