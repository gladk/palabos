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

#include "multiBlock/multiDataProcessorWrapper2D.h"
#include "multiBlock/multiBlockLattice2D.h"
#include "multiBlock/multiDataField2D.h"
#include "multiBlock/multiBlockOperations2D.h"
#include "atomicBlock/dataProcessor2D.h"
#include "core/plbDebug.h"

namespace plb {

/* *************** BoxProcessing2D, general case *************************** */

void applyProcessingFunctional(BoxProcessingFunctional2D* functional,
                               Box2D domain, std::vector<MultiBlock2D*> multiBlocks)
{
    executeDataProcessor( BoxProcessorGenerator2D(functional, domain),
                          multiBlocks );
}

void integrateProcessingFunctional(BoxProcessingFunctional2D* functional,
                                   Box2D domain,
                                   std::vector<MultiBlock2D*> multiBlocks,
                                   plint level)
{
    addInternalProcessor( BoxProcessorGenerator2D(functional, domain),
                          multiBlocks, level );
}


/* *************** DotProcessing, general case ***************************** */

void applyProcessingFunctional(DotProcessingFunctional2D* functional,
                               DotList2D const& dotList,
                               std::vector<MultiBlock2D*> multiBlocks)
{
    executeDataProcessor( DotProcessorGenerator2D(functional, dotList),
                          multiBlocks );
}

void integrateProcessingFunctional(DotProcessingFunctional2D* functional,
                                   DotList2D const& dotList,
                                   std::vector<MultiBlock2D*> multiBlocks,
                                   plint level)
{
    addInternalProcessor( DotProcessorGenerator2D(functional, dotList),
                          multiBlocks, level );
}

/* *************** BoundedBoxProcessing2D, general case *************************** */

void applyProcessingFunctional(BoundedBoxProcessingFunctional2D* functional,
                               Box2D domain, std::vector<MultiBlock2D*> multiBlocks,
                               plint boundaryWidth )
{
    std::vector<BoxProcessorGenerator2D*> generators;
    functional -> getGenerators(domain, boundaryWidth, generators);
    delete functional;
    for (pluint iGen=0; iGen<generators.size(); ++iGen) {
        executeDataProcessor( *generators[iGen], multiBlocks );
        delete generators[iGen];
    }
}

void integrateProcessingFunctional(BoundedBoxProcessingFunctional2D* functional,
                                   Box2D domain, std::vector<MultiBlock2D*> multiBlocks,
                                   plint boundaryWidth, plint level)
{
    std::vector<BoxProcessorGenerator2D*> generators;
    functional -> getGenerators(domain, boundaryWidth, generators);
    delete functional;
    for (pluint iGen=0; iGen<generators.size(); ++iGen) {
        addInternalProcessor( *generators[iGen], multiBlocks, level );
        delete generators[iGen];
    }
}

}  // namespace plb
