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

#include "dataProcessors/metaStuffWrapper2D.h"
#include "dataProcessors/metaStuffFunctional2D.h"
#include "dataProcessors/dataInitializerFunctional2D.h"
#include "dataProcessors/dataInitializerFunctional2D.hh"
#include "dataProcessors/dataInitializerWrapper2D.h"
#include "dataProcessors/dataInitializerWrapper2D.hh"
#include "atomicBlock/dataProcessorWrapper2D.h"
#include "atomicBlock/dataProcessorWrapper2D.hh"
#include "atomicBlock/dataProcessingFunctional2D.h"
#include "atomicBlock/dataProcessingFunctional2D.hh"
#include "multiBlock/multiDataProcessorWrapper2D.h"
#include "multiBlock/reductiveMultiDataProcessorWrapper2D.h"
#include "multiBlock/reductiveMultiDataProcessorWrapper2D.hh"
#include "multiBlock/multiDataProcessorWrapper2D.h"
#include "multiBlock/multiDataProcessorWrapper2D.hh"
#include "multiBlock/multiContainerBlock2D.h"


namespace plb {

bool allFlagsTrue(MultiBlock2D* multiBlock) {
    std::vector<MultiBlock2D*> singleBlockVector;
    singleBlockVector.push_back(multiBlock);
    AllFlagsTrueFunctional2D functional;
    applyProcessingFunctional (
            functional, multiBlock->getBoundingBox(),
            singleBlockVector );
    return functional.allTrue();
}

}  // namespace plb
