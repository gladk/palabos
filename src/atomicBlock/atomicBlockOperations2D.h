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
 * Operations on the 2D multiblock -- header file.
 *
 * The functions here implement algorithms for the execution of a
 * data-processor on a block, or the inclusion into a block. Although
 * these functions can be used by an end-user, they are inconvenient.
 * It is better to use the wrappers declared in dataProcessorWrapper2D.h .
 */

#ifndef ATOMIC_BLOCK_OPERATIONS_2D_H
#define ATOMIC_BLOCK_OPERATIONS_2D_H

#include "core/globalDefs.h"
#include "atomicBlock/dataProcessor2D.h"
#include "atomicBlock/atomicBlock2D.h"
#include <vector>

namespace plb {

void executeDataProcessor( DataProcessorGenerator2D const& generator,
                           std::vector<AtomicBlock2D*> objects );

void executeDataProcessor( DataProcessorGenerator2D const& generator,
                           AtomicBlock2D& object );

void executeDataProcessor( DataProcessorGenerator2D const& generator,
                           AtomicBlock2D& object1, AtomicBlock2D& object2 );


void executeDataProcessor( ReductiveDataProcessorGenerator2D& generator,
                           std::vector<AtomicBlock2D*> objects );

void executeDataProcessor( ReductiveDataProcessorGenerator2D& generator,
                           AtomicBlock2D& object );

void executeDataProcessor( ReductiveDataProcessorGenerator2D& generator,
                           AtomicBlock2D& object1, AtomicBlock2D& object2 );


void addInternalProcessor( DataProcessorGenerator2D const& generator,
                           std::vector<AtomicBlock2D*> objects, plint level=0 );

void addInternalProcessor( DataProcessorGenerator2D const& generator,
                           AtomicBlock2D& object, plint level=0 );

void addInternalProcessor( DataProcessorGenerator2D const& generator,
                           AtomicBlock2D& object1, AtomicBlock2D& object2,
                           plint level=0 );

} // namespace plb

#endif  // ATOMIC_BLOCK_OPERATIONS_2D_H
