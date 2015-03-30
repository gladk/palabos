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

/** \file
 * Operations on the 2D multiblock -- header file.
 */
#ifndef MULTI_BLOCK_OPERATIONS_2D_H
#define MULTI_BLOCK_OPERATIONS_2D_H

#include "core/globalDefs.h"
#include "core/blockStatistics.h"
#include "core/geometry2D.h"
#include <vector>

namespace plb {

class MultiBlock2D;
struct DataProcessorGenerator2D;
class ReductiveDataProcessorGenerator2D;

void executeDataProcessor( DataProcessorGenerator2D const& generator,
                           std::vector<MultiBlock2D*> multiBlocks );

void executeDataProcessor( DataProcessorGenerator2D const& generator,
                           MultiBlock2D& object );

void executeDataProcessor( DataProcessorGenerator2D const& generator,
                           MultiBlock2D& object1, MultiBlock2D& object2 );


void executeDataProcessor( ReductiveDataProcessorGenerator2D& generator,
                           std::vector<MultiBlock2D*> multiBlocks );

void executeDataProcessor( ReductiveDataProcessorGenerator2D& generator,
                           MultiBlock2D& object );

void executeDataProcessor( ReductiveDataProcessorGenerator2D& generator,
                           MultiBlock2D& object1, MultiBlock2D& object2 );


void addInternalProcessor( DataProcessorGenerator2D const& generator,
                           std::vector<MultiBlock2D*> multiBlocks, plint level=0 );

void addInternalProcessor( DataProcessorGenerator2D const& generator,
                           MultiBlock2D& object, plint level=0 );

void addInternalProcessor( DataProcessorGenerator2D const& generator,
                           MultiBlock2D& object1, MultiBlock2D& object2,
                           plint level=0 );

} // namespace plb

#endif  // MULTI_BLOCK_OPERATIONS_2D_H
