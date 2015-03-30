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
 * Operations on the 3D multiblock -- header file.
 */
#ifndef MULTI_BLOCK_OPERATIONS_3D_H
#define MULTI_BLOCK_OPERATIONS_3D_H

#include "core/globalDefs.h"
#include "core/blockStatistics.h"
#include "core/geometry3D.h"
#include <vector>

namespace plb {

class MultiBlock3D;
struct DataProcessorGenerator3D;
class ReductiveDataProcessorGenerator3D;

void executeDataProcessor( DataProcessorGenerator3D const& generator,
                           std::vector<MultiBlock3D*> multiBlocks );

void executeDataProcessor( DataProcessorGenerator3D const& generator,
                           MultiBlock3D& object );

void executeDataProcessor( DataProcessorGenerator3D const& generator,
                           MultiBlock3D& object1, MultiBlock3D& object2 );


void executeDataProcessor( ReductiveDataProcessorGenerator3D& generator,
                           std::vector<MultiBlock3D*> multiBlocks );

void executeDataProcessor( ReductiveDataProcessorGenerator3D& generator,
                           MultiBlock3D& object );

void executeDataProcessor( ReductiveDataProcessorGenerator3D& generator,
                           MultiBlock3D& object1, MultiBlock3D& object2 );


void addInternalProcessor( DataProcessorGenerator3D const& generator,
                           std::vector<MultiBlock3D*> multiBlocks, plint level=0 );

void addInternalProcessor( DataProcessorGenerator3D const& generator, MultiBlock3D& actor,
                           std::vector<MultiBlock3D*> multiBlockArgs, plint level=0 );

void addInternalProcessor( DataProcessorGenerator3D const& generator,
                           MultiBlock3D& object, plint level=0 );

void addInternalProcessor( DataProcessorGenerator3D const& generator,
                           MultiBlock3D& object1, MultiBlock3D& object2,
                           plint level=0 );

} // namespace plb

#endif  // MULTI_BLOCK_OPERATIONS_3D_H
