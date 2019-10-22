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

/* Main author: Daniel Lagrava
 **/

/** \file
 * Operations on the 2D multiGrid -- header file.
 */
#ifndef MULTI_GRID_OPERATIONS_2D_H
#define MULTI_GRID_OPERATIONS_2D_H

#include "core/globalDefs.h"
#include "core/blockStatistics.h"
#include "core/geometry2D.h"
#include "multiGrid/multiGrid2D.h"
#include <vector>

namespace plb {

class MultiBlock2D;
struct DataProcessorGenerator2D;
class ReductiveDataProcessorGenerator2D;

void executeDataProcessor( DataProcessorGenerator2D const& generator,
                           std::vector<MultiGrid2D*> multiGrids,
                           plint referenceLevel );


void executeDataProcessor( DataProcessorGenerator2D const& generator,
                           MultiGrid2D& object, plint referenceLevel );

void executeDataProcessor( DataProcessorGenerator2D const& generator,
                           MultiGrid2D& object1, MultiGrid2D& object2,
                           plint referenceLevel );


void executeDataProcessor( ReductiveDataProcessorGenerator2D& generator,
                           std::vector<MultiGrid2D*> multiGrids,
                           plint referenceLevel );

void executeDataProcessor( ReductiveDataProcessorGenerator2D& generator,
                           MultiGrid2D& object,
                           plint referenceLevel );

void executeDataProcessor( ReductiveDataProcessorGenerator2D& generator,
                           MultiGrid2D& object1, MultiGrid2D& object2,
                           plint referenceLevel );


void addInternalProcessor( DataProcessorGenerator2D const& generator,
                           std::vector<MultiGrid2D*> multiGrids,
                           plint referenceLevel,
                           plint processorLevel=0 );

void addInternalProcessor( DataProcessorGenerator2D const& generator,
                           MultiGrid2D& object,
                           plint referenceLevel,
                           plint processorLevel=0 );

void addInternalProcessor( DataProcessorGenerator2D const& generator,
                           MultiGrid2D& object1, MultiGrid2D& object2,
                           plint referenceLevel,
                           plint processorLevel=0 );

} // namespace plb

#endif  // MULTI_GRID_OPERATIONS_2D_H
