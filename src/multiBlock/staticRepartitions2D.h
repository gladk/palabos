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
 * Utilities for 2D multi data distributions -- header file.
 */

#ifndef STATIC_REPARTITIONS_2D_H
#define STATIC_REPARTITIONS_2D_H

#include "parallelism/mpiManager.h"
#include "core/globalDefs.h"
#include "atomicBlock/dataField2D.h"
#include "multiBlock/sparseBlockStructure2D.h"

namespace plb {

/// A 2D field of scalar values used to indicate the type of the cells.
/// Any positive value indicates an active (bulk, boundary) cell, 
/// while zero indicates a non-active (no-dynamics) cell
typedef ScalarField2D<unsigned char> CellTypeField2D;

/// Create a regular data distribution confined by domain.
SparseBlockStructure2D createRegularDistribution2D (
        plint nx, plint ny, plint numBlocksX, plint numBlocksY );

/// Create a nx-by-ny data distribution
SparseBlockStructure2D createRegularDistribution2D (
        Box2D const& domain, plint numBlocksX, plint numBlocksY );

/// Create a data distribution with regular blocks, as evenly distributed as possible
SparseBlockStructure2D
    createRegularDistribution2D(plint nx, plint ny,
                                int numProc = global::mpi().getSize());

/// Create a data distribution with regular blocks, as evenly distributed as possible
SparseBlockStructure2D createRegularDistribution2D(Box2D const& domain,
                                                   int numProc = global::mpi().getSize());

/// Re-create a distribution by covering it with regular blocks.
SparseBlockStructure2D reparallelize(SparseBlockStructure2D const& originalStructure,
                                     plint blockLx, plint blockLy);

/// Re-create a distribution by covering it with regular blocks.
SparseBlockStructure2D reparallelize(SparseBlockStructure2D const& originalStructure);

/// Create a data distribution by slicing the domain (a block of nX*nY cells
/// as defined by cellTypeField) into numBlocks blocks along the x-direction. 
/// The x-extent of the blocks is chosen such as to obtain an approximately 
/// equal number of active cells in each block.
SparseBlockStructure2D createXSlicedDistribution2D (
        CellTypeField2D const& cellTypeField,
        plint numBlocks );

/// cf above.
SparseBlockStructure2D createYSlicedDistribution2D (
        CellTypeField2D const& cellTypeField,
        plint numBlocks );
        
/// Create x-sliced data distribution, balancing the number of active cells between blocks,
/// implicitly setting numBlocks = #processors
SparseBlockStructure2D createXSlicedDistribution2D (
        CellTypeField2D const& cellTypeField );

/// cf above
SparseBlockStructure2D createYSlicedDistribution2D (
        CellTypeField2D const& cellTypeField );

}  // namespace plb

#endif  // STATIC_REPARTITIONS_2D_H
