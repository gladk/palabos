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
 * Utilities for 3D multi data distributions -- header file.
 */

#ifndef STATIC_REPARTITIONS_3D_H
#define STATIC_REPARTITIONS_3D_H

#include "parallelism/mpiManager.h"
#include "core/globalDefs.h"
#include "atomicBlock/dataField3D.h"
#include "multiBlock/sparseBlockStructure3D.h"

namespace plb {

/// A 3D field of scalar values used to indicate the type of the cells.
/// Any positive value indicates an active (bulk, boundary) cell, 
/// while zero indicates a non-active (no-dynamics) cell
typedef ScalarField3D<unsigned char> CellTypeField3D;

/// Create a regular data distribution confined by domain.
SparseBlockStructure3D createRegularDistribution3D (
        plint nx, plint ny, plint nz, plint numBlocksX, plint numBlocksY, plint numBlocksZ );

/// Create a nx-by-ny-by-nz data distribution
SparseBlockStructure3D createRegularDistribution3D (
        Box3D const& domain, plint numBlocksX, plint numBlocksY, plint numBlocksZ );

/// Create a data distribution with regular blocks, as evenly distributed as possible
SparseBlockStructure3D
    createRegularDistribution3D(plint nx, plint ny, plint nz,
                                int numProc = global::mpi().getSize());

/// Create a data distribution with regular blocks, as evenly distributed as possible
SparseBlockStructure3D createRegularDistribution3D(Box3D const& domain,
                                                   int numProc = global::mpi().getSize());

/// Create a data distribution with regular blocks only at the y and z directions, as evenly distributed as possible
SparseBlockStructure3D createRegularDistributionYZ3D (
        plint nx, plint ny, plint nz, int numProc = global::mpi().getSize() );

/// Create a data distribution with regular blocks only at the x and z directions, as evenly distributed as possible
SparseBlockStructure3D createRegularDistributionXZ3D (
        plint nx, plint ny, plint nz, int numProc = global::mpi().getSize() );

/// Create a data distribution with regular blocks only at the x and y directions, as evenly distributed as possible
SparseBlockStructure3D createRegularDistributionXY3D (
        plint nx, plint ny, plint nz, int numProc = global::mpi().getSize() );

/// Re-create a distribution by covering it with regular blocks.
SparseBlockStructure3D reparallelize(SparseBlockStructure3D const& originalStructure,
                                     plint blockLx, plint blockLy, plint blockLz);

/// Re-create a distribution by covering it with regular blocks.
SparseBlockStructure3D reparallelize(SparseBlockStructure3D const& originalStructure);

/// Create a data distribution by slicing the domain (a block of nX*nY*nZ cells
/// as defined by cellTypeField) into numBlocks blocks along the x-direction. 
/// The x-extent of the blocks is chosen such as to obtain an approximately 
/// equal number of active cells in each block.
SparseBlockStructure3D createXSlicedDistribution3D (
        CellTypeField3D const& cellTypeField,
        plint numBlocks );

/// cf above.
SparseBlockStructure3D createYSlicedDistribution3D (
        CellTypeField3D const& cellTypeField,
        plint numBlocks );
        
/// cf above.
SparseBlockStructure3D createZSlicedDistribution3D (
        CellTypeField3D const& cellTypeField,
        plint numBlocks );

/// Create x-sliced data distribution, balancing the number of active cells between blocks,
/// implicitly setting numBlocks = #processors
SparseBlockStructure3D createXSlicedDistribution3D (
        CellTypeField3D const& cellTypeField );

/// cf above
SparseBlockStructure3D createYSlicedDistribution3D (
        CellTypeField3D const& cellTypeField );

/// cf above
SparseBlockStructure3D createZSlicedDistribution3D (
        CellTypeField3D const& cellTypeField );

}  // namespace plb

#endif  // STATIC_REPARTITIONS_3D_H
