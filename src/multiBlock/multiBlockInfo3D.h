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
 * Helper functions for domain initialization -- header file.
 */
#include "core/globalDefs.h"
#include "multiBlock/multiBlock3D.h"
#include <string>


#ifndef MULTI_BLOCK_INFO_3D_H
#define MULTI_BLOCK_INFO_3D_H

namespace plb {

bool getMultiBlockInfo(MultiBlock3D const& multiBlock,
                       plint& nx, plint& ny, plint& nz, plint& numBlocks,
                       Box3D& smallest, Box3D& largest,
                       plint& numAllocatedCells);

std::string getMultiBlockInfo(MultiBlock3D const& multiBlock);

}  // namespace plb

#endif  // MULTI_BLOCK_INFO_3D_H
