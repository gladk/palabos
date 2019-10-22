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
 * I/O routines for 3D multiblock -- header file.
 */
#ifndef MULTI_BLOCK_WRITER_3D_H
#define MULTI_BLOCK_WRITER_3D_H

#include "core/globalDefs.h"
#include "atomicBlock/atomicContainerBlock3D.h"
#include "multiBlock/multiBlock3D.h"
#include "core/serializer.h"
#include "io/plbFiles.h"

namespace plb {


namespace parallelIO {

void save( MultiBlock3D& multiBlock, FileName fName,
           bool dynamicContent = true );

void saveFull( MultiBlock3D& multiBlock, FileName fName,
               IndexOrdering::OrderingT=IndexOrdering::forward, bool appendMode=false );

/** Convert the multi-block into raw serialized data, but preserve the
 *  parallel distribution. This is a preparation for parallel IO. The
 *  IDs of the original blocks are not preserved; they are renumbered
 *  contiguously, in increasing order of the original IDs.
 *  @var offset: Position, in bytes, at which the data corresponding
 *               to a given atomic-block must be written into the file.
 *               The size of the vector equals the number of atomic-blocks.
 *  @var myBlockIds: New, contiguously numbered IDs of the atomic-blocks
 *                   which are local to the current MPI thread.
 *  @var data: The serialized data of the local atomic-blocks.
 **/
void dumpData( MultiBlock3D& multiBlock, bool dynamicContent,
               std::vector<plint>& offset, std::vector<plint>& myBlockIds,
               std::vector<std::vector<char> >& data );

void writeXmlSpec( MultiBlock3D& multiBlock, FileName fName,
                   std::vector<plint> const& offset, bool dynamicContent );

}  // namespace parallelIO

}  // namespace plb

#endif  // MULTI_BLOCK_WRITER_3D_H

