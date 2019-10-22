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
#ifndef MPI_PARALLEL_IO
#define MPI_PARALLEL_IO

#include "core/globalDefs.h"
#include "io/plbFiles.h"
#include <string>
#include <vector>

namespace plb {

namespace parallelIO {

void writeRawData( FileName fName, std::vector<plint> const& myBlockIds,
                   std::vector<plint> const& offset, std::vector<std::vector<char> >& data,
                   bool appendMode=false );

void loadRawData( FileName fName,  std::vector<plint> const& myBlockIds,
                  std::vector<plint> const& offset, std::vector<std::vector<char> >& data );

}  // namespace parallelIO

}  // namespace plb

#endif  // MULTI_BLOCK_WRITER_3D_H
