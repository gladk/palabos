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

#ifndef DOMAIN_CLUSTERING_3D_H
#define DOMAIN_CLUSTERING_3D_H

#include "core/globalDefs.h"
#include "parallelism/mpiManager.h"
#include "multiBlock/sparseBlockStructure3D.h"
#include <vector>
#include <set>

namespace plb {

#ifdef PLB_MPI_PARALLEL

class BlockWithContact;

ExplicitThreadAttribution* optimalThreadAttribution(SparseBlockStructure3D const& sparseBlock, ThreadAttribution const& oldThreadAttribution);

std::vector<plint> optimalRepartition(SparseBlockStructure3D const& sparseBlock, std::vector<plint> localBlocks);

/// Provide explicitly a list of all blocks which, before re-attribution, are local (instead of providing
///   a thread-attribution. The body of this function is executed by all MPI threads except for the main one.
std::vector<plint> tryOptimalRepartition(SparseBlockStructure3D const& sparseBlock, std::vector<plint> localBlocks);

/// The body of this function is executed only by the main MPI thread.
std::vector<plint> mainProcessorRepartition(SparseBlockStructure3D const& sparseBlock, std::set<plint> localBlocks);

/// Execute one iteration step of the Las Vegas method.
void repartitionIter(SparseBlockStructure3D const& sparseBlock, std::set<plint>& localBlocks, int myPartner, bool exchangeAnyway);

/// On the current MPI thread compute the fitness of the current block repartition, i.e. the number
///   of contacts between local blocks, and returned detailed contact info for each block.
void evaluateNeighbors (
        SparseBlockStructure3D const& sparseBlock, std::set<plint> const& localBlocks,
        std::vector<BlockWithContact>& contacts, plint& totalContacts );

/// Count the number of contacts between local blocks. This could also be obtained through the function
///   evaluateNeighbors, which however would be slower.
plint getNumLocalNeighbors(SparseBlockStructure3D const& sparseBlock, std::set<plint> const& localBlocks);

/// Count the number of contacts for one given local block.
plint getNumLocalNeighbors(SparseBlockStructure3D const& sparseBlock, std::set<plint> const& localBlocks, plint blockId);

#endif // PLB_MPI_PARALLEL

}  // namespace plb

#endif  // DOMAIN_CLUSTERING_3D_H

