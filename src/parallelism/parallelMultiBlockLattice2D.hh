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
 * Helper classes for parallel 2D multiblock lattice -- generic implementation.
 */
#ifndef PARALLEL_MULTI_BLOCK_LATTICE_2D_HH
#define PARALLEL_MULTI_BLOCK_LATTICE_2D_HH

#include "parallelism/parallelMultiBlockLattice2D.h"
#include "parallelism/parallelDynamics.h"
#include "multiBlock/staticRepartitions2D.h"
#include "atomicBlock/blockLattice2D.h"


#ifdef PLB_MPI_PARALLEL

namespace plb {

/* *************** Class ParallelCellAccess2D ************************** */

template<typename T, template<typename U> class Descriptor>
ParallelCellAccess2D<T,Descriptor>::ParallelCellAccess2D()
    : parallelDynamics( 0 )
{ }

template<typename T, template<typename U> class Descriptor>
ParallelCellAccess2D<T,Descriptor>::~ParallelCellAccess2D() {
    delete parallelDynamics;
}

template<typename T, template<typename U> class Descriptor>
void ParallelCellAccess2D<T,Descriptor>::broadCastCell(Cell<T,Descriptor>& cell, plint fromBlock,
                                                    MultiBlockManagement2D const& multiBlockManagement) const
{
    const plint sizeOfCell = Descriptor<T>::q + Descriptor<T>::ExternalField::numScalars;
    char* cellData = new char[sizeOfCell*sizeof(T)];
    plint fromProc = multiBlockManagement.getThreadAttribution().getMpiProcess(fromBlock);
    if (global::mpi().getRank()==fromProc) {
        cell.serialize(cellData);
    }
    global::mpi().bCast(cellData, sizeOfCell, fromProc);
    cell.unSerialize(cellData);
    delete [] cellData;
}

template<typename T, template<typename U> class Descriptor>
Cell<T,Descriptor>& ParallelCellAccess2D<T,Descriptor>::getDistributedCell (
        plint iX, plint iY,
        MultiBlockManagement2D const& multiBlockManagement,
        std::map<plint,BlockLattice2D<T,Descriptor>*>& lattices )
{
    std::vector<plint> foundId;
    std::vector<plint> foundX, foundY;
    bool hasBulkCell = multiBlockManagement.findAllLocalRepresentations(iX,iY, foundId, foundX, foundY);
    baseCells.clear();
    for (pluint iBlock=0; iBlock<foundId.size(); ++iBlock) {
        plint foundBlock = foundId[iBlock];
        baseCells.push_back ( &lattices[foundBlock] -> get ( foundX[iBlock], foundY[iBlock] ) );
    }
    delete parallelDynamics;
    parallelDynamics = new ParallelDynamics<T,Descriptor>(baseCells, hasBulkCell);
    distributedCell.attributeDynamics(parallelDynamics);
    return distributedCell;
}

template<typename T, template<typename U> class Descriptor>
Cell<T,Descriptor> const& ParallelCellAccess2D<T,Descriptor>::getDistributedCell (
        plint iX, plint iY,
        MultiBlockManagement2D const& multiBlockManagement,
        std::map<plint,BlockLattice2D<T,Descriptor>*> const& lattices ) const
{
    std::vector<plint> foundId;
    std::vector<plint> foundX, foundY;
    bool hasBulkCell = multiBlockManagement.findAllLocalRepresentations(iX,iY, foundId, foundX, foundY);
    constBaseCells.clear();
    for (pluint iBlock=0; iBlock<foundId.size(); ++iBlock) {
        plint foundBlock = foundId[iBlock];
        typename std::map<plint,BlockLattice2D<T,Descriptor>*>::const_iterator it = lattices.find(foundBlock);
        constBaseCells.push_back ( &it->second -> get ( foundX[iBlock], foundY[iBlock] ) );
    }
    delete parallelDynamics;
    parallelDynamics = new ConstParallelDynamics<T,Descriptor>(constBaseCells, hasBulkCell);
    distributedCell.attributeDynamics(parallelDynamics);
    return distributedCell;
}

template<typename T, template<typename U> class Descriptor>
ParallelCellAccess2D<T,Descriptor>* ParallelCellAccess2D<T,Descriptor>::clone() const {
    return new ParallelCellAccess2D<T,Descriptor>;
}

}  // namespace plb

#endif  // PLB_MPI_PARALLEL

#endif  // PARALLEL_MULTI_BLOCK_LATTICE_2D_HH
