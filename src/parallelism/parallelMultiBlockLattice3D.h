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
 * Helper classes for parallel 2D multiblock lattice -- header file.
 */
#ifndef PARALLEL_MULTI_BLOCK_LATTICE_3D_H
#define PARALLEL_MULTI_BLOCK_LATTICE_3D_H

#include "core/globalDefs.h"
#include "parallelism/parallelBlockCommunicator3D.h"
#include "multiBlock/multiBlockLattice3D.h"

#ifdef PLB_MPI_PARALLEL

namespace plb {

template<typename T, template<typename U> class Descriptor>
class ParallelCellAccess3D : public MultiCellAccess3D<T,Descriptor> {
public:
    ParallelCellAccess3D();
    virtual ~ParallelCellAccess3D();
    virtual Cell<T,Descriptor>& getDistributedCell (
            plint iX, plint iY, plint iZ,
            MultiBlockManagement3D const& multiBlockManagement,
            std::map<plint,BlockLattice3D<T,Descriptor>*>& lattices );
    virtual Cell<T,Descriptor> const& getDistributedCell (
            plint iX, plint iY, plint iZ,
            MultiBlockManagement3D const& multiBlockManagement,
            std::map<plint,BlockLattice3D<T,Descriptor>*> const& lattices ) const;
    virtual void broadCastCell(Cell<T,Descriptor>& cell, plint fromBlock,
                               MultiBlockManagement3D const& multiBlockManagement) const;
    ParallelCellAccess3D<T,Descriptor>* clone() const;
private:
    mutable Cell<T,Descriptor> distributedCell;
    mutable std::vector<Cell<T,Descriptor>*> baseCells;
    mutable std::vector<Cell<T,Descriptor> const*> constBaseCells;
    mutable Dynamics<T,Descriptor>* parallelDynamics;
};

}  // namespace plb

#endif  // PLB_MPI_PARALLEL

#endif  // PARALLEL_MULTI_BLOCK_LATTICE_3D_H
