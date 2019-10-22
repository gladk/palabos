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
 * Helper classes for serial 2D multiblock lattice -- header file.
 */
#ifndef SERIAL_MULTI_BLOCK_LATTICE_2D_H
#define SERIAL_MULTI_BLOCK_LATTICE_2D_H

#include "core/globalDefs.h"
#include "multiBlock/multiBlockLattice2D.h"
#include <vector>

namespace plb {

template<typename T, template<typename U> class Descriptor>
class SerialCellAccess2D : public MultiCellAccess2D<T,Descriptor> {
public:
    SerialCellAccess2D();
    virtual void broadCastCell(Cell<T,Descriptor>& cell, plint fromBlock,
                               MultiBlockManagement2D const& multiBlockManagement) const;
    Cell<T,Descriptor>& getDistributedCell (
            plint iX, plint iY,
            MultiBlockManagement2D const& multiBlockManagement,
            std::map<plint,BlockLattice2D<T,Descriptor>*>& lattices );
    Cell<T,Descriptor> const& getDistributedCell (
            plint iX, plint iY,
            MultiBlockManagement2D const& multiBlockManagement,
            std::map<plint,BlockLattice2D<T,Descriptor>*> const& lattices ) const;
    virtual SerialCellAccess2D<T,Descriptor>* clone() const;
private:
    mutable plint locatedBlock;
};

}  // namespace plb

#endif  // SERIAL_MULTI_BLOCK_LATTICE_2D_H
