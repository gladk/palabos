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
 * Helper classes for serial 3D multiblock lattice -- generic implementation.
 */
#ifndef SERIAL_MULTI_BLOCK_LATTICE_3D_HH
#define SERIAL_MULTI_BLOCK_LATTICE_3D_HH

#include "multiBlock/serialMultiBlockLattice3D.h"
#include "multiBlock/multiBlockManagement3D.h"
#include "atomicBlock/blockLattice3D.h"

namespace plb {

////////////////////// Class SerialCellAccess3D /////////////////////

template<typename T, template<typename U> class Descriptor>
SerialCellAccess3D<T,Descriptor>::SerialCellAccess3D()
  : locatedBlock(0)
{ }

template<typename T, template<typename U> class Descriptor>
void SerialCellAccess3D<T,Descriptor>::broadCastCell(Cell<T,Descriptor>& cell, plint fromBlock,
                                                  MultiBlockManagement3D const& multiBlockManagement) const 
{
    // Nothing to do in the serial case
}

template<typename T, template<typename U> class Descriptor>
Cell<T,Descriptor>& SerialCellAccess3D<T,Descriptor>::getDistributedCell (
        plint iX, plint iY, plint iZ,
        MultiBlockManagement3D const& multiBlockManagement,
        std::map<plint,BlockLattice3D<T,Descriptor>*>& lattices )
{
    plint localX, localY, localZ;
#ifdef PLB_DEBUG
    bool ok =
#endif
        multiBlockManagement.findInLocalBulk
                  (iX,iY,iZ, locatedBlock, localX, localY, localZ);
    PLB_PRECONDITION( ok );
    return lattices[locatedBlock] -> get(localX,localY,localZ);
}

template<typename T, template<typename U> class Descriptor>
Cell<T,Descriptor> const& SerialCellAccess3D<T,Descriptor>::getDistributedCell (
        plint iX, plint iY, plint iZ,
        MultiBlockManagement3D const& multiBlockManagement,
        std::map<plint,BlockLattice3D<T,Descriptor>*> const& lattices ) const
{
    plint localX, localY, localZ;
#ifdef PLB_DEBUG
    bool ok =
#endif
        multiBlockManagement.findInLocalBulk
            (iX,iY,iZ, locatedBlock, localX, localY, localZ);
    PLB_PRECONDITION( ok );
    return lattices.find(locatedBlock)->second -> get(localX,localY,localZ);
}

template<typename T, template<typename U> class Descriptor>
SerialCellAccess3D<T,Descriptor>* SerialCellAccess3D<T,Descriptor>::clone() const {
    return new SerialCellAccess3D<T,Descriptor>;
}

}  // namespace plb

#endif  // SERIAL_MULTI_BLOCK_LATTICE_3D_HH
