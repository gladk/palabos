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

#ifndef SWIG_MULTI_BLOCK_LATTICE_INFO_2D_HH
#define SWIG_MULTI_BLOCK_LATTICE_INFO_2D_HH

#include "plbWrapper/lattice/multiBlockInfo2D.h"

namespace plb {

template<typename T, template<typename T> class Descriptor>
MultiBlockLatticeInfo<T,Descriptor>::MultiBlockLatticeInfo (
        MultiBlockLattice2D<T,Descriptor> const& multiBlock )
    : nx(), ny(),
      numBlocks(),
      numAllocatedCells()
{
    getMultiBlockInfo(multiBlock, nx, ny, numBlocks,
                      smallest, largest, numAllocatedCells);
}

template<typename T, template<typename T> class Descriptor>
plint MultiBlockLatticeInfo<T,Descriptor>::getNx() const {
    return nx;
}

template<typename T, template<typename T> class Descriptor>
plint MultiBlockLatticeInfo<T,Descriptor>::getNy() const {
    return ny;
}

template<typename T, template<typename T> class Descriptor>
plint MultiBlockLatticeInfo<T,Descriptor>::getNumBlocks() const {
    return numBlocks;
}

template<typename T, template<typename T> class Descriptor>
plint MultiBlockLatticeInfo<T,Descriptor>::getNumAllocatedCells() const {
    return numAllocatedCells;
}

template<typename T, template<typename T> class Descriptor>
Box2D MultiBlockLatticeInfo<T,Descriptor>::getSmallestBlock() const {
    return smallest;
}

template<typename T, template<typename T> class Descriptor>
Box2D MultiBlockLatticeInfo<T,Descriptor>::getLargestBlock() const {
    return largest;
}

}  // namespace plb

#endif  // SWIG_MULTI_BLOCK_LATTICE_INFO_2D_HH
