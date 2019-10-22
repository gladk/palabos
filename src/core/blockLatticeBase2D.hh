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
 * Base class for the 2D BlockLattice and MultiBlockLattice -- generic implementation.
 */
#ifndef BLOCK_LATTICE_BASE_2D_HH
#define BLOCK_LATTICE_BASE_2D_HH

#include "core/blockLatticeBase2D.h"
#include "core/latticeStatistics.h"
#include "atomicBlock/dataProcessorWrapper2D.h"
#include <cmath>

namespace plb {

/////////// class BlockLatticeBase2D //////////////////////////////

template<typename T, template<typename U> class Descriptor>
BlockLatticeBase2D<T,Descriptor>::BlockLatticeBase2D()
{ }

template<typename T, template<typename U> class Descriptor>
BlockLatticeBase2D<T,Descriptor>::~BlockLatticeBase2D()
{ }

template<typename T, template<typename U> class Descriptor>
void BlockLatticeBase2D<T,Descriptor>::swap(BlockLatticeBase2D<T,Descriptor>& rhs) {
    std::swap(timeCounter, rhs.timeCounter);
}

template<typename T, template<typename U> class Descriptor>
TimeCounter& BlockLatticeBase2D<T,Descriptor>::getTimeCounter() {
    return timeCounter;
}

template<typename T, template<typename U> class Descriptor>
TimeCounter const& BlockLatticeBase2D<T,Descriptor>::getTimeCounter() const {
    return timeCounter;
}

}  // namespace plb

#endif
