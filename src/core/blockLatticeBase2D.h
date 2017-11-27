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
 * Base class for the 2D BlockLattice and MultiBlockLattice -- header file.
 */
#ifndef BLOCK_LATTICE_BASE_2D_H
#define BLOCK_LATTICE_BASE_2D_H

#include "core/globalDefs.h"
#include "core/geometry2D.h"
#include "core/blockStatistics.h"

namespace plb {

template<typename T, template<typename U> class Descriptor> struct Dynamics;
template<typename T, template<typename U> class Descriptor> class Cell;

/// An interface to all the variants of (more or less) regular lattices.
template<typename T, template<typename U> class Descriptor>
class BlockLatticeBase2D {
public:
    BlockLatticeBase2D();
    virtual ~BlockLatticeBase2D();
    void swap(BlockLatticeBase2D<T,Descriptor>& rhs);
public:
    virtual Cell<T,Descriptor>& get(plint iX, plint iY) =0;
    virtual Cell<T,Descriptor> const& get(plint iX, plint iY) const =0;
    virtual void specifyStatisticsStatus (Box2D domain, bool status) =0;
    virtual void collide(Box2D domain) =0;
    virtual void collide() =0;
    virtual void stream(Box2D domain) =0;
    virtual void stream() =0;
    virtual void collideAndStream(Box2D domain) =0;
    virtual void collideAndStream() =0;
    virtual void incrementTime() =0;
    TimeCounter& getTimeCounter();
    TimeCounter const& getTimeCounter() const;
private:
    TimeCounter timeCounter;
};

}  // namespace plb

#endif
