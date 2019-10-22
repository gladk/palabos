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

#ifndef NON_LOCAL_DYNAMICS_2D_HH
#define NON_LOCAL_DYNAMICS_2D_HH

#include "core/globalDefs.h"
#include "core/nonLocalDynamics2D.h"

namespace plb {

template<typename T, template<typename U> class Descriptor>
NonLocalDynamics2D<T,Descriptor>::NonLocalDynamics2D(Dynamics<T,Descriptor>* baseDynamics_)
    : CompositeDynamics<T,Descriptor>(baseDynamics_, false)
{ }

template<typename T, template<typename U> class Descriptor>
bool NonLocalDynamics2D<T,Descriptor>::isNonLocal() const {
    return true;
}

template<typename T, template<typename U> class Descriptor>
void NonLocalDynamics2D<T,Descriptor>::prepareCollision(Cell<T,Descriptor>& cell) { }



template<typename T, template<typename U> class Descriptor>
NonLocalBoundaryDynamics2D<T,Descriptor>::NonLocalBoundaryDynamics2D(Dynamics<T,Descriptor>* baseDynamics_)
    : NonLocalDynamics2D<T,Descriptor>(baseDynamics_)
{ }

template<typename T, template<typename U> class Descriptor>
bool NonLocalBoundaryDynamics2D<T,Descriptor>::isBoundary() const {
    return true;
}

template<typename T, template<typename U> class Descriptor>
void NonLocalBoundaryDynamics2D<T,Descriptor>::nonLocalAction (
        plint iX, plint iY, BlockLattice2D<T,Descriptor>& lattice )
{
    boundaryCompletion(iX,iY, lattice);
}

}  // namespace plb

#endif  // NON_LOCAL_DYNAMICS_2D_HH

