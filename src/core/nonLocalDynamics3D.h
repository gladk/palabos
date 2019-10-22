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

#ifndef NON_LOCAL_DYNAMICS_3D_H
#define NON_LOCAL_DYNAMICS_3D_H

#include "core/globalDefs.h"
#include "core/dynamics.h"
#include "atomicBlock/blockLattice3D.h"

namespace plb {

template<typename T, template<typename U> class Descriptor>
class NonLocalDynamics3D : public CompositeDynamics<T,Descriptor> {
public:
    NonLocalDynamics3D(Dynamics<T,Descriptor>* baseDynamics_);
    virtual bool isNonLocal() const;
    virtual void prepareCollision(Cell<T,Descriptor>& cell);

    virtual void nonLocalAction(plint iX, plint iY, plint iZ, BlockLattice3D<T,Descriptor>& lattice) =0;
};

/*
template<typename T, template<typename U> class Descriptor>
class NonLocalBoundaryDynamics3D : public NonLocalDynamics3D<T,Descriptor> {
public:
    NonLocalBoundaryDynamics3D(Dynamics<T,Descriptor>* baseDynamics_);
    virtual bool isBoundary() const;

    virtual void nonLocalAction(plint iX, plint iY, plint iZ, BlockLattice3D<T,Descriptor>& lattice);
    virtual void boundaryCompletion(plint iX, plint iY, plint iZ, BlockLattice3D<T,Descriptor>& lattice) =0;
};
*/

}  // namespace plb

#endif  // NON_LOCAL_DYNAMICS_3D_H

