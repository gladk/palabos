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
 * Generator functions for all types of multi-blocks, to make them accessible to SWIG;
 * generic implementation.
 */
#ifndef MULTI_BLOCK_GENERATOR_3D_HH
#define MULTI_BLOCK_GENERATOR_3D_HH

#include "plbWrapper/lattice/multiBlockGenerator3D.h"
#include "multiBlock/defaultMultiBlockPolicy3D.h"

namespace plb {

template<typename T, template<typename U> class Descriptor>
MultiBlockLattice3D<T,Descriptor>*
    generateMultiBlockLattice3D(Box3D const& domain, Dynamics<T,Descriptor> const* backgroundDynamics)
{
    MultiBlockLattice3D<T,Descriptor>* lattice =
        new MultiBlockLattice3D<T,Descriptor> (
                defaultMultiBlockPolicy3D().getMultiBlockManagement(domain,Descriptor<T>::vicinity),
                defaultMultiBlockPolicy3D().getBlockCommunicator(),
                defaultMultiBlockPolicy3D().getCombinedStatistics(),
                defaultMultiBlockPolicy3D().getMultiCellAccess<T,Descriptor>(),
                backgroundDynamics->clone() );
    lattice->periodicity().toggleAll(true);
    return lattice;
}

template<typename T1, template<typename U> class Descriptor, typename T2>
MultiNTensorField3D<T2>*
    generateNTensorFieldFromLattice3D (
            MultiBlockLattice3D<T1,Descriptor> const& lattice,
            Box3D const& domain, plint ndim )
{
    MultiNTensorField3D<T2>* field = new MultiNTensorField3D<T2>(ndim, lattice, domain);
    field->periodicity() = lattice.periodicity();
    return field;
}

}  // namespace plb

#endif  // MULTI_BLOCK_GENERATOR_3D_HH
