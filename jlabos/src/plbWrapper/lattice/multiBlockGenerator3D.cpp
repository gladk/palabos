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

#ifdef COMPILE_3D

/** \file
 * Generator functions for all types of multi-blocks, to make them accessible to SWIG;
 * template instantiation.
 */
#include "plbWrapper/lattice/multiBlockGenerator3D.h"
#include "plbWrapper/lattice/multiBlockGenerator3D.hh"
#include "latticeBoltzmann/nearestNeighborLattices3D.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.hh"


namespace plb {

template
MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>*
    generateMultiBlockLattice3D (
            Box3D const& domain,
            Dynamics<FLOAT_T,descriptors::DESCRIPTOR_3D> const* backgroundDynamics );

template
MultiNTensorField3D<FLOAT_T>*
    generateNTensorFieldFromLattice3D (
            MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D> const& lattice,
            Box3D const& domain, plint ndim );

template
MultiNTensorField3D<int>*
    generateNTensorFieldFromLattice3D (
            MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D> const& lattice,
            Box3D const& domain, plint ndim);

}  // namespace plb

#endif  // COMPILE_3D
