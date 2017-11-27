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

#include "plbWrapper/lattice/latticeInitializerWrapper3D.h"
#include "plbWrapper/lattice/latticeInitializerWrapper3D.hh"
#include "latticeBoltzmann/nearestNeighborLattices3D.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.hh"

namespace plb {

template void pypalDefineDynamics<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice,
        Box3D domain, Dynamics<FLOAT_T,descriptors::DESCRIPTOR_3D>* dynamics);
template void maskedDefineDynamics<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice, MultiNTensorField3D<int>& mask,
        Box3D domain, Dynamics<FLOAT_T,descriptors::DESCRIPTOR_3D>* dynamics);
template void setBoundaryVelocity<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        MultiBlockLattice3D<FLOAT_T, descriptors::DESCRIPTOR_3D>& lattice,
        FLOAT_T* velocity, int numDimIs2, Box3D domain );
template void setBoundaryVelocity<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        MultiBlockLattice3D<FLOAT_T, descriptors::DESCRIPTOR_3D>& lattice,
        MultiNTensorField3D<FLOAT_T>& velocity, Box3D domain );
template void maskedSetBoundaryVelocity<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        MultiBlockLattice3D<FLOAT_T, descriptors::DESCRIPTOR_3D>& lattice,
        MultiNTensorField3D<FLOAT_T>& velocity,
        MultiNTensorField3D<int>& mask, Box3D domain );
template void initializeAtEquilibrium<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        MultiBlockLattice3D<FLOAT_T, descriptors::DESCRIPTOR_3D>& lattice, FLOAT_T rho,
        FLOAT_T* velocity, int numDimIs2, Box3D domain );
template void initializeAtEquilibrium<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        MultiBlockLattice3D<FLOAT_T, descriptors::DESCRIPTOR_3D>& lattice,
        MultiNTensorField3D<FLOAT_T>& density,
        MultiNTensorField3D<FLOAT_T>& velocity,
        Box3D domain );
template void maskedInitializeAtEquilibrium<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        MultiBlockLattice3D<FLOAT_T, descriptors::DESCRIPTOR_3D>& lattice,
        MultiNTensorField3D<FLOAT_T>& density,
        MultiNTensorField3D<FLOAT_T>& velocity,
        MultiNTensorField3D<int>& mask,
        Box3D domain );
template void setExternalVector<FLOAT_T, descriptors::DESCRIPTOR_3D> (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice,
        int vectorStartsAt, FLOAT_T* externalVector, int numDimIs2,
        Box3D domain);
template void setPopulations<FLOAT_T,descriptors::DESCRIPTOR_3D> (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice,
        FLOAT_T* populations, int numDimIsQ, Box3D domain );
template void setPopulations<FLOAT_T,descriptors::DESCRIPTOR_3D> (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice,
        MultiNTensorField3D<FLOAT_T>& populations, Box3D domain );
template void maskedSetPopulations<FLOAT_T,descriptors::DESCRIPTOR_3D> (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice,
        FLOAT_T* populations, int numDimIsQ,
        MultiNTensorField3D<int>& mask, Box3D domain );
template void maskedSetPopulations<FLOAT_T,descriptors::DESCRIPTOR_3D> (
        MultiBlockLattice3D<FLOAT_T,descriptors::DESCRIPTOR_3D>& lattice,
        MultiNTensorField3D<FLOAT_T>& velocity,
        MultiNTensorField3D<int>& mask, Box3D domain );

}  // namespace plb

#endif  // COMPILE_3D
