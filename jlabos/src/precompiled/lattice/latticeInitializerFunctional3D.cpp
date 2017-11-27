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

#include "dataProcessors/dataInitializerFunctional3D.h"
#include "dataProcessors/dataInitializerFunctional3D.hh"
#include "latticeBoltzmann/nearestNeighborLattices3D.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.hh"

#ifdef COMPILE_3D

namespace plb {

template class InstantiateDynamicsFunctional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class SetConstBoundaryVelocityFunctional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class SetConstBoundaryDensityFunctional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class IniConstEquilibriumFunctional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class InstantiateCompositeDynamicsFunctional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class SetExternalScalarFunctional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class SetExternalVectorFunctional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;

}  // namespace plb

#endif  // COMPILE_3D
