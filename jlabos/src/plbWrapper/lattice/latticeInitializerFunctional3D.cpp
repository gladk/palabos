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

#include "plbWrapper/lattice/latticeInitializerFunctional3D.h"
#include "plbWrapper/lattice/latticeInitializerFunctional3D.hh"
#include "latticeBoltzmann/nearestNeighborLattices3D.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.hh"

namespace plb {

template class MaskedIniDynamicsFunctional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class N_IniBoundaryVelocityFunctional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class Masked_N_IniBoundaryVelocityFunctional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class N_IniEquilibriumFunctional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class Masked_N_IniEquilibriumFunctional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class N_IniPopulationsFunctional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class Masked_N_IniPopulationsFunctional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class N_IniConstPopulationsFunctional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class Masked_N_IniConstPopulationsFunctional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;

}  // namespace plb

#endif  // COMPILE_3D
