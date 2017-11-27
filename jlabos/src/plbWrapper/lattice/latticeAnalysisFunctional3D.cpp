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
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- template instantiation.
 */
#include "plbWrapper/lattice/latticeAnalysisFunctional3D.h"
#include "plbWrapper/lattice/latticeAnalysisFunctional3D.hh"
#include "latticeBoltzmann/nearestNeighborLattices3D.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.hh"

namespace plb {
    
template class N_BoxDensityFunctional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class N_BoxKineticEnergyFunctional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class N_BoxVelocityComponentFunctional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class N_BoxVelocityNormFunctional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class N_BoxVelocityFunctional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class N_BoxPopulationFunctional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class N_BoxPopulationsFunctional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class N_BoxPiNeqFunctional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class N_BoxShearStressFunctional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class N_BoxStrainRateFromStressFunctional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class UPO_Rhs_Functional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class UPO_ApplyJ_Functional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class UPO_EnergyDerivative_Functional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;

template class Masked_N_BoxDensityFunctional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class Masked_N_BoxKineticEnergyFunctional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class Masked_N_BoxVelocityComponentFunctional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class Masked_N_BoxVelocityNormFunctional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class Masked_N_BoxVelocityFunctional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class Masked_N_BoxPopulationFunctional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class Masked_N_BoxPopulationsFunctional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class Masked_N_BoxPiNeqFunctional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class Masked_N_BoxShearStressFunctional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class Masked_N_BoxStrainRateFromStressFunctional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class Masked_UPO_Rhs_Functional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class Masked_UPO_ApplyJ_Functional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class Masked_UPO_EnergyDerivative_Functional3D<FLOAT_T, descriptors::DESCRIPTOR_3D>;

}

#endif  // COMPILE_3D
