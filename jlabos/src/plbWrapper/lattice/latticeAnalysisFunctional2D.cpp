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

#ifdef COMPILE_2D

/** \file
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- template instantiation.
 */
#include "plbWrapper/lattice/latticeAnalysisFunctional2D.h"
#include "plbWrapper/lattice/latticeAnalysisFunctional2D.hh"
#include "latticeBoltzmann/nearestNeighborLattices2D.h"
#include "latticeBoltzmann/nearestNeighborLattices2D.hh"

namespace plb {
    
template class N_BoxDensityFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D>;
template class N_BoxKineticEnergyFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D>;
template class N_BoxVelocityComponentFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D>;
template class N_BoxVelocityNormFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D>;
template class N_BoxVelocityFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D>;
template class N_BoxPopulationFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D>;
template class N_BoxPopulationsFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D>;
template class N_BoxPiNeqFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D>;
template class N_BoxShearStressFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D>;
template class N_BoxStrainRateFromStressFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D>;

template class Masked_N_BoxDensityFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D>;
template class Masked_N_BoxKineticEnergyFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D>;
template class Masked_N_BoxVelocityComponentFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D>;
template class Masked_N_BoxVelocityNormFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D>;
template class Masked_N_BoxVelocityFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D>;
template class Masked_N_BoxPopulationFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D>;
template class Masked_N_BoxPopulationsFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D>;
template class Masked_N_BoxPiNeqFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D>;
template class Masked_N_BoxShearStressFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D>;
template class Masked_N_BoxStrainRateFromStressFunctional2D<FLOAT_T, descriptors::DESCRIPTOR_2D>;

}

#endif  // COMPILE_2D
