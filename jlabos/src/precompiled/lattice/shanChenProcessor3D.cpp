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
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- template instantiation.
 */
#ifdef COMPILE_3D

#include "multiPhysics/shanChenProcessor3D.h"
#include "multiPhysics/shanChenProcessor3D.hh"
#include "multiPhysics/shanChenLattices3D.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.hh"

namespace plb {

#if NUMBIT_3D == 19
template class ShanChenMultiComponentProcessor3D<FLOAT_T, descriptors::ShanChenD3Q19Descriptor>;
template class ShanChenMultiComponentProcessor3D<FLOAT_T, descriptors::ForcedShanChenD3Q19Descriptor>;

template class ShanChenSingleComponentProcessor3D<FLOAT_T, descriptors::ShanChenD3Q19Descriptor>;
template class ShanChenSingleComponentProcessor3D<FLOAT_T, descriptors::ForcedShanChenD3Q19Descriptor>;
#endif

}

#endif  // COMPILE_3D
