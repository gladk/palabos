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

#include "atomicBlock/dataProcessingFunctional3D.h"
#include "atomicBlock/dataProcessingFunctional3D.hh"
#include "latticeBoltzmann/nearestNeighborLattices3D.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.hh"
#include "latticeBoltzmann/mrtLattices.h"
#include "latticeBoltzmann/mrtLattices.hh"
#include "multiPhysics/shanChenLattices3D.h"

namespace plb {

/* *************** Exotic lattices *********************************** */

template class LatticeBoxProcessingFunctional3D<FLOAT_T, descriptors::ShanChenD3Q19Descriptor>;
template class LatticeBoxProcessingFunctional3D<FLOAT_T, descriptors::ForcedShanChenD3Q19Descriptor>;
template class LatticeBoxProcessingFunctional3D<FLOAT_T, descriptors::MRTD3Q19Descriptor>;

template class BoxProcessingFunctional3D_L<FLOAT_T, descriptors::ShanChenD3Q19Descriptor>;
template class BoxProcessingFunctional3D_L<FLOAT_T, descriptors::ForcedShanChenD3Q19Descriptor>;
template class BoxProcessingFunctional3D_L<FLOAT_T, descriptors::MRTD3Q19Descriptor>;

/* *************** Boxed Data Processor functionals ****************** */

template class BoxProcessingFunctional3D_L<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class MaskedBoxProcessingFunctional3D_L<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class BoxProcessingFunctional3D_LL< FLOAT_T, descriptors::DESCRIPTOR_3D,
                                             FLOAT_T, descriptors::DESCRIPTOR_3D >;
template class BoxProcessingFunctional3D_LN<FLOAT_T, descriptors::DESCRIPTOR_3D, FLOAT_T>;
template class BoxProcessingFunctional3D_LN<FLOAT_T, descriptors::DESCRIPTOR_3D, int>;
template class MaskedBoxProcessingFunctional3D_LN<FLOAT_T, descriptors::DESCRIPTOR_3D, FLOAT_T>;
template class MaskedBoxProcessingFunctional3D_LN<FLOAT_T, descriptors::DESCRIPTOR_3D, int>;

template class BoxProcessingFunctional3D_LS<FLOAT_T, descriptors::DESCRIPTOR_3D, FLOAT_T>;
template class BoxProcessingFunctional3D_LS<FLOAT_T, descriptors::DESCRIPTOR_3D, int>;

/* *************** Bounded Boxed Data Processor functionals ****************** */

template class BoundedBoxProcessingFunctional3D_L<FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class BoundedBoxProcessingFunctional3D_LL< FLOAT_T, descriptors::DESCRIPTOR_3D,
                                                    FLOAT_T, descriptors::DESCRIPTOR_3D>;
template class BoundedBoxProcessingFunctional3D_LN< FLOAT_T, descriptors::DESCRIPTOR_3D,
                                                    FLOAT_T >;
template class BoundedMaskedBoxProcessingFunctional3D_LN< FLOAT_T, descriptors::DESCRIPTOR_3D,
                                                          FLOAT_T >;
template class BoundedLatticeBoxProcessingFunctional3D <
    FLOAT_T, descriptors::DESCRIPTOR_3D >;

}  // namespace plb

#endif  // COMPILE_3D
