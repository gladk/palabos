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

#include "atomicBlock/dataProcessingFunctional2D.h"
#include "atomicBlock/dataProcessingFunctional2D.hh"
#include "latticeBoltzmann/nearestNeighborLattices2D.h"
#include "latticeBoltzmann/nearestNeighborLattices2D.hh"
#include "latticeBoltzmann/mrtLattices.h"
#include "latticeBoltzmann/mrtLattices.hh"
#include "multiPhysics/shanChenLattices2D.h"

namespace plb {

/* *************** Exotic lattices *********************************** */

template class LatticeBoxProcessingFunctional2D<FLOAT_T, descriptors::ShanChenD2Q9Descriptor>;
template class LatticeBoxProcessingFunctional2D<FLOAT_T, descriptors::ForcedShanChenD2Q9Descriptor>;
template class LatticeBoxProcessingFunctional2D<FLOAT_T, descriptors::MRTD2Q9Descriptor>;

template class BoxProcessingFunctional2D_L<FLOAT_T, descriptors::ShanChenD2Q9Descriptor>;
template class BoxProcessingFunctional2D_L<FLOAT_T, descriptors::ForcedShanChenD2Q9Descriptor>;
template class BoxProcessingFunctional2D_L<FLOAT_T, descriptors::MRTD2Q9Descriptor>;

/* *************** Boxed Data Processor functionals ****************** */

template class BoxProcessingFunctional2D_L<FLOAT_T, descriptors::DESCRIPTOR_2D>;
template class MaskedBoxProcessingFunctional2D_L<FLOAT_T, descriptors::DESCRIPTOR_2D>;
template class BoxProcessingFunctional2D_LL< FLOAT_T, descriptors::DESCRIPTOR_2D,
                                             FLOAT_T, descriptors::DESCRIPTOR_2D >;
template class BoxProcessingFunctional2D_LN<FLOAT_T, descriptors::DESCRIPTOR_2D, FLOAT_T>;
template class BoxProcessingFunctional2D_LN<FLOAT_T, descriptors::DESCRIPTOR_2D, int>;
template class MaskedBoxProcessingFunctional2D_LN<FLOAT_T, descriptors::DESCRIPTOR_2D, FLOAT_T>;
template class MaskedBoxProcessingFunctional2D_LN<FLOAT_T, descriptors::DESCRIPTOR_2D, int>;

template class BoxProcessingFunctional2D_LS<FLOAT_T, descriptors::DESCRIPTOR_2D, FLOAT_T>;
template class BoxProcessingFunctional2D_LS<FLOAT_T, descriptors::DESCRIPTOR_2D, int>;

/* *************** Bounded Boxed Data Processor functionals ****************** */

template class BoundedBoxProcessingFunctional2D_L<FLOAT_T, descriptors::DESCRIPTOR_2D>;
template class BoundedBoxProcessingFunctional2D_LL< FLOAT_T, descriptors::DESCRIPTOR_2D,
                                                    FLOAT_T, descriptors::DESCRIPTOR_2D>;
template class BoundedBoxProcessingFunctional2D_LN< FLOAT_T, descriptors::DESCRIPTOR_2D,
                                                    FLOAT_T >;
template class BoundedMaskedBoxProcessingFunctional2D_LN< FLOAT_T, descriptors::DESCRIPTOR_2D,
                                                          FLOAT_T >;
template class BoundedLatticeBoxProcessingFunctional2D <
    FLOAT_T, descriptors::DESCRIPTOR_2D >;

}  // namespace plb

#endif  // COMPILE_2D
