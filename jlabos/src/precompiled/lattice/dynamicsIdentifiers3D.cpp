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

#include "core/dynamicsIdentifiers.h"
#include "core/dynamicsIdentifiers.hh"
#include "latticeBoltzmann/nearestNeighborLattices3D.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.hh"
#include "latticeBoltzmann/mrtLattices.h"
#include "latticeBoltzmann/mrtLattices.hh"
#include "multiPhysics/shanChenLattices3D.h"

namespace plb {

    namespace meta {

    template class DynamicsRegistration<FLOAT_T, descriptors::DESCRIPTOR_3D>;
    template DynamicsRegistration<FLOAT_T, descriptors::DESCRIPTOR_3D>&
        dynamicsRegistration<FLOAT_T, descriptors::DESCRIPTOR_3D>();
    template std::string
        constructIdNameChain<FLOAT_T, descriptors::DESCRIPTOR_3D>
            (std::vector<int> const& ids, std::string separator);

    template
        void createIdIndirection<FLOAT_T, descriptors::DESCRIPTOR_3D> (
                  std::map<int,std::string> const& foreignIdToName,
                  std::map<int,int>& idIndirect );


#if NUMBIT_3D == 19
    template class DynamicsRegistration<FLOAT_T, descriptors::ShanChenD3Q19Descriptor>;
    template DynamicsRegistration<FLOAT_T, descriptors::ShanChenD3Q19Descriptor>&
        dynamicsRegistration<FLOAT_T, descriptors::ShanChenD3Q19Descriptor>();

    template class DynamicsRegistration<FLOAT_T, descriptors::ForcedShanChenD3Q19Descriptor>;
    template DynamicsRegistration<FLOAT_T, descriptors::ForcedShanChenD3Q19Descriptor>&
        dynamicsRegistration<FLOAT_T, descriptors::ForcedShanChenD3Q19Descriptor>();

    template class DynamicsRegistration<FLOAT_T, descriptors::MRTD3Q19Descriptor>;
    template DynamicsRegistration<FLOAT_T, descriptors::MRTD3Q19Descriptor>&
        dynamicsRegistration<FLOAT_T, descriptors::MRTD3Q19Descriptor>();
#endif

    }  // namespace meta

}

#endif  // COMPILE_3D
