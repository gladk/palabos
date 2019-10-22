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

#include "core/dynamics.h"
#include "core/dynamics.hh"
#include "latticeBoltzmann/nearestNeighborLattices3D.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.hh"
#include "latticeBoltzmann/mrtLattices.h"
#include "latticeBoltzmann/mrtLattices.hh"

namespace plb {

    template class Dynamics<FLOAT_T, descriptors::DESCRIPTOR_3D>;
    template class BasicBulkDynamics<FLOAT_T, descriptors::DESCRIPTOR_3D>;
    template class CompositeDynamics<FLOAT_T, descriptors::DESCRIPTOR_3D>;
    template class PreparePopulationsDynamics<FLOAT_T, descriptors::DESCRIPTOR_3D>;
    template class BulkCompositeDynamics<FLOAT_T, descriptors::DESCRIPTOR_3D>;
    template class BounceBack<FLOAT_T, descriptors::DESCRIPTOR_3D>;
    template class NoDynamics<FLOAT_T, descriptors::DESCRIPTOR_3D>;

#if NUMBIT_3D == 19
    template class BasicBulkDynamics<FLOAT_T, descriptors::MRTD3Q19Descriptor>;
#endif

    template void constructIdChain<FLOAT_T, descriptors::DESCRIPTOR_3D> (
            Dynamics<FLOAT_T, descriptors::DESCRIPTOR_3D> const& dynamics, std::vector<int>& chain );

    template Dynamics<FLOAT_T, descriptors::DESCRIPTOR_3D> const&
        getBottomMostDynamics<FLOAT_T, descriptors::DESCRIPTOR_3D> (
                Dynamics<FLOAT_T, descriptors::DESCRIPTOR_3D> const& dynamics );

    template Dynamics<FLOAT_T, descriptors::DESCRIPTOR_3D>*
        cloneAndReplaceBottomDynamics<FLOAT_T, descriptors::DESCRIPTOR_3D> (
                Dynamics<FLOAT_T, descriptors::DESCRIPTOR_3D> const& dynamics,
                Dynamics<FLOAT_T, descriptors::DESCRIPTOR_3D>* newBottom );

    template Dynamics<FLOAT_T, descriptors::DESCRIPTOR_3D>*
        cloneAndInsertAtTopDynamics<FLOAT_T, descriptors::DESCRIPTOR_3D> (
                Dynamics<FLOAT_T, descriptors::DESCRIPTOR_3D> const& dynamics,
                CompositeDynamics<FLOAT_T, descriptors::DESCRIPTOR_3D>* newTop );

    template void serialize<FLOAT_T, descriptors::DESCRIPTOR_3D> (
            Dynamics<FLOAT_T, descriptors::DESCRIPTOR_3D> const& dynamics,
            std::vector<char>& data );

    template void serialize<FLOAT_T, descriptors::DESCRIPTOR_3D> (
            std::vector<Dynamics<FLOAT_T, descriptors::DESCRIPTOR_3D>*> const& dynamics,
            std::vector<char>& data );

    template pluint unserialize<FLOAT_T, descriptors::DESCRIPTOR_3D> (
            Dynamics<FLOAT_T, descriptors::DESCRIPTOR_3D>& dynamics,
            std::vector<char> const& data, pluint serializerPos );

    template void generateAndUnserializeDynamics<FLOAT_T, descriptors::DESCRIPTOR_3D> (
            std::vector<char> const& data,
            std::vector<Dynamics<FLOAT_T, descriptors::DESCRIPTOR_3D>*>& dynamics );

}

#endif  // COMPILE_3D
