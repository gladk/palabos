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

#include "core/blockIdentifiers.h"
#include "latticeBoltzmann/nearestNeighborLattices2D.h"
#include "latticeBoltzmann/nearestNeighborLattices2D.hh"
#include "latticeBoltzmann/nearestNeighborLattices3D.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.hh"

namespace plb {

namespace identifiers {

    template<>
    BlockId getScalarId<int>() {
        return IntScalarFieldId;
    }

    template<>
    BlockId getScalarId<float>() {
        return FloatScalarFieldId;
    }

    template<>
    BlockId getScalarId<double>() {
        return DoubleScalarFieldId;
    }


    template<>
    BlockId getNTensorId<int>() {
        return IntNTensorFieldId;
    }

    template<>
    BlockId getNTensorId<float>() {
        return FloatNTensorFieldId;
    }

    template<>
    BlockId getNTensorId<double>() {
        return DoubleNTensorFieldId;
    }


    template<>
    BlockId getTensorId<int,2>() {
        return IntTensorField2Id;
    }

    template<>
    BlockId getTensorId<float,2>() {
        return FloatTensorField2Id;
    }

    template<>
    BlockId getTensorId<double,2>() {
        return DoubleTensorField2Id;
    }

    template<>
    BlockId getTensorId<int,3>() {
        return IntTensorField3Id;
    }

    template<>
    BlockId getTensorId<float,3>() {
        return FloatTensorField3Id;
    }

    template<>
    BlockId getTensorId<double,3>() {
        return DoubleTensorField3Id;
    }

    template<>
    BlockId getTensorId<int,4>() {
        return IntTensorField4Id;
    }

    template<>
    BlockId getTensorId<float,4>() {
        return FloatTensorField4Id;
    }

    template<>
    BlockId getTensorId<double,4>() {
        return DoubleTensorField4Id;
    }

    template<>
    BlockId getTensorId<int,6>() {
        return IntTensorField6Id;
    }

    template<>
    BlockId getTensorId<float,6>() {
        return FloatTensorField6Id;
    }

    template<>
    BlockId getTensorId<double,6>() {
        return DoubleTensorField6Id;
    }

    template<>
    BlockId getTensorId<int,9>() {
        return IntTensorField9Id;
    }

    template<>
    BlockId getTensorId<float,9>() {
        return FloatTensorField9Id;
    }

    template<>
    BlockId getTensorId<double,9>() {
        return DoubleTensorField9Id;
    }

    template<>
    BlockId getLatticeId<int,descriptors::D2Q9Descriptor>() {
        return IntD2Q9BlockId;
    }

    template<>
    BlockId getLatticeId<float,descriptors::D2Q9Descriptor>() {
        return FloatD2Q9BlockId;
    }

    template<>
    BlockId getLatticeId<double,descriptors::D2Q9Descriptor>() {
        return DoubleD2Q9BlockId;
    }

    template<>
    BlockId getLatticeId<int,descriptors::D3Q13Descriptor>() {
        return IntD3Q13BlockId;
    }

    template<>
    BlockId getLatticeId<float,descriptors::D3Q13Descriptor>() {
        return FloatD3Q13BlockId;
    }

    template<>
    BlockId getLatticeId<double,descriptors::D3Q13Descriptor>() {
        return DoubleD3Q13BlockId;
    }

    template<>
    BlockId getLatticeId<int,descriptors::D3Q15Descriptor>() {
        return IntD3Q15BlockId;
    }

    template<>
    BlockId getLatticeId<float,descriptors::D3Q15Descriptor>() {
        return FloatD3Q15BlockId;
    }

    template<>
    BlockId getLatticeId<double,descriptors::D3Q15Descriptor>() {
        return DoubleD3Q15BlockId;
    }

    template<>
    BlockId getLatticeId<int,descriptors::D3Q19Descriptor>() {
        return IntD3Q19BlockId;
    }

    template<>
    BlockId getLatticeId<float,descriptors::D3Q19Descriptor>() {
        return FloatD3Q19BlockId;
    }

    template<>
    BlockId getLatticeId<double,descriptors::D3Q19Descriptor>() {
        return DoubleD3Q19BlockId;
    }

    template<>
    BlockId getLatticeId<int,descriptors::D3Q27Descriptor>() {
        return IntD3Q27BlockId;
    }

    template<>
    BlockId getLatticeId<float,descriptors::D3Q27Descriptor>() {
        return FloatD3Q27BlockId;
    }

    template<>
    BlockId getLatticeId<double,descriptors::D3Q27Descriptor>() {
        return DoubleD3Q27BlockId;
    }

    template<>
    BlockId getLatticeId<int,descriptors::ForcedD2Q9Descriptor>() {
        return IntD2Q9WithForceBlockId;
    }

    template<>
    BlockId getLatticeId<float,descriptors::ForcedD2Q9Descriptor>() {
        return FloatD2Q9WithForceBlockId;
    }

    template<>
    BlockId getLatticeId<double,descriptors::ForcedD2Q9Descriptor>() {
        return DoubleD2Q9WithForceBlockId;
    }

    template<>
    BlockId getLatticeId<int,descriptors::ForcedD3Q19Descriptor>() {
        return IntD3Q19WithForceBlockId;
    }

    template<>
    BlockId getLatticeId<float,descriptors::ForcedD3Q19Descriptor>() {
        return FloatD3Q19WithForceBlockId;
    }

    template<>
    BlockId getLatticeId<double,descriptors::ForcedD3Q19Descriptor>() {
        return DoubleD3Q19WithForceBlockId;
    }

}  // namespace identifiers

}  // namespace plb
