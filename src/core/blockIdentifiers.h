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


#ifndef BLOCK_IDENTIFIERS_H
#define BLOCK_IDENTIFIERS_H

#include "core/globalDefs.h"

namespace plb {

namespace identifiers {

    enum BlockId {
        UndefinedId           = 0,

        IntScalarFieldId      = 10000,
        FloatScalarFieldId    = 10001,
        DoubleScalarFieldId   = 10002,

        IntNTensorFieldId     = 20000,
        FloatNTensorFieldId   = 20001,
        DoubleNTensorFieldId  = 20002,

        IntTensorField2Id     = 30020,
        FloatTensorField2Id   = 30021,
        DoubleTensorField2Id  = 30022,

        IntTensorField3Id     = 30030,
        FloatTensorField3Id   = 30031,
        DoubleTensorField3Id  = 30032,

        IntTensorField4Id     = 30040,
        FloatTensorField4Id   = 30041,
        DoubleTensorField4Id  = 30042,

        IntTensorField6Id     = 30060,
        FloatTensorField6Id   = 30061,
        DoubleTensorField6Id  = 30062,

        IntTensorField9Id     = 30090,
        FloatTensorField9Id   = 30091,
        DoubleTensorField9Id  = 30092,

        IntD2Q5BlockId        = 42040,
        FloatD2Q5BlockId      = 42041,
        DoubleD2Q5BlockId     = 42042,

        IntD2Q9BlockId        = 42090,
        FloatD2Q9BlockId      = 42091,
        DoubleD2Q9BlockId     = 42092,

        IntD3Q7BlockId        = 43070,
        FloatD3Q7BlockId      = 43071,
        DoubleD3Q7BlockId     = 43072,

        IntD3Q13BlockId       = 43130,
        FloatD3Q13BlockId     = 43131,
        DoubleD3Q13BlockId    = 43132,

        IntD3Q15BlockId       = 43150,
        FloatD3Q15BlockId     = 43151,
        DoubleD3Q15BlockId    = 43152,

        IntD3Q19BlockId       = 43190,
        FloatD3Q19BlockId     = 43191,
        DoubleD3Q19BlockId    = 43192,

        IntD3Q27BlockId       = 43270,
        FloatD3Q27BlockId     = 43271,
        DoubleD3Q27BlockId    = 43272,

        IntD2Q5WithForceBlockId    = 52040,
        FloatD2Q5WithForceBlockId  = 52041,
        DoubleD2Q5WithForceBlockId = 52042,

        IntD2Q9WithForceBlockId    = 52090,
        FloatD2Q9WithForceBlockId  = 52091,
        DoubleD2Q9WithForceBlockId = 52092,

        IntD3Q7WithForceBlockId    = 53070,
        FloatD3Q7WithForceBlockId  = 53071,
        DoubleD3Q7WithForceBlockId = 53072,

        IntD3Q13WithForceBlockId    = 53130,
        FloatD3Q13WithForceBlockId  = 53131,
        DoubleD3Q13WithForceBlockId = 53132,

        IntD3Q15WithForceBlockId    = 53150,
        FloatD3Q15WithForceBlockId  = 53151,
        DoubleD3Q15WithForceBlockId = 53152,

        IntD3Q19WithForceBlockId    = 53190,
        FloatD3Q19WithForceBlockId  = 53191,
        DoubleD3Q19WithForceBlockId = 53192,

        IntD3Q27WithForceBlockId    = 53270,
        FloatD3Q27WithForceBlockId  = 53271,
        DoubleD3Q27WithForceBlockId = 53272,

        ContainerId                 = 60000,

        ParticleId                  = 70000,
    };

    template<typename T>
    BlockId getScalarId() {
        return UndefinedId;
    }

    template<typename T>
    BlockId getNTensorId() {
        return UndefinedId;
    }

    template<typename T, plint n>
    BlockId getTensorId() {
        return UndefinedId;
    }

    template<typename T, template<typename U> class Descriptor>
    BlockId getLatticeId() {
        return UndefinedId;
    }

    inline BlockId getContainerId() {
        return ContainerId;
    }

    inline BlockId getParticleId() {
        return ParticleId;
    }

}  // namespace identifiers

}  // namespace plb

#endif  // BLOCK_IDENTIFIERS_H
