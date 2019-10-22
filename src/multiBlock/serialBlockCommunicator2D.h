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
 * Serial version of the 2D block communicator -- header file.
 */
#ifndef SERIAL_BLOCK_COMMUNICATOR_2D_H
#define SERIAL_BLOCK_COMMUNICATOR_2D_H

#include "core/globalDefs.h"
#include "multiBlock/blockCommunicator2D.h"
#include "multiBlock/multiBlockManagement2D.h"

namespace plb {

class MultiBlockManagement2D;

class SerialBlockCommunicator2D : public BlockCommunicator2D {
public:
    SerialBlockCommunicator2D();
    virtual SerialBlockCommunicator2D* clone() const;
    virtual void communicate( std::vector<Overlap2D> const& overlaps,
                              MultiBlock2D const& originMultiBlock,
                              MultiBlock2D& destinationMultiBlock, modif::ModifT whichData ) const;
    virtual void duplicateOverlaps(MultiBlock2D& multiBlock, modif::ModifT whichData) const;
    virtual void signalPeriodicity() const;
private:
    void copyOverlap( Overlap2D const& overlap,
                      MultiBlock2D const& fromMultiBlock,
                      MultiBlock2D& toMultiBlock, modif::ModifT whichData ) const;
};

}  // namespace plb

#endif  // SERIAL_BLOCK_COMMUNICATOR_2D_H
