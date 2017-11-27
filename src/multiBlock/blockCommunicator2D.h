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
 * Block Communicator -- Abstract base class.
 */
#ifndef BLOCK_COMMUNICATOR_2D_H
#define BLOCK_COMMUNICATOR_2D_H

#include "core/globalDefs.h"
#include <vector>

namespace plb {

class MultiBlock2D;
class MultiBlockManagement2D;
class Overlap2D;

struct BlockCommunicator2D {
    virtual ~BlockCommunicator2D() { }
    virtual BlockCommunicator2D* clone() const =0;
    /// Fill the overlaps (the "envelopes") with data from the corresponding bulks.
    /** The variable whichData specifies which type of content (static/dynamic/full dynamics object)
     *  is being transmitted.
     **/
    virtual void duplicateOverlaps(MultiBlock2D& multiBlock, modif::ModifT whichData) const =0;
    /// Transmit data between two multi-blocks, according to a user-defined pattern.
    /** The variable whichData specifies which type of content (static/dynamic/full dynamics object)
     *  is being transmitted.
     **/
    virtual void communicate( std::vector<Overlap2D> const& overlaps,
                              MultiBlock2D const& originMultiBlock,
                              MultiBlock2D& destinationMultiBlock,
                              modif::ModifT whichData ) const =0;
    virtual void signalPeriodicity() const =0;
};

}  // namespace plb

#endif  // BLOCK_COMMUNICATOR_2D_H
