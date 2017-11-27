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
 * Copy 3D multiblocks on a new parallel distribution -- compiled code.
 */

#include "core/globalDefs.h"
#include "multiBlock/multiBlockGenerator3D.h"
#include "multiBlock/sparseBlockStructure3D.h"
#include "multiBlock/localMultiBlockInfo3D.h"
#include "multiBlock/nonLocalTransfer3D.h"
#include "multiBlock/defaultMultiBlockPolicy3D.h"
#include "multiBlock/multiBlockOperations3D.h"
#include "dataProcessors/dataAnalysisWrapper3D.h"
#include "dataProcessors/ntensorAnalysisWrapper3D.h"

namespace plb {

std::auto_ptr<MultiContainerBlock3D> generateMultiContainerBlock (
        MultiBlock3D& multiBlock, plint envelopeWidth )
{
    MultiBlockManagement3D sparseBlockManagement(multiBlock.getMultiBlockManagement());
    MultiContainerBlock3D* block = new MultiContainerBlock3D (
            MultiBlockManagement3D (
                sparseBlockManagement.getSparseBlockStructure(),
                sparseBlockManagement.getThreadAttribution().clone(),
                envelopeWidth, sparseBlockManagement.getRefinementLevel() ),
            defaultMultiBlockPolicy3D().getCombinedStatistics() );

    block->periodicity().toggle(0, multiBlock.periodicity().get(0));
    block->periodicity().toggle(1, multiBlock.periodicity().get(1));
    block->periodicity().toggle(2, multiBlock.periodicity().get(2));

    return std::auto_ptr<MultiContainerBlock3D>(block);
}

MultiContainerBlock3D* createMultiContainerBlock3D (
        MultiBlockManagement3D const& management,
        PeriodicitySwitch3D const& periodicity,
        plint envelopeWidth, plint gridLevel )
{
    MultiContainerBlock3D* block = new MultiContainerBlock3D (
            MultiBlockManagement3D (
                management.getSparseBlockStructure(),
                management.getThreadAttribution().clone(),
                envelopeWidth, gridLevel ),
            defaultMultiBlockPolicy3D().getCombinedStatistics() );

    block->periodicity().toggle(0, periodicity.get(0));
    block->periodicity().toggle(1, periodicity.get(1));
    block->periodicity().toggle(2, periodicity.get(2));

    return block;
}

MultiContainerBlock3D* createMultiContainerBlock3D (
        MultiBlockManagement3D const& management,
        plint envelopeWidth, plint gridLevel )
{
    MultiContainerBlock3D* block = new MultiContainerBlock3D (
            MultiBlockManagement3D (
                management.getSparseBlockStructure(),
                management.getThreadAttribution().clone(),
                envelopeWidth, gridLevel ),
            defaultMultiBlockPolicy3D().getCombinedStatistics() );

    return block;
}

}  // namespace plb

