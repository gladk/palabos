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
 * Copy 2D multiblocks on a new parallel distribution -- compiled code.
 */


#include "core/globalDefs.h"
#include "multiBlock/multiBlockGenerator2D.h"
#include "multiBlock/sparseBlockStructure2D.h"
#include "multiBlock/localMultiBlockInfo2D.h"
#include "multiBlock/nonLocalTransfer2D.h"
#include "multiBlock/defaultMultiBlockPolicy2D.h"
#include "multiBlock/multiBlockOperations2D.h"
#include "dataProcessors/dataAnalysisWrapper2D.h"
#include "dataProcessors/ntensorAnalysisWrapper2D.h"

namespace plb {

std::auto_ptr<MultiContainerBlock2D> generateMultiContainerBlock (
        MultiBlock2D& multiBlock, plint envelopeWidth )
{
    MultiBlockManagement2D sparseBlockManagement(multiBlock.getMultiBlockManagement());
    MultiContainerBlock2D* block = new MultiContainerBlock2D (
            MultiBlockManagement2D (
                sparseBlockManagement.getSparseBlockStructure(),
                sparseBlockManagement.getThreadAttribution().clone(),
                envelopeWidth, sparseBlockManagement.getRefinementLevel() ),
            defaultMultiBlockPolicy2D().getCombinedStatistics() );

    block->periodicity().toggle(0, multiBlock.periodicity().get(0));
    block->periodicity().toggle(1, multiBlock.periodicity().get(1));

    return std::auto_ptr<MultiContainerBlock2D>(block);
}

MultiContainerBlock2D* createMultiContainerBlock2D (
        MultiBlockManagement2D const& management,
        PeriodicitySwitch2D& periodicity,
        plint envelopeWidth, plint gridLevel )
{
    MultiContainerBlock2D* block = new MultiContainerBlock2D (
            MultiBlockManagement2D (
                management.getSparseBlockStructure(),
                management.getThreadAttribution().clone(),
                envelopeWidth, gridLevel ),
            defaultMultiBlockPolicy2D().getCombinedStatistics() );

    block->periodicity().toggle(0, periodicity.get(0));
    block->periodicity().toggle(1, periodicity.get(1));

    return block;
}

}  // namespace plb
