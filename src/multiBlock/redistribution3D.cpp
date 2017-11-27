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
 * Utilities for 3D multi data distributions -- header file.
 */

#include "core/globalDefs.h"
#include "multiBlock/redistribution3D.h"
#include <cstdlib>

namespace plb {

RandomRedistribute3D::RandomRedistribute3D(pluint rseed_)
    : rseed(rseed_)
{ }

MultiBlockManagement3D RandomRedistribute3D::redistribute (
        MultiBlockManagement3D const& original ) const
{
    ThreadAttribution const& originalAttribution = original.getThreadAttribution();
    SparseBlockStructure3D const& originalSparseBlock = original.getSparseBlockStructure();

    std::map<plint,Box3D> const& domains = originalSparseBlock.getBulks();
    std::vector<std::pair<plint,plint> > blockToProc(domains.size());
    std::map<plint,Box3D>::const_iterator it = domains.begin();
    for (plint pos=0; it != domains.end(); ++it, ++pos) {
        plint blockId = it->first;
        plint procId = originalAttribution.getMpiProcess(blockId);
        blockToProc[pos] = std::pair<plint,plint>(blockId,procId);
    }

    srand(rseed);
    plint numExchange = (plint)blockToProc.size()*10;
    for (plint iExch=0; iExch<numExchange; ++iExch) {
        plint id1 = rand() % blockToProc.size();
        plint id2 = rand() % blockToProc.size();
        std::swap(blockToProc[id1].second, blockToProc[id2].second);
    }

    ExplicitThreadAttribution* newAttribution = new ExplicitThreadAttribution;
    for (pluint i=0; i<blockToProc.size(); ++i) {
        newAttribution->addBlock(blockToProc[i].first, blockToProc[i].second);
    }

    return MultiBlockManagement3D (
            originalSparseBlock, newAttribution,
            original.getEnvelopeWidth(), original.getRefinementLevel() );
}

}  // namespace plb

