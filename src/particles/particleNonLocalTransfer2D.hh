/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2015 FlowKit Sarl
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

#ifndef PARTICLE_NON_LOCAL_TRANSFER_2D_HH
#define PARTICLE_NON_LOCAL_TRANSFER_2D_HH

#include "core/globalDefs.h"
#include "particles/particleNonLocalTransfer2D.h"
#include <vector>

namespace plb {

template<class ParticleFieldT>
void copy (                                                                   
        MultiParticleField2D<ParticleFieldT> const& from, Box2D const& fromDomain,
        MultiParticleField2D<ParticleFieldT>& to, Box2D const& toDomain )
{       
    Box2D fromDomain_(fromDomain);                                            
    Box2D toDomain_(toDomain);
    adjustEqualSize(fromDomain_, toDomain_);                                  
    std::vector<Overlap2D> dataTransfer = copyDomainDataTransfer (            
        from.getMultiBlockManagement().getSparseBlockStructure(), fromDomain_,
        to.getMultiBlockManagement().getSparseBlockStructure(), toDomain_ );
    to.getBlockCommunicator().communicate (                                   
        dataTransfer, from, to, modif::dynamicVariables );
    to.getBlockCommunicator().duplicateOverlaps(to, modif::dynamicVariables);             
}   

} // namespace plb

#endif  // PARTICLE_NON_LOCAL_TRANSFER_2D_HH
