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

#ifndef PARTICLE_NON_LOCAL_TRANSFER_3D_HH
#define PARTICLE_NON_LOCAL_TRANSFER_3D_HH

#include "core/globalDefs.h"
#include "particles/particleNonLocalTransfer3D.h"
#include "particles/particleProcessingFunctional3D.h"
#include "multiBlock/nonLocalTransfer3D.h"
#include <vector>

namespace plb {

template<class ParticleFieldT>
void copy (                                                                   
        MultiParticleField3D<ParticleFieldT> const& from, Box3D const& fromDomain,
        MultiParticleField3D<ParticleFieldT>& to, Box3D const& toDomain )
{       
    Box3D fromDomain_(fromDomain);                                            
    Box3D toDomain_(toDomain);
    adjustEqualSize(fromDomain_, toDomain_);                                  
    std::vector<Overlap3D> dataTransfer = copyDomainDataTransfer (            
        from.getMultiBlockManagement().getSparseBlockStructure(), fromDomain_,
        to.getMultiBlockManagement().getSparseBlockStructure(), toDomain_ );
    to.getBlockCommunicator().communicate (                                   
        dataTransfer, from, to, modif::dynamicVariables );
    to.getBlockCommunicator().duplicateOverlaps(to, modif::dynamicVariables);             
}   

template<typename T, template<typename U> class Descriptor, class ParticleFieldT>
void gatherParticles( MultiParticleField3D<ParticleFieldT>& particleField,
                      std::vector<Particle3D<T,Descriptor>*>& particles, Box3D domain )
{
    SparseBlockStructure3D blockStructure(domain);
    blockStructure.addBlock(domain, 0);
    plint envelopeWidth=1;
    MultiBlockManagement3D serialMultiBlockManagement (
            blockStructure, new OneToOneThreadAttribution, envelopeWidth );

    MultiParticleField3D<ParticleFieldT> multiSerialParticles (
            serialMultiBlockManagement,
            defaultMultiBlockPolicy3D().getCombinedStatistics() );

    copy ( particleField, domain, multiSerialParticles, domain );

    particles.clear();
    if (global::mpi().isMainProcessor()) {
        ParticleField3D<T,Descriptor>& atomicSerialParticles =
            dynamic_cast<ParticleField3D<T,Descriptor>&>(multiSerialParticles.getComponent(0));

        SmartBulk3D oneBlockBulk(serialMultiBlockManagement, 0);
        std::vector<Particle3D<T,Descriptor>*> found;
        atomicSerialParticles.findParticles(oneBlockBulk.toLocal(domain), found);
        for (pluint i=0; i<found.size(); ++i) {
            particles.push_back(found[i]->clone());
        }
    }
}

template<typename T, template<typename U> class Descriptor, class ParticleFieldT>
void injectParticlesAtMainProc( std::vector<Particle3D<T,Descriptor>*>& particles,
                                MultiParticleField3D<ParticleFieldT>& particleField, Box3D domain )
{
    SparseBlockStructure3D blockStructure(domain);
    blockStructure.addBlock(domain, 0);
    plint envelopeWidth=1;
    MultiBlockManagement3D serialMultiBlockManagement (
            blockStructure, new OneToOneThreadAttribution, envelopeWidth );

    MultiParticleField3D<ParticleFieldT> multiSerialParticles (
            serialMultiBlockManagement,
            defaultMultiBlockPolicy3D().getCombinedStatistics() );

    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(&multiSerialParticles);
    applyProcessingFunctional (
            new InjectParticlesFunctional3D<T,Descriptor>(particles),
            domain, particleArg );

    copy ( multiSerialParticles, domain, particleField, domain );
}


} // namespace plb

#endif  // PARTICLE_NON_LOCAL_TRANSFER_3D_HH
