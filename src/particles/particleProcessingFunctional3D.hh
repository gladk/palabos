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

#ifndef PARTICLE_PROCESSING_FUNCTIONAL_3D_HH
#define PARTICLE_PROCESSING_FUNCTIONAL_3D_HH

#include "particles/particleProcessingFunctional3D.h"
#include "particles/particleField3D.h"
#include "dataProcessors/metaStuffFunctional3D.h"
#include "core/geometry3D.h"
#include "core/plbDebug.h"
#include "core/blockStatistics.h"
#include "core/util.h"
#include "atomicBlock/atomicBlock3D.h"
#include "atomicBlock/blockLattice3D.h"
#include <algorithm>

namespace plb {

/* ******** CountParticlesFunctional3D *********************************** */

template<typename T, template<typename U> class Descriptor>
CountParticlesFunctional3D<T,Descriptor>::CountParticlesFunctional3D()
    : numParticlesId(this->getStatistics().subscribeIntSum())
{ }

template<typename T, template<typename U> class Descriptor>
void CountParticlesFunctional3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    ParticleField3D<T,Descriptor>& particleField
        = *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    std::vector<Particle3D<T,Descriptor>*> particles;
    particleField.findParticles(domain, particles);
    this->getStatistics().gatherIntSum(numParticlesId, (plint)particles.size());
}

template<typename T, template<typename U> class Descriptor>
CountParticlesFunctional3D<T,Descriptor>* CountParticlesFunctional3D<T,Descriptor>::clone() const {
    return new CountParticlesFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void CountParticlesFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
}

template<typename T, template<typename U> class Descriptor>
plint CountParticlesFunctional3D<T,Descriptor>::getNumParticles() const {
    return this->getStatistics().getIntSum(numParticlesId);
}


/* ******** CountParticlesSelectiveFunctional3D *********************************** */

template<typename T, template<typename U> class Descriptor>
CountParticlesSelectiveFunctional3D<T,Descriptor>::CountParticlesSelectiveFunctional3D(util::SelectInt* tags_)
    : numParticlesId(this->getStatistics().subscribeIntSum()),
      tags(tags_)
{ }

template<typename T, template<typename U> class Descriptor>
CountParticlesSelectiveFunctional3D<T,Descriptor>::~CountParticlesSelectiveFunctional3D()
{
    delete tags;
}

template<typename T, template<typename U> class Descriptor>
CountParticlesSelectiveFunctional3D<T,Descriptor>::CountParticlesSelectiveFunctional3D(CountParticlesSelectiveFunctional3D<T,Descriptor> const& rhs)
    : numParticlesId(this->getStatistics().subscribeIntSum()),
      tags(rhs.tags->clone())
{ }

template<typename T, template<typename U> class Descriptor>
CountParticlesSelectiveFunctional3D<T,Descriptor>&
    CountParticlesSelectiveFunctional3D<T,Descriptor>::operator=(CountParticlesSelectiveFunctional3D<T,Descriptor> const& rhs)
{
    CountParticlesSelectiveFunctional3D<T,Descriptor>(rhs).swap(*this);
    return *this;
}

template<typename T, template<typename U> class Descriptor>
void CountParticlesSelectiveFunctional3D<T,Descriptor>::swap(CountParticlesSelectiveFunctional3D<T,Descriptor>& rhs) {
    std::swap(numParticlesId, rhs.numParticlesId);
    std::swap(tags, rhs.tags);
    PlainReductiveBoxProcessingFunctional3D::swap(rhs);
}

template<typename T, template<typename U> class Descriptor>
void CountParticlesSelectiveFunctional3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    ParticleField3D<T,Descriptor>& particleField
        = *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    std::vector<Particle3D<T,Descriptor>*> particles;
    particleField.findParticles(domain, particles);
    for (pluint iParticle=0; iParticle<particles.size(); ++iParticle) {
        if ((*tags)(particles[iParticle]->getTag())) {
	    this->getStatistics().gatherIntSum(numParticlesId, 1);
        }
    }
}

template<typename T, template<typename U> class Descriptor>
CountParticlesSelectiveFunctional3D<T,Descriptor>* CountParticlesSelectiveFunctional3D<T,Descriptor>::clone() const {
    return new CountParticlesSelectiveFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void CountParticlesSelectiveFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
}

template<typename T, template<typename U> class Descriptor>
plint CountParticlesSelectiveFunctional3D<T,Descriptor>::getNumParticles() const {
    return this->getStatistics().getIntSum(numParticlesId);
}


/* ******** AverageParticleVelocityFunctional3D *********************************** */

template<typename T, template<typename U> class Descriptor>
AverageParticleVelocityFunctional3D<T,Descriptor>::AverageParticleVelocityFunctional3D()
    : averageVelocityId (
            Array<T,3>( this->getStatistics().subscribeAverage(),
                        this->getStatistics().subscribeAverage(),
                        this->getStatistics().subscribeAverage()) )
{ }

template<typename T, template<typename U> class Descriptor>
void AverageParticleVelocityFunctional3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    ParticleField3D<T,Descriptor>& particleField = *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    std::vector<Particle3D<T,Descriptor>*> particles;
    particleField.findParticles(domain, particles);
    for (pluint iParticle=0; iParticle<particles.size(); ++iParticle) {
        Array<T,3> velocity;
        particles[iParticle]->getVelocity(velocity);
        this->getStatistics().gatherAverage(averageVelocityId[0], velocity[0]);
        this->getStatistics().gatherAverage(averageVelocityId[1], velocity[1]);
        this->getStatistics().gatherAverage(averageVelocityId[2], velocity[2]);
        this->getStatistics().incrementStats();
    }
}

template<typename T, template<typename U> class Descriptor>
AverageParticleVelocityFunctional3D<T,Descriptor>* AverageParticleVelocityFunctional3D<T,Descriptor>::clone() const {
    return new AverageParticleVelocityFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void AverageParticleVelocityFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
}

template<typename T, template<typename U> class Descriptor>
Array<T,3> AverageParticleVelocityFunctional3D<T,Descriptor>::getAverageParticleVelocity() const {
    return Array<T,3> (
        this->getStatistics().getAverage(averageVelocityId[0]),
        this->getStatistics().getAverage(averageVelocityId[1]),
        this->getStatistics().getAverage(averageVelocityId[2]) );
}


/* ******** CopySelectParticles3D *********************************** */

template<typename T, template<typename U> class Descriptor>
CopySelectParticles3D<T,Descriptor>::CopySelectParticles3D(util::SelectInt* tags_)
    : tags(tags_)
{ }

template<typename T, template<typename U> class Descriptor>
CopySelectParticles3D<T,Descriptor>::~CopySelectParticles3D()
{
    delete tags;
}

template<typename T, template<typename U> class Descriptor>
CopySelectParticles3D<T,Descriptor>::CopySelectParticles3D(CopySelectParticles3D<T,Descriptor> const& rhs)
    : tags(rhs.tags->clone())
{ }

template<typename T, template<typename U> class Descriptor>
CopySelectParticles3D<T,Descriptor>&
    CopySelectParticles3D<T,Descriptor>::operator=(CopySelectParticles3D<T,Descriptor> const& rhs)
{
    CopySelectParticles3D<T,Descriptor>(rhs).swap(*this);
    return *this;
}

template<typename T, template<typename U> class Descriptor>
void CopySelectParticles3D<T,Descriptor>::swap(CopySelectParticles3D<T,Descriptor>& rhs) {
    std::swap(tags, rhs.tags);
}

template<typename T, template<typename U> class Descriptor>
void CopySelectParticles3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==2 );
    ParticleField3D<T,Descriptor>& fromParticles = *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    ParticleField3D<T,Descriptor>& toParticles = *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[1]);

    std::vector<Particle3D<T,Descriptor>*> particles;
    fromParticles.findParticles(domain, particles);
    for (pluint iParticle=0; iParticle<particles.size(); ++iParticle) {
        if ((*tags)(particles[iParticle]->getTag())) {
            toParticles.addParticle(domain, particles[iParticle]->clone());
        }
    }
}

template<typename T, template<typename U> class Descriptor>
CopySelectParticles3D<T,Descriptor>* CopySelectParticles3D<T,Descriptor>::clone() const {
    return new CopySelectParticles3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void CopySelectParticles3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing;          // From Particle field.
    modified[1] = modif::dynamicVariables; // To Particle field.
}


/* ******** InjectParticlesFunctional3D *********************************** */

template<typename T, template<typename U> class Descriptor>
InjectParticlesFunctional3D<T,Descriptor>::InjectParticlesFunctional3D(std::vector<Particle3D<T,Descriptor>*>& particles_)
    : particles(particles_)
{ }

template<typename T, template<typename U> class Descriptor>
void InjectParticlesFunctional3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    ParticleField3D<T,Descriptor>& particleField = *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    for (pluint iParticle=0; iParticle<particles.size(); ++iParticle) {
        // The function addParticle tests automatically if the particle corresponds to the domain
        //   and needs to be added. In either case, the particle is consumed by the function.
        particleField.addParticle(domain, particles[iParticle]->clone());
    }
}

template<typename T, template<typename U> class Descriptor>
InjectParticlesFunctional3D<T,Descriptor>* InjectParticlesFunctional3D<T,Descriptor>::clone() const {
    return new InjectParticlesFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void InjectParticlesFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables;  // Particle field.
}


/* ******** InjectRandomParticlesFunctional3D *********************************** */

template<typename T, template<typename U> class Descriptor>
InjectRandomParticlesFunctional3D<T,Descriptor>::InjectRandomParticlesFunctional3D (
        Particle3D<T,Descriptor>* particleTemplate_, T probabilityPerCell_ )
    : particleTemplate(particleTemplate_),
      probabilityPerCell(probabilityPerCell_)
{ }

template<typename T, template<typename U> class Descriptor>
InjectRandomParticlesFunctional3D<T,Descriptor>::InjectRandomParticlesFunctional3D (
        InjectRandomParticlesFunctional3D<T,Descriptor> const& rhs)
    : particleTemplate(rhs.particleTemplate->clone()),
      probabilityPerCell(rhs.probabilityPerCell)
{ }

template<typename T, template<typename U> class Descriptor>
InjectRandomParticlesFunctional3D<T,Descriptor>&
    InjectRandomParticlesFunctional3D<T,Descriptor>::operator= (
        InjectRandomParticlesFunctional3D<T,Descriptor> const& rhs )
{
    InjectRandomParticlesFunctional3D<T,Descriptor>(rhs).swap(*this);
    return *this;
}

template<typename T, template<typename U> class Descriptor>
void InjectRandomParticlesFunctional3D<T,Descriptor>::swap (
        InjectRandomParticlesFunctional3D<T,Descriptor>& rhs )
{
    std::swap(particleTemplate, rhs.particleTemplate);
    std::swap(probabilityPerCell, rhs.probabilityPerCell);
}

template<typename T, template<typename U> class Descriptor>
InjectRandomParticlesFunctional3D<T,Descriptor>::~InjectRandomParticlesFunctional3D()
{
    delete particleTemplate;
}

template<typename T, template<typename U> class Descriptor>
void InjectRandomParticlesFunctional3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    ParticleField3D<T,Descriptor>& particleField = *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                T randNumber = (T)rand() / (T)RAND_MAX;
                if(randNumber<probabilityPerCell) {
                    T randX = (T)rand() / (T)RAND_MAX - (T)1;
                    T randY = (T)rand() / (T)RAND_MAX - (T)1;
                    T randZ = (T)rand() / (T)RAND_MAX - (T)1;
                    Particle3D<T,Descriptor>* newparticle = particleTemplate->clone();
                    newparticle->getPosition() =
                            Array<T,3> (
                                    particleField.getLocation().x + iX + randX,
                                    particleField.getLocation().y + iY + randY,
                                    particleField.getLocation().z + iZ + randZ );
                    particleField.addParticle(domain, newparticle);
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
InjectRandomParticlesFunctional3D<T,Descriptor>* InjectRandomParticlesFunctional3D<T,Descriptor>::clone() const {
    return new InjectRandomParticlesFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void InjectRandomParticlesFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables;  // Particle field.
}


/* ******** AnalyticalInjectRandomParticlesFunctional3D *********************************** */

template<typename T, template<typename U> class Descriptor, class DomainFunctional>
AnalyticalInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional>::
    AnalyticalInjectRandomParticlesFunctional3D (
         Particle3D<T,Descriptor>* particleTemplate_, T probabilityPerCell_, DomainFunctional functional_ )
    : particleTemplate(particleTemplate_),
      probabilityPerCell(probabilityPerCell_),
      functional(functional_)
{ }

template<typename T, template<typename U> class Descriptor, class DomainFunctional>
AnalyticalInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional>::AnalyticalInjectRandomParticlesFunctional3D (
        AnalyticalInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional> const& rhs)
    : particleTemplate(rhs.particleTemplate->clone()),
      probabilityPerCell(rhs.probabilityPerCell),
      functional(rhs.functional)
{ }

template<typename T, template<typename U> class Descriptor, class DomainFunctional>
AnalyticalInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional>&
    AnalyticalInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional>::operator= (
        AnalyticalInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional> const& rhs )
{
    AnalyticalInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional>(rhs).swap(*this);
    return *this;
}

template<typename T, template<typename U> class Descriptor, class DomainFunctional>
void AnalyticalInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional>::swap (
        AnalyticalInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional>& rhs )
{
    std::swap(particleTemplate, rhs.particleTemplate);
    std::swap(probabilityPerCell, rhs.probabilityPerCell);
    std::swap(functional, rhs.functional);
}

template<typename T, template<typename U> class Descriptor, class DomainFunctional>
AnalyticalInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional>::~AnalyticalInjectRandomParticlesFunctional3D()
{
    delete particleTemplate;
}

template<typename T, template<typename U> class Descriptor, class DomainFunctional>
void AnalyticalInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    ParticleField3D<T,Descriptor>& particleField = *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                T randNumber = (T)rand() / (T)RAND_MAX;
                if(randNumber<probabilityPerCell) {
                    T randX = (T)rand() / (T)RAND_MAX - (T)1;
                    T randY = (T)rand() / (T)RAND_MAX - (T)1;
                    T randZ = (T)rand() / (T)RAND_MAX - (T)1;
                    Array<T,3> position (
                            particleField.getLocation().x + iX + randX,
                            particleField.getLocation().y + iY + randY,
                            particleField.getLocation().z + iZ + randZ );
                    if (functional(position)) {
                        Particle3D<T,Descriptor>* newparticle = particleTemplate->clone();
                        newparticle->getPosition() = position;
                        particleField.addParticle(domain, newparticle);
                    }
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor, class DomainFunctional>
AnalyticalInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional>*
    AnalyticalInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional>::clone() const {
    return new AnalyticalInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional>(*this);
}

template<typename T, template<typename U> class Descriptor, class DomainFunctional>
void AnalyticalInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables;  // Particle field.
}


/* ******** MaskedInjectRandomParticlesFunctional3D *********************************** */

template<typename T, template<typename U> class Descriptor, class DomainFunctional>
MaskedInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional>::MaskedInjectRandomParticlesFunctional3D (
        Particle3D<T,Descriptor>* particleTemplate_, T probabilityPerCell_, DomainFunctional functional_, int flag_ )
    : particleTemplate(particleTemplate_),
      probabilityPerCell(probabilityPerCell_),
      functional(functional_),
      flag(flag_)
{ }

template<typename T, template<typename U> class Descriptor, class DomainFunctional>
MaskedInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional>::MaskedInjectRandomParticlesFunctional3D (
        MaskedInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional> const& rhs)
    : particleTemplate(rhs.particleTemplate->clone()),
      probabilityPerCell(rhs.probabilityPerCell),
      functional(rhs.functional),
      flag(rhs.flag)
{ }

template<typename T, template<typename U> class Descriptor, class DomainFunctional>
MaskedInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional>&
    MaskedInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional>::operator= (
        MaskedInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional> const& rhs )
{
    MaskedInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional>(rhs).swap(*this);
    return *this;
}

template<typename T, template<typename U> class Descriptor, class DomainFunctional>
void MaskedInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional>::swap (
        MaskedInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional>& rhs )
{
    std::swap(particleTemplate, rhs.particleTemplate);
    std::swap(probabilityPerCell, rhs.probabilityPerCell);
    std::swap(functional, rhs.functional);
    std::swap(flag, rhs.flag);
}

template<typename T, template<typename U> class Descriptor, class DomainFunctional>
MaskedInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional>::~MaskedInjectRandomParticlesFunctional3D()
{
    delete particleTemplate;
}

template<typename T, template<typename U> class Descriptor, class DomainFunctional>
void MaskedInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==2 );
    ParticleField3D<T,Descriptor>& particleField = *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    ScalarField3D<int>& flagMatrix = *dynamic_cast<ScalarField3D<int>*>(blocks[1]);
    Dot3D offset = computeRelativeDisplacement(particleField, flagMatrix);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                T randNumber = (T)rand() / (T)RAND_MAX;
                if(randNumber<probabilityPerCell) {
                    T randX = (T)rand() / (T)RAND_MAX - (T)0.5;
                    T randY = (T)rand() / (T)RAND_MAX - (T)0.5;
                    T randZ = (T)rand() / (T)RAND_MAX - (T)0.5;
                    Array<T,3> position (
                            particleField.getLocation().x + iX + randX,
                            particleField.getLocation().y + iY + randY,
                            particleField.getLocation().z + iZ + randZ );
                    
                    if (functional(position) && flagMatrix.get(iX+offset.x,iY+offset.y,iZ+offset.z)==flag) {
                        Particle3D<T,Descriptor>* newparticle = particleTemplate->clone();
                        newparticle->getPosition() = position;
                        particleField.addParticle(domain, newparticle);
                    }
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor, class DomainFunctional>
MaskedInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional>*
    MaskedInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional>::clone() const
{
    return new MaskedInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional>(*this);
}

template<typename T, template<typename U> class Descriptor, class DomainFunctional>
void MaskedInjectRandomParticlesFunctional3D<T,Descriptor,DomainFunctional>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables;  // Particle field.
    modified[1] = modif::nothing;  // Mask.
}


/* ******** InjectEquallySpacedParticlesFunctional3D *********************************** */

template<typename T, template<typename U> class Descriptor>
InjectEquallySpacedParticlesFunctional3D<T,Descriptor>::InjectEquallySpacedParticlesFunctional3D (
        Particle3D<T,Descriptor>* particleTemplate_, plint nx_, plint ny_, plint nz_ )
    : particleTemplate(particleTemplate_),
      nx(nx_),
      ny(ny_),
      nz(nz_)
{
    PLB_ASSERT(nx > 0);
    PLB_ASSERT(ny > 0);
    PLB_ASSERT(nz > 0);
}

template<typename T, template<typename U> class Descriptor>
InjectEquallySpacedParticlesFunctional3D<T,Descriptor>::InjectEquallySpacedParticlesFunctional3D (
        InjectEquallySpacedParticlesFunctional3D<T,Descriptor> const& rhs)
    : particleTemplate(rhs.particleTemplate->clone()),
      nx(rhs.nx),
      ny(rhs.ny),
      nz(rhs.nz)
{ }

template<typename T, template<typename U> class Descriptor>
InjectEquallySpacedParticlesFunctional3D<T,Descriptor>&
    InjectEquallySpacedParticlesFunctional3D<T,Descriptor>::operator= (
        InjectEquallySpacedParticlesFunctional3D<T,Descriptor> const& rhs )
{
    InjectEquallySpacedParticlesFunctional3D<T,Descriptor>(rhs).swap(*this);
    return *this;
}

template<typename T, template<typename U> class Descriptor>
void InjectEquallySpacedParticlesFunctional3D<T,Descriptor>::swap (
        InjectEquallySpacedParticlesFunctional3D<T,Descriptor>& rhs )
{
    std::swap(particleTemplate, rhs.particleTemplate);
    std::swap(nx, rhs.nx);
    std::swap(ny, rhs.ny);
    std::swap(nz, rhs.nz);
}

template<typename T, template<typename U> class Descriptor>
InjectEquallySpacedParticlesFunctional3D<T,Descriptor>::~InjectEquallySpacedParticlesFunctional3D()
{
    delete particleTemplate;
}

template<typename T, template<typename U> class Descriptor>
void InjectEquallySpacedParticlesFunctional3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    ParticleField3D<T,Descriptor> *particleField = dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    PLB_ASSERT(particleField);

    Dot3D loc = particleField->getLocation();

    T dx = 1.0 / (T) nx;
    T dy = 1.0 / (T) ny;
    T dz = 1.0 / (T) nz;

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        T x0 = iX + loc.x - 0.5;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            T y0 = iY + loc.y - 0.5;
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                T z0 = iZ + loc.z - 0.5;
                for (plint i = 0; i < nx; i++) {
                    T posX = x0 + i * dx;
                    for (plint j = 0; j < ny; j++) {
                        T posY = y0 + j * dy;
                        for (plint k = 0; k < nz; k++) {
                            T posZ = z0 + k * dz;

                            Particle3D<T,Descriptor>* newparticle = particleTemplate->clone();
                            newparticle->getPosition() = Array<T,3>(posX, posY, posZ);
                            particleField->addParticle(domain, newparticle);
                        }
                    }
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
InjectEquallySpacedParticlesFunctional3D<T,Descriptor>* InjectEquallySpacedParticlesFunctional3D<T,Descriptor>::clone() const {
    return new InjectEquallySpacedParticlesFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void InjectEquallySpacedParticlesFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables;  // Particle field.
}


/* ******** MaskedInjectEquallySpacedParticlesFunctional3D *********************************** */

template<typename T, template<typename U> class Descriptor>
MaskedInjectEquallySpacedParticlesFunctional3D<T,Descriptor>::MaskedInjectEquallySpacedParticlesFunctional3D (
        Particle3D<T,Descriptor>* particleTemplate_, plint nx_, plint ny_, plint nz_, int flag_ )
    : particleTemplate(particleTemplate_),
      nx(nx_),
      ny(ny_),
      nz(nz_),
      flag(flag_)
{
    PLB_ASSERT(nx > 0);
    PLB_ASSERT(ny > 0);
    PLB_ASSERT(nz > 0);
}

template<typename T, template<typename U> class Descriptor>
MaskedInjectEquallySpacedParticlesFunctional3D<T,Descriptor>::MaskedInjectEquallySpacedParticlesFunctional3D (
        MaskedInjectEquallySpacedParticlesFunctional3D<T,Descriptor> const& rhs)
    : particleTemplate(rhs.particleTemplate->clone()),
      nx(rhs.nx),
      ny(rhs.ny),
      nz(rhs.nz),
      flag(rhs.flag)
{ }

template<typename T, template<typename U> class Descriptor>
MaskedInjectEquallySpacedParticlesFunctional3D<T,Descriptor>&
    MaskedInjectEquallySpacedParticlesFunctional3D<T,Descriptor>::operator= (
        MaskedInjectEquallySpacedParticlesFunctional3D<T,Descriptor> const& rhs )
{
    MaskedInjectEquallySpacedParticlesFunctional3D<T,Descriptor>(rhs).swap(*this);
    return *this;
}

template<typename T, template<typename U> class Descriptor>
void MaskedInjectEquallySpacedParticlesFunctional3D<T,Descriptor>::swap (
        MaskedInjectEquallySpacedParticlesFunctional3D<T,Descriptor>& rhs )
{
    std::swap(particleTemplate, rhs.particleTemplate);
    std::swap(nx, rhs.nx);
    std::swap(ny, rhs.ny);
    std::swap(nz, rhs.nz);
    std::swap(flag, rhs.flag);
}

template<typename T, template<typename U> class Descriptor>
MaskedInjectEquallySpacedParticlesFunctional3D<T,Descriptor>::~MaskedInjectEquallySpacedParticlesFunctional3D()
{
    delete particleTemplate;
}

template<typename T, template<typename U> class Descriptor>
void MaskedInjectEquallySpacedParticlesFunctional3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==2 );
    ParticleField3D<T,Descriptor> *particleField = dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    PLB_ASSERT(particleField);
    ScalarField3D<int> *flagMatrix = dynamic_cast<ScalarField3D<int>*>(blocks[1]);
    PLB_ASSERT(flagMatrix);

    Dot3D loc = particleField->getLocation();
    Dot3D offset = computeRelativeDisplacement(*particleField, *flagMatrix);

    T dx = 1.0 / (T) nx;
    T dy = 1.0 / (T) ny;
    T dz = 1.0 / (T) nz;

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        T x0 = iX + loc.x - 0.5;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            T y0 = iY + loc.y - 0.5;
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                T z0 = iZ + loc.z - 0.5;
                if (flagMatrix->get(iX+offset.x,iY+offset.y,iZ+offset.z)==flag) {
                    for (plint i = 0; i < nx; i++) {
                        T posX = x0 + i * dx;
                        for (plint j = 0; j < ny; j++) {
                            T posY = y0 + j * dy;
                            for (plint k = 0; k < nz; k++) {
                                T posZ = z0 + k * dz;

                                Particle3D<T,Descriptor>* newparticle = particleTemplate->clone();
                                newparticle->getPosition() = Array<T,3>(posX, posY, posZ);
                                particleField->addParticle(domain, newparticle);
                            }
                        }
                    }
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
MaskedInjectEquallySpacedParticlesFunctional3D<T,Descriptor>* MaskedInjectEquallySpacedParticlesFunctional3D<T,Descriptor>::clone() const {
    return new MaskedInjectEquallySpacedParticlesFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void MaskedInjectEquallySpacedParticlesFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables;  // Particle field.
    modified[1] = modif::nothing;  // Mask.
}


/* ******** AbsorbParticlesFunctional3D *********************************** */

template<typename T, template<typename U> class Descriptor>
void AbsorbParticlesFunctional3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    particleField.removeParticles(domain);
}

template<typename T, template<typename U> class Descriptor>
AbsorbParticlesFunctional3D<T,Descriptor>* AbsorbParticlesFunctional3D<T,Descriptor>::clone() const {
    return new AbsorbParticlesFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void AbsorbParticlesFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables;  // Particle field.
}


/* ******** AbsorbParticlesFunctionalSelective3D *********************************** */

template<typename T, template<typename U> class Descriptor>
AbsorbParticlesFunctionalSelective3D<T,Descriptor>::AbsorbParticlesFunctionalSelective3D(util::SelectInt* tags_)
    : tags(tags_)
{ }

template<typename T, template<typename U> class Descriptor>
AbsorbParticlesFunctionalSelective3D<T,Descriptor>::~AbsorbParticlesFunctionalSelective3D()
{
    delete tags;
}

template<typename T, template<typename U> class Descriptor>
AbsorbParticlesFunctionalSelective3D<T,Descriptor>::AbsorbParticlesFunctionalSelective3D(AbsorbParticlesFunctionalSelective3D<T,Descriptor> const& rhs)
    : tags(rhs.tags->clone())
{ }

template<typename T, template<typename U> class Descriptor>
AbsorbParticlesFunctionalSelective3D<T,Descriptor>&
    AbsorbParticlesFunctionalSelective3D<T,Descriptor>::operator=(AbsorbParticlesFunctionalSelective3D<T,Descriptor> const& rhs)
{
    AbsorbParticlesFunctionalSelective3D<T,Descriptor>(rhs).swap(*this);
    return *this;
}

template<typename T, template<typename U> class Descriptor>
void AbsorbParticlesFunctionalSelective3D<T,Descriptor>::swap(AbsorbParticlesFunctionalSelective3D<T,Descriptor>& rhs) {
    std::swap(tags, rhs.tags);
}

template<typename T, template<typename U> class Descriptor>
void AbsorbParticlesFunctionalSelective3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    std::vector<Particle3D<T,Descriptor>*> found;
    particleField.findParticles(domain, found);
    std::vector<Particle3D<T,Descriptor>*> remaining;
    for (pluint i=0; i<found.size(); ++i) {
        if(!(*tags)(found[i]->getTag())) {
            remaining.push_back(found[i]->clone());
        }
    }
    particleField.removeParticles(domain);
    for (pluint i=0; i<remaining.size(); ++i) {
        particleField.addParticle(domain, remaining[i]);
    }
}

template<typename T, template<typename U> class Descriptor>
AbsorbParticlesFunctionalSelective3D<T,Descriptor>* AbsorbParticlesFunctionalSelective3D<T,Descriptor>::clone() const {
    return new AbsorbParticlesFunctionalSelective3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void AbsorbParticlesFunctionalSelective3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables;  // Particle field.
}


/* ******** RemoveParticlesFromWall3D *********************************** */

template<typename T, template<typename U> class Descriptor>
RemoveParticlesFromWall3D<T,Descriptor>::RemoveParticlesFromWall3D(int wallFlag_)
    : wallFlag(wallFlag_)
{ }

template<typename T, template<typename U> class Descriptor>
void RemoveParticlesFromWall3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_ASSERT(blocks.size()==2);
    ParticleField3D<T,Descriptor> *particleField = dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    PLB_ASSERT(particleField);
    ScalarField3D<int> *flag = dynamic_cast<ScalarField3D<int>*>(blocks[1]);
    PLB_ASSERT(flag);

    Dot3D ofs = computeRelativeDisplacement(*particleField, *flag);

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (flag->get(iX+ofs.x, iY+ofs.y, iZ+ofs.z) == wallFlag) {
                    Box3D node(iX, iX, iY, iY, iZ, iZ);
                    particleField->removeParticles(node);
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
RemoveParticlesFromWall3D<T,Descriptor>* RemoveParticlesFromWall3D<T,Descriptor>::clone() const
{
    return new RemoveParticlesFromWall3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void RemoveParticlesFromWall3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables; // Particle field.
    modified[1] = modif::nothing;          // Flag matrix.
}


/* ******** PushParticlesAwayFromWall3D *********************************** */

template<typename T, template<typename U> class Descriptor>
PushParticlesAwayFromWall3D<T,Descriptor>::PushParticlesAwayFromWall3D(T cutOffValue_, T movingDistance_,
        int wallFlag_, int fluidFlag_)
    : cutOffValue(cutOffValue_),
      movingDistance(movingDistance_),
      wallFlag(wallFlag_),
      fluidFlag(fluidFlag_)
{ }

template<typename T, template<typename U> class Descriptor>
void PushParticlesAwayFromWall3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_ASSERT(blocks.size()==2);
    ParticleField3D<T,Descriptor> *particleField = dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    PLB_ASSERT(particleField);
    ScalarField3D<int> *flag = dynamic_cast<ScalarField3D<int>*>(blocks[1]);
    PLB_ASSERT(flag);

    Dot3D ofs = computeRelativeDisplacement(*particleField, *flag);

    // The next block is for optimization purposes.
    {
        std::vector<Particle3D<T,Descriptor>*> particles;
        particleField->findParticles(domain, particles);
        if (particles.size() == 0) {
            return;
        }
    }

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Box3D node(iX, iX, iY, iY, iZ, iZ);
                std::vector<Particle3D<T,Descriptor>*> particles;
                particleField->findParticles(node, particles);
                if (particles.size() == 0) {
                    continue;
                }

                bool hasFluidNeighbour = false;
                if (flag->get(iX+ofs.x, iY+ofs.y, iZ+ofs.z) == wallFlag) {
                    for (int i = -1; i < 2; i++) {
                        for (int j = -1; j < 2; j++) {
                            for (int k = -1; k < 2; k++) {
                                if (flag->get(iX+i+ofs.x, iY+j+ofs.y, iZ+k+ofs.z) == fluidFlag) {
                                    hasFluidNeighbour = true;
                                    goto hasFluidNeighbourLabel;
                                }
                            }
                        }
                    }
                }
hasFluidNeighbourLabel:
                if (flag->get(iX+ofs.x, iY+ofs.y, iZ+ofs.z) == wallFlag && !hasFluidNeighbour) {
                    continue;
                }

                int numWallCells = 0;
                Array<int,3> intWallNormal(0, 0, 0);
                for (int i = -1; i < 2; i++) {
                    for (int j = -1; j < 2; j++) {
                        for (int k = -1; k < 2; k++) {
                            if (flag->get(iX+i+ofs.x, iY+j+ofs.y, iZ+k+ofs.z) == wallFlag) {
                                numWallCells++;
                                intWallNormal += Array<int,3>(-i, -j, -k);
                            }
                        }
                    }
                }
                if (numWallCells == 0) {
                    continue;
                } else {
                    int norm2intWallNormal = intWallNormal[0] * intWallNormal[0] + 
                                             intWallNormal[1] * intWallNormal[1] + 
                                             intWallNormal[2] * intWallNormal[2];
                    if (norm2intWallNormal != 0) {
                        T normWallNormal = std::sqrt((T) norm2intWallNormal);
                        Array<T,3> wallNormal;
                        wallNormal[0] = (T) intWallNormal[0] / normWallNormal;
                        wallNormal[1] = (T) intWallNormal[1] / normWallNormal;
                        wallNormal[2] = (T) intWallNormal[2] / normWallNormal;

                        for (pluint iParticle=0; iParticle<particles.size(); ++iParticle) {
                            PointParticle3D<T,Descriptor> *pp =
                                dynamic_cast<PointParticle3D<T,Descriptor>*>(particles[iParticle]);
                            PLB_ASSERT(pp);
                            if (normSqr(pp->getVelocity()) <= cutOffValue) {
                                pp->getVelocity() = movingDistance * wallNormal;
                            }
                        }
                    }
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
PushParticlesAwayFromWall3D<T,Descriptor>* PushParticlesAwayFromWall3D<T,Descriptor>::clone() const
{
    return new PushParticlesAwayFromWall3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void PushParticlesAwayFromWall3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables; // Particle field.
    modified[1] = modif::nothing;          // Flag matrix.
}


/* ******** FluidToParticleCoupling3D *********************************** */

template<typename T, template<typename U> class Descriptor>
FluidToParticleCoupling3D<T,Descriptor>::FluidToParticleCoupling3D(T scaling_)
    : scaling(scaling_)
{ }

template<typename T, template<typename U> class Descriptor>
void FluidToParticleCoupling3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==2 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    BlockLattice3D<T,Descriptor>& fluid =
        *dynamic_cast<BlockLattice3D<T,Descriptor>*>(blocks[1]);
    particleField.fluidToParticleCoupling(domain, fluid, scaling);
}

template<typename T, template<typename U> class Descriptor>
FluidToParticleCoupling3D<T,Descriptor>* FluidToParticleCoupling3D<T,Descriptor>::clone() const {
    return new FluidToParticleCoupling3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void FluidToParticleCoupling3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables; // Particle field.
    modified[1] = modif::nothing;  // Fluid.
}

/* ******** VelocityToParticleCoupling3D *********************************** */

template<typename T, template<typename U> class Descriptor>
VelocityToParticleCoupling3D<T,Descriptor>::VelocityToParticleCoupling3D(T scaling_)
    : scaling(scaling_)
{ }

template<typename T, template<typename U> class Descriptor>
void VelocityToParticleCoupling3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==2 );
    ParticleField3D<T,Descriptor>* particleFieldPtr =
        dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    PLB_ASSERT( particleFieldPtr );
    ParticleField3D<T,Descriptor>& particleField = *particleFieldPtr;

    TensorField3D<T,3>* velocityPtr =
        dynamic_cast<TensorField3D<T,3>*>(blocks[1]);
    PLB_ASSERT( velocityPtr );
    TensorField3D<T,3>& velocity = *velocityPtr;

    particleField.velocityToParticleCoupling(domain, velocity, scaling);
}

template<typename T, template<typename U> class Descriptor>
VelocityToParticleCoupling3D<T,Descriptor>* VelocityToParticleCoupling3D<T,Descriptor>::clone() const {
    return new VelocityToParticleCoupling3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void VelocityToParticleCoupling3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables; // Particle field.
    modified[1] = modif::nothing;  // Fluid.
}


/* ******** RhoBarJtoParticleCoupling3D *********************************** */

template<typename T, template<typename U> class Descriptor>
RhoBarJtoParticleCoupling3D<T,Descriptor>::RhoBarJtoParticleCoupling3D (
        bool velIsJ_, T scaling_)
    : velIsJ(velIsJ_),
      scaling(scaling_)
{ }

template<typename T, template<typename U> class Descriptor>
void RhoBarJtoParticleCoupling3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==2 );
    ParticleField3D<T,Descriptor>* particleFieldPtr =
        dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    PLB_ASSERT( particleFieldPtr );
    ParticleField3D<T,Descriptor>& particleField = *particleFieldPtr;

    NTensorField3D<T>* rhoBarJptr =
        dynamic_cast<NTensorField3D<T>*>(blocks[1]);
    PLB_ASSERT( rhoBarJptr );
    NTensorField3D<T>& rhoBarJ = *rhoBarJptr;

    particleField.rhoBarJtoParticleCoupling(domain, rhoBarJ, velIsJ, scaling);
}

template<typename T, template<typename U> class Descriptor>
RhoBarJtoParticleCoupling3D<T,Descriptor>* RhoBarJtoParticleCoupling3D<T,Descriptor>::clone() const {
    return new RhoBarJtoParticleCoupling3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void RhoBarJtoParticleCoupling3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables; // Particle field.
    modified[1] = modif::nothing;  // RhoBarJ;
}


/* ******** AdvanceParticlesFunctional3D *********************************** */

template<typename T, template<typename U> class Descriptor>
AdvanceParticlesFunctional3D<T,Descriptor>::AdvanceParticlesFunctional3D (
        T cutOffValue_ )
  : cutOffValue(cutOffValue_)
{ }

template<typename T, template<typename U> class Descriptor>
void AdvanceParticlesFunctional3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    particleField.advanceParticles(domain, cutOffValue);
}

template<typename T, template<typename U> class Descriptor>
AdvanceParticlesFunctional3D<T,Descriptor>* AdvanceParticlesFunctional3D<T,Descriptor>::clone() const {
    return new AdvanceParticlesFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT AdvanceParticlesFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;  // Important: access envelope as well,
                                          // because particles are streamed from the
                                          // envelope into the bulk.
}

template<typename T, template<typename U> class Descriptor>
void AdvanceParticlesFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables;  // Particle field.
}


/* ******** AdvanceParticlesEveryWhereFunctional3D *********************************** */

template<typename T, template<typename U> class Descriptor>
AdvanceParticlesEveryWhereFunctional3D<T,Descriptor>::AdvanceParticlesEveryWhereFunctional3D (
        T cutOffValue_ )
  : cutOffValue(cutOffValue_)
{ }

template<typename T, template<typename U> class Descriptor>
void AdvanceParticlesEveryWhereFunctional3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    // Important: acts on the full domain, not only the provided one.
    particleField.advanceParticles(particleField.getBoundingBox(), cutOffValue);
}

template<typename T, template<typename U> class Descriptor>
AdvanceParticlesEveryWhereFunctional3D<T,Descriptor>* AdvanceParticlesEveryWhereFunctional3D<T,Descriptor>::clone() const {
    return new AdvanceParticlesEveryWhereFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void AdvanceParticlesEveryWhereFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables;  // Particle field.
}


/* ******** VerletUpdateVelocity3D *********************************** */

template<typename T, template<typename U> class Descriptor>
void VerletUpdateVelocity3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);

    std::vector<Particle3D<T,Descriptor>*> found;
    particleField.findParticles(domain, found);

    for (pluint iParticle=0; iParticle<found.size(); ++iParticle) {
        Particle3D<T,Descriptor>* nonTypeParticle = found[iParticle];
        VerletParticle3D<T,Descriptor>* particle =
            dynamic_cast<VerletParticle3D<T,Descriptor>*>(nonTypeParticle);
        PLB_ASSERT( particle );

        Array<T,3> a(particle->get_a());
        particle->set_v(particle->get_vHalfTime() + (T)0.5*a);
        PLB_ASSERT( norm(particle->get_v()) < 1. );
    }
}

template<typename T, template<typename U> class Descriptor>
VerletUpdateVelocity3D<T,Descriptor>* VerletUpdateVelocity3D<T,Descriptor>::clone() const {
    return new VerletUpdateVelocity3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void VerletUpdateVelocity3D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = true;  // Particle field.
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT VerletUpdateVelocity3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}

template<typename T, template<typename U> class Descriptor>
void VerletUpdateVelocity3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables; // Particle field.
}



/* ******** VerletUpdateVelocitySelective3D *********************************** */

template<typename T, template<typename U> class Descriptor>
VerletUpdateVelocitySelective3D<T,Descriptor>::VerletUpdateVelocitySelective3D(util::SelectInt* tags_)
    : tags(tags_)
{ }

template<typename T, template<typename U> class Descriptor>
VerletUpdateVelocitySelective3D<T,Descriptor>::~VerletUpdateVelocitySelective3D()
{
    delete tags;
}

template<typename T, template<typename U> class Descriptor>
VerletUpdateVelocitySelective3D<T,Descriptor>::VerletUpdateVelocitySelective3D(VerletUpdateVelocitySelective3D<T,Descriptor> const& rhs)
    : tags(rhs.tags->clone())
{ }

template<typename T, template<typename U> class Descriptor>
VerletUpdateVelocitySelective3D<T,Descriptor>&
    VerletUpdateVelocitySelective3D<T,Descriptor>::operator=(VerletUpdateVelocitySelective3D<T,Descriptor> const& rhs)
{
    VerletUpdateVelocitySelective3D<T,Descriptor>(rhs).swap(*this);
    return *this;
}

template<typename T, template<typename U> class Descriptor>
void VerletUpdateVelocitySelective3D<T,Descriptor>::swap(VerletUpdateVelocitySelective3D<T,Descriptor>& rhs) {
    std::swap(tags, rhs.tags);
}

template<typename T, template<typename U> class Descriptor>
void VerletUpdateVelocitySelective3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);

    std::vector<Particle3D<T,Descriptor>*> found;
    particleField.findParticles(domain, found);

    for (pluint iParticle=0; iParticle<found.size(); ++iParticle) {
        Particle3D<T,Descriptor>* nonTypeParticle = found[iParticle];
        if ((*tags)(nonTypeParticle->getTag())) {
            VerletParticle3D<T,Descriptor>* particle =
                dynamic_cast<VerletParticle3D<T,Descriptor>*>(nonTypeParticle);
            PLB_ASSERT( particle );

            Array<T,3> a(particle->get_a());
            particle->set_v(particle->get_vHalfTime() + (T)0.5*a);
            PLB_ASSERT( norm(particle->get_v()) < 1. );
        }
    }
}

template<typename T, template<typename U> class Descriptor>
VerletUpdateVelocitySelective3D<T,Descriptor>* VerletUpdateVelocitySelective3D<T,Descriptor>::clone() const {
    return new VerletUpdateVelocitySelective3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void VerletUpdateVelocitySelective3D<T,Descriptor>::getModificationPattern(std::vector<bool>& isWritten) const {
    isWritten[0] = true;  // Particle field.
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT VerletUpdateVelocitySelective3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}

template<typename T, template<typename U> class Descriptor>
void VerletUpdateVelocitySelective3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables; // Particle field.
}


template<typename T, template<typename U> class Descriptor>
void addWallParticles (
    MultiParticleField3D<DenseParticleField3D<T,Descriptor> >& particles,
    TriangleBoundary3D<T>& boundary )
{
    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(&particles);
    applyProcessingFunctional (
            new CreateParticleFromVertexNoTagging3D<T,Descriptor> (
                    boundary.getMesh(),
                    new RestParticle3D<T,Descriptor>(0, Array<T,3>((T)0.,(T)0.,(T)0.))),
            particles.getBoundingBox(), particleArg );
}

template<typename T, template<typename U> class Descriptor, class ParticleFieldT>
void addWallParticlesGeneric (
    MultiParticleField3D<ParticleFieldT>& particles,
    TriangleBoundary3D<T>& boundary )
{
    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(&particles);
    applyProcessingFunctional (
            new CreateParticleFromVertexNoTagging3D<T,Descriptor> (
                    boundary.getMesh(),
                    new RestParticle3D<T,Descriptor>(0, Array<T,3>((T)0.,(T)0.,(T)0.))),
            particles.getBoundingBox(), particleArg );
}


/* ******** CountAndAccumulateParticles3D *********************************** */

template<typename T, template<typename U> class Descriptor>
void CountAndAccumulateParticles3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==2 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    ScalarField3D<plint>& numParticlefield =
        *dynamic_cast<ScalarField3D<plint>*>(blocks[1]);
    Dot3D offset = computeRelativeDisplacement(particleField, numParticlefield);
    std::vector<Particle3D<T,Descriptor>*> particles;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                particleField.findParticles(Box3D(iX,iX,iY,iY,iZ,iZ), particles);
                numParticlefield.get(iX+offset.x,iY+offset.y,iZ+offset.z) += particles.size();
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
CountAndAccumulateParticles3D<T,Descriptor>* CountAndAccumulateParticles3D<T,Descriptor>::clone() const {
    return new CountAndAccumulateParticles3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void CountAndAccumulateParticles3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing;          // Particle field.
    modified[1] = modif::staticVariables;  // Scalar field.
}

/* ******** CountAndAccumulateTaggedParticles3D *********************************** */

template<typename T, template<typename U> class Descriptor>
CountAndAccumulateTaggedParticles3D<T,Descriptor>::CountAndAccumulateTaggedParticles3D(plint tag_)
    : tag(tag_)
{ }

template<typename T, template<typename U> class Descriptor>
void CountAndAccumulateTaggedParticles3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==2 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    ScalarField3D<plint>& numParticlefield =
        *dynamic_cast<ScalarField3D<plint>*>(blocks[1]);
    Dot3D offset = computeRelativeDisplacement(particleField, numParticlefield);
    std::vector<Particle3D<T,Descriptor>*> particles;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                particleField.findParticles(Box3D(iX,iX,iY,iY,iZ,iZ), particles);
                for (pluint iParticle=0; iParticle<particles.size(); ++iParticle) {
                    if (particles[iParticle]->getTag() == tag) {
                        ++numParticlefield.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                    }
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
CountAndAccumulateTaggedParticles3D<T,Descriptor>* CountAndAccumulateTaggedParticles3D<T,Descriptor>::clone() const {
    return new CountAndAccumulateTaggedParticles3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void CountAndAccumulateTaggedParticles3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing;          // Particle field.
    modified[1] = modif::staticVariables;  // Scalar field.
}

/* ******** CountAndAccumulateTaggedParticlesRefined3D *********************************** */

template<typename T, template<typename U> class Descriptor>
CountAndAccumulateTaggedParticlesRefined3D<T,Descriptor>::CountAndAccumulateTaggedParticlesRefined3D(
        plint tag_, plint dxScale_)
    : tag(tag_),
      dxScale(dxScale_)
{
    PLB_PRECONDITION(dxScale <= -1);
}

template<typename T, template<typename U> class Descriptor>
void CountAndAccumulateTaggedParticlesRefined3D<T,Descriptor>::processGenericBlocks(
        Box3D coarseDomain, std::vector<AtomicBlock3D*> blocks)
{
    PLB_PRECONDITION(blocks.size() == 2);

    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    ScalarField3D<plint>& numParticlefield =
        *dynamic_cast<ScalarField3D<plint>*>(blocks[1]);

    T stretch = util::twoToThePower(-dxScale);
    Dot3D refinedLocation = numParticlefield.getLocation();

    std::vector<Particle3D<T,Descriptor>*> particles;
    // Loop over the coarse particle grid and gather the particles.
    for (plint iX=coarseDomain.x0; iX<=coarseDomain.x1; ++iX) {
        for (plint iY=coarseDomain.y0; iY<=coarseDomain.y1; ++iY) {
            for (plint iZ=coarseDomain.z0; iZ<=coarseDomain.z1; ++iZ) {
                particleField.findParticles(Box3D(iX,iX,iY,iY,iZ,iZ), particles);
                for (pluint iParticle=0; iParticle<particles.size(); ++iParticle) {
                    // If the tag is acceptable, transform to the refined coordinates and count.
                    if (particles[iParticle]->getTag() == tag) {
                        Array<T,3> refinedParticlePosition = stretch *
                            particles[iParticle]->getPosition();
                        plint refinedX = -1;
                        plint refinedY = -1;
                        plint refinedZ = -1;
                        computeGridPosition(refinedParticlePosition, refinedLocation,
                                refinedX, refinedY, refinedZ);
                        if (contained(refinedX, refinedY, refinedZ, numParticlefield.getBoundingBox())) {
                            ++numParticlefield.get(refinedX, refinedY, refinedZ);
                        }
                    }
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
CountAndAccumulateTaggedParticlesRefined3D<T,Descriptor>*
CountAndAccumulateTaggedParticlesRefined3D<T,Descriptor>::clone() const
{
    return new CountAndAccumulateTaggedParticlesRefined3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void CountAndAccumulateTaggedParticlesRefined3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing;          // Particle field.
    modified[1] = modif::staticVariables;  // Scalar field.
}

template<typename T, template<typename U> class Descriptor>
plint CountAndAccumulateTaggedParticlesRefined3D<T,Descriptor>::nearestCell(T pos)
{
    T afterComma = pos-floor(pos);
    if (pos >= (T)0.) {
        if (afterComma > (T)0.5) {
            return (plint)(pos+(T)0.75);
        }
        else {
            return (plint)(pos);
        }
    }
    else {
        if (afterComma > (T)0.5) {
            return (plint)(pos);
        }
        else {
            return (plint)(pos-(T)0.75);
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void CountAndAccumulateTaggedParticlesRefined3D<T,Descriptor>::computeGridPosition(
        Array<T,3> const& position, Dot3D const& location, plint& iX, plint& iY, plint& iZ)
{
    iX = nearestCell(position[0]) - location.x;
    iY = nearestCell(position[1]) - location.y;
    iZ = nearestCell(position[2]) - location.z;
}

/* ******** CountTaggedParticles3D *********************************** */

template<typename T, template<typename U> class Descriptor>
CountTaggedParticles3D<T,Descriptor>::CountTaggedParticles3D(util::SelectInt* tags_)
    : tags(tags_)
{ }

template<typename T, template<typename U> class Descriptor>
CountTaggedParticles3D<T,Descriptor>::~CountTaggedParticles3D()
{
    delete tags;
}

template<typename T, template<typename U> class Descriptor>
CountTaggedParticles3D<T,Descriptor>::CountTaggedParticles3D(CountTaggedParticles3D<T,Descriptor> const& rhs)
    : tags(rhs.tags->clone())
{ }

template<typename T, template<typename U> class Descriptor>
CountTaggedParticles3D<T,Descriptor>&
    CountTaggedParticles3D<T,Descriptor>::operator=(CountTaggedParticles3D<T,Descriptor> const& rhs)
{
    CountTaggedParticles3D<T,Descriptor>(rhs).swap(*this);
    return *this;
}

template<typename T, template<typename U> class Descriptor>
void CountTaggedParticles3D<T,Descriptor>::swap(CountTaggedParticles3D<T,Descriptor>& rhs) {
    std::swap(tags, rhs.tags);
}

template<typename T, template<typename U> class Descriptor>
void CountTaggedParticles3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==2 );
    ParticleField3D<T,Descriptor>& particleField =
        *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
    ScalarField3D<plint>& numParticlefield =
        *dynamic_cast<ScalarField3D<plint>*>(blocks[1]);
    Dot3D offset = computeRelativeDisplacement(particleField, numParticlefield);
    std::vector<Particle3D<T,Descriptor>*> particles;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                particleField.findParticles(Box3D(iX,iX,iY,iY,iZ,iZ), particles);
                numParticlefield.get(iX+offset.x,iY+offset.y,iZ+offset.z) = 0;
                for (pluint iParticle=0; iParticle<particles.size(); ++iParticle) {
                    if ((*tags)(particles[iParticle]->getTag())) {
                        ++numParticlefield.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                    }
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
CountTaggedParticles3D<T,Descriptor>* CountTaggedParticles3D<T,Descriptor>::clone() const {
    return new CountTaggedParticles3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void CountTaggedParticles3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing;          // Particle field.
    modified[1] = modif::staticVariables;  // Scalar field.
}


template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
plint countParticles (
                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, Box3D const& domain )
{
    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(&particles);

    CountParticlesFunctional3D<T,Descriptor> functional;
    applyProcessingFunctional(functional, domain, particleArg);
    return functional.getNumParticles();
}

template< typename T, template<typename U> class Descriptor,
          template<typename T_, template<typename U_> class Descriptor_> class ParticleFieldT >
plint countParticles (
                MultiParticleField3D<ParticleFieldT<T,Descriptor> >& particles, Box3D const& domain, util::SelectInt* tags )
{
    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(&particles);

    CountParticlesSelectiveFunctional3D<T,Descriptor> functional(tags);
    applyProcessingFunctional(functional, domain, particleArg);
    return functional.getNumParticles();
}

template<typename T, template<typename U> class Descriptor>
void injectParticles(std::vector<Particle3D<T,Descriptor>*>& injectedParticles,
                     MultiParticleField3D<DenseParticleField3D<T,Descriptor> >& particles, Box3D domain)
{
    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(&particles);
    applyProcessingFunctional (
            new InjectParticlesFunctional3D<T,Descriptor>(injectedParticles), domain, particleArg );
}


}  // namespace plb

#endif  // PARTICLE_PROCESSING_FUNCTIONAL_3D_HH
