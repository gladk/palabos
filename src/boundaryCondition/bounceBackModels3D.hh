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

/** \file
 * BounceBack dynamics models in 3D -- generic implementation.
 */
#ifndef BOUNCE_BACK_MODELS_3D_HH
#define BOUNCE_BACK_MODELS_3D_HH

#include "boundaryCondition/bounceBackModels3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "core/cell.h"
#include "atomicBlock/dataProcessorWrapper3D.h"
#include "multiBlock/multiBlockGenerator3D.h"
#include <typeinfo>

namespace plb {

/* *************** Class InitializeMomentumExchangeFunctional3D ************* */

template<typename T, template<typename U> class Descriptor>
class InitializeMomentumExchangeFunctional3D : public BoxProcessingFunctional3D_L<T,Descriptor> {
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice)
    {
        static const int bounceBackId = BounceBack<T,Descriptor>().getId();
        static const int noDynamicsId = NoDynamics<T,Descriptor>().getId();
        static const Array<plint,3> tmp((plint) 0, (plint) 0, (plint) 0);
        static const int mEBounceBackId = MomentumExchangeBounceBack<T,Descriptor>(tmp).getId();

        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    Dynamics<T,Descriptor>& dynamics = lattice.get(iX,iY,iZ).getDynamics();
                    if (dynamics.getId() == mEBounceBackId) {
                        std::vector<plint> fluidDirections;
                        for (plint iPop=1; iPop<Descriptor<T>::q; ++iPop) {
                            plint nextX = iX + Descriptor<T>::c[iPop][0];
                            plint nextY = iY + Descriptor<T>::c[iPop][1];
                            plint nextZ = iZ + Descriptor<T>::c[iPop][2];
                            Dynamics<T,Descriptor>& partner = lattice.get(nextX,nextY,nextZ).getDynamics();
                            int partnerId = partner.getId();
                            if (partnerId != mEBounceBackId && partnerId != bounceBackId && partnerId != noDynamicsId) {
                                fluidDirections.push_back(iPop);
                            }
                        }

                        if (!fluidDirections.empty()) {
                            MomentumExchangeBounceBack<T,Descriptor>& bounceBackDynamics =
                                dynamic_cast<MomentumExchangeBounceBack<T,Descriptor>&>(dynamics);
                            bounceBackDynamics.setFluidDirections(fluidDirections);
                        }
                    }
                }
            }
        }
    }
    virtual BlockDomain::DomainT appliesTo() const {
        // Apply only to bulk, because envelopes should not contribute
        //   to the statistics.
        return BlockDomain::bulk;
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
    }
    virtual InitializeMomentumExchangeFunctional3D<T,Descriptor>* clone() const 
    {
        return new InitializeMomentumExchangeFunctional3D<T,Descriptor>(*this);
    }
};

/* ************* Class MomentumExchangeComplexDomainFunctional3D ** */

template<typename T, template<typename U> class Descriptor>
class MomentumExchangeComplexDomainFunctional3D
    : public BoxProcessingFunctional3D_L<T,Descriptor>
{
public:
    MomentumExchangeComplexDomainFunctional3D( DomainFunctional3D* domain_ )
        : domain(domain_)
    { }
    MomentumExchangeComplexDomainFunctional3D (
            MomentumExchangeComplexDomainFunctional3D<T,Descriptor> const& rhs )
        : domain(rhs.domain->clone())
    { }
    MomentumExchangeComplexDomainFunctional3D<T,Descriptor>& operator= (
            MomentumExchangeComplexDomainFunctional3D<T,Descriptor> const& rhs )
    {
        delete domain; domain = rhs.domain->clone();
        return *this;
    }

    virtual ~MomentumExchangeComplexDomainFunctional3D() {
        delete domain;
    }
    virtual void process(Box3D boundingBox, BlockLattice3D<T,Descriptor>& lattice)
    {
        static const int bounceBackId = BounceBack<T,Descriptor>().getId();
        static const int noDynamicsId = NoDynamics<T,Descriptor>().getId();
        static const Array<plint,3> tmp((plint) 0, (plint) 0, (plint) 0);
        static const int mEBounceBackId = MomentumExchangeBounceBack<T,Descriptor>(tmp).getId();

        Dot3D relativeOffset = lattice.getLocation();
        for (plint iX=boundingBox.x0; iX<=boundingBox.x1; ++iX) {
            for (plint iY=boundingBox.y0; iY<=boundingBox.y1; ++iY) {
                for (plint iZ=boundingBox.z0; iZ<=boundingBox.z1; ++iZ) {
                    if ((*domain)(iX+relativeOffset.x,iY+relativeOffset.y,iZ+relativeOffset.z)) {
                        Dynamics<T,Descriptor>& dynamics = lattice.get(iX,iY,iZ).getDynamics();
                        if (dynamics.getId() == mEBounceBackId) {
                            std::vector<plint> fluidDirections;
                            for (plint iPop=1; iPop<Descriptor<T>::q; ++iPop) {
                                plint nextX = iX + Descriptor<T>::c[iPop][0];
                                plint nextY = iY + Descriptor<T>::c[iPop][1];
                                plint nextZ = iZ + Descriptor<T>::c[iPop][2];
                                Dynamics<T,Descriptor>& partner = lattice.get(nextX,nextY,nextZ).getDynamics();
                                int partnerId = partner.getId();
                                if (partnerId != mEBounceBackId && partnerId != bounceBackId && partnerId != noDynamicsId) {
                                    fluidDirections.push_back(iPop);
                                }
                            }

                            if (!fluidDirections.empty()) {
                                MomentumExchangeBounceBack<T,Descriptor>& bounceBackDynamics =
                                    dynamic_cast<MomentumExchangeBounceBack<T,Descriptor>&>(dynamics);
                                bounceBackDynamics.setFluidDirections(fluidDirections);
                            }
                        }
                    }
                }
            }
        }
    }
    virtual BlockDomain::DomainT appliesTo() const {
        // Apply only to bulk, because envelopes should not contribute
        //   to the statistics.
        return BlockDomain::bulk;
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
    }
    virtual MomentumExchangeComplexDomainFunctional3D<T,Descriptor>* clone() const 
    {
        return new MomentumExchangeComplexDomainFunctional3D<T,Descriptor>(*this);
    }
private:
    DomainFunctional3D* domain;
};

template<typename T, template<typename U> class Descriptor>
class InitializeDotMomentumExchangeFunctional3D : public DotProcessingFunctional3D_L<T,Descriptor> {
public:
    virtual void process(DotList3D const& dotList, BlockLattice3D<T,Descriptor>& lattice)
    {
        static const int bounceBackId = BounceBack<T,Descriptor>().getId();
        static const int noDynamicsId = NoDynamics<T,Descriptor>().getId();
        static const Array<plint,3> tmp((plint) 0, (plint) 0, (plint) 0);
        static const int mEBounceBackId = MomentumExchangeBounceBack<T,Descriptor>(tmp).getId();

        for (plint iDot=0; iDot<dotList.getN(); ++iDot) {
            Dot3D const& dot = dotList.getDot(iDot);
            plint iX = dot.x;
            plint iY = dot.y;
            plint iZ = dot.z;
            Dynamics<T,Descriptor>& dynamics = lattice.get(iX,iY,iZ).getDynamics();
            if (dynamics.getId() == mEBounceBackId) {
                std::vector<plint> fluidDirections;
                for (plint iPop=1; iPop<Descriptor<T>::q; ++iPop) {
                    plint nextX = iX + Descriptor<T>::c[iPop][0];
                    plint nextY = iY + Descriptor<T>::c[iPop][1];
                    plint nextZ = iZ + Descriptor<T>::c[iPop][2];
                    Dynamics<T,Descriptor>& partner = lattice.get(nextX,nextY,nextZ).getDynamics();
                    int partnerId = partner.getId();
                    if (partnerId != mEBounceBackId && partnerId != bounceBackId && partnerId != noDynamicsId) {
                        fluidDirections.push_back(iPop);
                    }
                }

                if (!fluidDirections.empty()) {
                    MomentumExchangeBounceBack<T,Descriptor>& bounceBackDynamics =
                        dynamic_cast<MomentumExchangeBounceBack<T,Descriptor>&>(dynamics);
                    bounceBackDynamics.setFluidDirections(fluidDirections);
                }
            }
        }
    }
    virtual BlockDomain::DomainT appliesTo() const {
        // Apply only to bulk, because envelopes should not contribute
        //   to the statistics.
        return BlockDomain::bulk;
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
    }
    virtual InitializeDotMomentumExchangeFunctional3D<T,Descriptor>* clone() const 
    {
        return new InitializeDotMomentumExchangeFunctional3D<T,Descriptor>(*this);
    }
};


template<typename T, template<class U> class Descriptor>
void initializeMomentumExchange (
        BlockLattice3D<T,Descriptor>& lattice, Box3D domain )
{
    applyProcessingFunctional (
        new InitializeMomentumExchangeFunctional3D<T,Descriptor>(), domain, lattice );
}

template<typename T, template<class U> class Descriptor>
void initializeMomentumExchange (
        BlockLattice3D<T,Descriptor>& lattice, Box3D boundingBox,
        DomainFunctional3D* domain )
{
    applyProcessingFunctional (
            new MomentumExchangeComplexDomainFunctional3D<T,Descriptor>(domain), boundingBox, lattice );
}

template<typename T, template<class U> class Descriptor>
void initializeMomentumExchange (
        BlockLattice3D<T,Descriptor>& lattice, DotList3D const& dotList )
{
    applyProcessingFunctional (
        new InitializeDotMomentumExchangeFunctional3D<T,Descriptor>(), dotList, lattice );
}


template<typename T, template<class U> class Descriptor>
void initializeMomentumExchange (
        MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain )
{
    applyProcessingFunctional (
        new InitializeMomentumExchangeFunctional3D<T,Descriptor>(), domain, lattice );
}

template<typename T, template<class U> class Descriptor>
void initializeMomentumExchange (
        MultiBlockLattice3D<T,Descriptor>& lattice, Box3D boundingBox,
        DomainFunctional3D* domain )
{
    applyProcessingFunctional (
            new MomentumExchangeComplexDomainFunctional3D<T,Descriptor>(domain), boundingBox, lattice );
}

template<typename T, template<class U> class Descriptor>
void initializeMomentumExchange (
        MultiBlockLattice3D<T,Descriptor>& lattice, DotList3D const& dotList )
{
    applyProcessingFunctional (
        new InitializeDotMomentumExchangeFunctional3D<T,Descriptor>(), dotList, lattice );
}


template<typename T, template<typename U> class Descriptor> 
void CountBBNeighborsFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<int>& neighbors)
{
    static plint bbid = BounceBack<T,Descriptor>().getId();
    Dot3D offset = computeRelativeDisplacement(lattice, neighbors);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint numNeighbors=0;
                if (lattice.get(iX,iY,iZ).getDynamics().getId()!=bbid) {
                    for (plint dx=-1; dx<=+1; ++dx) {
                        for (plint dy=-1; dy<=+1; ++dy) {
                            for (plint dz=-1; dz<=+1; ++dz)  {
                                if (lattice.get(iX+dx,iY+dy,iZ+dz).getDynamics().getId()==bbid) {
                                    ++numNeighbors;
                                }
                            }
                        }
                    }
                }
                neighbors.get(iX+offset.x,iY+offset.y,iZ+offset.z) = numNeighbors;
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
CountBBNeighborsFunctional3D<T,Descriptor>* CountBBNeighborsFunctional3D<T,Descriptor>::clone() const
{
    return new CountBBNeighborsFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void CountBBNeighborsFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT CountBBNeighborsFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}



template<typename T, template<typename U> class Descriptor> 
void CountBBNeighbors_NTensor3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice, NTensorField3D<int>& neighbors )
{
    PLB_PRECONDITION( neighbors.getNdim() == 1 );
    static plint bbid = BounceBack<T,Descriptor>().getId();
    Dot3D offset = computeRelativeDisplacement(lattice, neighbors);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint numNeighbors=0;
                if (lattice.get(iX,iY,iZ).getDynamics().getId()!=bbid) {
                    for (plint dx=-1; dx<=+1; ++dx) {
                        for (plint dy=-1; dy<=+1; ++dy) {
                            for (plint dz=-1; dz<=+1; ++dz)  {
                                if (lattice.get(iX+dx,iY+dy,iZ+dz).getDynamics().getId()==bbid) {
                                    ++numNeighbors;
                                }
                            }
                        }
                    }
                }
                *neighbors.get(iX+offset.x,iY+offset.y,iZ+offset.z) = numNeighbors;
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
CountBBNeighbors_NTensor3D<T,Descriptor>* CountBBNeighbors_NTensor3D<T,Descriptor>::clone() const
{
    return new CountBBNeighbors_NTensor3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void CountBBNeighbors_NTensor3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT CountBBNeighbors_NTensor3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T, template<typename U> class Descriptor> 
void MaskedCountBBNeighbors_NTensor3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
        NTensorField3D<int>& neighbors, NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( neighbors.getNdim() == 1 );
    static plint bbid = BounceBack<T,Descriptor>().getId();
    Dot3D offset = computeRelativeDisplacement(lattice, neighbors);
    Dot3D maskOfs = computeRelativeDisplacement(lattice, mask);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                int numNeighbors=0;
                if ( lattice.get(iX,iY,iZ).getDynamics().getId()!=bbid &&
                     *mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z) )
                {
                    for (plint dx=-1; dx<=+1; ++dx) {
                        for (plint dy=-1; dy<=+1; ++dy) {
                            for (plint dz=-1; dz<=+1; ++dz)  {
                                if (lattice.get(iX+dx,iY+dy,iZ+dz).getDynamics().getId()==bbid) {
                                    ++numNeighbors;
                                }
                            }
                        }
                    }
                }
                *neighbors.get(iX+offset.x,iY+offset.y,iZ+offset.z) = numNeighbors;
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
MaskedCountBBNeighbors_NTensor3D<T,Descriptor>* MaskedCountBBNeighbors_NTensor3D<T,Descriptor>::clone() const
{
    return new MaskedCountBBNeighbors_NTensor3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void MaskedCountBBNeighbors_NTensor3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT MaskedCountBBNeighbors_NTensor3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T, template<typename U> class Descriptor> 
void countBBNeighbors( MultiBlockLattice3D<T,Descriptor>& lattice,
                       MultiNTensorField3D<int>& neighbors, Box3D domain )
{
    applyProcessingFunctional (
            new CountBBNeighbors_NTensor3D<T,Descriptor>, domain, lattice, neighbors );
}

template<typename T, template<typename U> class Descriptor> 
void maskedCountBBNeighbors( MultiBlockLattice3D<T,Descriptor>& lattice,
                       MultiNTensorField3D<int>& neighbors,
                       MultiNTensorField3D<int>& mask, Box3D domain )
{
    applyProcessingFunctional (
            new MaskedCountBBNeighbors_NTensor3D<T,Descriptor>, domain, lattice, neighbors, mask );
}

template<typename T, template<typename U> class Descriptor> 
MultiNTensorField3D<int>* countBBNeighbors(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain) {
    MultiNTensorField3D<int>* neighbors
        = generateMultiNTensorField<int>(lattice, domain, 1);
    countBBNeighbors(lattice, *neighbors, domain);
    return neighbors;
}

template<typename T, template<typename U> class Descriptor> 
MultiNTensorField3D<int>* maskedCountBBNeighbors( MultiBlockLattice3D<T,Descriptor>& lattice,
                                                  MultiNTensorField3D<int>& mask, Box3D domain )
{
    MultiNTensorField3D<int>* neighbors
        = generateMultiNTensorField<int>(lattice, domain, 1);
    maskedCountBBNeighbors(lattice, *neighbors, mask, domain);
    return neighbors;
}

}  // namespace plb

#endif  // BOUNCE_BACK_MODELS_3D_H
