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
 * BounceBack dynamics models in 2D -- generic implementation.
 */

#ifndef BOUNCE_BACK_MODELS_2D_HH
#define BOUNCE_BACK_MODELS_2D_HH

#include "boundaryCondition/bounceBackModels2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "core/cell.h"
#include "atomicBlock/dataProcessorWrapper2D.h"

namespace plb {

/* *************** Class InitializeMomentumExchangeFunctional2D ************* */

template<typename T, template<typename U> class Descriptor>
class InitializeMomentumExchangeFunctional2D : public BoxProcessingFunctional2D_L<T,Descriptor> {
public:
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice)
    {
        static const int mEBounceBackId = MomentumExchangeBounceBack<T,Descriptor>(Array<plint,2>::zero()).getId();

        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                Dynamics<T,Descriptor>& dynamics = lattice.get(iX,iY).getDynamics();
                if (dynamics.getId() == mEBounceBackId) {
                    std::vector<plint> fluidDirections;
                    for (plint iPop=1; iPop<Descriptor<T>::q; ++iPop) {
                        plint nextX = iX + Descriptor<T>::c[iPop][0];
                        plint nextY = iY + Descriptor<T>::c[iPop][1];
                        Dynamics<T,Descriptor>& partner = lattice.get(nextX,nextY).getDynamics();
                        if (partner.hasMoments()) {
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
    virtual BlockDomain::DomainT appliesTo() const {
        // Apply only to bulk, because envelopes should not contribute
        //   to the statistics.
        return BlockDomain::bulk;
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::dynamicVariables;
    }
    virtual InitializeMomentumExchangeFunctional2D<T,Descriptor>* clone() const
    {
        return new InitializeMomentumExchangeFunctional2D<T,Descriptor>(*this);
    }
};

/* ************* Class MomentumExchangeComplexDomainFunctional2D ** */

template<typename T, template<typename U> class Descriptor>
class MomentumExchangeComplexDomainFunctional2D
    : public BoxProcessingFunctional2D_L<T,Descriptor>
{
public:
    MomentumExchangeComplexDomainFunctional2D( DomainFunctional2D* domain_ )
        : domain(domain_)
    { }
    MomentumExchangeComplexDomainFunctional2D (
            MomentumExchangeComplexDomainFunctional2D<T,Descriptor> const& rhs )
        : domain(rhs.domain->clone())
    { }
    MomentumExchangeComplexDomainFunctional2D<T,Descriptor>& operator= (
            MomentumExchangeComplexDomainFunctional2D<T,Descriptor> const& rhs )
    {
        delete domain; domain = rhs.domain->clone();
        return *this;
    }

    virtual ~MomentumExchangeComplexDomainFunctional2D() {
        delete domain;
    }
    virtual void process(Box2D boundingBox, BlockLattice2D<T,Descriptor>& lattice)
    {
        static const int mEBounceBackId = MomentumExchangeBounceBack<T,Descriptor>(Array<plint,2>::zero()).getId();

        Dot2D relativeOffset = lattice.getLocation();
        for (plint iX=boundingBox.x0; iX<=boundingBox.x1; ++iX) {
            for (plint iY=boundingBox.y0; iY<=boundingBox.y1; ++iY) {
                if ((*domain)(iX+relativeOffset.x,iY+relativeOffset.y)) {
                    Dynamics<T,Descriptor>& dynamics = lattice.get(iX,iY).getDynamics();
                    if (dynamics.getId() == mEBounceBackId) {
                        std::vector<plint> fluidDirections;
                        for (plint iPop=1; iPop<Descriptor<T>::q; ++iPop) {
                            plint nextX = iX + Descriptor<T>::c[iPop][0];
                            plint nextY = iY + Descriptor<T>::c[iPop][1];
                            Dynamics<T,Descriptor>& partner = lattice.get(nextX,nextY).getDynamics();
                            if (partner.hasMoments()) {
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
        modified[0] = modif::dynamicVariables;
    }
    virtual MomentumExchangeComplexDomainFunctional2D<T,Descriptor>* clone() const
    {
        return new MomentumExchangeComplexDomainFunctional2D<T,Descriptor>(*this);
    }
private:
    DomainFunctional2D* domain;
};

/* ************* Class InitializeDotMomentumExchangeFunctional2D ** */

template<typename T, template<typename U> class Descriptor>
class InitializeDotMomentumExchangeFunctional2D : public DotProcessingFunctional2D_L<T,Descriptor> {
public:
    virtual void process(DotList2D const& dotList, BlockLattice2D<T,Descriptor>& lattice)
    {
        static const int mEBounceBackId = MomentumExchangeBounceBack<T,Descriptor>(Array<plint,2>::zero()).getId();

        for (plint iDot=0; iDot<dotList.getN(); ++iDot) {
            Dot2D const& dot = dotList.getDot(iDot);
            plint iX = dot.x;
            plint iY = dot.y;
            Dynamics<T,Descriptor>& dynamics = lattice.get(iX,iY).getDynamics();
            if (dynamics.getId() == mEBounceBackId) {
                std::vector<plint> fluidDirections;
                for (plint iPop=1; iPop<Descriptor<T>::q; ++iPop) {
                    plint nextX = iX + Descriptor<T>::c[iPop][0];
                    plint nextY = iY + Descriptor<T>::c[iPop][1];
                    Dynamics<T,Descriptor>& partner = lattice.get(nextX,nextY).getDynamics();
                    if (partner.hasMoments()) {
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
        modified[0] = modif::dynamicVariables;
    }
    virtual InitializeDotMomentumExchangeFunctional2D<T,Descriptor>* clone() const
    {
        return new InitializeDotMomentumExchangeFunctional2D<T,Descriptor>(*this);
    }
};

template<typename T, template<class U> class Descriptor>
void initializeMomentumExchange (
        BlockLattice2D<T,Descriptor>& lattice, Box2D domain )
{
    applyProcessingFunctional (
        new InitializeMomentumExchangeFunctional2D<T,Descriptor>(), domain, lattice );
}

template<typename T, template<class U> class Descriptor>
void initializeMomentumExchange (
        BlockLattice2D<T,Descriptor>& lattice, Box2D boundingBox,
        DomainFunctional2D* domain )
{
    applyProcessingFunctional (
            new MomentumExchangeComplexDomainFunctional2D<T,Descriptor>(domain), boundingBox, lattice );
}

template<typename T, template<class U> class Descriptor>
void initializeMomentumExchange (
        BlockLattice2D<T,Descriptor>& lattice, DotList2D const& dotList )
{
    applyProcessingFunctional (
        new InitializeDotMomentumExchangeFunctional2D<T,Descriptor>(), dotList, lattice );
}

template<typename T, template<class U> class Descriptor>
void initializeMomentumExchange (
        MultiBlockLattice2D<T,Descriptor>& lattice, Box2D domain )
{
    applyProcessingFunctional (
        new InitializeMomentumExchangeFunctional2D<T,Descriptor>(), domain, lattice );
}

template<typename T, template<class U> class Descriptor>
void initializeMomentumExchange (
        MultiBlockLattice2D<T,Descriptor>& lattice, Box2D boundingBox,
        DomainFunctional2D* domain )
{
    applyProcessingFunctional (
            new MomentumExchangeComplexDomainFunctional2D<T,Descriptor>(domain), boundingBox, lattice );
}

template<typename T, template<class U> class Descriptor>
void initializeMomentumExchange (
        MultiBlockLattice2D<T,Descriptor>& lattice, DotList2D const& dotList )
{
    applyProcessingFunctional (
        new InitializeDotMomentumExchangeFunctional2D<T,Descriptor>(), dotList, lattice );
}

}  // namespace plb

#endif  // BOUNCE_BACK_MODELS_2D_HH
