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
 * BounceBack dynamics models in 3D -- header file.
 */

#ifndef BOUNCE_BACK_MODELS_3D_H
#define BOUNCE_BACK_MODELS_3D_H

#include "boundaryCondition/bounceBackModels.h"
#include "dataProcessors/dataInitializerFunctional3D.h"

namespace plb {

template<typename T, template<class U> class Descriptor>
void initializeMomentumExchange (
        BlockLattice3D<T,Descriptor>& lattice, Box3D domain );

template<typename T, template<class U> class Descriptor>
void initializeMomentumExchange (
        BlockLattice3D<T,Descriptor>& lattice, Box3D boundingBox,
        DomainFunctional3D* domain );

template<typename T, template<class U> class Descriptor>
void initializeMomentumExchange (
        BlockLattice3D<T,Descriptor>& lattice, DotList3D const& dotList );


template<typename T, template<class U> class Descriptor>
void initializeMomentumExchange (
        MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain );

template<typename T, template<class U> class Descriptor>
void initializeMomentumExchange (
        MultiBlockLattice3D<T,Descriptor>& lattice, Box3D boundingBox,
        DomainFunctional3D* domain );

template<typename T, template<class U> class Descriptor>
void initializeMomentumExchange (
        MultiBlockLattice3D<T,Descriptor>& lattice, DotList3D const& dotList );



template<typename T, template<typename U> class Descriptor> 
class CountBBNeighborsFunctional3D : public BoxProcessingFunctional3D_LS<T,Descriptor,int>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       ScalarField3D<int>& neighbors);
    virtual CountBBNeighborsFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};


template<typename T, template<typename U> class Descriptor> 
class CountBBNeighbors_NTensor3D : public BoxProcessingFunctional3D_LN<T,Descriptor,int>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       NTensorField3D<int>& neighbors);
    virtual CountBBNeighbors_NTensor3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};


template<typename T, template<typename U> class Descriptor> 
class MaskedCountBBNeighbors_NTensor3D : public MaskedBoxProcessingFunctional3D_LN<T,Descriptor,int>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       NTensorField3D<int>& neighbors,
                                       NTensorField3D<int>& mask);
    virtual MaskedCountBBNeighbors_NTensor3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
void countBBNeighbors( MultiBlockLattice3D<T,Descriptor>& lattice,
                       MultiNTensorField3D<int>& neighbors, Box3D domain );

template<typename T, template<typename U> class Descriptor> 
void maskedCountBBNeighbors( MultiBlockLattice3D<T,Descriptor>& lattice,
                       MultiNTensorField3D<int>& neighbors,
                       MultiNTensorField3D<int>& mask, Box3D domain );

template<typename T, template<typename U> class Descriptor> 
MultiNTensorField3D<int>* countBBNeighbors(MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain);

template<typename T, template<typename U> class Descriptor> 
MultiNTensorField3D<int>* maskedCountBBNeighbors( MultiBlockLattice3D<T,Descriptor>& lattice,
                                                MultiNTensorField3D<int>& mask, Box3D domain);


/* *************** Class ComputeMomentumExchangeFunctional3D ************* */

template<typename T, template<typename U> class Descriptor>
class ComputeMomentumExchangeFunctional3D : public ReductiveBoxProcessingFunctional3D_L<T,Descriptor> {
public:
    ComputeMomentumExchangeFunctional3D();
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice);
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual ComputeMomentumExchangeFunctional3D<T,Descriptor>* clone() const;
    Array<T,3> getForce() const;
private:
    Array<plint,3> forceIds;
};

template<typename T, template<typename U> class Descriptor>
class MaskedComputeMomentumExchangeFunctional3D : public ReductiveBoxProcessingFunctional3D_LS<T,Descriptor,int> {
public:
    MaskedComputeMomentumExchangeFunctional3D();
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<int>& mask);
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual MaskedComputeMomentumExchangeFunctional3D<T,Descriptor>* clone() const;
    Array<T,3> getForce() const;
private:
    Array<plint,3> forceIds;
};



}  // namespace plb

#endif  // BOUNCE_BACK_MODELS_3D_H

