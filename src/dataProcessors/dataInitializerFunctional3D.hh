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
 * Functionals for domain initialization -- generic implementation.
 */
#ifndef DATA_INITIALIZER_FUNCTIONAL_3D_HH
#define DATA_INITIALIZER_FUNCTIONAL_3D_HH

#include "dataProcessors/dataInitializerFunctional3D.h"
#include "core/cell.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "atomicBlock/blockLattice3D.h"
#include "multiGrid/multiGridUtil.h"

namespace plb {

/* *************** PART I ******************************************** */
/* *************** Initialization of the block-lattice *************** */
/* ******************************************************************* */

template<typename T, template<class U> class Descriptor>
OneCellFunctional3D<T,Descriptor>::~OneCellFunctional3D()
{ }

template<typename T, template<class U> class Descriptor>
BlockDomain::DomainT OneCellFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<class U> class Descriptor>
void OneCellFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
}


template<typename T, template<class U> class Descriptor>
void OneCellFunctional3D<T,Descriptor>::setscale(int dxScale, int dtScale)
{ }


template<typename T, template<class U> class Descriptor>
OneCellIndexedFunctional3D<T,Descriptor>::~OneCellIndexedFunctional3D()
{ }

template<typename T, template<class U> class Descriptor>
BlockDomain::DomainT OneCellIndexedFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<class U> class Descriptor>
void OneCellIndexedFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::staticVariables;
}

template<typename T, template<class U> class Descriptor>
void OneCellIndexedFunctional3D<T,Descriptor>::setscale(int dxScale, int dtScale)
{ }


template<typename T, template<class U> class Descriptor>
GenericLatticeFunctional3D<T,Descriptor>::GenericLatticeFunctional3D (
        OneCellFunctional3D<T,Descriptor>* f_ )
    : f(f_)
{ }

template<typename T, template<class U> class Descriptor>
GenericLatticeFunctional3D<T,Descriptor>::GenericLatticeFunctional3D (
        GenericLatticeFunctional3D<T,Descriptor> const& rhs )
    : f(rhs.f->clone())
{ }

template<typename T, template<class U> class Descriptor>
GenericLatticeFunctional3D<T,Descriptor>::~GenericLatticeFunctional3D() {
    delete f;
}

template<typename T, template<class U> class Descriptor>
GenericLatticeFunctional3D<T,Descriptor>& GenericLatticeFunctional3D<T,Descriptor>::operator= (
        GenericLatticeFunctional3D<T,Descriptor> const& rhs )
{
    delete f;
    f = rhs.f->clone();
    return *this;
}

template<typename T, template<class U> class Descriptor>
void GenericLatticeFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                f->execute(lattice.get(iX,iY,iZ));
            }
        }
    }
}

template<typename T, template<class U> class Descriptor>
GenericLatticeFunctional3D<T,Descriptor>*
    GenericLatticeFunctional3D<T,Descriptor>::clone() const
{
    return new GenericLatticeFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<class U> class Descriptor>
void GenericLatticeFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    f->getTypeOfModification(modified);
}

template<typename T, template<class U> class Descriptor>
BlockDomain::DomainT GenericLatticeFunctional3D<T,Descriptor>::appliesTo() const {
    return f->appliesTo();
}

template<typename T, template<class U> class Descriptor>
void GenericLatticeFunctional3D<T,Descriptor>::setscale(int dxScale, int dtScale) {
    f->setscale(dxScale, dtScale);
}



template<typename T, template<class U> class Descriptor>
GenericIndexedLatticeFunctional3D<T,Descriptor>::GenericIndexedLatticeFunctional3D (
        OneCellIndexedFunctional3D<T,Descriptor>* f_ )
    : f(f_)
{ }

template<typename T, template<class U> class Descriptor>
GenericIndexedLatticeFunctional3D<T,Descriptor>::GenericIndexedLatticeFunctional3D (
        GenericIndexedLatticeFunctional3D<T,Descriptor> const& rhs )
    : f(rhs.f->clone())
{ }

template<typename T, template<class U> class Descriptor>
GenericIndexedLatticeFunctional3D<T,Descriptor>::~GenericIndexedLatticeFunctional3D() {
    delete f;
}

template<typename T, template<class U> class Descriptor>
GenericIndexedLatticeFunctional3D<T,Descriptor>&
    GenericIndexedLatticeFunctional3D<T,Descriptor>::operator= (
        GenericIndexedLatticeFunctional3D<T,Descriptor> const& rhs )
{
    delete f;
    f = rhs.f->clone();
    return *this;
}

template<typename T, template<class U> class Descriptor>
void GenericIndexedLatticeFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    Dot3D relativeOffset = lattice.getLocation();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                f->execute ( iX+relativeOffset.x,
                             iY+relativeOffset.y,
                             iZ+relativeOffset.z,
                             lattice.get(iX,iY,iZ) );
            }
        }
    }
}

template<typename T, template<class U> class Descriptor>
GenericIndexedLatticeFunctional3D<T,Descriptor>*
    GenericIndexedLatticeFunctional3D<T,Descriptor>::clone() const
{
    return new GenericIndexedLatticeFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<class U> class Descriptor>
void GenericIndexedLatticeFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    f->getTypeOfModification(modified);
}

template<typename T, template<class U> class Descriptor>
BlockDomain::DomainT GenericIndexedLatticeFunctional3D<T,Descriptor>::appliesTo() const
{
    return f->appliesTo();
}

template<typename T, template<class U> class Descriptor>
void GenericIndexedLatticeFunctional3D<T,Descriptor>::setscale(int dxScale, int dtScale)
{
    f->setscale(dxScale, dtScale);
}


/* *************** Class InstantiateDynamicsFunctional3D ************* */

template<typename T, template<typename U> class Descriptor>
InstantiateDynamicsFunctional3D<T,Descriptor>::InstantiateDynamicsFunctional3D (
        Dynamics<T,Descriptor>* dynamics_ )
    : dynamics(dynamics_)
{ }

template<typename T, template<typename U> class Descriptor>
InstantiateDynamicsFunctional3D<T,Descriptor>::InstantiateDynamicsFunctional3D (
        InstantiateDynamicsFunctional3D<T,Descriptor> const& rhs )
    : dynamics(rhs.dynamics->clone())
{ }

template<typename T, template<typename U> class Descriptor>
InstantiateDynamicsFunctional3D<T,Descriptor>&
    InstantiateDynamicsFunctional3D<T,Descriptor>::operator= (
        InstantiateDynamicsFunctional3D<T,Descriptor> const& rhs )
{
    delete dynamics; dynamics = rhs.dynamics->clone();
    return *this;
}


template<typename T, template<typename U> class Descriptor>
InstantiateDynamicsFunctional3D<T,Descriptor>::~InstantiateDynamicsFunctional3D() {
    delete dynamics;
}

template<typename T, template<typename U> class Descriptor>
void InstantiateDynamicsFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                lattice.attributeDynamics(iX,iY,iZ, dynamics->clone());
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT InstantiateDynamicsFunctional3D<T,Descriptor>::appliesTo() const {
    // Dynamics needs to be instantiated everywhere, including envelope.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void InstantiateDynamicsFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor>
InstantiateDynamicsFunctional3D<T,Descriptor>*
    InstantiateDynamicsFunctional3D<T,Descriptor>::clone() const 
{
    return new InstantiateDynamicsFunctional3D<T,Descriptor>(*this);
}


/* ************* Class InstantiateComplexDomainDynamicsFunctional3D ** */

template<typename T, template<typename U> class Descriptor>
InstantiateComplexDomainDynamicsFunctional3D<T,Descriptor>::InstantiateComplexDomainDynamicsFunctional3D (
        Dynamics<T,Descriptor>* dynamics_, DomainFunctional3D* domain_ )
    : dynamics(dynamics_),
      domain(domain_)
{ }

template<typename T, template<typename U> class Descriptor>
InstantiateComplexDomainDynamicsFunctional3D<T,Descriptor>::InstantiateComplexDomainDynamicsFunctional3D (
        InstantiateComplexDomainDynamicsFunctional3D<T,Descriptor> const& rhs )
    : dynamics(rhs.dynamics->clone()),
      domain(rhs.domain->clone())
{ }

template<typename T, template<typename U> class Descriptor>
InstantiateComplexDomainDynamicsFunctional3D<T,Descriptor>&
    InstantiateComplexDomainDynamicsFunctional3D<T,Descriptor>::operator= (
        InstantiateComplexDomainDynamicsFunctional3D<T,Descriptor> const& rhs )
{
    delete dynamics; dynamics = rhs.dynamics->clone();
    delete domain; domain = rhs.domain->clone();
    return *this;
}

template<typename T, template<typename U> class Descriptor>
InstantiateComplexDomainDynamicsFunctional3D<T,Descriptor>::~InstantiateComplexDomainDynamicsFunctional3D()
{
    delete dynamics;
    delete domain;
}

template<typename T, template<typename U> class Descriptor>
void InstantiateComplexDomainDynamicsFunctional3D<T,Descriptor>::process (
        Box3D boundingBox, BlockLattice3D<T,Descriptor>& lattice )
{
    Dot3D relativeOffset = lattice.getLocation();
    for (plint iX=boundingBox.x0; iX<=boundingBox.x1; ++iX) {
        for (plint iY=boundingBox.y0; iY<=boundingBox.y1; ++iY) {
            for (plint iZ=boundingBox.z0; iZ<=boundingBox.z1; ++iZ) {
                if ((*domain)(iX+relativeOffset.x,iY+relativeOffset.y,iZ+relativeOffset.z)) {
                    lattice.attributeDynamics(iX,iY,iZ, dynamics->clone());
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT InstantiateComplexDomainDynamicsFunctional3D<T,Descriptor>::appliesTo() const
{
    // Dynamics needs to be instantiated everywhere, including envelope.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void InstantiateComplexDomainDynamicsFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor>
InstantiateComplexDomainDynamicsFunctional3D<T,Descriptor>*
    InstantiateComplexDomainDynamicsFunctional3D<T,Descriptor>::clone() const 
{
    return new InstantiateComplexDomainDynamicsFunctional3D<T,Descriptor>(*this);
}


/* ************* Class InstantiateDotDynamicsFunctional3D ******************* */

template<typename T, template<typename U> class Descriptor>
InstantiateDotDynamicsFunctional3D<T,Descriptor>::InstantiateDotDynamicsFunctional3D (
        Dynamics<T,Descriptor>* dynamics_ )
    : dynamics(dynamics_)
{ }

template<typename T, template<typename U> class Descriptor>
InstantiateDotDynamicsFunctional3D<T,Descriptor>::InstantiateDotDynamicsFunctional3D (
        InstantiateDotDynamicsFunctional3D<T,Descriptor> const& rhs )
    : dynamics(rhs.dynamics->clone())
{ }

template<typename T, template<typename U> class Descriptor>
InstantiateDotDynamicsFunctional3D<T,Descriptor>&
    InstantiateDotDynamicsFunctional3D<T,Descriptor>::operator= (
        InstantiateDotDynamicsFunctional3D<T,Descriptor> const& rhs )
{
    delete dynamics; dynamics = rhs.dynamics->clone();
    return *this;
}

template<typename T, template<typename U> class Descriptor>
InstantiateDotDynamicsFunctional3D<T,Descriptor>::~InstantiateDotDynamicsFunctional3D()
{
    delete dynamics;
}

template<typename T, template<typename U> class Descriptor>
void InstantiateDotDynamicsFunctional3D<T,Descriptor>::process (
        DotList3D const& dotList, BlockLattice3D<T,Descriptor>& lattice )
{
    for (plint iDot=0; iDot<dotList.getN(); ++iDot) {
        Dot3D const& dot = dotList.getDot(iDot);
        lattice.attributeDynamics(dot.x, dot.y, dot.z, dynamics->clone());
    }
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT InstantiateDotDynamicsFunctional3D<T,Descriptor>::appliesTo() const
{
    // Dynamics needs to be instantiated everywhere, including envelope.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void InstantiateDotDynamicsFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor>
InstantiateDotDynamicsFunctional3D<T,Descriptor>*
    InstantiateDotDynamicsFunctional3D<T,Descriptor>::clone() const 
{
    return new InstantiateDotDynamicsFunctional3D<T,Descriptor>(*this);
}


/* ************* Class DynamicsFromMaskFunctional3D ************************ */

template<typename T, template<typename U> class Descriptor>
DynamicsFromMaskFunctional3D<T,Descriptor>::DynamicsFromMaskFunctional3D (
        Dynamics<T,Descriptor>* dynamics_, bool whichFlag_ )
    : dynamics(dynamics_), whichFlag(whichFlag_)
{ }

template<typename T, template<typename U> class Descriptor>
DynamicsFromMaskFunctional3D<T,Descriptor>::DynamicsFromMaskFunctional3D (
        DynamicsFromMaskFunctional3D<T,Descriptor> const& rhs )
    : dynamics(rhs.dynamics->clone()),
      whichFlag(rhs.whichFlag)
{ }

template<typename T, template<typename U> class Descriptor>
DynamicsFromMaskFunctional3D<T,Descriptor>&
    DynamicsFromMaskFunctional3D<T,Descriptor>::operator= (
        DynamicsFromMaskFunctional3D<T,Descriptor> const& rhs )
{
    delete dynamics; dynamics = rhs.dynamics->clone();
    whichFlag = rhs.whichFlag;
    return *this;
}

template<typename T, template<typename U> class Descriptor>
DynamicsFromMaskFunctional3D<T,Descriptor>::~DynamicsFromMaskFunctional3D() {
    delete dynamics;
}

template<typename T, template<typename U> class Descriptor>
void DynamicsFromMaskFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                      ScalarField3D<bool>& mask )
{
    Dot3D offset = computeRelativeDisplacement(lattice, mask);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                bool flag = mask.get(iX+offset.x, iY+offset.y, iZ+offset.z);
                if ( util::boolIsEqual(flag, whichFlag) ) {
                    lattice.attributeDynamics(iX,iY,iZ, dynamics->clone());
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT DynamicsFromMaskFunctional3D<T,Descriptor>::appliesTo() const {
    // Dynamics needs to be instantiated everywhere, including envelope.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void DynamicsFromMaskFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T, template<typename U> class Descriptor>
DynamicsFromMaskFunctional3D<T,Descriptor>*
    DynamicsFromMaskFunctional3D<T,Descriptor>::clone() const 
{
    return new DynamicsFromMaskFunctional3D<T,Descriptor>(*this);
}


/* ************* Class DynamicsFromIntMaskFunctional3D ************************ */

template<typename T, template<typename U> class Descriptor>
DynamicsFromIntMaskFunctional3D<T,Descriptor>::DynamicsFromIntMaskFunctional3D (
        Dynamics<T,Descriptor>* dynamics_, int whichFlag_ )
    : dynamics(dynamics_), whichFlag(whichFlag_)
{ }

template<typename T, template<typename U> class Descriptor>
DynamicsFromIntMaskFunctional3D<T,Descriptor>::DynamicsFromIntMaskFunctional3D (
        DynamicsFromIntMaskFunctional3D<T,Descriptor> const& rhs )
    : dynamics(rhs.dynamics->clone()),
      whichFlag(rhs.whichFlag)
{ }

template<typename T, template<typename U> class Descriptor>
DynamicsFromIntMaskFunctional3D<T,Descriptor>&
    DynamicsFromIntMaskFunctional3D<T,Descriptor>::operator= (
        DynamicsFromIntMaskFunctional3D<T,Descriptor> const& rhs )
{
    delete dynamics; dynamics = rhs.dynamics->clone();
    whichFlag = rhs.whichFlag;
    return *this;
}

template<typename T, template<typename U> class Descriptor>
DynamicsFromIntMaskFunctional3D<T,Descriptor>::~DynamicsFromIntMaskFunctional3D() {
    delete dynamics;
}

template<typename T, template<typename U> class Descriptor>
void DynamicsFromIntMaskFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                      ScalarField3D<int>& mask )
{
    Dot3D offset = computeRelativeDisplacement(lattice, mask);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                int flag = mask.get(iX+offset.x, iY+offset.y, iZ+offset.z);
                if ( flag == whichFlag ) {
                    lattice.attributeDynamics(iX,iY,iZ, dynamics->clone());
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT DynamicsFromIntMaskFunctional3D<T,Descriptor>::appliesTo() const {
    // Dynamics needs to be instantiated everywhere, including envelope.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void DynamicsFromIntMaskFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T, template<typename U> class Descriptor>
DynamicsFromIntMaskFunctional3D<T,Descriptor>*
    DynamicsFromIntMaskFunctional3D<T,Descriptor>::clone() const 
{
    return new DynamicsFromIntMaskFunctional3D<T,Descriptor>(*this);
}

/* ************* Class RecomposeFromOrderZeroVariablesFunctional3D ******************* */

template<typename T, template<typename U> class Descriptor>
void RecomposeFromOrderZeroVariablesFunctional3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks )
{
    BlockLattice3D<T,Descriptor>& lattice =
        dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]);
    ScalarField3D<T> const& rhoField =
        dynamic_cast<ScalarField3D<T> const&>(*atomicBlocks[1]);
    TensorField3D<T,3> const& uField =
        dynamic_cast<TensorField3D<T,3> const&>(*atomicBlocks[2]);
    TensorField3D<T,Descriptor<T>::q> const& fNeqField =
        dynamic_cast<TensorField3D<T,Descriptor<T>::q> const&>(*atomicBlocks[3]);

    Dot3D offset1 = computeRelativeDisplacement(lattice, rhoField);
    Dot3D offset2 = computeRelativeDisplacement(lattice, uField);
    Dot3D offset3 = computeRelativeDisplacement(lattice, fNeqField);

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Cell<T,Descriptor>& cell = lattice.get(iX,iY,iZ);
                T rho               = rhoField.get(iX+offset1.x, iY+offset1.y, iZ+offset1.z);
                Array<T,3> const& u = uField.get(iX+offset2.x, iY+offset2.y, iZ+offset2.z);
                Array<T,Descriptor<T>::q> const& fNeq = fNeqField.get(iX+offset3.x, iY+offset3.y, iZ+offset3.z);

                plint recomposeOrder = 0;
                std::vector<T> rawData(cell.getDynamics().numDecomposedVariables(recomposeOrder));

                // Convert rho --> rhoBar.
                rawData[0] = Descriptor<T>::rhoBar(rho);

                // Convert u --> j
                rawData[1] = rho*u[0];
                rawData[2] = rho*u[1];
                rawData[3] = rho*u[2];
                
                for (pluint iPop = 1+Descriptor<T>::d; iPop < 1+Descriptor<T>::d+Descriptor<T>::q; ++iPop) {
                    rawData[iPop] = fNeq[iPop-(1+Descriptor<T>::d)];
                }
                
                for (pluint iPop = 1+Descriptor<T>::d+Descriptor<T>::q; iPop < rawData.size(); ++iPop) {
                    rawData[iPop] = *cell.getExternal(iPop-(1+Descriptor<T>::d+Descriptor<T>::q));
                }

                // Recompose the cell.
                cell.getDynamics().recompose(cell, rawData, recomposeOrder);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
RecomposeFromOrderZeroVariablesFunctional3D<T,Descriptor>*
    RecomposeFromOrderZeroVariablesFunctional3D<T,Descriptor>::clone() const
{
    return new RecomposeFromOrderZeroVariablesFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT RecomposeFromOrderZeroVariablesFunctional3D<T,Descriptor>::appliesTo() const
{
    // We could directly apply to the envelope too, but let's keep it
    //   bulk-only for future compatibility.
    return BlockDomain::bulk;
}

template<typename T, template<typename U> class Descriptor>
void RecomposeFromOrderZeroVariablesFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
    modified[3] = modif::nothing;
}

/* ************* Class RecomposeFromFlowVariablesFunctional3D ******************* */

template<typename T, template<typename U> class Descriptor>
void RecomposeFromFlowVariablesFunctional3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks )
{
    BlockLattice3D<T,Descriptor>& lattice =
        dynamic_cast<BlockLattice3D<T,Descriptor>&>(*atomicBlocks[0]);
    ScalarField3D<T> const& rhoField =
        dynamic_cast<ScalarField3D<T> const&>(*atomicBlocks[1]);
    TensorField3D<T,3> const& uField =
        dynamic_cast<TensorField3D<T,3> const&>(*atomicBlocks[2]);
    TensorField3D<T,6> const& SField =
        dynamic_cast<TensorField3D<T,6> const&>(*atomicBlocks[3]);

    Dot3D offset1 = computeRelativeDisplacement(lattice, rhoField);
    Dot3D offset2 = computeRelativeDisplacement(lattice, uField);
    Dot3D offset3 = computeRelativeDisplacement(lattice, SField);

    std::vector<T> rawData(10+Descriptor<T>::ExternalField::numScalars);

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Cell<T,Descriptor>& cell = lattice.get(iX,iY,iZ);
                T rho               = rhoField.get(iX+offset1.x, iY+offset1.y, iZ+offset1.z);
                Array<T,3> const& u = uField.get(iX+offset2.x, iY+offset2.y, iZ+offset2.z);
                Array<T,6> const& S = SField.get(iX+offset3.x, iY+offset3.y, iZ+offset3.z);

                // Convert rho --> rhoBar.
                rawData[0] = Descriptor<T>::rhoBar(rho);

                // Convert u --> j
                rawData[1] = rho*u[0];
                rawData[2] = rho*u[1];
                rawData[3] = rho*u[2];

                // Convert S --> PiNeq.
                T omega = cell.getDynamics().getOmega();
                T prefactor = - Descriptor<T>::cs2 * rho * (T)2 / omega;
                rawData[4] = S[0] * prefactor;
                rawData[5] = S[1] * prefactor;
                rawData[6] = S[2] * prefactor;
                rawData[7] = S[3] * prefactor;
                rawData[8] = S[4] * prefactor;
                rawData[9] = S[5] * prefactor;

                // Recompose the cell.
                plint recomposeOrder = 1;
                cell.getDynamics().recompose(cell, rawData, recomposeOrder);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
RecomposeFromFlowVariablesFunctional3D<T,Descriptor>*
    RecomposeFromFlowVariablesFunctional3D<T,Descriptor>::clone() const
{
    return new RecomposeFromFlowVariablesFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT RecomposeFromFlowVariablesFunctional3D<T,Descriptor>::appliesTo() const
{
    // We could directly apply to the envelope too, but let's keep it
    //   bulk-only for future compatibility.
    return BlockDomain::bulk;
}

template<typename T, template<typename U> class Descriptor>
void RecomposeFromFlowVariablesFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
    modified[3] = modif::nothing;
}


/* ************* Class AssignOmegaFunctional3D ******************* */

template<typename T, template<typename U> class Descriptor>
AssignOmegaFunctional3D<T,Descriptor>::AssignOmegaFunctional3D(T omega_)
    : omega(omega_)
{ }

template<typename T, template<typename U> class Descriptor>
void AssignOmegaFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    // Define dimensions of a viscosity.
    int dimDx = 2;
    int dimDt = -1;
    T scaleFactor = scaleFromReference(this->getDxScale(), dimDx,
                                       this->getDtScale(), dimDt);
    T nu_cs2 = (T)1/omega - (T)1/(T)2;
    T scaledOmega = (T)1/(scaleFactor*nu_cs2 + (T)1/(T)2);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                lattice.get(iX,iY,iZ).getDynamics().setOmega(scaledOmega);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
AssignOmegaFunctional3D<T,Descriptor>*
    AssignOmegaFunctional3D<T,Descriptor>::clone() const
{
    return new AssignOmegaFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT AssignOmegaFunctional3D<T,Descriptor>::appliesTo() const
{
    // Omega needs to be set on envelope nodes as well, because the dynamics object
    //   is being modified.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void AssignOmegaFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::staticVariables;
}


/* ************* Class AssignScalarFieldOmegaFunctional3D ******************* */

template<typename T, template<typename U> class Descriptor>
AssignScalarFieldOmegaFunctional3D<T,Descriptor>::AssignScalarFieldOmegaFunctional3D()
{ }

template<typename T, template<typename U> class Descriptor>
void AssignScalarFieldOmegaFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T> &omega )
{
    // Define dimensions of a viscosity.
    int dimDx = 2;
    int dimDt = -1;
    T scaleFactor = scaleFromReference(this->getDxScale(), dimDx,
                                       this->getDtScale(), dimDt);

    Dot3D off = computeRelativeDisplacement(lattice,omega);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        plint oX = iX + off.x;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            plint oY = iY + off.y;
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint oZ = iZ + off.z;
            
                T nu_cs2 = (T)1/omega.get(oX,oY,oZ) - (T)1/(T)2;
                T scaledOmega = (T)1/(scaleFactor*nu_cs2 + (T)1/(T)2);
                
                lattice.get(iX,iY,iZ).getDynamics().setOmega(scaledOmega);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
AssignScalarFieldOmegaFunctional3D<T,Descriptor>*
    AssignScalarFieldOmegaFunctional3D<T,Descriptor>::clone() const
{
    return new AssignScalarFieldOmegaFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT AssignScalarFieldOmegaFunctional3D<T,Descriptor>::appliesTo() const
{
    // Omega needs to be set on envelope nodes as well, because the dynamics object
    //   is being modified.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void AssignScalarFieldOmegaFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}


/* ************* Class SetConstBoundaryVelocityFunctional3D ******************* */

template<typename T, template<typename U> class Descriptor>
SetConstBoundaryVelocityFunctional3D<T,Descriptor>::SetConstBoundaryVelocityFunctional3D (
        Array<T,Descriptor<T>::d> velocity )
    : u(velocity)
{ }

template<typename T, template<typename U> class Descriptor>
void SetConstBoundaryVelocityFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    int dimDx = 1;
    int dimDt = -1;
    T scaleFactor = scaleFromReference(this->getDxScale(), dimDx,
                                       this->getDtScale(), dimDt);
    Array<T,3> scaledU = u*scaleFactor;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                lattice.get(iX,iY,iZ).defineVelocity(scaledU);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
SetConstBoundaryVelocityFunctional3D<T,Descriptor>*
    SetConstBoundaryVelocityFunctional3D<T,Descriptor>::clone() const
{
    return new SetConstBoundaryVelocityFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT SetConstBoundaryVelocityFunctional3D<T,Descriptor>::appliesTo() const
{
    // Boundary condition needs to be set on envelope nodes as well to ensure
    //   proper behavior.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void SetConstBoundaryVelocityFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::staticVariables;
}


/* ************* Class SetConstBoundaryDensityFunctional3D ******************* */

template<typename T, template<typename U> class Descriptor>
SetConstBoundaryDensityFunctional3D<T,Descriptor>::SetConstBoundaryDensityFunctional3D(T rho_)
    : rho(rho_)
{ }

template<typename T, template<typename U> class Descriptor>
void SetConstBoundaryDensityFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                lattice.get(iX,iY,iZ).defineDensity(rho);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
SetConstBoundaryDensityFunctional3D<T,Descriptor>*
    SetConstBoundaryDensityFunctional3D<T,Descriptor>::clone() const
{
    return new SetConstBoundaryDensityFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT SetConstBoundaryDensityFunctional3D<T,Descriptor>::appliesTo() const
{
    // Boundary condition needs to be set on envelope nodes as well to ensure
    //   proper behavior.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void SetConstBoundaryDensityFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
}


/* ************* Class IniConstEquilibriumFunctional3D ******************* */

template<typename T, template<typename U> class Descriptor>
IniConstEquilibriumFunctional3D<T,Descriptor>::IniConstEquilibriumFunctional3D (
        T density_, Array<T,Descriptor<T>::d> velocity, T temperature )
    : rhoBar(Descriptor<T>::rhoBar(density_)),
      j     (density_*velocity[0], density_*velocity[1], density_*velocity[2]),
      jSqr  (VectorTemplate<T,Descriptor>::normSqr(j)),
      thetaBar(temperature-(T)1)
{ }

template<typename T, template<typename U> class Descriptor>
void IniConstEquilibriumFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    int dimDx = 1;
    int dimDt = -1;
    T scaleFactor = scaleFromReference(this->getDxScale(), dimDx,
                                       this->getDtScale(), dimDt);
    Array<T,3> scaledJ = j*scaleFactor;
    T scaledJsqr = jSqr*scaleFactor*scaleFactor;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
                    lattice.get(iX,iY,iZ)[iPop] =
                        lattice.get(iX,iY,iZ).computeEquilibrium(iPop, rhoBar, scaledJ, scaledJsqr, thetaBar);
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
IniConstEquilibriumFunctional3D<T,Descriptor>*
    IniConstEquilibriumFunctional3D<T,Descriptor>::clone() const
{
    return new IniConstEquilibriumFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT IniConstEquilibriumFunctional3D<T,Descriptor>::appliesTo() const
{
    // Include boundary right away, to avoid need for envelope update.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void IniConstEquilibriumFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
}


/* ************* Class StripeOffDensityOffsetFunctional3D ******************* */

template<typename T, template<typename U> class Descriptor>
StripeOffDensityOffsetFunctional3D<T,Descriptor>::StripeOffDensityOffsetFunctional3D(T deltaRho_)
    : deltaRho(deltaRho_)
{ }

template<typename T, template<typename U> class Descriptor>
void StripeOffDensityOffsetFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Cell<T,Descriptor>& cell         = lattice.get(iX,iY,iZ);
                Dynamics<T,Descriptor>& dynamics = cell.getDynamics();
                plint orderOfDecomposition = 0;
                std::vector<T> rawData;
                dynamics.decompose(cell, rawData, orderOfDecomposition);
                T& rhoBar = rawData[0];
                rhoBar -= deltaRho;
                dynamics.recompose(cell, rawData, orderOfDecomposition);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
StripeOffDensityOffsetFunctional3D<T,Descriptor>*
    StripeOffDensityOffsetFunctional3D<T,Descriptor>::clone() const
{
    return new StripeOffDensityOffsetFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT StripeOffDensityOffsetFunctional3D<T,Descriptor>::appliesTo() const
{
    // Include boundary right away, to avoid need for envelope update.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void StripeOffDensityOffsetFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::staticVariables;
}

/* ************* Class InstantiateCompositeDynamicsFunctional3D ******************* */

template<typename T, template<typename U> class Descriptor>
InstantiateCompositeDynamicsFunctional3D<T,Descriptor>::InstantiateCompositeDynamicsFunctional3D (
        CompositeDynamics<T,Descriptor>* compositeDynamics_ )
    : compositeDynamics(compositeDynamics_)
{ }

template<typename T, template<typename U> class Descriptor>
InstantiateCompositeDynamicsFunctional3D<T,Descriptor>::InstantiateCompositeDynamicsFunctional3D (
        InstantiateCompositeDynamicsFunctional3D<T,Descriptor> const& rhs )
    : compositeDynamics(rhs.compositeDynamics->clone())
{ }

template<typename T, template<typename U> class Descriptor>
InstantiateCompositeDynamicsFunctional3D<T,Descriptor>&
    InstantiateCompositeDynamicsFunctional3D<T,Descriptor>::operator= (
        InstantiateCompositeDynamicsFunctional3D<T,Descriptor> const& rhs )
{
    delete compositeDynamics; compositeDynamics = rhs.compositeDynamics->clone();
    return *this;
}


template<typename T, template<typename U> class Descriptor>
InstantiateCompositeDynamicsFunctional3D<T,Descriptor>::~InstantiateCompositeDynamicsFunctional3D()
{
    delete compositeDynamics;
}

template<typename T, template<typename U> class Descriptor>
void InstantiateCompositeDynamicsFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                lattice.attributeDynamics(iX,iY,iZ,
                        cloneAndInsertAtTopDynamics (
                            lattice.get(iX,iY,iZ).getDynamics(),
                            compositeDynamics->clone() ) );
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void InstantiateCompositeDynamicsFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dataStructure;
}

template<typename T, template<typename U> class Descriptor>
InstantiateCompositeDynamicsFunctional3D<T,Descriptor>*
    InstantiateCompositeDynamicsFunctional3D<T,Descriptor>::clone() const 
{
    return new InstantiateCompositeDynamicsFunctional3D<T,Descriptor>(*this);
}


/* ************* Class SetExternalScalarFunctional3D ******************* */

template<typename T, template<typename U> class Descriptor>
SetExternalScalarFunctional3D<T,Descriptor>::SetExternalScalarFunctional3D (
        int whichScalar_, T externalScalar_)
    : whichScalar(whichScalar_),
      externalScalar(externalScalar_)
{
    PLB_ASSERT(whichScalar < Descriptor<T>::ExternalField::numScalars);
}

template<typename T, template<typename U> class Descriptor>
void SetExternalScalarFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                *lattice.get(iX,iY,iZ).getExternal(whichScalar) = externalScalar;
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT SetExternalScalarFunctional3D<T,Descriptor>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void SetExternalScalarFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor>
SetExternalScalarFunctional3D<T,Descriptor>*
    SetExternalScalarFunctional3D<T,Descriptor>::clone() const 
{
    return new SetExternalScalarFunctional3D<T,Descriptor>(*this);
}


/* ************* Class SetExternalScalarFromScalarFieldFunctional3D ******************* */

template<typename T, template<typename U> class Descriptor>
SetExternalScalarFromScalarFieldFunctional3D<T,Descriptor>::SetExternalScalarFromScalarFieldFunctional3D (
        int whichScalar_)
    : whichScalar(whichScalar_)
{
    PLB_ASSERT(whichScalar < Descriptor<T>::ExternalField::numScalars);
}

template<typename T, template<typename U> class Descriptor>
void SetExternalScalarFromScalarFieldFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T> &scalar)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalar);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        plint oX = iX + offset.x;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            plint oY = iY + offset.y;
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint oZ = iZ + offset.z;
                *lattice.get(iX,iY,iZ).getExternal(whichScalar) = scalar.get(oX,oY,oZ);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT SetExternalScalarFromScalarFieldFunctional3D<T,Descriptor>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void SetExternalScalarFromScalarFieldFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor>
SetExternalScalarFromScalarFieldFunctional3D<T,Descriptor>*
    SetExternalScalarFromScalarFieldFunctional3D<T,Descriptor>::clone() const 
{
    return new SetExternalScalarFromScalarFieldFunctional3D<T,Descriptor>(*this);
}


/* ************* Class MaskedSetExternalScalarFunctional3D ******************* */

template<typename T, template<typename U> class Descriptor>
MaskedSetExternalScalarFunctional3D<T,Descriptor>::MaskedSetExternalScalarFunctional3D (
        int flag_, int whichScalar_, T externalScalar_)
    : flag(flag_),
      whichScalar(whichScalar_),
      externalScalar(externalScalar_)
{
    PLB_ASSERT(whichScalar < Descriptor<T>::ExternalField::numScalars);
}

template<typename T, template<typename U> class Descriptor>
void MaskedSetExternalScalarFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<int>& mask )
{
    Dot3D offset = computeRelativeDisplacement(lattice, mask);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        plint oX = iX + offset.x;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            plint oY = iY + offset.y;
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint oZ = iZ + offset.z;
                if (mask.get(oX,oY,oZ) == flag) {
                    *lattice.get(iX,iY,iZ).getExternal(whichScalar) = externalScalar;
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT MaskedSetExternalScalarFunctional3D<T,Descriptor>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void MaskedSetExternalScalarFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;   // Lattice
    modified[1] = modif::nothing;   // Mask
}

template<typename T, template<typename U> class Descriptor>
MaskedSetExternalScalarFunctional3D<T,Descriptor>*
    MaskedSetExternalScalarFunctional3D<T,Descriptor>::clone() const 
{
    return new MaskedSetExternalScalarFunctional3D<T,Descriptor>(*this);
}


/* ************* Class SetGenericExternalScalarFunctional3D ******************* */

template<typename T, template<typename U> class Descriptor, class Functional>
SetGenericExternalScalarFunctional3D<T,Descriptor,Functional>::SetGenericExternalScalarFunctional3D (
        int whichScalar_, Functional const& functional_)
    : whichScalar(whichScalar_),
      functional(functional_)
{
    PLB_ASSERT(whichScalar < Descriptor<T>::ExternalField::numScalars);
}

template<typename T, template<typename U> class Descriptor, class Functional>
void SetGenericExternalScalarFunctional3D<T,Descriptor,Functional>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    Dot3D absOffset = lattice.getLocation();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                *lattice.get(iX,iY,iZ).getExternal(whichScalar)
                    = functional(iX+absOffset.x,iY+absOffset.y,iZ+absOffset.z);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor, class Functional>
void SetGenericExternalScalarFunctional3D<T,Descriptor,Functional>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor, class Functional>
SetGenericExternalScalarFunctional3D<T,Descriptor,Functional>*
    SetGenericExternalScalarFunctional3D<T,Descriptor,Functional>::clone() const 
{
    return new SetGenericExternalScalarFunctional3D<T,Descriptor,Functional>(*this);
}


/* ************* Class MaskedSetGenericExternalScalarFunctional3D ******************* */

template<typename T, template<typename U> class Descriptor, class Functional>
MaskedSetGenericExternalScalarFunctional3D<T,Descriptor,Functional>::MaskedSetGenericExternalScalarFunctional3D (
        int flag_, int whichScalar_, Functional const& functional_)
    : flag(flag_),
      whichScalar(whichScalar_),
      functional(functional_)
{
    PLB_ASSERT(whichScalar < Descriptor<T>::ExternalField::numScalars);
}

template<typename T, template<typename U> class Descriptor, class Functional>
void MaskedSetGenericExternalScalarFunctional3D<T,Descriptor,Functional>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<int>& mask )
{
    Dot3D offset = computeRelativeDisplacement(lattice, mask);
    Dot3D absOffset = lattice.getLocation();

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        plint oX = iX + offset.x;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            plint oY = iY + offset.y;
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint oZ = iZ + offset.z;

                if (mask.get(oX,oY,oZ) == flag) {
                    *lattice.get(iX,iY,iZ).getExternal(whichScalar)
                        = functional(iX+absOffset.x,iY+absOffset.y,iZ+absOffset.z);
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor, class Functional>
void MaskedSetGenericExternalScalarFunctional3D<T,Descriptor,Functional>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;   // Lattice.
    modified[1] = modif::nothing;           // Mask.
}

template<typename T, template<typename U> class Descriptor, class Functional>
MaskedSetGenericExternalScalarFunctional3D<T,Descriptor,Functional>*
    MaskedSetGenericExternalScalarFunctional3D<T,Descriptor,Functional>::clone() const 
{
    return new MaskedSetGenericExternalScalarFunctional3D<T,Descriptor,Functional>(*this);
}


/* ************* Class AddToExternalScalarFunctional3D ******************* */

template<typename T, template<typename U> class Descriptor>
AddToExternalScalarFunctional3D<T,Descriptor>::AddToExternalScalarFunctional3D (
        int whichScalar_, T externalScalar_)
    : whichScalar(whichScalar_),
      externalScalar(externalScalar_)
{
    PLB_ASSERT(whichScalar < Descriptor<T>::ExternalField::numScalars);
}

template<typename T, template<typename U> class Descriptor>
void AddToExternalScalarFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                *lattice.get(iX,iY,iZ).getExternal(whichScalar) += externalScalar;
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT AddToExternalScalarFunctional3D<T,Descriptor>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void AddToExternalScalarFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor>
AddToExternalScalarFunctional3D<T,Descriptor>*
    AddToExternalScalarFunctional3D<T,Descriptor>::clone() const 
{
    return new AddToExternalScalarFunctional3D<T,Descriptor>(*this);
}


/* ************* Class SetExternalVectorFunctional3D ******************* */

template<typename T, template<typename U> class Descriptor>
SetExternalVectorFunctional3D<T,Descriptor>::SetExternalVectorFunctional3D (
        int vectorStartsAt_, Array<T,Descriptor<T>::d> & externalVector_)
    : vectorStartsAt(vectorStartsAt_),
      externalVector(externalVector_)
{
    PLB_ASSERT( vectorStartsAt+Descriptor<T>::d <=
                Descriptor<T>::ExternalField::numScalars );
}

template<typename T, template<typename U> class Descriptor>
void SetExternalVectorFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                externalVector.to_cArray(lattice.get(iX,iY,iZ).getExternal(vectorStartsAt));
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT SetExternalVectorFunctional3D<T,Descriptor>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void SetExternalVectorFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor>
SetExternalVectorFunctional3D<T,Descriptor>*
    SetExternalVectorFunctional3D<T,Descriptor>::clone() const 
{
    return new SetExternalVectorFunctional3D<T,Descriptor>(*this);
}


/* ************* Class MaskedSetExternalVectorFunctional3D ******************* */

template<typename T, template<typename U> class Descriptor>
MaskedSetExternalVectorFunctional3D<T,Descriptor>::MaskedSetExternalVectorFunctional3D (
        int flag_, int vectorStartsAt_, Array<T,Descriptor<T>::d> & externalVector_)
    : flag(flag_),
      vectorStartsAt(vectorStartsAt_),
      externalVector(externalVector_)
{
    PLB_ASSERT( vectorStartsAt+Descriptor<T>::d <=
                Descriptor<T>::ExternalField::numScalars );
}

template<typename T, template<typename U> class Descriptor>
void MaskedSetExternalVectorFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<int>& mask )
{
    Dot3D offset = computeRelativeDisplacement(lattice, mask);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        plint oX = iX + offset.x;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            plint oY = iY + offset.y;
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint oZ = iZ + offset.z;
                if (mask.get(oX,oY,oZ) == flag) {
                    externalVector.to_cArray(lattice.get(iX,iY,iZ).getExternal(vectorStartsAt));
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void MaskedSetExternalVectorFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;   // Lattice.
    modified[1] = modif::nothing;           // Mask.
}

template<typename T, template<typename U> class Descriptor>
MaskedSetExternalVectorFunctional3D<T,Descriptor>*
    MaskedSetExternalVectorFunctional3D<T,Descriptor>::clone() const 
{
    return new MaskedSetExternalVectorFunctional3D<T,Descriptor>(*this);
}


/* ************* Class SetGenericExternalVectorFunctional3D ******************* */

template<typename T, template<typename U> class Descriptor, class Functional>
SetGenericExternalVectorFunctional3D<T,Descriptor,Functional>::SetGenericExternalVectorFunctional3D (
        int vectorBeginsAt_, Functional const& functional_)
    : vectorBeginsAt(vectorBeginsAt_),
      functional(functional_)
{
    PLB_ASSERT(vectorBeginsAt < Descriptor<T>::ExternalField::numScalars);
}

template<typename T, template<typename U> class Descriptor, class Functional>
void SetGenericExternalVectorFunctional3D<T,Descriptor,Functional>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    Dot3D absOffset = lattice.getLocation();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Array<T,Descriptor<T>::d> u;
                functional(iX+absOffset.x,iY+absOffset.y,iZ+absOffset.z,u);
                u.to_cArray(lattice.get(iX,iY,iZ).getExternal(vectorBeginsAt));
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor, class Functional>
void SetGenericExternalVectorFunctional3D<T,Descriptor,Functional>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor, class Functional>
SetGenericExternalVectorFunctional3D<T,Descriptor,Functional>*
    SetGenericExternalVectorFunctional3D<T,Descriptor,Functional>::clone() const 
{
    return new SetGenericExternalVectorFunctional3D<T,Descriptor,Functional>(*this);
}


/* ************* Class MaskedSetGenericExternalVectorFunctional3D ******************* */

template<typename T, template<typename U> class Descriptor, class Functional>
MaskedSetGenericExternalVectorFunctional3D<T,Descriptor,Functional>::MaskedSetGenericExternalVectorFunctional3D (
        int flag_, int vectorBeginsAt_, Functional const& functional_)
    : flag(flag_),
      vectorBeginsAt(vectorBeginsAt_),
      functional(functional_)
{
    PLB_ASSERT(vectorBeginsAt < Descriptor<T>::ExternalField::numScalars);
}

template<typename T, template<typename U> class Descriptor, class Functional>
void MaskedSetGenericExternalVectorFunctional3D<T,Descriptor,Functional>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<int>& mask )
{
    Dot3D offset = computeRelativeDisplacement(lattice, mask);
    Dot3D absOffset = lattice.getLocation();

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        plint oX = iX + offset.x;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            plint oY = iY + offset.y;
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint oZ = iZ + offset.z;

                if (mask.get(oX,oY,oZ) == flag) {
                    Array<T,Descriptor<T>::d> u;
                    functional(iX+absOffset.x,iY+absOffset.y,iZ+absOffset.z,u);
                    u.to_cArray(lattice.get(iX,iY,iZ).getExternal(vectorBeginsAt));
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor, class Functional>
void MaskedSetGenericExternalVectorFunctional3D<T,Descriptor,Functional>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;   // Lattice.
    modified[1] = modif::nothing;           // Mask.
}

template<typename T, template<typename U> class Descriptor, class Functional>
MaskedSetGenericExternalVectorFunctional3D<T,Descriptor,Functional>*
    MaskedSetGenericExternalVectorFunctional3D<T,Descriptor,Functional>::clone() const 
{
    return new MaskedSetGenericExternalVectorFunctional3D<T,Descriptor,Functional>(*this);
}


/* ************* Class SetExternalVectorFromTensorFieldFunctional3D ******************* */

template<typename T, template<typename U> class Descriptor, int nDim>
SetExternalVectorFromTensorFieldFunctional3D<T,Descriptor,nDim>::SetExternalVectorFromTensorFieldFunctional3D (
        int vectorStartsAt_)
    : vectorStartsAt(vectorStartsAt_)
{
    PLB_ASSERT( vectorStartsAt+nDim <=
                Descriptor<T>::ExternalField::numScalars );
}

template<typename T, template<typename U> class Descriptor, int nDim>
void SetExternalVectorFromTensorFieldFunctional3D<T,Descriptor,nDim>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice, TensorField3D<T,nDim>& tensor )
{
    Dot3D offset = computeRelativeDisplacement(lattice, tensor);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        plint oX = iX + offset.x;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            plint oY = iY + offset.y;
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint oZ = iZ + offset.z;
                Cell<T,Descriptor>& cell = lattice.get(iX,iY,iZ);
                Array<T,nDim> externalVector = tensor.get(oX,oY,oZ);
                
                for (plint iD=0; iD<nDim; ++iD) {
                    *cell.getExternal(vectorStartsAt+iD) = externalVector[iD];
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor, int nDim>
BlockDomain::DomainT SetExternalVectorFromTensorFieldFunctional3D<T,Descriptor,nDim>::appliesTo() const
{
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor, int nDim>
void SetExternalVectorFromTensorFieldFunctional3D<T,Descriptor,nDim>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T, template<typename U> class Descriptor, int nDim>
SetExternalVectorFromTensorFieldFunctional3D<T,Descriptor,nDim>*
    SetExternalVectorFromTensorFieldFunctional3D<T,Descriptor,nDim>::clone() const 
{
    return new SetExternalVectorFromTensorFieldFunctional3D<T,Descriptor,nDim>(*this);
}



/* ************* Class InterpolatePopulationsFunctional3D ******************* */

template<typename T, template<typename U> class Descriptor>
InterpolatePopulationsFunctional3D<T,Descriptor>::InterpolatePopulationsFunctional3D (
        plint minIter_, plint maxIter_ )
    : minIter(minIter_),
      maxIter(maxIter_)
{ }

template<typename T, template<typename U> class Descriptor>
void InterpolatePopulationsFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice1,
                      BlockLattice3D<T,Descriptor>& lattice2 )
{
    plint currentIter = (plint)lattice1.getTimeCounter().getTime();
    T lambda = (T)0.5 + 0.5*(T)(currentIter-minIter)/(T)(maxIter-minIter);
    if (lambda<(T)0.5) lambda = (T)0.5;
    if (lambda>T(1.)) lambda = T(1.);
    Dot3D offset = computeRelativeDisplacement(lattice1, lattice2);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        plint oX = iX + offset.x;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            plint oY = iY + offset.y;
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint oZ = iZ + offset.z;
                Cell<T,Descriptor>& cell1 = lattice1.get(iX,iY,iZ);
                Cell<T,Descriptor>& cell2 = lattice2.get(oX,oY,oZ);
                for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
                    cell1[iPop] = lambda*cell1[iPop] + (1.-lambda)*cell2[iPop];
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT InterpolatePopulationsFunctional3D<T,Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template<typename T, template<typename U> class Descriptor>
void InterpolatePopulationsFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T, template<typename U> class Descriptor>
InterpolatePopulationsFunctional3D<T,Descriptor>*
    InterpolatePopulationsFunctional3D<T,Descriptor>::clone() const 
{
    return new InterpolatePopulationsFunctional3D<T,Descriptor>(*this);
}


/* *************** PART II ******************************************* */
/* *************** Initialization of scalar- and tensor-fields ******* */
/* ******************************************************************* */


/* ************** Class IniConstScalarFunctional3D ******************* */

template<typename T>
IniConstScalarFunctional3D<T>::IniConstScalarFunctional3D(T value_)
    : value(value_)
{ }

template<typename T>
void IniConstScalarFunctional3D<T>::process(Box3D domain, ScalarField3D<T>& field) {
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                field.get(iX,iY,iZ) = value;
            }
        }
    }
}

template<typename T>
IniConstScalarFunctional3D<T>* IniConstScalarFunctional3D<T>::clone() const {
    return new IniConstScalarFunctional3D<T>(*this);
}

template<typename T>
BlockDomain::DomainT IniConstScalarFunctional3D<T>::appliesTo() const {
    // Include boundary right away, to avoid need for envelope update.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T>
void IniConstScalarFunctional3D<T>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
}


/* ************** Class MaskedIniConstScalarFunctional3D ******************* */

template<typename T>
MaskedIniConstScalarFunctional3D<T>::MaskedIniConstScalarFunctional3D (
        int flag_, T value_ )
    : flag(flag_),
      value(value_)
{ }

template<typename T>
void MaskedIniConstScalarFunctional3D<T>::process (
        Box3D domain, ScalarField3D<T>& field, ScalarField3D<int>& mask)
{
    Dot3D offset = computeRelativeDisplacement(field, mask);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (mask.get(iX+offset.x,iY+offset.y,iZ+offset.z)==flag) {
                    field.get(iX,iY,iZ) = value;
                }
            }
        }
    }
}

template<typename T>
MaskedIniConstScalarFunctional3D<T>* MaskedIniConstScalarFunctional3D<T>::clone() const {
    return new MaskedIniConstScalarFunctional3D<T>(*this);
}

template<typename T>
BlockDomain::DomainT MaskedIniConstScalarFunctional3D<T>::appliesTo() const {
    // Include boundary right away, to avoid need for envelope update.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T>
void MaskedIniConstScalarFunctional3D<T>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;  // Scalar-Field.
    modified[1] = modif::nothing; // Mask.
}


/* ************** Class IniConstTensorFunctional3D ******************* */

template<typename T, int nDim>
IniConstTensorFunctional3D<T,nDim>::IniConstTensorFunctional3D (
        Array<T,nDim> const& value_ )
    : value(value_)
{ }

template<typename T, int nDim>
void IniConstTensorFunctional3D<T,nDim>::process (
        Box3D domain, TensorField3D<T,nDim>& field )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                field.get(iX,iY,iZ) = value;
            }
        }
    }
}

template<typename T, int nDim>
IniConstTensorFunctional3D<T,nDim>* IniConstTensorFunctional3D<T,nDim>::clone() const {
    return new IniConstTensorFunctional3D<T,nDim>(*this);
}

template<typename T, int nDim>
BlockDomain::DomainT IniConstTensorFunctional3D<T,nDim>::appliesTo() const {
    // Include boundary right away, to avoid need for envelope update.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, int nDim>
void IniConstTensorFunctional3D<T,nDim>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
}


/* ************** Class MaskedIniConstTensorFunctional3D ******************* */

template<typename T, int nDim>
MaskedIniConstTensorFunctional3D<T,nDim>::MaskedIniConstTensorFunctional3D (
        int flag_, Array<T,nDim> const& value_ )
    : flag(flag_),
      value(value_)
{ }

template<typename T, int nDim>
void MaskedIniConstTensorFunctional3D<T,nDim>::process (
        Box3D domain,
        ScalarField3D<int>& mask, TensorField3D<T,nDim>& field )
{
    Dot3D offset = computeRelativeDisplacement(mask, field);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (mask.get(iX,iY,iZ)==flag) {
                    field.get(iX+offset.x,iY+offset.y,iZ+offset.z) = value;
                }
            }
        }
    }
}

template<typename T, int nDim>
MaskedIniConstTensorFunctional3D<T,nDim>* MaskedIniConstTensorFunctional3D<T,nDim>::clone() const {
    return new MaskedIniConstTensorFunctional3D<T,nDim>(*this);
}

template<typename T, int nDim>
BlockDomain::DomainT MaskedIniConstTensorFunctional3D<T,nDim>::appliesTo() const {
    // Include boundary right away, to avoid need for envelope update.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, int nDim>
void MaskedIniConstTensorFunctional3D<T,nDim>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing;  // Mask.
    modified[1] = modif::staticVariables;   // Tensor-Field.
}


/* ************** Class SetToCoordinateFunctional3D ****************** */

template<typename T>
SetToCoordinateFunctional3D<T>::SetToCoordinateFunctional3D(plint index_)
    : index(index_)
{
    PLB_ASSERT( index >= 0 && index <=2 );
}

template<typename T>
void SetToCoordinateFunctional3D<T>::process(Box3D domain, ScalarField3D<T>& field) {
    Dot3D relativeOffset = field.getLocation();
    Array<plint,3> ofs(relativeOffset.x, relativeOffset.y, relativeOffset.z);
    Array<plint,3> pos;
    for ( pos[0]=domain.x0; pos[0]<=domain.x1; ++pos[0] ) {
        for ( pos[1]=domain.y0; pos[1]<=domain.y1; ++pos[1] ) {
            for ( pos[2]=domain.z0; pos[2]<=domain.z1; ++pos[2] ) {
                field.get(pos[0],pos[1],pos[2]) = (T) (pos[index]+ofs[index]);
            }
        }
    }
}

template<typename T>
SetToCoordinateFunctional3D<T>* SetToCoordinateFunctional3D<T>::clone() const {
    return new SetToCoordinateFunctional3D<T>(*this);
}

template<typename T>
BlockDomain::DomainT SetToCoordinateFunctional3D<T>::appliesTo() const {
    // Boundary cannot be included, because periodic boundaries
    //   would get the wrong value.
    return BlockDomain::bulk;
}

template<typename T>
void SetToCoordinateFunctional3D<T>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
}


/* ************** Class SetToCoordinatesFunctional3D ***************** */

template<typename T>
SetToCoordinatesFunctional3D<T>::SetToCoordinatesFunctional3D()
{ }

template<typename T>
void SetToCoordinatesFunctional3D<T>::process(Box3D domain, TensorField3D<T,3>& field) {
    Dot3D relativeOffset = field.getLocation();
    Array<plint,3> ofs(relativeOffset.x, relativeOffset.y, relativeOffset.z);
    Array<plint,3> pos;
    for ( pos[0]=domain.x0; pos[0]<=domain.x1; ++pos[0] ) {
        for ( pos[1]=domain.y0; pos[1]<=domain.y1; ++pos[1] ) {
            for ( pos[2]=domain.z0; pos[2]<=domain.z1; ++pos[2] ) {
                Array<T,3>& cell = field.get(pos[0], pos[1], pos[2]);
                cell[0] = (T) (pos[0]+ofs[0]);
                cell[1] = (T) (pos[1]+ofs[1]);
                cell[2] = (T) (pos[2]+ofs[2]);
            }
        }
    }
}

template<typename T>
SetToCoordinatesFunctional3D<T>* SetToCoordinatesFunctional3D<T>::clone() const {
    return new SetToCoordinatesFunctional3D<T>(*this);
}

template<typename T>
BlockDomain::DomainT SetToCoordinatesFunctional3D<T>::appliesTo() const {
    // Boundary cannot be included, because periodic boundaries
    //   would get the wrong value.
    return BlockDomain::bulk;
}

template<typename T>
void SetToCoordinatesFunctional3D<T>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
}

/* ************** Class SetTensorComponentFunctional3D ***************** */

template<typename T, int nDim>
SetTensorComponentFunctional3D<T,nDim>::SetTensorComponentFunctional3D(int whichDim_)
    : whichDim(whichDim_)
{ }

template<typename T, int nDim>
void SetTensorComponentFunctional3D<T,nDim>::process (
        Box3D domain, ScalarField3D<T>& scalarField,
                      TensorField3D<T,nDim>& tensorField )
{
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                tensorField.get(iX+offset.x, iY+offset.y, iZ+offset.z)[whichDim] =
                    scalarField.get(iX,iY,iZ);
            }
        }
    }
}

template<typename T, int nDim>
SetTensorComponentFunctional3D<T,nDim>* SetTensorComponentFunctional3D<T,nDim>::clone() const {
    return new SetTensorComponentFunctional3D<T,nDim>(*this);
}

template<typename T, int nDim>
BlockDomain::DomainT SetTensorComponentFunctional3D<T,nDim>::appliesTo() const {
    // Boundary cannot be included, because periodic boundaries
    //   would get the wrong value.
    return BlockDomain::bulk;
}

template<typename T, int nDim>
void SetTensorComponentFunctional3D<T,nDim>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}
}  // namespace plb

#endif  // DATA_INITIALIZER_FUNCTIONAL_3D_HH

