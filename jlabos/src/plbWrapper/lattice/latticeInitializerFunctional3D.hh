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

#ifdef COMPILE_3D

/** \file
 * Functionals for domain initialization -- generic implementation.
 */
#ifndef SWIG_LATTICE_INITIALIZER_FUNCTIONAL_3D_HH
#define SWIG_LATTICE_INITIALIZER_FUNCTIONAL_3D_HH

#include "plbWrapper/block/dataInitializerFunctional3D.h"
#include "core/cell.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/roundOffPolicy.h"
#include "atomicBlock/blockLattice3D.h"

namespace plb {

/* ************* Class MaskedIniDynamicsFunctional3D ************************ */

template<typename T, template<typename U> class Descriptor>
MaskedIniDynamicsFunctional3D<T,Descriptor>::MaskedIniDynamicsFunctional3D (
        Dynamics<T,Descriptor>* dynamics_)
    : dynamics(dynamics_)
{ }

template<typename T, template<typename U> class Descriptor>
MaskedIniDynamicsFunctional3D<T,Descriptor>::MaskedIniDynamicsFunctional3D (
        MaskedIniDynamicsFunctional3D<T,Descriptor> const& rhs )
    : dynamics(rhs.dynamics->clone())
{ }

template<typename T, template<typename U> class Descriptor>
MaskedIniDynamicsFunctional3D<T,Descriptor>&
    MaskedIniDynamicsFunctional3D<T,Descriptor>::operator= (
        MaskedIniDynamicsFunctional3D<T,Descriptor> const& rhs )
{
    delete dynamics; dynamics = rhs.dynamics->clone();
    return *this;
}

template<typename T, template<typename U> class Descriptor>
MaskedIniDynamicsFunctional3D<T,Descriptor>::~MaskedIniDynamicsFunctional3D() {
    delete dynamics;
}

template<typename T, template<typename U> class Descriptor>
void MaskedIniDynamicsFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                      NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( mask.getNdim()==1 );
    Dot3D offset = computeRelativeDisplacement(lattice, mask);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+offset.x, iY+offset.y,iZ+offset.z)) {
                    lattice.attributeDynamics(iX,iY,iZ, dynamics->clone());
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT MaskedIniDynamicsFunctional3D<T,Descriptor>::appliesTo() const {
    // Dynamics needs to be instantiated everywhere, including envelope.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void MaskedIniDynamicsFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T, template<typename U> class Descriptor>
MaskedIniDynamicsFunctional3D<T,Descriptor>*
    MaskedIniDynamicsFunctional3D<T,Descriptor>::clone() const 
{
    return new MaskedIniDynamicsFunctional3D<T,Descriptor>(*this);
}

/* ************* Class N_IniBoundaryVelocityFunctional3D ************************ */

template<typename T, template<typename U> class Descriptor>
void N_IniBoundaryVelocityFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                      NTensorField3D<T>& velocity )
{
    PLB_PRECONDITION( velocity.getNdim()==3 );
    Dot3D offset = computeRelativeDisplacement(lattice, velocity);
    Array<T,3> velocityArray;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                velocityArray.from_cArray(velocity.get(iX+offset.x, iY+offset.y,iZ+offset.z));
                lattice.get(iX,iY,iZ).defineVelocity(velocityArray);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT N_IniBoundaryVelocityFunctional3D<T,Descriptor>::appliesTo() const {
    // Dynamics needs to be instantiated everywhere, including envelope.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void N_IniBoundaryVelocityFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T, template<typename U> class Descriptor>
N_IniBoundaryVelocityFunctional3D<T,Descriptor>*
    N_IniBoundaryVelocityFunctional3D<T,Descriptor>::clone() const 
{
    return new N_IniBoundaryVelocityFunctional3D<T,Descriptor>(*this);
}


/* ************* Class Masked_N_IniBoundaryVelocityFunctional3D ************************ */

template<typename T, template<typename U> class Descriptor>
void Masked_N_IniBoundaryVelocityFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                      NTensorField3D<T>& velocity,
                      NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( velocity.getNdim()==3 );
    PLB_PRECONDITION( mask.getNdim()==1 );
    Dot3D offset = computeRelativeDisplacement(lattice, velocity);
    Dot3D maskOfs = computeRelativeDisplacement(lattice, mask);
    Array<T,3> velocityArray;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x,iY+maskOfs.y, iZ+maskOfs.z)) {
                    velocityArray.from_cArray(velocity.get(iX+offset.x, iY+offset.y,iZ+offset.z));
                    lattice.get(iX,iY,iZ).defineVelocity(velocityArray);
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT Masked_N_IniBoundaryVelocityFunctional3D<T,Descriptor>::appliesTo() const {
    // Dynamics needs to be instantiated everywhere, including envelope.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void Masked_N_IniBoundaryVelocityFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
}

template<typename T, template<typename U> class Descriptor>
Masked_N_IniBoundaryVelocityFunctional3D<T,Descriptor>*
    Masked_N_IniBoundaryVelocityFunctional3D<T,Descriptor>::clone() const 
{
    return new Masked_N_IniBoundaryVelocityFunctional3D<T,Descriptor>(*this);
}

/* ************* Class N_IniEquilibriumFunctional3D ************************ */

template<typename T, template<typename U> class Descriptor>
void N_IniEquilibriumFunctional3D<T,Descriptor>::processGenericBlocks (
            Box3D domain, std::vector<AtomicBlock3D*> blocks)
{
    BlockLattice3D<T,Descriptor>& lattice =
        dynamic_cast<BlockLattice3D<T,Descriptor>&>(*blocks[0]);
    NTensorField3D<T>& density =
        dynamic_cast<NTensorField3D<T>&>(*blocks[1]);
    PLB_PRECONDITION( density.getNdim()==1 );
    NTensorField3D<T>& velocity =
        dynamic_cast<NTensorField3D<T>&>(*blocks[2]);
    PLB_PRECONDITION( velocity.getNdim()==3 );
    Dot3D densityOfs = computeRelativeDisplacement(lattice, density);
    Dot3D velocityOfs = computeRelativeDisplacement(lattice, velocity);
    Array<T,3> j;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                T rho = *density.get(iX+densityOfs.x, iY+densityOfs.y, iZ+densityOfs.z);
                T rhoBar = Descriptor<T>::rhoBar(rho);
                j.from_cArray(velocity.get(iX+velocityOfs.x, iY+velocityOfs.y, iZ+velocityOfs.z));
                j *= rho;
                T jSqr = VectorTemplate<T,Descriptor>::normSqr(j);
                for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
                    lattice.get(iX,iY,iZ)[iPop] =
                        lattice.get(iX,iY,iZ).computeEquilibrium(iPop, rhoBar, j, jSqr);
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT N_IniEquilibriumFunctional3D<T,Descriptor>::appliesTo() const {
    // Spare the effort of communication by writing directly to envelope.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void N_IniEquilibriumFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;   // Lattice
    modified[1] = modif::nothing;  // Density
    modified[2] = modif::nothing;  // Velocity
}

template<typename T, template<typename U> class Descriptor>
N_IniEquilibriumFunctional3D<T,Descriptor>*
    N_IniEquilibriumFunctional3D<T,Descriptor>::clone() const 
{
    return new N_IniEquilibriumFunctional3D<T,Descriptor>(*this);
}


/* ************* Class Masked_N_IniEquilibriumFunctional3D ************************ */

template<typename T, template<typename U> class Descriptor>
void Masked_N_IniEquilibriumFunctional3D<T,Descriptor>::processGenericBlocks (
            Box3D domain, std::vector<AtomicBlock3D*> blocks)
{
    BlockLattice3D<T,Descriptor>& lattice =
        dynamic_cast<BlockLattice3D<T,Descriptor>&>(*blocks[0]);
    NTensorField3D<T>& density =
        dynamic_cast<NTensorField3D<T>&>(*blocks[1]);
    PLB_PRECONDITION( density.getNdim()==1 );
    NTensorField3D<T>& velocity =
        dynamic_cast<NTensorField3D<T>&>(*blocks[2]);
    PLB_PRECONDITION( velocity.getNdim()==3 );
    NTensorField3D<int>& mask =
        dynamic_cast<NTensorField3D<int>&>(*blocks[3]);
    PLB_PRECONDITION( mask.getNdim()==1 );
    Dot3D densityOfs = computeRelativeDisplacement(lattice, density);
    Dot3D velocityOfs = computeRelativeDisplacement(lattice, velocity);
    Dot3D maskOfs = computeRelativeDisplacement(lattice, mask);
    Array<T,3> j;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                    T rho = *density.get(iX+densityOfs.x, iY+densityOfs.y, iZ+densityOfs.z);
                    T rhoBar = Descriptor<T>::rhoBar(rho);
                    j.from_cArray(velocity.get(iX+velocityOfs.x, iY+velocityOfs.y, iZ+velocityOfs.z));
                    j *= rho;
                    T jSqr = VectorTemplate<T,Descriptor>::normSqr(j);
                    for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
                        lattice.get(iX,iY,iZ)[iPop] =
                            lattice.get(iX,iY,iZ).computeEquilibrium(iPop, rhoBar, j, jSqr);
                    }
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT Masked_N_IniEquilibriumFunctional3D<T,Descriptor>::appliesTo() const {
    // Spare the effort of communication by writing directly to envelope.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void Masked_N_IniEquilibriumFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;   // Lattice
    modified[1] = modif::nothing;  // Density
    modified[2] = modif::nothing;  // Velocity
    modified[3] = modif::nothing;  // Mask
}

template<typename T, template<typename U> class Descriptor>
Masked_N_IniEquilibriumFunctional3D<T,Descriptor>*
    Masked_N_IniEquilibriumFunctional3D<T,Descriptor>::clone() const 
{
    return new Masked_N_IniEquilibriumFunctional3D<T,Descriptor>(*this);
}

/* ************* Class N_IniPopulationsFunctional3D ************************ */

template<typename T, template<typename U> class Descriptor>
void N_IniPopulationsFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                      NTensorField3D<T>& populations )
{
    PLB_PRECONDITION( populations.getNdim()==Descriptor<T>::q );
    Dot3D offset = computeRelativeDisplacement(lattice, populations);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Cell<T,Descriptor>& cell = lattice.get(iX,iY,iZ);
                T const* pop = populations.get(iX+offset.x, iY+offset.y,iZ+offset.z) ;
                for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
                    cell[iPop] = fBar<T,Descriptor>(pop[iPop],iPop);
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT N_IniPopulationsFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void N_IniPopulationsFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T, template<typename U> class Descriptor>
N_IniPopulationsFunctional3D<T,Descriptor>*
    N_IniPopulationsFunctional3D<T,Descriptor>::clone() const 
{
    return new N_IniPopulationsFunctional3D<T,Descriptor>(*this);
}


/* ************* Class Masked_N_IniPopulationsFunctional3D ************************ */

template<typename T, template<typename U> class Descriptor>
void Masked_N_IniPopulationsFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                      NTensorField3D<T>& populations,
                      NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( populations.getNdim()==Descriptor<T>::q );
    PLB_PRECONDITION( mask.getNdim()==1 );
    Dot3D offset = computeRelativeDisplacement(lattice, populations);
    Dot3D maskOfs= computeRelativeDisplacement(lattice, mask);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x,iY+maskOfs.y, iZ+maskOfs.z)) {
                    Cell<T,Descriptor>& cell = lattice.get(iX,iY,iZ);
                    T const* pop = populations.get(iX+offset.x, iY+offset.y,iZ+offset.z) ;
                    for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
                        cell[iPop] = fBar<T,Descriptor>(pop[iPop],iPop);
                    }
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT Masked_N_IniPopulationsFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void Masked_N_IniPopulationsFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
}

template<typename T, template<typename U> class Descriptor>
Masked_N_IniPopulationsFunctional3D<T,Descriptor>*
    Masked_N_IniPopulationsFunctional3D<T,Descriptor>::clone() const 
{
    return new Masked_N_IniPopulationsFunctional3D<T,Descriptor>(*this);
}


/* ************* Class N_IniConstPopulationsFunctional3D ************************ */

template<typename T, template<typename U> class Descriptor>
N_IniConstPopulationsFunctional3D<T,Descriptor>::N_IniConstPopulationsFunctional3D (
        std::vector<T> const& pop_ )
    : pop(pop_)
{
    PLB_ASSERT( pop.size() == (pluint) Descriptor<T>::q );
}

template<typename T, template<typename U> class Descriptor>
void N_IniConstPopulationsFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Cell<T,Descriptor>& cell = lattice.get(iX,iY,iZ);
                for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
                    cell[iPop] = fBar<T,Descriptor>(pop[iPop],iPop);
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT N_IniConstPopulationsFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void N_IniConstPopulationsFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor>
N_IniConstPopulationsFunctional3D<T,Descriptor>*
    N_IniConstPopulationsFunctional3D<T,Descriptor>::clone() const 
{
    return new N_IniConstPopulationsFunctional3D<T,Descriptor>(*this);
}


/* ************* Class Masked_N_IniConstPopulationsFunctional3D ************************ */

template<typename T, template<typename U> class Descriptor>
Masked_N_IniConstPopulationsFunctional3D<T,Descriptor>::Masked_N_IniConstPopulationsFunctional3D (
        std::vector<T> const& pop_ )
    : pop(pop_)
{
    PLB_ASSERT( pop.size() == (pluint) Descriptor<T>::q );
}

template<typename T, template<typename U> class Descriptor>
void Masked_N_IniConstPopulationsFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                      NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( mask.getNdim()==1 );
    Dot3D maskOfs= computeRelativeDisplacement(lattice, mask);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x,iY+maskOfs.y, iZ+maskOfs.z)) {
                    Cell<T,Descriptor>& cell = lattice.get(iX,iY,iZ);
                    for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
                        cell[iPop] = fBar<T,Descriptor>(pop[iPop],iPop);
                    }
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT Masked_N_IniConstPopulationsFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void Masked_N_IniConstPopulationsFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T, template<typename U> class Descriptor>
Masked_N_IniConstPopulationsFunctional3D<T,Descriptor>*
    Masked_N_IniConstPopulationsFunctional3D<T,Descriptor>::clone() const 
{
    return new Masked_N_IniConstPopulationsFunctional3D<T,Descriptor>(*this);
}

}  // namespace plb

#endif  // SWIG_LATTICE_INITIALIZER_FUNCTIONAL_3D_HH

#endif  // CoMPILE_3D
