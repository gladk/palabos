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

#ifndef ADVECTION_DIFFUSION_3D_HH
#define ADVECTION_DIFFUSION_3D_HH

#include "multiPhysics/advectionDiffusion3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "core/util.h"
#include "finiteDifference/finiteDifference3D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "multiBlock/multiDataProcessorWrapper3D.h"

namespace plb {

template< typename T, template<typename U> class TemperatureDescriptor >
void VelocityToPassiveAdvDiff3D<T,TemperatureDescriptor>::process (
        Box3D domain,
        BlockLattice3D<T,TemperatureDescriptor>& temperature,
        TensorField3D<T,3>& velocity )
{
    const int velOffset = TemperatureDescriptor<T>::ExternalField::velocityBeginsAt;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                T *u = temperature.get(iX,iY,iZ).getExternal(velOffset);
                velocity.get(iX,iY,iZ).to_cArray(u);
            }
        }
    }
}

template< typename T, template<typename U> class TemperatureDescriptor >
VelocityToPassiveAdvDiff3D<T,TemperatureDescriptor>*
    VelocityToPassiveAdvDiff3D<T,TemperatureDescriptor>::clone() const
{
    return new VelocityToPassiveAdvDiff3D<T,TemperatureDescriptor>(*this);
}

template< typename T, template<typename U> class TemperatureDescriptor >
BlockDomain::DomainT VelocityToPassiveAdvDiff3D<T,TemperatureDescriptor>::appliesTo() const {
    return BlockDomain::bulk;
}

template< typename T, template<typename U> class TemperatureDescriptor >
void VelocityToPassiveAdvDiff3D<T,TemperatureDescriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}


template< typename T, template<typename U> class TemperatureDescriptor >
void N_VelocityToPassiveAdvDiff3D<T,TemperatureDescriptor>::process (
        Box3D domain,
        BlockLattice3D<T,TemperatureDescriptor>& temperature,
        NTensorField3D<T>& velocity )
{
    PLB_PRECONDITION( velocity.getNdim()==3 );
    Dot3D offset = computeRelativeDisplacement(temperature, velocity);
    const int velOffset = TemperatureDescriptor<T>::ExternalField::velocityBeginsAt;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                T* u_to = temperature.get(iX,iY,iZ).getExternal(velOffset);
                T* u_from = velocity.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                u_to[0] = u_from[0];
                u_to[1] = u_from[1];
                u_to[2] = u_from[2];
            }
        }
    }
}

template< typename T, template<typename U> class TemperatureDescriptor >
N_VelocityToPassiveAdvDiff3D<T,TemperatureDescriptor>*
    N_VelocityToPassiveAdvDiff3D<T,TemperatureDescriptor>::clone() const
{
    return new N_VelocityToPassiveAdvDiff3D<T,TemperatureDescriptor>(*this);
}

template< typename T, template<typename U> class TemperatureDescriptor >
BlockDomain::DomainT N_VelocityToPassiveAdvDiff3D<T,TemperatureDescriptor>::appliesTo() const {
    return BlockDomain::bulk;
}

template< typename T, template<typename U> class TemperatureDescriptor >
void N_VelocityToPassiveAdvDiff3D<T,TemperatureDescriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}


template< typename T,
          template<typename U1> class FluidDescriptor, 
          template<typename U2> class ScalarDescriptor
        >
LatticeToPassiveAdvDiff3D<T,FluidDescriptor,ScalarDescriptor>::LatticeToPassiveAdvDiff3D(T scaling_)
    : scaling(scaling_)
{ }

template< typename T,
          template<typename U1> class FluidDescriptor, 
          template<typename U2> class ScalarDescriptor
        >
void LatticeToPassiveAdvDiff3D<T,FluidDescriptor,ScalarDescriptor>::process (
              Box3D domain,
              BlockLattice3D<T,FluidDescriptor>& fluid,
              BlockLattice3D<T,ScalarDescriptor>& scalar )
{
    const int velOffset = ScalarDescriptor<T>::ExternalField::velocityBeginsAt;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                T *u = scalar.get(iX,iY,iZ).getExternal(velOffset);
                Array<T,3> velocity;
                fluid.get(iX,iY,iZ).computeVelocity(velocity);
                velocity *= scaling;
                velocity.to_cArray(u);
            }
        }
    }
}

template< typename T,
          template<typename U1> class FluidDescriptor, 
          template<typename U2> class ScalarDescriptor
        >
LatticeToPassiveAdvDiff3D<T,FluidDescriptor,ScalarDescriptor>*
    LatticeToPassiveAdvDiff3D<T,FluidDescriptor,ScalarDescriptor>::clone() const
{
    return new LatticeToPassiveAdvDiff3D<T,FluidDescriptor,ScalarDescriptor>(*this);
}

template< typename T,
          template<typename U1> class FluidDescriptor, 
          template<typename U2> class ScalarDescriptor
        >
BlockDomain::DomainT LatticeToPassiveAdvDiff3D<T,FluidDescriptor,ScalarDescriptor>::appliesTo() const {
    return BlockDomain::bulk;
}

template< typename T,
          template<typename U1> class FluidDescriptor, 
          template<typename U2> class ScalarDescriptor
        >
void LatticeToPassiveAdvDiff3D<T,FluidDescriptor,ScalarDescriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}


template< typename T, template<typename U> class Descriptor >
CrystallizeAndAggregate<T,Descriptor>::CrystallizeAndAggregate(T Ncr_, T Nag_)
    : Ncr(Ncr_),
      Nag(Nag_)
{ }

template< typename T, template<typename U> class Descriptor >
void CrystallizeAndAggregate<T,Descriptor>::process (
        Box3D domain,
        BlockLattice3D<T,Descriptor>& lattice )
{
    plint bbId = BounceBack<T,Descriptor>().getId();
    Box3D extendedDomain(domain.enlarge(1));
    ScalarField3D<int> bbNodes(lattice.getNx(),lattice.getNy(),lattice.getNz());
    for (plint iX=extendedDomain.x0; iX<=extendedDomain.x1; ++iX) {
        for (plint iY=extendedDomain.y0; iY<=extendedDomain.y1; ++iY) {
            for (plint iZ=extendedDomain.z0; iZ<=extendedDomain.z1; ++iZ) {
                bbNodes.get(iX,iY,iZ) = lattice.get(iX,iY,iZ).getDynamics().getId()==bbId ? 1:0;
            }
        }
    }
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (!bbNodes.get(iX,iY,iZ)) {
                    T N = lattice.get(iX,iY,iZ).computeDensity();
                    if (N>=Ncr) {
                        lattice.attributeDynamics(iX,iY,iZ, new BounceBack<T,Descriptor>());
                    }
                    else if (N>=Nag) {
                        int numNeighbors = 0;
                        for (plint dx=-1; dx<=+1; ++dx) {
                            for (plint dy=-1; dy<=+1; ++dy) {
                                for (plint dz=-1; dz<=+1; ++dz) {
                                    numNeighbors += bbNodes.get(iX,iY,iZ);
                                }
                            }
                        }
                        if (numNeighbors>0) {
                            lattice.attributeDynamics(iX,iY,iZ, new BounceBack<T,Descriptor>());
                        }
                    }
                }
            }
        }
    }
}

template< typename T, template<typename U> class Descriptor >
CrystallizeAndAggregate<T,Descriptor>*
    CrystallizeAndAggregate<T,Descriptor>::clone() const
{
    return new CrystallizeAndAggregate<T,Descriptor>(*this);
}

template< typename T, template<typename U> class Descriptor >
BlockDomain::DomainT CrystallizeAndAggregate<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}

template< typename T, template<typename U> class Descriptor >
void CrystallizeAndAggregate<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::dynamicVariables;
}

template< typename T, template<typename U> class Descriptor >
void crystallizeAndAggregate(MultiBlockLattice3D<T,Descriptor>& lattice, T Ncr, T Nag, Box3D domain)
{
    applyProcessingFunctional (
            new CrystallizeAndAggregate<T,Descriptor>(Ncr, Nag), domain, lattice );

}

template< typename T,
          template<typename U1> class FluidDescriptor, 
          template<typename U2> class ScalarDescriptor
        >
void latticeToPassiveAdvDiff (
        MultiBlockLattice3D<T,FluidDescriptor>& fluid, MultiBlockLattice3D<T,ScalarDescriptor>& scalar, Box3D domain )
{
    applyProcessingFunctional (
            new LatticeToPassiveAdvDiff3D<T,FluidDescriptor,ScalarDescriptor>(), domain, fluid, scalar);
}


template< typename T, template<typename U> class TemperatureDescriptor >
void NVelocityToPassiveAdvDiff(MultiBlockLattice3D<T,TemperatureDescriptor>& temperature, MultiNTensorField3D<T>& velocity, Box3D domain)
{
    PLB_PRECONDITION( velocity.getNdim()==3 );
    applyProcessingFunctional (
            new N_VelocityToPassiveAdvDiff3D<T,TemperatureDescriptor>(), domain, temperature, velocity);
}

}  // namespace plb

#endif  // ADVECTION_DIFFUSION_3D_HH

