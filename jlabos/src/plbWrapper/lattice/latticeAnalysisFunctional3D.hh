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
 * Helper functions for domain initialization -- header file.
 */
#ifndef SWIG_LATTICE_ANALYSIS_FUNCTIONAL_3D_HH
#define SWIG_LATTICE_ANALYSIS_FUNCTIONAL_3D_HH

#include "plbWrapper/lattice/latticeAnalysisFunctional3D.h"
#include "core/blockStatistics.h"
#include "core/plbDebug.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "latticeBoltzmann/dynamicsTemplates.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "finiteDifference/fdStencils1D.h"
#include "atomicBlock/atomicBlock3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataField3D.h"
#include <cmath>

namespace plb {

template<typename T, template<typename U> class Descriptor> 
void N_BoxDensityFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice, NTensorField3D<T>& scalarField)
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                *scalarField.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                    = lattice.get(iX,iY,iZ).computeDensity();
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
N_BoxDensityFunctional3D<T,Descriptor>* N_BoxDensityFunctional3D<T,Descriptor>::clone() const
{
    return new N_BoxDensityFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void N_BoxDensityFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT N_BoxDensityFunctional3D<T,Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}


template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxDensityFunctional3D<T,Descriptor>::process (
        Box3D domain,
        BlockLattice3D<T,Descriptor>& lattice, NTensorField3D<T>& scalarField,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    Dot3D maskOfs = computeRelativeDisplacement(lattice, mask);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x,iY+maskOfs.y, iZ+maskOfs.z)) {
                    *scalarField.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                        = lattice.get(iX,iY,iZ).computeDensity();
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
Masked_N_BoxDensityFunctional3D<T,Descriptor>* Masked_N_BoxDensityFunctional3D<T,Descriptor>::clone() const
{
    return new Masked_N_BoxDensityFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxDensityFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT Masked_N_BoxDensityFunctional3D<T,Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}


template<typename T, template<typename U> class Descriptor> 
void N_BoxKineticEnergyFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice, NTensorField3D<T>& scalarField)
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Array<T,Descriptor<T>::d> velocity;
                lattice.get(iX,iY,iZ).computeVelocity(velocity);
                *scalarField.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                    = VectorTemplate<T,Descriptor>::normSqr(velocity) / (T)2;
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
N_BoxKineticEnergyFunctional3D<T,Descriptor>* N_BoxKineticEnergyFunctional3D<T,Descriptor>::clone() const
{
    return new N_BoxKineticEnergyFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void N_BoxKineticEnergyFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT N_BoxKineticEnergyFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}



template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxKineticEnergyFunctional3D<T,Descriptor>::process (
        Box3D domain,
        BlockLattice3D<T,Descriptor>& lattice, NTensorField3D<T>& scalarField,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    Dot3D maskOfs = computeRelativeDisplacement(lattice, mask);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x,iY+maskOfs.y, iZ+maskOfs.z)) {
                    Array<T,Descriptor<T>::d> velocity;
                    lattice.get(iX,iY,iZ).computeVelocity(velocity);
                    *scalarField.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                        = VectorTemplate<T,Descriptor>::normSqr(velocity) / (T)2;
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
Masked_N_BoxKineticEnergyFunctional3D<T,Descriptor>* Masked_N_BoxKineticEnergyFunctional3D<T,Descriptor>::clone() const
{
    return new Masked_N_BoxKineticEnergyFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxKineticEnergyFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT Masked_N_BoxKineticEnergyFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T, template<typename U> class Descriptor> 
void N_BoxVelocityNormFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice, NTensorField3D<T>& scalarField)
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Array<T,Descriptor<T>::d> velocity;
                lattice.get(iX,iY,iZ).computeVelocity(velocity);
                *scalarField.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                    = sqrt( VectorTemplate<T,Descriptor>::normSqr(velocity) );
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
N_BoxVelocityNormFunctional3D<T,Descriptor>* N_BoxVelocityNormFunctional3D<T,Descriptor>::clone() const
{
    return new N_BoxVelocityNormFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void N_BoxVelocityNormFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT N_BoxVelocityNormFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}




template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxVelocityNormFunctional3D<T,Descriptor>::process (
        Box3D domain,
        BlockLattice3D<T,Descriptor>& lattice, NTensorField3D<T>& scalarField,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    Dot3D maskOfs = computeRelativeDisplacement(lattice, mask);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x,iY+maskOfs.y, iZ+maskOfs.z)) {
                    Array<T,Descriptor<T>::d> velocity;
                    lattice.get(iX,iY,iZ).computeVelocity(velocity);
                    *scalarField.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                        = sqrt( VectorTemplate<T,Descriptor>::normSqr(velocity) );
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
Masked_N_BoxVelocityNormFunctional3D<T,Descriptor>* Masked_N_BoxVelocityNormFunctional3D<T,Descriptor>::clone() const
{
    return new Masked_N_BoxVelocityNormFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxVelocityNormFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT Masked_N_BoxVelocityNormFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}



template<typename T, template<typename U> class Descriptor> 
N_BoxVelocityComponentFunctional3D<T,Descriptor>::N_BoxVelocityComponentFunctional3D(int iComponent_)
    : iComponent(iComponent_)
{ }

template<typename T, template<typename U> class Descriptor> 
void N_BoxVelocityComponentFunctional3D<T,Descriptor>::process (
        Box3D domain,
        BlockLattice3D<T,Descriptor>& lattice,
        NTensorField3D<T>& scalarField )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    Array<T,Descriptor<T>::d> velocity;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                lattice.get(iX,iY,iZ).computeVelocity(velocity);
                *scalarField.get(iX+offset.x,iY+offset.y,iZ+offset.z) = velocity[iComponent];
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
N_BoxVelocityComponentFunctional3D<T,Descriptor>*
    N_BoxVelocityComponentFunctional3D<T,Descriptor>::clone() const
{
    return new N_BoxVelocityComponentFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void N_BoxVelocityComponentFunctional3D<T,Descriptor>::
    getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT N_BoxVelocityComponentFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T, template<typename U> class Descriptor> 
Masked_N_BoxVelocityComponentFunctional3D<T,Descriptor>::Masked_N_BoxVelocityComponentFunctional3D(int iComponent_)
    : iComponent(iComponent_)
{ }

template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxVelocityComponentFunctional3D<T,Descriptor>::process (
        Box3D domain,
        BlockLattice3D<T,Descriptor>& lattice,
        NTensorField3D<T>& scalarField,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    Dot3D maskOfs = computeRelativeDisplacement(lattice, mask);
    Array<T,Descriptor<T>::d> velocity;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x,iY+maskOfs.y, iZ+maskOfs.z)) {
                    lattice.get(iX,iY,iZ).computeVelocity(velocity);
                    *scalarField.get(iX+offset.x,iY+offset.y,iZ+offset.z) = velocity[iComponent];
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
Masked_N_BoxVelocityComponentFunctional3D<T,Descriptor>*
    Masked_N_BoxVelocityComponentFunctional3D<T,Descriptor>::clone() const
{
    return new Masked_N_BoxVelocityComponentFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxVelocityComponentFunctional3D<T,Descriptor>::
    getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT Masked_N_BoxVelocityComponentFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T, template<typename U> class Descriptor> 
void N_BoxVelocityFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
        NTensorField3D<T>& velocity)
{
    PLB_PRECONDITION( velocity.getNdim()==3 );
    Dot3D offset = computeRelativeDisplacement(lattice, velocity);
    Array<T,Descriptor<T>::d> velocity_value;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                lattice.get(iX,iY,iZ).computeVelocity(velocity_value);
                velocity_value.to_cArray(velocity.get(iX+offset.x,iY+offset.y,iZ+offset.z));
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
N_BoxVelocityFunctional3D<T,Descriptor>* N_BoxVelocityFunctional3D<T,Descriptor>::clone() const
{
    return new N_BoxVelocityFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void N_BoxVelocityFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT N_BoxVelocityFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}



template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxVelocityFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
        NTensorField3D<T>& velocity,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( velocity.getNdim()==3 );
    Dot3D offset = computeRelativeDisplacement(lattice, velocity);
    Dot3D maskOfs = computeRelativeDisplacement(lattice, mask);
    Array<T,Descriptor<T>::d> velocity_value;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x,iY+maskOfs.y, iZ+maskOfs.z)) {
                    lattice.get(iX,iY,iZ).computeVelocity(velocity_value);
                    velocity_value.to_cArray(velocity.get(iX+offset.x,iY+offset.y,iZ+offset.z));
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
Masked_N_BoxVelocityFunctional3D<T,Descriptor>* Masked_N_BoxVelocityFunctional3D<T,Descriptor>::clone() const
{
    return new Masked_N_BoxVelocityFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxVelocityFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT Masked_N_BoxVelocityFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T, template<typename U> class Descriptor> 
void N_BoxPiNeqFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
        NTensorField3D<T>& PiNeq )
{
    PLB_PRECONDITION( PiNeq.getNdim()==6 );
    Dot3D offset = computeRelativeDisplacement(lattice, PiNeq);
    Array<T,6> PiNeq_value;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                lattice.get(iX,iY,iZ).computePiNeq(PiNeq_value);
                PiNeq_value.to_cArray(PiNeq.get(iX+offset.x,iY+offset.y,iZ+offset.z));
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
N_BoxPiNeqFunctional3D<T,Descriptor>* N_BoxPiNeqFunctional3D<T,Descriptor>::clone() const
{
    return new N_BoxPiNeqFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void N_BoxPiNeqFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT N_BoxPiNeqFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}



template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxPiNeqFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
        NTensorField3D<T>& PiNeq,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( PiNeq.getNdim()==6 );
    Dot3D offset = computeRelativeDisplacement(lattice, PiNeq);
    Dot3D maskOfs = computeRelativeDisplacement(lattice, mask);
    Array<T,6> PiNeq_value;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x,iY+maskOfs.y, iZ+maskOfs.z)) {
                    lattice.get(iX,iY,iZ).computePiNeq(PiNeq_value);
                    PiNeq_value.to_cArray(PiNeq.get(iX+offset.x,iY+offset.y,iZ+offset.z));
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
Masked_N_BoxPiNeqFunctional3D<T,Descriptor>* Masked_N_BoxPiNeqFunctional3D<T,Descriptor>::clone() const
{
    return new Masked_N_BoxPiNeqFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxPiNeqFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT Masked_N_BoxPiNeqFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T, template<typename U> class Descriptor> 
void N_BoxShearStressFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
        NTensorField3D<T>& ShearStress )
{
    PLB_PRECONDITION( ShearStress.getNdim()==6 );
    Dot3D offset = computeRelativeDisplacement(lattice, ShearStress);
    Array<T,6> ShearStress_value;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                lattice.get(iX,iY,iZ).computeShearStress(ShearStress_value);
                ShearStress_value.to_cArray(ShearStress.get(iX+offset.x,iY+offset.y,iZ+offset.z));
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
N_BoxShearStressFunctional3D<T,Descriptor>* N_BoxShearStressFunctional3D<T,Descriptor>::clone() const
{
    return new N_BoxShearStressFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void N_BoxShearStressFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT N_BoxShearStressFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}



template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxShearStressFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
        NTensorField3D<T>& ShearStress,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( ShearStress.getNdim()==6 );
    Dot3D offset = computeRelativeDisplacement(lattice, ShearStress);
    Dot3D maskOfs = computeRelativeDisplacement(lattice, mask);
    Array<T,6> ShearStress_value;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x,iY+maskOfs.y, iZ+maskOfs.z)) {
                    lattice.get(iX,iY,iZ).computeShearStress(ShearStress_value);
                    ShearStress_value.to_cArray(ShearStress.get(iX+offset.x,iY+offset.y,iZ+offset.z));
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
Masked_N_BoxShearStressFunctional3D<T,Descriptor>* Masked_N_BoxShearStressFunctional3D<T,Descriptor>::clone() const
{
    return new Masked_N_BoxShearStressFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxShearStressFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT Masked_N_BoxShearStressFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T, template<typename U> class Descriptor> 
void N_BoxStrainRateFromStressFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
        NTensorField3D<T>& S )
{
    PLB_PRECONDITION( S.getNdim()==6 );
    Dot3D offset = computeRelativeDisplacement(lattice, S);
    Array<T,6> element;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Cell<T,Descriptor> const& cell = lattice.get(iX,iY,iZ);
                cell.computePiNeq(element);
                T omega     = cell.getDynamics().getOmega();
                T rhoBar    = cell.getDynamics().computeRhoBar(cell);
                T prefactor = - omega * Descriptor<T>::invCs2 *
                                Descriptor<T>::invRho(rhoBar) / (T)2;
                for (int iTensor=0; iTensor<SymmetricTensor<T,Descriptor>::n; ++iTensor) {
                    element[iTensor] *= prefactor;
                }
                element.to_cArray(S.get(iX+offset.x,iY+offset.y,iZ+offset.z));
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
N_BoxStrainRateFromStressFunctional3D<T,Descriptor>*
    N_BoxStrainRateFromStressFunctional3D<T,Descriptor>::clone() const
{
    return new N_BoxStrainRateFromStressFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void N_BoxStrainRateFromStressFunctional3D<T,Descriptor>::
    getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT N_BoxStrainRateFromStressFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}



template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxStrainRateFromStressFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
        NTensorField3D<T>& S,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( S.getNdim()==6 );
    Dot3D offset = computeRelativeDisplacement(lattice, S);
    Dot3D maskOfs = computeRelativeDisplacement(lattice, mask);
    Array<T,6> element;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x,iY+maskOfs.y, iZ+maskOfs.z)) {
                    Cell<T,Descriptor> const& cell = lattice.get(iX,iY,iZ);
                    cell.computePiNeq(element);
                    T omega     = cell.getDynamics().getOmega();
                    T rhoBar    = cell.getDynamics().computeRhoBar(cell);
                    T prefactor = - omega * Descriptor<T>::invCs2 *
                                    Descriptor<T>::invRho(rhoBar) / (T)2;
                    for (int iTensor=0; iTensor<SymmetricTensor<T,Descriptor>::n; ++iTensor) {
                        element[iTensor] *= prefactor;
                    }
                    element.to_cArray(S.get(iX+offset.x,iY+offset.y,iZ+offset.z));
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
Masked_N_BoxStrainRateFromStressFunctional3D<T,Descriptor>*
    Masked_N_BoxStrainRateFromStressFunctional3D<T,Descriptor>::clone() const
{
    return new Masked_N_BoxStrainRateFromStressFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxStrainRateFromStressFunctional3D<T,Descriptor>::
    getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT Masked_N_BoxStrainRateFromStressFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T, template<typename U> class Descriptor> 
N_BoxPopulationFunctional3D<T,Descriptor>::N_BoxPopulationFunctional3D(plint iComponent_)
    : iComponent(iComponent_)
{ }

template<typename T, template<typename U> class Descriptor> 
void N_BoxPopulationFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice, NTensorField3D<T>& scalarField)
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                *scalarField.get(iX+offset.x,iY+offset.y,iZ+offset.z) =
                    fullF<T,Descriptor>(lattice.get(iX,iY,iZ)[iComponent], iComponent);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
N_BoxPopulationFunctional3D<T,Descriptor>* N_BoxPopulationFunctional3D<T,Descriptor>::clone() const
{
    return new N_BoxPopulationFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void N_BoxPopulationFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT N_BoxPopulationFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T, template<typename U> class Descriptor> 
Masked_N_BoxPopulationFunctional3D<T,Descriptor>::Masked_N_BoxPopulationFunctional3D(plint iComponent_)
    : iComponent(iComponent_)
{ }

template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxPopulationFunctional3D<T,Descriptor>::process (
        Box3D domain,
        BlockLattice3D<T,Descriptor>& lattice, NTensorField3D<T>& scalarField,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    Dot3D maskOfs = computeRelativeDisplacement(lattice, mask);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x,iY+maskOfs.y, iZ+maskOfs.z)) {
                    *scalarField.get(iX+offset.x,iY+offset.y,iZ+offset.z) =
                        fullF<T,Descriptor>(lattice.get(iX,iY,iZ)[iComponent], iComponent);
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
Masked_N_BoxPopulationFunctional3D<T,Descriptor>* Masked_N_BoxPopulationFunctional3D<T,Descriptor>::clone() const
{
    return new Masked_N_BoxPopulationFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxPopulationFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT Masked_N_BoxPopulationFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T, template<typename U> class Descriptor> 
void N_BoxPopulationsFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                      NTensorField3D<T>& populations )
{
    PLB_PRECONDITION( populations.getNdim()==Descriptor<T>::q );
    Dot3D offset = computeRelativeDisplacement(lattice, populations);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Cell<T,Descriptor> const& cell = lattice.get(iX,iY,iZ);
                T* pop = populations.get(iX+offset.x,iY+offset.y,iZ+offset.z) ;
                for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
                    pop[iPop] = fullF<T,Descriptor>(cell[iPop], iPop);
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
N_BoxPopulationsFunctional3D<T,Descriptor>* N_BoxPopulationsFunctional3D<T,Descriptor>::clone() const
{
    return new N_BoxPopulationsFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void N_BoxPopulationsFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT N_BoxPopulationsFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxPopulationsFunctional3D<T,Descriptor>::process (
        Box3D domain,
        BlockLattice3D<T,Descriptor>& lattice, NTensorField3D<T>& populations,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( populations.getNdim()==Descriptor<T>::q );
    Dot3D offset = computeRelativeDisplacement(lattice, populations);
    Dot3D maskOfs = computeRelativeDisplacement(lattice, mask);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x,iY+maskOfs.y, iZ+maskOfs.z)) {
                    Cell<T,Descriptor> const& cell = lattice.get(iX,iY,iZ);
                    T* pop = populations.get(iX+offset.x,iY+offset.y,iZ+offset.z) ;
                    for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
                        pop[iPop] = fullF<T,Descriptor>(cell[iPop], iPop);
                    }
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
Masked_N_BoxPopulationsFunctional3D<T,Descriptor>* Masked_N_BoxPopulationsFunctional3D<T,Descriptor>::clone() const
{
    return new Masked_N_BoxPopulationsFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxPopulationsFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT Masked_N_BoxPopulationsFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T, template<typename U> class Descriptor> 
UPO_Rhs_Functional3D<T,Descriptor>::UPO_Rhs_Functional3D(T omega_)
    : omega(omega_)
{ }

template<typename T, template<typename U> class Descriptor> 
void UPO_Rhs_Functional3D<T,Descriptor>::process (
        Box3D domain, NTensorField3D<T>& lattice,
                      NTensorField3D<T>& result )
{
    PLB_PRECONDITION( lattice.getNdim()==Descriptor<T>::q );
    PLB_PRECONDITION( result.getNdim()==Descriptor<T>::q );
    Dot3D offset = computeRelativeDisplacement(lattice, result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                T rho;
                Array<T,3> j;
                T const* cell = lattice.get(iX,iY,iZ);
                T* resultCell = result.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                Array<T,Descriptor<T>::q> latticePop;
                latticePop.from_cArray(cell);
                // This is going to yield rho, because the populations here are without round-off optimization.
                momentTemplatesImpl<T,Descriptor<T> >::get_rhoBar_j(latticePop, rho, j);
                T rhoBar = Descriptor<T>::rhoBar(rho);
                T invRho = 1./rho;
                T jSqr = VectorTemplate<T,Descriptor>::normSqr(j);
                for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
                    plint cix = Descriptor<T>::c[iPop][0];
                    plint ciy = Descriptor<T>::c[iPop][1];
                    plint ciz = Descriptor<T>::c[iPop][2];
                    T const* behind = lattice.get(iX-cix, iY-ciy, iZ-ciz);
                    T const* infront = lattice.get(iX+cix, iY+ciy, iZ+ciz);
                    resultCell[iPop] =
                        (T)0.5 * (behind[iPop]-infront[iPop]) +
                        omega*( dynamicsTemplatesImpl<T,Descriptor<T> >::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr)
                                // Add t_i to translate from optimized to non-optimized.
                                + Descriptor<T>::t[iPop]-cell[iPop] );
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
UPO_Rhs_Functional3D<T,Descriptor>* UPO_Rhs_Functional3D<T,Descriptor>::clone() const
{
    return new UPO_Rhs_Functional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void UPO_Rhs_Functional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT UPO_Rhs_Functional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T, template<typename U> class Descriptor> 
Masked_UPO_Rhs_Functional3D<T,Descriptor>::Masked_UPO_Rhs_Functional3D(T omega_)
    : omega(omega_)
{ }

template<typename T, template<typename U> class Descriptor> 
void Masked_UPO_Rhs_Functional3D<T,Descriptor>::process (
        Box3D domain, NTensorField3D<T>& lattice,
                      NTensorField3D<T>& result,
                      NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( lattice.getNdim()==Descriptor<T>::q );
    PLB_PRECONDITION( result.getNdim()==Descriptor<T>::q );
    Dot3D offset = computeRelativeDisplacement(lattice, result);
    Dot3D maskOfs = computeRelativeDisplacement(lattice, mask);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x,iY+maskOfs.y, iZ+maskOfs.z)) {
                    T rho;
                    Array<T,3> j;
                    T const* cell = lattice.get(iX,iY,iZ);
                    T* resultCell = result.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                    Array<T,Descriptor<T>::q> latticePop;
                    latticePop.from_cArray(cell);
                    // This is going to yield rho, because the populations here are without round-off optimization.
                    momentTemplatesImpl<T,Descriptor<T> >::get_rhoBar_j(latticePop, rho, j);
                    T rhoBar = Descriptor<T>::rhoBar(rho);
                    T invRho = 1./rho;
                    T jSqr = VectorTemplate<T,Descriptor>::normSqr(j);
                    for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
                        plint cix = Descriptor<T>::c[iPop][0];
                        plint ciy = Descriptor<T>::c[iPop][1];
                        plint ciz = Descriptor<T>::c[iPop][2];
                        T const* behind = lattice.get(iX-cix, iY-ciy, iZ-ciz);
                        T const* infront = lattice.get(iX+cix, iY+ciy, iZ+ciz);
                        resultCell[iPop] =
                            (T)0.5 * (behind[iPop]-infront[iPop]) +
                            omega*( dynamicsTemplatesImpl<T,Descriptor<T> >::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr)
                                    // Add t_i to translate from optimized to non-optimized.
                                    + Descriptor<T>::t[iPop]-cell[iPop] );
                    }
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
Masked_UPO_Rhs_Functional3D<T,Descriptor>* Masked_UPO_Rhs_Functional3D<T,Descriptor>::clone() const
{
    return new Masked_UPO_Rhs_Functional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void Masked_UPO_Rhs_Functional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing; // Mask
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT Masked_UPO_Rhs_Functional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T, template<typename U> class Descriptor> 
UPO_ApplyJ_Functional3D<T,Descriptor>::UPO_ApplyJ_Functional3D(T omega_)
    : omega(omega_)
{ }

template<typename T, template<typename U> class Descriptor> 
void UPO_ApplyJ_Functional3D<T,Descriptor>::process (
        Box3D domain, std::vector<NTensorField3D<T>*> lattices )
{
    typedef Descriptor<T> D;
    typedef VectorTemplate<T,Descriptor> V;
    NTensorField3D<T> const& f = *lattices[0];
    NTensorField3D<T> const& g = *lattices[1];
    NTensorField3D<T>& h = *lattices[2];
    Dot3D gOfs = computeRelativeDisplacement(f, g);
    Dot3D hOfs = computeRelativeDisplacement(f, h);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                T const* fCell = f.get(iX,iY,iZ);
                T const* gCell = g.get(iX+gOfs.x,iY+gOfs.y,iZ+gOfs.z);
                T* hCell = h.get(iX+hOfs.x,iY+hOfs.y,iZ+hOfs.z);
                Array<T,Descriptor<T>::q> fPop, gPop;
                fPop.from_cArray(&fCell[0]);
                gPop.from_cArray(&gCell[0]);
                T fRho, gRho;
                Array<T,3> fJ, gJ;
                // This is going to yield rho, because the populations here are without round-off optimization.
                momentTemplatesImpl<T,Descriptor<T> >::get_rhoBar_j(fPop, fRho, fJ);
                momentTemplatesImpl<T,Descriptor<T> >::get_rhoBar_j(gPop, gRho, gJ);
                Array<T,3> fU = fJ/fRho;
                for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
                    plint cix = D::c[iPop][0];
                    plint ciy = D::c[iPop][1];
                    plint ciz = D::c[iPop][2];
                    T ci_rhou_g = cix*gJ[0]+ciy*gJ[1]+ciz*gJ[2];
                    T ci_u_f = cix*fU[0]+ciy*fU[1]+ciz*fU[2];
                    T rhou_g_u_f = V::scalarProduct(gJ, fU);
                    T uSqr_f = V::normSqr(fU);
                    T const* behind = g.get(iX+gOfs.x-cix, iY+gOfs.y-ciy, iZ+gOfs.z-ciz);
                    T const* infront = g.get(iX+gOfs.x+cix, iY+gOfs.y+ciy, iZ+gOfs.z+ciz);
                    T ti = D::t[iPop];
                    hCell[iPop] =
                        (T)0.5 * (behind[iPop]-infront[iPop])
                        + omega* (
                            ti*( gRho+D::invCs2*ci_rhou_g+D::invCs2*D::invCs2*(ci_rhou_g*ci_u_f-D::cs2*rhou_g_u_f)
                                 - (T)0.5*D::invCs2*D::invCs2*ci_u_f*ci_u_f*gRho
                                 + (T)0.5*D::invCs2*uSqr_f*gRho )
                            - fullF<T,Descriptor>(gCell[iPop],iPop) );
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
UPO_ApplyJ_Functional3D<T,Descriptor>* UPO_ApplyJ_Functional3D<T,Descriptor>::clone() const
{
    return new UPO_ApplyJ_Functional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void UPO_ApplyJ_Functional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;  // f
    modified[1] = modif::nothing;  // g
    modified[2] = modif::staticVariables;   // h
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT UPO_ApplyJ_Functional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T, template<typename U> class Descriptor> 
Masked_UPO_ApplyJ_Functional3D<T,Descriptor>::Masked_UPO_ApplyJ_Functional3D(T omega_)
    : omega(omega_)
{ }

template<typename T, template<typename U> class Descriptor> 
void Masked_UPO_ApplyJ_Functional3D<T,Descriptor>::process (
        Box3D domain,
        std::vector<NTensorField3D<T>*> lattices,
        NTensorField3D<int>& mask )
{
    typedef Descriptor<T> D;
    typedef VectorTemplate<T,Descriptor> V;
    NTensorField3D<T> const& f = *lattices[0];
    NTensorField3D<T> const& g = *lattices[1];
    NTensorField3D<T>& h = *lattices[2];
    Dot3D gOfs = computeRelativeDisplacement(f, g);
    Dot3D hOfs = computeRelativeDisplacement(f, h);
    Dot3D maskOfs = computeRelativeDisplacement(f, mask);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x,iY+maskOfs.y, iZ+maskOfs.z)) {
                    T const* fCell = f.get(iX,iY,iZ);
                    T const* gCell = g.get(iX+gOfs.x,iY+gOfs.y,iZ+gOfs.z);
                    T* hCell = h.get(iX+hOfs.x,iY+hOfs.y,iZ+hOfs.z);
                    Array<T,Descriptor<T>::q> fPop, gPop;
                    fPop.from_cArray(&fCell[0]);
                    gPop.from_cArray(&gCell[0]);
                    T fRho, gRho;
                    Array<T,3> fJ, gJ;
                    // This is going to yield rho, because the populations here are without round-off optimization.
                    momentTemplatesImpl<T,Descriptor<T> >::get_rhoBar_j(fPop, fRho, fJ);
                    momentTemplatesImpl<T,Descriptor<T> >::get_rhoBar_j(gPop, gRho, gJ);
                    Array<T,3> fU = fJ/fRho;
                    for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
                        plint cix = D::c[iPop][0];
                        plint ciy = D::c[iPop][1];
                        plint ciz = D::c[iPop][2];
                        T ci_rhou_g = cix*gJ[0]+ciy*gJ[1]+ciz*gJ[2];
                        T ci_u_f = cix*fU[0]+ciy*fU[1]+ciz*fU[2];
                        T rhou_g_u_f = V::scalarProduct(gJ, fU);
                        T uSqr_f = V::normSqr(fU);
                        T const* behind = g.get(iX+gOfs.x-cix, iY+gOfs.y-ciy, iZ+gOfs.z-ciz);
                        T const* infront = g.get(iX+gOfs.x+cix, iY+gOfs.y+ciy, iZ+gOfs.z+ciz);
                        T ti = D::t[iPop];
                        hCell[iPop] =
                            (T)0.5 * (behind[iPop]-infront[iPop])
                            + omega* (
                                ti*( gRho+D::invCs2*ci_rhou_g+D::invCs2*D::invCs2*(ci_rhou_g*ci_u_f-D::cs2*rhou_g_u_f)
                                     - (T)0.5*D::invCs2*D::invCs2*ci_u_f*ci_u_f*gRho
                                     + (T)0.5*D::invCs2*uSqr_f*gRho )
                                - fullF<T,Descriptor>(gCell[iPop],iPop) );
                    }
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
Masked_UPO_ApplyJ_Functional3D<T,Descriptor>* Masked_UPO_ApplyJ_Functional3D<T,Descriptor>::clone() const
{
    return new Masked_UPO_ApplyJ_Functional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void Masked_UPO_ApplyJ_Functional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;  // f
    modified[1] = modif::nothing;  // g
    modified[2] = modif::staticVariables;   // h
    modified[3] = modif::nothing;  // mask
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT Masked_UPO_ApplyJ_Functional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T, template<typename U> class Descriptor> 
void UPO_EnergyDerivative_Functional3D<T,Descriptor>::process (
        Box3D domain, NTensorField3D<T>& lattice,
                      NTensorField3D<T>& derivative )
{
    PLB_PRECONDITION( lattice.getNdim()==Descriptor<T>::q );
    PLB_PRECONDITION( derivative.getNdim()==Descriptor<T>::q );
    typedef Descriptor<T> D;
    typedef VectorTemplate<T,Descriptor> V;
    Dot3D offset = computeRelativeDisplacement(lattice, derivative);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                T rho;
                Array<T,3> j;
                T const* cell = lattice.get(iX,iY,iZ);
                T* derivativeCell = derivative.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                Array<T,Descriptor<T>::q> latticePop;
                latticePop.from_cArray(&cell[0]);
                // This is going to yield rho, because the populations here are without round-off optimization.
                momentTemplatesImpl<T,Descriptor<T> >::get_rhoBar_j(latticePop, rho, j);
                Array<T,3> u = j/rho;
                for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
                    Array<T,3> c_m_u( D::c[iPop][0] - u[0],
                                      D::c[iPop][1] - u[1],
                                      D::c[iPop][2] - u[2] );
                    derivativeCell[iPop] = V::scalarProduct(u,c_m_u);
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
UPO_EnergyDerivative_Functional3D<T,Descriptor>* UPO_EnergyDerivative_Functional3D<T,Descriptor>::clone() const
{
    return new UPO_EnergyDerivative_Functional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void UPO_EnergyDerivative_Functional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT UPO_EnergyDerivative_Functional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}



template<typename T, template<typename U> class Descriptor> 
void Masked_UPO_EnergyDerivative_Functional3D<T,Descriptor>::process (
        Box3D domain, NTensorField3D<T>& lattice,
                      NTensorField3D<T>& derivative,
                      NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( lattice.getNdim()==Descriptor<T>::q );
    PLB_PRECONDITION( derivative.getNdim()==Descriptor<T>::q );
    typedef Descriptor<T> D;
    typedef VectorTemplate<T,Descriptor> V;
    Dot3D offset = computeRelativeDisplacement(lattice, derivative);
    Dot3D maskOfs = computeRelativeDisplacement(lattice, mask);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x,iY+maskOfs.y, iZ+maskOfs.z)) {
                    T rho;
                    Array<T,3> j;
                    T const* cell = lattice.get(iX,iY,iZ);
                    T* derivativeCell = derivative.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                    Array<T,Descriptor<T>::q> latticePop;
                    latticePop.from_cArray(&cell[0]);
                    // This is going to yield rho, because the populations here are without round-off optimization.
                    momentTemplatesImpl<T,Descriptor<T> >::get_rhoBar_j(latticePop, rho, j);
                    Array<T,3> u = j/rho;
                    for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
                        Array<T,3> c_m_u( D::c[iPop][0] - u[0],
                                          D::c[iPop][1] - u[1],
                                          D::c[iPop][2] - u[2] );
                        derivativeCell[iPop] = V::scalarProduct(u,c_m_u);
                    }
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
Masked_UPO_EnergyDerivative_Functional3D<T,Descriptor>* Masked_UPO_EnergyDerivative_Functional3D<T,Descriptor>::clone() const
{
    return new Masked_UPO_EnergyDerivative_Functional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void Masked_UPO_EnergyDerivative_Functional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;  // Mask.
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT Masked_UPO_EnergyDerivative_Functional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}


}  // namespace plb

#endif  // SWIG_LATTICE_ANALYSIS_FUNCTIONAL_3D_HH
