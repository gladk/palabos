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
#ifndef SWIG_LATTICE_ANALYSIS_FUNCTIONAL_2D_HH
#define SWIG_LATTICE_ANALYSIS_FUNCTIONAL_2D_HH

#include "plbWrapper/lattice/latticeAnalysisFunctional2D.h"
#include "core/blockStatistics.h"
#include "core/plbDebug.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "finiteDifference/fdStencils1D.h"
#include "atomicBlock/atomicBlock2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataField2D.h"
#include <cmath>

namespace plb {

template<typename T, template<typename U> class Descriptor> 
void N_BoxDensityFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice, NTensorField2D<T>& scalarField)
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    Dot2D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            *scalarField.get(iX+offset.x,iY+offset.y)
                = lattice.get(iX,iY).computeDensity();
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
N_BoxDensityFunctional2D<T,Descriptor>* N_BoxDensityFunctional2D<T,Descriptor>::clone() const
{
    return new N_BoxDensityFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void N_BoxDensityFunctional2D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT N_BoxDensityFunctional2D<T,Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}


template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxDensityFunctional2D<T,Descriptor>::process (
        Box2D domain,
        BlockLattice2D<T,Descriptor>& lattice, NTensorField2D<T>& scalarField,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    Dot2D offset = computeRelativeDisplacement(lattice, scalarField);
    Dot2D maskOfs = computeRelativeDisplacement(lattice, mask);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            if (*mask.get(iX+maskOfs.x,iY+maskOfs.y)) {
                *scalarField.get(iX+offset.x,iY+offset.y)
                    = lattice.get(iX,iY).computeDensity();
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
Masked_N_BoxDensityFunctional2D<T,Descriptor>* Masked_N_BoxDensityFunctional2D<T,Descriptor>::clone() const
{
    return new Masked_N_BoxDensityFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxDensityFunctional2D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT Masked_N_BoxDensityFunctional2D<T,Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}


template<typename T, template<typename U> class Descriptor> 
void N_BoxKineticEnergyFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice, NTensorField2D<T>& scalarField)
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    Dot2D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            Array<T,Descriptor<T>::d> velocity;
            lattice.get(iX,iY).computeVelocity(velocity);
            *scalarField.get(iX+offset.x,iY+offset.y)
                = VectorTemplate<T,Descriptor>::normSqr(velocity) / (T)2;
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
N_BoxKineticEnergyFunctional2D<T,Descriptor>* N_BoxKineticEnergyFunctional2D<T,Descriptor>::clone() const
{
    return new N_BoxKineticEnergyFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void N_BoxKineticEnergyFunctional2D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT N_BoxKineticEnergyFunctional2D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}



template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxKineticEnergyFunctional2D<T,Descriptor>::process (
        Box2D domain,
        BlockLattice2D<T,Descriptor>& lattice, NTensorField2D<T>& scalarField,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    Dot2D offset = computeRelativeDisplacement(lattice, scalarField);
    Dot2D maskOfs = computeRelativeDisplacement(lattice, mask);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            if (*mask.get(iX+maskOfs.x,iY+maskOfs.y)) {
                Array<T,Descriptor<T>::d> velocity;
                lattice.get(iX,iY).computeVelocity(velocity);
                *scalarField.get(iX+offset.x,iY+offset.y)
                    = VectorTemplate<T,Descriptor>::normSqr(velocity) / (T)2;
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
Masked_N_BoxKineticEnergyFunctional2D<T,Descriptor>* Masked_N_BoxKineticEnergyFunctional2D<T,Descriptor>::clone() const
{
    return new Masked_N_BoxKineticEnergyFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxKineticEnergyFunctional2D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT Masked_N_BoxKineticEnergyFunctional2D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T, template<typename U> class Descriptor> 
void N_BoxVelocityNormFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice, NTensorField2D<T>& scalarField)
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    Dot2D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            Array<T,Descriptor<T>::d> velocity;
            lattice.get(iX,iY).computeVelocity(velocity);
            *scalarField.get(iX+offset.x,iY+offset.y)
                = sqrt( VectorTemplate<T,Descriptor>::normSqr(velocity) );
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
N_BoxVelocityNormFunctional2D<T,Descriptor>* N_BoxVelocityNormFunctional2D<T,Descriptor>::clone() const
{
    return new N_BoxVelocityNormFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void N_BoxVelocityNormFunctional2D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT N_BoxVelocityNormFunctional2D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}




template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxVelocityNormFunctional2D<T,Descriptor>::process (
        Box2D domain,
        BlockLattice2D<T,Descriptor>& lattice, NTensorField2D<T>& scalarField,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    Dot2D offset = computeRelativeDisplacement(lattice, scalarField);
    Dot2D maskOfs = computeRelativeDisplacement(lattice, mask);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            if (*mask.get(iX+maskOfs.x,iY+maskOfs.y)) {
                Array<T,Descriptor<T>::d> velocity;
                lattice.get(iX,iY).computeVelocity(velocity);
                *scalarField.get(iX+offset.x,iY+offset.y)
                    = sqrt( VectorTemplate<T,Descriptor>::normSqr(velocity) );
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
Masked_N_BoxVelocityNormFunctional2D<T,Descriptor>* Masked_N_BoxVelocityNormFunctional2D<T,Descriptor>::clone() const
{
    return new Masked_N_BoxVelocityNormFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxVelocityNormFunctional2D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT Masked_N_BoxVelocityNormFunctional2D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}



template<typename T, template<typename U> class Descriptor> 
N_BoxVelocityComponentFunctional2D<T,Descriptor>::N_BoxVelocityComponentFunctional2D(int iComponent_)
    : iComponent(iComponent_)
{ }

template<typename T, template<typename U> class Descriptor> 
void N_BoxVelocityComponentFunctional2D<T,Descriptor>::process (
        Box2D domain,
        BlockLattice2D<T,Descriptor>& lattice,
        NTensorField2D<T>& scalarField )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    Dot2D offset = computeRelativeDisplacement(lattice, scalarField);
    Array<T,Descriptor<T>::d> velocity;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            lattice.get(iX,iY).computeVelocity(velocity);
            *scalarField.get(iX+offset.x,iY+offset.y) = velocity[iComponent];
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
N_BoxVelocityComponentFunctional2D<T,Descriptor>*
    N_BoxVelocityComponentFunctional2D<T,Descriptor>::clone() const
{
    return new N_BoxVelocityComponentFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void N_BoxVelocityComponentFunctional2D<T,Descriptor>::
    getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT N_BoxVelocityComponentFunctional2D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T, template<typename U> class Descriptor> 
Masked_N_BoxVelocityComponentFunctional2D<T,Descriptor>::Masked_N_BoxVelocityComponentFunctional2D(int iComponent_)
    : iComponent(iComponent_)
{ }

template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxVelocityComponentFunctional2D<T,Descriptor>::process (
        Box2D domain,
        BlockLattice2D<T,Descriptor>& lattice,
        NTensorField2D<T>& scalarField,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    Dot2D offset = computeRelativeDisplacement(lattice, scalarField);
    Dot2D maskOfs = computeRelativeDisplacement(lattice, mask);
    Array<T,Descriptor<T>::d> velocity;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            if (*mask.get(iX+maskOfs.x,iY+maskOfs.y)) {
                lattice.get(iX,iY).computeVelocity(velocity);
                *scalarField.get(iX+offset.x,iY+offset.y) = velocity[iComponent];
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
Masked_N_BoxVelocityComponentFunctional2D<T,Descriptor>*
    Masked_N_BoxVelocityComponentFunctional2D<T,Descriptor>::clone() const
{
    return new Masked_N_BoxVelocityComponentFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxVelocityComponentFunctional2D<T,Descriptor>::
    getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT Masked_N_BoxVelocityComponentFunctional2D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T, template<typename U> class Descriptor> 
void N_BoxVelocityFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
        NTensorField2D<T>& velocity)
{
    PLB_PRECONDITION( velocity.getNdim()==2 );
    Dot2D offset = computeRelativeDisplacement(lattice, velocity);
    Array<T,Descriptor<T>::d> velocity_value;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            lattice.get(iX,iY).computeVelocity(velocity_value);
            velocity_value.to_cArray(velocity.get(iX+offset.x,iY+offset.y));
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
N_BoxVelocityFunctional2D<T,Descriptor>* N_BoxVelocityFunctional2D<T,Descriptor>::clone() const
{
    return new N_BoxVelocityFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void N_BoxVelocityFunctional2D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT N_BoxVelocityFunctional2D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}



template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxVelocityFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
        NTensorField2D<T>& velocity,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( velocity.getNdim()==2 );
    Dot2D offset = computeRelativeDisplacement(lattice, velocity);
    Dot2D maskOfs = computeRelativeDisplacement(lattice, mask);
    Array<T,Descriptor<T>::d> velocity_value;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            if (*mask.get(iX+maskOfs.x,iY+maskOfs.y)) {
                lattice.get(iX,iY).computeVelocity(velocity_value);
                velocity_value.to_cArray(velocity.get(iX+offset.x,iY+offset.y));
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
Masked_N_BoxVelocityFunctional2D<T,Descriptor>* Masked_N_BoxVelocityFunctional2D<T,Descriptor>::clone() const
{
    return new Masked_N_BoxVelocityFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxVelocityFunctional2D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT Masked_N_BoxVelocityFunctional2D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T, template<typename U> class Descriptor> 
void N_BoxPiNeqFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
        NTensorField2D<T>& PiNeq )
{
    PLB_PRECONDITION( PiNeq.getNdim()==3 );
    Dot2D offset = computeRelativeDisplacement(lattice, PiNeq);
    Array<T,3> PiNeq_value;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            lattice.get(iX,iY).computePiNeq(PiNeq_value);
            PiNeq_value.to_cArray(PiNeq.get(iX+offset.x,iY+offset.y));
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
N_BoxPiNeqFunctional2D<T,Descriptor>* N_BoxPiNeqFunctional2D<T,Descriptor>::clone() const
{
    return new N_BoxPiNeqFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void N_BoxPiNeqFunctional2D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT N_BoxPiNeqFunctional2D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}



template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxPiNeqFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
        NTensorField2D<T>& PiNeq,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( PiNeq.getNdim()==3 );
    Dot2D offset = computeRelativeDisplacement(lattice, PiNeq);
    Dot2D maskOfs = computeRelativeDisplacement(lattice, mask);
    Array<T,3> PiNeq_value;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            if (*mask.get(iX+maskOfs.x,iY+maskOfs.y)) {
                lattice.get(iX,iY).computePiNeq(PiNeq_value);
                PiNeq_value.to_cArray(PiNeq.get(iX+offset.x,iY+offset.y));
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
Masked_N_BoxPiNeqFunctional2D<T,Descriptor>* Masked_N_BoxPiNeqFunctional2D<T,Descriptor>::clone() const
{
    return new Masked_N_BoxPiNeqFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxPiNeqFunctional2D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT Masked_N_BoxPiNeqFunctional2D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}

template<typename T, template<typename U> class Descriptor> 
void N_BoxShearStressFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
        NTensorField2D<T>& ShearStress )
{
    PLB_PRECONDITION( ShearStress.getNdim()==3 );
    Dot2D offset = computeRelativeDisplacement(lattice, ShearStress);
    Array<T,3> ShearStress_value;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            lattice.get(iX,iY).computeShearStress(ShearStress_value);
            ShearStress_value.to_cArray(ShearStress.get(iX+offset.x,iY+offset.y));
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
N_BoxShearStressFunctional2D<T,Descriptor>* N_BoxShearStressFunctional2D<T,Descriptor>::clone() const
{
    return new N_BoxShearStressFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void N_BoxShearStressFunctional2D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT N_BoxShearStressFunctional2D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}



template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxShearStressFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
        NTensorField2D<T>& ShearStress,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( ShearStress.getNdim()==3 );
    Dot2D offset = computeRelativeDisplacement(lattice, ShearStress);
    Dot2D maskOfs = computeRelativeDisplacement(lattice, mask);
    Array<T,3> ShearStress_value;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            if (*mask.get(iX+maskOfs.x,iY+maskOfs.y)) {
                lattice.get(iX,iY).computeShearStress(ShearStress_value);
                ShearStress_value.to_cArray(ShearStress.get(iX+offset.x,iY+offset.y));
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
Masked_N_BoxShearStressFunctional2D<T,Descriptor>* Masked_N_BoxShearStressFunctional2D<T,Descriptor>::clone() const
{
    return new Masked_N_BoxShearStressFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxShearStressFunctional2D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT Masked_N_BoxShearStressFunctional2D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}



template<typename T, template<typename U> class Descriptor> 
void N_BoxStrainRateFromStressFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
        NTensorField2D<T>& S )
{
    PLB_PRECONDITION( S.getNdim()==3 );
    Dot2D offset = computeRelativeDisplacement(lattice, S);
    Array<T,3> element;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            Cell<T,Descriptor> const& cell = lattice.get(iX,iY);
            cell.computePiNeq(element);
            T omega     = cell.getDynamics().getOmega();
            T rhoBar    = cell.getDynamics().computeRhoBar(cell);
            T prefactor = - omega * Descriptor<T>::invCs2 *
                            Descriptor<T>::invRho(rhoBar) / (T)2;
            for (int iTensor=0; iTensor<SymmetricTensor<T,Descriptor>::n; ++iTensor) {
                element[iTensor] *= prefactor;
            }
            element.to_cArray(S.get(iX+offset.x,iY+offset.y));
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
N_BoxStrainRateFromStressFunctional2D<T,Descriptor>*
    N_BoxStrainRateFromStressFunctional2D<T,Descriptor>::clone() const
{
    return new N_BoxStrainRateFromStressFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void N_BoxStrainRateFromStressFunctional2D<T,Descriptor>::
    getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT N_BoxStrainRateFromStressFunctional2D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}



template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxStrainRateFromStressFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
        NTensorField2D<T>& S,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( S.getNdim()==3 );
    Dot2D offset = computeRelativeDisplacement(lattice, S);
    Dot2D maskOfs = computeRelativeDisplacement(lattice, mask);
    Array<T,3> element;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            if (*mask.get(iX+maskOfs.x,iY+maskOfs.y)) {
                Cell<T,Descriptor> const& cell = lattice.get(iX,iY);
                cell.computePiNeq(element);
                T omega     = cell.getDynamics().getOmega();
                T rhoBar    = cell.getDynamics().computeRhoBar(cell);
                T prefactor = - omega * Descriptor<T>::invCs2 *
                                Descriptor<T>::invRho(rhoBar) / (T)2;
                for (int iTensor=0; iTensor<SymmetricTensor<T,Descriptor>::n; ++iTensor) {
                    element[iTensor] *= prefactor;
                }
                element.to_cArray(S.get(iX+offset.x,iY+offset.y));
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
Masked_N_BoxStrainRateFromStressFunctional2D<T,Descriptor>*
    Masked_N_BoxStrainRateFromStressFunctional2D<T,Descriptor>::clone() const
{
    return new Masked_N_BoxStrainRateFromStressFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxStrainRateFromStressFunctional2D<T,Descriptor>::
    getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT Masked_N_BoxStrainRateFromStressFunctional2D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T, template<typename U> class Descriptor> 
N_BoxPopulationFunctional2D<T,Descriptor>::N_BoxPopulationFunctional2D(plint iComponent_)
    : iComponent(iComponent_)
{ }

template<typename T, template<typename U> class Descriptor> 
void N_BoxPopulationFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice, NTensorField2D<T>& scalarField)
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    Dot2D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            *scalarField.get(iX+offset.x,iY+offset.y) =
                fullF<T,Descriptor>(lattice.get(iX,iY)[iComponent], iComponent);
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
N_BoxPopulationFunctional2D<T,Descriptor>* N_BoxPopulationFunctional2D<T,Descriptor>::clone() const
{
    return new N_BoxPopulationFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void N_BoxPopulationFunctional2D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT N_BoxPopulationFunctional2D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T, template<typename U> class Descriptor> 
Masked_N_BoxPopulationFunctional2D<T,Descriptor>::Masked_N_BoxPopulationFunctional2D(plint iComponent_)
    : iComponent(iComponent_)
{ }

template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxPopulationFunctional2D<T,Descriptor>::process (
        Box2D domain,
        BlockLattice2D<T,Descriptor>& lattice, NTensorField2D<T>& scalarField,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    Dot2D offset = computeRelativeDisplacement(lattice, scalarField);
    Dot2D maskOfs = computeRelativeDisplacement(lattice, mask);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            if (*mask.get(iX+maskOfs.x,iY+maskOfs.y)) {
                *scalarField.get(iX+offset.x,iY+offset.y) =
                    fullF<T,Descriptor>(lattice.get(iX,iY)[iComponent], iComponent);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
Masked_N_BoxPopulationFunctional2D<T,Descriptor>* Masked_N_BoxPopulationFunctional2D<T,Descriptor>::clone() const
{
    return new Masked_N_BoxPopulationFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxPopulationFunctional2D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT Masked_N_BoxPopulationFunctional2D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T, template<typename U> class Descriptor> 
void N_BoxPopulationsFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
                      NTensorField2D<T>& populations )
{
    PLB_PRECONDITION( populations.getNdim()==Descriptor<T>::q );
    Dot2D offset = computeRelativeDisplacement(lattice, populations);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            Cell<T,Descriptor> const& cell = lattice.get(iX,iY);
            T* pop = populations.get(iX+offset.x,iY+offset.y) ;
            for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
                pop[iPop] = fullF<T,Descriptor>(cell[iPop], iPop);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
N_BoxPopulationsFunctional2D<T,Descriptor>* N_BoxPopulationsFunctional2D<T,Descriptor>::clone() const
{
    return new N_BoxPopulationsFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void N_BoxPopulationsFunctional2D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT N_BoxPopulationsFunctional2D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxPopulationsFunctional2D<T,Descriptor>::process (
        Box2D domain,
        BlockLattice2D<T,Descriptor>& lattice, NTensorField2D<T>& populations,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( populations.getNdim()==Descriptor<T>::q );
    Dot2D offset = computeRelativeDisplacement(lattice, populations);
    Dot2D maskOfs = computeRelativeDisplacement(lattice, mask);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            if (*mask.get(iX+maskOfs.x,iY+maskOfs.y)) {
                Cell<T,Descriptor> const& cell = lattice.get(iX,iY);
                T* pop = populations.get(iX+offset.x,iY+offset.y) ;
                for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
                    pop[iPop] = fullF<T,Descriptor>(cell[iPop], iPop);
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
Masked_N_BoxPopulationsFunctional2D<T,Descriptor>* Masked_N_BoxPopulationsFunctional2D<T,Descriptor>::clone() const
{
    return new Masked_N_BoxPopulationsFunctional2D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void Masked_N_BoxPopulationsFunctional2D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT Masked_N_BoxPopulationsFunctional2D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulk;
}

}  // namespace plb

#endif  // SWIG_LATTICE_ANALYSIS_FUNCTIONAL_2D_HH
