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
 * Functionals for domain initialization -- generic implementation.
 */
#ifndef SWIG_DATA_INITIALIZER_FUNCTIONAL_3D_HH
#define SWIG_DATA_INITIALIZER_FUNCTIONAL_3D_HH

#include "plbWrapper/block/dataInitializerFunctional3D.h"
#include "core/cell.h"
#include "atomicBlock/dataField3D.h"

namespace plb {

/* ************** Class IniConstNTensorFunctional3D ******************* */

template<typename T>
IniConstNTensorFunctional3D<T>::IniConstNTensorFunctional3D (
        std::vector<T> const& value_ )
    : value(value_)
{ }

template<typename T>
void IniConstNTensorFunctional3D<T>::process (
        Box3D domain, NTensorField3D<T>& field )
{
    PLB_PRECONDITION( field.getNdim() == (plint) value.size() );
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                for (pluint iDim=0; iDim<value.size(); ++iDim) {
                    field.get(iX,iY,iZ)[iDim] = value[iDim];
                }
            }
        }
    }
}

template<typename T>
IniConstNTensorFunctional3D<T>* IniConstNTensorFunctional3D<T>::clone() const {
    return new IniConstNTensorFunctional3D<T>(*this);
}

template<typename T>
BlockDomain::DomainT IniConstNTensorFunctional3D<T>::appliesTo() const {
    // Include boundary right away, to avoid need for envelope update.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T>
void IniConstNTensorFunctional3D<T>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
}

template<typename T>
void IniConstNTensorFunctional3D<T>::rescale(double dxScale, double dtScale)
{ }


/* ************** Class MaskedIniConstNTensorFunctional3D ******************* */

template<typename T>
MaskedIniConstNTensorFunctional3D<T>::MaskedIniConstNTensorFunctional3D (
        std::vector<T> const& value_ )
    : value(value_)
{ }

template<typename T>
void MaskedIniConstNTensorFunctional3D<T>::process (
        Box3D domain,
        NTensorField3D<T>& field,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( field.getNdim() == (plint) value.size() );
    Dot3D maskOfs = computeRelativeDisplacement(field, mask);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x,iY+maskOfs.y, iZ+maskOfs.z)) {
                    for (pluint iDim=0; iDim<value.size(); ++iDim) {
                        field.get(iX,iY,iZ)[iDim] = value[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
MaskedIniConstNTensorFunctional3D<T>* MaskedIniConstNTensorFunctional3D<T>::clone() const {
    return new MaskedIniConstNTensorFunctional3D<T>(*this);
}

template<typename T>
BlockDomain::DomainT MaskedIniConstNTensorFunctional3D<T>::appliesTo() const {
    // Include boundary right away, to avoid need for envelope update.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T>
void MaskedIniConstNTensorFunctional3D<T>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
void MaskedIniConstNTensorFunctional3D<T>::rescale(double dxScale, double dtScale)
{ }


/* ************** Class SetToCoordinateFunctional3D ****************** */

template<typename T>
SetToNCoordinateFunctional3D<T>::SetToNCoordinateFunctional3D(plint index_)
    : index(index_)
{
    PLB_ASSERT( index >= 0 && index <=2 );
}

template<typename T>
void SetToNCoordinateFunctional3D<T>::process(Box3D domain, NTensorField3D<T>& field) {
    PLB_PRECONDITION( field.getNdim()==1 );
    Dot3D relativeOffset = field.getLocation();
    Array<plint,3> ofs(relativeOffset.x, relativeOffset.y, relativeOffset.z);
    Array<plint,3> pos;
    for ( pos[0]=domain.x0; pos[0]<=domain.x1; ++pos[0] ) {
        for ( pos[1]=domain.y0; pos[1]<=domain.y1; ++pos[1] ) {
            for ( pos[2]=domain.z0; pos[2]<=domain.z1; ++pos[2] ) {
                *field.get(pos[0],pos[1],pos[2]) = (T) (pos[index]+ofs[index]);
            }
        }
    }
}

template<typename T>
SetToNCoordinateFunctional3D<T>* SetToNCoordinateFunctional3D<T>::clone() const {
    return new SetToNCoordinateFunctional3D<T>(*this);
}

template<typename T>
BlockDomain::DomainT SetToNCoordinateFunctional3D<T>::appliesTo() const {
    // Boundary cannot be included, because periodic boundaries
    //   would get the wrong value.
    return BlockDomain::bulk;
}

template<typename T>
void SetToNCoordinateFunctional3D<T>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
}

template<typename T>
void SetToNCoordinateFunctional3D<T>::rescale(double dxScale, double dtScale)
{ }


/* ************** Class MaskedSetToCoordinateFunctional3D ****************** */

template<typename T>
MaskedSetToNCoordinateFunctional3D<T>::MaskedSetToNCoordinateFunctional3D(plint index_)
    : index(index_)
{
    PLB_ASSERT( index >= 0 && index <=1 );
}

template<typename T>
void MaskedSetToNCoordinateFunctional3D<T>::process (
        Box3D domain, NTensorField3D<T>& field,
        NTensorField3D<int>& mask)
{
    PLB_PRECONDITION( field.getNdim()==1 );
    Dot3D relativeOffset = field.getLocation();
    Dot3D maskOfs = computeRelativeDisplacement(field, mask);
    Array<plint,3> ofs(relativeOffset.x, relativeOffset.y, relativeOffset.z);
    Array<plint,3> pos;
    for ( pos[0]=domain.x0; pos[0]<=domain.x1; ++pos[0] ) {
        for ( pos[1]=domain.y0; pos[1]<=domain.y1; ++pos[1] ) {
            for ( pos[2]=domain.z0; pos[2]<=domain.z1; ++pos[2] ) {
                if (*mask.get(pos[0]+maskOfs.x,pos[1]+maskOfs.y,pos[2]+maskOfs.z)) {
                    *field.get(pos[0],pos[1],pos[2]) = (T) (pos[index]+ofs[index]);
                }
            }
        }
    }
}

template<typename T>
MaskedSetToNCoordinateFunctional3D<T>* MaskedSetToNCoordinateFunctional3D<T>::clone() const {
    return new MaskedSetToNCoordinateFunctional3D<T>(*this);
}

template<typename T>
BlockDomain::DomainT MaskedSetToNCoordinateFunctional3D<T>::appliesTo() const {
    // Boundary cannot be included, because periodic boundaries
    //   would get the wrong value.
    return BlockDomain::bulk;
}

template<typename T>
void MaskedSetToNCoordinateFunctional3D<T>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
void MaskedSetToNCoordinateFunctional3D<T>::rescale(double dxScale, double dtScale)
{ }


/* ************** Class SetToNCoordinatesFunctional3D ***************** */

template<typename T>
SetToNCoordinatesFunctional3D<T>::SetToNCoordinatesFunctional3D()
{ }

template<typename T>
void SetToNCoordinatesFunctional3D<T>::process(Box3D domain, NTensorField3D<T>& field) {
    PLB_PRECONDITION( field.getNdim()==3 );
    Dot3D relativeOffset = field.getLocation();
    Array<plint,3> ofs(relativeOffset.x, relativeOffset.y, relativeOffset.z);
    Array<plint,3> pos;
    for ( pos[0]=domain.x0; pos[0]<=domain.x1; ++pos[0] ) {
        for ( pos[1]=domain.y0; pos[1]<=domain.y1; ++pos[1] ) {
            for ( pos[2]=domain.z0; pos[2]<=domain.z1; ++pos[2] ) {
                T* cell = field.get(pos[0], pos[1], pos[2]);
                cell[0] = (T) (pos[0]+ofs[0]);
                cell[1] = (T) (pos[1]+ofs[1]);
                cell[2] = (T) (pos[2]+ofs[2]);
            }
        }
    }
}

template<typename T>
SetToNCoordinatesFunctional3D<T>* SetToNCoordinatesFunctional3D<T>::clone() const {
    return new SetToNCoordinatesFunctional3D<T>(*this);
}

template<typename T>
BlockDomain::DomainT SetToNCoordinatesFunctional3D<T>::appliesTo() const {
    // Boundary cannot be included, because periodic boundaries
    //   would get the wrong value.
    return BlockDomain::bulk;
}

template<typename T>
void SetToNCoordinatesFunctional3D<T>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
}

template<typename T>
void SetToNCoordinatesFunctional3D<T>::rescale(double dxScale, double dtScale)
{ }


/* ************** Class MaskedSetToNCoordinatesFunctional3D ***************** */

template<typename T>
MaskedSetToNCoordinatesFunctional3D<T>::MaskedSetToNCoordinatesFunctional3D()
{ }

template<typename T>
void MaskedSetToNCoordinatesFunctional3D<T>::process (
        Box3D domain, NTensorField3D<T>& field,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( field.getNdim()==3 );
    Dot3D maskOfs = computeRelativeDisplacement(field, mask);
    Dot3D relativeOffset = field.getLocation();
    Array<plint,3> ofs(relativeOffset.x, relativeOffset.y, relativeOffset.z);
    Array<plint,3> pos;
    for ( pos[0]=domain.x0; pos[0]<=domain.x1; ++pos[0] ) {
        for ( pos[1]=domain.y0; pos[1]<=domain.y1; ++pos[1] ) {
            for ( pos[2]=domain.z0; pos[2]<=domain.z1; ++pos[2] ) {
                if (*mask.get(pos[0]+maskOfs.x,pos[1]+maskOfs.y,pos[2]+maskOfs.z)) {
                    T* cell = field.get(pos[0], pos[1], pos[2]);
                    cell[0] = (T) (pos[0]+ofs[0]);
                    cell[1] = (T) (pos[1]+ofs[1]);
                    cell[2] = (T) (pos[2]+ofs[2]);
                }
            }
        }
    }
}

template<typename T>
MaskedSetToNCoordinatesFunctional3D<T>* MaskedSetToNCoordinatesFunctional3D<T>::clone() const {
    return new MaskedSetToNCoordinatesFunctional3D<T>(*this);
}

template<typename T>
BlockDomain::DomainT MaskedSetToNCoordinatesFunctional3D<T>::appliesTo() const {
    // Boundary cannot be included, because periodic boundaries
    //   would get the wrong value.
    return BlockDomain::bulk;
}

template<typename T>
void MaskedSetToNCoordinatesFunctional3D<T>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
void MaskedSetToNCoordinatesFunctional3D<T>::rescale(double dxScale, double dtScale)
{ }


/* ************** Class SetTensorComponentFunctional3D ***************** */

template<typename T>
SetNTensorComponentFunctional3D<T>::SetNTensorComponentFunctional3D(int whichDim_)
    : whichDim(whichDim_)
{ }

template<typename T>
void SetNTensorComponentFunctional3D<T>::process (
        Box3D domain, NTensorField3D<T>& scalarField,
                      NTensorField3D<T>& tensorField )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    PLB_PRECONDITION( tensorField.getNdim()==3 );
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                tensorField.get(iX+offset.x, iY+offset.y,iZ+offset.z)[whichDim] =
                    *scalarField.get(iX,iY,iZ);
            }
        }
    }
}

template<typename T>
SetNTensorComponentFunctional3D<T>* SetNTensorComponentFunctional3D<T>::clone() const {
    return new SetNTensorComponentFunctional3D<T>(*this);
}

template<typename T>
BlockDomain::DomainT SetNTensorComponentFunctional3D<T>::appliesTo() const {
    // Boundary cannot be included, because periodic boundaries
    //   would get the wrong value.
    return BlockDomain::bulk;
}

template<typename T>
void SetNTensorComponentFunctional3D<T>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
void SetNTensorComponentFunctional3D<T>::rescale(double dxScale, double dtScale)
{ }

/* ************** Class MaskedSetTensorComponentFunctional3D ***************** */

template<typename T>
MaskedSetNTensorComponentFunctional3D<T>::MaskedSetNTensorComponentFunctional3D(int whichDim_)
    : whichDim(whichDim_)
{ }

template<typename T>
void MaskedSetNTensorComponentFunctional3D<T>::process (
        Box3D domain, NTensorField3D<T>& scalarField,
                      NTensorField3D<T>& tensorField,
                      NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    PLB_PRECONDITION( tensorField.getNdim()==3 );
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    Dot3D maskOfs = computeRelativeDisplacement(scalarField, mask);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x,iY+maskOfs.y, iZ+maskOfs.z)) {
                    tensorField.get(iX+offset.x, iY+offset.y, iZ+offset.z)[whichDim] =
                        *scalarField.get(iX,iY,iZ);
                }
            }
        }
    }
}

template<typename T>
MaskedSetNTensorComponentFunctional3D<T>* MaskedSetNTensorComponentFunctional3D<T>::clone() const {
    return new MaskedSetNTensorComponentFunctional3D<T>(*this);
}

template<typename T>
BlockDomain::DomainT MaskedSetNTensorComponentFunctional3D<T>::appliesTo() const {
    // Boundary cannot be included, because periodic boundaries
    //   would get the wrong value.
    return BlockDomain::bulk;
}

template<typename T>
void MaskedSetNTensorComponentFunctional3D<T>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T>
void MaskedSetNTensorComponentFunctional3D<T>::rescale(double dxScale, double dtScale)
{ }


/* ************** Class AssignNTensorFunctional3D ***************** */

template<typename T>
AssignNTensorFunctional3D<T>::AssignNTensorFunctional3D()
{ }

template<typename T>
void AssignNTensorFunctional3D<T>::process (
        Box3D domain, NTensorField3D<T>& A,
                      NTensorField3D<T>& B )
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    plint ndim = A.getNdim();

    Dot3D offset = computeRelativeDisplacement(A,B);
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    *A.get(iX,iY,iZ) = *B.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        A.get(iX,iY,iZ)[iDim] = B.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
AssignNTensorFunctional3D<T>* AssignNTensorFunctional3D<T>::clone() const {
    return new AssignNTensorFunctional3D<T>(*this);
}

template<typename T>
BlockDomain::DomainT AssignNTensorFunctional3D<T>::appliesTo() const {
    // Boundary cannot be included, because periodic boundaries
    //   would get the wrong value.
    return BlockDomain::bulk;
}

template<typename T>
void AssignNTensorFunctional3D<T>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
void AssignNTensorFunctional3D<T>::rescale(double dxScale, double dtScale)
{ }


/* ************** Class MaskedAssignNTensorFunctional3D ***************** */

template<typename T>
MaskedAssignNTensorFunctional3D<T>::MaskedAssignNTensorFunctional3D()
{ }

template<typename T>
void MaskedAssignNTensorFunctional3D<T>::process (
        Box3D domain, NTensorField3D<T>& A,
                      NTensorField3D<T>& B,
                      NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    plint ndim = A.getNdim();

    Dot3D offset = computeRelativeDisplacement(A,B);
    Dot3D maskOfs = computeRelativeDisplacement(A, mask);
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x,iY+maskOfs.y, iZ+maskOfs.z)) {
                        *A.get(iX,iY,iZ) = *B.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x,iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            A.get(iX,iY,iZ)[iDim] = B.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim];
                        }
                    }
                }
            }
        }
    }
}

template<typename T>
MaskedAssignNTensorFunctional3D<T>* MaskedAssignNTensorFunctional3D<T>::clone() const {
    return new MaskedAssignNTensorFunctional3D<T>(*this);
}

template<typename T>
BlockDomain::DomainT MaskedAssignNTensorFunctional3D<T>::appliesTo() const {
    // Boundary cannot be included, because periodic boundaries
    //   would get the wrong value.
    return BlockDomain::bulk;
}

template<typename T>
void MaskedAssignNTensorFunctional3D<T>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
}

template<typename T>
void MaskedAssignNTensorFunctional3D<T>::rescale(double dxScale, double dtScale)
{ }

}  // namespace plb

#endif  // SWIG_DATA_INITIALIZER_FUNCTIONAL_3D_HH
