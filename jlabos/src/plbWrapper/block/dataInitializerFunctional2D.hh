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
#ifndef SWIG_DATA_INITIALIZER_FUNCTIONAL_2D_HH
#define SWIG_DATA_INITIALIZER_FUNCTIONAL_2D_HH

#include "plbWrapper/block/dataInitializerFunctional2D.h"
#include "core/cell.h"
#include "atomicBlock/dataField2D.h"

namespace plb {

/* ************** Class IniConstNTensorFunctional2D ******************* */

template<typename T>
IniConstNTensorFunctional2D<T>::IniConstNTensorFunctional2D (
        std::vector<T> const& value_ )
    : value(value_)
{ }

template<typename T>
void IniConstNTensorFunctional2D<T>::process (
        Box2D domain, NTensorField2D<T>& field )
{
    PLB_PRECONDITION( field.getNdim() == (plint) value.size() );
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (pluint iDim=0; iDim<value.size(); ++iDim) {
                field.get(iX,iY)[iDim] = value[iDim];
            }
        }
    }
}

template<typename T>
IniConstNTensorFunctional2D<T>* IniConstNTensorFunctional2D<T>::clone() const {
    return new IniConstNTensorFunctional2D<T>(*this);
}

template<typename T>
BlockDomain::DomainT IniConstNTensorFunctional2D<T>::appliesTo() const {
    // Include boundary right away, to avoid need for envelope update.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T>
void IniConstNTensorFunctional2D<T>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
}

template<typename T>
void IniConstNTensorFunctional2D<T>::rescale(double dxScale, double dtScale)
{ }


/* ************** Class MaskedIniConstNTensorFunctional2D ******************* */

template<typename T>
MaskedIniConstNTensorFunctional2D<T>::MaskedIniConstNTensorFunctional2D (
        std::vector<T> const& value_ )
    : value(value_)
{ }

template<typename T>
void MaskedIniConstNTensorFunctional2D<T>::process (
        Box2D domain,
        NTensorField2D<T>& field,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( field.getNdim() == (plint) value.size() );
    Dot2D maskOfs = computeRelativeDisplacement(field, mask);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            if (*mask.get(iX+maskOfs.x,iY+maskOfs.y)) {
                for (pluint iDim=0; iDim<value.size(); ++iDim) {
                    field.get(iX,iY)[iDim] = value[iDim];
                }
            }
        }
    }
}

template<typename T>
MaskedIniConstNTensorFunctional2D<T>* MaskedIniConstNTensorFunctional2D<T>::clone() const {
    return new MaskedIniConstNTensorFunctional2D<T>(*this);
}

template<typename T>
BlockDomain::DomainT MaskedIniConstNTensorFunctional2D<T>::appliesTo() const {
    // Include boundary right away, to avoid need for envelope update.
    return BlockDomain::bulkAndEnvelope;
}

template<typename T>
void MaskedIniConstNTensorFunctional2D<T>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
void MaskedIniConstNTensorFunctional2D<T>::rescale(double dxScale, double dtScale)
{ }


/* ************** Class SetToCoordinateFunctional2D ****************** */

template<typename T>
SetToNCoordinateFunctional2D<T>::SetToNCoordinateFunctional2D(plint index_)
    : index(index_)
{
    PLB_ASSERT( index >= 0 && index <=1 );
}

template<typename T>
void SetToNCoordinateFunctional2D<T>::process(Box2D domain, NTensorField2D<T>& field) {
    PLB_PRECONDITION( field.getNdim()==1 );
    Dot2D relativeOffset = field.getLocation();
    Array<plint,2> ofs(relativeOffset.x, relativeOffset.y);
    Array<plint,2> pos;
    for ( pos[0]=domain.x0; pos[0]<=domain.x1; ++pos[0] ) {
        for ( pos[1]=domain.y0; pos[1]<=domain.y1; ++pos[1] ) {
            *field.get(pos[0],pos[1]) = (T) (pos[index]+ofs[index]);
        }
    }
}

template<typename T>
SetToNCoordinateFunctional2D<T>* SetToNCoordinateFunctional2D<T>::clone() const {
    return new SetToNCoordinateFunctional2D<T>(*this);
}

template<typename T>
BlockDomain::DomainT SetToNCoordinateFunctional2D<T>::appliesTo() const {
    // Boundary cannot be included, because periodic boundaries
    //   would get the wrong value.
    return BlockDomain::bulk;
}

template<typename T>
void SetToNCoordinateFunctional2D<T>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
}

template<typename T>
void SetToNCoordinateFunctional2D<T>::rescale(double dxScale, double dtScale)
{ }


/* ************** Class MaskedSetToCoordinateFunctional2D ****************** */

template<typename T>
MaskedSetToNCoordinateFunctional2D<T>::MaskedSetToNCoordinateFunctional2D(plint index_)
    : index(index_)
{
    PLB_ASSERT( index >= 0 && index <=1 );
}

template<typename T>
void MaskedSetToNCoordinateFunctional2D<T>::process (
        Box2D domain, NTensorField2D<T>& field,
        NTensorField2D<int>& mask)
{
    PLB_PRECONDITION( field.getNdim()==1 );
    Dot2D relativeOffset = field.getLocation();
    Dot2D maskOfs = computeRelativeDisplacement(field, mask);
    Array<plint,2> ofs(relativeOffset.x, relativeOffset.y);
    Array<plint,2> pos;
    for ( pos[0]=domain.x0; pos[0]<=domain.x1; ++pos[0] ) {
        for ( pos[1]=domain.y0; pos[1]<=domain.y1; ++pos[1] ) {
            if (*mask.get(pos[0]+maskOfs.x,pos[1]+maskOfs.y)) {
                *field.get(pos[0],pos[1]) = (T) (pos[index]+ofs[index]);
            }
        }
    }
}

template<typename T>
MaskedSetToNCoordinateFunctional2D<T>* MaskedSetToNCoordinateFunctional2D<T>::clone() const {
    return new MaskedSetToNCoordinateFunctional2D<T>(*this);
}

template<typename T>
BlockDomain::DomainT MaskedSetToNCoordinateFunctional2D<T>::appliesTo() const {
    // Boundary cannot be included, because periodic boundaries
    //   would get the wrong value.
    return BlockDomain::bulk;
}

template<typename T>
void MaskedSetToNCoordinateFunctional2D<T>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
void MaskedSetToNCoordinateFunctional2D<T>::rescale(double dxScale, double dtScale)
{ }


/* ************** Class SetToNCoordinatesFunctional2D ***************** */

template<typename T>
SetToNCoordinatesFunctional2D<T>::SetToNCoordinatesFunctional2D()
{ }

template<typename T>
void SetToNCoordinatesFunctional2D<T>::process(Box2D domain, NTensorField2D<T>& field) {
    PLB_PRECONDITION( field.getNdim()==2 );
    Dot2D relativeOffset = field.getLocation();
    Array<plint,2> ofs(relativeOffset.x, relativeOffset.y);
    Array<plint,2> pos;
    for ( pos[0]=domain.x0; pos[0]<=domain.x1; ++pos[0] ) {
        for ( pos[1]=domain.y0; pos[1]<=domain.y1; ++pos[1] ) {
            T* cell = field.get(pos[0], pos[1]);
            cell[0] = (T) (pos[0]+ofs[0]);
            cell[1] = (T) (pos[1]+ofs[1]);
        }
    }
}

template<typename T>
SetToNCoordinatesFunctional2D<T>* SetToNCoordinatesFunctional2D<T>::clone() const {
    return new SetToNCoordinatesFunctional2D<T>(*this);
}

template<typename T>
BlockDomain::DomainT SetToNCoordinatesFunctional2D<T>::appliesTo() const {
    // Boundary cannot be included, because periodic boundaries
    //   would get the wrong value.
    return BlockDomain::bulk;
}

template<typename T>
void SetToNCoordinatesFunctional2D<T>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
}

template<typename T>
void SetToNCoordinatesFunctional2D<T>::rescale(double dxScale, double dtScale)
{ }


/* ************** Class MaskedSetToNCoordinatesFunctional2D ***************** */

template<typename T>
MaskedSetToNCoordinatesFunctional2D<T>::MaskedSetToNCoordinatesFunctional2D()
{ }

template<typename T>
void MaskedSetToNCoordinatesFunctional2D<T>::process (
        Box2D domain, NTensorField2D<T>& field,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( field.getNdim()==2 );
    Dot2D maskOfs = computeRelativeDisplacement(field, mask);
    Dot2D relativeOffset = field.getLocation();
    Array<plint,2> ofs(relativeOffset.x, relativeOffset.y);
    Array<plint,2> pos;
    for ( pos[0]=domain.x0; pos[0]<=domain.x1; ++pos[0] ) {
        for ( pos[1]=domain.y0; pos[1]<=domain.y1; ++pos[1] ) {
            if (*mask.get(pos[0]+maskOfs.x,pos[1]+maskOfs.y)) {
                T* cell = field.get(pos[0], pos[1]);
                cell[0] = (T) (pos[0]+ofs[0]);
                cell[1] = (T) (pos[1]+ofs[1]);
            }
        }
    }
}

template<typename T>
MaskedSetToNCoordinatesFunctional2D<T>* MaskedSetToNCoordinatesFunctional2D<T>::clone() const {
    return new MaskedSetToNCoordinatesFunctional2D<T>(*this);
}

template<typename T>
BlockDomain::DomainT MaskedSetToNCoordinatesFunctional2D<T>::appliesTo() const {
    // Boundary cannot be included, because periodic boundaries
    //   would get the wrong value.
    return BlockDomain::bulk;
}

template<typename T>
void MaskedSetToNCoordinatesFunctional2D<T>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
void MaskedSetToNCoordinatesFunctional2D<T>::rescale(double dxScale, double dtScale)
{ }


/* ************** Class SetTensorComponentFunctional2D ***************** */

template<typename T>
SetNTensorComponentFunctional2D<T>::SetNTensorComponentFunctional2D(int whichDim_)
    : whichDim(whichDim_)
{ }

template<typename T>
void SetNTensorComponentFunctional2D<T>::process (
        Box2D domain, NTensorField2D<T>& scalarField,
                      NTensorField2D<T>& tensorField )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    PLB_PRECONDITION( tensorField.getNdim()==2 );
    Dot2D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            tensorField.get(iX+offset.x, iY+offset.y)[whichDim] =
                *scalarField.get(iX,iY);
        }
    }
}

template<typename T>
SetNTensorComponentFunctional2D<T>* SetNTensorComponentFunctional2D<T>::clone() const {
    return new SetNTensorComponentFunctional2D<T>(*this);
}

template<typename T>
BlockDomain::DomainT SetNTensorComponentFunctional2D<T>::appliesTo() const {
    // Boundary cannot be included, because periodic boundaries
    //   would get the wrong value.
    return BlockDomain::bulk;
}

template<typename T>
void SetNTensorComponentFunctional2D<T>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
void SetNTensorComponentFunctional2D<T>::rescale(double dxScale, double dtScale)
{ }

/* ************** Class MaskedSetTensorComponentFunctional2D ***************** */

template<typename T>
MaskedSetNTensorComponentFunctional2D<T>::MaskedSetNTensorComponentFunctional2D(int whichDim_)
    : whichDim(whichDim_)
{ }

template<typename T>
void MaskedSetNTensorComponentFunctional2D<T>::process (
        Box2D domain, NTensorField2D<T>& scalarField,
                      NTensorField2D<T>& tensorField,
                      NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    PLB_PRECONDITION( tensorField.getNdim()==2 );
    Dot2D offset = computeRelativeDisplacement(scalarField, tensorField);
    Dot2D maskOfs = computeRelativeDisplacement(scalarField, mask);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            if (*mask.get(iX+maskOfs.x,iY+maskOfs.y)) {
                tensorField.get(iX+offset.x, iY+offset.y)[whichDim] =
                    *scalarField.get(iX,iY);
            }
        }
    }
}

template<typename T>
MaskedSetNTensorComponentFunctional2D<T>* MaskedSetNTensorComponentFunctional2D<T>::clone() const {
    return new MaskedSetNTensorComponentFunctional2D<T>(*this);
}

template<typename T>
BlockDomain::DomainT MaskedSetNTensorComponentFunctional2D<T>::appliesTo() const {
    // Boundary cannot be included, because periodic boundaries
    //   would get the wrong value.
    return BlockDomain::bulk;
}

template<typename T>
void MaskedSetNTensorComponentFunctional2D<T>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T>
void MaskedSetNTensorComponentFunctional2D<T>::rescale(double dxScale, double dtScale)
{ }


/* ************** Class AssignNTensorFunctional2D ***************** */

template<typename T>
AssignNTensorFunctional2D<T>::AssignNTensorFunctional2D()
{ }

template<typename T>
void AssignNTensorFunctional2D<T>::process (
        Box2D domain, NTensorField2D<T>& A,
                      NTensorField2D<T>& B )
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    plint ndim = A.getNdim();

    Dot2D offset = computeRelativeDisplacement(A,B);
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                *A.get(iX,iY) = *B.get(iX+offset.x,iY+offset.y);
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    A.get(iX,iY)[iDim] = B.get(iX+offset.x,iY+offset.y)[iDim];
                }
            }
        }
    }
}

template<typename T>
AssignNTensorFunctional2D<T>* AssignNTensorFunctional2D<T>::clone() const {
    return new AssignNTensorFunctional2D<T>(*this);
}

template<typename T>
BlockDomain::DomainT AssignNTensorFunctional2D<T>::appliesTo() const {
    // Boundary cannot be included, because periodic boundaries
    //   would get the wrong value.
    return BlockDomain::bulk;
}

template<typename T>
void AssignNTensorFunctional2D<T>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
void AssignNTensorFunctional2D<T>::rescale(double dxScale, double dtScale)
{ }


/* ************** Class MaskedAssignNTensorFunctional2D ***************** */

template<typename T>
MaskedAssignNTensorFunctional2D<T>::MaskedAssignNTensorFunctional2D()
{ }

template<typename T>
void MaskedAssignNTensorFunctional2D<T>::process (
        Box2D domain, NTensorField2D<T>& A,
                      NTensorField2D<T>& B,
                      NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    plint ndim = A.getNdim();

    Dot2D offset = computeRelativeDisplacement(A,B);
    Dot2D maskOfs = computeRelativeDisplacement(A, mask);
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x,iY+maskOfs.y)) {
                    *A.get(iX,iY) = *B.get(iX+offset.x,iY+offset.y);
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x,iY+maskOfs.y)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        A.get(iX,iY)[iDim] = B.get(iX+offset.x,iY+offset.y)[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
MaskedAssignNTensorFunctional2D<T>* MaskedAssignNTensorFunctional2D<T>::clone() const {
    return new MaskedAssignNTensorFunctional2D<T>(*this);
}

template<typename T>
BlockDomain::DomainT MaskedAssignNTensorFunctional2D<T>::appliesTo() const {
    // Boundary cannot be included, because periodic boundaries
    //   would get the wrong value.
    return BlockDomain::bulk;
}

template<typename T>
void MaskedAssignNTensorFunctional2D<T>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
}

template<typename T>
void MaskedAssignNTensorFunctional2D<T>::rescale(double dxScale, double dtScale)
{ }

}  // namespace plb

#endif  // SWIG_DATA_INITIALIZER_FUNCTIONAL_2D_HH
