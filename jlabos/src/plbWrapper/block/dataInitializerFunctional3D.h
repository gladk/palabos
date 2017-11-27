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
 * Functionals for domain initialization -- header file.
 */
#ifndef SWIG_DATA_INITIALIZER_FUNCTIONAL_3D_H
#define SWIG_DATA_INITIALIZER_FUNCTIONAL_3D_H

#include "core/globalDefs.h"
#include "atomicBlock/dataProcessingFunctional3D.h"
#include "core/dynamics.h"

namespace plb {

template<typename T>
class IniConstNTensorFunctional3D : public BoxProcessingFunctional3D_N<T>
{
public:
    IniConstNTensorFunctional3D(std::vector<T> const& value_);
    virtual void process(Box3D domain, NTensorField3D<T>& field);
    virtual IniConstNTensorFunctional3D<T>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual void rescale(double dxScale, double dtScale);
private:
    std::vector<T> value;
};

template<typename T>
class MaskedIniConstNTensorFunctional3D : public MaskedBoxProcessingFunctional3D_N<T>
{
public:
    MaskedIniConstNTensorFunctional3D(std::vector<T> const& value_);
    virtual void process(Box3D domain,
                         NTensorField3D<T>& field,
                         NTensorField3D<int>& mask );
    virtual MaskedIniConstNTensorFunctional3D<T>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual void rescale(double dxScale, double dtScale);
private:
    std::vector<T> value;
};

template<typename T>
class SetToNCoordinateFunctional3D : public BoxProcessingFunctional3D_N<T>
{
public:
    SetToNCoordinateFunctional3D(plint index_);
    virtual void process(Box3D domain, NTensorField3D<T>& field);
    virtual SetToNCoordinateFunctional3D<T>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual void rescale(double dxScale, double dtScale);
private:
    plint index;
};

template<typename T>
class MaskedSetToNCoordinateFunctional3D : public MaskedBoxProcessingFunctional3D_N<T>
{
public:
    MaskedSetToNCoordinateFunctional3D(plint index_);
    virtual void process(Box3D domain,
                         NTensorField3D<T>& field,
                         NTensorField3D<int>& mask );
    virtual MaskedSetToNCoordinateFunctional3D<T>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual void rescale(double dxScale, double dtScale);
private:
    plint index;
};

template<typename T>
class SetToNCoordinatesFunctional3D : public BoxProcessingFunctional3D_N<T>
{
public:
    SetToNCoordinatesFunctional3D();
    virtual void process(Box3D domain, NTensorField3D<T>& field);
    virtual SetToNCoordinatesFunctional3D<T>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual void rescale(double dxScale, double dtScale);
};

template<typename T>
class MaskedSetToNCoordinatesFunctional3D : public MaskedBoxProcessingFunctional3D_N<T>
{
public:
    MaskedSetToNCoordinatesFunctional3D();
    virtual void process(Box3D domain,
                         NTensorField3D<T>& field,
                         NTensorField3D<int>& mask);
    virtual MaskedSetToNCoordinatesFunctional3D<T>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual void rescale(double dxScale, double dtScale);
};

template<typename T>
class SetNTensorComponentFunctional3D : public BoxProcessingFunctional3D_NN<T,T>
{
public:
    SetNTensorComponentFunctional3D(int whichDim_);
    virtual void process(Box3D domain, NTensorField3D<T>& scalarField,
                                       NTensorField3D<T>& tensorField);
    virtual SetNTensorComponentFunctional3D<T>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual void rescale(double dxScale, double dtScale);
private:
    int whichDim;
};

template<typename T>
class MaskedSetNTensorComponentFunctional3D : public MaskedBoxProcessingFunctional3D_NN<T,T>
{
public:
    MaskedSetNTensorComponentFunctional3D(int whichDim_);
    virtual void process(Box3D domain, NTensorField3D<T>& scalarField,
                                       NTensorField3D<T>& tensorField,
                                       NTensorField3D<int>& mask);
    virtual MaskedSetNTensorComponentFunctional3D<T>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual void rescale(double dxScale, double dtScale);
private:
    int whichDim;
};

template<typename T>
class AssignNTensorFunctional3D : public BoxProcessingFunctional3D_NN<T,T>
{
public:
    AssignNTensorFunctional3D();
    virtual void process(Box3D domain, NTensorField3D<T>& A,
                                       NTensorField3D<T>& B);
    virtual AssignNTensorFunctional3D<T>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual void rescale(double dxScale, double dtScale);
};

template<typename T>
class MaskedAssignNTensorFunctional3D : public MaskedBoxProcessingFunctional3D_NN<T,T>
{
public:
    MaskedAssignNTensorFunctional3D();
    virtual void process(Box3D domain, NTensorField3D<T>& A,
                                       NTensorField3D<T>& B,
                                       NTensorField3D<int>& mask);
    virtual MaskedAssignNTensorFunctional3D<T>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual void rescale(double dxScale, double dtScale);
};

}  // namespace plb

#endif  // SWIG_DATA_INITIALIZER_FUNCTIONAL_3D_H
