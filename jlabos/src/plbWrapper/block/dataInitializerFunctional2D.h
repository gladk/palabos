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
#ifndef SWIG_DATA_INITIALIZER_FUNCTIONAL_2D_H
#define SWIG_DATA_INITIALIZER_FUNCTIONAL_2D_H

#include "core/globalDefs.h"
#include "atomicBlock/dataProcessingFunctional2D.h"
#include "core/dynamics.h"

namespace plb {

template<typename T>
class IniConstNTensorFunctional2D : public BoxProcessingFunctional2D_N<T>
{
public:
    IniConstNTensorFunctional2D(std::vector<T> const& value_);
    virtual void process(Box2D domain, NTensorField2D<T>& field);
    virtual IniConstNTensorFunctional2D<T>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual void rescale(double dxScale, double dtScale);
private:
    std::vector<T> value;
};

template<typename T>
class MaskedIniConstNTensorFunctional2D : public MaskedBoxProcessingFunctional2D_N<T>
{
public:
    MaskedIniConstNTensorFunctional2D(std::vector<T> const& value_);
    virtual void process(Box2D domain,
                         NTensorField2D<T>& field,
                         NTensorField2D<int>& mask );
    virtual MaskedIniConstNTensorFunctional2D<T>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual void rescale(double dxScale, double dtScale);
private:
    std::vector<T> value;
};

template<typename T>
class SetToNCoordinateFunctional2D : public BoxProcessingFunctional2D_N<T>
{
public:
    SetToNCoordinateFunctional2D(plint index_);
    virtual void process(Box2D domain, NTensorField2D<T>& field);
    virtual SetToNCoordinateFunctional2D<T>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual void rescale(double dxScale, double dtScale);
private:
    plint index;
};

template<typename T>
class MaskedSetToNCoordinateFunctional2D : public MaskedBoxProcessingFunctional2D_N<T>
{
public:
    MaskedSetToNCoordinateFunctional2D(plint index_);
    virtual void process(Box2D domain,
                         NTensorField2D<T>& field,
                         NTensorField2D<int>& mask );
    virtual MaskedSetToNCoordinateFunctional2D<T>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual void rescale(double dxScale, double dtScale);
private:
    plint index;
};

template<typename T>
class SetToNCoordinatesFunctional2D : public BoxProcessingFunctional2D_N<T>
{
public:
    SetToNCoordinatesFunctional2D();
    virtual void process(Box2D domain, NTensorField2D<T>& field);
    virtual SetToNCoordinatesFunctional2D<T>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual void rescale(double dxScale, double dtScale);
};

template<typename T>
class MaskedSetToNCoordinatesFunctional2D : public MaskedBoxProcessingFunctional2D_N<T>
{
public:
    MaskedSetToNCoordinatesFunctional2D();
    virtual void process(Box2D domain,
                         NTensorField2D<T>& field,
                         NTensorField2D<int>& mask);
    virtual MaskedSetToNCoordinatesFunctional2D<T>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual void rescale(double dxScale, double dtScale);
};

template<typename T>
class SetNTensorComponentFunctional2D : public BoxProcessingFunctional2D_NN<T,T>
{
public:
    SetNTensorComponentFunctional2D(int whichDim_);
    virtual void process(Box2D domain, NTensorField2D<T>& scalarField,
                                       NTensorField2D<T>& tensorField);
    virtual SetNTensorComponentFunctional2D<T>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual void rescale(double dxScale, double dtScale);
private:
    int whichDim;
};

template<typename T>
class MaskedSetNTensorComponentFunctional2D : public MaskedBoxProcessingFunctional2D_NN<T,T>
{
public:
    MaskedSetNTensorComponentFunctional2D(int whichDim_);
    virtual void process(Box2D domain, NTensorField2D<T>& scalarField,
                                       NTensorField2D<T>& tensorField,
                                       NTensorField2D<int>& mask);
    virtual MaskedSetNTensorComponentFunctional2D<T>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual void rescale(double dxScale, double dtScale);
private:
    int whichDim;
};

template<typename T>
class AssignNTensorFunctional2D : public BoxProcessingFunctional2D_NN<T,T>
{
public:
    AssignNTensorFunctional2D();
    virtual void process(Box2D domain, NTensorField2D<T>& A,
                                       NTensorField2D<T>& B);
    virtual AssignNTensorFunctional2D<T>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual void rescale(double dxScale, double dtScale);
};

template<typename T>
class MaskedAssignNTensorFunctional2D : public MaskedBoxProcessingFunctional2D_NN<T,T>
{
public:
    MaskedAssignNTensorFunctional2D();
    virtual void process(Box2D domain, NTensorField2D<T>& A,
                                       NTensorField2D<T>& B,
                                       NTensorField2D<int>& mask);
    virtual MaskedAssignNTensorFunctional2D<T>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual void rescale(double dxScale, double dtScale);
};

}  // namespace plb

#endif  // SWIG_DATA_INITIALIZER_FUNCTIONAL_2D_H
