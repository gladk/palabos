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
#ifndef TYPE_CONVERTER_FUNCTIONAL_2D_H
#define TYPE_CONVERTER_FUNCTIONAL_2D_H

#include "core/globalDefs.h"
#include "atomicBlock/dataProcessingFunctional2D.h"

namespace plb {
	
/* *************** Data Functionals for scalar-fields **************** */

template<typename T, typename U>
class FromComplexToRealScalarFieldFunctional2D: public BoxProcessingFunctional2D_SS<T,U>
{
public:
	virtual void process(Box2D domain, ScalarField2D<T>& field1, ScalarField2D<U>& field2);
	virtual FromComplexToRealScalarFieldFunctional2D<T,U>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
};
template<typename T, typename U>
class FromComplexToImaginaryScalarFieldFunctional2D: public BoxProcessingFunctional2D_SS<T,U>
{
public:
	virtual void process(Box2D domain, ScalarField2D<T>& field1, ScalarField2D<U>& field2);
	virtual FromComplexToImaginaryScalarFieldFunctional2D<T,U>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
};
	
/* *************** Data Functionals for Tensor-fields **************** */

template<typename T, typename U,int d>
class FromComplexToRealTensorFieldFunctional2D: public BoxProcessingFunctional2D_TT<T,d,U,d>
{
public:
	virtual void process(Box2D domain, TensorField2D<T,d>& field1, TensorField2D<U,d>& field2);
	virtual FromComplexToRealTensorFieldFunctional2D<T,U,d>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
};
template<typename T, typename U,int d>
class FromComplexToImaginaryTensorFieldFunctional2D: public BoxProcessingFunctional2D_TT<T,d,U,d>
{
public:
	virtual void process(Box2D domain, TensorField2D<T,d>& field1, TensorField2D<U,d>& field2);
	virtual FromComplexToImaginaryTensorFieldFunctional2D<T,U,d>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
};

	
}  // namespace plb

#endif  // TYPE_CONVERTER_FUNCTIONAL_2D_H
