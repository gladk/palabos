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
#ifndef TYPE_CONVERTER_FUNCTIONAL_2D_HH
#define TYPE_CONVERTER_FUNCTIONAL_2D_HH

#include "dataProcessors/dataAnalysisFunctional2D.h"
#include "core/plbDebug.h"
#include "atomicBlock/atomicBlock2D.h"
#include "atomicBlock/dataField2D.h"
#include "typeConverterFunctional2D.h"
#include <cmath>

namespace plb {

/* *************** Data Functionals for scalar-fields **************** */

template<typename T, typename U>
void FromComplexToRealScalarFieldFunctional2D<T,U>::process (Box2D domain, ScalarField2D<T>& field1,
											    ScalarField2D<U>& field2 )
{
    Dot2D offset = computeRelativeDisplacement(field1, field2);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            field2.get(iX+offset.x,iY+offset.y) = field1.get(iX,iY).real();
        }
    }
}

template<typename T, typename U>
	FromComplexToRealScalarFieldFunctional2D<T,U>* FromComplexToRealScalarFieldFunctional2D<T,U>::clone() const
{
    return new FromComplexToRealScalarFieldFunctional2D<T,U>(*this);
}

template<typename T, typename U>
void FromComplexToRealScalarFieldFunctional2D<T,U>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}


template<typename T, typename U>
void FromComplexToImaginaryScalarFieldFunctional2D<T,U>::process (Box2D domain, ScalarField2D<T>& field1,
												  ScalarField2D<U>& field2 )
{
	Dot2D offset = computeRelativeDisplacement(field1, field2);
	for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
		for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
			field2.get(iX+offset.x,iY+offset.y) = field1.get(iX,iY).imaginary();
		}
	}
}

template<typename T, typename U>
FromComplexToImaginaryScalarFieldFunctional2D<T,U>* FromComplexToImaginaryScalarFieldFunctional2D<T,U>::clone() const
{
	return new FromComplexToImaginaryScalarFieldFunctional2D<T,U>(*this);
}

template<typename T, typename U>
void FromComplexToImaginaryScalarFieldFunctional2D<T,U>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}


/* *************** Data Functionals for Tensor-fields **************** */

template<typename T, typename U,int d>
void FromComplexToRealTensorFieldFunctional2D<T,U,d>::process (Box2D domain, TensorField2D<T,d>& field1,
												  TensorField2D<U,d>& field2 )
{
	Dot2D offset = computeRelativeDisplacement(field1, field2);
	for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
		for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
			for (plint id=0; id<d; id++ ) {
				field2.get(iX+offset.x,iY+offset.y)[id] = field1.get(iX,iY)[id].real();
			}
		}
	}
}

template<typename T, typename U,int d>
FromComplexToRealTensorFieldFunctional2D<T,U,d>* FromComplexToRealTensorFieldFunctional2D<T,U,d>::clone() const
{
	return new FromComplexToRealTensorFieldFunctional2D<T,U,d>(*this);
}

template<typename T, typename U,int d>
void FromComplexToRealTensorFieldFunctional2D<T,U,d>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}


template<typename T, typename U,int d>
void FromComplexToImaginaryTensorFieldFunctional2D<T,U,d>::process (Box2D domain, TensorField2D<T,d>& field1,
													   TensorField2D<U,d>& field2 )
{
	Dot2D offset = computeRelativeDisplacement(field1, field2);
	for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
		for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
			for (plint id=0; id<d; id++ ) {
				field2.get(iX+offset.x,iY+offset.y)[id] = field1.get(iX,iY)[id].imaginary();
			}
		}
	}
}

template<typename T, typename U,int d>
FromComplexToImaginaryTensorFieldFunctional2D<T,U,d>* FromComplexToImaginaryTensorFieldFunctional2D<T,U,d>::clone() const
{
	return new FromComplexToImaginaryTensorFieldFunctional2D<T,U,d>(*this);
}

template<typename T, typename U,int d>
void FromComplexToImaginaryTensorFieldFunctional2D<T,U,d>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

}//namespace plb

#endif  // TYPE_CONVERTER_FUNCTIONAL_2D_HH
