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

#include "plbWrapper/block/dataInitializerWrapper3D.h"
#include "plbWrapper/block/dataInitializerWrapper3D.hh"

namespace plb {

template void setToConstant(MultiNTensorField3D<PRECOMP_T>& field, Box3D domain, PRECOMP_T value);
template void maskedSetToConstant(MultiNTensorField3D<PRECOMP_T>& field, MultiNTensorField3D<int>& mask, Box3D domain, PRECOMP_T value);
template void setToConstant( MultiNTensorField3D<PRECOMP_T>& field, Box3D domain,
                             PRECOMP_T* in_value, int size );
template void maskedSetToConstant( MultiNTensorField3D<PRECOMP_T>& field, MultiNTensorField3D<int>& mask, Box3D domain,
                                   PRECOMP_T* in_value, int size );
template void setToCoordinate(MultiNTensorField3D<PRECOMP_T>& field, Box3D domain, plint index);
template void maskedSetToCoordinate(MultiNTensorField3D<PRECOMP_T>& field, MultiNTensorField3D<int>& mask, Box3D domain, plint index);
template void setToCoordinates(MultiNTensorField3D<PRECOMP_T>& field, Box3D domain);
template void maskedSetToCoordinates(MultiNTensorField3D<PRECOMP_T>& field, MultiNTensorField3D<int>& mask, Box3D domain);

template void assignComponent( MultiNTensorField3D<PRECOMP_T>& tensorField, int whichComponent,
                               MultiNTensorField3D<PRECOMP_T>& scalarField, Box3D domain );
template void maskedAssignComponent( MultiNTensorField3D<PRECOMP_T>& tensorField, int whichComponent,
                                     MultiNTensorField3D<PRECOMP_T>& scalarField, MultiNTensorField3D<int>& mask, Box3D domain );

}  // namespace plb
