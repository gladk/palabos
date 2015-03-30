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

/** \file
 * Helper functions for domain initialization -- header file.
 */
#ifndef NTENSOR_ANALYSIS_FUNCTIONAL_3D_HH
#define NTENSOR_ANALYSIS_FUNCTIONAL_3D_HH

#include "dataProcessors/ntensorAnalysisFunctional3D.h"
#include "core/plbDebug.h"
#include "core/util.h"
#include "core/blockStatistics.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "finiteDifference/fdStencils1D.h"
#include "atomicBlock/atomicBlock3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataField3D.h"
#include <cmath>
#include <limits>

namespace plb {

template<typename T1, typename T2>
void CopyConvertNTensorFunctional3D<T1,T2>::process (
        Box3D domain, NTensorField3D<T1>& field1,
                      NTensorField3D<T2>& field2 )
{

    PLB_PRECONDITION( field1.getNdim() == field2.getNdim());
    plint ndim = field1.getNdim();
    Dot3D offset = computeRelativeDisplacement(field1, field2);

    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    for (int iDim=0; iDim<ndim; ++iDim) {
                        *field2.get(iX+offset.x,iY+offset.y,iZ+offset.z) =
                            (T2) *field1.get(iX,iY,iZ);
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
                    for (int iDim=0; iDim<ndim; ++iDim) {
                        field2.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim] =
                            (T2) field1.get(iX,iY,iZ)[iDim];
                    }
                }
            }
        }
    }
}

template<typename T1, typename T2>
CopyConvertNTensorFunctional3D<T1,T2>* CopyConvertNTensorFunctional3D<T1,T2>::clone() const
{
    return new CopyConvertNTensorFunctional3D<T1,T2>(*this);
}

template<typename T1, typename T2>
void CopyConvertNTensorFunctional3D<T1,T2>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T1, typename T2>
BlockDomain::DomainT CopyConvertNTensorFunctional3D<T1,T2>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

}  // namespace plb

#endif  // NTENSOR_ANALYSIS_FUNCTIONAL_3D_HH
