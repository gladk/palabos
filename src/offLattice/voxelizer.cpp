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


#include "core/globalDefs.h"
#include "offLattice/voxelizer.h"
#include "offLattice/voxelizer.hh"
#include "atomicBlock/dataField3D.h"
#include "multiBlock/multiBlockGenerator3D.h"
#include "multiBlock/multiDataProcessorWrapper3D.h"
#include "multiBlock/multiDataProcessorWrapper3D.hh"
#include "multiBlock/reductiveMultiDataProcessorWrapper3D.h"
#include "dataProcessors/dataInitializerWrapper3D.h"
#include "dataProcessors/dataInitializerWrapper3D.hh"
#include "dataProcessors/dataAnalysisWrapper3D.h"
#include "dataProcessors/dataAnalysisWrapper3D.hh"
#include "dataProcessors/dataAnalysisFunctional3D.h"
#include "dataProcessors/dataAnalysisFunctional3D.hh"
#include "dataProcessors/metaStuffWrapper3D.h"

namespace plb {

void convertUndeterminedToFlag(MultiScalarField3D<int>& voxels, int flag) {
    applyProcessingFunctional(
            new UndeterminedToFlagFunctional3D(flag), voxels.getBoundingBox(), voxels);
}

UndeterminedToFlagFunctional3D::UndeterminedToFlagFunctional3D(int flag_)
    : flag(flag_)
{ }

UndeterminedToFlagFunctional3D* UndeterminedToFlagFunctional3D::clone() const {
    return new UndeterminedToFlagFunctional3D(*this);
}

void UndeterminedToFlagFunctional3D::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
}

void UndeterminedToFlagFunctional3D::process(Box3D domain, ScalarField3D<int>& voxels) {
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (voxels.get(iX,iY,iZ)==voxelFlag::undetermined) {
                    voxels.get(iX,iY,iZ) = flag;
                }
            }
        }
    }
}

} // namespace plb

