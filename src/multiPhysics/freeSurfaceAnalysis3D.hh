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

#ifndef FREE_SURFACE_ANALYSIS_3D_HH
#define FREE_SURFACE_ANALYSIS_3D_HH

#include "multiPhysics/freeSurfaceAnalysis3D.h"

namespace plb {

template<typename T, template<typename U> class Descriptor>
FS_AverageMassFunctional3D<T,Descriptor>::FS_AverageMassFunctional3D()
    : averageMassId(this->getStatistics().subscribeAverage())
{ }

template<typename T, template<typename U> class Descriptor>
void FS_AverageMassFunctional3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks )
{
    using namespace twoPhaseFlag;
    FreeSurfaceProcessorParam3D<T,Descriptor> param(atomicBlocks);
    BlockStatistics& statistics = this->getStatistics();

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                statistics.gatherAverage(averageMassId, param.mass(iX,iY,iZ));
                statistics.incrementStats();
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
FS_AverageMassFunctional3D<T,Descriptor>* FS_AverageMassFunctional3D<T,Descriptor>::clone() const
{
    return new FS_AverageMassFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
T FS_AverageMassFunctional3D<T,Descriptor>::getAverageMass() const {
    return this->getStatistics().getAverage(averageMassId);
}

template<typename T, template<typename U> class Descriptor>
T freeSurfaceAverageMass(std::vector<MultiBlock3D*> twoPhaseArgs, Box3D domain)
{
    FS_AverageMassFunctional3D<T,Descriptor> functional;
    applyProcessingFunctional(functional, domain, twoPhaseArgs);
    return functional.getAverageMass();
}


template<typename T, template<typename U> class Descriptor>
FS_TotalMassFunctional3D<T,Descriptor>::FS_TotalMassFunctional3D()
    : totalMassId(this->getStatistics().subscribeSum())
{ }

template<typename T, template<typename U> class Descriptor>
void FS_TotalMassFunctional3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks )
{
    using namespace twoPhaseFlag;
    FreeSurfaceProcessorParam3D<T,Descriptor> param(atomicBlocks);
    BlockStatistics& statistics = this->getStatistics();

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                statistics.gatherSum(totalMassId, param.mass(iX,iY,iZ));
                statistics.incrementStats();
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
FS_TotalMassFunctional3D<T,Descriptor>* FS_TotalMassFunctional3D<T,Descriptor>::clone() const
{
    return new FS_TotalMassFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
T FS_TotalMassFunctional3D<T,Descriptor>::getTotalMass() const {
    return this->getStatistics().getSum(totalMassId);
}


template<typename T, template<typename U> class Descriptor>
T freeSurfaceTotalMass(std::vector<MultiBlock3D*> twoPhaseArgs, Box3D domain)
{
    FS_TotalMassFunctional3D<T,Descriptor> functional;
    applyProcessingFunctional(functional, domain, twoPhaseArgs);
    return functional.getTotalMass();
}


template<typename T, template<typename U> class Descriptor>
FS_AverageDensityFunctional3D<T,Descriptor>::FS_AverageDensityFunctional3D()
    : averageDensityId(this->getStatistics().subscribeAverage())
{ }

template<typename T, template<typename U> class Descriptor>
void FS_AverageDensityFunctional3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks )
{
    using namespace twoPhaseFlag;
    FreeSurfaceProcessorParam3D<T,Descriptor> param(atomicBlocks);
    BlockStatistics& statistics = this->getStatistics();
    BounceBack<T,Descriptor> BBdynamics;

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (param.cell(iX,iY,iZ).getDynamics().getId() != BBdynamics.getId() ) {
                    statistics.gatherAverage(averageDensityId, param.getDensity(iX,iY,iZ));
                    statistics.incrementStats();
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
FS_AverageDensityFunctional3D<T,Descriptor>* FS_AverageDensityFunctional3D<T,Descriptor>::clone() const
{
    return new FS_AverageDensityFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
T FS_AverageDensityFunctional3D<T,Descriptor>::getAverageDensity() const {
    return this->getStatistics().getAverage(averageDensityId);
}


template<typename T, template<typename U> class Descriptor>
T freeSurfaceAverageDensity(std::vector<MultiBlock3D*> twoPhaseArgs, Box3D domain)
{
    FS_AverageDensityFunctional3D<T,Descriptor> functional;
    applyProcessingFunctional(functional, domain, twoPhaseArgs);
    return functional.getAverageDensity();
}


template<typename T, template<typename U> class Descriptor>
FS_AverageVolumeFractionFunctional3D<T,Descriptor>::FS_AverageVolumeFractionFunctional3D()
    : averageVfId(this->getStatistics().subscribeAverage())
{ }

template<typename T, template<typename U> class Descriptor>
void FS_AverageVolumeFractionFunctional3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks )
{
    using namespace twoPhaseFlag;
    FreeSurfaceProcessorParam3D<T,Descriptor> param(atomicBlocks);
    BlockStatistics& statistics = this->getStatistics();

    BounceBack<T,Descriptor> BBdynamics;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (param.cell(iX,iY,iZ).getDynamics().getId() != BBdynamics.getId() ) {
                            statistics.gatherAverage(averageVfId, param.volumeFraction(iX,iY,iZ));
                            statistics.incrementStats();
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
FS_AverageVolumeFractionFunctional3D<T,Descriptor>* FS_AverageVolumeFractionFunctional3D<T,Descriptor>::clone() const
{
    return new FS_AverageVolumeFractionFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
T FS_AverageVolumeFractionFunctional3D<T,Descriptor>::getAverageVolumeFraction() const {
    return this->getStatistics().getAverage(averageVfId);
}


template<typename T, template<typename U> class Descriptor>
T freeSurfaceAverageVolumeFraction(std::vector<MultiBlock3D*> twoPhaseArgs, Box3D domain) {
    FS_AverageVolumeFractionFunctional3D<T,Descriptor> functional;
    applyProcessingFunctional(functional, domain, twoPhaseArgs);
    return functional.getAverageVolumeFraction();
}



template<typename T, template<typename U> class Descriptor>
CountFreeSurfaceElementsFunctional3D<T,Descriptor>::CountFreeSurfaceElementsFunctional3D(plint flagToLookFor_)
    : numCellsId(this->getStatistics().subscribeIntSum()),
      flagToLookFor(flagToLookFor_)
{ }

template<typename T, template<typename U> class Descriptor>
void CountFreeSurfaceElementsFunctional3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks )
{
    using namespace twoPhaseFlag;
    FreeSurfaceProcessorParam3D<T,Descriptor> param(atomicBlocks);
    BlockStatistics& statistics = this->getStatistics();

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                int materialIndex = param.flag(iX,iY,iZ);
                if (materialIndex==flagToLookFor) {  // Fluid Cell
                    statistics.gatherIntSum(numCellsId, 1);
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
CountFreeSurfaceElementsFunctional3D<T,Descriptor>* CountFreeSurfaceElementsFunctional3D<T,Descriptor>::clone() const
{
    return new CountFreeSurfaceElementsFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
plint CountFreeSurfaceElementsFunctional3D<T,Descriptor>::getNumInterfaceCells() const {
    return this->getStatistics().getIntSum(numCellsId);
}

template<typename T, template<typename U> class Descriptor>
plint countFreeSurfaceElements(std::vector<MultiBlock3D*> twoPhaseArgs, plint flagToLookFor, Box3D domain) {
    CountFreeSurfaceElementsFunctional3D<T,Descriptor> functional(flagToLookFor);
    applyProcessingFunctional(functional, domain, twoPhaseArgs);
    return functional.getNumInterfaceCells();
}


template<typename T, template<typename U> class Descriptor>
FS_AverageMomentumFunctional3D<T,Descriptor>::FS_AverageMomentumFunctional3D()
    : averageMomentumId (
            this->getStatistics().subscribeAverage(),
            this->getStatistics().subscribeAverage(),
            this->getStatistics().subscribeAverage() )
{ }

template<typename T, template<typename U> class Descriptor>
void FS_AverageMomentumFunctional3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks )
{
    using namespace twoPhaseFlag;
    FreeSurfaceProcessorParam3D<T,Descriptor> param(atomicBlocks);
    BlockStatistics& statistics = this->getStatistics();
    BounceBack<T,Descriptor> BBdynamics;
    NoDynamics<T,Descriptor> NNdynamics;

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if ( param.cell(iX,iY,iZ).getDynamics().getId() != BBdynamics.getId() &&
                             param.cell(iX,iY,iZ).getDynamics().getId() != NNdynamics.getId()  )
                {
                    Array<T,Descriptor<T>::d> j = param.getMomentum(iX,iY,iZ);
                            statistics.gatherAverage(averageMomentumId[0], j[0]);
                            statistics.gatherAverage(averageMomentumId[1], j[1]);
                            statistics.gatherAverage(averageMomentumId[2], j[2]);
                            statistics.incrementStats();
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
FS_AverageMomentumFunctional3D<T,Descriptor>* FS_AverageMomentumFunctional3D<T,Descriptor>::clone() const
{
    return new FS_AverageMomentumFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
Array<T,3> FS_AverageMomentumFunctional3D<T,Descriptor>::getAverageMomentum() const {
    return Array<T,3> (
            this->getStatistics().getAverage(averageMomentumId[0]),
            this->getStatistics().getAverage(averageMomentumId[1]),
            this->getStatistics().getAverage(averageMomentumId[2]) );
}

template<typename T, template<typename U> class Descriptor>
Array<T,3> freeSurfaceAverageMomentum(std::vector<MultiBlock3D*> twoPhaseArgs, Box3D domain)
{
    FS_AverageMomentumFunctional3D<T,Descriptor> functional;
    applyProcessingFunctional(functional, domain, twoPhaseArgs);
    return functional.getAverageMomentum();
}


template<typename T, template<typename U> class Descriptor>
FS_AverageHeightFunctional3D<T,Descriptor>::FS_AverageHeightFunctional3D()
    : averageHeightId (this->getStatistics().subscribeAverage())
{ }

template<typename T, template<typename U> class Descriptor>
void FS_AverageHeightFunctional3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks )
{
    using namespace twoPhaseFlag;
    FreeSurfaceProcessorParam3D<T,Descriptor> param(atomicBlocks);
    BlockStatistics& statistics = this->getStatistics();
    BounceBack<T,Descriptor> BBdynamics;
    Dot3D absOffset = param.absOffset();
    
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (param.cell(iX,iY,iZ).getDynamics().getId() != BBdynamics.getId() )
                {
                    T localHeight = T(0);
                    if (param.volumeFraction(iX,iY,iZ) == 1) {
                        localHeight = T(absOffset.z);
                    }
                    statistics.gatherAverage(averageHeightId, localHeight);
                    statistics.incrementStats();
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
FS_AverageHeightFunctional3D<T,Descriptor>* FS_AverageHeightFunctional3D<T,Descriptor>::clone() const
{
    return new FS_AverageHeightFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
T FS_AverageHeightFunctional3D<T,Descriptor>::getAverageHeight() const {
    return this->getStatistics().getAverage(averageHeightId);
}

template<typename T, template<typename U> class Descriptor>
T freeSurfaceAverageHeight(std::vector<MultiBlock3D*> twoPhaseArgs, Box3D domain)
{
    FS_AverageHeightFunctional3D<T,Descriptor> functional;
    applyProcessingFunctional(functional, domain, twoPhaseArgs);
    return functional.getAverageHeight();
}


template<typename T, template<typename U> class Descriptor>
GetWaterLevelAtxyFunctional3D<T,Descriptor>::GetWaterLevelAtxyFunctional3D()
    : numFluidOccupiedCellId (this->getStatistics().subscribeIntSum())
{ }

template<typename T, template<typename U> class Descriptor>
void GetWaterLevelAtxyFunctional3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks )
{
    using namespace twoPhaseFlag;
    FreeSurfaceProcessorParam3D<T,Descriptor> param(atomicBlocks);
    BlockStatistics& statistics = this->getStatistics();
    plint bbDynamics = BounceBack<T,Descriptor>().getId();

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (param.cell(iX,iY,iZ).getDynamics().getId() != bbDynamics) {
                    if (param.volumeFraction(iX,iY,iZ) >= 0.5 ) {
                        statistics.gatherIntSum(numFluidOccupiedCellId, 1);
                    }
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
GetWaterLevelAtxyFunctional3D<T,Descriptor>* GetWaterLevelAtxyFunctional3D<T,Descriptor>::clone() const
{
    return new GetWaterLevelAtxyFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
plint GetWaterLevelAtxyFunctional3D<T,Descriptor>::getNumFluidCellsAtXY() const {
    return this->getStatistics().getIntSum(numFluidOccupiedCellId);
}

template<typename T, template<typename U> class Descriptor>
T getAverageHeightAtXY(std::vector<MultiBlock3D*> twoPhaseArgs, plint N, Box3D domain)
{
    GetWaterLevelAtxyFunctional3D<T,Descriptor> functional;
    applyProcessingFunctional(functional, domain, twoPhaseArgs);
    plint length_domain = domain.x1-domain.x0 ; // number of cell along y direction
    if (length_domain==0)
        length_domain =1;
    plint width_domain = domain.y1-domain.y0 ; // number of cell along y direction
    if (width_domain==0)
        width_domain =1;
    T heightAtXY = functional.getNumFluidCellsAtXY()/(T(N)*length_domain*width_domain);
    return heightAtXY;
}

}  // namespace plb

#endif  // FREE_SURFACE_ANALYSIS_3D_HH

