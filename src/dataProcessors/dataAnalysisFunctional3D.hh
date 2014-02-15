/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2013 FlowKit Sarl
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
 * Data processors for data analysis -- header file.
 */
#ifndef DATA_ANALYSIS_FUNCTIONAL_3D_HH
#define DATA_ANALYSIS_FUNCTIONAL_3D_HH

#include "dataProcessors/dataAnalysisFunctional3D.h"
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

/* *************** PART I ******************************************** */
/* *************** Analysis of the block-lattice ********************* */
/* ******************************************************************* */

/* *************** Reductive Data Functionals for BlockLattice ******* */

template<typename T, template<typename U> class Descriptor> 
BoxSumRhoBarFunctional3D<T,Descriptor>::BoxSumRhoBarFunctional3D()
    : sumRhoBarId(this->getStatistics().subscribeSum())
{ }

template<typename T, template<typename U> class Descriptor> 
void BoxSumRhoBarFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Cell<T,Descriptor> const& cell = lattice.get(iX,iY,iZ);
                statistics.gatherSum (
                        sumRhoBarId, cell.getDynamics().computeRhoBar(cell)
                );
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxSumRhoBarFunctional3D<T,Descriptor>*
    BoxSumRhoBarFunctional3D<T,Descriptor>::clone() const
{
    return new BoxSumRhoBarFunctional3D(*this);
}

template<typename T, template<typename U> class Descriptor> 
T BoxSumRhoBarFunctional3D<T,Descriptor>::getSumRhoBar() const {
    return this->getStatistics().getSum(sumRhoBarId);
}


template<typename T, template<typename U> class Descriptor> 
BoxSumEnergyFunctional3D<T,Descriptor>::BoxSumEnergyFunctional3D()
    : sumEnergyId(this->getStatistics().subscribeSum())
{ }

template<typename T, template<typename U> class Descriptor> 
void BoxSumEnergyFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Array<T,Descriptor<T>::d> velocity;
                lattice.get(iX,iY,iZ).computeVelocity(velocity);
                T uNormSqr = VectorTemplate<T,Descriptor>::normSqr(velocity);
                statistics.gatherSum(sumEnergyId, uNormSqr);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxSumEnergyFunctional3D<T,Descriptor>*
    BoxSumEnergyFunctional3D<T,Descriptor>::clone() const
{
    return new BoxSumEnergyFunctional3D(*this);
}

template<typename T, template<typename U> class Descriptor> 
T BoxSumEnergyFunctional3D<T,Descriptor>::getSumEnergy() const {
    return this->getStatistics().getSum(sumEnergyId) / (T)2;
}

template<typename T, template<typename U> class Descriptor> 
void BoxSumEnergyFunctional3D<T,Descriptor>::getDimensionsX(std::vector<int>& dimensions) const
{
    dimensions.resize(1);
    dimensions[0] = 2;
}

template<typename T, template<typename U> class Descriptor> 
void BoxSumEnergyFunctional3D<T,Descriptor>::getDimensionsT(std::vector<int>& dimensions) const
{
    dimensions.resize(1);
    dimensions[0] = -2;
}


/* ******** DensitySingleProbe3D *********************************** */

template<typename T, template<typename U> class Descriptor>
DensitySingleProbe3D<T,Descriptor>::DensitySingleProbe3D (
        std::vector<Array<T,3> > const& positions_ )
    : positions(positions_)
{
    densityIds.resize(positions.size());
    for (pluint i=0; i<positions.size(); ++i) {
        densityIds[i] = this->getStatistics().subscribeSum();
    }
}

template<typename T, template<typename U> class Descriptor>
void DensitySingleProbe3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    std::vector<Dot3D> cellPos(8);
    std::vector<T> weights(8);
    T density;
    for (pluint i=0; i<positions.size(); ++i) {
        density= (T) 0;
        Array<T,3> position(positions[i]);
        Dot3D referenceCellPos((plint)position[0], (plint)position[1], (plint)position[2]);
        referenceCellPos -= lattice.getLocation();
        if (contained(referenceCellPos, domain)) {
            linearInterpolationCoefficients(lattice, position, cellPos, weights);
            for (plint iCell=0; iCell<8; ++iCell) {
                Cell<T,Descriptor> const& cell = lattice.get(cellPos[iCell].x,cellPos[iCell].y,cellPos[iCell].z);
                T cellDensity = cell.computeDensity();
                density+=weights[iCell]*cellDensity;
            }
        }
       this->getStatistics().gatherSum(densityIds[i], density);
    }
}

template<typename T, template<typename U> class Descriptor>
DensitySingleProbe3D<T,Descriptor>*
    DensitySingleProbe3D<T,Descriptor>::clone() const
{
    return new DensitySingleProbe3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
std::vector<T> DensitySingleProbe3D<T,Descriptor>::getDensities() const {
    std::vector<T> densities(positions.size());
    for (pluint i=0; i<positions.size(); ++i) {
        densities[i] = this->getStatistics().getSum(densityIds[i]);
    }
    return densities;
}


/* ******** VelocitySingleProbe3D *********************************** */

template<typename T, template<typename U> class Descriptor>
VelocitySingleProbe3D<T,Descriptor>::VelocitySingleProbe3D (
        std::vector<Array<T,3> > const& positions_ )
    : positions(positions_)
{
    velIds.resize(positions.size());
    for (pluint iVel=0; iVel<positions.size(); ++iVel) {
        velIds[iVel][0] = this->getStatistics().subscribeSum();
        velIds[iVel][1] = this->getStatistics().subscribeSum();
        velIds[iVel][2] = this->getStatistics().subscribeSum();
    }
}

template<typename T, template<typename U> class Descriptor>
void VelocitySingleProbe3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    std::vector<Dot3D> cellPos(8);
    std::vector<T> weights(8);
    Array<T,3> velocity;
    for (pluint iVel=0; iVel<positions.size(); ++iVel) {
        velocity.resetToZero();
        Array<T,3> position(positions[iVel]);
        Dot3D referenceCellPos((plint)position[0], (plint)position[1], (plint)position[2]);
        referenceCellPos -= lattice.getLocation();
        if (contained(referenceCellPos, domain)) {
            linearInterpolationCoefficients(lattice, position, cellPos, weights);
            for (plint iCell=0; iCell<8; ++iCell) {
                Cell<T,Descriptor> const& cell = lattice.get(cellPos[iCell].x,cellPos[iCell].y,cellPos[iCell].z);
                Array<T,3> cellVelocity;
                cell.computeVelocity(cellVelocity);
                velocity+=weights[iCell]*cellVelocity;
            }
        }
       this->getStatistics().gatherSum(velIds[iVel][0], velocity[0]);
       this->getStatistics().gatherSum(velIds[iVel][1], velocity[1]);
       this->getStatistics().gatherSum(velIds[iVel][2], velocity[2]);
    }
}

template<typename T, template<typename U> class Descriptor>
VelocitySingleProbe3D<T,Descriptor>*
    VelocitySingleProbe3D<T,Descriptor>::clone() const
{
    return new VelocitySingleProbe3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
std::vector<Array<T,3> > VelocitySingleProbe3D<T,Descriptor>::getVelocities() const {
    std::vector<Array<T,3> > velocities(positions.size());
    for (pluint iVel=0; iVel<positions.size(); ++iVel) {
        velocities[iVel][0] = this->getStatistics().getSum(velIds[iVel][0]);
        velocities[iVel][1] = this->getStatistics().getSum(velIds[iVel][1]);
        velocities[iVel][2] = this->getStatistics().getSum(velIds[iVel][2]);
    }
    return velocities;
}


/* ******** VorticitySingleProbe3D *********************************** */

template<typename T, template<typename U> class Descriptor>
VorticitySingleProbe3D<T,Descriptor>::VorticitySingleProbe3D (
        std::vector<Array<T,3> > const& positions_ )
    : positions(positions_)
{
    vorticityIds.resize(positions.size());
    for (pluint iVel=0; iVel<positions.size(); ++iVel) {
        vorticityIds[iVel][0] = this->getStatistics().subscribeSum();
        vorticityIds[iVel][1] = this->getStatistics().subscribeSum();
        vorticityIds[iVel][2] = this->getStatistics().subscribeSum();
    }
}

template<typename T, template<typename U> class Descriptor>
void VorticitySingleProbe3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    Array<T,3> vorticity;
    for (pluint iVort=0; iVort<positions.size(); ++iVort) {
        vorticity.resetToZero();
        Array<T,3> position(positions[iVort]);
        Dot3D referenceCellPos((plint)position[0], (plint)position[1], (plint)position[2]);
        referenceCellPos -= lattice.getLocation();
        if (contained(referenceCellPos, domain)) {
            vorticity = fdLattice::firstOrderBulkVorticity(lattice, referenceCellPos.x, referenceCellPos.y, referenceCellPos.z);
        }
       this->getStatistics().gatherSum(vorticityIds[iVort][0], vorticity[0]);
       this->getStatistics().gatherSum(vorticityIds[iVort][1], vorticity[1]);
       this->getStatistics().gatherSum(vorticityIds[iVort][2], vorticity[2]);
    }
}

template<typename T, template<typename U> class Descriptor>
VorticitySingleProbe3D<T,Descriptor>*
    VorticitySingleProbe3D<T,Descriptor>::clone() const
{
    return new VorticitySingleProbe3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
std::vector<Array<T,3> > VorticitySingleProbe3D<T,Descriptor>::getVorticities() const {
    std::vector<Array<T,3> > vorticities(positions.size());
    for (pluint iVort=0; iVort<positions.size(); ++iVort) {
        vorticities[iVort][0] = this->getStatistics().getSum(vorticityIds[iVort][0]);
        vorticities[iVort][1] = this->getStatistics().getSum(vorticityIds[iVort][1]);
        vorticities[iVort][2] = this->getStatistics().getSum(vorticityIds[iVort][2]);
    }
    return vorticities;
}


/* *************** Data Functionals for BlockLattice ***************** */

template<typename T, template<typename U> class Descriptor> 
void CopyPopulationsFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& latticeFrom,
                      BlockLattice3D<T,Descriptor>& latticeTo )
{
    Dot3D offset = computeRelativeDisplacement(latticeFrom, latticeTo);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        plint ofX = iX + offset.x;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            plint ofY = iY + offset.y;
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                latticeTo.get(ofX,ofY,iZ+offset.z).
                    attributeValues(latticeFrom.get(iX,iY,iZ));
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
CopyPopulationsFunctional3D<T,Descriptor>* CopyPopulationsFunctional3D<T,Descriptor>::clone() const
{
    return new CopyPopulationsFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void CopyPopulationsFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT CopyPopulationsFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

template<typename T1, template<typename U1> class Descriptor1, typename T2, template<typename U2> class Descriptor2> 
void CopyConvertPopulationsFunctional3D<T1,Descriptor1, T2,Descriptor2>::process (
    Box3D domain, BlockLattice3D<T1,Descriptor1>& latticeFrom,
    BlockLattice3D<T2,Descriptor2>& latticeTo )
{
    Dot3D offset = computeRelativeDisplacement(latticeFrom, latticeTo);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        plint ofX = iX + offset.x;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            plint ofY = iY + offset.y;
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    latticeTo.get(ofX,ofY,iZ+offset.z).
                    attributeValues(latticeFrom.get(iX,iY,iZ));
            }
        }
    }
}

template<typename T1, template<typename U1> class Descriptor1, typename T2, template<typename U2> class Descriptor2> 
CopyConvertPopulationsFunctional3D<T1,Descriptor1,T2,Descriptor2>* CopyConvertPopulationsFunctional3D<T1,Descriptor1,T2,Descriptor2>::clone() const
{
    return new CopyConvertPopulationsFunctional3D<T1,Descriptor1,T2,Descriptor2>(*this);
}

template<typename T1, template<typename U1> class Descriptor1, typename T2, template<typename U2> class Descriptor2> 
void CopyConvertPopulationsFunctional3D<T1,Descriptor1,T2,Descriptor2>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T1, template<typename U1> class Descriptor1, typename T2, template<typename U2> class Descriptor2> 
BlockDomain::DomainT CopyConvertPopulationsFunctional3D<T1,Descriptor1,T2,Descriptor2>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, template<typename U> class Descriptor> 
void LatticeCopyAllFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& latticeFrom,
                      BlockLattice3D<T,Descriptor>& latticeTo )
{
    std::vector<char> data;
    Dot3D offset = computeRelativeDisplacement(latticeFrom, latticeTo);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Cell<T,Descriptor> const& fromCell = latticeFrom.get(iX,iY,iZ);
                Cell<T,Descriptor>& toCell = latticeTo.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                toCell.attributeValues(fromCell);

                data.clear();
                HierarchicSerializer serializer(data, fromCell.getDynamics().getId());
                fromCell.getDynamics().serialize(serializer);

                HierarchicUnserializer unSerializer(data, 0);
                toCell.getDynamics().unserialize(unSerializer);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
LatticeCopyAllFunctional3D<T,Descriptor>* LatticeCopyAllFunctional3D<T,Descriptor>::clone() const
{
    return new LatticeCopyAllFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void LatticeCopyAllFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT LatticeCopyAllFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}



template<typename T, template<typename U> class Descriptor> 
void LatticeRegenerateFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& latticeFrom,
                      BlockLattice3D<T,Descriptor>& latticeTo )
{
    Dot3D offset = computeRelativeDisplacement(latticeFrom, latticeTo);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                latticeTo.get(iX+offset.x,iY+offset.y,iZ+offset.z).attributeValues(latticeFrom.get(iX,iY,iZ));
                latticeTo.attributeDynamics(iX+offset.x,iY+offset.y,iZ+offset.z,
                                           latticeFrom.get(iX,iY,iZ).getDynamics().clone());
                latticeTo.get(iX+offset.x,iY+offset.y,iZ+offset.z).
                    attributeValues(latticeFrom.get(iX,iY,iZ));
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
LatticeRegenerateFunctional3D<T,Descriptor>* LatticeRegenerateFunctional3D<T,Descriptor>::clone() const
{
    return new LatticeRegenerateFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void LatticeRegenerateFunctional3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::nothing;
    // Full dynamics object must be recreated, because this data processor
    //   re-attributes a new dynamics and acts on the bulk only.
    modified[1] = modif::dataStructure;
}

template<typename T, template<typename U> class Descriptor> 
void BoxDensityFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& scalarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                scalarField.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                    = lattice.get(iX,iY,iZ).computeDensity();
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxDensityFunctional3D<T,Descriptor>* BoxDensityFunctional3D<T,Descriptor>::clone() const
{
    return new BoxDensityFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void BoxDensityFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT BoxDensityFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, template<typename U> class Descriptor> 
void BoxRhoBarFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& scalarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
       for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Cell<T,Descriptor> const& cell = lattice.get(iX,iY,iZ);
                scalarField.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                    = cell.getDynamics().computeRhoBar(cell);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxRhoBarFunctional3D<T,Descriptor>* BoxRhoBarFunctional3D<T,Descriptor>::clone() const
{
    return new BoxRhoBarFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void BoxRhoBarFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT BoxRhoBarFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, template<typename U> class Descriptor> 
void BoxRhoBarJfunctional3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> fields )
{
    BlockLattice3D<T,Descriptor>& lattice = *dynamic_cast<BlockLattice3D<T,Descriptor>*>(fields[0]);
    ScalarField3D<T>& rhoBarField = *dynamic_cast<ScalarField3D<T>*>(fields[1]);
    TensorField3D<T,3>& jField = *dynamic_cast<TensorField3D<T,3>*>(fields[2]);
    Dot3D offset1 = computeRelativeDisplacement(lattice, rhoBarField);
    Dot3D offset2 = computeRelativeDisplacement(lattice, jField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
       for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Cell<T,Descriptor> const& cell = lattice.get(iX,iY,iZ);
                momentTemplates<T,Descriptor>::get_rhoBar_j (
                        cell,
                        rhoBarField.get(iX+offset1.x,iY+offset1.y,iZ+offset1.z),
                        jField.get(iX+offset2.x,iY+offset2.y,iZ+offset2.z) );
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxRhoBarJfunctional3D<T,Descriptor>* BoxRhoBarJfunctional3D<T,Descriptor>::clone() const
{
    return new BoxRhoBarJfunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void BoxRhoBarJfunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;  // lattice
    modified[1] = modif::staticVariables;   // rhoBar
    modified[2] = modif::staticVariables;   // j
}

template<typename T, template<typename U> class Descriptor> 
void PackedRhoBarJfunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                      NTensorField3D<T>& rhoBarJField )
{
    Dot3D offset = computeRelativeDisplacement(lattice, rhoBarJField);
    T rhoBar;
    Array<T,3> j;
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
       for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Cell<T,Descriptor> const& cell = lattice.get(iX,iY,iZ);
                momentTemplates<T,Descriptor>::get_rhoBar_j(cell, rhoBar, j);
                T* rhoBarJ = rhoBarJField.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                *rhoBarJ = rhoBar;
                j.to_cArray(rhoBarJ+1);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
PackedRhoBarJfunctional3D<T,Descriptor>* PackedRhoBarJfunctional3D<T,Descriptor>::clone() const
{
    return new PackedRhoBarJfunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void PackedRhoBarJfunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;  // lattice
    modified[1] = modif::staticVariables;   // rhoBarJ
}


template<typename T>
void DensityFromRhoBarJfunctional3D<T>::process (
        Box3D domain, ScalarField3D<T>& density,
                      NTensorField3D<T>& rhoBarJField )
{
    Dot3D offset = computeRelativeDisplacement(density, rhoBarJField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
       for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                T rhoBar = *rhoBarJField.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                density.get(iX,iY,iZ) = rhoBar+(T)1;
            }
        }
    }
}

template<typename T>
DensityFromRhoBarJfunctional3D<T>* DensityFromRhoBarJfunctional3D<T>::clone() const
{
    return new DensityFromRhoBarJfunctional3D<T>(*this);
}

template<typename T>
void DensityFromRhoBarJfunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;  // density
    modified[1] = modif::nothing;  // rhoBarJ
}

template<typename T>
VelocityFromRhoBarJfunctional3D<T>::VelocityFromRhoBarJfunctional3D(bool velIsJ_)
    : velIsJ(velIsJ_)
{ } 

template<typename T>
void VelocityFromRhoBarJfunctional3D<T>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> fields)
{
    PLB_ASSERT(fields.size()==2);
    TensorField3D<T,3>* velocity = dynamic_cast<TensorField3D<T,3>*>(fields[0]);
    NTensorField3D<T>* rhoBarJfield = dynamic_cast<NTensorField3D<T>*>(fields[1]);
    PLB_ASSERT(velocity);
    PLB_ASSERT(rhoBarJfield);
    Dot3D offset = computeRelativeDisplacement(*velocity, *rhoBarJfield);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
       for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Array<T,3> nextVelocity;
                nextVelocity.from_cArray((*rhoBarJfield).get(iX+offset.x,iY+offset.y,iZ+offset.z)+1);
                if (!velIsJ) {
                    nextVelocity /= (*(*rhoBarJfield).get(iX+offset.x,iY+offset.y,iZ+offset.z))+(T)1;
                }
                (*velocity).get(iX,iY,iZ) = nextVelocity;
            }
        }
    }
}

template<typename T>
VelocityFromRhoBarJfunctional3D<T>* VelocityFromRhoBarJfunctional3D<T>::clone() const
{
    return new VelocityFromRhoBarJfunctional3D<T>(*this);
}

template<typename T>
void VelocityFromRhoBarJfunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;  // velocity
    modified[1] = modif::nothing;  // rhoBarJ
}


template<typename T, template<typename U> class Descriptor> 
void BoxKineticEnergyFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& scalarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Array<T,Descriptor<T>::d> velocity;
                lattice.get(iX,iY,iZ).computeVelocity(velocity);
                scalarField.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                    = VectorTemplate<T,Descriptor>::normSqr(velocity) / (T)2;
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxKineticEnergyFunctional3D<T,Descriptor>* BoxKineticEnergyFunctional3D<T,Descriptor>::clone() const
{
    return new BoxKineticEnergyFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void BoxKineticEnergyFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT BoxKineticEnergyFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}



template<typename T, template<typename U> class Descriptor> 
void BoxVelocityNormFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& scalarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Array<T,Descriptor<T>::d> velocity;
                lattice.get(iX,iY,iZ).computeVelocity(velocity);
                scalarField.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                    = sqrt( VectorTemplate<T,Descriptor>::normSqr(velocity) );
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxVelocityNormFunctional3D<T,Descriptor>* BoxVelocityNormFunctional3D<T,Descriptor>::clone() const
{
    return new BoxVelocityNormFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void BoxVelocityNormFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT BoxVelocityNormFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor> 
BoxVelocityComponentFunctional3D<T,Descriptor>::BoxVelocityComponentFunctional3D(int iComponent_)
    : iComponent(iComponent_)
{ }

template<typename T, template<typename U> class Descriptor> 
void BoxVelocityComponentFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& scalarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Array<T,Descriptor<T>::d> velocity;
                lattice.get(iX,iY,iZ).computeVelocity(velocity);
                scalarField.get(iX+offset.x,iY+offset.y,iZ+offset.z) = velocity[iComponent];
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxVelocityComponentFunctional3D<T,Descriptor>* BoxVelocityComponentFunctional3D<T,Descriptor>::clone() const
{
    return new BoxVelocityComponentFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void BoxVelocityComponentFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT BoxVelocityComponentFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, template<typename U> class Descriptor> 
void BoxVelocityFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice, TensorField3D<T,Descriptor<T>::d>& tensorField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                lattice.get(iX,iY,iZ).computeVelocity (
                        tensorField.get(iX+offset.x,iY+offset.y,iZ+offset.z) );
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxVelocityFunctional3D<T,Descriptor>* BoxVelocityFunctional3D<T,Descriptor>::clone() const
{
    return new BoxVelocityFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void BoxVelocityFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}


template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT BoxVelocityFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, template<typename U> class Descriptor> 
void BoxTemperatureFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& scalarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                scalarField.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                    = lattice.get(iX,iY,iZ).computeTemperature();
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxTemperatureFunctional3D<T,Descriptor>* BoxTemperatureFunctional3D<T,Descriptor>::clone() const
{
    return new BoxTemperatureFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void BoxTemperatureFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT BoxTemperatureFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, template<typename U> class Descriptor> 
void BoxPiNeqFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
        TensorField3D<T,SymmetricTensor<T,Descriptor>::n>& PiNeq )
{
    Dot3D offset = computeRelativeDisplacement(lattice, PiNeq);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                lattice.get(iX,iY,iZ).computePiNeq (
                        PiNeq.get(iX+offset.x,iY+offset.y,iZ+offset.z) );
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxPiNeqFunctional3D<T,Descriptor>* BoxPiNeqFunctional3D<T,Descriptor>::clone() const
{
    return new BoxPiNeqFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void BoxPiNeqFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT BoxPiNeqFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, template<typename U> class Descriptor> 
void BoxShearStressFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
        TensorField3D<T,SymmetricTensor<T,Descriptor>::n>& PiNeq )
{
    Dot3D offset = computeRelativeDisplacement(lattice, PiNeq);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                lattice.get(iX,iY,iZ).computeShearStress (
                        PiNeq.get(iX+offset.x,iY+offset.y,iZ+offset.z) );
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxShearStressFunctional3D<T,Descriptor>* BoxShearStressFunctional3D<T,Descriptor>::clone() const
{
    return new BoxShearStressFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void BoxShearStressFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT BoxShearStressFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, template<typename U> class Descriptor> 
void BoxStrainRateFromStressFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
        TensorField3D<T,SymmetricTensor<T,Descriptor>::n>& S )
{
    Dot3D offset = computeRelativeDisplacement(lattice, S);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Cell<T,Descriptor> const& cell = lattice.get(iX,iY,iZ);
                Array<T,SymmetricTensor<T,Descriptor>::n>& element
                    = S.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                cell.computePiNeq(element);
                T omega     = cell.getDynamics().getOmega();
                T rhoBar    = cell.getDynamics().computeRhoBar(cell);
                T prefactor = - omega * Descriptor<T>::invCs2 *
                                Descriptor<T>::invRho(rhoBar) / (T)2;
                for (int iTensor=0; iTensor<SymmetricTensor<T,Descriptor>::n; ++iTensor) {
                    element[iTensor] *= prefactor;
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxStrainRateFromStressFunctional3D<T,Descriptor>* BoxStrainRateFromStressFunctional3D<T,Descriptor>::clone() const
{
    return new BoxStrainRateFromStressFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void BoxStrainRateFromStressFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT BoxStrainRateFromStressFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

template<typename T>
void BoxQcriterionFunctional3D<T>::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields) {
    TensorField3D<T,3>& vorticity = *dynamic_cast<TensorField3D<T,3>*>(fields[0]);
    TensorField3D<T,6>& strain = *dynamic_cast<TensorField3D<T,6>*>(fields[1]);
    ScalarField3D<T>& qCriterion = *dynamic_cast<ScalarField3D<T>*>(fields[2]);
    Dot3D offset1 = computeRelativeDisplacement(vorticity, strain);
    Dot3D offset2 = computeRelativeDisplacement(vorticity, qCriterion);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        plint oX1 = iX + offset1.x;
        plint oX2 = iX + offset2.x;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            plint oY1 = iY + offset1.y;
            plint oY2 = iY + offset2.y;
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint oZ1 = iZ + offset1.z;
                plint oZ2 = iZ + offset2.z;
                
                T vortNorm = VectorTemplateImpl<T,3>::normSqr(vorticity.get(iX,iY,iZ));
                T normStrain = SymmetricTensorImpl<T,3>::tensorNormSqr(strain.get(oX1,oY1,oZ1));
                
                qCriterion.get(oX2,oY2,oZ2) = (vortNorm-(T)2*normStrain)/(T)4;
            }
        }
    }
}

template<typename T>
BoxQcriterionFunctional3D<T>* BoxQcriterionFunctional3D<T>::clone() const {
    return new BoxQcriterionFunctional3D<T>(*this);
}

template<typename T>
void BoxQcriterionFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT BoxQcriterionFunctional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, template<typename U> class Descriptor> 
BoxPopulationFunctional3D<T,Descriptor>::BoxPopulationFunctional3D(plint iComponent_)
    : iComponent(iComponent_)
{ }

template<typename T, template<typename U> class Descriptor> 
void BoxPopulationFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& scalarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                scalarField.get(iX+offset.x,iY+offset.y,iZ+offset.z) = lattice.get(iX,iY,iZ)[iComponent];
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxPopulationFunctional3D<T,Descriptor>* BoxPopulationFunctional3D<T,Descriptor>::clone() const
{
    return new BoxPopulationFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void BoxPopulationFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT BoxPopulationFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, template<typename U> class Descriptor> 
BoxEquilibriumFunctional3D<T,Descriptor>::BoxEquilibriumFunctional3D(plint iComponent_)
    : iComponent(iComponent_)
{ }

template<typename T, template<typename U> class Descriptor> 
void BoxEquilibriumFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& equilibrium)
{
    Dot3D offset = computeRelativeDisplacement(lattice, equilibrium);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                T rhoBar;
                Array<T,Descriptor<T>::d> j;
                Cell<T,Descriptor> const& cell = lattice.get(iX,iY,iZ);
                cell.getDynamics().computeRhoBarJ(cell, rhoBar, j);
                T jSqr = normSqr(j);
                equilibrium.get(iX+offset.x,iY+offset.y,iZ+offset.z) =
                    cell.computeEquilibrium(iComponent, rhoBar, j, jSqr);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxEquilibriumFunctional3D<T,Descriptor>* BoxEquilibriumFunctional3D<T,Descriptor>::clone() const
{
    return new BoxEquilibriumFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void BoxEquilibriumFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT BoxEquilibriumFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, template<typename U> class Descriptor> 
BoxAllPopulationsFunctional3D<T,Descriptor>::BoxAllPopulationsFunctional3D()
{ }

template<typename T, template<typename U> class Descriptor> 
void BoxAllPopulationsFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice, TensorField3D<T,Descriptor<T>::q>& tensorField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        plint oX = iX+offset.x;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            plint oY = iY+offset.y;
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                    tensorField.get(oX,oY,iZ+offset.z)[iPop] = lattice.get(iX,iY,iZ)[iPop];
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxAllPopulationsFunctional3D<T,Descriptor>* BoxAllPopulationsFunctional3D<T,Descriptor>::clone() const
{
    return new BoxAllPopulationsFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void BoxAllPopulationsFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT BoxAllPopulationsFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor> 
void BoxAllEquilibriumFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice, TensorField3D<T,Descriptor<T>::q>& equilibrium)
{
    Dot3D offset = computeRelativeDisplacement(lattice, equilibrium);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                T rhoBar;
                Array<T,Descriptor<T>::d> j;
                Cell<T,Descriptor> const& cell = lattice.get(iX,iY,iZ);
                cell.getDynamics().computeRhoBarJ(cell, rhoBar, j);
                T jSqr = normSqr(j);
                for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                    equilibrium.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iPop] =
                        cell.computeEquilibrium(iPop, rhoBar, j, jSqr);
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxAllEquilibriumFunctional3D<T,Descriptor>* BoxAllEquilibriumFunctional3D<T,Descriptor>::clone() const
{
    return new BoxAllEquilibriumFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void BoxAllEquilibriumFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT BoxAllEquilibriumFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor> 
BoxAllPopulationsToLatticeFunctional3D<T,Descriptor>::BoxAllPopulationsToLatticeFunctional3D()
{ }

template<typename T, template<typename U> class Descriptor> 
void BoxAllPopulationsToLatticeFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice, TensorField3D<T,Descriptor<T>::q>& tensorField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                    lattice.get(iX,iY,iZ)[iPop] = tensorField.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iPop];
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxAllPopulationsToLatticeFunctional3D<T,Descriptor>* BoxAllPopulationsToLatticeFunctional3D<T,Descriptor>::clone() const
{
    return new BoxAllPopulationsToLatticeFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void BoxAllPopulationsToLatticeFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT BoxAllPopulationsToLatticeFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor> 
void BoxOmegaFunctional3D<T,Descriptor>::process (
    Box3D domain, BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& scalarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Cell<T,Descriptor> const& cell = lattice.get(iX,iY,iZ);
                scalarField.get(iX+offset.x,iY+offset.y,iZ+offset.z) = cell.getDynamics().getOmega();
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
BoxOmegaFunctional3D<T,Descriptor>* BoxOmegaFunctional3D<T,Descriptor>::clone() const
{
    return new BoxOmegaFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void BoxOmegaFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT BoxOmegaFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, template<typename U> class Descriptor>
void BoxKinematicViscosityFunctional3D<T,Descriptor>::process (
    Box3D domain, BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& scalarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Cell<T,Descriptor> const& cell = lattice.get(iX,iY,iZ);
                scalarField.get(iX+offset.x,iY+offset.y,iZ+offset.z) = Descriptor<T>::cs2*((T)1/cell.getDynamics().getOmega()-(T)0.5);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
BoxKinematicViscosityFunctional3D<T,Descriptor>* BoxKinematicViscosityFunctional3D<T,Descriptor>::clone() const
{
    return new BoxKinematicViscosityFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void BoxKinematicViscosityFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT BoxKinematicViscosityFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor>
void BoxExternalForceFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice, TensorField3D<T,Descriptor<T>::d>& tensorField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        plint oX = iX + offset.x;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            plint oY = iY + offset.y;
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                T *force = lattice.get(iX,iY,iZ).getExternal(Descriptor<T>::ExternalField::forceBeginsAt);
                tensorField.get(oX,oY,iZ+offset.z).from_cArray(force);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
BoxExternalForceFunctional3D<T,Descriptor>* BoxExternalForceFunctional3D<T,Descriptor>::clone() const
{
    return new BoxExternalForceFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void BoxExternalForceFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}


template<typename T, template<typename U> class Descriptor>
BoxExternalScalarFunctional3D<T,Descriptor>::BoxExternalScalarFunctional3D(int whichScalar_)
    : whichScalar(whichScalar_)
{
    PLB_ASSERT(whichScalar < Descriptor<T>::ExternalField::numScalars);
}

template<typename T, template<typename U> class Descriptor>
void BoxExternalScalarFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& scalarField)
{
    Dot3D offset = computeRelativeDisplacement(lattice, scalarField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        plint oX = iX + offset.x;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            plint oY = iY + offset.y;
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                scalarField.get(oX,oY,iZ+offset.z)
                    = *lattice.get(iX,iY,iZ).getExternal(whichScalar);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
BoxExternalScalarFunctional3D<T,Descriptor>* BoxExternalScalarFunctional3D<T,Descriptor>::clone() const
{
    return new BoxExternalScalarFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void BoxExternalScalarFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

/* *************** PART II ******************************************* */
/* *************** Analysis of the scalar-field ********************** */
/* ******************************************************************* */

/* *************** Reductive Data Functionals for scalar-field ******* */

template<typename T>
BoxScalarSumFunctional3D<T>::BoxScalarSumFunctional3D()
    : sumScalarId(this->getStatistics().subscribeSum())
{ }

template<typename T>
void BoxScalarSumFunctional3D<T>::process (
        Box3D domain, ScalarField3D<T>& scalarField )
{
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                statistics.gatherSum(sumScalarId, (double)scalarField.get(iX,iY,iZ));
            }
        }
    }
}

template<typename T>
BoxScalarSumFunctional3D<T>* BoxScalarSumFunctional3D<T>::clone() const
{
    return new BoxScalarSumFunctional3D<T>(*this);
}

template<typename T>
T BoxScalarSumFunctional3D<T>::getSumScalar() const {
    double doubleSum = this->getStatistics().getSum(sumScalarId);
    // The sum is internally computed on floating-point values. If T is
    //   integer, the value must be rounded at the end.
    if (std::numeric_limits<T>::is_integer) {
        return (T) util::roundToInt(doubleSum);
    }
    return (T) doubleSum;
}


template<typename T>
MaskedBoxScalarAverageFunctional3D<T>::MaskedBoxScalarAverageFunctional3D(int flag_)
    : averageScalarId(this->getStatistics().subscribeAverage()),
      flag(flag_)
{ }

template<typename T>
void MaskedBoxScalarAverageFunctional3D<T>::process (
        Box3D domain,
        ScalarField3D<T>& scalarField,
        ScalarField3D<int>& mask )
{
    Dot3D offset = computeRelativeDisplacement(scalarField, mask);
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (mask.get(iX+offset.x, iY+offset.y, iZ+offset.z)==flag) {
                    statistics.gatherAverage(averageScalarId, (double)scalarField.get(iX,iY,iZ));
                    statistics.incrementStats();
                }
            }
        }
    }
}

template<typename T>
MaskedBoxScalarAverageFunctional3D<T>* MaskedBoxScalarAverageFunctional3D<T>::clone() const
{
    return new MaskedBoxScalarAverageFunctional3D<T>(*this);
}

template<typename T>
T MaskedBoxScalarAverageFunctional3D<T>::getAverageScalar() const {
    double doubleAverage = this->getStatistics().getAverage(averageScalarId);
    // The average is internally computed on floating-point values. If T is
    //   integer, the value must be rounded at the end.
    if (std::numeric_limits<T>::is_integer) {
        return (T) util::roundToInt(doubleAverage);
    }
    return (T) doubleAverage;
}


template<typename T>
BoxScalarMinFunctional3D<T>::BoxScalarMinFunctional3D()
    : maxScalarId(this->getStatistics().subscribeMax())
{ }

template<typename T>
void BoxScalarMinFunctional3D<T>::process (
        Box3D domain, ScalarField3D<T>& scalarField )
{
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                // BlockStatistics computes only maximum, no minimum. Therefore,
                //   the relation min(x) = -max(-x) is used.
                statistics.gatherMax(maxScalarId, -(double)scalarField.get(iX,iY,iZ));
            }
        }
    }
}

template<typename T>
BoxScalarMinFunctional3D<T>* BoxScalarMinFunctional3D<T>::clone() const
{
    return new BoxScalarMinFunctional3D<T>(*this);
}

template<typename T>
T BoxScalarMinFunctional3D<T>::getMinScalar() const {
    // The minus sign accounts for the relation min(x) = -max(-x).
    // The sum is internally computed on floating-point values. If T is
    //   integer, the value must be rounded at the end.
    double doubleMin = - this->getStatistics().getMax(maxScalarId);
    if (std::numeric_limits<T>::is_integer) {
        return (T) util::roundToInt(doubleMin);
    }
    return (T) doubleMin;
}


template<typename T>
BoxScalarMaxFunctional3D<T>::BoxScalarMaxFunctional3D()
    : maxScalarId(this->getStatistics().subscribeMax())
{ }

template<typename T>
void BoxScalarMaxFunctional3D<T>::process (
        Box3D domain, ScalarField3D<T>& scalarField )
{
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                statistics.gatherMax(maxScalarId, (double)scalarField.get(iX,iY,iZ));
            }
        }
    }
}

template<typename T>
BoxScalarMaxFunctional3D<T>* BoxScalarMaxFunctional3D<T>::clone() const
{
    return new BoxScalarMaxFunctional3D<T>(*this);
}

template<typename T>
T BoxScalarMaxFunctional3D<T>::getMaxScalar() const {
    // The sum is internally computed on floating-point values. If T is
    //   integer, the value must be rounded at the end.
    double doubleMax = this->getStatistics().getMax(maxScalarId);
    if (std::numeric_limits<T>::is_integer) {
        return (T) util::roundToInt(doubleMax);
    }
    return (T) doubleMax;
}


template<typename T>
BoundedBoxScalarSumFunctional3D<T>::BoundedBoxScalarSumFunctional3D()
    : sumScalarId(this->getStatistics().subscribeSum())
{ }

template<typename T>
void BoundedBoxScalarSumFunctional3D<T>::processBulk (
        Box3D domain, ScalarField3D<T>& scalarField )
{
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                statistics.gatherSum(sumScalarId, (double)scalarField.get(iX,iY,iZ));
            }
        }
    }
}

template<typename T>
void BoundedBoxScalarSumFunctional3D<T>::processPlane (
        int direction, int orientation,
        Box3D domain, ScalarField3D<T>& scalarField )
{
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                // Plane boundary nodes have a weight of 0.5, because only 50% of the
                //   cell centered at the node is inside the computational domain.
                statistics.gatherSum(sumScalarId, (double)scalarField.get(iX,iY,iZ) / 2.);
            }
        }
    }
}

template<typename T>
void BoundedBoxScalarSumFunctional3D<T>::processEdge (
        int plane, int normal1, int normal2,
        Box3D domain, ScalarField3D<T>& scalarField )
{
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                // Edge nodes have a weight of 0.25, because only 25% of the
                //   cell centered at the node is inside the computational domain.
                statistics.gatherSum(sumScalarId, (double)scalarField.get(iX,iY,iZ) / 4.);
            }
        }
    }
}

template<typename T>
void BoundedBoxScalarSumFunctional3D<T>::processCorner (
        int normalX, int normalY, int normalZ,
        Box3D domain, ScalarField3D<T>& scalarField )
{
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                // Corner nodes have a weight of 0.125, because only 1/8 of the
                //   cell centered at the node is inside the computational domain.
                statistics.gatherSum(sumScalarId, (double)scalarField.get(iX,iY,iZ) / 8.);
            }
        }
    }
}

template<typename T>
BoundedBoxScalarSumFunctional3D<T>* BoundedBoxScalarSumFunctional3D<T>::clone() const
{
    return new BoundedBoxScalarSumFunctional3D<T>(*this);
}

template<typename T>
T BoundedBoxScalarSumFunctional3D<T>::getSumScalar() const {
    double doubleSum = this->getStatistics().getSum(sumScalarId);
    // The sum is internally computed on floating-point values. If T is
    //   integer, the value must be rounded at the end.
    if (std::numeric_limits<T>::is_integer) {
        return (T) util::roundToInt(doubleSum);
    }
    return (T) doubleSum;
}

template<typename T1, typename T2>
void CopyConvertScalarFunctional3D<T1,T2>::process (
            Box3D domain, ScalarField3D<T1>& field1,
                      ScalarField3D<T2>& field2 )
{
    Dot3D offset = computeRelativeDisplacement(field1, field2);

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                field2.get(iX+offset.x,iY+offset.y,iZ+offset.z) =
                    (T2) field1.get(iX,iY,iZ);
            }
        }
    }
}

template<typename T1, typename T2>
CopyConvertScalarFunctional3D<T1,T2>* CopyConvertScalarFunctional3D<T1,T2>::clone() const
{
    return new CopyConvertScalarFunctional3D<T1,T2>(*this);
}

template<typename T1, typename T2>
void CopyConvertScalarFunctional3D<T1,T2>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T1, typename T2>
BlockDomain::DomainT CopyConvertScalarFunctional3D<T1,T2>::appliesTo() const {
    return BlockDomain::bulk;
}



/* *************** Data Functionals for scalar-fields **************** */

template<typename T>
void ExtractScalarSubDomainFunctional3D<T>::process (
        Box3D domain, ScalarField3D<T>& field1,
                      ScalarField3D<T>& field2 )
{
    Dot3D offset = computeRelativeDisplacement(field1, field2);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                field2.get(iX+offset.x,iY+offset.y,iZ+offset.z) = field1.get(iX,iY,iZ);
            }
        }
    }
}

template<typename T>
ExtractScalarSubDomainFunctional3D<T>* ExtractScalarSubDomainFunctional3D<T>::clone() const
{
    return new ExtractScalarSubDomainFunctional3D<T>(*this);
}

template<typename T>
void ExtractScalarSubDomainFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT ExtractScalarSubDomainFunctional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

/* ******** compute sqrt functional 3D ************************************* */

template<typename T>
void ComputeScalarSqrtFunctional3D<T>::process (
        Box3D domain, ScalarField3D<T>& A,
                      ScalarField3D<T>& B )
{
    Dot3D offset = computeRelativeDisplacement(A, B);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                PLB_ASSERT( A.get(iX,iY,iZ) >= T());
                B.get(iX+offset.x,iY+offset.y,iZ+offset.z) = std::sqrt( A.get(iX,iY,iZ) );
            }
        }
    }
}

template<typename T>
ComputeScalarSqrtFunctional3D<T>* ComputeScalarSqrtFunctional3D<T>::clone() const
{
    return new ComputeScalarSqrtFunctional3D<T>(*this);
}

template<typename T>
void ComputeScalarSqrtFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT ComputeScalarSqrtFunctional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** compute sqrt functional 3D ************************************* */

template<typename T>
void ComputeAbsoluteValueFunctional3D<T>::process (
        Box3D domain, ScalarField3D<T>& A,
                      ScalarField3D<T>& B )
{
    Dot3D offset = computeRelativeDisplacement(A, B);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                B.get(iX+offset.x,iY+offset.y,iZ+offset.z) = fabs( A.get(iX,iY,iZ) );
            }
        }
    }
}

template<typename T>
ComputeAbsoluteValueFunctional3D<T>* ComputeAbsoluteValueFunctional3D<T>::clone() const
{
    return new ComputeAbsoluteValueFunctional3D<T>(*this);
}

template<typename T>
void ComputeAbsoluteValueFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT ComputeAbsoluteValueFunctional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

/* ******** compute sqrt functional 3D ************************************* */

template<typename T, int nDim>
void ComputeTensorSqrtFunctional3D<T,nDim>::process (
        Box3D domain, TensorField3D<T,nDim>& A,
                      TensorField3D<T,nDim>& B )
{
    Dot3D offset = computeRelativeDisplacement(A, B);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                for (plint iD = 0; iD < nDim; ++iD) {
                    PLB_ASSERT( A.get(iX,iY,iZ)[iD] >= T());
                    B.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iD] = std::sqrt( A.get(iX,iY,iZ)[iD] );
                }
            }
        }
    }
}

template<typename T, int nDim>
ComputeTensorSqrtFunctional3D<T,nDim>* ComputeTensorSqrtFunctional3D<T,nDim>::clone() const
{
    return new ComputeTensorSqrtFunctional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void ComputeTensorSqrtFunctional3D<T,nDim>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, int nDim>
BlockDomain::DomainT ComputeTensorSqrtFunctional3D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_lt_alpha_functional3D ************************************* */

template<typename T>
A_lt_alpha_functional3D<T>::A_lt_alpha_functional3D(T alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_lt_alpha_functional3D<T>::process (
        Box3D domain, ScalarField3D<T>& A, ScalarField3D<int>& result )
{
    Dot3D offset = computeRelativeDisplacement(A, result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                result.get(iX+offset.x,iY+offset.y,iZ+offset.z) =
                    A.get(iX,iY,iZ) < alpha ? 1 : 0;
            }
        }
    }
}

template<typename T>
A_lt_alpha_functional3D<T>* A_lt_alpha_functional3D<T>::clone() const {
    return new A_lt_alpha_functional3D<T>(*this);
}

template<typename T>
void A_lt_alpha_functional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_lt_alpha_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_gt_alpha_functional3D ************************************* */

template<typename T>
A_gt_alpha_functional3D<T>::A_gt_alpha_functional3D(T alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_gt_alpha_functional3D<T>::process (
        Box3D domain, ScalarField3D<T>& A, ScalarField3D<int>& result )
{
    Dot3D offset = computeRelativeDisplacement(A, result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                result.get(iX+offset.x,iY+offset.y,iZ+offset.z) =
                    A.get(iX,iY,iZ) > alpha ? 1 : 0;
            }
        }
    }
}

template<typename T>
A_gt_alpha_functional3D<T>* A_gt_alpha_functional3D<T>::clone() const {
    return new A_gt_alpha_functional3D<T>(*this);
}

template<typename T>
void A_gt_alpha_functional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_gt_alpha_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_plus_alpha_functional3D ************************************* */

template<typename T>
A_plus_alpha_functional3D<T>::A_plus_alpha_functional3D(T alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_plus_alpha_functional3D<T>::process (
        Box3D domain, ScalarField3D<T>& A, ScalarField3D<T>& result )
{
    Dot3D offset = computeRelativeDisplacement(A, result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                result.get(iX+offset.x,iY+offset.y,iZ+offset.z) = A.get(iX,iY,iZ) + alpha;
            }
        }
    }
}

template<typename T>
A_plus_alpha_functional3D<T>* A_plus_alpha_functional3D<T>::clone() const {
    return new A_plus_alpha_functional3D<T>(*this);
}

template<typename T>
void A_plus_alpha_functional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_plus_alpha_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_minus_alpha_functional3D ************************************** */

template<typename T>
A_minus_alpha_functional3D<T>::A_minus_alpha_functional3D(T alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_minus_alpha_functional3D<T>::process (
        Box3D domain, ScalarField3D<T>& A, ScalarField3D<T>& result )
{
    Dot3D offset = computeRelativeDisplacement(A, result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                result.get(iX+offset.x,iY+offset.y,iZ+offset.z) = A.get(iX,iY,iZ) - alpha;
            }
        }
    }
}

template<typename T>
A_minus_alpha_functional3D<T>* A_minus_alpha_functional3D<T>::clone() const {
    return new A_minus_alpha_functional3D<T>(*this);
}

template<typename T>
void A_minus_alpha_functional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_minus_alpha_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Alpha_minus_A_functional3D ************************************* */

template<typename T>
Alpha_minus_A_functional3D<T>::Alpha_minus_A_functional3D(T alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void Alpha_minus_A_functional3D<T>::process (
        Box3D domain, ScalarField3D<T>& A, ScalarField3D<T>& result )
{
    Dot3D offset = computeRelativeDisplacement(A, result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                result.get(iX+offset.x,iY+offset.y,iZ+offset.z) = alpha - A.get(iX,iY,iZ);
            }
        }
    }
}

template<typename T>
Alpha_minus_A_functional3D<T>* Alpha_minus_A_functional3D<T>::clone() const {
    return new Alpha_minus_A_functional3D<T>(*this);
}

template<typename T>
void Alpha_minus_A_functional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT Alpha_minus_A_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}



/* ******** A_times_alpha_functional3D ************************************* */

template<typename T>
A_times_alpha_functional3D<T>::A_times_alpha_functional3D(T alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_times_alpha_functional3D<T>::process (
        Box3D domain, ScalarField3D<T>& A, ScalarField3D<T>& result )
{
    Dot3D offset = computeRelativeDisplacement(A, result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                result.get(iX+offset.x,iY+offset.y,iZ+offset.z) = A.get(iX,iY,iZ) * alpha;
            }
        }
    }
}

template<typename T>
A_times_alpha_functional3D<T>* A_times_alpha_functional3D<T>::clone() const {
    return new A_times_alpha_functional3D<T>(*this);
}

template<typename T>
void A_times_alpha_functional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_times_alpha_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_dividedBy_alpha_functional3D ************************************* */

template<typename T>
A_dividedBy_alpha_functional3D<T>::A_dividedBy_alpha_functional3D(T alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_dividedBy_alpha_functional3D<T>::process (
        Box3D domain, ScalarField3D<T>& A, ScalarField3D<T>& result )
{
    Dot3D offset = computeRelativeDisplacement(A, result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                result.get(iX+offset.x,iY+offset.y,iZ+offset.z) = A.get(iX,iY,iZ) / alpha;
            }
        }
    }
}

template<typename T>
A_dividedBy_alpha_functional3D<T>* A_dividedBy_alpha_functional3D<T>::clone() const {
    return new A_dividedBy_alpha_functional3D<T>(*this);
}

template<typename T>
void A_dividedBy_alpha_functional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_dividedBy_alpha_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Alpha_dividedBy_A_functional3D ************************************* */

template<typename T>
Alpha_dividedBy_A_functional3D<T>::Alpha_dividedBy_A_functional3D(T alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void Alpha_dividedBy_A_functional3D<T>::process (
        Box3D domain, ScalarField3D<T>& A, ScalarField3D<T>& result )
{
    Dot3D offset = computeRelativeDisplacement(A, result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                result.get(iX+offset.x,iY+offset.y,iZ+offset.z) = alpha / A.get(iX,iY,iZ);
            }
        }
    }
}

template<typename T>
Alpha_dividedBy_A_functional3D<T>* Alpha_dividedBy_A_functional3D<T>::clone() const {
    return new Alpha_dividedBy_A_functional3D<T>(*this);
}

template<typename T>
void Alpha_dividedBy_A_functional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT Alpha_dividedBy_A_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}



/* ******** A_plus_alpha_inplace_functional3D ************************************* */

template<typename T>
A_plus_alpha_inplace_functional3D<T>::A_plus_alpha_inplace_functional3D(T alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_plus_alpha_inplace_functional3D<T>::process (
        Box3D domain, ScalarField3D<T>& A)
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                A.get(iX,iY,iZ) += alpha;
            }
        }
    }
}

template<typename T>
A_plus_alpha_inplace_functional3D<T>* A_plus_alpha_inplace_functional3D<T>::clone() const {
    return new A_plus_alpha_inplace_functional3D<T>(*this);
}

template<typename T>
void A_plus_alpha_inplace_functional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_plus_alpha_inplace_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_minus_alpha_inplace_functional3D ************************************** */

template<typename T>
A_minus_alpha_inplace_functional3D<T>::A_minus_alpha_inplace_functional3D(T alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_minus_alpha_inplace_functional3D<T>::process (
        Box3D domain, ScalarField3D<T>& A)
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                A.get(iX,iY,iZ) -= alpha;
            }
        }
    }
}

template<typename T>
A_minus_alpha_inplace_functional3D<T>* A_minus_alpha_inplace_functional3D<T>::clone() const {
    return new A_minus_alpha_inplace_functional3D<T>(*this);
}

template<typename T>
void A_minus_alpha_inplace_functional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_minus_alpha_inplace_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_times_alpha_inplace_functional3D ************************************* */

template<typename T>
A_times_alpha_inplace_functional3D<T>::A_times_alpha_inplace_functional3D(T alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_times_alpha_inplace_functional3D<T>::process (
        Box3D domain, ScalarField3D<T>& A )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                A.get(iX,iY,iZ) *= alpha;
            }
        }
    }
}

template<typename T>
A_times_alpha_inplace_functional3D<T>* A_times_alpha_inplace_functional3D<T>::clone() const {
    return new A_times_alpha_inplace_functional3D<T>(*this);
}

template<typename T>
void A_times_alpha_inplace_functional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_times_alpha_inplace_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_dividedBy_alpha_inplace_functional3D ************************************* */

template<typename T>
A_dividedBy_alpha_inplace_functional3D<T>::A_dividedBy_alpha_inplace_functional3D(T alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_dividedBy_alpha_inplace_functional3D<T>::process (
        Box3D domain, ScalarField3D<T>& A )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                A.get(iX,iY,iZ) /= alpha;
            }
        }
    }
}

template<typename T>
A_dividedBy_alpha_inplace_functional3D<T>* A_dividedBy_alpha_inplace_functional3D<T>::clone() const {
    return new A_dividedBy_alpha_inplace_functional3D<T>(*this);
}

template<typename T>
void A_dividedBy_alpha_inplace_functional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_dividedBy_alpha_inplace_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_lt_B_functional3D ****************************************** */

template<typename T>
void A_lt_B_functional3D<T>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> fields )
{
    PLB_PRECONDITION( fields.size()==3 );
    ScalarField3D<T>& A = *dynamic_cast<ScalarField3D<T>*>(fields[0]);
    ScalarField3D<T>& B = *dynamic_cast<ScalarField3D<T>*>(fields[1]);
    ScalarField3D<int>& result = *dynamic_cast<ScalarField3D<int>*>(fields[2]);
    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                    = A.get(iX,iY,iZ) < B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z) ? 1:0;
            }
        }
    }
}

template<typename T>
A_lt_B_functional3D<T>* A_lt_B_functional3D<T>::clone() const {
    return new A_lt_B_functional3D<T>(*this);
}

template<typename T>
void A_lt_B_functional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_lt_B_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_gt_B_functional3D ****************************************** */

template<typename T>
void A_gt_B_functional3D<T>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> fields )
{
    PLB_PRECONDITION( fields.size()==3 );
    ScalarField3D<T>& A = *dynamic_cast<ScalarField3D<T>*>(fields[0]);
    ScalarField3D<T>& B = *dynamic_cast<ScalarField3D<T>*>(fields[1]);
    ScalarField3D<int>& result = *dynamic_cast<ScalarField3D<int>*>(fields[2]);
    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                    = A.get(iX,iY,iZ) > B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z) ? 1:0;
            }
        }
    }
}

template<typename T>
A_gt_B_functional3D<T>* A_gt_B_functional3D<T>::clone() const {
    return new A_gt_B_functional3D<T>(*this);
}

template<typename T>
void A_gt_B_functional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_gt_B_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_plus_B_functional3D ****************************************** */

template<typename T>
void A_plus_B_functional3D<T>::process (
        Box3D domain, std::vector<ScalarField3D<T>*> fields )
{
    PLB_PRECONDITION( fields.size()==3 );
    ScalarField3D<T>& A = *fields[0];
    ScalarField3D<T>& B = *fields[1];
    ScalarField3D<T>& result = *fields[2];
    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                    = A.get(iX,iY,iZ) + B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z);
            }
        }
    }
}

template<typename T>
A_plus_B_functional3D<T>* A_plus_B_functional3D<T>::clone() const {
    return new A_plus_B_functional3D<T>(*this);
}

template<typename T>
void A_plus_B_functional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_plus_B_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_minus_B_functional3D ****************************************** */

template<typename T>
void A_minus_B_functional3D<T>::process (
        Box3D domain, std::vector<ScalarField3D<T>*> fields )
{
    PLB_PRECONDITION( fields.size()==3 );
    ScalarField3D<T>& A = *fields[0];
    ScalarField3D<T>& B = *fields[1];
    ScalarField3D<T>& result = *fields[2];
    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                    = A.get(iX,iY,iZ) - B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z);
            }
        }
    }
}

template<typename T>
A_minus_B_functional3D<T>* A_minus_B_functional3D<T>::clone() const {
    return new A_minus_B_functional3D<T>(*this);
}

template<typename T>
void A_minus_B_functional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_minus_B_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_times_B_functional3D ****************************************** */

template<typename T>
void A_times_B_functional3D<T>::process (
        Box3D domain, std::vector<ScalarField3D<T>*> fields )
{
    PLB_PRECONDITION( fields.size()==3 );
    ScalarField3D<T>& A = *fields[0];
    ScalarField3D<T>& B = *fields[1];
    ScalarField3D<T>& result = *fields[2];
    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                    = A.get(iX,iY,iZ) * B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z);
            }
        }
    }
}

template<typename T>
A_times_B_functional3D<T>* A_times_B_functional3D<T>::clone() const {
    return new A_times_B_functional3D<T>(*this);
}

template<typename T>
void A_times_B_functional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_times_B_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_dividedBy_B_functional3D ****************************************** */

template<typename T>
void A_dividedBy_B_functional3D<T>::process (
        Box3D domain, std::vector<ScalarField3D<T>*> fields )
{
    PLB_PRECONDITION( fields.size()==3 );
    ScalarField3D<T>& A = *fields[0];
    ScalarField3D<T>& B = *fields[1];
    ScalarField3D<T>& result = *fields[2];
    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                    = A.get(iX,iY,iZ) / B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z);
            }
        }
    }
}

template<typename T>
A_dividedBy_B_functional3D<T>* A_dividedBy_B_functional3D<T>::clone() const {
    return new A_dividedBy_B_functional3D<T>(*this);
}

template<typename T>
void A_dividedBy_B_functional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_dividedBy_B_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_plus_B_inplace_functional3D ****************************************** */

template<typename T>
void A_plus_B_inplace_functional3D<T>::process (
        Box3D domain, ScalarField3D<T>& A, ScalarField3D<T>& B)
{
    Dot3D offset = computeRelativeDisplacement(A,B);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                A.get(iX,iY,iZ) += B.get(iX+offset.x,iY+offset.y,iZ+offset.z);
            }
        }
    }
}

template<typename T>
A_plus_B_inplace_functional3D<T>* A_plus_B_inplace_functional3D<T>::clone() const {
    return new A_plus_B_inplace_functional3D<T>(*this);
}

template<typename T>
void A_plus_B_inplace_functional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT A_plus_B_inplace_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_minus_B_inplace_functional3D ****************************************** */

template<typename T>
void A_minus_B_inplace_functional3D<T>::process (
        Box3D domain, ScalarField3D<T>& A, ScalarField3D<T>& B)
{
    Dot3D offset = computeRelativeDisplacement(A,B);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                A.get(iX,iY,iZ) -= B.get(iX+offset.x,iY+offset.y,iZ+offset.z);
            }
        }
    }
}

template<typename T>
A_minus_B_inplace_functional3D<T>* A_minus_B_inplace_functional3D<T>::clone() const {
    return new A_minus_B_inplace_functional3D<T>(*this);
}

template<typename T>
void A_minus_B_inplace_functional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT A_minus_B_inplace_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_times_B_inplace_functional3D ****************************************** */

template<typename T>
void A_times_B_inplace_functional3D<T>::process (
        Box3D domain, ScalarField3D<T>& A, ScalarField3D<T>& B)
{
    Dot3D offset = computeRelativeDisplacement(A,B);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                A.get(iX,iY,iZ) *= B.get(iX+offset.x,iY+offset.y,iZ+offset.z);
            }
        }
    }
}

template<typename T>
A_times_B_inplace_functional3D<T>* A_times_B_inplace_functional3D<T>::clone() const {
    return new A_times_B_inplace_functional3D<T>(*this);
}

template<typename T>
void A_times_B_inplace_functional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT A_times_B_inplace_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_dividedBy_B_inplace_functional3D ****************************************** */

template<typename T>
void A_dividedBy_B_inplace_functional3D<T>::process (
        Box3D domain, ScalarField3D<T>& A, ScalarField3D<T>& B)
{
    Dot3D offset = computeRelativeDisplacement(A,B);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                A.get(iX,iY,iZ) -= B.get(iX+offset.x,iY+offset.y,iZ+offset.z);
            }
        }
    }
}

template<typename T>
A_dividedBy_B_inplace_functional3D<T>* A_dividedBy_B_inplace_functional3D<T>::clone() const {
    return new A_dividedBy_B_inplace_functional3D<T>(*this);
}

template<typename T>
void A_dividedBy_B_inplace_functional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT A_dividedBy_B_inplace_functional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

template<typename T>
void UniformlyBoundScalarField3D<T>::process (
        Box3D domain, ScalarField3D<T>& data )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                T& val = data.get(iX,iY,iZ);
                if (fabs(val) > bound)
                    val = val > T() ? bound : -bound;
            }
        }
    }
}


/* *************** PART III ****************************************** */
/* *************** Analysis of the tensor-field ********************** */
/* ******************************************************************* */

template<typename T1, typename T2, int nDim>
void CopyConvertTensorFunctional3D<T1,T2,nDim>::process (
        Box3D domain, TensorField3D<T1,nDim>& field1,
                      TensorField3D<T2,nDim>& field2 )
{
    Dot3D offset = computeRelativeDisplacement(field1, field2);

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                for (int iDim=0; iDim<nDim; ++iDim) {
                    field2.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim] =
                        (T2) field1.get(iX,iY,iZ)[iDim];
                }
            }
        }
    }
}

template<typename T1, typename T2, int nDim>
CopyConvertTensorFunctional3D<T1,T2,nDim>* CopyConvertTensorFunctional3D<T1,T2,nDim>::clone() const
{
    return new CopyConvertTensorFunctional3D<T1,T2,nDim>(*this);
}

template<typename T1, typename T2, int nDim>
void CopyConvertTensorFunctional3D<T1,T2,nDim>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T1, typename T2, int nDim>
BlockDomain::DomainT CopyConvertTensorFunctional3D<T1,T2,nDim>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T, int nDim>
void ExtractTensorSubDomainFunctional3D<T,nDim>::process (
        Box3D domain, TensorField3D<T,nDim>& field1,
                      TensorField3D<T,nDim>& field2 )
{
    Dot3D offset = computeRelativeDisplacement(field1, field2);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                for (int iDim=0; iDim<nDim; ++iDim) {
                    field2.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim]
                        = field1.get(iX,iY,iZ)[iDim];
                }
            }
        }
    }
}

template<typename T, int nDim>
ExtractTensorSubDomainFunctional3D<T,nDim>* ExtractTensorSubDomainFunctional3D<T,nDim>::clone() const
{
    return new ExtractTensorSubDomainFunctional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void ExtractTensorSubDomainFunctional3D<T,nDim>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, int nDim>
BlockDomain::DomainT ExtractTensorSubDomainFunctional3D<T,nDim>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T, int nDim>
ExtractTensorComponentFunctional3D<T,nDim>::ExtractTensorComponentFunctional3D(int iComponent_)
    : iComponent(iComponent_)
{ }

template<typename T, int nDim>
void ExtractTensorComponentFunctional3D<T,nDim>::process (
        Box3D domain, ScalarField3D<T>& scalarField,
                      TensorField3D<T,nDim>& tensorField )
{
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                scalarField.get(iX,iY,iZ)
                    = tensorField.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iComponent];
            }
        }
    }
}

template<typename T, int nDim>
ExtractTensorComponentFunctional3D<T,nDim>* ExtractTensorComponentFunctional3D<T,nDim>::clone() const
{
    return new ExtractTensorComponentFunctional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void ExtractTensorComponentFunctional3D<T,nDim>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T, int nDim>
BlockDomain::DomainT ExtractTensorComponentFunctional3D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, int nDim>
void ComputeNormFunctional3D<T,nDim>::process (
        Box3D domain, ScalarField3D<T>& scalarField,
                      TensorField3D<T,nDim>& tensorField )
{
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                scalarField.get(iX,iY,iZ) = std::sqrt( VectorTemplateImpl<T,nDim>::normSqr (
                                                           tensorField.get(iX+offset.x,iY+offset.y,iZ+offset.z) ) );
            }
        }
    }
}

template<typename T, int nDim>
ComputeNormFunctional3D<T,nDim>* ComputeNormFunctional3D<T,nDim>::clone() const
{
    return new ComputeNormFunctional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void ComputeNormFunctional3D<T,nDim>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T, int nDim>
BlockDomain::DomainT ComputeNormFunctional3D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T, int nDim>
void ComputeNormSqrFunctional3D<T,nDim>::process (
        Box3D domain, ScalarField3D<T>& scalarField,
                      TensorField3D<T,nDim>& tensorField )
{
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                scalarField.get(iX,iY,iZ) = VectorTemplateImpl<T,nDim>::normSqr (
                                            tensorField.get(iX+offset.x,iY+offset.y,iZ+offset.z) );
            }
        }
    }
}

template<typename T, int nDim>
ComputeNormSqrFunctional3D<T,nDim>* ComputeNormSqrFunctional3D<T,nDim>::clone() const
{
    return new ComputeNormSqrFunctional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void ComputeNormSqrFunctional3D<T,nDim>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T, int nDim>
BlockDomain::DomainT ComputeNormSqrFunctional3D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T>
void ComputeSymmetricTensorNormFunctional3D<T>::process (
        Box3D domain, ScalarField3D<T>& scalarField,
                      TensorField3D<T,6>& tensorField )
{
    typedef SymmetricTensorImpl<T,3> tensor;
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Array<T,6>& el = tensorField.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                scalarField.get(iX,iY,iZ) = std::sqrt ( 
                        // Count diagonal components once ...
                                util::sqr(el[tensor::xx]) + util::sqr(el[tensor::yy]) + util::sqr(el[tensor::zz]) +
                        // .. and off-diagonal components twice, due to symmetry.
                        (T)2 * (util::sqr(el[tensor::xy]) + util::sqr(el[tensor::xz]) +util::sqr(el[tensor::yz])) );
            }
        }
    }
}

template<typename T>
ComputeSymmetricTensorNormFunctional3D<T>* ComputeSymmetricTensorNormFunctional3D<T>::clone() const
{
    return new ComputeSymmetricTensorNormFunctional3D<T>(*this);
}

template<typename T>
void ComputeSymmetricTensorNormFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT ComputeSymmetricTensorNormFunctional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T>
void ComputeSymmetricTensorNormSqrFunctional3D<T>::process (
        Box3D domain, ScalarField3D<T>& scalarField,
                      TensorField3D<T,6>& tensorField )
{
    typedef SymmetricTensorImpl<T,3> tensor;
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Array<T,6>& el = tensorField.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                scalarField.get(iX,iY,iZ) = 
                        // Count diagonal components once ...
                                util::sqr(el[tensor::xx]) + util::sqr(el[tensor::yy]) + util::sqr(el[tensor::zz]) +
                        // .. and off-diagonal components twice, due to symmetry.
                        (T)2 * (util::sqr(el[tensor::xy]) + util::sqr(el[tensor::xz]) +util::sqr(el[tensor::yz]));
            }
        }
    }
}

template<typename T>
ComputeSymmetricTensorNormSqrFunctional3D<T>* ComputeSymmetricTensorNormSqrFunctional3D<T>::clone() const
{
    return new ComputeSymmetricTensorNormSqrFunctional3D<T>(*this);
}

template<typename T>
void ComputeSymmetricTensorNormSqrFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT ComputeSymmetricTensorNormSqrFunctional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}



template<typename T>
void ComputeSymmetricTensorTraceFunctional3D<T>::process (
        Box3D domain, ScalarField3D<T>& scalarField,
                      TensorField3D<T,6>& tensorField )
{
    typedef SymmetricTensorImpl<T,3> tensor;
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Array<T,6>& el = tensorField.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                scalarField.get(iX,iY,iZ) = el[tensor::xx] + el[tensor::yy] + el[tensor::zz];
            }
        }
    }
}

template<typename T>
ComputeSymmetricTensorTraceFunctional3D<T>* ComputeSymmetricTensorTraceFunctional3D<T>::clone() const
{
    return new ComputeSymmetricTensorTraceFunctional3D<T>(*this);
}

template<typename T>
void ComputeSymmetricTensorTraceFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT ComputeSymmetricTensorTraceFunctional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T>
void BoxBulkGradientFunctional3D<T>::process (
        Box3D domain, ScalarField3D<T>& phi,
                      TensorField3D<T,3>& gradient )
{
    Dot3D offset = computeRelativeDisplacement(phi, gradient);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                plint iZ2 = iZ+offset.z;
                gradient.get(iX2,iY2,iZ2)[0] =
                    fdDataField::bulkXderiv(phi, iX,iY,iZ);
                gradient.get(iX2,iY2,iZ2)[1] =
                    fdDataField::bulkYderiv(phi, iX,iY,iZ);
                gradient.get(iX2,iY2,iZ2)[2] =
                    fdDataField::bulkZderiv(phi, iX,iY,iZ);
            }
        }
    }
}

template<typename T>
BoxBulkGradientFunctional3D<T>* BoxBulkGradientFunctional3D<T>::clone() const
{
    return new BoxBulkGradientFunctional3D<T>(*this);
}

template<typename T>
void BoxBulkGradientFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT BoxBulkGradientFunctional3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}



template<typename T>
void BoxGradientFunctional3D<T>::processBulk (
        Box3D domain, ScalarField3D<T>& phi, TensorField3D<T,3>& gradient )
{
    Dot3D offset = computeRelativeDisplacement(phi, gradient);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                plint iZ2 = iZ+offset.z;
                gradient.get(iX2,iY2,iZ2)[0] =
                    fdDataField::bulkXderiv(phi, iX,iY,iZ);
                gradient.get(iX2,iY2,iZ2)[1] =
                    fdDataField::bulkYderiv(phi, iX,iY,iZ);
                gradient.get(iX2,iY2,iZ2)[2] =
                    fdDataField::bulkZderiv(phi, iX,iY,iZ);
            }
        }
    }
}

template<typename T>
void BoxGradientFunctional3D<T>::processPlane (
        int direction, int orientation, Box3D domain,
        ScalarField3D<T>& phi, TensorField3D<T,3>& gradient )
{
    Dot3D offset = computeRelativeDisplacement(phi, gradient);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                plint iZ2 = iZ+offset.z;
                gradient.get(iX2,iY2,iZ2)[0] =
                    fdDataField::planeXderiv(phi,direction,orientation, iX,iY,iZ);
                gradient.get(iX2,iY2,iZ2)[1] =
                    fdDataField::planeYderiv(phi,direction,orientation, iX,iY,iZ);
                gradient.get(iX2,iY2,iZ2)[2] =
                    fdDataField::planeZderiv(phi,direction,orientation, iX,iY,iZ);
            }
        }
    }
}

template<typename T>
void BoxGradientFunctional3D<T>::processEdge (
        int plane, int normal1, int normal2, Box3D domain,
        ScalarField3D<T>& phi, TensorField3D<T,3>& gradient )
{
    Dot3D offset = computeRelativeDisplacement(phi, gradient);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                plint iZ2 = iZ+offset.z;
                gradient.get(iX2,iY2,iZ2)[0] =
                    fdDataField::edgeXderiv(phi,plane,normal1,normal2, iX,iY,iZ);
                gradient.get(iX2,iY2,iZ2)[1] =
                    fdDataField::edgeYderiv(phi,plane,normal1,normal2, iX,iY,iZ);
                gradient.get(iX2,iY2,iZ2)[2] =
                    fdDataField::edgeZderiv(phi,plane,normal1,normal2, iX,iY,iZ);
            }
        }
    }
}

template<typename T>
void BoxGradientFunctional3D<T>::processCorner (
        int normalX, int normalY, int normalZ, Box3D domain,
        ScalarField3D<T>& phi, TensorField3D<T,3>& gradient )
{

    Dot3D offset = computeRelativeDisplacement(phi, gradient);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                plint iZ2 = iZ+offset.z;
                gradient.get(iX2,iY2,iZ2)[0] =
                    fdDataField::cornerXderiv(phi,normalX,normalY,normalZ, iX,iY,iZ);
                gradient.get(iX2,iY2,iZ2)[1] =
                    fdDataField::cornerYderiv(phi,normalX,normalY,normalZ, iX,iY,iZ);
                gradient.get(iX2,iY2,iZ2)[2] =
                    fdDataField::cornerZderiv(phi,normalX,normalY,normalZ, iX,iY,iZ);
            }
        }
    }
}


template<typename T>
BoxGradientFunctional3D<T>* BoxGradientFunctional3D<T>::clone() const
{
    return new BoxGradientFunctional3D<T>(*this);
}

template<typename T>
void BoxGradientFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}


template<typename T>
BlockDomain::DomainT BoxGradientFunctional3D<T>::appliesTo() const {
    // Don't apply to envelope, because nearest neighbors need to be accessed.
    return BlockDomain::bulk;
}


template<typename T, int nDim>
void BoxBulkVorticityFunctional3D<T,nDim>::process (
        Box3D domain, TensorField3D<T,nDim>& velocity,
                      TensorField3D<T,nDim>& vorticity )
{
    Dot3D offset = computeRelativeDisplacement(velocity, vorticity);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        plint iX2 = iX+offset.x;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            plint iY2 = iY+offset.y;
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iZ2 = iZ+offset.z;
                vorticity.get(iX2,iY2,iZ2)[0] =
                    fdDataField::bulkVorticityX(velocity, iX,iY,iZ);
                vorticity.get(iX2,iY2,iZ2)[1] = 
                    fdDataField::bulkVorticityY(velocity, iX,iY,iZ);
                vorticity.get(iX2,iY2,iZ2)[2] = 
                    fdDataField::bulkVorticityZ(velocity, iX,iY,iZ);
            }
        }
    }
}

template<typename T, int nDim>
BoxBulkVorticityFunctional3D<T,nDim>* BoxBulkVorticityFunctional3D<T,nDim>::clone() const
{
    return new BoxBulkVorticityFunctional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void BoxBulkVorticityFunctional3D<T,nDim>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, int nDim>
BlockDomain::DomainT BoxBulkVorticityFunctional3D<T,nDim>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T, int nDim>
void BoxVorticityFunctional3D<T,nDim>::processBulk (
        Box3D domain, TensorField3D<T,nDim>& velocity, TensorField3D<T,nDim>& vorticity )
{
    Dot3D offset = computeRelativeDisplacement(velocity, vorticity);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                plint iZ2 = iZ+offset.z;
                vorticity.get(iX2,iY2,iZ2)[0] =
                    fdDataField::bulkVorticityX(velocity, iX,iY,iZ);
                vorticity.get(iX2,iY2,iZ2)[1] = 
                    fdDataField::bulkVorticityY(velocity, iX,iY,iZ);
                vorticity.get(iX2,iY2,iZ2)[2] = 
                    fdDataField::bulkVorticityZ(velocity, iX,iY,iZ);
            }
        }
    }
}

template<typename T, int nDim>
void BoxVorticityFunctional3D<T,nDim>::processPlane (
        int direction, int orientation, Box3D domain,
        TensorField3D<T,nDim>& velocity, TensorField3D<T,nDim>& vorticity )
{
    Dot3D offset = computeRelativeDisplacement(velocity, vorticity);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                plint iZ2 = iZ+offset.z;
                vorticity.get(iX2,iY2,iZ2)[0] = 
                    fdDataField::planeVorticityX(velocity,direction,orientation, iX,iY,iZ);
                vorticity.get(iX2,iY2,iZ2)[1] = 
                    fdDataField::planeVorticityY(velocity,direction,orientation, iX,iY,iZ);
                vorticity.get(iX2,iY2,iZ2)[2] = 
                    fdDataField::planeVorticityZ(velocity,direction,orientation, iX,iY,iZ);
            }
        }
    }
}

template<typename T, int nDim>
void BoxVorticityFunctional3D<T,nDim>::processEdge (
        int plane, int normal1, int normal2, Box3D domain,
        TensorField3D<T,nDim>& velocity, TensorField3D<T,nDim>& vorticity )
{
    Dot3D offset = computeRelativeDisplacement(velocity, vorticity);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                plint iZ2 = iZ+offset.z;
                vorticity.get(iX2,iY2,iZ2)[0] = 
                    fdDataField::edgeVorticityX(velocity,plane,normal1,normal2, iX,iY,iZ);
                vorticity.get(iX2,iY2,iZ2)[1] = 
                    fdDataField::edgeVorticityY(velocity,plane,normal1,normal2, iX,iY,iZ);
                vorticity.get(iX2,iY2,iZ2)[2] = 
                    fdDataField::edgeVorticityZ(velocity,plane,normal1,normal2, iX,iY,iZ);
            }
        }
    }
}

template<typename T, int nDim>
void BoxVorticityFunctional3D<T,nDim>::processCorner (
        int normalX, int normalY, int normalZ, Box3D domain,
        TensorField3D<T,nDim>& velocity, TensorField3D<T,nDim>& vorticity )
{

    Dot3D offset = computeRelativeDisplacement(velocity, vorticity);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                plint iZ2 = iZ+offset.z;
                vorticity.get(iX2,iY2,iZ2)[0] = 
                    fdDataField::cornerVorticityX(velocity,normalX,normalY,normalZ, iX,iY,iZ);
                vorticity.get(iX2,iY2,iZ2)[1] = 
                    fdDataField::cornerVorticityY(velocity,normalX,normalY,normalZ, iX,iY,iZ);
                vorticity.get(iX2,iY2,iZ2)[2] = 
                    fdDataField::cornerVorticityZ(velocity,normalX,normalY,normalZ, iX,iY,iZ);
            }
        }
    }
}


template<typename T, int nDim>
BoxVorticityFunctional3D<T,nDim>* BoxVorticityFunctional3D<T,nDim>::clone() const
{
    return new BoxVorticityFunctional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void BoxVorticityFunctional3D<T,nDim>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}


template<typename T, int nDim>
BlockDomain::DomainT BoxVorticityFunctional3D<T,nDim>::appliesTo() const {
    // Don't apply to envelope, because nearest neighbors need to be accessed.
    return BlockDomain::bulk;
}


template<typename T, int nDim>
void BoxBulkStrainRateFunctional3D<T,nDim>::process (
        Box3D domain, TensorField3D<T,nDim>& velocity,
                      TensorField3D<T,SymmetricTensorImpl<T,nDim>::n>& S )
{
    typedef SymmetricTensorImpl<T,nDim> tensor;
    Dot3D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                plint iZ2 = iZ+offset.z;
                Array<T,SymmetricTensorImpl<T,nDim>::n>& el = S.get(iX2,iY2,iZ2);
                el[tensor::xx] = fdDataField::bulkXderiv(velocity, iX, iY, iZ, 0);
                el[tensor::xy] = ( fdDataField::bulkXderiv(velocity, iX, iY, iZ, 1) +
                                   fdDataField::bulkYderiv(velocity, iX, iY, iZ, 0) ) / (T)2;
                el[tensor::xz] = ( fdDataField::bulkXderiv(velocity, iX, iY, iZ, 2) +
                                   fdDataField::bulkZderiv(velocity, iX, iY, iZ, 0) ) / (T)2;
                el[tensor::yy] = fdDataField::bulkYderiv(velocity, iX, iY, iZ, 1);
                el[tensor::yz] = ( fdDataField::bulkYderiv(velocity, iX, iY, iZ, 2) +
                                   fdDataField::bulkZderiv(velocity, iX, iY, iZ, 1) ) / (T)2;
                el[tensor::zz] = fdDataField::bulkZderiv(velocity, iX, iY, iZ, 2);
            }
        }
    }
}

template<typename T, int nDim>
BoxBulkStrainRateFunctional3D<T,nDim>* BoxBulkStrainRateFunctional3D<T,nDim>::clone() const
{
    return new BoxBulkStrainRateFunctional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void BoxBulkStrainRateFunctional3D<T,nDim>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, int nDim>
BlockDomain::DomainT BoxBulkStrainRateFunctional3D<T,nDim>::appliesTo() const {
//     return BlockDomain::bulkAndEnvelope;
    return BlockDomain::bulk;
}


template<typename T, int nDim>
void BoxStrainRateFunctional3D<T,nDim>::processBulk (
        Box3D domain, TensorField3D<T,nDim>& velocity, TensorField3D<T,SymmetricTensorImpl<T,nDim>::n>& S )
{
    typedef SymmetricTensorImpl<T,nDim> tensor;
    Dot3D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                plint iZ2 = iZ+offset.z;
                Array<T,SymmetricTensorImpl<T,nDim>::n>& el = S.get(iX2,iY2,iZ2);
                el[tensor::xx] = fdDataField::bulkXderiv(velocity, iX, iY, iZ, 0);
                el[tensor::xy] = ( fdDataField::bulkXderiv(velocity, iX, iY, iZ, 1) +
                                   fdDataField::bulkYderiv(velocity, iX, iY, iZ, 0) ) / (T)2;
                el[tensor::xz] = ( fdDataField::bulkXderiv(velocity, iX, iY, iZ, 2) +
                                   fdDataField::bulkZderiv(velocity, iX, iY, iZ, 0) ) / (T)2;
                el[tensor::yy] = fdDataField::bulkYderiv(velocity, iX, iY, iZ, 1);
                el[tensor::yz] = ( fdDataField::bulkYderiv(velocity, iX, iY, iZ, 2) +
                                   fdDataField::bulkZderiv(velocity, iX, iY, iZ, 1) ) / (T)2;
                el[tensor::zz] = fdDataField::bulkZderiv(velocity, iX, iY, iZ, 2);
            }
        }
    }
}

template<typename T, int nDim>
void BoxStrainRateFunctional3D<T,nDim>::processPlane (
        int direction, int orientation, Box3D domain,
        TensorField3D<T,nDim>& velocity, TensorField3D<T,SymmetricTensorImpl<T,nDim>::n>& S )
{
    typedef SymmetricTensorImpl<T,nDim> tensor;
    Dot3D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                plint iZ2 = iZ+offset.z;
                Array<T,SymmetricTensorImpl<T,nDim>::n>& el = S.get(iX2,iY2,iZ2);
                el[tensor::xx] = fdDataField::planeXderiv(velocity, direction,orientation, iX, iY, iZ, 0);
                el[tensor::xy] = ( fdDataField::planeXderiv(velocity, direction,orientation, iX, iY, iZ, 1) +
                               fdDataField::planeYderiv(velocity, direction,orientation, iX, iY, iZ, 0) ) / (T)2;
                el[tensor::xz] = ( fdDataField::planeXderiv(velocity, direction,orientation, iX, iY, iZ, 2) +
                               fdDataField::planeZderiv(velocity, direction,orientation, iX, iY, iZ, 0) ) / (T)2;
                el[tensor::yy] = fdDataField::planeYderiv(velocity, direction,orientation, iX, iY, iZ, 1);
                el[tensor::yz] = ( fdDataField::planeYderiv(velocity, direction,orientation, iX, iY, iZ, 2) +
                               fdDataField::planeZderiv(velocity, direction,orientation, iX, iY, iZ, 1) ) / (T)2;
                el[tensor::zz] = fdDataField::planeZderiv(velocity, direction,orientation, iX, iY, iZ, 2);
            }
        }
    }
}

template<typename T, int nDim>
void BoxStrainRateFunctional3D<T,nDim>::processEdge (
        int plane, int normal1, int normal2, Box3D domain,
        TensorField3D<T,nDim>& velocity, TensorField3D<T,SymmetricTensorImpl<T,nDim>::n>& S )
{
    typedef SymmetricTensorImpl<T,nDim> tensor;
    Dot3D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                plint iZ2 = iZ+offset.z;
                Array<T,SymmetricTensorImpl<T,nDim>::n>& el = S.get(iX2,iY2,iZ2);
                el[tensor::xx] = fdDataField::edgeXderiv(velocity, plane,normal1,normal2, iX, iY, iZ, 0);
                el[tensor::xy] = ( fdDataField::edgeXderiv(velocity, plane,normal1,normal2, iX, iY, iZ, 1) +
                               fdDataField::edgeYderiv(velocity, plane,normal1,normal2, iX, iY, iZ, 0) ) / (T)2;
                el[tensor::xz] = ( fdDataField::edgeXderiv(velocity, plane,normal1,normal2, iX, iY, iZ, 2) +
                               fdDataField::edgeZderiv(velocity, plane,normal1,normal2, iX, iY, iZ, 0) ) / (T)2;
                el[tensor::yy] = fdDataField::edgeYderiv(velocity, plane,normal1,normal2, iX, iY, iZ, 1);
                el[tensor::yz] = ( fdDataField::edgeYderiv(velocity, plane,normal1,normal2, iX, iY, iZ, 2) +
                               fdDataField::edgeZderiv(velocity, plane,normal1,normal2, iX, iY, iZ, 1) ) / (T)2;
                el[tensor::zz] = fdDataField::edgeZderiv(velocity, plane,normal1,normal2, iX, iY, iZ, 2);
            }
        }
    }
}

template<typename T, int nDim>
void BoxStrainRateFunctional3D<T,nDim>::processCorner (
        int normalX, int normalY, int normalZ, Box3D domain,
        TensorField3D<T,nDim>& velocity, TensorField3D<T,SymmetricTensorImpl<T,nDim>::n>& S )
{

    typedef SymmetricTensorImpl<T,nDim> tensor;
    Dot3D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                plint iZ2 = iZ+offset.z;
                Array<T,SymmetricTensorImpl<T,nDim>::n>& el = S.get(iX2,iY2,iZ2);
                el[tensor::xx] = fdDataField::cornerXderiv(velocity, normalX,normalY,normalZ, iX, iY, iZ, 0);
                el[tensor::xy] = ( fdDataField::cornerXderiv(velocity, normalX,normalY,normalZ, iX, iY, iZ, 1) +
                                   fdDataField::cornerYderiv(velocity, normalX,normalY,normalZ, iX, iY, iZ, 0) ) / (T)2;
                el[tensor::xz] = ( fdDataField::cornerXderiv(velocity, normalX,normalY,normalZ, iX, iY, iZ, 2) +
                                   fdDataField::cornerZderiv(velocity, normalX,normalY,normalZ, iX, iY, iZ, 0) ) / (T)2;
                el[tensor::yy] = fdDataField::cornerYderiv(velocity, normalX,normalY,normalZ, iX, iY, iZ, 1);
                el[tensor::yz] = ( fdDataField::cornerYderiv(velocity, normalX,normalY,normalZ, iX, iY, iZ, 2) +
                                   fdDataField::cornerZderiv(velocity, normalX,normalY,normalZ, iX, iY, iZ, 1) ) / (T)2;
                el[tensor::zz] = fdDataField::cornerZderiv(velocity, normalX,normalY,normalZ, iX, iY, iZ, 2);
            }
        }
    }
}


template<typename T, int nDim>
BoxStrainRateFunctional3D<T,nDim>* BoxStrainRateFunctional3D<T,nDim>::clone() const
{
    return new BoxStrainRateFunctional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void BoxStrainRateFunctional3D<T,nDim>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, int nDim>
BlockDomain::DomainT BoxStrainRateFunctional3D<T,nDim>::appliesTo() const {
    // Don't apply to envelope, because nearest neighbors need to be accessed.
    return BlockDomain::bulk;
}


/* ******** BoxBulkDivergenceFunctional3D ****************************************** */

template<typename T, int nDim>
void BoxBulkDivergenceFunctional3D<T,nDim>::process (
        Box3D domain, ScalarField3D<T>& divergence,
                      TensorField3D<T,nDim>& velocity )
{
    Dot3D offset = computeRelativeDisplacement(divergence, velocity);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                plint iZ2 = iZ+offset.z;
                divergence.get(iX,iY,iZ) =
                    fdDataField::bulkXderiv(velocity, iX2,iY2,iZ2, 0) +
                    fdDataField::bulkYderiv(velocity, iX2,iY2,iZ2, 1) +
                    fdDataField::bulkZderiv(velocity, iX2,iY2,iZ2, 2);
            }
        }
    }
}

template<typename T, int nDim>
BoxBulkDivergenceFunctional3D<T,nDim>* BoxBulkDivergenceFunctional3D<T,nDim>::clone() const
{
    return new BoxBulkDivergenceFunctional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void BoxBulkDivergenceFunctional3D<T,nDim>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;   // Divergence.
    modified[1] = modif::nothing;  // Velocity.
}

template<typename T, int nDim>
BlockDomain::DomainT BoxBulkDivergenceFunctional3D<T,nDim>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Tensor_A_plus_B_functional3D ****************************************** */

template<typename T, int nDim>
void Tensor_A_plus_B_functional3D<T,nDim>::process (
        Box3D domain, std::vector<TensorField3D<T,nDim>*> fields )
{
    PLB_PRECONDITION( fields.size()==3 );
    TensorField3D<T,nDim>& A = *fields[0];
    TensorField3D<T,nDim>& B = *fields[1];
    TensorField3D<T,nDim>& result = *fields[2];
    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                    = A.get(iX,iY,iZ) + B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z);
            }
        }
    }
}

template<typename T, int nDim>
Tensor_A_plus_B_functional3D<T,nDim>* Tensor_A_plus_B_functional3D<T,nDim>::clone() const {
    return new Tensor_A_plus_B_functional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void Tensor_A_plus_B_functional3D<T,nDim>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template<typename T, int nDim>
BlockDomain::DomainT Tensor_A_plus_B_functional3D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Tensor_A_minus_B_functional3D ************************************ */

template<typename T, int nDim>
void Tensor_A_minus_B_functional3D<T,nDim>::process (
        Box3D domain, std::vector<TensorField3D<T,nDim>*> fields )
{
    PLB_PRECONDITION( fields.size()==3 );
    TensorField3D<T,nDim>& A = *fields[0];
    TensorField3D<T,nDim>& B = *fields[1];
    TensorField3D<T,nDim>& result = *fields[2];
    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                    = A.get(iX,iY,iZ) - B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z);
            }
        }
    }
}

template<typename T, int nDim>
Tensor_A_minus_B_functional3D<T,nDim>* Tensor_A_minus_B_functional3D<T,nDim>::clone() const {
    return new Tensor_A_minus_B_functional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void Tensor_A_minus_B_functional3D<T,nDim>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template<typename T, int nDim>
BlockDomain::DomainT Tensor_A_minus_B_functional3D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Tensor_A_times_B_functional3D ************************************ */

template<typename T, int nDim>
void Tensor_A_times_B_functional3D<T,nDim>::process (
        Box3D domain, std::vector<TensorField3D<T,nDim>*> fields )
{
    PLB_PRECONDITION( fields.size()==3 );
    TensorField3D<T,nDim>& A = *fields[0];
    TensorField3D<T,nDim>& B = *fields[1];
    TensorField3D<T,nDim>& result = *fields[2];
    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                    = A.get(iX,iY,iZ) * B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z);
            }
        }
    }
}

template<typename T, int nDim>
Tensor_A_times_B_functional3D<T,nDim>* Tensor_A_times_B_functional3D<T,nDim>::clone() const {
    return new Tensor_A_times_B_functional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void Tensor_A_times_B_functional3D<T,nDim>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template<typename T, int nDim>
BlockDomain::DomainT Tensor_A_times_B_functional3D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

/* ******** IndexContraction_SymmetricTensor_A_SymmetricTensor_B_functional3D ************************************ */

template<typename T, int nDim>
void IndexContraction_SymmetricTensor_A_SymmetricTensor_B_functional3D<T,nDim>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> fields )
{
    PLB_PRECONDITION( fields.size()==3 );
    TensorField3D<T,SymmetricTensorImpl<T,nDim>::n>& A = *dynamic_cast<TensorField3D<T,SymmetricTensorImpl<T,nDim>::n>*>(fields[0]);
    TensorField3D<T,SymmetricTensorImpl<T,nDim>::n>& B = *dynamic_cast<TensorField3D<T,SymmetricTensorImpl<T,nDim>::n>*>(fields[1]);
    ScalarField3D<T>& result = *dynamic_cast<ScalarField3D<T>*>(fields[2]);

    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z) =
                        SymmetricTensorImpl<T,nDim>::contractIndexes(A.get(iX,iY,iZ),B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z));
            }
        }
    }
}

template<typename T, int nDim>
IndexContraction_SymmetricTensor_A_SymmetricTensor_B_functional3D<T,nDim>* IndexContraction_SymmetricTensor_A_SymmetricTensor_B_functional3D<T,nDim>::clone() const {
    return new IndexContraction_SymmetricTensor_A_SymmetricTensor_B_functional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void IndexContraction_SymmetricTensor_A_SymmetricTensor_B_functional3D<T,nDim>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template<typename T, int nDim>
BlockDomain::DomainT IndexContraction_SymmetricTensor_A_SymmetricTensor_B_functional3D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

/* ******** TensorProduct_A_A_functional3D ************************************ */

template<typename T, int nDim>
void TensorProduct_A_A_functional3D<T,nDim>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> fields )
{
    PLB_PRECONDITION( fields.size()==2 );
    TensorField3D<T,nDim>& A = *dynamic_cast<TensorField3D<T,nDim>*>(fields[0]);
    TensorField3D<T,SymmetricTensorImpl<T,nDim>::n>& result = *dynamic_cast<TensorField3D<T,SymmetricTensorImpl<T,nDim>::n>*>(fields[1]);

    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iPi = 0;
                for (plint iA = 0; iA < nDim; ++iA) {
                    for (plint iB = iA; iB < nDim; ++iB) {
                        result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)[iPi]
                            = A.get(iX,iY,iZ)[iA] * A.get(iX,iY,iZ)[iB];
                        ++iPi;
                    }
                }
            }
        }
    }
}

template<typename T, int nDim>
TensorProduct_A_A_functional3D<T,nDim>* TensorProduct_A_A_functional3D<T,nDim>::clone() const {
    return new TensorProduct_A_A_functional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void TensorProduct_A_A_functional3D<T,nDim>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, int nDim>
BlockDomain::DomainT TensorProduct_A_A_functional3D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

/* ******** Scalar_A_times_Tensor_B_functional3D ************************************ */

template<typename T, int nDim>
void Scalar_A_times_Tensor_B_functional3D<T,nDim>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> fields )
{
    PLB_PRECONDITION( fields.size()==3 );
    ScalarField3D<T>& A = *dynamic_cast<ScalarField3D<T>*>(fields[0]);
    TensorField3D<T,nDim>& B = *dynamic_cast<TensorField3D<T,nDim>*>(fields[1]);
    TensorField3D<T,nDim>& result = *dynamic_cast<TensorField3D<T,nDim>*>(fields[2]);
    
    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                    = A.get(iX,iY,iZ) * B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z);
            }
        }
    }
}

template<typename T, int nDim>
Scalar_A_times_Tensor_B_functional3D<T,nDim>* Scalar_A_times_Tensor_B_functional3D<T,nDim>::clone() const {
    return new Scalar_A_times_Tensor_B_functional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void Scalar_A_times_Tensor_B_functional3D<T,nDim>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template<typename T, int nDim>
BlockDomain::DomainT Scalar_A_times_Tensor_B_functional3D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

/* ******** Tensor_A_dividedBy_B_functional3D ************************************ */

template<typename T, int nDim>
void Tensor_A_dividedBy_B_functional3D<T,nDim>::process (
        Box3D domain, std::vector<TensorField3D<T,nDim>*> fields )
{
    PLB_PRECONDITION( fields.size()==3 );
    TensorField3D<T,nDim>& A = *fields[0];
    TensorField3D<T,nDim>& B = *fields[1];
    TensorField3D<T,nDim>& result = *fields[2];
    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                    = A.get(iX,iY,iZ) / B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z);
            }
        }
    }
}

template<typename T, int nDim>
Tensor_A_dividedBy_B_functional3D<T,nDim>* Tensor_A_dividedBy_B_functional3D<T,nDim>::clone() const {
    return new Tensor_A_dividedBy_B_functional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void Tensor_A_dividedBy_B_functional3D<T,nDim>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template<typename T, int nDim>
BlockDomain::DomainT Tensor_A_dividedBy_B_functional3D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Tensor_A_plus_B_inplace_functional3D ************************************ */

template<typename T, int nDim>
void Tensor_A_plus_B_inplace_functional3D<T,nDim>::process (
        Box3D domain, TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B)
{
    Dot3D offset = computeRelativeDisplacement(A,B);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                A.get(iX,iY,iZ) += B.get(iX+offset.x,iY+offset.y,iZ+offset.z);
            }
        }
    }
}

template<typename T, int nDim>
Tensor_A_plus_B_inplace_functional3D<T,nDim>* Tensor_A_plus_B_inplace_functional3D<T,nDim>::clone() const {
    return new Tensor_A_plus_B_inplace_functional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void Tensor_A_plus_B_inplace_functional3D<T,nDim>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T, int nDim>
BlockDomain::DomainT Tensor_A_plus_B_inplace_functional3D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Tensor_A_minus_B_inplace_functional3D ************************************ */

template<typename T, int nDim>
void Tensor_A_minus_B_inplace_functional3D<T,nDim>::process (
        Box3D domain, TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B)
{
    Dot3D offset = computeRelativeDisplacement(A,B);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                A.get(iX,iY,iZ) -= B.get(iX+offset.x,iY+offset.y,iZ+offset.z);
            }
        }
    }
}

template<typename T, int nDim>
Tensor_A_minus_B_inplace_functional3D<T,nDim>* Tensor_A_minus_B_inplace_functional3D<T,nDim>::clone() const {
    return new Tensor_A_minus_B_inplace_functional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void Tensor_A_minus_B_inplace_functional3D<T,nDim>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T, int nDim>
BlockDomain::DomainT Tensor_A_minus_B_inplace_functional3D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Tensor_A_times_B_inplace_functional3D ************************************ */

template<typename T, int nDim>
void Tensor_A_times_B_inplace_functional3D<T,nDim>::process (
        Box3D domain, TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B)
{
    Dot3D offset = computeRelativeDisplacement(A,B);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                A.get(iX,iY,iZ) *= B.get(iX+offset.x,iY+offset.y,iZ+offset.z);
            }
        }
    }
}

template<typename T, int nDim>
Tensor_A_times_B_inplace_functional3D<T,nDim>* Tensor_A_times_B_inplace_functional3D<T,nDim>::clone() const {
    return new Tensor_A_times_B_inplace_functional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void Tensor_A_times_B_inplace_functional3D<T,nDim>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T, int nDim>
BlockDomain::DomainT Tensor_A_times_B_inplace_functional3D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Tensor_A_times_alpha_inplace_functional3D ************************************ */

template<typename T, int nDim>
Tensor_A_times_alpha_inplace_functional3D<T,nDim>::
    Tensor_A_times_alpha_inplace_functional3D(T alpha_)
        : alpha(alpha_)
{ }

template<typename T, int nDim>
void Tensor_A_times_alpha_inplace_functional3D<T,nDim>::process (
        Box3D domain, TensorField3D<T,nDim>& A)
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                A.get(iX,iY,iZ) *= alpha;
            }
        }
    }
}

template<typename T, int nDim>
Tensor_A_times_alpha_inplace_functional3D<T,nDim>*
    Tensor_A_times_alpha_inplace_functional3D<T,nDim>::clone() const
{
    return new Tensor_A_times_alpha_inplace_functional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void Tensor_A_times_alpha_inplace_functional3D<T,nDim>::
    getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::staticVariables;
}

template<typename T, int nDim>
BlockDomain::DomainT Tensor_A_times_alpha_inplace_functional3D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

/* ******** Tensor_A_times_alpha_functional3D ************************************ */

template<typename T, int nDim>
Tensor_A_times_alpha_functional3D<T,nDim>::
    Tensor_A_times_alpha_functional3D(T alpha_)
        : alpha(alpha_)
{ }

template<typename T, int nDim>
void Tensor_A_times_alpha_functional3D<T,nDim>::process (
    Box3D domain, TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& result)
{
    Dot3D offset = computeRelativeDisplacement(A, result);
    
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        plint bX = iX + offset.x;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            plint bY = iY + offset.y;
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint bZ = iZ + offset.z;
                result.get(bX,bY,bZ) = alpha*A.get(iX,iY,iZ);
            }
        }
    }
}

template<typename T, int nDim>
Tensor_A_times_alpha_functional3D<T,nDim>*
    Tensor_A_times_alpha_functional3D<T,nDim>::clone() const
{
    return new Tensor_A_times_alpha_functional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void Tensor_A_times_alpha_functional3D<T,nDim>::
    getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T, int nDim>
BlockDomain::DomainT Tensor_A_times_alpha_functional3D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Tensor_A_dividedBy_B_inplace_functional3D ************************************ */

template<typename T, int nDim>
void Tensor_A_dividedBy_B_inplace_functional3D<T,nDim>::process (
        Box3D domain, TensorField3D<T,nDim>& A, TensorField3D<T,nDim>& B)
{
    Dot3D offset = computeRelativeDisplacement(A,B);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                A.get(iX,iY,iZ) /= B.get(iX+offset.x,iY+offset.y,iZ+offset.z);
            }
        }
    }
}

template<typename T, int nDim>
Tensor_A_dividedBy_B_inplace_functional3D<T,nDim>* Tensor_A_dividedBy_B_inplace_functional3D<T,nDim>::clone() const {
    return new Tensor_A_dividedBy_B_inplace_functional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void Tensor_A_dividedBy_B_inplace_functional3D<T,nDim>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T, int nDim>
BlockDomain::DomainT Tensor_A_dividedBy_B_inplace_functional3D<T,nDim>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Normalize_Tensor_functional3D ****************************************** */

template<typename T, int nDim>
void Normalize_Tensor_functional3D<T,nDim>::process (
        Box3D domain, TensorField3D<T,nDim>& data, TensorField3D<T,nDim>& result )
{
    T eps = getEpsilon<T>(precision);
    Dot3D offset = computeRelativeDisplacement(data, result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Array<T,nDim>& a = data.get(iX,iY,iZ);
                T norm = T();
                for (int i = 0; i < nDim; norm += a[i]*a[i], i++)
                    ;
                norm = sqrt(norm);
                if (norm <= eps)
                    result.get(iX+offset.x,iY+offset.y,iZ+offset.z).resetToZero();
                else
                    result.get(iX+offset.x,iY+offset.y,iZ+offset.z) = a / norm;
            }
        }
    }
}

template<typename T, int nDim>
Normalize_Tensor_functional3D<T,nDim>* Normalize_Tensor_functional3D<T,nDim>::clone() const {
    return new Normalize_Tensor_functional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void Normalize_Tensor_functional3D<T,nDim>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}


/* ******** Normalize_Tensor_inplace_functional3D ************************************ */

template<typename T, int nDim>
void Normalize_Tensor_inplace_functional3D<T,nDim>::process (
        Box3D domain, TensorField3D<T,nDim>& data)
{
    T eps = getEpsilon<T>(precision);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Array<T,nDim>& a = data.get(iX,iY,iZ);
                T norm = T();
                for (int i = 0; i < nDim; norm += a[i]*a[i], i++)
                    ;
                norm = sqrt(norm);
                if (norm <= eps)
                    a.resetToZero();
                else
                    a /= norm;
            }
        }
    }
}

template<typename T, int nDim>
Normalize_Tensor_inplace_functional3D<T,nDim>*
    Normalize_Tensor_inplace_functional3D<T,nDim>::clone() const
{
    return new Normalize_Tensor_inplace_functional3D<T,nDim>(*this);
}

template<typename T, int nDim>
void Normalize_Tensor_inplace_functional3D<T,nDim>::
    getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::staticVariables;
}


/* ************* Class LBMsmoothen3D ******************* */

template<typename T, template<typename U> class Descriptor>
void LBMsmoothen3D<T,Descriptor>::process (
        Box3D domain, ScalarField3D<T>& data, ScalarField3D<T>& result )
{
    typedef Descriptor<T> D;
    Dot3D offset = computeRelativeDisplacement(data, result);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                result.get(iX+offset.x,iY+offset.y,iZ+offset.z) = (T)0;
                T sum = (T) 0;
                for (plint iPop=1; iPop<D::q; ++iPop) {
                    plint nextX = iX+D::c[iPop][0];
                    plint nextY = iY+D::c[iPop][1];
                    plint nextZ = iZ+D::c[iPop][2];
                    sum += D::t[iPop];
                    result.get(iX+offset.x,iY+offset.y,iZ+offset.z) +=
                        D::t[iPop] * data.get(nextX,nextY,nextZ);
                }
                result.get(iX+offset.x,iY+offset.y,iZ+offset.z) /= sum;
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
LBMsmoothen3D<T,Descriptor>*
    LBMsmoothen3D<T,Descriptor>::clone() const
{
    return new LBMsmoothen3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void LBMsmoothen3D<T,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}


/* ************* Class Smoothen3D ******************* */

template<typename T>
void Smoothen3D<T>::process(Box3D domain, ScalarField3D<T>& data, ScalarField3D<T>& result)
{
    Dot3D offset = computeRelativeDisplacement(data, result);

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                T *res = &result.get(iX+offset.x, iY+offset.y, iZ+offset.z);
                *res = (T) 0;
                int n = 0;
                for (int i = -1; i < 2; i++) {
                    plint nextX = iX + i;
                    for (int j = -1; j < 2; j++) {
                        plint nextY = iY + j;
                        for (int k = -1; k < 2; k++) {
                            plint nextZ = iZ + k;
                            if (!(i == 0 && j == 0 && k == 0)) {
                                n++;
                                *res += data.get(nextX, nextY, nextZ);
                            }
                        }
                    }
                }
                *res /= (T) n;
            }
        }
    }
}

template<typename T>
Smoothen3D<T>* Smoothen3D<T>::clone() const
{
    return new Smoothen3D<T>(*this);
}

template<typename T>
void Smoothen3D<T>::getTypeOfModification (std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;         // Data.
    modified[1] = modif::nothing;         // Flags.
    modified[2] = modif::staticVariables; // Result.
}


/* ************* Class LBMsmoothenTensor3D ******************* */

template<typename T, int nDim, template<typename U> class Descriptor>
void LBMsmoothenTensor3D<T,nDim,Descriptor>::process (
        Box3D domain, TensorField3D<T,nDim>& data, TensorField3D<T,nDim>& result )
{
    typedef Descriptor<T> D;
    Dot3D offset = computeRelativeDisplacement(data, result);

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                result.get(iX+offset.x,iY+offset.y,iZ+offset.z).resetToZero();
                T sum = (T) 0;
                for (plint iPop=1; iPop<D::q; ++iPop) {
                    plint nextX = iX+D::c[iPop][0];
                    plint nextY = iY+D::c[iPop][1];
                    plint nextZ = iZ+D::c[iPop][2];
                    sum += D::t[iPop];
                    result.get(iX+offset.x,iY+offset.y,iZ+offset.z) +=
                        D::t[iPop] * data.get(nextX,nextY,nextZ);
                }
                result.get(iX+offset.x,iY+offset.y,iZ+offset.z) /= sum;
            }
        }
    }
}

template<typename T, int nDim, template<typename U> class Descriptor>
LBMsmoothenTensor3D<T,nDim,Descriptor>*
    LBMsmoothenTensor3D<T,nDim,Descriptor>::clone() const
{
    return new LBMsmoothenTensor3D<T,nDim,Descriptor>(*this);
}

template<typename T, int nDim, template<typename U> class Descriptor>
void LBMsmoothenTensor3D<T,nDim,Descriptor>::getTypeOfModification (
        std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}


/* ************* Class SmoothenTensor3D ******************* */

template<typename T, int nDim>
void SmoothenTensor3D<T,nDim>::process (
        Box3D domain, TensorField3D<T,nDim>& data, TensorField3D<T,nDim>& result )
{
    Dot3D offset = computeRelativeDisplacement(data, result);

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Array<T,3> *res = &result.get(iX+offset.x, iY+offset.y, iZ+offset.z);
                res->resetToZero();
                int n = 0;
                for (int i = -1; i < 2; i++) {
                    plint nextX = iX + i;
                    for (int j = -1; j < 2; j++) {
                        plint nextY = iY + j;
                        for (int k = -1; k < 2; k++) {
                            plint nextZ = iZ + k;
                            if (!(i == 0 && j == 0 && k == 0)) {
                                n++;
                                *res += data.get(nextX, nextY, nextZ);
                            }
                        }
                    }
                }
                *res /= (T) n;
            }
        }
    }
}

template<typename T, int nDim>
SmoothenTensor3D<T,nDim>* SmoothenTensor3D<T,nDim>::clone() const
{
    return new SmoothenTensor3D<T,nDim>(*this);
}

template<typename T, int nDim>
void SmoothenTensor3D<T,nDim>::getTypeOfModification (
        std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;         // Data.
    modified[1] = modif::nothing;         // Flags.
    modified[2] = modif::staticVariables; // Result.
}


/* ************* Class MollifyScalar3D ******************* */

template<typename T>
MollifyScalar3D<T>::MollifyScalar3D(T l_, plint d_, Box3D globalDomain_, int exclusionFlag_)
    : l(l_),
      d(d_),
      globalDomain(globalDomain_),
      exclusionFlag(exclusionFlag_)
{ }

template<typename T>
void MollifyScalar3D<T>::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields)
{
    PLB_ASSERT(fields.size() == 3);
    ScalarField3D<T> *scalar = dynamic_cast<ScalarField3D<T>*>(fields[0]);
    PLB_ASSERT(scalar);
    ScalarField3D<int> *flag = dynamic_cast<ScalarField3D<int>*>(fields[1]);
    PLB_ASSERT(flag);
    ScalarField3D<T> *result = dynamic_cast<ScalarField3D<T>*>(fields[2]);
    PLB_ASSERT(result);

    Dot3D absOfs = scalar->getLocation();

    Dot3D ofsF = computeRelativeDisplacement(*scalar, *flag);
    Dot3D ofsR = computeRelativeDisplacement(*scalar, *result);

    static T pi = acos((T) -1);

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (flag->get(iX + ofsF.x, iY + ofsF.y, iZ + ofsF.z) == exclusionFlag) {
                    result->get(iX + ofsR.x, iY + ofsR.y, iZ + ofsR.z) = scalar->get(iX, iY, iZ);
                    continue;
                }

                T sum = 0.0;
                for (plint dx = -d; dx <= d; dx++) {
                    plint i = iX + dx;
                    for (plint dy = -d; dy <= d; dy++) {
                        plint j = iY + dy;
                        for (plint dz = -d; dz <= d; dz++) {
                            plint k = iZ + dz;
                            if (contained(i + absOfs.x, j + absOfs.y, k + absOfs.z, globalDomain)) {
                                if (flag->get(i + ofsF.x, j + ofsF.y, k + ofsF.z) != exclusionFlag) {
                                    T r = sqrt(dx * dx + dy * dy + dz * dz);
                                    T integrand = scalar->get(i, j, k);
                                    T mollifier = 0.0;
                                    if (r < l) {
                                        mollifier = (1.0 + cos(pi * r / l));
                                    }
                                    sum += integrand * mollifier;
                                }
                            }
                        }
                    }
                }
                sum *= (1.0 / (2.0 * l));
                result->get(iX + ofsR.x, iY + ofsR.y, iZ + ofsR.z) = sum;
            }
        }
    }
}

template<typename T>
MollifyScalar3D<T>* MollifyScalar3D<T>::clone() const
{
    return new MollifyScalar3D<T>(*this);
}

template<typename T>
void MollifyScalar3D<T>::getTypeOfModification (std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;         // Data.
    modified[1] = modif::nothing;         // Flags.
    modified[2] = modif::staticVariables; // Result.
}


/* ************* Class MollifyTensor3D ******************* */

template<typename T, int nDim>
MollifyTensor3D<T,nDim>::MollifyTensor3D(T l_, plint d_, Box3D globalDomain_, int exclusionFlag_)
    : l(l_),
      d(d_),
      globalDomain(globalDomain_),
      exclusionFlag(exclusionFlag_)
{ }

template<typename T, int nDim>
void MollifyTensor3D<T,nDim>::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields)
{
    PLB_ASSERT(fields.size() == 3);
    TensorField3D<T,nDim> *tensor = dynamic_cast<TensorField3D<T,nDim>*>(fields[0]);
    PLB_ASSERT(tensor);
    ScalarField3D<int> *flag = dynamic_cast<ScalarField3D<int>*>(fields[1]);
    PLB_ASSERT(flag);
    TensorField3D<T,nDim> *result = dynamic_cast<TensorField3D<T,nDim>*>(fields[2]);
    PLB_ASSERT(result);

    Dot3D absOfs = tensor->getLocation();

    Dot3D ofsF = computeRelativeDisplacement(*tensor, *flag);
    Dot3D ofsR = computeRelativeDisplacement(*tensor, *result);

    static T pi = acos((T) -1);

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (flag->get(iX + ofsF.x, iY + ofsF.y, iZ + ofsF.z) == exclusionFlag) {
                    result->get(iX + ofsR.x, iY + ofsR.y, iZ + ofsR.z) = tensor->get(iX, iY, iZ);
                    continue;
                }

                Array<T,nDim> sum;
                sum.resetToZero();
                for (plint dx = -d; dx <= d; dx++) {
                    plint i = iX + dx;
                    for (plint dy = -d; dy <= d; dy++) {
                        plint j = iY + dy;
                        for (plint dz = -d; dz <= d; dz++) {
                            plint k = iZ + dz;
                            if (contained(i + absOfs.x, j + absOfs.y, k + absOfs.z, globalDomain)) {
                                if (flag->get(i + ofsF.x, j + ofsF.y, k + ofsF.z) != exclusionFlag) {
                                    T r = sqrt(dx * dx + dy * dy + dz * dz);
                                    Array<T,nDim> integrand = tensor->get(i, j, k);
                                    T mollifier = 0.0;
                                    if (r < l) {
                                        mollifier = (1.0 + cos(pi * r / l));
                                    }
                                    sum += integrand * mollifier;
                                }
                            }
                        }
                    }
                }
                sum *= (1.0 / (2.0 * l));
                result->get(iX + ofsR.x, iY + ofsR.y, iZ + ofsR.z) = sum;
            }
        }
    }
}

template<typename T, int nDim>
MollifyTensor3D<T,nDim>* MollifyTensor3D<T,nDim>::clone() const
{
    return new MollifyTensor3D<T,nDim>(*this);
}

template<typename T, int nDim>
void MollifyTensor3D<T,nDim>::getTypeOfModification (std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;         // Data.
    modified[1] = modif::nothing;         // Flags.
    modified[2] = modif::staticVariables; // Result.
}


/* ************* Class LBMcomputeGradient3D ******************* */

template<typename T,template<typename U> class Descriptor>
void LBMcomputeGradient3D<T,Descriptor>::process (
        Box3D domain, ScalarField3D<T>& scalarField, TensorField3D<T,3>& gradient )
{
    typedef Descriptor<T> D;
    Dot3D ofs = computeRelativeDisplacement(scalarField, gradient);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Array<T,3>& gradientOfScalar = gradient.get(iX+ofs.x,iY+ofs.y,iZ+ofs.z);                    
                gradientOfScalar.resetToZero();
                // Compute the gradient of a scalar function "the lattice Boltzmann way".
                for (plint iPop=1; iPop < D::q; ++iPop ) {
                    plint nextX = iX + D::c[iPop][0];
                    plint nextY = iY + D::c[iPop][1];
                    plint nextZ = iZ + D::c[iPop][2];

                    gradientOfScalar[0] += D::t[iPop]*D::c[iPop][0]*scalarField.get(nextX,nextY,nextZ);
                    gradientOfScalar[1] += D::t[iPop]*D::c[iPop][1]*scalarField.get(nextX,nextY,nextZ);
                    gradientOfScalar[2] += D::t[iPop]*D::c[iPop][2]*scalarField.get(nextX,nextY,nextZ);
                }
                gradientOfScalar *= D::invCs2;                     
            }
        }
    }
}


/* ************* Class LBMcomputeDivergence3D ******************* */

template<typename T,template<typename U> class Descriptor>
void LBMcomputeDivergence3D<T,Descriptor>::process (
        Box3D domain, ScalarField3D<T>& divergence, TensorField3D<T,3>& vectorField )
{
    typedef Descriptor<T> D;
    Dot3D ofs = computeRelativeDisplacement(divergence, vectorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                // Compute the divergence of a vector function "the lattice Boltzmann way".
                T sum = T();
                for (plint iPop=1; iPop < D::q; ++iPop ) {
                    plint nextX = iX + D::c[iPop][0];
                    plint nextY = iY + D::c[iPop][1];
                    plint nextZ = iZ + D::c[iPop][2];

                    Array<T,3>& vectorFunction = vectorField.get(nextX+ofs.x,nextY+ofs.y,nextZ+ofs.z);                    

                    T tmp = D::c[iPop][0]*vectorFunction[0] +
                            D::c[iPop][1]*vectorFunction[1] +
                            D::c[iPop][2]*vectorFunction[2];

                    sum += D::t[iPop]*tmp;
                }
                divergence.get(iX,iY,iZ) = D::invCs2 * sum;                     
            }
        }
    }
}

}  // namespace plb

#endif  // DATA_ANALYSIS_FUNCTIONAL_3D_HH
