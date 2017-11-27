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

#ifndef OFF_LATTICE_BOUNDARY_CONDITION_3D_H
#define OFF_LATTICE_BOUNDARY_CONDITION_3D_H

#include "core/globalDefs.h"
#include "offLattice/triangularSurfaceMesh.h"
#include "offLattice/offLatticeBoundaryProfiles3D.h"
#include "offLattice/triangleBoundary3D.h"
#include "offLattice/triangleToDef.h"
#include "offLattice/guoOffLatticeModel3D.h"
#include "multiBlock/multiBlockLattice3D.h"

namespace plb {

template< typename T,
          template<typename U> class Descriptor,
          class BoundaryType >
class OffLatticeBoundaryCondition3D {
public:
    OffLatticeBoundaryCondition3D (
            OffLatticeModel3D<T,BoundaryType>* offLatticeModel_,
            VoxelizedDomain3D<T>& voxelizedDomain_,
            MultiBlockLattice3D<T,Descriptor>& lattice_ );
    OffLatticeBoundaryCondition3D (
            OffLatticeModel3D<T,BoundaryType>* offLatticeModel_,
            VoxelizedDomain3D<T>& voxelizedDomain_,
            MultiBlockLattice3D<T,Descriptor>& lattice_,
            MultiParticleField3D<DenseParticleField3D<T,Descriptor> >& particleField_ );
    OffLatticeBoundaryCondition3D (
            OffLatticeBoundaryCondition3D<T,Descriptor,BoundaryType> const& rhs );
    ~OffLatticeBoundaryCondition3D();
    MultiBlockLattice3D<T,Descriptor> const& getLattice() const { return lattice; }
    VoxelizedDomain3D<T> const& getVoxelizedDomain() const { return voxelizedDomain; }
    VoxelizedDomain3D<T>& getVoxelizedDomain() { return voxelizedDomain; }
    void apply();
    void insert(plint processorLevel = 1);
    void apply(std::vector<MultiBlock3D*> const& completionArg);
    void insert(std::vector<MultiBlock3D*> const& completionArg, plint processorLevel = 1);
    Array<T,3> getForceOnObject();
    std::auto_ptr<MultiTensorField3D<T,3> > computeVelocity(Box3D domain);
    std::auto_ptr<MultiTensorField3D<T,3> > computeVelocity();
    std::auto_ptr<MultiTensorField3D<T,3> > computeVorticity(Box3D domain);
    std::auto_ptr<MultiTensorField3D<T,3> > computeVorticity();
    std::auto_ptr<MultiScalarField3D<T> > computeVelocityNorm(Box3D domain);
    std::auto_ptr<MultiScalarField3D<T> > computeVelocityNorm();
    std::auto_ptr<MultiScalarField3D<T> > computeVorticityNorm(Box3D domain);
    std::auto_ptr<MultiScalarField3D<T> > computeVorticityNorm();
    std::auto_ptr<MultiScalarField3D<T> > computeVelocityComponent(Box3D domain, plint iComp);
    std::auto_ptr<MultiScalarField3D<T> > computeVelocityComponent(plint iComp);
    std::auto_ptr<MultiScalarField3D<T> > computePressure(Box3D domain);
    std::auto_ptr<MultiScalarField3D<T> > computePressure();
    std::auto_ptr<MultiScalarField3D<T> > computeDensity(Box3D domain, T solidDensity=T());
    std::auto_ptr<MultiScalarField3D<T> > computeDensity(T solidDensity=T());
    std::auto_ptr<MultiScalarField3D<T> > computeStrainRateNorm();
    std::auto_ptr<MultiScalarField3D<T> > computeStrainRateNorm(Box3D domain);
    std::auto_ptr<MultiTensorField3D<T,SymmetricTensor<T,Descriptor>::n> > computeStrainRate();
    std::auto_ptr<MultiTensorField3D<T,SymmetricTensor<T,Descriptor>::n> > computeStrainRate(Box3D domain);
    std::auto_ptr<MultiScalarField3D<T> > computeShearStressNorm();
    std::auto_ptr<MultiScalarField3D<T> > computeShearStressNorm(Box3D domain);
    std::auto_ptr<MultiTensorField3D<T,SymmetricTensor<T,Descriptor>::n> > computeShearStress();
    std::auto_ptr<MultiTensorField3D<T,SymmetricTensor<T,Descriptor>::n> > computeShearStress(Box3D domain);
    T computeAverageVelocityComponent(Box3D domain, plint iComponent);
    Array<T,3> computeAverageVelocity(Box3D domain);
    T computeAverageDensity(Box3D domain);
    T computeAverageDensity();
    T computeAverageEnergy(Box3D domain);
    T computeAverageEnergy();
    T computeRMSvorticity(Box3D domain);
    T computeRMSvorticity();
    T computeAverageShearStressNorm(Box3D domain);
    T computeAverageShearStressNorm();
    T computeRMSshearStressNorm(Box3D domain);
    T computeRMSshearStressNorm();
private:
    VoxelizedDomain3D<T>& voxelizedDomain;
    MultiBlockLattice3D<T,Descriptor>& lattice;
    MultiBlock3D& boundaryShapeArg;
    OffLatticeModel3D<T,BoundaryType>* offLatticeModel;
    MultiContainerBlock3D offLatticePattern;
};

}  // namespace plb

#endif  // OFF_LATTICE_BOUNDARY_CONDITION_3D_H
