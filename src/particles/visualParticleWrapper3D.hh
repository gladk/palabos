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

#ifndef VISUAL_PARTICLE_WRAPPER_3D_HH
#define VISUAL_PARTICLE_WRAPPER_3D_HH

#include "core/globalDefs.h"
#include "core/array.h"
#include "core/geometry3D.h"
#include "particles/visualParticleWrapper3D.h"
#include "particles/visualParticleFunctional3D.h"
#include "particles/visualParticle3D.h"

#include <vector>
#include <string>
#include <cmath>

namespace plb {

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiParticleField3D<DenseParticleField3D<T,Descriptor> > >
    computeSurfaceForce (
            TriangleBoundary3D<T>& boundary,
            VoxelizedDomain3D<T>& voxelizedDomain,
            MultiBlockLattice3D<T,Descriptor>& lattice,
            bool incompressibleModel, bool dynamicMesh )
{
    return computeSurfaceForce(boundary, voxelizedDomain, lattice,
                               lattice.getBoundingBox(), incompressibleModel, dynamicMesh);
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiParticleField3D<DenseParticleField3D<T,Descriptor> > >
    computeSurfaceForce (
            TriangleBoundary3D<T>& boundary,
            VoxelizedDomain3D<T>& voxelizedDomain,
            MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain,
            bool incompressibleModel, bool dynamicMesh )
{
    if (dynamicMesh) {
        // Select Dynamic, Open mesh.
        boundary.pushSelect(0,1);
    }
    else {
        // Select Static, Open mesh.
        boundary.pushSelect(0,0);
    }
    std::auto_ptr<MultiParticleField3D<DenseParticleField3D<T,Descriptor> > >
        particles( new MultiParticleField3D<DenseParticleField3D<T,Descriptor> > (
                       lattice.getMultiBlockManagement(),
                       defaultMultiBlockPolicy3D().getCombinedStatistics() ) );

    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(particles.get());
    applyProcessingFunctional ( 
            new CreateParticleFromVertex3D<T,Descriptor,VisualParticle3D<T,Descriptor> > (
                    boundary.getMesh() ),
            domain, particleArg );

    std::vector<MultiBlock3D*> particleFluidVoxelArg;
    particleFluidVoxelArg.push_back(particles.get());
    particleFluidVoxelArg.push_back(&lattice);
    particleFluidVoxelArg.push_back(&voxelizedDomain.getVoxelMatrix());
    applyProcessingFunctional (
            new ComputeParticleForce3D<T,Descriptor> (
                    boundary, voxelizedDomain.getFlowType(), incompressibleModel ),
            domain, particleFluidVoxelArg);
    // Restore original setting.
    boundary.popSelect();
    return particles;
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiParticleField3D<DenseParticleField3D<T,Descriptor> > >
    computeSurfaceForce (
            TriangleBoundary3D<T>& boundary,
            MultiParticleField3D<DenseParticleField3D<T,Descriptor> >& parallelWall,
            VoxelizedDomain3D<T>& voxelizedDomain,
            MultiBlockLattice3D<T,Descriptor>& lattice,
            bool incompressibleModel, bool dynamicMesh )
{
    if (dynamicMesh) {
        // Select Dynamic, Open mesh.
        boundary.pushSelect(0,1);
    }
    else {
        // Select Static, Open mesh.
        boundary.pushSelect(0,0);
    }
    std::auto_ptr<MultiParticleField3D<DenseParticleField3D<T,Descriptor> > >
        particles( new MultiParticleField3D<DenseParticleField3D<T,Descriptor> > (
                       lattice.getMultiBlockManagement(),
                       defaultMultiBlockPolicy3D().getCombinedStatistics() ) );

    std::vector<MultiBlock3D*> particleParticleArg;
    particleParticleArg.push_back(&parallelWall);
    particleParticleArg.push_back(particles.get());
    applyProcessingFunctional ( 
            new CloneParticles3D<T,Descriptor>(new VisualParticle3D<T,Descriptor>),
            lattice.getBoundingBox(), particleParticleArg );

    std::vector<MultiBlock3D*> particleFluidVoxelArg;
    particleFluidVoxelArg.push_back(particles.get());
    particleFluidVoxelArg.push_back(&lattice);
    particleFluidVoxelArg.push_back(&voxelizedDomain.getVoxelMatrix());
    applyProcessingFunctional (
            new ComputeParticleForce3D<T,Descriptor> (
                    boundary, voxelizedDomain.getFlowType(), incompressibleModel ),
            lattice.getBoundingBox(), particleFluidVoxelArg);
    // Restore original setting.
    boundary.popSelect();
    return particles;
}


template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiParticleField3D<DenseParticleField3D<T,Descriptor> > >
    computeSurfaceScalar (
            TriangleBoundary3D<T>& boundary,
            VoxelizedDomain3D<T>& voxelizedDomain,
            MultiScalarField3D<T>& scalarField,
            const std::vector<int> &usableFlowTypes,
            plint whichScalar, bool dynamicMesh )
{
    return computeSurfaceScalar<T,Descriptor>(boundary, voxelizedDomain, scalarField,
                               scalarField.getBoundingBox(), usableFlowTypes, whichScalar, dynamicMesh);
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiParticleField3D<DenseParticleField3D<T,Descriptor> > >
    computeSurfaceScalar (
            TriangleBoundary3D<T>& boundary,
            VoxelizedDomain3D<T>& voxelizedDomain,
            MultiScalarField3D<T>& scalarField, Box3D domain,
            const std::vector<int> &usableFlowTypes,
            plint whichScalar, bool dynamicMesh )
{
    if (dynamicMesh) {
        // Select Dynamic, Open mesh.
        boundary.pushSelect(0,1);
    }
    else {
        // Select Static, Open mesh.
        boundary.pushSelect(0,0);
    }

    std::auto_ptr<MultiScalarField3D<T> > extendedScalarField 
        = extendEnvelopeWidth(scalarField, 2);

    std::auto_ptr<MultiParticleField3D<DenseParticleField3D<T,Descriptor> > >
        particles( new MultiParticleField3D<DenseParticleField3D<T,Descriptor> > (
                       extendedScalarField->getMultiBlockManagement(),
                       defaultMultiBlockPolicy3D().getCombinedStatistics() ) );

    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(particles.get());
    applyProcessingFunctional ( 
            new CreateParticleFromVertex3D<T,Descriptor,VisualParticle3D<T,Descriptor> > (
                    boundary.getMesh() ),
            domain, particleArg );


    std::vector<MultiBlock3D*> particleFluidVoxelArg;
    particleFluidVoxelArg.push_back(particles.get());
    particleFluidVoxelArg.push_back(extendedScalarField.get());
    particleFluidVoxelArg.push_back(&voxelizedDomain.getVoxelMatrix());
    applyProcessingFunctional (
            new LinearScalarFieldToParticle3D<T,Descriptor> (
                    usableFlowTypes ),
            domain, particleFluidVoxelArg);
    // Restore original setting.
    boundary.popSelect();
    return particles;
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiParticleField3D<DenseParticleField3D<T,Descriptor> > >
    getParticlePosAndVelocity (
            MultiParticleField3D<DenseParticleField3D<T,Descriptor> >& originalParticles )
{
    std::auto_ptr<MultiParticleField3D<DenseParticleField3D<T,Descriptor> > >
        particles( new MultiParticleField3D<DenseParticleField3D<T,Descriptor> > (
                       originalParticles.getMultiBlockManagement(),
                       defaultMultiBlockPolicy3D().getCombinedStatistics() ) );

    std::vector<MultiBlock3D*> particleParticleArg;
    particleParticleArg.push_back(&originalParticles);
    particleParticleArg.push_back(particles.get());
    applyProcessingFunctional ( 
            new GetParticleVelocity3D<T,Descriptor>(), particles->getBoundingBox(), particleParticleArg );

    return particles;
}

template<typename T, template<typename U> class Descriptor>
void injectVisualScalarFieldParticles (
        MultiScalarField3D<T>& scalarField,
        MultiParticleField3D<DenseParticleField3D<T,Descriptor> >& particleField,
        Box3D domain, Particle3D<T,Descriptor>* particleTemplate,
        T mostUnlikely, T mostLikely,
        T probability, plint numShotsPerCell )
{
    plint ny = scalarField.getNy();
    plint nz = scalarField.getNz();

    std::vector<MultiBlock3D*> particleScalarArg;
    particleScalarArg.push_back(&particleField);
    particleScalarArg.push_back(&scalarField);

    applyProcessingFunctional (
            new VisualScalarFieldParticles3D<T,Descriptor> (
                mostUnlikely, mostLikely, numShotsPerCell,
                probability, ny, nz, particleTemplate ),
            domain, particleScalarArg );
}

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiParticleField3D<DenseParticleField3D<T,Descriptor> > >
    injectVisualScalarFieldParticles (
        MultiScalarField3D<T>& scalarField,
        Box3D domain, Particle3D<T,Descriptor>* particleTemplate,
        T mostUnlikely, T mostLikely,
        T probability, plint numShotsPerCell )
{
    std::auto_ptr<MultiParticleField3D<DenseParticleField3D<T,Descriptor> > >
        particleField( new MultiParticleField3D<DenseParticleField3D<T,Descriptor> > (
                       scalarField.getMultiBlockManagement(),
                       defaultMultiBlockPolicy3D().getCombinedStatistics() ) );
    injectVisualScalarFieldParticles (
            scalarField, *particleField, domain, particleTemplate,
            mostUnlikely, mostLikely, probability, numShotsPerCell );
    return particleField;
}

template<typename T, template<typename U> class Descriptor, class ParticleFieldT>
void scalarFieldToParticles(MultiScalarField3D<T>& scalar, MultiParticleField3D<ParticleFieldT>& particleField, Box3D domain)
{
    std::vector<MultiBlock3D*> particleScalarArg;
    particleScalarArg.push_back(&particleField);
    particleScalarArg.push_back(&scalar);
    applyProcessingFunctional (
            new ScalarFieldToParticle3D<T,Descriptor>(), domain, particleScalarArg );
}

template<typename T, int nDim, template<typename U> class Descriptor, class ParticleFieldT>
void tensorFieldToParticles(MultiTensorField3D<T,nDim>& tensor, MultiParticleField3D<ParticleFieldT>& particleField, Box3D domain)
{
    std::vector<MultiBlock3D*> arg;
    arg.push_back(&particleField);
    arg.push_back(&tensor);
    applyProcessingFunctional (
            new TensorFieldToParticle3D<T,nDim,Descriptor>(), domain, arg );
}


template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiParticleField3D<DenseParticleField3D<T,Descriptor> > >
    scalarFieldToSurface(TriangularSurfaceMesh<T>& mesh, MultiScalarField3D<T>& scalar, Box3D domain)
{
    std::auto_ptr<MultiParticleField3D<DenseParticleField3D<T,Descriptor> > >
        particles( new MultiParticleField3D<DenseParticleField3D<T,Descriptor> > (
                       scalar.getMultiBlockManagement(),
                       defaultMultiBlockPolicy3D().getCombinedStatistics() ) );

    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(particles.get());
    applyProcessingFunctional ( 
      new CreateParticleFromVertex3D<T,Descriptor,VisualParticle3D<T,Descriptor> >(mesh),
      domain, particleArg );

    std::vector<MultiBlock3D*> particleScalarArg;
    particleScalarArg.push_back(particles.get());
    particleScalarArg.push_back(&scalar);
    applyProcessingFunctional (
            new ScalarFieldToParticle3D<T,Descriptor>(),
            domain, particleScalarArg);

    return particles;
}

}  // namespace plb

#endif  // VISUAL_PARTICLE_WRAPPER_3D_HH

