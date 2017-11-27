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
  * This program solves the 3D Poiseuille flow inside a tube. Particles
  * are injected and their positions are saved to disk. The code
  * demonstrates the use of Guo off lattice boundary conditions,
  * particle injection and absorbtion as well as domain voxelization.
  **/

#include "palabos3D.h"
#include "palabos3D.hh"

using namespace plb;

typedef double T;
typedef Array<T,3> Velocity;
#define DESCRIPTOR descriptors::D3Q19Descriptor

plint xDirection  = 0;
plint yDirection  = 1;
plint borderWidth = 1;  // Because the Guo boundary condition acts in a one-cell layer.
                        // Requirement: margin>=borderWidth.
plint margin      = 1;  // Extra margin of allocated cells around the obstacle.
plint extraLayer  = 0;  // Make the bounding box larger; for visualization purposes
                        //   only. For the simulation, it is OK to have extraLayer=0.
const plint blockSize = 20; // Zero means: no sparse representation.
const plint extendedEnvelopeWidth = 2; // Because the Guo boundary condition needs 2-cell neighbor access.

plint n0 = 15; // Reference resolution.
plint ny = 80; // Resolution (cylinder diameter).
plint nz = ny;
plint nx = ny;
T radius = (T)ny/2.;
Array<T,3> originalCenter(0.0, 0.0, 0.0);
Array<T,3> inletCenter(0.0, 0.0, 0.0), outletCenter(0.0, 0.0, 0.0);
T length = (T)nx; // Cylinder length (another good choice: 4*nx).
plint nAxial = nx/2; // Two parameters needed for the creation of the triangularized cylinder surface.
plint nCirc  = 3*ny/2;

T Reynolds = 10.;
T uRef = 0.02;

T uAverage = uRef * n0 / ny;

plint maxIter  = 100000; // Maximum number of iterations for the simulation.
plint outIter  = 100; // Number of iterations for printing the average kinetic energy on the screen.
plint saveIter = 100; // Number of iterations for saving data in the disk.

bool useAllDirections = true; // Extrapolation scheme for the off lattice boundary condition.
bool useRegularized = true; // Use an off lattice boundary condition which is closer in spirit to
                            //   regularized boundary conditions.

plint particleTimeFactor = 1;  // If this variable has value 2, this means that the particles have
                               //   a time step two times bigger than the fluid.
T particleProbabilityPerCell = 2.5e-3;   // Probability of injecting a particle into an injection cell at each time step.
T cutOffSpeedSqr = 1.e-8; // Criterion to eliminate particles with very small velocity.

// This function object is a domain functional which is used for injecting particles from
//   an area close to the inlet of the tube.
class CircularInjection {
public:
    CircularInjection ( T radius_, Array<T,3> center_ )
        : radius(radius_),
          center(center_)
    { }
    bool operator()(Array<T,3> const& pos) const {
        return (util::sqr(pos[1]-center[1]) + util::sqr(pos[2]-center[2]))
               < util::sqr(radius);
    }
private:
    T radius;
    Array<T,3> center;
};

void iniLattice( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                 VoxelizedDomain3D<T>& voxelizedDomain )
{
    // Switch all remaining outer cells to no-dynamics, except the outer
    //   boundary layer, and keep the rest as BGKdynamics.
    defineDynamics(lattice, voxelizedDomain.getVoxelMatrix(), lattice.getBoundingBox(),
                   new NoDynamics<T,DESCRIPTOR>, voxelFlag::outside);
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), (T) 1., Array<T,3>((T)0.,(T)0.,(T)0.));
    lattice.initialize();
}

int main(int argc, char* argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./");

    // Create the cylinder surface as a set of triangles.
    TriangleSet<T> triangleSet;
    triangleSet = constructCylinder<T>(originalCenter, radius, radius, length, nAxial, nCirc);

    // The next few lines of code are typical. They transform the surface geometry of the
    //   tube to more efficient data structures that are internally used by palabos.
    //   The TriangleBoundary3D structure will be later used to assign proper boundary conditions.
    DEFscaledMesh<T> defMesh(triangleSet, ny, yDirection, margin, extraLayer);
    defMesh.getMesh().inflate();
    TriangleBoundary3D<T> boundary(defMesh);

    pcout << "Number of inlets/outlets which were closed: "
          << boundary.getInletOutlet(xDirection).size() << std::endl;
    inletCenter = computeBaryCenter(boundary.getMesh(),boundary.getInletOutlet(xDirection)[0]);
    outletCenter = computeBaryCenter(boundary.getMesh(),boundary.getInletOutlet(xDirection)[1]);
    radius = computeInnerRadius(boundary.getMesh(),boundary.getInletOutlet(xDirection)[0]);

    pcout << "Inlet center: " << inletCenter[0] << "," << inletCenter[1] << "," << inletCenter[2] << std::endl;
    pcout << "Outlet center: " << outletCenter[0] << "," << outletCenter[1] << "," << outletCenter[2] << std::endl;

    T nu = uAverage * 2.*radius / Reynolds;
    T omega = 1./(3.*nu+0.5);
    pcout << "omega=" << omega << std::endl;

    boundary.getMesh().writeAsciiSTL("cylinder.stl");
    pcout << "Number of triangles: " << boundary.getMesh().getNumTriangles() << std::endl;

    // The tube simulation is an interior (as opposed to exterior) flow problem. For
    //   this reason, the lattice nodes that lie inside the computational domain must
    //   be identified and distinguished from the ones that lie outside of it. This is
    //   handled by the following voxelization process.
    pcout << std::endl << "Voxelizing the domain." << std::endl;
    const int flowType = voxelFlag::inside;
    VoxelizedDomain3D<T> voxelizedDomain (
            boundary, flowType, extraLayer, borderWidth, extendedEnvelopeWidth, blockSize);
    pcout << getMultiBlockInfo(voxelizedDomain.getVoxelMatrix()) << std::endl;

    std::auto_ptr<MultiBlockLattice3D<T,DESCRIPTOR> > lattice 
        = generateMultiBlockLattice<T,DESCRIPTOR> (
                voxelizedDomain.getVoxelMatrix(), extendedEnvelopeWidth, new BGKdynamics<T,DESCRIPTOR>(omega) );

    // The Guo off lattice boundary condition is set up.
    pcout << "Creating boundary condition." << std::endl;
    BoundaryProfiles3D<T,Velocity> profiles;
    profiles.defineInletOutletTags(boundary, xDirection);
    profiles.setInletOutlet (
            new PoiseuilleProfile3D<T>(uAverage),
            new PoiseuilleProfile3D<T>(-uAverage) );
    GuoOffLatticeModel3D<T,DESCRIPTOR>* model =
            new GuoOffLatticeModel3D<T,DESCRIPTOR> (
                new TriangleFlowShape3D<T,Array<T,3> > (
                    voxelizedDomain.getBoundary(), profiles),
                flowType, useAllDirections );
    model->selectUseRegularizedModel(useRegularized);
    OffLatticeBoundaryCondition3D<T,DESCRIPTOR,Velocity> boundaryCondition (
            model, voxelizedDomain, *lattice);
    boundaryCondition.insert();

    pcout << std::endl << "Initializing lattice." << std::endl;
    iniLattice(*lattice, voxelizedDomain);

    // Particles

    // Definition of a particle field.
    MultiParticleField3D<DenseParticleField3D<T,DESCRIPTOR> >* particles=0;
    particles = new MultiParticleField3D<DenseParticleField3D<T,DESCRIPTOR> > (
        lattice->getMultiBlockManagement(),
        defaultMultiBlockPolicy3D().getCombinedStatistics() );

    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(particles);

    std::vector<MultiBlock3D*> particleFluidArg;
    particleFluidArg.push_back(particles);
    particleFluidArg.push_back(lattice.get());

    // Functional that advances the particles to their new position at each
    //   predefined time step.
    integrateProcessingFunctional (
            new AdvanceParticlesEveryWhereFunctional3D<T,DESCRIPTOR>(cutOffSpeedSqr),
            lattice->getBoundingBox(), particleArg, 0);
    // Functional that assigns the particle velocity according to the particle's
    //   position in the fluid.
    integrateProcessingFunctional (
            new FluidToParticleCoupling3D<T,DESCRIPTOR>((T)particleTimeFactor),
            lattice->getBoundingBox(), particleFluidArg, 1 );

    // Definition of a domain from which particles will be injected in the flow field.
    //   The specific domain is close to the inlet of the tube.
    Box3D injectionDomain(lattice->getBoundingBox());
    injectionDomain.x0 = inletCenter[0] +2;
    injectionDomain.x1 = inletCenter[0] +4;

    // Definition of simple mass-less particles.
    Particle3D<T,DESCRIPTOR>* particleTemplate=0;
    particleTemplate = new PointParticle3D<T,DESCRIPTOR>(0, Array<T,3>(0.,0.,0.), Array<T,3>(0.,0.,0.));

    // For the sake of illustration, particles are being injected in an area close to the
    // center of the tube.
    T injectionRadius = radius/4.;
    // Functional which injects particles with predefined probability from the specified injection domain.
    integrateProcessingFunctional (
            new AnalyticalInjectRandomParticlesFunctional3D<T,DESCRIPTOR,CircularInjection> (
                particleTemplate, particleProbabilityPerCell, CircularInjection(injectionRadius, inletCenter) ),
            injectionDomain, particleArg, 0 );

    // Definition of an absorbtion domain for the particles. The specific domain is very close to the
    //   exit of the tube.
    Box3D absorbtionDomain(lattice->getBoundingBox());
    absorbtionDomain.x0 = outletCenter[0] -2;
    absorbtionDomain.x1 = outletCenter[0] -2;

    // Functional which absorbs the particles which reach the specified absorbtion domain.
    integrateProcessingFunctional (
            new AbsorbParticlesFunctional3D<T,DESCRIPTOR>,
            absorbtionDomain, particleArg, 0 );

    particles->executeInternalProcessors();

    bool checkForErrors = true;
    pcout << std::endl << "Starting simulation." << std::endl;
    for (plint i=0; i<maxIter; ++i) {
        if (i%outIter==0) {
            pcout << "Iteration= " << i << "; "
                  << "Average energy: "
                  << boundaryCondition.computeAverageEnergy() << std::endl;
            pcout << "Number of particles in the tube: "
                  << countParticles(*particles, particles->getBoundingBox()) << std::endl;
                  
        }
        if (i%saveIter==0 && i>0) {
            pcout << "Write visualization files." << std::endl;
            VtkImageOutput3D<T> vtkOut("volume", 1.); 
            vtkOut.writeData<float>(*boundaryCondition.computePressure(), "p", 1.);
            vtkOut.writeData<float>(*boundaryCondition.computeVelocityNorm(), "u", 1.);

            pcout << "Write particle output file." << std::endl;
            writeAsciiParticlePos(*particles, "particle_positions.dat");
            writeParticleVtk(*particles, "particles.vtk");
        }

        lattice->collideAndStream();

        if (i%particleTimeFactor==0) {
            particles->executeInternalProcessors();
        }

        if (checkForErrors) {
            abortIfErrorsOccurred();
            checkForErrors = false;
        }
    }

    delete particles;
    return 0;
}
