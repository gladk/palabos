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
  * This program solves the 3D Poiseuille flow inside a cone. Particles
  * are injected and their positions are saved in a VTK file. The code
  * demonstrates the use of Guo off lattice boundary conditions,
  * particle injection and absorbtion as well as domain voxelization.
  **/

#include "palabos3D.h"
#include "palabos3D.hh"

using namespace plb;

typedef double T;
typedef Array<T,3> Velocity;
#define DESCRIPTOR descriptors::D3Q19Descriptor

// The following are low-level constants needed to voxelize the domain and
// to create the boundary condition.
plint xDirection  = 0;
plint yDirection  = 1;
plint borderWidth = 1;
plint margin      = 1;
plint extraLayer  = 0;
// The computational domain is covered by blocks of size blockSize x blockSize x blockSize
plint blockSize = 0; // Zero means: no sparse representation.
// The parallel envelope of the fluid (2 is needed for the off-lattice boundary condition).
plint extendedEnvelopeWidth = 2;
int flowType = voxelFlag::inside;  // This is an inner-flow problem.
bool useAllDirections = true;      // Extrapolation scheme for the off lattice boundary condition.
bool useRegularized = true;        // Use an off lattice boundary condition which is closer in spirit to
                                   //   regularized boundary conditions.

plint coneFactor  = 3;  // Ratio between cone outlet and inlet.
plint lengthRatio = 3;  // Ratio between cone length and inlet diameter.
plint n0 = 20;          // Resolution of the cone inlet.

// Geometrical parameters of the cone, in lattice units, to be computed.
T radius = (T)0.;
Array<T,3> inletCenter(0.0, 0.0, 0.0), outletCenter(0.0, 0.0, 0.0);
Box3D injectionDomain;

T Reynolds = 50.;    // Reynolds number wrt. inlet diameter.
T uInject = 0.02;    // Injection velocity in lattice units.
T omega    = 0.;     // Relaxation parameter.

plint maxIter  = 1000000; // Maximum number of iterations for the simulation.
plint outIter  = 100;     // Iteration intervals for printing the average kinetic energy on the screen.
plint saveIter = 100;     // Iteration intervals for saving data in the disk.
plint maxFluidIter = 1200;      // Stop iterating the fluid after this iteration.
plint startParticleIter = 1200; // Start iterating the particles after this iteration.

T particleProbabilityPerCell = 4.e-6; // Per-cell rate of particle injection.
T fluidCompliance = 6.e-5;            // The fluid->particle force is equal to the velocity difference times the fluid-compliance.
Array<T,3> gravity(-3.0e-7, 0., 0.);  // A body force acting on each particle.
T forceAmplitude = 1.e-4;             // Amplitude for the pairwise particle interaction force.
T cutOffLength = 3.5;                 // Cutoff length for particle interaction in lattice units (make sure the interaction force
                                      // is negligible at this distance).
plint commEnvelope = (plint)cutOffLength +1;  // Width of communication envelope needed to cope with neighborhoods within the cutoff length.
                                              // Must be bigger than blockSize/2.

// The containers for the fluid and the particles, and the data structures for the
// boundary-condition.
typedef DenseParticleField3D<T,DESCRIPTOR> ParticleFieldT;
MultiBlockLattice3D<T,DESCRIPTOR>* fluidLattice                         =0;
MultiParticleField3D<ParticleFieldT>* particles    =0;
TriangleBoundary3D<T>* triangleBd                                       =0;
VoxelizedDomain3D<T>* voxelizedDomain                                   =0;
BoundaryProfiles3D<T,Velocity> profiles;
OffLatticeBoundaryCondition3D<T,DESCRIPTOR,Velocity>* boundaryCondition =0;


/* ******** ParticleInteractionForce *********************************** */

// This is the most crucial part of this code: it defines the pair-wise interaction force
// between particles. This is a data processor: it executes algorithm for pair-wise
// interaction on a specific sub-domain of the simulation, defined by the parameter
// "domain". The rule is: you must modify properly all particles inside the domain, and
// you may read from all particles in the domain plus an envelope of size "commEnvelope",
// because Palabos has implemented parallelization in such a way that for every domain
// on which a data processor is executed, an envelope of size "commEnvelope" at most contains
// data which is properly synchronized with other processors.
template<typename T, template<typename U> class Descriptor>
class ParticleInteractionForce : public BoxProcessingFunctional3D
{
public:
    ParticleInteractionForce(T forceAmplitude_, Array<T,3> gravity_)
        : forceAmplitude(forceAmplitude_),
          gravity(gravity_)
    { }
    // This is the central part of the code, where the particle interaction is defined.
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks)
    {
        PLB_PRECONDITION( blocks.size()==1 );
        ParticleField3D<T,Descriptor>& particleField = *dynamic_cast<ParticleField3D<T,Descriptor>*>(blocks[0]);
        Dot3D offset = particleField.getLocation();

        std::vector<Particle3D<T,Descriptor>*> particles;
        std::vector<Particle3D<T,Descriptor>*> neighbors;
        particleField.findParticles(domain, particles);
        // Loop over all particles assigned to this data processor.
        for (pluint iParticle=0; iParticle<particles.size(); ++iParticle) {
            Particle3D<T,Descriptor>* nonTypeParticle = particles[iParticle];
            if (nonTypeParticle->getTag()!=0) {  // Exclude the wall particles.
                VerletParticle3D<T,Descriptor>* particle =
                    dynamic_cast<VerletParticle3D<T,Descriptor>*>(nonTypeParticle);
                // Compute a neighborhood, in integer coordinates, which contains at least
                // all particles inside a radius of cutOffLength.
                Array<T,3> position = particle->getPosition();
                // Convert global particle coordinates to local coordinates for the domain
                // treated by this data processor.
                plint x = util::roundToInt(position[0])-offset.x;
                plint y = util::roundToInt(position[1])-offset.y;
                plint z = util::roundToInt(position[2])-offset.z;
                Box3D neighborhood(x-commEnvelope,x+commEnvelope,
                                   y-commEnvelope,y+commEnvelope, z-commEnvelope,z+commEnvelope);
                neighbors.clear();
                T rCritical = 2.0;
                T rCriticalSqr = util::sqr(rCritical);
                // Use the particle hash to find neighboring particles efficiently.
                particleField.findParticles(neighborhood, neighbors);
                Array<T,3> force((T)0.,(T)0.,(T)0.);
                for (pluint iNeighbor=0; iNeighbor<neighbors.size(); ++iNeighbor) {
                    Array<T,3> r = particle->getPosition() - neighbors[iNeighbor]->getPosition();
                    T rSqr = normSqr(r);
                    if (rSqr<util::sqr(cutOffLength) && rSqr>1.e-6) {
                        // The potential is amplitude*r^{-6} for r>rCritical, and r*amplitude for r<=rCritical.
                        // This is to aovid numerical instability due to huge forces when two
                        // particles come too close (for example in a badly designed initial condition).
                        if (rSqr>rCriticalSqr) {
                            // The potential is r^{-6}, so the force is r^{-7}. An additional r^{-1}
                            // is needed to normalize the vector r.
                            force += r/(rSqr*rSqr*rSqr*rSqr) * forceAmplitude;
                        }
                        else {
                            // Inside a radius of rCritical (in lattice units), the force is constant.
                            T rNorm = std::sqrt(rSqr);
                            force += r/(rNorm*rCritical*rCriticalSqr*rCriticalSqr*rCriticalSqr) * forceAmplitude;
                        }
                    }
                }
                // Particle acceleration according to Newton's law.
                particle->set_a (
                        particle->get_a() + (gravity+force) * particle->get_invRho() );
            }
        }
    }
    virtual ParticleInteractionForce<T,Descriptor>* clone() const
    {
        return new ParticleInteractionForce<T,Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::dynamicVariables;
    }
private:
    T forceAmplitude;
    Array<T,3> gravity;
};


// This function object defines the area (here a circle) in which the particles are injected.
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

// In this part of the code, the triangular surface mesh of the cone is constructed, and the
// domain is voxelized: a flag matrix is created which tags cells that are inside resp. outside
// the cone.
void buildConeGeometry() {
    // Create the cone surface as a set of triangles.
    TriangleSet<T> triangleSet;

    plint nAxial = n0*lengthRatio/2;  // Triangle resolution along the length.
    plint nCirc  = 3*n0*coneFactor/2; // Triangle resolution along the circumference.
    // Create the cone in arbitrary units, as it will be rescaled afterwards.
    triangleSet = constructCylinder<T>(Array<T,3>(0.,0.,0.), 1., (T)coneFactor, (T)lengthRatio*2., nAxial, nCirc);

    // Transform the surface geometry of the tube to more efficient data structures that are internally used by palabos.
    DEFscaledMesh<T> defMesh(triangleSet, n0*coneFactor, yDirection, margin, extraLayer);
    defMesh.getMesh().inflate();
    triangleBd = new TriangleBoundary3D<T>(defMesh);

    pcout << "Number of inlets/outlets which were closed: "
          << triangleBd->getInletOutlet(xDirection).size() << std::endl;

    // Write a surface mesh that can be visualized in Paraview, along with the particles.
    triangleBd->getMesh().writeBinarySTL("cone.stl");
    pcout << "Number of triangles: " << triangleBd->getMesh().getNumTriangles() << std::endl;
    pcout << "Number of vertices: " << triangleBd->getMesh().getNumVertices() << std::endl;

    // The tube simulation is an interior (as opposed to exterior) flow problem. For
    //   this reason, the lattice nodes that lie inside the computational domain must
    //   be identified and distinguished from the ones that lie outside of it. This is
    //   handled by the following voxelization process.
    pcout << std::endl << "Voxelizing the domain." << std::endl;
    voxelizedDomain = new VoxelizedDomain3D<T> (
            *triangleBd, flowType, extraLayer, borderWidth, extendedEnvelopeWidth, blockSize);
    pcout << getMultiBlockInfo(voxelizedDomain->getVoxelMatrix()) << std::endl;
    // At this point, the MultiScalarField voxelizedDomain->getVoxelMatrix() contains the cell tags.
}

// Compute the dimensions of the cone in lattice units, and the relaxation parameter omega for the fluid.
void computePhysicalParameters()
{
    inletCenter = computeBaryCenter(triangleBd->getMesh(),triangleBd->getInletOutlet(xDirection)[0]);
    outletCenter = computeBaryCenter(triangleBd->getMesh(),triangleBd->getInletOutlet(xDirection)[1]);
    radius = computeInnerRadius(triangleBd->getMesh(),triangleBd->getInletOutlet(xDirection)[0]);

    pcout << "Inlet center: " << inletCenter[0] << "," << inletCenter[1] << "," << inletCenter[2] << std::endl;
    pcout << "Outlet center: " << outletCenter[0] << "," << outletCenter[1] << "," << outletCenter[2] << std::endl;

    T nu = uInject * 2.*radius / Reynolds;
    omega = 1./(3.*nu+0.5);
    pcout << "omega=" << omega << std::endl;
}

// Instantiate the multi-block lattices for the fluid and the particles, and create
// the fluid boundary condition on the cone.
void createLattices() {
    // Allocate the data for the fluid.
    fluidLattice = generateMultiBlockLattice<T,DESCRIPTOR> (
                  voxelizedDomain->getVoxelMatrix(), extendedEnvelopeWidth,
                  new BGKdynamics<T,DESCRIPTOR>(omega) ).release();

    // Set up the Guo off-lattice boundary condition on the walls of the cone.
    pcout << "Creating boundary condition." << std::endl;
    // In order to distinguish the inlet from the outlet, say that they are sorted
    // in x-direction.
    profiles.defineInletOutletTags(*triangleBd, xDirection);
    // As previously said, things are considered in increasing x-direction. In this
    // direction, the first opening (the inlet) shall implement a velocity plug-profile,
    // and the second opening (the outlet), a Neumann velocity profile with fixed
    // density.
    profiles.setInletOutlet (
            new VelocityPlugProfile3D<T>(uInject),
            new DensityNeumannBoundaryProfile3D<T>() );
    // The algorithm for the off-lattice boundary condition is the Guo one.
    GuoOffLatticeModel3D<T,DESCRIPTOR>* model =
            new GuoOffLatticeModel3D<T,DESCRIPTOR> (
                new TriangleFlowShape3D<T,Array<T,3> > (
                    voxelizedDomain->getBoundary(), profiles),
                flowType, useAllDirections );
    model->selectUseRegularizedModel(useRegularized);
    boundaryCondition =
        new OffLatticeBoundaryCondition3D<T,DESCRIPTOR,Velocity> (
                model, *voxelizedDomain, *fluidLattice );
    // Insert the boundary condition into the lattice (by means of a data processor)
    // so it is executed at the end of every collision-streaming cycle.
    boundaryCondition->insert();

    pcout << std::endl << "Initializing lattice." << std::endl;
    // Switch all remaining outer cells to NoDynamics (for computational efficiency),
    // except the outer boundary layer which is needed to implement the boundary condition,
    // and keep the rest as BGKdynamics.
    defineDynamics(*fluidLattice, voxelizedDomain->getVoxelMatrix(), fluidLattice->getBoundingBox(),
                   new NoDynamics<T,DESCRIPTOR>, voxelFlag::outside);
    // Default-initialize at equilibrium with zero-velocity and fixed pressure.
    initializeAtEquilibrium(*fluidLattice, fluidLattice->getBoundingBox(), (T)1., Array<T,3>((T)0.,(T)0.,(T)0.));
    // Execute all data processors once to start the simulation off with well-defined initial values.
    fluidLattice->initialize();

    // Create a multi-block distribution for the particles which is the same as the one for
    // the lattice. This is necessary, because afterwards a fluid->particle coupling is
    // implemented.
    MultiBlockManagement3D const& latticeManagement(fluidLattice->getMultiBlockManagement());
    MultiBlockManagement3D particleManagement (
            latticeManagement.getSparseBlockStructure(),
            latticeManagement.getThreadAttribution().clone(),
            commEnvelope, latticeManagement.getRefinementLevel() );

    // Allocation of the data for the particles. This only allocates the particle-hash
    // (a grid on which the particles are stored). The particles themselves are injected
    // later on.
    particles = new MultiParticleField3D<ParticleFieldT> (
        particleManagement, defaultMultiBlockPolicy3D().getCombinedStatistics() );
}

// Put in place the particle-particle, particle-wall, and fluid-particle interactions,
// and the algorithms for the Verlet integration of the particle motion.
void setupCouplings()
{
    // Add a particle on each vertex of the cone wall. These particles don't move, but
    // they interact with the injected particles, so that these rebounce from the wall.
    addWallParticlesGeneric<T,DESCRIPTOR,ParticleFieldT>(*particles, *triangleBd);

    // In the following the data processors for the equations of motion and the
    // interaction terms are manually added to the particle field. In other situations,
    // data processors are added implicitly, such as for example the data processors
    // for the implementation of a boundary condition in the function call
    // boundaryCondition->insert() above.

    // Prepare the arguments to be provided to the data processors: the particle field
    // for Verlet integration, and particles+fluid for the coupling terms.
    std::vector<MultiBlock3D*> particleArg;
    particleArg.push_back(particles);

    std::vector<MultiBlock3D*> particleFluidArg;
    particleFluidArg.push_back(particles);
    particleFluidArg.push_back(fluidLattice);

    // Calls the function advance() on each particle which, for Verlet particles,
    // updates the position according to v(t+0.5)
    integrateProcessingFunctional (
            new AdvanceParticlesEveryWhereFunctional3D<T,DESCRIPTOR>(-1.0),
            fluidLattice->getBoundingBox(), particleArg, 0);
    // Compute fluid->particle force (which includes friction).
    integrateProcessingFunctional (
            new FluidToParticleCoupling3D<T,DESCRIPTOR>(1.),
            fluidLattice->getBoundingBox(), particleFluidArg, 1 );
    // Compute particle->particle and gravity forces, with help of the data processor
    // ParticleInteractionForce defined above.
    integrateProcessingFunctional (
            new ParticleInteractionForce<T,DESCRIPTOR>(forceAmplitude, gravity),
            fluidLattice->getBoundingBox(), particleArg, 1 );
    // Integrate the velocity according to the Verlet algorithm.
    integrateProcessingFunctional (
            new VerletUpdateVelocitySelective3D<T,DESCRIPTOR>(new util::SelectLargerEqualInt(1)),
            fluidLattice->getBoundingBox(), particleArg, 1 );

    // Definition of a domain from which particles will be injected in the flow field.
    //   The specific domain is close to the inlet of the tube.
    injectionDomain = Box3D(fluidLattice->getBoundingBox());
    injectionDomain.x0 = inletCenter[0] +3;
    injectionDomain.x1 = inletCenter[0] +5;

    // Definition of an absorbtion domains for the particles, one at the bottom and one at the top.
    Box3D absorbtionDomain1(fluidLattice->getBoundingBox());
    absorbtionDomain1.x0 = outletCenter[0] -2;
    absorbtionDomain1.x1 = outletCenter[0] -2;
    Box3D absorbtionDomain2(fluidLattice->getBoundingBox());
    absorbtionDomain2.x0 = inletCenter[0] +1;
    absorbtionDomain2.x1 = inletCenter[0] +1;

    // Functional which absorbs the particles in the specified domains, to avoid that they leave the
    // computational domain.
    integrateProcessingFunctional (
            new AbsorbParticlesFunctionalSelective3D<T,DESCRIPTOR>(new util::SelectLargerEqualInt(1)),
            absorbtionDomain1, particleArg, 0 );
    integrateProcessingFunctional (
            new AbsorbParticlesFunctionalSelective3D<T,DESCRIPTOR>(new util::SelectLargerEqualInt(1)),
            absorbtionDomain2, particleArg, 0 );

    // Execute all data processors once to start the simulation off with well-defined initial values.
    particles->executeInternalProcessors();
}


int main(int argc, char* argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./");
    srand(global::mpi().getRank()+10);

    // Set up the problem.
    buildConeGeometry();
    computePhysicalParameters();
    createLattices();
    setupCouplings();

    // Create a template for the particles to be injected.
    VerletParticle3D<T,DESCRIPTOR> *particleTemplate=0;
    particleTemplate = new VerletParticle3D<T,DESCRIPTOR>(1, Array<T,3>(0.,0.,0.));
    particleTemplate->setFluidCompliance(fluidCompliance);

    pcout << std::endl << "Starting simulation." << std::endl;
    bool checkForErrors = true;
    for (plint i=0; i<maxIter; ++i) {
        // Start injecting particles, and iterating them, once the fluid is stationary (it
        // would also be OK to start injecting particles earlier, if you prefer).
        if (i>=startParticleIter) {
            std::vector<MultiBlock3D*> particleArg;
            particleArg.push_back(particles);
            // Progressively increase the injection rate.
            T probability = particleProbabilityPerCell;
            //T probability = particleProbabilityPerCell*(1.+(T)i/6000.);
            Particle3D<T,DESCRIPTOR>* particle = particleTemplate->clone();
            // Attribute a color to the first particles by changing their tags (the tags
            // are written in the VTK file).
            if (i<=startParticleIter+20000) {
                particle->setTag(i-startParticleIter);
            }
            applyProcessingFunctional (
                    new AnalyticalInjectRandomParticlesFunctional3D<T,DESCRIPTOR,CircularInjection> (
                        particle, probability, CircularInjection(radius/2.0, inletCenter) ),
                    injectionDomain, particleArg );
        }
        // For efficiency reasons, stop iterating the fluid once it is assumed to have reached a stationary state.
        if (i<maxFluidIter) {
            // Iterate the fluid.
            fluidLattice->collideAndStream();
        }
        if (i>=startParticleIter) {
            // Iterate the particles..
            particles->executeInternalProcessors();
        }

        if (checkForErrors) {
            abortIfErrorsOccurred();
            checkForErrors = false;
        }

        if (i%outIter==0) {
            pcout << "Iteration= " << i << "; "
                  << "Average energy: "
                  << boundaryCondition->computeAverageEnergy() << std::endl;
            pcout << "Number of particles in the tube: "
                  << countParticles(*particles, particles->getBoundingBox()) << std::endl;
        }
        if (i == maxFluidIter) {
            pcout << "Writing fluid data.." << std::endl;
            VtkImageOutput3D<T> vtkOut("volume", 1.); 
            vtkOut.writeData<float>(*boundaryCondition->computePressure(), "p", 1.);
            vtkOut.writeData<float>(*boundaryCondition->computeVelocityNorm(), "u", 1.);
        }

        if (i%saveIter==0 && i > startParticleIter) {
            pcout << "Writing particle data.." << std::endl;
            std::string fName = createFileName("particles_",i,8)+".vtk";
            writeSelectedParticleVtk<T,DESCRIPTOR>(*particles, fName, particles->getBoundingBox(), util::SelectLargerEqualInt(1));
        }
    }
}

