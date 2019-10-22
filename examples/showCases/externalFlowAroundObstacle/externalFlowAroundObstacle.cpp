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

/* \file
 * External flow around a 3D obstacle.
 * This example demonstrates many features of Palabos:
 * Loading geometries from STL files.
 * Using the voxelizer.
 * Using off-lattice boundary conditions.
 * Imposing sophisticated outflow boundary conditions.
 * Using sponge zones.
 * Imposing time dependent inlet boundary conditions.
 * Integrating particles for the computation of streamlines.
 * Computing the force on objects.
 * */

#include "palabos3D.h"
#include "palabos3D.hh"

using namespace plb;
using namespace std;

typedef double T;
typedef Array<T,3> Velocity;
#define DESCRIPTOR descriptors::D3Q19Descriptor
typedef DenseParticleField3D<T,DESCRIPTOR> ParticleFieldT;

#define PADDING 8

static std::string outputDir("./tmp/");


// This structure holds all the user-defined parameters, and some
// derived values needed for the simulation.
struct Param
{
    T nu;                               // Kinematic viscosity.
    T lx, ly, lz;                       // Size of computational domain, in physical units.
    T cx, cy, cz;                       // Position of the center of the obstacle, in physical units.
    plint cxLB, cyLB, czLB;             // Position of the center of the obstacle, in lattice units.
    bool freeSlipWall;                  // Use free-slip condition on obstacle, as opposed to no-slip?
    bool lateralFreeSlip;               // Use free-slip lateral boundaries or periodic ones?
    T maxT, statT, imageT, vtkT;        // Time, in physical units, at which events occur.
    plint resolution;                   // Number of lattice nodes along a reference length.
    T inletVelocity;                    // Inlet x-velocity in physical units.
    T uLB;                              // Velocity in lattice units (numerical parameters).
    bool useSmago;                      // Use a Smagorinsky LES model or not.
    T cSmago;                           // Parameter for the Smagorinsky LES model.
    plint nx, ny, nz;                   // Grid resolution of bounding box.
    T omega;                            // Relaxation parameter.
    T dx, dt;                           // Discrete space and time steps.
    plint maxIter, statIter;            // Time for events in lattice units.
    plint imageIter, vtkIter;
    bool useParticles;                  // Simulate particles or not.
    int particleTimeFactor;             // If the particle time factor is 2, then the integration time step
                                        //   for the particles is twice that of the fluid.
    T particleProbabilityPerCell;       // Probability of injection of a particle at an injection cell at each time step.
    T cutOffSpeedSqr;                   // Criterion to eliminate particles with very small velocity.
    int maxNumParticlesToWrite;         // Maximum number of particles in the output VTK files.

    T outletSpongeZoneWidth;            // Width of the outlet sponge zone.
    plint numOutletSpongeCells;         // Number of the lattice nodes contained in the outlet sponge zone.
    int outletSpongeZoneType;           // Type of the outlet sponge zone (Viscosity or Smagorinsky).
    T targetSpongeCSmago;               // Target Smagorinsky parameter at the end of the Smagorinsky sponge Zone.
    plint initialIter;                  // Number of initial iterations until the inlet velocity reaches its final value.

    Box3D inlet, outlet, lateral1;      // Outer domain boundaries in lattice units.
    Box3D lateral2, lateral3, lateral4;

    std::string geometry_fname;

    Param()
    { }

    Param(std::string xmlFname)
    {
        XMLreader document(xmlFname);
        document["geometry"]["filename"].read(geometry_fname);
        document["geometry"]["center"]["x"].read(cx);
        document["geometry"]["center"]["y"].read(cy);
        document["geometry"]["center"]["z"].read(cz);
        document["geometry"]["freeSlipWall"].read(freeSlipWall);
        document["geometry"]["lateralFreeSlip"].read(lateralFreeSlip);
        document["geometry"]["domain"]["x"].read(lx);
        document["geometry"]["domain"]["y"].read(ly);
        document["geometry"]["domain"]["z"].read(lz);

        document["numerics"]["nu"].read(nu);
        document["numerics"]["inletVelocity"].read(inletVelocity);
        document["numerics"]["resolution"].read(resolution);
        document["numerics"]["uLB"].read(uLB);
        document["numerics"]["useSmago"].read(useSmago);
        if (useSmago) {
            document["numerics"]["cSmago"].read(cSmago);
        }

        document["numerics"]["useParticles"].read(useParticles);
        if (useParticles) {
            document["numerics"]["particleTimeFactor"].read(particleTimeFactor);
            document["numerics"]["particleProbabilityPerCell"].read(particleProbabilityPerCell);
            document["numerics"]["cutOffSpeedSqr"].read(cutOffSpeedSqr);
            document["numerics"]["maxNumParticlesToWrite"].read(maxNumParticlesToWrite);
        }

        document["numerics"]["outletSpongeZoneWidth"].read(outletSpongeZoneWidth);
        std::string zoneType;
        document["numerics"]["outletSpongeZoneType"].read(zoneType);
        if ((util::tolower(zoneType)).compare("viscosity") == 0) {
            outletSpongeZoneType = 0;
        } else if ((util::tolower(zoneType)).compare("smagorinsky") == 0) {
            outletSpongeZoneType = 1;
        } else {
            pcout << "The sponge zone type must be either \"Viscosity\" or \"Smagorinsky\"." << std::endl;
            exit(-1);
        }
        document["numerics"]["targetSpongeCSmago"].read(targetSpongeCSmago);

        document["numerics"]["initialIter"].read(initialIter);

        document["output"]["maxT"].read(maxT);
        document["output"]["statT"].read(statT);
        document["output"]["imageT"].read(imageT);
        document["output"]["vtkT"].read(vtkT);

        computeLBparameters();
    }

    void computeLBparameters()
    {
        dx = ly / (resolution - 1.0);
        dt = (uLB/inletVelocity) * dx;
        T nuLB = nu * dt/(dx*dx);
        omega = 1.0/(DESCRIPTOR<T>::invCs2*nuLB+0.5);
        if (lateralFreeSlip) {
            nx = util::roundToInt(lx/dx) + 1;
            ny = util::roundToInt(ly/dx) + 1;
            nz = util::roundToInt(lz/dx) + 1;
        } else {
            nx = util::roundToInt(lx/dx) + 1;
            ny = util::roundToInt(ly/dx);
            nz = util::roundToInt(lz/dx);
        }
        cxLB = util::roundToInt(cx/dx);
        cyLB = util::roundToInt(cy/dx);
        czLB = util::roundToInt(cz/dx);
        maxIter   = util::roundToInt(maxT/dt);
        statIter  = util::roundToInt(statT/dt);
        imageIter = util::roundToInt(imageT/dt);
        vtkIter   = util::roundToInt(vtkT/dt);
        numOutletSpongeCells = util::roundToInt(outletSpongeZoneWidth/dx);

        inlet    = Box3D(0,      0,      0,      ny-1,   0,      nz-1);
        outlet   = Box3D(nx-1,   nx-1,   0,      ny-1,   0,      nz-1);
        lateral1 = Box3D(1,      nx-2,   0,      0,      0,      nz-1);
        lateral2 = Box3D(1,      nx-2,   ny-1,   ny-1,   0,      nz-1);
        lateral3 = Box3D(1,      nx-2,   1,      ny-2,   0,      0);
        lateral4 = Box3D(1,      nx-2,   1,      ny-2,   nz-1,   nz-1);
    }

    Box3D boundingBox() const
    {
        return Box3D(0, nx-1, 0, ny-1, 0, nz-1);
    }

    T getInletVelocity(plint iIter)
    {
        static T pi = std::acos((T) -1.0);

        if (iIter >= initialIter) {
            return uLB;
        }

        if (iIter < 0) {
            iIter = 0;
        }

        return uLB * std::sin(pi * iIter / (2.0 * initialIter));
    }
};

Param param;

// Instantiate the boundary conditions for the outer domain.
void outerDomainBoundaries(MultiBlockLattice3D<T,DESCRIPTOR> *lattice,
                           MultiScalarField3D<T> *rhoBar,
                           MultiTensorField3D<T,3> *j,
                           OnLatticeBoundaryCondition3D<T,DESCRIPTOR> *bc)
{
    Array<T,3> uBoundary(param.getInletVelocity(0), 0.0, 0.0);

    if (param.lateralFreeSlip) {
        pcout << "Free-slip lateral boundaries." << std::endl;

        lattice->periodicity().toggleAll(false);
        rhoBar->periodicity().toggleAll(false);
        j->periodicity().toggleAll(false);

        bc->setVelocityConditionOnBlockBoundaries(*lattice, param.inlet, boundary::dirichlet);
        setBoundaryVelocity(*lattice, param.inlet, uBoundary);

        bc->setVelocityConditionOnBlockBoundaries(*lattice, param.lateral1, boundary::freeslip);
        bc->setVelocityConditionOnBlockBoundaries(*lattice, param.lateral2, boundary::freeslip);
        bc->setVelocityConditionOnBlockBoundaries(*lattice, param.lateral3, boundary::freeslip);
        bc->setVelocityConditionOnBlockBoundaries(*lattice, param.lateral4, boundary::freeslip);
        setBoundaryVelocity(*lattice, param.lateral1, uBoundary);
        setBoundaryVelocity(*lattice, param.lateral2, uBoundary);
        setBoundaryVelocity(*lattice, param.lateral3, uBoundary);
        setBoundaryVelocity(*lattice, param.lateral4, uBoundary);

        // The VirtualOutlet is a sophisticated outflow boundary condition.
        Box3D globalDomain(lattice->getBoundingBox());
        std::vector<MultiBlock3D*> bcargs;
        bcargs.push_back(lattice);
        bcargs.push_back(rhoBar);
        bcargs.push_back(j);
        T outsideDensity = 1.0;
        int bcType = 1;
        integrateProcessingFunctional(new VirtualOutlet<T,DESCRIPTOR>(outsideDensity, globalDomain, bcType),
                param.outlet, bcargs, 2);
        setBoundaryVelocity(*lattice, param.outlet, uBoundary);
    } else {
        pcout << "Periodic lateral boundaries." << std::endl;

        lattice->periodicity().toggleAll(true);
        rhoBar->periodicity().toggleAll(true);
        j->periodicity().toggleAll(true);

        lattice->periodicity().toggle(0, false);
        rhoBar->periodicity().toggle(0, false);
        j->periodicity().toggle(0, false);

        bc->addVelocityBoundary0N(param.inlet, *lattice);
        setBoundaryVelocity(*lattice, param.inlet, uBoundary);

        // The VirtualOutlet is a sophisticated outflow boundary condition.
        // The "globalDomain" argument for the boundary condition must be
        // bigger than the actual bounding box of the simulation for
        // the directions which are periodic.
        Box3D globalDomain(lattice->getBoundingBox());
        globalDomain.y0 -= 2; // y-periodicity
        globalDomain.y1 += 2;
        globalDomain.z0 -= 2; // z-periodicity
        globalDomain.z1 += 2;
        std::vector<MultiBlock3D*> bcargs;
        bcargs.push_back(lattice);
        bcargs.push_back(rhoBar);
        bcargs.push_back(j);
        T outsideDensity = 1.0;
        int bcType = 1;
        integrateProcessingFunctional(new VirtualOutlet<T,DESCRIPTOR>(outsideDensity, globalDomain, bcType),
                param.outlet, bcargs, 2);
        setBoundaryVelocity(*lattice, param.outlet, uBoundary);
    }
}

// Write VTK file for the flow around the obstacle, to be viewed with Paraview.
void writeVTK(OffLatticeBoundaryCondition3D<T,DESCRIPTOR,Velocity>& bc, plint iT)
{
    VtkImageOutput3D<T> vtkOut(createFileName("volume", iT, PADDING));
    vtkOut.writeData<float>( *bc.computeVelocityNorm(param.boundingBox()),
                             "velocityNorm", param.dx/param.dt );
    vtkOut.writeData<3,float>(*bc.computeVelocity(param.boundingBox()), "velocity", param.dx/param.dt);
    vtkOut.writeData<float>( *bc.computePressure(param.boundingBox()),
                             "pressure", param.dx*param.dx/(param.dt*param.dt) );
}

// Write PPM images on slices.
void writePPM(OffLatticeBoundaryCondition3D<T,DESCRIPTOR,Velocity>& bc, plint iT)
{
    Box3D xSlice(param.cxLB, param.cxLB, 0,          param.ny-1, 0,          param.nz-1);
    Box3D ySlice(0,          param.nx-1, param.cyLB, param.cyLB, 0,          param.nz-1);
    Box3D zSlice(0,          param.nx-1, 0,          param.ny-1, param.czLB, param.czLB);

    ImageWriter<T> writer("leeloo");
    writer.writeScaledPpm(createFileName("vnorm_xslice", iT, PADDING), *bc.computeVelocityNorm(xSlice));
    writer.writeScaledPpm(createFileName("vnorm_yslice", iT, PADDING), *bc.computeVelocityNorm(ySlice));
    writer.writeScaledPpm(createFileName("vnorm_zslice", iT, PADDING), *bc.computeVelocityNorm(zSlice));
}

void runProgram()
{
    /*
     * Read the obstacle geometry.
     */

    pcout << std::endl << "Reading STL data for the obstacle geometry." << std::endl;
    Array<T,3> center(param.cx, param.cy, param.cz);
    Array<T,3> centerLB(param.cxLB, param.cyLB, param.czLB);
    // The triangle-set defines the surface of the geometry.
    TriangleSet<T> triangleSet(param.geometry_fname, DBL);

    // Place the obstacle in the correct place in the simulation domain.
    // Here the "geometric center" of the obstacle is computed manually,
    // by computing first its bounding cuboid. In cases that the STL
    // file with the geometry of the obstacle contains its center as
    // the point, say (0, 0, 0), then the following variable
    // "obstacleCenter" must be set to (0, 0, 0) manually.
    Cuboid<T> bCuboid = triangleSet.getBoundingCuboid();
    Array<T,3> obstacleCenter = (T) 0.5 * (bCuboid.lowerLeftCorner + bCuboid.upperRightCorner);
    triangleSet.translate(-obstacleCenter);
    triangleSet.scale(1.0/param.dx); // In lattice units from now on...
    triangleSet.translate(centerLB);
    triangleSet.writeBinarySTL(outputDir+"obstacle_LB.stl");

    // The DEFscaledMesh, and the triangle-boundary are more sophisticated data
    // structures used internally by Palabos to treat the boundary.
    plint xDirection = 0;
    plint borderWidth = 1;      // Because Guo acts in a one-cell layer.
                                // Requirement: margin>=borderWidth.
    plint margin = 1;           // Extra margin of allocated cells around the obstacle, for the case of moving walls.
    plint blockSize = 0;        // Size of blocks in the sparse/parallel representation.
                                // Zero means: don't use sparse representation.
    DEFscaledMesh<T> defMesh(triangleSet, 0, xDirection, margin, Dot3D(0, 0, 0));
    TriangleBoundary3D<T> boundary(defMesh);
    //boundary.getMesh().inflate();

    pcout << "tau = " << 1.0/param.omega << std::endl;
    pcout << "dx = " << param.dx << std::endl;
    pcout << "dt = " << param.dt << std::endl;
    pcout << "Number of iterations in an integral time scale: " << (plint) (1.0/param.dt) << std::endl;

    /*
     * Voxelize the domain.
     */

    // Voxelize the domain means: decide which lattice nodes are inside the obstacle and which are outside.
    pcout << std::endl << "Voxelizing the domain." << std::endl;
    plint extendedEnvelopeWidth = 2;   // Extrapolated off-lattice BCs.
    const int flowType = voxelFlag::outside;
    VoxelizedDomain3D<T> voxelizedDomain (
            boundary, flowType, param.boundingBox(), borderWidth, extendedEnvelopeWidth, blockSize );
    pcout << getMultiBlockInfo(voxelizedDomain.getVoxelMatrix()) << std::endl;

    /*
     * Generate the lattice, the density and momentum blocks.
     */

    pcout << "Generating the lattice, the rhoBar and j fields." << std::endl;
    MultiBlockLattice3D<T,DESCRIPTOR> *lattice = new MultiBlockLattice3D<T,DESCRIPTOR>(voxelizedDomain.getVoxelMatrix());
    if (param.useSmago) {
        defineDynamics(*lattice, lattice->getBoundingBox(),
                new SmagorinskyBGKdynamics<T,DESCRIPTOR>(param.omega, param.cSmago));
        pcout << "Using Smagorinsky BGK dynamics." << std::endl;
    } else {
        defineDynamics(*lattice, lattice->getBoundingBox(),
                new BGKdynamics<T,DESCRIPTOR>(param.omega));
        pcout << "Using BGK dynamics." << std::endl;
    }
    bool velIsJ = false;
    defineDynamics(*lattice, voxelizedDomain.getVoxelMatrix(), lattice->getBoundingBox(),
            new NoDynamics<T,DESCRIPTOR>(), voxelFlag::inside);
    lattice->toggleInternalStatistics(false);

    // The rhoBar and j fields are used at both the collision and at the implementation of the
    // outflow boundary condition.
    plint envelopeWidth = 1;
    MultiScalarField3D<T> *rhoBar = generateMultiScalarField<T>((MultiBlock3D&) *lattice, envelopeWidth).release();
    rhoBar->toggleInternalStatistics(false);

    MultiTensorField3D<T,3> *j = generateMultiTensorField<T,3>((MultiBlock3D&) *lattice, envelopeWidth).release();
    j->toggleInternalStatistics(false);

    std::vector<MultiBlock3D*> lattice_rho_bar_j_arg;
    lattice_rho_bar_j_arg.push_back(lattice);
    lattice_rho_bar_j_arg.push_back(rhoBar);
    lattice_rho_bar_j_arg.push_back(j);
    integrateProcessingFunctional(
            new ExternalRhoJcollideAndStream3D<T,DESCRIPTOR>(),
            lattice->getBoundingBox(), lattice_rho_bar_j_arg, 0);
    integrateProcessingFunctional(
            new BoxRhoBarJfunctional3D<T,DESCRIPTOR>(),
            lattice->getBoundingBox(), lattice_rho_bar_j_arg, 3); // rhoBar and j are computed at level 3 because
                                                                  // the boundary conditions are on levels 0, 1 and 2.

    /*
     * Generate the off-lattice boundary condition on the obstacle and the outer-domain boundary conditions.
     */

    pcout << "Generating boundary conditions." << std::endl;

    OffLatticeBoundaryCondition3D<T,DESCRIPTOR,Velocity> *boundaryCondition;
    
    BoundaryProfiles3D<T,Velocity> profiles;
    bool useAllDirections=true;
    OffLatticeModel3D<T,Velocity>* offLatticeModel=0;
    if (param.freeSlipWall) {
        profiles.setWallProfile(new FreeSlipProfile3D<T>);
    }
    else {
        profiles.setWallProfile(new NoSlipProfile3D<T>);
    }
    offLatticeModel =
         new GuoOffLatticeModel3D<T,DESCRIPTOR> (
            new TriangleFlowShape3D<T,Array<T,3> >(voxelizedDomain.getBoundary(), profiles),
            flowType, useAllDirections );
    offLatticeModel->setVelIsJ(velIsJ);
    boundaryCondition = new OffLatticeBoundaryCondition3D<T,DESCRIPTOR,Velocity>(
            offLatticeModel, voxelizedDomain, *lattice);
    
    boundaryCondition->insert();

    // The boundary condition algorithm or the outer domain.
    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* outerBoundaryCondition
        = createLocalBoundaryCondition3D<T,DESCRIPTOR>();
    outerDomainBoundaries(lattice, rhoBar, j, outerBoundaryCondition);

    /*
     * Implement the outlet sponge zone.
     */

    if (param.numOutletSpongeCells > 0) {
        T bulkValue;
        Array<plint,6> numSpongeCells;

        if (param.outletSpongeZoneType == 0) {
            pcout << "Generating an outlet viscosity sponge zone." << std::endl;
            bulkValue = param.omega;
        } else if (param.outletSpongeZoneType == 1) {
            pcout << "Generating an outlet Smagorinsky sponge zone." << std::endl;
            bulkValue = param.cSmago;
        } else {
            pcout << "Error: unknown type of sponge zone." << std::endl;
            exit(-1);
        }

        // Number of sponge zone lattice nodes at all the outer domain boundaries.
        // So: 0 means the boundary at x = 0
        //     1 means the boundary at x = nx-1
        //     2 means the boundary at y = 0
        //     and so on...
        numSpongeCells[0] = 0;
        numSpongeCells[1] = param.numOutletSpongeCells;
        numSpongeCells[2] = 0;
        numSpongeCells[3] = 0;
        numSpongeCells[4] = 0;
        numSpongeCells[5] = 0;

        std::vector<MultiBlock3D*> args;
        args.push_back(lattice);

        if (param.outletSpongeZoneType == 0) {
            applyProcessingFunctional(new ViscositySpongeZone3D<T,DESCRIPTOR>(
                        param.nx, param.ny, param.nz, bulkValue, numSpongeCells),
                    lattice->getBoundingBox(), args);
        } else {
            applyProcessingFunctional(new SmagorinskySpongeZone3D<T,DESCRIPTOR>(
                        param.nx, param.ny, param.nz, bulkValue, param.targetSpongeCSmago, numSpongeCells),
                    lattice->getBoundingBox(), args);
        }
    }

    /*
     * Setting the initial conditions.
     */

    // Initial condition: Constant pressure and velocity-at-infinity everywhere.
    Array<T,3> uBoundary(param.getInletVelocity(0), (T)0.0, (T)0.0);
    initializeAtEquilibrium(*lattice, lattice->getBoundingBox(), (T)1.0, uBoundary);
    applyProcessingFunctional(
            new BoxRhoBarJfunctional3D<T,DESCRIPTOR>(),
            lattice->getBoundingBox(), lattice_rho_bar_j_arg); // Compute rhoBar and j before VirtualOutlet is executed.
    //lattice->executeInternalProcessors(1); // Execute all processors except the ones at level 0.
    //lattice->executeInternalProcessors(2);
    //lattice->executeInternalProcessors(3);

    /*
     * Particles (streamlines).
     */

    // This part of the code that relates to particles, is purely for visualization
    // purposes. Particles are used to compute streamlines essentially.

    // Definition of a particle field.
    MultiParticleField3D<ParticleFieldT>* particles = 0;

    if (param.useParticles) {
        particles = new MultiParticleField3D<ParticleFieldT> (
            lattice->getMultiBlockManagement(),
            defaultMultiBlockPolicy3D().getCombinedStatistics() );

        std::vector<MultiBlock3D*> particleArg;
        particleArg.push_back(particles);

        std::vector<MultiBlock3D*> particleFluidArg;
        particleFluidArg.push_back(particles);
        particleFluidArg.push_back(lattice);

        // Functional that advances the particles to their new position at each predefined time step.
        integrateProcessingFunctional (
                new AdvanceParticlesEveryWhereFunctional3D<T,DESCRIPTOR>(param.cutOffSpeedSqr),
                lattice->getBoundingBox(), particleArg, 0);
        // Functional that assigns the particle velocity according to the particle's position in the fluid.
        integrateProcessingFunctional (
                new FluidToParticleCoupling3D<T,DESCRIPTOR>((T) param.particleTimeFactor),
                lattice->getBoundingBox(), particleFluidArg, 1 );

        // Definition of a domain from which particles will be injected in the flow field.
        Box3D injectionDomain(0, 0, centerLB[1]-0.25*param.ny, centerLB[1]+0.25*param.ny,
                centerLB[2]-0.25*param.nz, centerLB[2]+0.25*param.nz);

        // Definition of simple mass-less particles.
        Particle3D<T,DESCRIPTOR>* particleTemplate=0;
        particleTemplate = new PointParticle3D<T,DESCRIPTOR>(0, Array<T,3>(0.,0.,0.), Array<T,3>(0.,0.,0.));

        // Functional which injects particles with predefined probability from the specified injection domain.
        std::vector<MultiBlock3D*> particleInjectionArg;
        particleInjectionArg.push_back(particles);

        integrateProcessingFunctional (
                new InjectRandomParticlesFunctional3D<T,DESCRIPTOR>(particleTemplate, param.particleProbabilityPerCell),
                injectionDomain, particleInjectionArg, 0 );

        // Definition of an absorbtion domain for the particles.
        Box3D absorbtionDomain(param.outlet);

        // Functional which absorbs the particles which reach the specified absorbtion domain.
        integrateProcessingFunctional (
                new AbsorbParticlesFunctional3D<T,DESCRIPTOR>, absorbtionDomain, particleArg, 0 );

        particles->executeInternalProcessors();
    }

    /*
     * Starting the simulation.
     */

    plb_ofstream energyFile((outputDir+"average_energy.dat").c_str());

    pcout << std::endl;
    pcout << "Starting simulation." << std::endl;
    bool checkForErrors = true;
    for (plint i = 0; i < param.maxIter; ++i) {
        if (i <= param.initialIter) {
            Array<T,3> uBoundary(param.getInletVelocity(i), 0.0, 0.0);
            setBoundaryVelocity(*lattice, param.inlet, uBoundary);
        }

        if (i % param.statIter == 0) {
             pcout << "At iteration " << i << ", t = " << i*param.dt << std::endl;
             if (i != 0) {
                 Array<T,3> force(boundaryCondition->getForceOnObject());
                 T factor = util::sqr(util::sqr(param.dx)) / util::sqr(param.dt);
                 pcout << "Force on object over fluid density: F[x] = " << force[0]*factor << ", F[y] = "
                       << force[1]*factor << ", F[z] = " << force[2]*factor << std::endl;
             }
             T avEnergy = boundaryCondition->computeAverageEnergy() * util::sqr(param.dx) / util::sqr(param.dt);
             pcout << "Average kinetic energy over fluid density: E = " << avEnergy << std::endl;
             energyFile << i*param.dt << "  " << avEnergy << std::endl;
             pcout << std::endl;
        }

        if (i % param.vtkIter == 0) {
            pcout << "Writing VTK at time t = " << i*param.dt << endl;
            writeVTK(*boundaryCondition, i);
            if (param.useParticles) {
                writeParticleVtk<T,DESCRIPTOR> (
                        *particles, createFileName(outputDir+"particles_", i, PADDING) + ".vtk",
                        param.maxNumParticlesToWrite );
            }
        }

        if (i % param.imageIter == 0) {
            pcout << "Writing PPM image at time t = " << i*param.dt << endl;
            writePPM(*boundaryCondition, i);
        }

        lattice->executeInternalProcessors();
        lattice->incrementTime();
        if (param.useParticles && i % param.particleTimeFactor == 0) {
            particles->executeInternalProcessors();
        }
        if (checkForErrors) {
            abortIfErrorsOccurred();
            checkForErrors = false;
        }
    }
    
    energyFile.close();
    delete outerBoundaryCondition;
    delete boundaryCondition;
    if (param.useParticles) {
        delete particles;
    }
    delete j;
    delete rhoBar;
    delete lattice;
} 

int main(int argc, char* argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir(outputDir);

    // The try-catch blocks catch exceptions in case an error occurs,
    // and terminate the program properly with a nice error message.

    // 1. Read command-line parameter: the input file name.
    string xmlFileName;
    try {
        global::argv(1).read(xmlFileName);
    }
    catch (PlbIOException& exception) {
        pcout << "Wrong parameters; the syntax is: " 
              << (std::string)global::argv(0) << " input-file.xml" << std::endl;
        return -1;
    }

    // 2. Read input parameters from the XML file.
    try {
        param = Param(xmlFileName);
    }
    catch (PlbIOException& exception) {
        pcout << exception.what() << std::endl;
        return -1;
    }

    // 3. Execute the main program.
    try {
        runProgram();
    }
    catch (PlbIOException& exception) {
        pcout << exception.what() << std::endl;
        return -1;
    }
}

