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
 * This is an example to show how to use the immersed boundary method in
 * Palabos. An immersed rectangular surface oscillates in the main
 * stream of a fluid. The angle of oscillation obeys a sinusoidal
 * function in time.
 */

#include "palabos3D.h"
#include "palabos3D.hh"

#include <cstdlib>
#include <cmath>

using namespace plb;

typedef double T;

#define PADDING 8

#define DESCRIPTOR descriptors::D3Q19Descriptor

plint smallEnvelopeWidth  = 1;
plint largeEnvelopeWidth  = 4; // Because of immersed walls.

std::string outputDir("./tmp/");

T pi = std::acos((T) -1);

// This structure holds all the user-defined parameters, and some
// derived values needed for the simulation.
//
// The values used in this example lead to a large mesh of
// dimensions 708 x 216 x 216. If a desktop computer is
// used to run this code, and not a parallel cluster,
// be sure to reduce the size of the grid.
struct Param {
    // User defined parameters.

    T lx, ly, lz;                                   // Dimensions of the simulation domain in physical units.
    T xSide, ySide;                                 // Lengths of the x and y sides of the immersed rectangle in physical units.
    plint resolution;                               // Number of lattice nodes along a length equal to ySide.
    T dt;                                           // Discrete time step.
    plint maxIter, statIter, outIter;               // Maximum number of iterations, and number of iterations for output.
    plint startIter;
    T nu;                                           // Fluid kinematic viscosity in physical units.
    T inletVelocity;                                // Inlet velocity in physical units.
    Array<T,3> mountPoint;                          // Mount point of the surface in the physical simulation domain.
    int ibIter;                                     // Iterations for the immersed boundary method.
    bool periodicBoundaries;                        // Periodic or freeSlip lateral boundaries.
    T ampl;                                         // Angle amplitude of the rectangle surface oscillation.
    T freq;                                         // Angle frequency of the rectangle surface oscillation.

    // Derived parameters.

    plint nx, ny, nz;                               // Grid dimensions of the simulation domain.
    T inletVelocityLB;                              // Parameters in lattice units.
    T xSideLB, ySideLB;
    Array<T,3> mountPointLB;
    T amplLB, freqLB;
    T omega;                                        // Relaxation parameter for the fluid.
    T dx;                                           // Discrete spatial step.
} param;


void setParam()
{
    // User input data (all in physical units).

    param.lx = 29.5;
    param.ly = 9.0;
    param.lz = 9.0;
    param.xSide = 2.0;
    param.ySide = 1.0;
    param.resolution = 25;
    param.dt = 0.00082;
    param.maxIter = 500000;
    param.statIter = 200;
    param.outIter = 500;
    param.startIter = 24000;
    param.inletVelocity = 1.0;
    T Re = 50.0;
    param.nu = param.inletVelocity * param.ySide / Re;
    param.mountPoint[0] = 4.5;
    param.mountPoint[1] = 4.5;
    param.mountPoint[2] = 4.5;
    param.ibIter = 4;
    param.periodicBoundaries = true;
    param.ampl = 0.25 * pi;
    param.freq = 0.1;
}


void computeLbParam()
{
    // Derived quantities.

    param.dx = param.ySide / (param.resolution - 1.0);
    if (param.periodicBoundaries) {
        param.nx = util::roundToInt(param.lx / param.dx) + 1;
        param.ny = util::roundToInt(param.ly / param.dx);
        param.nz = util::roundToInt(param.lz / param.dx);
    } else {
        param.nx = util::roundToInt(param.lx / param.dx) + 1;
        param.ny = util::roundToInt(param.ly / param.dx) + 1;
        param.nz = util::roundToInt(param.lz / param.dx) + 1;
    }
    param.xSideLB = param.xSide / param.dx;
    param.ySideLB = param.ySide / param.dx;
    param.inletVelocityLB = param.inletVelocity * param.dt / param.dx;
    param.mountPointLB = ((T) 1.0 / param.dx) * param.mountPoint;
    T nuLB = param.nu * param.dt / (param.dx * param.dx);
    param.omega = 1.0 / (DESCRIPTOR<T>::invCs2 * nuLB + 0.5);
    param.amplLB = param.ampl;
    param.freqLB = param.freq * param.dt;
}


void printParam(void)
{
    pcout << "lx = " << param.lx << std::endl;
    pcout << "ly = " << param.ly << std::endl;
    pcout << "lz = " << param.lz << std::endl;
    pcout << "resolution = " << param.resolution << std::endl;
    pcout << "xSide = " << param.xSide << std::endl;
    pcout << "ySide = " << param.ySide << std::endl;
    pcout << "ampl = " << param.ampl << std::endl;
    pcout << "freq = " << param.freq << std::endl;
    pcout << "dt = " << param.dt << std::endl;
    pcout << "maxIter = " << param.maxIter << std::endl;
    pcout << "statIter = " << param.statIter << std::endl;
    pcout << "outIter = " << param.outIter << std::endl;
    pcout << "inletVelocity = " << param.inletVelocity << std::endl;
    pcout << "nu = " << param.nu << std::endl;
    pcout << "Re = " << param.inletVelocity * param.ySide / param.nu << std::endl;
    pcout << "mountPoint = (" << param.mountPoint[0] << ", " << param.mountPoint[1] << ", " << param.mountPoint[2] << ")" << std::endl;
    pcout << "ibIter = " << param.ibIter << std::endl;
    pcout << "nx = " << param.nx << std::endl;
    pcout << "ny = " << param.ny << std::endl;
    pcout << "nz = " << param.nz << std::endl;
    pcout << "inletVelocityLB = " << param.inletVelocityLB << std::endl;
    pcout << "mountPointLB = (" << param.mountPointLB[0] << ", " << param.mountPointLB[1] << ", " << param.mountPointLB[2] << ")" << std::endl;
    pcout << "omega = " << param.omega << std::endl;
    pcout << "tau = " << 1.0 / param.omega << std::endl;
    pcout << "dx = " << param.dx << std::endl;
    pcout << "dt / dx = " << param.dt / param.dx << std::endl;
    pcout << "dt / (dx * dx) = " << param.dt / (param.dx * param.dx) << std::endl;
    pcout << "xSideLB = " << param.xSideLB << std::endl;
    pcout << "ySideLB = " << param.ySideLB << std::endl;
    pcout << "amplLB = " << param.amplLB << std::endl;
    pcout << "freqLB = " << param.freqLB << std::endl;
    pcout << std::endl;
}


// The inlet velocity starts from zero and increases in time until
// it reaches a target value specified by the user.
T getInletVelocity(T targetValue, plint iT)
{
    return targetValue * util::sinIncreasingFunction<T>(iT, param.startIter);
}


// The angle of oscillation of the immersed rectangle is given as
// a sinusoidal function of time.
T getAngle(T t)
{
    return 0.5 * pi + param.amplLB * std::sin(2.0 * pi * param.freqLB * t);
}


// The angular velocity of the immersed rectangle is the time derivative
// of the angle function.
T getAngularVelocity(T t)
{
    return 2.0 * pi * param.freqLB * param.amplLB * std::cos(2.0 * pi * param.freqLB * t);
}


// This class computes the velocity of the immersed surface at every
// time step. The immersed surface has a defined position given by 
// an angle of the form:
//   phi = phi_0 + A * sin(2 * pi * f * t)
// So its velocity is given by the cross product between its 
// angular velocity (time derivative of the angle) and its position.
// The axis of rotation is the y-axis.
class SurfaceVelocity {
public:
    SurfaceVelocity(T t_)
        : t(t_)
    { }
    Array<T,3> operator()(Array<T,3> const& pos)
    {
        T xp = pos[0] - param.mountPointLB[0];
        T zp = pos[2] - param.mountPointLB[2];
        T angularVelocity = getAngularVelocity(t);
        Array<T,3> velocity;
        velocity[0] =  angularVelocity * zp;
        velocity[1] =  0.0;
        velocity[2] = -angularVelocity * xp;
        return velocity;
    }
private:
    T t;
};


// Write an STL file of the instantaneous geometry of the
// moving immersed surface.
void writeSurface(TriangleSet<T> triangleSet, plint iT)
{
    static Array<T,3> normedAxis((T) 0, (T) 1, (T) 0);

    T initialAngle = 0.5 * pi;
    T angle = getAngle(iT) - initialAngle;
    triangleSet.translate(-param.mountPoint);
    triangleSet.rotateAtOrigin(normedAxis, angle);
    triangleSet.translate( param.mountPoint);
    triangleSet.writeBinarySTL(createFileName(outputDir+"surface_", iT, PADDING) + ".stl");
}


void writeVTK(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, Array<T,3> center, plint iT)
{
    plint cx = util::roundToInt(center[0]);
    plint cy = util::roundToInt(center[1]);
    plint cz = util::roundToInt(center[2]);

    Box3D box_x(cx - 1, cx + 1, 0, param.ny - 1, 0, param.nz - 1);
    Box3D box_y(0, param.nx - 1, cy - 1, cy + 1, 0, param.nz - 1);
    Box3D box_z(0, param.nx - 1, 0, param.ny - 1, cz - 1, cz + 1);

    VtkImageOutput3D<T> vtkOut_x(createFileName("fluid_x_", iT, PADDING), param.dx);
    std::auto_ptr<MultiTensorField3D<T,3> > vx = computeVelocity(lattice, box_x);

    vtkOut_x.writeData<3,float>(*vx, "v", param.dx / param.dt);
    std::auto_ptr<MultiScalarField3D<T> > rhox = computeDensity(lattice, box_x);
    vtkOut_x.writeData<float>(*rhox, "p", (param.dx * param.dx) / (param.dt * param.dt));

    VtkImageOutput3D<T> vtkOut_y(createFileName("fluid_y_", iT, PADDING), param.dx);
    std::auto_ptr<MultiTensorField3D<T,3> > vy = computeVelocity(lattice, box_y);

    vtkOut_y.writeData<3,float>(*vy, "v", param.dx / param.dt);
    std::auto_ptr<MultiScalarField3D<T> > rhoy = computeDensity(lattice, box_y);
    vtkOut_y.writeData<float>(*rhoy, "p", (param.dx * param.dx) / (param.dt * param.dt));

    VtkImageOutput3D<T> vtkOut_z(createFileName("fluid_z_", iT, PADDING), param.dx);
    std::auto_ptr<MultiTensorField3D<T,3> > vz = computeVelocity(lattice, box_z);

    vtkOut_z.writeData<3,float>(*vz, "v", param.dx / param.dt);
    std::auto_ptr<MultiScalarField3D<T> > rhoz = computeDensity(lattice, box_z);
    vtkOut_z.writeData<float>(*rhoz, "p", (param.dx * param.dx) / (param.dt * param.dt));
}


int main(int argc, char* argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir(outputDir);

    setParam();
    computeLbParam();
    printParam();

    // Multi blocks.

    pcout << "Generating multi-blocks." << std::endl;

    Dynamics<T,DESCRIPTOR> *dynamics = 0;
    dynamics = new IncBGKdynamics<T,DESCRIPTOR>(param.omega);
    bool incompressibleModel = true;
    pcout << "Dynamics: Incompressible BGK." << std::endl;

    MultiBlockLattice3D<T,DESCRIPTOR> *lattice = new MultiBlockLattice3D<T,DESCRIPTOR>(
            param.nx, param.ny, param.nz, dynamics->clone());
    defineDynamics(*lattice, lattice->getBoundingBox(), dynamics->clone());
    delete dynamics;
    lattice->toggleInternalStatistics(false);

    MultiScalarField3D<T> *rhoBar = generateMultiScalarField<T>((MultiBlock3D&) *lattice, largeEnvelopeWidth).release();
    rhoBar->toggleInternalStatistics(false);

    MultiTensorField3D<T,3> *j = generateMultiTensorField<T,3>((MultiBlock3D&) *lattice, largeEnvelopeWidth).release();
    j->toggleInternalStatistics(false);

    std::vector<MultiBlock3D*> rhoBarJarg;
    rhoBarJarg.push_back(lattice);
    rhoBarJarg.push_back(rhoBar);
    rhoBarJarg.push_back(j);
    integrateProcessingFunctional(
            new ExternalRhoJcollideAndStream3D<T,DESCRIPTOR>(),
            lattice->getBoundingBox(), rhoBarJarg, 0);
    integrateProcessingFunctional(
            new BoxRhoBarJfunctional3D<T,DESCRIPTOR>(),
            lattice->getBoundingBox(), rhoBarJarg, 2);

    // Boundary conditions.

    pcout << "Generating boundary conditions." << std::endl;

    if (param.periodicBoundaries) {
        lattice->periodicity().toggle(0, false);
        lattice->periodicity().toggle(1, true);
        lattice->periodicity().toggle(2, true);

        rhoBar->periodicity().toggle(0, false);
        rhoBar->periodicity().toggle(1, true);
        rhoBar->periodicity().toggle(2, true);

        j->periodicity().toggle(0, false);
        j->periodicity().toggle(1, true);
        j->periodicity().toggle(2, true);

        pcout << "Periodic lateral boundaries." << std::endl;
    } else {
        lattice->periodicity().toggleAll(false);
        rhoBar->periodicity().toggleAll(false);
        j->periodicity().toggleAll(false);

        pcout << "Free-slip lateral boundaries." << std::endl;
    }

    OnLatticeBoundaryCondition3D<T,DESCRIPTOR> *bc = createLocalBoundaryCondition3D<T,DESCRIPTOR>();

    Box3D inlet(0, 0, 0, param.ny - 1, 0, param.nz - 1);
    Box3D outlet(param.nx - 1, param.nx - 1, 0, param.ny - 1, 0, param.nz - 1);
    if (param.periodicBoundaries) {
        // Inlet velocity boundary condition.
        bc->addVelocityBoundary0N(inlet, *lattice);
        T vx = getInletVelocity(param.inletVelocityLB, 0);
        setBoundaryVelocity(*lattice, inlet, Array<T,3>(vx, (T) 0, (T) 0));

        // Outflow boundary condition.
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
                outlet, bcargs, 1);
        setBoundaryVelocity(*lattice, outlet, Array<T,3>(vx, (T) 0, (T) 0));
    } else {
        // Inlet velocity boundary condition.
        bc->setVelocityConditionOnBlockBoundaries(*lattice, inlet, boundary::dirichlet);
        T vx = getInletVelocity(param.inletVelocityLB, 0);
        setBoundaryVelocity(*lattice, inlet, Array<T,3>(vx, (T) 0, (T) 0));

        // Free-slip condition on lateral boundaries.
        Box3D lateral1(1, param.nx - 2, 0, 0, 0, param.nz - 1);
        Box3D lateral2(1, param.nx - 2, param.ny - 1, param.ny - 1, 0, param.nz - 1);
        Box3D lateral3(1, param.nx - 2, 1, param.ny - 2, 0, 0);
        Box3D lateral4(1, param.nx - 2, 1, param.ny - 2, param.nz - 1, param.nz - 1);

        bc->setVelocityConditionOnBlockBoundaries(*lattice, lateral1, boundary::freeslip);
        bc->setVelocityConditionOnBlockBoundaries(*lattice, lateral2, boundary::freeslip);
        bc->setVelocityConditionOnBlockBoundaries(*lattice, lateral3, boundary::freeslip);
        bc->setVelocityConditionOnBlockBoundaries(*lattice, lateral4, boundary::freeslip);
        setBoundaryVelocity(*lattice, lateral1, Array<T,3>(vx, (T) 0, (T) 0));
        setBoundaryVelocity(*lattice, lateral2, Array<T,3>(vx, (T) 0, (T) 0));
        setBoundaryVelocity(*lattice, lateral3, Array<T,3>(vx, (T) 0, (T) 0));
        setBoundaryVelocity(*lattice, lateral4, Array<T,3>(vx, (T) 0, (T) 0));

        // Outflow boundary condition.
        // The VirtualOutlet is a sophisticated outflow boundary condition.
        Box3D globalDomain(lattice->getBoundingBox());
        std::vector<MultiBlock3D*> bcargs;
        bcargs.push_back(lattice);
        bcargs.push_back(rhoBar);
        bcargs.push_back(j);
        T outsideDensity = 1.0;
        int bcType = 1;
        integrateProcessingFunctional(new VirtualOutlet<T,DESCRIPTOR>(outsideDensity, globalDomain, bcType),
                outlet, bcargs, 1);
        setBoundaryVelocity(*lattice, outlet, Array<T,3>(vx, (T) 0, (T) 0));
    }

    delete bc;

    // Initialization.

    T vx = getInletVelocity(param.inletVelocityLB, 0);
    initializeAtEquilibrium(*lattice, lattice->getBoundingBox(),
            (T)1.0, Array<T,3>(vx, (T) 0, (T) 0));
    applyProcessingFunctional(
            new BoxRhoBarJfunctional3D<T,DESCRIPTOR>(),
            lattice->getBoundingBox(), rhoBarJarg);
    T energy = computeAverageEnergy(*lattice) * (param.dx * param.dx) / (param.dt * param.dt);
    pcout << "Initial average kinetic energy: " << energy << std::endl;
    pcout << std::endl;

    // Immersed rectangle surface.

    pcout << "Creating the immersed rectangle surface." << std::endl;

    // The immersed boundary method needs a set of vertices and a set of areas
    // that correspond to each vertex. These vertices and areas describe the
    // time dependent geometry of the surface at each time step.
    std::vector<Array<T,3> > vertices;
    std::vector<T> areas;

    // The initial geometry of the immersed surface can be created analyticaly or, alternatively,
    // it can be loaded by reading a user provided STL file.
    TriangleSet<T> rectangleTriangleSet = constructRectangle<T>(param.xSideLB, param.ySideLB,
            util::roundToInt(2.0 * param.xSideLB), util::roundToInt(2.0 * param.ySideLB));
    Array<T,3> mount = Array<T,3>((T) 0, (T) -0.5 * param.ySideLB, (T) 0) + param.mountPointLB;
    rectangleTriangleSet.translate(mount);

    DEFscaledMesh<T> *rectangleDef = new DEFscaledMesh<T>(rectangleTriangleSet, 0, 0, 0, Dot3D(0, 0, 0));
    pcout << "The rectangle has " << rectangleDef->getMesh().getNumVertices() << " vertices and " <<
        rectangleDef->getMesh().getNumTriangles() << " triangles." << std::endl;
    for (pluint iVertex = 0; iVertex < (pluint) rectangleDef->getMesh().getNumVertices(); iVertex++) {
        vertices.push_back(rectangleDef->getMesh().getVertex(iVertex));
        areas.push_back(rectangleDef->getMesh().computeVertexArea(iVertex));
    }
    delete rectangleDef;

    // The object "rectangleTriangleSet" from now on will be used  only for visualization purposes,
    // so it is changed back to physical units.
    rectangleTriangleSet.scale(param.dx);

    // The next container block is necessary for the immersed-wall algorithm.
    MultiContainerBlock3D container(*rhoBar);

    pcout << std::endl;

    plb_ofstream energyFile;
    std::string energyFname = outputDir + "averageEnergy.dat";
    energyFile.open(energyFname.c_str());

    pcout << "Starting simulation." << std::endl;

    for (plint iT = 0; iT < param.maxIter; iT++) {
        // Set the time dependent inlet boundary velocity.
        if (iT <= param.startIter) {
            T vx = getInletVelocity(param.inletVelocityLB, iT);
            setBoundaryVelocity(*lattice, inlet, Array<T,3>(vx, (T) 0, (T) 0));
        }

        if ((iT % param.statIter == 0 && iT != 0) || iT == param.maxIter - 1) {
            pcout << "At iteration " << iT << ", t = " << iT * param.dt << std::endl;
            T energy = computeAverageEnergy(*lattice) * (param.dx * param.dx) / (param.dt * param.dt);
            pcout << "Average kinetic energy: " << energy << std::endl;
            energyFile << iT * param.dt << " " << energy << std::endl;
            pcout << std::endl;
        }

        if (iT % param.outIter == 0 || iT == param.maxIter - 1) {
            pcout << "Output to disk at iteration: " << iT << std::endl;
            writeVTK(*lattice, param.mountPointLB, iT);
            writeSurface(rectangleTriangleSet, iT);
            pcout << std::endl;
        }

        lattice->executeInternalProcessors(); // Execute all processors and communicate appropriately.

        // Immersed walls algorithm.

        T timeLB = iT + 1.0;

        // New position of the immersed surface at the next time iteration.
        // The angle of rotation is given by a function of the form:
        //   phi = phi_0 + A * sin(2 * pi * f * t)
        // The rotation axis is the y-axis.
        for (plint iVertex = 0; iVertex < (plint) vertices.size(); iVertex++) {
            T z = vertices[iVertex][2] - param.mountPointLB[2];
            T x = vertices[iVertex][0] - param.mountPointLB[0];
            T r = std::sqrt(z * z + x * x);
            T phi = getAngle(timeLB);
            vertices[iVertex][2] = r * std::cos(phi) + param.mountPointLB[2];
            vertices[iVertex][0] = r * std::sin(phi) + param.mountPointLB[0];
        }

        // Instantiate the immersed wall data and performed the immersed boundary iterations.
        instantiateImmersedWallData(vertices, areas, container);
        for (int i = 0; i < param.ibIter; i++) {
            inamuroIteration(SurfaceVelocity(timeLB),
                    *rhoBar, *j, container, (T) 1.0 / param.omega, incompressibleModel);
        }
    }

    energyFile.close();

    delete j;
    delete rhoBar;
    delete lattice;

    return 0;
}

