
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
 * External flow around many 3D objects (moving and static).
 * This example demonstrates many features of Palabos:
 * Configuration with an input XML file.
 * Loading geometries from STL files.
 * Automatically refining 3D surface meshes.
 * Using the Immersed Boundary Method simultaneously for moving and static objects.
 * Imposing sophisticated outflow boundary conditions.
 * Using sponge zones.
 * Imposing time dependent inlet boundary conditions.
 * Computing the force on objects.
 * Computing the torque on moving objects.
 * Checkpointing (saving the state of the simulation and restarting).
 * Parallel I/O.
 * */

/* 
 * Description and usage:
 *
 * This code solves the flow around a rotating cube and a static sphere.
 *
 * The code needs to be compiled only once, and then subsequent runs are configured
 * completely through the input XML file. The user can provide his own geometries
 * as STL files.
 *
 * To run the code with, say 4, processes use:
 *
 * mpirun -np 4 ./generalExternalFlow generalExternalFlow.xml
 *
 * The input XML file contains many comments on how to configure the simulation.
 * To stop and save the state of the simulation, just create a file called "abort"
 * in the working directory. To restart the simulation use:
 *
 * mpirun -np 4 ./generalExternalFlow generalExternalFlow.xml continue.xml
 *
 * The file "continue.xml" is generated automatically.
 *
 */


#include "palabos3D.h"
#include "palabos3D.hh"

#include <cstdlib>
#include <cmath>
#include <memory>

using namespace plb;

typedef double T;

#define DESCRIPTOR descriptors::D3Q19Descriptor

std::string outDir("./tmp/");

struct SimulationParameters {
    /*
     * Parameters set by the user.
     * All user input variables and all data in external input files must be in the same system of units.
     */

    std::vector<T> xDomain;                         // Extent in the x-direction of the physical simulation domain.
    std::vector<T> yDomain;                         // Extent in the y-direction of the physical simulation domain.
    std::vector<T> zDomain;                         // Extent in the z-direction of the physical simulation domain.

    std::vector<std::string> movingSurfaceFileNames;    // Files with the moving immersed surface geometries.
    std::vector<std::string> staticSurfaceFileNames;    // Files with the static immersed surface geometries.
    plint numMovingSurfaces;                            // Number of moving immersed surfaces.
    plint numStaticSurfaces;                            // Number of static immersed surfaces.
    plint numSurfaces;                                  // Number of all immersed surfaces.
    bool refineSurfaceMeshes;                           // Immersed surfaces are represented as triangle surface meshes.
                                                        // Because of the immersed boundary method used, the size of
                                                        // the triangles must approximately be the same as the lattice cell size.
                                                        // The value of this parameter provides the option to refine the
                                                        // triangular representation of the immersed surfaces.
    T targetMaxEdgeLength;                              // If surface mesh refinement has been selected, this parameter
                                                        // defines the target value of the maximum edge length of the new
                                                        // refined surface mesh.
    plint maxNumSurfaceRefinements;                     // This is the maximum number of surface refinements to be performed,
                                                        // whether or not the target maximum triangle edge length has been
                                                        // achieved.

    std::vector<ConnectedTriangleSet<T> > allSurfaces;

    int flowDirection;                              // Direction of the flow (0 -> positive x, 1 -> positive y, 2 -> positive z).

    T characteristicLength;                         // Length to define the resolution and the Reynolds number.
    plint resolution;                               // Total number of lattice nodes in the characteristic length.

    T dt;                                           // Discrete time step.

    plint maxIter;                                  // Maximum number of iterations.
    plint statIter;                                 // Number of iterations for terminal output.
    plint outIter;                                  // Number of iterations for disk output.
    plint cpIter;                                   // Number of iterations for checkpointing.
    plint abIter;                                   // Number of iterations for checking for user-driven program abortion.
    plint startIter;                                // Number of initial iterations for smoothly increasing the inlet velocity.

    int ibIter;                                     // Iterations for the immersed boundary method.

    std::string abortFileName;                      // File for signaling program abortion.
    std::string xmlContinueFileName;                // XML file for restarting.
    std::string baseFileName;                       // Basename of the checkpoint files.

    T rho;                                          // Fluid density in physical units
    T nu;                                           // Fluid kinematic viscosity in physical units.
    T ambientPressure;                              // Absolute stagnation pressure in physical units.
    T cSmago;                                       // Smagorinsky parameter.
    Array<T,3> inletVelocity;                       // Inlet velocity vector in physical units.
    std::vector<Array<T,3> > angularVelocities;     // Angular velocity vectors for all moving surfaces in physical units.
    std::vector<Array<T,3> > rotationAxisPoints;    // Points on the rotation axis of every moving surface in physical units.

    bool lateralPeriodic;                           // Use periodic lateral boundaries, or free-slip ones?
    int outflowBcType;                              // Type of the outflow boundary condition.
                                                    // If 0, then a constant pressure is imposed at the outlet.
                                                    // If 1, then FluidPressureOutlet3D is used.
                                                    // If 2, then VirtualOutlet is used.
    Array<T,6> spongeWidths;                        // Sponge zone widths (0.0 for no sponge zone for the corresponding
                                                    // lattice boundary).
    bool useSmagorinskySponges;                     // Type of the sponge zones: viscosity or Smagorinsky.
    T targetCSmago;                                 // Smagorinsky parameter at the end of the sponges, when Smagorinsky
                                                    // sponge zones are used.

    bool useParallelIO;                             // For a desktop PC this should be "false", for a cluster "true".

    Precision precision;                            // Precision for geometric operations.

    std::vector<Array<T,3> > torqueAxisPoints;      // Axes with respect to which the torques on moving objects is computed.
    std::vector<Array<T,3> > torqueAxisDirections;

    bool outputInDomain;                            // Save data on disk in a volume domain or not?
    Cuboid<T> outputCuboid;                         // Volume domain for disk output.

    bool outputOnSlices;                            // Save data on disk on a set of slices or not?
    std::vector<T> xPositions;                      // Positions of the x-slices for output.
    std::vector<T> xyRange;                         // y range of the x-slices.
    std::vector<T> xzRange;                         // z range of the x-slices.
    std::vector<T> yPositions;                      // Positions of the y-slices for output.
    std::vector<T> yzRange;                         // z range of the y-slices.
    std::vector<T> yxRange;                         // x range of the y-slices.
    std::vector<T> zPositions;                      // Positions of the z-slices for output.
    std::vector<T> zxRange;                         // x range of the z-slices.
    std::vector<T> zyRange;                         // y range of the z-slices.

    /*
     * Parameters NOT set by the user.
     */

    T lx, ly, lz;                                   // Dimensions of the physical simulation domain.
    plint nx, ny, nz;                               // Grid dimensions of the simulation domain.
    T rho_LB;
    T ambientPressure_LB;
    Array<T,3> inletVelocity_LB;
    std::vector<Array<T,3> > angularVelocities_LB;
    std::vector<Array<T,3> > rotationAxisPoints_LB;
    std::vector<Array<T,3> > rotationAxisUnitVectors;
    std::vector<plint> startIds;
    std::vector<Array<T,3> > vertices;
    std::vector<T> areas;
    std::vector<int> flags;
    plint nextIter;
    Array<plint,6> numSpongeCells;                  // Number of sponge zone cells per lattice boundary.
    T omega;                                        // Relaxation parameter for the fluid.
    T dx;                                           // Discrete spatial step.
    Array<T,3> physicalLocation;                    // Location of the physical domain.
    Box3D inlet, outlet;
    Box3D lateral1, lateral2;
    Box3D lateral3, lateral4;
    plint smallEnvelopeWidth;                       // Standard width.
    plint largeEnvelopeWidth;                       // For velocity because of immersed walls.
    bool saveDynamicContent;
    plint fileNamePadding;
    bool incompressibleModel;

    std::vector<Array<T,3> > torqueAxisPoints_LB;

    Box3D outputDomain;                             // Domains for disk output.
    std::vector<Box3D> xSlices;
    std::vector<Box3D> ySlices;
    std::vector<Box3D> zSlices;
};

T toPhys(T lbVal, plint direction, T dx, Array<T,3> const& location)
{
    PLB_ASSERT(direction >= 0 && direction <= 2);
    return (lbVal * dx + location[direction]);
}

Array<T,3> toPhys(Array<T,3> const& lbVal, T dx, Array<T,3> const& location)
{
    return (lbVal * dx + location);
}

T toLB(T physVal, plint direction, T dx, Array<T,3> const& location)
{
    PLB_ASSERT(direction >= 0 && direction <= 2);
    return (physVal - location[direction]) / dx;
}

Array<T,3> toLB(Array<T,3> const& physVal, T dx, Array<T,3> const& location)
{
    return (physVal - location) / dx;
}

void readUserDefinedSimulationParameters(std::string xmlInputFileName, SimulationParameters& param)
{
    XMLreader document(xmlInputFileName);

    document["geometry"]["simulationDomain"]["x"].read(param.xDomain);
    PLB_ASSERT(param.xDomain.size() == 2 && param.xDomain[1] > param.xDomain[0]);
    document["geometry"]["simulationDomain"]["y"].read(param.yDomain);
    PLB_ASSERT(param.yDomain.size() == 2 && param.yDomain[1] > param.yDomain[0]);
    document["geometry"]["simulationDomain"]["z"].read(param.zDomain);
    PLB_ASSERT(param.zDomain.size() == 2 && param.zDomain[1] > param.zDomain[0]);
    document["geometry"]["movingSurfaceFileNames"].read(param.movingSurfaceFileNames);
    param.numMovingSurfaces = param.movingSurfaceFileNames.size();
    document["geometry"]["staticSurfaceFileNames"].read(param.staticSurfaceFileNames);
    param.numStaticSurfaces = param.staticSurfaceFileNames.size();
    param.numSurfaces = param.numStaticSurfaces + param.numMovingSurfaces;
    PLB_ASSERT(param.numSurfaces > 0);
    document["geometry"]["flowDirection"].read(param.flowDirection);
    PLB_ASSERT(param.flowDirection == 0 || param.flowDirection == 1 || param.flowDirection == 2);

    document["numerics"]["refineSurfaceMeshes"].read(param.refineSurfaceMeshes);
    if (param.refineSurfaceMeshes) {
        document["numerics"]["targetMaxEdgeLength"].read(param.targetMaxEdgeLength);
        document["numerics"]["maxNumSurfaceRefinements"].read(param.maxNumSurfaceRefinements);
        PLB_ASSERT(param.maxNumSurfaceRefinements > 0);
    }
    std::string precision;
    document["numerics"]["precision"].read(precision);
    PLB_ASSERT(precision == "FLT" || precision == "DBL" || precision == "LDBL" || precision == "INF");
    if (precision == "FLT") {
        param.precision = FLT;
    } else if (precision == "DBL") {
        param.precision = DBL;
    } else if (precision == "LDBL") {
        param.precision = LDBL;
    } else {
        param.precision = INF;
    }

    document["numerics"]["characteristicLength"].read(param.characteristicLength);
    document["numerics"]["resolution"].read(param.resolution);
    PLB_ASSERT(param.resolution > 2);
    document["numerics"]["dt"].read(param.dt);
    document["numerics"]["maxIter"].read(param.maxIter);
    PLB_ASSERT(param.maxIter > 0);
    document["numerics"]["ibIter"].read(param.ibIter);
    PLB_ASSERT(param.ibIter >= 0);

    document["numerics"]["startIter"].read(param.startIter);
    PLB_ASSERT(param.startIter >= 0);
    document["numerics"]["cSmago"].read(param.cSmago);
    document["numerics"]["inletVelocity"].read<T,3>(param.inletVelocity);
    document["numerics"]["lateralPeriodic"].read(param.lateralPeriodic);

    {
        std::vector<T> x, y, z;
        document["numerics"]["angularVelocities"]["x"].read(x);
        document["numerics"]["angularVelocities"]["y"].read(y);
        document["numerics"]["angularVelocities"]["z"].read(z);
        PLB_ASSERT(x.size() == y.size() && y.size() == z.size());
        PLB_ASSERT(x.size() == param.movingSurfaceFileNames.size());
        plint sz = x.size();
        param.angularVelocities.resize(sz);
        for (plint iSurface = 0; iSurface < sz; iSurface++) {
            param.angularVelocities[iSurface][0] = x[iSurface];
            param.angularVelocities[iSurface][1] = y[iSurface];
            param.angularVelocities[iSurface][2] = z[iSurface];
        }
    }
    {
        std::vector<T> x, y, z;
        document["numerics"]["rotationAxisPoints"]["x"].read(x);
        document["numerics"]["rotationAxisPoints"]["y"].read(y);
        document["numerics"]["rotationAxisPoints"]["z"].read(z);
        PLB_ASSERT(x.size() == y.size() && y.size() == z.size());
        PLB_ASSERT(x.size() == param.movingSurfaceFileNames.size());
        plint sz = x.size();
        param.rotationAxisPoints.resize(sz);
        for (plint iSurface = 0; iSurface < sz; iSurface++) {
            param.rotationAxisPoints[iSurface][0] = x[iSurface];
            param.rotationAxisPoints[iSurface][1] = y[iSurface];
            param.rotationAxisPoints[iSurface][2] = z[iSurface];
        }
    }

    std::vector<T> sWidths;
    document["numerics"]["spongeZones"]["xWidths"].read(sWidths);
    PLB_ASSERT(sWidths.size() == 2);
    param.spongeWidths[0] = sWidths[0];
    param.spongeWidths[1] = sWidths[1];
    sWidths.clear();
    document["numerics"]["spongeZones"]["yWidths"].read(sWidths);
    PLB_ASSERT(sWidths.size() == 2);
    param.spongeWidths[2] = sWidths[0];
    param.spongeWidths[3] = sWidths[1];
    sWidths.clear();
    document["numerics"]["spongeZones"]["zWidths"].read(sWidths);
    PLB_ASSERT(sWidths.size() == 2);
    param.spongeWidths[4] = sWidths[0];
    param.spongeWidths[5] = sWidths[1];
    sWidths.clear();
    document["numerics"]["spongeZones"]["useSmagorinskySponges"].read(param.useSmagorinskySponges);
    if (param.useSmagorinskySponges) {
        document["numerics"]["spongeZones"]["targetCSmago"].read(param.targetCSmago);
    }

    document["numerics"]["outflowBcType"].read(param.outflowBcType);
    PLB_ASSERT(param.outflowBcType >= 0 && param.outflowBcType <= 2);
    document["numerics"]["abortFileName"].read(param.abortFileName);
    document["numerics"]["xmlContinueFileName"].read(param.xmlContinueFileName);
    document["numerics"]["baseFileName"].read(param.baseFileName);
    document["numerics"]["useParallelIO"].read(param.useParallelIO);

    document["fluid"]["rho"].read(param.rho);
    document["fluid"]["nu"].read(param.nu);
    document["fluid"]["ambientPressure"].read(param.ambientPressure);

    document["output"]["statIter"].read(param.statIter);
    PLB_ASSERT(param.statIter > 0);
    document["output"]["outIter"].read(param.outIter);
    PLB_ASSERT(param.outIter > 0);
    document["output"]["cpIter"].read(param.cpIter);
    document["output"]["abIter"].read(param.abIter);
    PLB_ASSERT(param.abIter > 0);

    {
        std::vector<T> x, y, z;
        document["output"]["torques"]["axisPoints"]["x"].read(x);
        document["output"]["torques"]["axisPoints"]["y"].read(y);
        document["output"]["torques"]["axisPoints"]["z"].read(z);
        PLB_ASSERT(x.size() == y.size() && y.size() == z.size());
        PLB_ASSERT(x.size() == param.movingSurfaceFileNames.size());
        plint sz = x.size();
        param.torqueAxisPoints.resize(sz);
        for (plint iSurface = 0; iSurface < sz; iSurface++) {
            param.torqueAxisPoints[iSurface][0] = x[iSurface];
            param.torqueAxisPoints[iSurface][1] = y[iSurface];
            param.torqueAxisPoints[iSurface][2] = z[iSurface];
        }
    }
    {
        std::vector<T> x, y, z;
        document["output"]["torques"]["axisDirections"]["x"].read(x);
        document["output"]["torques"]["axisDirections"]["y"].read(y);
        document["output"]["torques"]["axisDirections"]["z"].read(z);
        PLB_ASSERT(x.size() == y.size() && y.size() == z.size());
        PLB_ASSERT(x.size() == param.movingSurfaceFileNames.size());
        plint sz = x.size();
        param.torqueAxisDirections.resize(sz);
        for (plint iSurface = 0; iSurface < sz; iSurface++) {
            param.torqueAxisDirections[iSurface][0] = x[iSurface];
            param.torqueAxisDirections[iSurface][1] = y[iSurface];
            param.torqueAxisDirections[iSurface][2] = z[iSurface];

            T axisNorm = norm(param.torqueAxisDirections[iSurface]);
            PLB_ASSERT(!util::isZero(axisNorm));
            param.torqueAxisDirections[iSurface] /= axisNorm;
        }
    }

    document["output"]["outputInDomain"].read(param.outputInDomain);
    if (param.outputInDomain) {
        std::vector<T> x, y, z;
        document["output"]["outputDomain"]["x"].read(x);
        PLB_ASSERT(x.size() == 2 && x[1] > x[0]);
        document["output"]["outputDomain"]["y"].read(y);
        PLB_ASSERT(y.size() == 2 && y[1] > y[0]);
        document["output"]["outputDomain"]["z"].read(z);
        PLB_ASSERT(z.size() == 2 && z[1] > z[0]);
        param.outputCuboid.lowerLeftCorner[0] = x[0];
        param.outputCuboid.lowerLeftCorner[1] = y[0];
        param.outputCuboid.lowerLeftCorner[2] = z[0];
        param.outputCuboid.upperRightCorner[0] = x[1];
        param.outputCuboid.upperRightCorner[1] = y[1];
        param.outputCuboid.upperRightCorner[2] = z[1];
    }

    document["output"]["outputOnSlices"].read(param.outputOnSlices);
    if (param.outputOnSlices) {
        document["output"]["outputSlices"]["xSlices"]["xPositions"].read(param.xPositions);
        document["output"]["outputSlices"]["xSlices"]["yRange"].read(param.xyRange);
        PLB_ASSERT(param.xyRange.size() == 2 && param.xyRange[1] > param.xyRange[0]);
        document["output"]["outputSlices"]["xSlices"]["zRange"].read(param.xzRange);
        PLB_ASSERT(param.xzRange.size() == 2 && param.xzRange[1] > param.xzRange[0]);

        document["output"]["outputSlices"]["ySlices"]["yPositions"].read(param.yPositions);
        document["output"]["outputSlices"]["ySlices"]["zRange"].read(param.yzRange);
        PLB_ASSERT(param.yzRange.size() == 2 && param.yzRange[1] > param.yzRange[0]);
        document["output"]["outputSlices"]["ySlices"]["xRange"].read(param.yxRange);
        PLB_ASSERT(param.yxRange.size() == 2 && param.yxRange[1] > param.yxRange[0]);

        document["output"]["outputSlices"]["zSlices"]["zPositions"].read(param.zPositions);
        document["output"]["outputSlices"]["zSlices"]["xRange"].read(param.zxRange);
        PLB_ASSERT(param.zxRange.size() == 2 && param.zxRange[1] > param.zxRange[0]);
        document["output"]["outputSlices"]["zSlices"]["yRange"].read(param.zyRange);
        PLB_ASSERT(param.zyRange.size() == 2 && param.zyRange[1] > param.zyRange[0]);
    }
}

void computeOutputDomain(SimulationParameters& param)
{
    if (!param.outputInDomain) {
        return;
    }

    Array<T,3> llc = param.outputCuboid.lowerLeftCorner;
    Array<T,3> urc = param.outputCuboid.upperRightCorner;

    plint x0 = util::roundToInt(toLB(llc[0], 0, param.dx, param.physicalLocation));
    plint y0 = util::roundToInt(toLB(llc[1], 1, param.dx, param.physicalLocation));
    plint z0 = util::roundToInt(toLB(llc[2], 2, param.dx, param.physicalLocation));

    plint x1 = util::roundToInt(toLB(urc[0], 0, param.dx, param.physicalLocation));
    plint y1 = util::roundToInt(toLB(urc[1], 1, param.dx, param.physicalLocation));
    plint z1 = util::roundToInt(toLB(urc[2], 2, param.dx, param.physicalLocation));

    PLB_ASSERT(x1 >= x0 && y1 >= y0 && z1 >= z0);

    param.outputDomain = Box3D(x0, x1, y0, y1, z0, z1);
}

void computeOutputSlices(SimulationParameters& param)
{
    if (!param.outputOnSlices) {
        return;
    }

    {
        param.xSlices.clear();

        plint y0 = util::roundToInt(toLB(param.xyRange[0], 1, param.dx, param.physicalLocation));
        plint y1 = util::roundToInt(toLB(param.xyRange[1], 1, param.dx, param.physicalLocation));
        plint z0 = util::roundToInt(toLB(param.xzRange[0], 2, param.dx, param.physicalLocation));
        plint z1 = util::roundToInt(toLB(param.xzRange[1], 2, param.dx, param.physicalLocation));
        PLB_ASSERT(y1 >= y0 && z1 >= z0);

        for (plint i = 0; i < (plint) param.xPositions.size(); i++) {
            plint xPos = util::roundToInt(toLB(param.xPositions[i], 0, param.dx, param.physicalLocation));
            plint x0 = xPos - 1;
            plint x1 = xPos + 1;
            param.xSlices.push_back(Box3D(x0, x1, y0, y1, z0, z1));
        }
    }

    {
        param.ySlices.clear();

        plint z0 = util::roundToInt(toLB(param.yzRange[0], 2, param.dx, param.physicalLocation));
        plint z1 = util::roundToInt(toLB(param.yzRange[1], 2, param.dx, param.physicalLocation));
        plint x0 = util::roundToInt(toLB(param.yxRange[0], 0, param.dx, param.physicalLocation));
        plint x1 = util::roundToInt(toLB(param.yxRange[1], 0, param.dx, param.physicalLocation));
        PLB_ASSERT(z1 >= z0 && x1 >= x0);

        for (plint i = 0; i < (plint) param.yPositions.size(); i++) {
            plint yPos = util::roundToInt(toLB(param.yPositions[i], 1, param.dx, param.physicalLocation));
            plint y0 = yPos - 1;
            plint y1 = yPos + 1;
            param.ySlices.push_back(Box3D(x0, x1, y0, y1, z0, z1));
        }
    }

    {
        param.zSlices.clear();

        plint x0 = util::roundToInt(toLB(param.zxRange[0], 0, param.dx, param.physicalLocation));
        plint x1 = util::roundToInt(toLB(param.zxRange[1], 0, param.dx, param.physicalLocation));
        plint y0 = util::roundToInt(toLB(param.zyRange[0], 1, param.dx, param.physicalLocation));
        plint y1 = util::roundToInt(toLB(param.zyRange[1], 1, param.dx, param.physicalLocation));
        PLB_ASSERT(x1 >= x0 && y1 >= y0);

        for (plint i = 0; i < (plint) param.zPositions.size(); i++) {
            plint zPos = util::roundToInt(toLB(param.zPositions[i], 2, param.dx, param.physicalLocation));
            plint z0 = zPos - 1;
            plint z1 = zPos + 1;
            param.zSlices.push_back(Box3D(x0, x1, y0, y1, z0, z1));
        }
    }
}

void defineOuterDomain(SimulationParameters& param)
{
    if (param.flowDirection == 0) {
        param.inlet       = Box3D(0,          0,          0,          param.ny-1, 0,          param.nz-1);
        param.outlet      = Box3D(param.nx-1, param.nx-1, 0,          param.ny-1, 0,          param.nz-1);
        param.lateral1    = Box3D(1,          param.nx-2, 0,          0,          0,          param.nz-1);
        param.lateral2    = Box3D(1,          param.nx-2, param.ny-1, param.ny-1, 0,          param.nz-1);
        param.lateral3    = Box3D(1,          param.nx-2, 1,          param.ny-2, 0,          0);
        param.lateral4    = Box3D(1,          param.nx-2, 1,          param.ny-2, param.nz-1, param.nz-1);
    } else if (param.flowDirection == 1) {
        param.inlet       = Box3D(0,          param.nx-1, 0,          0,          0,          param.nz-1);
        param.outlet      = Box3D(0,          param.nx-1, param.ny-1, param.ny-1, 0,          param.nz-1);
        param.lateral1    = Box3D(0,          param.nx-1, 1,          param.ny-2, 0,          0);
        param.lateral2    = Box3D(0,          param.nx-1, 1,          param.ny-2, param.nz-1, param.nz-1);
        param.lateral3    = Box3D(0,          0,          1,          param.ny-2, 1,          param.nz-2);
        param.lateral4    = Box3D(param.nx-1, param.nx-1, 1,          param.ny-2, 1,          param.nz-2);
    } else {
        param.inlet       = Box3D(0,          param.nx-1, 0,          param.ny-1, 0,          0);
        param.outlet      = Box3D(0,          param.nx-1, 0,          param.ny-1, param.nz-1, param.nz-1);
        param.lateral1    = Box3D(0,          0,          0,          param.ny-1, 1,          param.nz-2);
        param.lateral2    = Box3D(param.nx-1, param.nx-1, 0,          param.ny-1, 1,          param.nz-2);
        param.lateral3    = Box3D(1,          param.nx-2, 0,          0,          1,          param.nz-2);
        param.lateral4    = Box3D(1,          param.nx-2, param.ny-1, param.ny-1, 1,          param.nz-2);
    }
}

void calculateDerivedSimulationParameters(SimulationParameters& param)
{
    // Derived quantities.

    param.smallEnvelopeWidth = 1;
    param.largeEnvelopeWidth = 4;
    param.saveDynamicContent = true;
    param.fileNamePadding = 8;

    param.lx = param.xDomain[1] - param.xDomain[0];
    param.ly = param.yDomain[1] - param.yDomain[0];
    param.lz = param.zDomain[1] - param.zDomain[0];
    param.physicalLocation = Array<T,3>(param.xDomain[0], param.yDomain[0], param.zDomain[0]);

    param.dx = param.characteristicLength / (param.resolution - 1.0);

    if (param.lateralPeriodic) {
        if (param.flowDirection == 0) {
            param.nx = util::roundToInt(param.lx / param.dx) + 1;
            param.ny = util::roundToInt(param.ly / param.dx);
            param.nz = util::roundToInt(param.lz / param.dx);
        } else if (param.flowDirection == 1) {
            param.nx = util::roundToInt(param.lx / param.dx);
            param.ny = util::roundToInt(param.ly / param.dx) + 1;
            param.nz = util::roundToInt(param.lz / param.dx);
        } else {
            param.nx = util::roundToInt(param.lx / param.dx);
            param.ny = util::roundToInt(param.ly / param.dx);
            param.nz = util::roundToInt(param.lz / param.dx) + 1;
        }
    } else {
        param.nx = util::roundToInt(param.lx / param.dx) + 1;
        param.ny = util::roundToInt(param.ly / param.dx) + 1;
        param.nz = util::roundToInt(param.lz / param.dx) + 1;
    }

    param.rho_LB = 1.0;
    param.ambientPressure_LB = (1.0/param.rho) * (param.dt*param.dt/(param.dx*param.dx)) * param.ambientPressure;
    param.inletVelocity_LB = param.inletVelocity * (param.dt / param.dx);

    param.angularVelocities_LB.resize(param.angularVelocities.size());
    for (plint iSurface = 0; iSurface < (plint) param.angularVelocities_LB.size(); iSurface++) {
        param.angularVelocities_LB[iSurface] = param.angularVelocities[iSurface] * param.dt;
    }
    param.rotationAxisPoints_LB.resize(param.rotationAxisPoints.size());
    for (plint iSurface = 0; iSurface < (plint) param.rotationAxisPoints_LB.size(); iSurface++) {
        param.rotationAxisPoints_LB[iSurface] = toLB(param.rotationAxisPoints[iSurface], param.dx, param.physicalLocation);
    }
    param.rotationAxisUnitVectors.resize(param.angularVelocities.size());
    for (plint iSurface = 0; iSurface < (plint) param.rotationAxisUnitVectors.size(); iSurface++) {
        T angVelNorm = norm(param.angularVelocities[iSurface]);
        if (util::isZero(angVelNorm)) {
            param.rotationAxisUnitVectors[iSurface] = Array<T,3>((T) 1, (T) 0, (T) 0);
        } else {
            param.rotationAxisUnitVectors[iSurface] = param.angularVelocities[iSurface] / angVelNorm;
        }
    }

    if (param.refineSurfaceMeshes) {
        if (param.targetMaxEdgeLength < 0.0) {
            param.targetMaxEdgeLength = param.dx;
        }
    }

    for (int iSponge = 0; iSponge < 6; iSponge++) {
        param.numSpongeCells[iSponge] = util::roundToInt(param.spongeWidths[iSponge] / param.dx);
    }

    T nu_LB = param.nu * param.dt / (param.dx * param.dx);
    param.omega = 1.0 / (DESCRIPTOR<T>::invCs2 * nu_LB + 0.5);

    param.torqueAxisPoints_LB.resize(param.torqueAxisPoints.size());
    for (plint iSurface = 0; iSurface < (plint) param.torqueAxisPoints_LB.size(); iSurface++) {
        param.torqueAxisPoints_LB[iSurface] = toLB(param.torqueAxisPoints[iSurface], param.dx, param.physicalLocation);
    }

    computeOutputDomain(param);
    computeOutputSlices(param);
    defineOuterDomain(param);
}

void printSimulationParameters(SimulationParameters const& param)
{
    pcout << "xDomain = [" << param.xDomain[0] << ", " << param.xDomain[1] << "]" << std::endl;
    pcout << "yDomain = [" << param.yDomain[0] << ", " << param.yDomain[1] << "]" << std::endl;
    pcout << "zDomain = [" << param.zDomain[0] << ", " << param.zDomain[1] << "]" << std::endl;

    for (plint i = 0; i < (plint) param.movingSurfaceFileNames.size(); i++) {
        pcout << "movingSurfaceFileNames[" << i << "] = " << param.movingSurfaceFileNames[i] << std::endl;
    }
    for (plint i = 0; i < (plint) param.staticSurfaceFileNames.size(); i++) {
        pcout << "staticSurfaceFileNames[" << i << "] = " << param.staticSurfaceFileNames[i] << std::endl;
    }
    pcout << "refineSurfaceMeshes = " << (param.refineSurfaceMeshes ? "true" : "false") << std::endl;
    if (param.refineSurfaceMeshes) {
        pcout << "targetMaxEdgeLength = " << param.targetMaxEdgeLength << std::endl;
        pcout << "maxNumSurfaceRefinements = " << param.maxNumSurfaceRefinements << std::endl;
    }
    pcout << "flowDirection = " << param.flowDirection << std::endl;
    pcout << "characteristicLength = " << param.characteristicLength << std::endl;
    pcout << "resolution = " << param.resolution << std::endl;
    pcout << "dt = " << param.dt << std::endl;
    pcout << "maxIter = " << param.maxIter << std::endl;
    pcout << "statIter = " << param.statIter << std::endl;
    pcout << "outIter = " << param.outIter << std::endl;
    pcout << "cpIter = " << param.cpIter << std::endl;
    pcout << "abIter = " << param.abIter << std::endl;
    pcout << "startIter = " << param.startIter << std::endl;
    pcout << "ibIter = " << param.ibIter << std::endl;
    pcout << "abortFileName = " << param.abortFileName << std::endl;
    pcout << "xmlContinueFileName = " << param.xmlContinueFileName << std::endl;
    pcout << "baseFileName = " << param.baseFileName << std::endl;
    pcout << "rho = " << param.rho << std::endl;
    pcout << "nu = " << param.nu << std::endl;
    pcout << "cSmago = " << param.cSmago << std::endl;
    pcout << "inletVelocity = [" << param.inletVelocity[0] << ", "
                                 << param.inletVelocity[1] << ", "
                                 << param.inletVelocity[2] << "]" << std::endl;
    pcout << "lateralPeriodic = " << (param.lateralPeriodic ? "true" : "false") << std::endl;
    for (plint iSurface = 0; iSurface < (plint) param.angularVelocities.size(); iSurface++) {
        pcout << "angularVelocities[" << iSurface << "] = [" << param.angularVelocities[iSurface][0] << ", "
                                                             << param.angularVelocities[iSurface][1] << ", "
                                                             << param.angularVelocities[iSurface][2] << "]" << std::endl;
    }
    for (plint iSurface = 0; iSurface < (plint) param.rotationAxisPoints.size(); iSurface++) {
        pcout << "rotationAxisPoints[" << iSurface << "] = [" << param.rotationAxisPoints[iSurface][0] << ", "
                                                              << param.rotationAxisPoints[iSurface][1] << ", "
                                                              << param.rotationAxisPoints[iSurface][2] << "]" << std::endl;
    }
    pcout << "outflowBcType = " << param.outflowBcType << std::endl;

    for (int iSponge = 0; iSponge < 6; iSponge++) {
        pcout << "spongeWidths[" << iSponge << "] = " << param.spongeWidths[iSponge] << std::endl;
    }
    pcout << "useSmagorinskySponges = " << (param.useSmagorinskySponges ? "true" : "false") << std::endl;
    if (param.useSmagorinskySponges) {
        pcout << "targetCSmago = " << param.targetCSmago << std::endl;
    }
    for (plint iSurface = 0; iSurface < (plint) param.torqueAxisPoints.size(); iSurface++) {
        pcout << "torqueAxisPoints[" << iSurface << "] = [" << param.torqueAxisPoints[iSurface][0] << ", "
                                                            << param.torqueAxisPoints[iSurface][1] << ", "
                                                            << param.torqueAxisPoints[iSurface][2] << "]" << std::endl;
    }
    for (plint iSurface = 0; iSurface < (plint) param.torqueAxisDirections.size(); iSurface++) {
        pcout << "torqueAxisDirections[" << iSurface << "] = [" << param.torqueAxisDirections[iSurface][0] << ", "
                                                                << param.torqueAxisDirections[iSurface][1] << ", "
                                                                << param.torqueAxisDirections[iSurface][2] << "]" << std::endl;
    }

    pcout << "useParallelIO = " << (param.useParallelIO ? "true" : "false") << std::endl;
    pcout << "precision = " << (param.precision == FLT ? "FLT" :
            (param.precision == DBL ? "DBL" :
             (param.precision == LDBL ? "LDBL" :
              "INF"))) << std::endl;

    pcout << "lx = " << param.lx << std::endl;
    pcout << "ly = " << param.ly << std::endl;
    pcout << "lz = " << param.lz << std::endl;
    pcout << "nx = " << param.nx << std::endl;
    pcout << "ny = " << param.ny << std::endl;
    pcout << "nz = " << param.nz << std::endl;
    pcout << "rho_LB = " << param.rho_LB << std::endl;
    pcout << "inletVelocity_LB = [" << param.inletVelocity_LB[0] << ", "
                                    << param.inletVelocity_LB[1] << ", "
                                    << param.inletVelocity_LB[2] << "]" << std::endl;
    for (plint iSurface = 0; iSurface < (plint) param.angularVelocities_LB.size(); iSurface++) {
        pcout << "angularVelocities_LB[" << iSurface << "] = [" << param.angularVelocities_LB[iSurface][0] << ", "
                                                                << param.angularVelocities_LB[iSurface][1] << ", "
                                                                << param.angularVelocities_LB[iSurface][2] << "]" << std::endl;
    }
    for (plint iSurface = 0; iSurface < (plint) param.rotationAxisPoints_LB.size(); iSurface++) {
        pcout << "rotationAxisPoints_LB[" << iSurface << "] = [" << param.rotationAxisPoints_LB[iSurface][0] << ", "
                                                                 << param.rotationAxisPoints_LB[iSurface][1] << ", "
                                                                 << param.rotationAxisPoints_LB[iSurface][2] << "]" << std::endl;
    }
    for (int iSponge = 0; iSponge < 6; iSponge++) {
        pcout << "numSpongeCells[" << iSponge << "] = " << param.numSpongeCells[iSponge] << std::endl;
    }
    for (plint iSurface = 0; iSurface < (plint) param.torqueAxisPoints_LB.size(); iSurface++) {
        pcout << "torqueAxisPoints_LB[" << iSurface << "] = [" << param.torqueAxisPoints_LB[iSurface][0] << ", "
                                                               << param.torqueAxisPoints_LB[iSurface][1] << ", "
                                                               << param.torqueAxisPoints_LB[iSurface][2] << "]" << std::endl;
    }
    pcout << "Re = " << param.inletVelocity[param.flowDirection] * param.characteristicLength / param.nu << std::endl;
    pcout << "omega = " << param.omega << std::endl;
    pcout << "tau = " << 1.0 / param.omega << std::endl;
    pcout << "dx = " << param.dx << std::endl;
    pcout << "dt / dx = " << param.dt / param.dx << std::endl;
    pcout << "dt / (dx * dx) = " << param.dt / (param.dx * param.dx) << std::endl;
    pcout << "uLB = " << param.dt / param.dx * param.inletVelocity[param.flowDirection] << std::endl;
    pcout << "physicalLocation = (" << param.physicalLocation[0] << ", " << param.physicalLocation[1] << ", "
          << param.physicalLocation[2] << ")" << std::endl;
    pcout << "inlet = [" << param.inlet.x0 << ", " << param.inlet.x1 << ", "
                         << param.inlet.y0 << ", " << param.inlet.y1 << ", "
                         << param.inlet.z0 << ", " << param.inlet.z1 << "]" << std::endl;
    pcout << "outlet = [" << param.outlet.x0 << ", " << param.outlet.x1 << ", "
                          << param.outlet.y0 << ", " << param.outlet.y1 << ", "
                          << param.outlet.z0 << ", " << param.outlet.z1 << "]" << std::endl;
    pcout << "lateral1 = [" << param.lateral1.x0 << ", " << param.lateral1.x1 << ", "
                            << param.lateral1.y0 << ", " << param.lateral1.y1 << ", "
                            << param.lateral1.z0 << ", " << param.lateral1.z1 << "]" << std::endl;
    pcout << "lateral2 = [" << param.lateral2.x0 << ", " << param.lateral2.x1 << ", "
                            << param.lateral2.y0 << ", " << param.lateral2.y1 << ", "
                            << param.lateral2.z0 << ", " << param.lateral2.z1 << "]" << std::endl;
    pcout << "lateral3 = [" << param.lateral3.x0 << ", " << param.lateral3.x1 << ", "
                            << param.lateral3.y0 << ", " << param.lateral3.y1 << ", "
                            << param.lateral3.z0 << ", " << param.lateral3.z1 << "]" << std::endl;
    pcout << "lateral4 = [" << param.lateral4.x0 << ", " << param.lateral4.x1 << ", "
                            << param.lateral4.y0 << ", " << param.lateral4.y1 << ", "
                            << param.lateral4.z0 << ", " << param.lateral4.z1 << "]" << std::endl;
    pcout << "smallEnvelopeWidth = " << param.smallEnvelopeWidth << std::endl;
    pcout << "largeEnvelopeWidth = " << param.largeEnvelopeWidth << std::endl;
    pcout << "saveDynamicContent = " << (param.saveDynamicContent ? "true" : "false") << std::endl;
    pcout << "fileNamePadding = " << param.fileNamePadding << std::endl;
    pcout << "incompressibleModel = " << (param.incompressibleModel ? "true" : "false") << std::endl;
    pcout << std::endl;
}

Array<T,3> getVelocity(Array<T,3> targetValue, plint iIter, plint startIter)
{
    return (targetValue * util::sinIncreasingFunction<T>(iIter, startIter));
}

void initializeImmersedSurfaceData(SimulationParameters& param)
{
    std::vector<std::string> allSurfaceFileNames;
    // Later we count on the fact that moving surfaces are included in the global vectors in the beginning.
    allSurfaceFileNames.insert(allSurfaceFileNames.end(), param.movingSurfaceFileNames.begin(), param.movingSurfaceFileNames.end());
    allSurfaceFileNames.insert(allSurfaceFileNames.end(), param.staticSurfaceFileNames.begin(), param.staticSurfaceFileNames.end());

    param.vertices.resize(0);
    param.areas.resize(0);
    param.flags.resize(0);
    param.startIds.resize(0);

    for (plint iSurface = 0; iSurface < param.numSurfaces; iSurface++) {
        TriangleSet<T> *surfaceTriangleSet = new TriangleSet<T>(allSurfaceFileNames[iSurface], param.precision);
        surfaceTriangleSet->writeBinarySTL(outDir + allSurfaceFileNames[iSurface]);

        if (param.refineSurfaceMeshes) {
            bool succeeded = surfaceTriangleSet->refineRecursively(param.targetMaxEdgeLength,
                    param.maxNumSurfaceRefinements);
            if (!succeeded) {
                pcout << std::endl;
                pcout << "WARNING: The target maximum triangle edge length " << param.targetMaxEdgeLength
                      << " for the immersed surface " << iSurface << std::endl
                      << "         was not reached after " << param.maxNumSurfaceRefinements << " refinement iterations."
                      << std::endl;
                pcout << std::endl;
            }
            FileName newSurfaceFileName(allSurfaceFileNames[iSurface]);
            newSurfaceFileName.setName(outDir + newSurfaceFileName.getName() + "_refined");
            surfaceTriangleSet->writeBinarySTL(newSurfaceFileName.get());
        }

        T maxEdgeLength = surfaceTriangleSet->getMaxEdgeLength();
        surfaceTriangleSet->scale(1.0 / param.dx);
        surfaceTriangleSet->translate(-param.physicalLocation / param.dx);

        ConnectedTriangleSet<T> connectedTriangleSet(*surfaceTriangleSet);
        delete surfaceTriangleSet;
        plint numVertices = connectedTriangleSet.getNumVertices();
        plint numTriangles = connectedTriangleSet.getNumTriangles();

        pcout << "The immersed surface " << iSurface <<" has " << numVertices
              << " vertices and " << numTriangles << " triangles." << std::endl;
        pcout << "The immersed surface " << iSurface <<" has a maximum triangle edge length of " << maxEdgeLength << std::endl;
        if (maxEdgeLength >= 4.0 * param.dx) {
            pcout << std::endl;
            pcout << "CAUTION: The maximum triangle edge length for the immersed surface " << iSurface << " is greater than "
                  << " 4 times dx."
                  << std::endl;
            pcout << "         The immersed boundary method will not work correctly. Surface refinement is necessary."
                  << std::endl;
            pcout << std::endl;
            exit(1);
        } else if (maxEdgeLength > param.dx) {
            pcout << std::endl;
            pcout << "WARNING: The maximum triangle edge length for the immersed surface " << iSurface << " is greater than dx."
                  << std::endl;
            pcout << "         The immersed boundary method might not work in an optimal way. Surface refinement is recommended."
                  << std::endl;
            pcout << std::endl;
        }

        param.startIds.push_back(param.vertices.size());
        for (plint iVertex = 0; iVertex < numVertices; iVertex++) {
            param.vertices.push_back(connectedTriangleSet.getVertex(iVertex));
            T area;
            Array<T,3> unitNormal;
            connectedTriangleSet.computeVertexAreaAndUnitNormal(iVertex, area, unitNormal);
            param.areas.push_back(area);
            param.flags.push_back(iSurface);
        }

        param.allSurfaces.push_back(connectedTriangleSet);  // Keep the initial surface topology.
    }
}

class VelFunction {
public:
    VelFunction(SimulationParameters& param_)
        : param(param_)
    { }

    Array<T,3> operator()(pluint id)
    {
        plint iSurface = -1;
        for (plint iMovingSurface = 0; iMovingSurface < param.numMovingSurfaces; iMovingSurface++) {
            plint startId = param.startIds[iMovingSurface];
            plint numVertices = param.allSurfaces[iMovingSurface].getNumVertices();
            if ((plint) id >= startId && (plint) id < startId + numVertices) {
                iSurface = iMovingSurface;
                break;
            }
        }

        if (iSurface == -1) {   // We are on a static surface.
            return Array<T,3>((T) 0, (T) 0, (T) 0);
        }

        Array<T,3> const& p0 = param.rotationAxisPoints_LB[iSurface];
        Array<T,3> const& omegaMax = param.angularVelocities_LB[iSurface];
        Array<T,3> const& position = param.vertices[id];

        Array<T,3> omega = getVelocity(omegaMax, param.nextIter, param.startIter);

        return(getExactRotationalVelocity(position, omega, p0));
    }

private:
    SimulationParameters& param;
};

class SurfaceNormalFunction {
public:
    SurfaceNormalFunction(SimulationParameters& param_)
        : param(param_)
    { }

    Array<T,3> operator()(pluint id)
    {
        plint iSurface = -1;
        for (plint i = 0; i < param.numSurfaces - 1; i++) {
            if ((plint) id >= param.startIds[i] && (plint) id < param.startIds[i+1]) {
                iSurface = i;
                break;
            }
        }
        if (iSurface == -1) {
            iSurface = param.numSurfaces - 1;
        }
        PLB_ASSERT(iSurface >= 0);

        plint localId = id - param.startIds[iSurface];
        T area;
        Array<T,3> unitNormal;
        param.allSurfaces[iSurface].computeVertexAreaAndUnitNormal(localId, area, unitNormal, &param.vertices,
                param.startIds[iSurface]);

        return unitNormal;
    }

private:
    SimulationParameters& param;
};

void updateMovingSurfaces(SimulationParameters& param, plint iIter)
{
    for (plint iMovingSurface = 0; iMovingSurface < param.numMovingSurfaces; iMovingSurface++) {
        Array<T,3> const& n = param.rotationAxisUnitVectors[iMovingSurface];
        Array<T,3> const& p0 = param.rotationAxisPoints_LB[iMovingSurface];
        Array<T,3> const& omegaMax = param.angularVelocities_LB[iMovingSurface];
        Array<T,3> omega = getVelocity(omegaMax, iIter, param.startIter);

        plint numVertices = param.allSurfaces[iMovingSurface].getNumVertices();
        plint startId = param.startIds[iMovingSurface];
        for (plint iVertex = 0; iVertex < numVertices; iVertex++) {
            plint id = iVertex + startId;
            Array<T,3>& position = param.vertices[id];
            position = getRotatedPosition(position, omega, n, p0);
        }
    }
}

void outputMovingSurfaces(SimulationParameters& param, std::string baseName, plint iIter)
{
    plint numDigits = util::val2str(param.numMovingSurfaces).length();
    for (plint iMovingSurface = 0; iMovingSurface < param.numMovingSurfaces; iMovingSurface++) {
        TriangleSet<T>* triangleSet = param.allSurfaces[iMovingSurface].toTriangleSet(param.precision,
                &param.vertices, param.startIds[iMovingSurface]);
        PLB_ASSERT(triangleSet != 0);
        triangleSet->scale(param.dx);
        triangleSet->translate(param.physicalLocation);

        std::string fname = createFileName(
                createFileName(baseName + "moving_surface_", iMovingSurface, numDigits+1)+"_", iIter, param.fileNamePadding)
            + ".stl";
        triangleSet->writeBinarySTL(fname);

        delete triangleSet;
    }
}

void saveMovingSurfaces(SimulationParameters& param, std::string baseName, plint iIter)
{
    // Checkpointing of moving surfaces is kept to a very basic level for simplicity.
    // We assume that if one process exits, then all the others will exit as well.
    if (global::mpi().isMainProcessor()) {
        std::string fname = createFileName(baseName + "surfaces_", iIter, param.fileNamePadding) + ".dat";
        FILE* fp = fopen(fname.c_str(), "wb");
        PLB_ASSERT(fp != 0);

        plint numVertices = param.vertices.size();
        if ((plint) fwrite(&param.vertices[0], sizeof(Array<T,3>), numVertices, fp) != numVertices) {
            fclose(fp);
            remove(fname.c_str());
            std::cout << "Error in saving surface data." << std::endl;
            exit(1);
        }
    }

    // If all the processes do not exit in the case that one does, then the code
    // will deliberately hang here.
    global::mpi().barrier();
}

void readMovingSurfaces(SimulationParameters& param, std::string baseName, plint iIter)
{
    // Checkpointing of moving surfaces is kept to a very basic level for simplicity.
    // We assume that if one process exits, then all the others will exit as well.
    // We also assume that all processes can read a file from the filesystem.
    std::string fname = createFileName(baseName + "surfaces_", iIter, param.fileNamePadding) + ".dat";
    FILE* fp = fopen(fname.c_str(), "rb");
    PLB_ASSERT(fp != 0);

    plint numVertices = param.vertices.size();
    if ((plint) fread(&param.vertices[0], sizeof(Array<T,3>), numVertices, fp) != numVertices) {
        std::cout << "Error in loading surface data." << std::endl;
        exit(1);
    }

    // If all the processes do not exit in the case that one or more do, then the code
    // will deliberately hang here.
    global::mpi().barrier();
}

Box3D boundAllSurfaces(SimulationParameters& param)
{
    Box3D totalBound;
    for (plint iSurface = 0; iSurface < param.numSurfaces; iSurface++) {
        TriangleSet<T>* triangleSet = param.allSurfaces[iSurface].toTriangleSet(param.precision,
                &param.vertices, param.startIds[iSurface]);
        PLB_ASSERT(triangleSet != 0);
        Cuboid<T> bCuboid = triangleSet->getBoundingCuboid();
        delete triangleSet;

        Array<T,3> const& llc = bCuboid.lowerLeftCorner;
        Array<T,3> const& urc = bCuboid.upperRightCorner;
        Box3D localBound;
        localBound.x0 = llc[0];
        localBound.y0 = llc[1];
        localBound.z0 = llc[2];
        localBound.x1 = urc[0] + 1;
        localBound.y1 = urc[1] + 1;
        localBound.z1 = urc[2] + 1;
        localBound = localBound.enlarge(2);

        if (iSurface == 0) {
            totalBound = localBound;
        } else {
            totalBound = bound(totalBound, localBound);
        }
    }

    return totalBound;
}

void createFluidBlocks(SimulationParameters& param, MultiBlockLattice3D<T,DESCRIPTOR>*& lattice,
        MultiScalarField3D<T>*& rhoBar, MultiTensorField3D<T,3>*& j, MultiContainerBlock3D*& container, 
        std::vector<MultiBlock3D*>& lattice_rho_bar_j_arg)
{
    Dynamics<T,DESCRIPTOR> *dynamics = new SmagorinskyBGKdynamics<T,DESCRIPTOR>(param.omega, param.cSmago);
    param.incompressibleModel = false;
    pcout << "Dynamics: Smagorinsky BGK." << std::endl;

    Box3D fullDomain(0, param.nx-1, 0, param.ny-1, 0, param.nz-1);
    lattice = generateMultiBlockLattice<T,DESCRIPTOR>(fullDomain, dynamics->clone(),
            param.smallEnvelopeWidth).release();
    defineDynamics(*lattice, lattice->getBoundingBox(), dynamics->clone());
    delete dynamics;
    lattice->toggleInternalStatistics(false);

    rhoBar = generateMultiScalarField<T>((MultiBlock3D&) *lattice, param.largeEnvelopeWidth).release();
    rhoBar->toggleInternalStatistics(false);

    j = generateMultiTensorField<T,3>((MultiBlock3D&) *lattice, param.largeEnvelopeWidth).release();
    j->toggleInternalStatistics(false);

    container = new MultiContainerBlock3D((MultiBlock3D&) *rhoBar);

    lattice_rho_bar_j_arg.clear();
    lattice_rho_bar_j_arg.push_back(lattice);
    lattice_rho_bar_j_arg.push_back(rhoBar);
    lattice_rho_bar_j_arg.push_back(j);
    integrateProcessingFunctional(
            new ExternalRhoJcollideAndStream3D<T,DESCRIPTOR>(),
            lattice->getBoundingBox(), lattice_rho_bar_j_arg, 0);
    integrateProcessingFunctional(
            new BoxRhoBarJfunctional3D<T,DESCRIPTOR>(),
            lattice->getBoundingBox(), lattice_rho_bar_j_arg, 3); // Boundary conditions are executed at levels 0, 1 and 2.

    // Integrate the immersed boundary processors in the lattice multi-block.

    std::vector<MultiBlock3D*> args;

    plint pl = 4;

    args.resize(0);
    args.push_back(container);
    integrateProcessingFunctional(
            new InstantiateImmersedWallDataWithIndexedTagging3D<T>(param.vertices, param.areas, param.flags),
            container->getBoundingBox(), *lattice, args, pl);
    pl++;

    for (plint i = 0; i < param.ibIter; i++) {
        args.resize(0);
        args.push_back(rhoBar);
        args.push_back(j);
        args.push_back(container);
        integrateProcessingFunctional(
            new IndexedInamuroIteration3D<T,VelFunction>(
                VelFunction(param), 1.0 / param.omega, param.incompressibleModel),
            rhoBar->getBoundingBox(), *lattice, args, pl);
        pl++;
    }
}

void createSpongeZones(SimulationParameters const& param, bool continueSimulation,
        MultiBlockLattice3D<T,DESCRIPTOR> *lattice)
{
    plint totalNumSpongeCells = 0;
    for (int iSponge = 0; iSponge < 6; iSponge++) {
        totalNumSpongeCells += param.numSpongeCells[iSponge];
    }

    if (!continueSimulation && totalNumSpongeCells > 0) {
        if (param.useSmagorinskySponges) {
            T bulkValue, targetValue;
            pcout << "Generating Smagorinsky sponge zones." << std::endl;
            bulkValue = param.cSmago;
            targetValue = param.targetCSmago;

            std::vector<MultiBlock3D*> args;
            args.push_back(lattice);
            applyProcessingFunctional(new SmagorinskySpongeZone3D<T,DESCRIPTOR>(
                        param.nx, param.ny, param.nz, bulkValue, targetValue, param.numSpongeCells),
                    lattice->getBoundingBox(), args);
        } else {
            T bulkValue;
            pcout << "Generating viscosity sponge zones." << std::endl;
            bulkValue = param.omega;

            std::vector<MultiBlock3D*> args;
            args.push_back(lattice);
            applyProcessingFunctional(new ViscositySpongeZone3D<T,DESCRIPTOR>(
                        param.nx, param.ny, param.nz, bulkValue, param.numSpongeCells),
                    lattice->getBoundingBox(), args);
        }
    }
}

void outerDomainBoundaryConditions(SimulationParameters const& param,
        MultiBlockLattice3D<T,DESCRIPTOR> *lattice,
        MultiScalarField3D<T> *rhoBar, MultiTensorField3D<T,3> *j,
        OnLatticeBoundaryCondition3D<T,DESCRIPTOR> *bc)
{
    Array<T,3> velocity = getVelocity(param.inletVelocity_LB, 0, param.startIter);

    if (param.lateralPeriodic) {
        pcout << "Periodic lateral boundaries." << std::endl;

        lattice->periodicity().toggleAll(true);
        rhoBar->periodicity().toggleAll(true);
        j->periodicity().toggleAll(true);

        lattice->periodicity().toggle(param.flowDirection, false);
        rhoBar->periodicity().toggle(param.flowDirection, false);
        j->periodicity().toggle(param.flowDirection, false);

        // Inlet boundary condition.

        if (param.flowDirection == 0) {
            bc->addVelocityBoundary0N(param.inlet, *lattice);
        } else if (param.flowDirection == 1) {
            bc->addVelocityBoundary1N(param.inlet, *lattice);
        } else {
            bc->addVelocityBoundary2N(param.inlet, *lattice);
        }

        setBoundaryVelocity(*lattice, param.inlet, velocity);
    } else {
        pcout << "Free-slip lateral boundaries." << std::endl;

        lattice->periodicity().toggleAll(false);
        rhoBar->periodicity().toggleAll(false);
        j->periodicity().toggleAll(false);

        bc->setVelocityConditionOnBlockBoundaries(*lattice, param.inlet, boundary::dirichlet);
        setBoundaryVelocity(*lattice, param.inlet, velocity);

        bc->setVelocityConditionOnBlockBoundaries(*lattice, param.lateral1, boundary::freeslip);
        bc->setVelocityConditionOnBlockBoundaries(*lattice, param.lateral2, boundary::freeslip);
        bc->setVelocityConditionOnBlockBoundaries(*lattice, param.lateral3, boundary::freeslip);
        bc->setVelocityConditionOnBlockBoundaries(*lattice, param.lateral4, boundary::freeslip);
        setBoundaryVelocity(*lattice, param.lateral1, velocity);
        setBoundaryVelocity(*lattice, param.lateral2, velocity);
        setBoundaryVelocity(*lattice, param.lateral3, velocity);
        setBoundaryVelocity(*lattice, param.lateral4, velocity);
    }

    // Outlet boundary condition.

    if (param.outflowBcType == 0) {
        if (param.flowDirection == 0) {
            bc->addPressureBoundary0P(param.outlet, *lattice);
        } else if (param.flowDirection == 1) {
            bc->addPressureBoundary1P(param.outlet, *lattice);
        } else {
            bc->addPressureBoundary2P(param.outlet, *lattice);
        }
        setBoundaryDensity(*lattice, param.outlet, param.rho_LB);
        setBoundaryVelocity(*lattice, param.outlet, velocity);
    } else if (param.outflowBcType == 1) {
        if (param.flowDirection == 0) {
            integrateProcessingFunctional(new FluidPressureOutlet3D<T,DESCRIPTOR,0,+1>(), param.outlet, *lattice, 2);
        } else if (param.flowDirection == 1) {
            integrateProcessingFunctional(new FluidPressureOutlet3D<T,DESCRIPTOR,1,+1>(), param.outlet, *lattice, 2);
        } else {
            integrateProcessingFunctional(new FluidPressureOutlet3D<T,DESCRIPTOR,2,+1>(), param.outlet, *lattice, 2);
        }
        setBoundaryVelocity(*lattice, param.outlet, velocity);
    } else {
        if (param.lateralPeriodic) {
            Box3D globalDomain(lattice->getBoundingBox());
            if (param.flowDirection == 0) {
                globalDomain.y0 -= 2; // y-periodicity
                globalDomain.y1 += 2;
                globalDomain.z0 -= 2; // z-periodicity
                globalDomain.z1 += 2;
            } else if (param.flowDirection == 1) {
                globalDomain.z0 -= 2; // z-periodicity
                globalDomain.z1 += 2;
                globalDomain.x0 -= 2; // x-periodicity
                globalDomain.x1 += 2;
            } else {
                globalDomain.x0 -= 2; // x-periodicity
                globalDomain.x1 += 2;
                globalDomain.y0 -= 2; // y-periodicity
                globalDomain.y1 += 2;
            }

            std::vector<MultiBlock3D*> bcargs;
            bcargs.push_back(lattice);
            bcargs.push_back(rhoBar);
            bcargs.push_back(j);
            int type = 1;
            integrateProcessingFunctional(
                    new VirtualOutlet<T,DESCRIPTOR>(param.rho_LB, globalDomain, type),
                    param.outlet, bcargs, 2);
            setBoundaryVelocity(*lattice, param.outlet, velocity);
        } else {
            Box3D globalDomain(lattice->getBoundingBox());
            std::vector<MultiBlock3D*> bcargs;
            bcargs.push_back(lattice);
            bcargs.push_back(rhoBar);
            bcargs.push_back(j);
            int type = 1;
            integrateProcessingFunctional(
                    new VirtualOutlet<T,DESCRIPTOR>(param.rho_LB, globalDomain, type),
                    param.outlet, bcargs, 2);
            setBoundaryVelocity(*lattice, param.outlet, velocity);
        }
    }
}

void initializeSimulation(SimulationParameters& param, bool continueSimulation, std::string xmlRestartFileName,
        plint& iniIter, MultiBlockLattice3D<T,DESCRIPTOR> *lattice, std::vector<MultiBlock3D*>& lattice_rho_bar_j_arg,
        std::vector<MultiBlock3D*>& checkpointBlocks)
{
    if (!continueSimulation) {
        Array<T,3> velocity = getVelocity(param.inletVelocity_LB, 0, param.startIter);
        initializeAtEquilibrium(*lattice, lattice->getBoundingBox(), param.rho_LB, velocity);
        applyProcessingFunctional(
                new BoxRhoBarJfunctional3D<T,DESCRIPTOR>(),
                lattice->getBoundingBox(), lattice_rho_bar_j_arg);
        T energy = computeAverageEnergy(*lattice) * param.rho * (param.dx * param.dx) / (param.dt * param.dt);
        pcout << "Initial average kinetic energy: " << energy << std::endl;
    }
    pcout << std::endl;

    if (continueSimulation) {
        pcout << "Reading state of the simulation from file: " << xmlRestartFileName << std::endl;
        loadState(checkpointBlocks, iniIter, param.saveDynamicContent, xmlRestartFileName);
        lattice->resetTime(iniIter);
        readMovingSurfaces(param, param.baseFileName, iniIter);
        pcout << std::endl;
    }

    pcout << std::endl;
}

void writeVTK(SimulationParameters const& param, MultiBlockLattice3D<T,DESCRIPTOR> *lattice,
        plint iIter)
{
    T outletPressure = param.rho_LB;
    if (param.outflowBcType != 0) {
        outletPressure = computeAverageDensity(*lattice, param.outlet);
    }
    outletPressure *= param.rho * (param.dx * param.dx) / (param.dt * param.dt) * DESCRIPTOR<T>::cs2;

    T pressureScale = param.rho * (param.dx * param.dx) / (param.dt * param.dt) * DESCRIPTOR<T>::cs2;
    T pressureOffset = param.ambientPressure - outletPressure;

    if (param.outputInDomain) {
        std::string fname = createFileName(outDir + "domain_", iIter, param.fileNamePadding);
        VtkImageOutput3D<T> vtkOut(fname, param.dx, param.physicalLocation);

        std::auto_ptr<MultiScalarField3D<T> > p = computeDensity(*lattice, param.outputDomain);
        vtkOut.writeData<float>(*p, "pressure", pressureScale, pressureOffset);
        p.reset();

        std::auto_ptr<MultiTensorField3D<T,3> > v = computeVelocity(*lattice, param.outputDomain);
        vtkOut.writeData<float>(*computeNorm(*v), "velocityNorm", param.dx / param.dt);
        vtkOut.writeData<3,float>(*v, "velocity", param.dx / param.dt);
        std::auto_ptr<MultiTensorField3D<T,3> > vort = computeVorticity(*v);
        vtkOut.writeData<float>(*computeNorm(*vort), "vorticityNorm", 1.0 / param.dt);
        vtkOut.writeData<3,float>(*vort, "vorticity", 1.0 / param.dt);
    }

    if (param.outputOnSlices) {
        plint numXdigits = util::val2str(param.xSlices.size()).length();
        for (plint i = 0; i < (plint) param.xSlices.size(); i++) {
            Box3D box_x = param.xSlices[i];

            std::string fname = createFileName(
                    createFileName(outDir + "x_slice_", i, numXdigits+1)+"_", iIter, param.fileNamePadding);
            VtkImageOutput3D<T> vtkOut_x(fname, param.dx, param.physicalLocation);

            std::auto_ptr<MultiScalarField3D<T> > px = computeDensity(*lattice, box_x);
            vtkOut_x.writeData<float>(*px, "pressure", pressureScale, pressureOffset);
            px.reset();

            std::auto_ptr<MultiTensorField3D<T,3> > vx = computeVelocity(*lattice, box_x);
            vtkOut_x.writeData<float>(*computeNorm(*vx), "velocityNorm", param.dx / param.dt);
            vtkOut_x.writeData<3,float>(*vx, "velocity", param.dx / param.dt);
            std::auto_ptr<MultiTensorField3D<T,3> > vortx = computeVorticity(*vx);
            vtkOut_x.writeData<float>(*computeNorm(*vortx), "vorticityNorm", 1.0 / param.dt);
            vtkOut_x.writeData<3,float>(*vortx, "vorticity", 1.0 / param.dt);
        }

        plint numYdigits = util::val2str(param.ySlices.size()).length();
        for (plint i = 0; i < (plint) param.ySlices.size(); i++) {
            Box3D box_y = param.ySlices[i];

            std::string fname = createFileName(
                    createFileName(outDir + "y_slice_", i, numYdigits+1)+"_", iIter, param.fileNamePadding);
            VtkImageOutput3D<T> vtkOut_y(fname, param.dx, param.physicalLocation);

            std::auto_ptr<MultiScalarField3D<T> > py = computeDensity(*lattice, box_y);
            vtkOut_y.writeData<float>(*py, "pressure", pressureScale, pressureOffset);
            py.reset();

            std::auto_ptr<MultiTensorField3D<T,3> > vy = computeVelocity(*lattice, box_y);
            vtkOut_y.writeData<float>(*computeNorm(*vy), "velocityNorm", param.dx / param.dt);
            vtkOut_y.writeData<3,float>(*vy, "velocity", param.dx / param.dt);
            std::auto_ptr<MultiTensorField3D<T,3> > vorty = computeVorticity(*vy);
            vtkOut_y.writeData<float>(*computeNorm(*vorty), "vorticityNorm", 1.0 / param.dt);
            vtkOut_y.writeData<3,float>(*vorty, "vorticity", 1.0 / param.dt);
        }

        plint numZdigits = util::val2str(param.zSlices.size()).length();
        for (plint i = 0; i < (plint) param.zSlices.size(); i++) {
            Box3D box_z = param.zSlices[i];

            std::string fname = createFileName(
                    createFileName(outDir + "z_slice_", i, numZdigits+1)+"_", iIter, param.fileNamePadding);
            VtkImageOutput3D<T> vtkOut_z(fname, param.dx, param.physicalLocation);

            std::auto_ptr<MultiScalarField3D<T> > pz = computeDensity(*lattice, box_z);
            vtkOut_z.writeData<float>(*pz, "pressure", pressureScale, pressureOffset);
            pz.reset();

            std::auto_ptr<MultiTensorField3D<T,3> > vz = computeVelocity(*lattice, box_z);
            vtkOut_z.writeData<float>(*computeNorm(*vz), "velocityNorm", param.dx / param.dt);
            vtkOut_z.writeData<3,float>(*vz, "velocity", param.dx / param.dt);
            std::auto_ptr<MultiTensorField3D<T,3> > vortz = computeVorticity(*vz);
            vtkOut_z.writeData<float>(*computeNorm(*vortz), "vorticityNorm", 1.0 / param.dt);
            vtkOut_z.writeData<3,float>(*vortz, "vorticity", 1.0 / param.dt);
        }
    }
}


int main(int argc, char* argv[])
{
    plbInit(&argc, &argv);

    std::cout.precision(10);
    std::scientific(std::cout);

    // Command-line arguments

    if (argc != 2 && argc != 3) {
        pcout << "Usage: " << argv[0] << " xml-input-file-name [xml-continue-file-name]" << std::endl;
        exit(1);
    }

    std::string xmlInputFileName;
    xmlInputFileName = std::string(argv[1]);

    std::string xmlRestartFileName;
    bool continueSimulation = false;
    if (argc == 3) {
        xmlRestartFileName = std::string(argv[2]);
        continueSimulation = true;
    }

    // Set the simulation parameters.

    SimulationParameters param;

    readUserDefinedSimulationParameters(xmlInputFileName, param);
    calculateDerivedSimulationParameters(param);

    global::IOpolicy().activateParallelIO(param.useParallelIO);

    // Immersed surfaces.

    pcout << "Processing immersed surface geometries." << std::endl;

    initializeImmersedSurfaceData(param);

    // Fluid.

    pcout << "Generating fluid blocks." << std::endl;

    MultiBlockLattice3D<T,DESCRIPTOR> *lattice = 0;
    MultiScalarField3D<T> *rhoBar = 0;
    MultiTensorField3D<T,3> *j = 0;
    MultiContainerBlock3D *container = 0;
    std::vector<MultiBlock3D*> lattice_rho_bar_j_arg;

    createFluidBlocks(param, lattice, rhoBar, j, container, lattice_rho_bar_j_arg);
    printSimulationParameters(param);

    // Create boundary sponges.

    createSpongeZones(param, continueSimulation, lattice);

    // Boundary conditions.

    pcout << "Generating outer domain boundary conditions." << std::endl;

    OnLatticeBoundaryCondition3D<T,DESCRIPTOR> *bc = createLocalBoundaryCondition3D<T,DESCRIPTOR>();
    outerDomainBoundaryConditions(param, lattice, rhoBar, j, bc);
    delete bc;

    // Initialization.

    std::vector<MultiBlock3D*> checkpointBlocks;
    checkpointBlocks.push_back(lattice);
    checkpointBlocks.push_back(rhoBar);
    checkpointBlocks.push_back(j);

    plint iniIter = 0;

    initializeSimulation(param, continueSimulation, xmlRestartFileName, iniIter, lattice, lattice_rho_bar_j_arg,
            checkpointBlocks);

    FILE *fpEnergy = 0;
    if (global::mpi().isMainProcessor()) {
        std::string fileName = outDir + "average_energy.dat";
        fpEnergy = fopen(fileName.c_str(), continueSimulation ? "a" : "w");
        PLB_ASSERT(fpEnergy != 0);
    }

    std::vector<FILE *> fpForces(param.numSurfaces);
    if (global::mpi().isMainProcessor()) {
        plint numDigits = util::val2str(param.numSurfaces).length();
        for (plint iSurface = 0; iSurface < param.numSurfaces; iSurface++) {
            std::string fileName;
            if (param.numSurfaces != 1) {
                fileName = createFileName(outDir + "total_force_on_surface_", iSurface, numDigits+1) + ".dat";
            } else {
                fileName = outDir + "total_force_on_surface.dat";
            }
            fpForces[iSurface] = fopen(fileName.c_str(), continueSimulation ? "a" : "w");
            PLB_ASSERT(fpForces[iSurface] != 0);
        }
    }

    std::vector<FILE *> fpTorques(param.numMovingSurfaces);
    if (global::mpi().isMainProcessor()) {
        plint numDigits = util::val2str(param.numMovingSurfaces).length();
        for (plint iSurface = 0; iSurface < param.numMovingSurfaces; iSurface++) {
            std::string fileName;
            if (param.numMovingSurfaces != 1) {
                fileName = createFileName(outDir + "total_torque_on_surface_", iSurface, numDigits+1) + ".dat";
            } else {
                fileName = outDir + "total_torque_on_surface.dat";
            }
            fpTorques[iSurface] = fopen(fileName.c_str(), continueSimulation ? "a" : "w");
            PLB_ASSERT(fpTorques[iSurface] != 0);
        }
    }

    // Starting iterations.

    pcout << "Starting simulation." << std::endl;
    bool stopExecution = false;
    param.nextIter = iniIter + 1;
    for (plint iIter = iniIter; iIter < param.maxIter && !stopExecution; iIter++) {
        param.nextIter = iIter + 1;
        if (iIter <= param.startIter) {
            Array<T,3> velocity = getVelocity(param.inletVelocity_LB, iIter, param.startIter);
            setBoundaryVelocity(*lattice, param.inlet, velocity);
        }

        if (iIter != iniIter && (iIter % param.statIter == 0 || iIter == param.maxIter - 1)) {
            pcout << "At iteration " << iIter << ", t = " << iIter * param.dt << std::endl;

            // Average kinetic energy.

            T energy = computeAverageEnergy(*lattice) * param.rho * (param.dx * param.dx) / (param.dt * param.dt);
            pcout << "Average kinetic energy: " << energy << std::endl;
            if (global::mpi().isMainProcessor()) {
                double t = (double) (iIter * param.dt);
                double ed = (double) energy;
                fprintf(fpEnergy, "% .8e\t% .8e\n", t, ed);
                fflush(fpEnergy);
            }

            // Forces on immersed surfaces.

            Box3D totalBound = boundAllSurfaces(param);
            T averageDensity = computeAverageDensity(*lattice, totalBound);
            T forceConversion = 2.0 * param.rho * (param.dx * param.dx * param.dx * param.dx) / (param.dt * param.dt);
            recomputeImmersedForce(SurfaceNormalFunction(param), param.omega, averageDensity, *lattice,
                    *container, param.largeEnvelopeWidth, lattice->getBoundingBox(), param.incompressibleModel);
            for (plint iSurface = 0; iSurface < param.numSurfaces; iSurface++) {
                Array<T,3> force = -reduceImmersedForce<T>(*container, iSurface) * forceConversion;
                if (param.numSurfaces != 1) {
                    pcout << "Force on surface " << iSurface << ": ";
                } else {
                    pcout << "Force on surface: ";
                }
                pcout << "(" << force[0] << ", " << force[1] << ", " << force[2] << ")" << std::endl;
                if (global::mpi().isMainProcessor()) {
                    double t = (double) (iIter * param.dt);
                    double f0 = (double) force[0];
                    double f1 = (double) force[1];
                    double f2 = (double) force[2];
                    fprintf(fpForces[iSurface], "% .8e\t% .8e\t% .8e\t% .8e\n", t, f0, f1, f2);
                    fflush(fpForces[iSurface]);
                }
            }

            // Torques on immersed moving surfaces.

            T torqueConversion = 2.0 * param.rho * (param.dx * param.dx * param.dx * param.dx * param.dx) / (param.dt * param.dt);
            for (plint iSurface = 0; iSurface < param.numMovingSurfaces; iSurface++) {
                Array<T,3> torque = -reduceAxialTorqueImmersed<T>(*container, param.torqueAxisPoints_LB[iSurface],
                        param.torqueAxisDirections[iSurface], iSurface) * torqueConversion;
                if (param.numMovingSurfaces != 1) {
                    pcout << "Torque on moving surface " << iSurface << ": ";
                } else {
                    pcout << "Torque on moving surface: ";
                }
                pcout << "(" << torque[0] << ", " << torque[1] << ", " << torque[2] << ")" << std::endl;
                if (global::mpi().isMainProcessor()) {
                    double t = (double) (iIter * param.dt);
                    double t0 = (double) torque[0];
                    double t1 = (double) torque[1];
                    double t2 = (double) torque[2];
                    fprintf(fpTorques[iSurface], "% .8e\t% .8e\t% .8e\t% .8e\n", t, t0, t1, t2);
                    fflush(fpTorques[iSurface]);
                }
            }

            pcout << "Time for one fluid iteration: " << global::timer("lb-iter").getTime() / (T) param.statIter << std::endl;
            global::timer("lb-iter").reset();
            pcout << std::endl;
        }

        if (iIter % param.outIter == 0 || iIter == param.maxIter - 1) {
            pcout << "Output to disk at iteration: " << iIter << std::endl;
            writeVTK(param, lattice, iIter);
            outputMovingSurfaces(param, outDir, iIter);
            pcout << std::endl;
        }

        if ((param.cpIter > 0 && iIter % param.cpIter == 0 && iIter != iniIter) ||
             iIter == param.maxIter - 1) {
            pcout << "Saving the state of the simulation at iteration: " << iIter << std::endl;
            saveState(checkpointBlocks, iIter, param.saveDynamicContent, param.xmlContinueFileName,
                    param.baseFileName, param.fileNamePadding);
            saveMovingSurfaces(param, param.baseFileName, iIter);
            pcout << std::endl;
        }

        //if (iIter % param.abIter == 0 || iIter == param.maxIter - 1) {
        if (iIter % param.abIter == 0) {
            stopExecution = abortExecution(param.abortFileName, checkpointBlocks, iIter,
                    param.saveDynamicContent, param.xmlContinueFileName,
                    param.baseFileName, param.fileNamePadding);

            if (stopExecution) {
                saveMovingSurfaces(param, param.baseFileName, iIter);
                pcout << "Aborting execution at iteration: " << iIter << std::endl;
                pcout << std::endl;
            }
        }

        // Here we could use either "param.nextIter" or "iIter".
        updateMovingSurfaces(param, param.nextIter);

        global::timer("lb-iter").start();
        lattice->executeInternalProcessors(); // Execute all processors and communicate appropriately.
        global::timer("lb-iter").stop();
        lattice->incrementTime();
    }

    if (global::mpi().isMainProcessor()) {
        fclose(fpEnergy);
        for (plint iSurface = 0; iSurface < param.numSurfaces; iSurface++) {
            fclose(fpForces[iSurface]);
        }
        for (plint iSurface = 0; iSurface < param.numMovingSurfaces; iSurface++) {
            fclose(fpTorques[iSurface]);
        }
    }

    delete container;
    delete j;
    delete rhoBar;
    delete lattice;

    exit(0);
}

