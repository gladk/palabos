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


#include "palabos3D.h"
#include "palabos3D.hh"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <map>
#include <sstream>
#include <vector>

using namespace plb;

typedef double T;
typedef Array<T,3> Velocity;

#define RESCALER ConvectiveNoForceRescaler
#define DESCRIPTOR descriptors::D3Q19Descriptor

struct SimulationParameters {
    /*
     * Parameters set by the user.
     */

    // Geometry.

    std::string staticSurfaceFileName;    // Files with the static immersed surface geometries.

    // Grid.
    std::string gridDensityFunctionFile;            // File with the discrete grid density function.
    plint minLeafLevel;                             // Octree leaf levels.
    plint maxLeafLevel;
    plint nBlock;                                   // Block size.

    // Boundary Conditions.
    T inletVelocity;                                // Inlet velocity x-component in physical units.

    int outflowBcType;                              // Type of the outflow boundary condition.
                                                    // If 0, then a constant velocity is imposed at the outlet.
                                                    // If 1, then a constant pressure is imposed at the outlet.
                                                    // If 2, then a velocity Neumann condition is imposed at the outlet.

    Array<T,6> spongeWidths;                        // Sponge zone widths (0.0 for no sponge zone for the corresponding
                                                    // lattice boundary).
    // Numerics.
    Precision precision;                    // Precision for geometric operations.
    T L_Ref;                                // Length to define the Reynolds number.
    T u_Ref;                                // Characteristic velocity.
    T u_LB;                                 // Lattice velocity (w.r.t the characteristic velocity).
    plint maxIter;                          // Maximum number of iterations at the coarsest level.

    // Fluid.

    T rho;                  // Fluid density in physical units
    T Re;                   // Reynolds number.
    T nu;                   // Fluid kinematic viscosity in physical units.
    T ambientPressure;      // Absolute stagnation pressure in physical units.

    // Output.

    std::string outDir;                 // Output directory.
    plint statIter;                     // Number of iterations for terminal output at the coarsest level.
    plint outIter;                      // Number of iterations for disk output at the coarsest level.
    bool computeAverages;               // Compute average and RMS values or not?
    plint avgIter;                      // Number of iterations to start averaging at the coarsest level.
    plint minOutputLevel;               // Minimum grid refinement level for output on disk.
    plint maxOutputLevel;               // Maximum grid refinement level for output on disk.

    bool outputInDomain;                // Save data on disk in a volume domain or not?
    Cuboid<T> outputCuboid;             // Volume domain for disk output.

    bool outputOnSlices;                // Save data on disk on a set of slices or not?
    std::vector<T> xPositions;          // Positions of the x-slices for output.
    std::vector<T> xyRange;             // y range of the x-slices.
    std::vector<T> xzRange;             // z range of the x-slices.
    std::vector<T> yPositions;          // Positions of the y-slices for output.
    std::vector<T> yzRange;             // z range of the y-slices.
    std::vector<T> yxRange;             // x range of the y-slices.
    std::vector<T> zPositions;          // Positions of the z-slices for output.
    std::vector<T> zxRange;             // x range of the z-slices.
    std::vector<T> zyRange;             // y range of the z-slices.

    plint cpIter;                       // Number of iterations for checkpointing.
    plint abIter;                       // Number of iterations for checking for user-driven program abortion.
    std::string abortFileName;          // File for signaling program abortion (inside the outDir).
    std::string xmlContinueFileName;    // XML file for restarting (inside the outDir).
    std::string baseFileName;           // Basename of the checkpoint files (inside the outDir).
    bool useParallelIO;                 // For a desktop PC this should be "false", for a cluster "true".

    /*
     * Parameters NOT set by the user.
     */

    std::string surfaceName;

    OctreeGridStructure ogs;
    Cuboid<T> fullDomain;
    plint finestLevel;
    T dxCoarsest, dtCoarsest;
    T dxFinest, dtFinest;

    T rho_LB;
    Array<T,3> inletVelocity_LB;
    std::vector<T> omega;
    Array<T,3> physicalLocation;
    plint smallEnvelopeWidth;
    plint mediumEnvelopeWidth;
    plint largeEnvelopeWidth;
    bool incompressibleModel;
    plint fileNamePadding;
    bool saveDynamicContent;

    std::map<plint, std::vector<Box3D> > outputDomains;     // Output domains per output level.
    std::vector<std::string> outputDomainNames;
};

T toLB(T physVal, plint direction, T dx, Array<T,3> const& location)
{
    PLB_ASSERT(direction >= 0 && direction <= 2);
    return((physVal - location[direction]) / dx);
}

Array<T,3> toLB(Array<T,3> const& physVal, T dx, Array<T,3> const& location)
{
    return((physVal - location) / dx);
}

void readUserDefinedSimulationParameters(std::string xmlInputFileName, SimulationParameters& param)
{
    XMLreader document(xmlInputFileName);

    std::vector<T> fullDomainX, fullDomainY, fullDomainZ;
    document["geometry"]["simulationDomain"]["x"].read(fullDomainX);
    PLB_ASSERT(fullDomainX.size() == 2 && fullDomainX[1] > fullDomainX[0]);
    document["geometry"]["simulationDomain"]["y"].read(fullDomainY);
    PLB_ASSERT(fullDomainY.size() == 2 && fullDomainY[1] > fullDomainY[0]);
    document["geometry"]["simulationDomain"]["z"].read(fullDomainZ);
    PLB_ASSERT(fullDomainZ.size() == 2 && fullDomainZ[1] > fullDomainZ[0]);
    param.fullDomain.lowerLeftCorner[0] = fullDomainX[0];
    param.fullDomain.lowerLeftCorner[1] = fullDomainY[0];
    param.fullDomain.lowerLeftCorner[2] = fullDomainZ[0];
    param.fullDomain.upperRightCorner[0] = fullDomainX[1];
    param.fullDomain.upperRightCorner[1] = fullDomainY[1];
    param.fullDomain.upperRightCorner[2] = fullDomainZ[1];

    document["geometry"]["staticSurfaceFileName"].read(param.staticSurfaceFileName);

    document["grid"]["gridDensityFunctionFile"].read(param.gridDensityFunctionFile);
    abortIfCannotOpenFileForReading(param.gridDensityFunctionFile);
    plint numLevels = 0;
    document["grid"]["numLevels"].read(numLevels);
    PLB_ASSERT(numLevels >= 1);
    plint maxOctreeLevel = 0;
    document["grid"]["maxOctreeLevel"].read(maxOctreeLevel);
    if (maxOctreeLevel - numLevels + 1 < 0) {
        maxOctreeLevel = numLevels - 1;
    }
    param.maxLeafLevel = maxOctreeLevel;
    param.minLeafLevel = maxOctreeLevel - numLevels + 1;
    document["grid"]["nBlock"].read(param.nBlock);
    PLB_ASSERT(param.nBlock >= 6);

    document["boundaryConditions"]["inletVelocity"].read(param.inletVelocity);

    document["boundaryConditions"]["outflowBcType"].read(param.outflowBcType);
    PLB_ASSERT(param.outflowBcType >= 0 && param.outflowBcType <= 2);

    std::vector<T> zoneWidths;
    document["boundaryConditions"]["spongeZones"]["xWidths"].read(zoneWidths);
    PLB_ASSERT(zoneWidths.size() == 2);
    param.spongeWidths[0] = zoneWidths[0];
    param.spongeWidths[1] = zoneWidths[1];
    zoneWidths.clear();
    document["boundaryConditions"]["spongeZones"]["yWidths"].read(zoneWidths);
    PLB_ASSERT(zoneWidths.size() == 2);
    param.spongeWidths[2] = zoneWidths[0];
    param.spongeWidths[3] = zoneWidths[1];
    zoneWidths.clear();
    document["boundaryConditions"]["spongeZones"]["zWidths"].read(zoneWidths);
    PLB_ASSERT(zoneWidths.size() == 2);
    param.spongeWidths[4] = zoneWidths[0];
    param.spongeWidths[5] = zoneWidths[1];
    zoneWidths.clear();

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
    document["numerics"]["characteristicLength"].read(param.L_Ref);
    PLB_ASSERT(util::greaterThan_abs(param.L_Ref, (T) 0));
    document["numerics"]["characteristicVelocity"].read(param.u_Ref);
    PLB_ASSERT(util::greaterThan_abs(param.u_Ref, (T) 0));
    document["numerics"]["uLB"].read(param.u_LB);
    document["numerics"]["maxIter"].read(param.maxIter);
    PLB_ASSERT(param.maxIter > 0);

    document["fluid"]["rho"].read(param.rho);
    document["fluid"]["Re"].read(param.Re);
    param.nu = param.u_Ref * param.L_Ref / param.Re;
    document["fluid"]["ambientPressure"].read(param.ambientPressure);

    std::string outDir;
    document["output"]["outDir"].read(outDir);
    if (outDir[outDir.size() - 1] != '/') {
        outDir += '/';
    }
    param.outDir = outDir;
    abortIfCannotCreateFileInDir(param.outDir, "plb-checkfile.txt");

    document["output"]["statIter"].read(param.statIter);
    PLB_ASSERT(param.statIter > 0);
    document["output"]["outIter"].read(param.outIter);
    PLB_ASSERT(param.outIter > 0);
    document["output"]["computeAverages"].read(param.computeAverages);
    
    if (param.computeAverages) {
        document["output"]["avgIter"].read(param.avgIter);
        PLB_ASSERT(param.avgIter >= 0);
    } else {
        param.avgIter = 0;
    }

    document["output"]["minOutputLevel"].read(param.minOutputLevel);
    document["output"]["maxOutputLevel"].read(param.maxOutputLevel);
    PLB_ASSERT(param.minOutputLevel <= param.maxOutputLevel);

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

    document["output"]["cpIter"].read(param.cpIter);
    document["output"]["abIter"].read(param.abIter);
    PLB_ASSERT(param.abIter > 0);
    document["output"]["abortFileName"].read(param.abortFileName);
    document["output"]["xmlContinueFileName"].read(param.xmlContinueFileName);
    document["output"]["baseFileName"].read(param.baseFileName);
    document["output"]["useParallelIO"].read(param.useParallelIO);

    param.abortFileName = param.outDir + param.abortFileName;
    param.xmlContinueFileName = param.outDir + param.xmlContinueFileName;
    param.baseFileName = param.outDir + param.baseFileName;
}

void createOctreeGridStructure(SimulationParameters &param)
{
    // Default parameters for the octree grid generation.

    bool useSamples = false;
    plint numSamples = -1;
    plint maxIter = 100;
    bool removeBlocks = true;
    bool fineToCoarse = true;
    int numLevelsToGroupBlocks = 0;
    int numLevelsToGroupOverlaps = -1;
    bool strongGrouping = false;
    bool verbose = true;
    bool stlOutput = true;
    std::string stlBaseName = "octree";
    T gridDensityScaleFactor = (T)1;

    OctreeGridGenerator<T> octreeGridGenerator(
            param.fullDomain,
            param.gridDensityFunctionFile,
            param.minLeafLevel,
            param.maxLeafLevel,
            param.nBlock,
            global::mpi().getSize(),
            gridDensityScaleFactor,
            useSamples,
            numSamples,
            maxIter,
            removeBlocks,
            fineToCoarse,
            numLevelsToGroupBlocks,
            numLevelsToGroupOverlaps,
            strongGrouping,
            param.outDir,
            verbose,
            stlOutput,
            stlBaseName);

    param.ogs = octreeGridGenerator.generateOctreeGridStructure();
    param.fullDomain = octreeGridGenerator.getFullDomain();
    param.dxFinest = octreeGridGenerator.getDxFinestLevel();

    param.ogs.writeXML(octreeGridGenerator.getOutDir() + "octreeGridStructure.xml");

    {
        XMLwriter gridInfoXML;
        Array<T,6> fullDomain(octreeGridGenerator.getFullDomain().x0(), octreeGridGenerator.getFullDomain().x1(),
                              octreeGridGenerator.getFullDomain().y0(), octreeGridGenerator.getFullDomain().y1(),
                              octreeGridGenerator.getFullDomain().z0(), octreeGridGenerator.getFullDomain().z1());
        gridInfoXML["fullDomain"].set<T,6>(fullDomain);
        gridInfoXML["dxFinestLevel"].set(octreeGridGenerator.getDxFinestLevel());
        gridInfoXML.print(octreeGridGenerator.getOutDir() + "octreeGridInfo.xml");
    }
}

Box3D computeOutputDomain(SimulationParameters const& param, Cuboid<T> const& outputCuboid, plint level)
{
    T dx = param.dxFinest * (T) util::intTwoToThePower(param.finestLevel - level);
    Box3D outputBox;
    outputBox.x0 = (plint) toLB(outputCuboid.x0(), 0, dx, param.physicalLocation);
    outputBox.x1 = (plint) toLB(outputCuboid.x1(), 0, dx, param.physicalLocation) + (plint) 1;
    outputBox.y0 = (plint) toLB(outputCuboid.y0(), 1, dx, param.physicalLocation);
    outputBox.y1 = (plint) toLB(outputCuboid.y1(), 1, dx, param.physicalLocation) + (plint) 1;
    outputBox.z0 = (plint) toLB(outputCuboid.z0(), 2, dx, param.physicalLocation);
    outputBox.z1 = (plint) toLB(outputCuboid.z1(), 2, dx, param.physicalLocation) + (plint) 1;
    return(outputBox);
}

void computeAllOutputDomains(SimulationParameters& param)
{
    std::vector<Cuboid<T> > outputCuboids;
    if (param.outputInDomain) {
        outputCuboids.push_back(param.outputCuboid);
        param.outputDomainNames.push_back("domain");
    }

    if (param.outputOnSlices) {
        plint numXdigits = util::val2str(param.xPositions.size()).length();
        for (plint i = 0; i < (plint) param.xPositions.size(); i++) {
            Array<T,3> llc(param.xPositions[i], param.xyRange[0], param.xzRange[0]);
            Array<T,3> urc(param.xPositions[i], param.xyRange[1], param.xzRange[1]);
            outputCuboids.push_back(Cuboid<T>(llc, urc));
            param.outputDomainNames.push_back(createFileName("slice_x_", i, numXdigits + 1));
        }

        plint numYdigits = util::val2str(param.yPositions.size()).length();
        for (plint i = 0; i < (plint) param.yPositions.size(); i++) {
            Array<T,3> llc(param.yxRange[0], param.yPositions[i], param.yzRange[0]);
            Array<T,3> urc(param.yxRange[1], param.yPositions[i], param.yzRange[1]);
            outputCuboids.push_back(Cuboid<T>(llc, urc));
            param.outputDomainNames.push_back(createFileName("slice_y_", i, numYdigits + 1));
        }

        plint numZdigits = util::val2str(param.zPositions.size()).length();
        for (plint i = 0; i < (plint) param.zPositions.size(); i++) {
            Array<T,3> llc(param.zxRange[0], param.zyRange[0], param.zPositions[i]);
            Array<T,3> urc(param.zxRange[1], param.zyRange[1], param.zPositions[i]);
            outputCuboids.push_back(Cuboid<T>(llc, urc));
            param.outputDomainNames.push_back(createFileName("slice_z_", i, numZdigits + 1));
        }
    }

    std::vector<plint> domainExistsInNumLevels(outputCuboids.size(), 0);
    for (plint iLevel = param.minOutputLevel; iLevel <= param.maxOutputLevel; iLevel++) {
        std::vector<Box3D> outputDomainsAtLevel;
        for (plint iCuboid = 0; iCuboid < (plint) outputCuboids.size(); iCuboid++) {
            Box3D outputDomain = computeOutputDomain(param, outputCuboids[iCuboid], param.minOutputLevel);
            if (iLevel != param.minOutputLevel) {
                outputDomain = outputDomain.multiply(util::intTwoToThePower(iLevel - param.minOutputLevel));
                //outputDomain.x1--;
                //outputDomain.y1--;
                //outputDomain.z1--;
            }
            if (outputDomain.getNx() <= 0 || outputDomain.getNy() <= 0 || outputDomain.getNz() <= 0) {
                outputDomain = Box3D(-1, -1, -1, -1, -1, -1);
            } else {
                domainExistsInNumLevels[iCuboid] += 1;
            }
            outputDomainsAtLevel.push_back(outputDomain);
        }
        PLB_ASSERT(param.outputDomainNames.size() == outputDomainsAtLevel.size());
        param.outputDomains[iLevel] = outputDomainsAtLevel;
    }
    for (plint iCuboid = 0; iCuboid < (plint) outputCuboids.size(); iCuboid++) {
        PLB_ASSERT(domainExistsInNumLevels[iCuboid] != 0);
    }
}

void calculateDerivedSimulationParameters(SimulationParameters& param)
{
    // Derived quantities.

    param.smallEnvelopeWidth = 1;
    param.mediumEnvelopeWidth = 2;
    param.largeEnvelopeWidth = 4;
    param.fileNamePadding = 8;
    param.saveDynamicContent = true;

    plint numLevels = param.ogs.getNumLevels();
    PLB_ASSERT(numLevels >= 1);
    param.finestLevel = numLevels - 1;
    if (param.minOutputLevel < 0) {
        param.minOutputLevel = 0;
    }
    if (param.maxOutputLevel > param.finestLevel) {
        param.maxOutputLevel = param.finestLevel;
    }

    param.physicalLocation = Array<T,3>(param.fullDomain.x0(), param.fullDomain.y0(), param.fullDomain.z0());

    param.dxCoarsest = param.dxFinest * (T) util::intTwoToThePower(param.finestLevel);
    param.dtFinest = (param.u_LB / param.u_Ref) * param.dxFinest;
    param.dtCoarsest = param.dtFinest * (T) util::intTwoToThePower(param.finestLevel);

    param.rho_LB = 1.0;
    param.inletVelocity_LB = Array<T,3>(param.inletVelocity * param.dtFinest / param.dxFinest, (T) 0, (T) 0);

    param.omega.resize(numLevels);
    for (plint iLevel = 0; iLevel < numLevels; iLevel++) {
        T dx = param.dxFinest * (T) util::intTwoToThePower(param.finestLevel - iLevel);
        T dt = param.dtFinest * (T) util::intTwoToThePower(param.finestLevel - iLevel);
        T nu_LB = param.nu * dt / (dx * dx);
        param.omega[iLevel] = (T) 1 / (DESCRIPTOR<T>::invCs2 * nu_LB + (T) 0.5);
    }

    FileName fileName(param.staticSurfaceFileName);
    param.surfaceName = fileName.getName();

    computeAllOutputDomains(param);
}

void printSimulationParameters(SimulationParameters const& param)
{
    pcout << "inletVelocity = " << param.inletVelocity << std::endl;
    pcout << "outflowBcType = " << param.outflowBcType << std::endl;

    for (int iZone = 0; iZone < 6; iZone++) {
        pcout << "spongeWidths[" << iZone << "] = " << param.spongeWidths[iZone] << std::endl;
    }

    pcout << "precision = " << (param.precision == FLT ? "FLT" :
            (param.precision == DBL ? "DBL" :
             (param.precision == LDBL ? "LDBL" :
              "INF"))) << std::endl;
    pcout << "L_Ref = " << param.L_Ref << std::endl;
    pcout << "u_Ref = " << param.u_Ref << std::endl;
    pcout << "u_LB = " << param.u_LB << std::endl;
    pcout << "maxIter = " << param.maxIter << std::endl;

    pcout << "rho = " << param.rho << std::endl;
    pcout << "nu = " << param.nu << std::endl;
    pcout << "ambientPressure = " << param.ambientPressure << std::endl;

    pcout << "outDir = " << param.outDir << std::endl;
    pcout << "statIter = " << param.statIter << std::endl;
    pcout << "outIter = " << param.outIter << std::endl;

    pcout << "computeAverages = " << (param.computeAverages ? "true" : "false") << std::endl;
    if (param.computeAverages) {
        pcout << "avgIter = " << param.avgIter << std::endl;
    }

    pcout << "cpIter = " << param.cpIter << std::endl;
    pcout << "abIter = " << param.abIter << std::endl;
    pcout << "abortFileName = " << param.abortFileName << std::endl;
    pcout << "xmlContinueFileName = " << param.xmlContinueFileName << std::endl;
    pcout << "baseFileName = " << param.baseFileName << std::endl;
    pcout << "useParallelIO = " << (param.useParallelIO ? "true" : "false") << std::endl;

    pcout << "finestLevel = " << param.finestLevel << std::endl;

    pcout << "inletVelocity_LB = [" << param.inletVelocity_LB[0] << ", "
                                    << param.inletVelocity_LB[1] << ", "
                                    << param.inletVelocity_LB[2] << "]" << std::endl;
    pcout << "Re = " << param.u_Ref * param.L_Ref / param.nu << std::endl;
    pcout << "omegaFinest = " << param.omega[param.finestLevel] << std::endl;
    pcout << "tauFinest = " << (T) 1 / param.omega[param.finestLevel] << std::endl;
    pcout << "omegaCoarsest = " << param.omega[0] << std::endl;
    pcout << "tauCoarsest = " << (T) 1 / param.omega[0] << std::endl;
    pcout << "dxFinest = " << param.dxFinest << std::endl;
    pcout << "dtFinest / dxFinest = " << param.dtFinest / param.dxFinest << std::endl;
    pcout << "dtFinest / (dxFinest * dxFinest) = " << param.dtFinest / (param.dxFinest * param.dxFinest) << std::endl;
    pcout << "physicalLocation = (" << param.physicalLocation[0] << ", " << param.physicalLocation[1] << ", "
          << param.physicalLocation[2] << ")" << std::endl;
    pcout << "incompressibleModel = " << (param.incompressibleModel ? "true" : "false") << std::endl;
    pcout << std::endl;
}

void createZones(SimulationParameters const& param, MultiBlockLattice3D<T,DESCRIPTOR>& lattice, plint level)
{
    T dx = param.dxFinest * (T) util::intTwoToThePower(param.finestLevel - level);

    Array<plint,6> numSpongeCells;
    plint totalNumSpongeCells = 0;
    for (plint iZone = 0; iZone < 6; iZone++) {
        numSpongeCells[iZone] = util::roundToInt(param.spongeWidths[iZone] / dx);
        totalNumSpongeCells += numSpongeCells[iZone];
    }
    plint nx = util::roundToInt(toLB(param.fullDomain.x1(), 0, dx, param.physicalLocation)) + 1; 
    plint ny = util::roundToInt(toLB(param.fullDomain.y1(), 1, dx, param.physicalLocation)) + 1; 
    plint nz = util::roundToInt(toLB(param.fullDomain.z1(), 2, dx, param.physicalLocation)) + 1; 
    Box3D fullBox(0, nx-1, 0, ny-1, 0, nz-1);

    if (totalNumSpongeCells > 0) {
        pcout << "Generating viscosity sponge zone at level: " << level << std::endl;
        T bulkValue = param.omega[level];

        std::vector<MultiBlock3D*> args;
        args.push_back(&lattice);
        applyProcessingFunctional(new ViscositySpongeZone3D<T,DESCRIPTOR>(nx, ny, nz, bulkValue, numSpongeCells),
                lattice.getBoundingBox(), args);
    }
}

void applyOuterBoundaryConditions(SimulationParameters const& param, 
    MultiLevelCoupling3D<T,DESCRIPTOR,RESCALER>& lattices,
    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* bc)
{
    Box3D coarsestBoundingBox = lattices.getOgs().getClosedCover(0);
    for (plint iLevel = 0; iLevel <= param.finestLevel; iLevel++) {
        pcout << "Generating outer domain boundary conditions at level: " << iLevel << std::endl;
        MultiBlockLattice3D<T,DESCRIPTOR>& lattice = lattices.getLevel(iLevel);

        lattice.periodicity().toggleAll(false);

        Box3D boundingBox = coarsestBoundingBox.multiply(util::intTwoToThePower(iLevel));
        Box3D box = boundingBox; 

        Box3D inlet   (box.x0,   box.x0, box.y0,   box.y1,   box.z0,   box.z1);
        Box3D outlet  (box.x1,   box.x1, box.y0+1, box.y1-1, box.z0+1, box.z1-1);
        Box3D yBottom (box.x0+1, box.x1, box.y0,   box.y0,   box.z0,   box.z1);
        Box3D yTop    (box.x0+1, box.x1, box.y1,   box.y1,   box.z0,   box.z1);
        Box3D zBottom (box.x0+1, box.x1, box.y0+1, box.y1-1, box.z0,   box.z0);
        Box3D zTop    (box.x0+1, box.x1, box.y0+1, box.y1-1, box.z1,   box.z1);

        // Inlet boundary condition.

        bc->setVelocityConditionOnBlockBoundaries(lattice, boundingBox, inlet, boundary::dirichlet);

        Array<T,3> velocity(param.inletVelocity_LB);
        setBoundaryVelocity(lattice, inlet, velocity);

        Array<T,3> zero((T) 0, (T) 0, (T) 0);

        // Lateral boundary conditions.
        bc->setVelocityConditionOnBlockBoundaries(lattice, boundingBox, yBottom, boundary::dirichlet);
        setBoundaryVelocity(lattice, yBottom, velocity);
        
        bc->setVelocityConditionOnBlockBoundaries(lattice, boundingBox, yTop, boundary::dirichlet);
        setBoundaryVelocity(lattice, yTop, velocity);
        
        bc->setVelocityConditionOnBlockBoundaries(lattice, boundingBox, zBottom, boundary::dirichlet);
        setBoundaryVelocity(lattice, zBottom, velocity);
        
        bc->setVelocityConditionOnBlockBoundaries(lattice, boundingBox, zTop, boundary::dirichlet);
        setBoundaryVelocity(lattice, zTop, velocity);
    
        // Outlet boundary condition.

        if (param.outflowBcType == 0) {
            bc->setVelocityConditionOnBlockBoundaries(lattice, boundingBox, outlet, boundary::dirichlet);
            setBoundaryVelocity(lattice, outlet, velocity);
        } else if (param.outflowBcType == 1) {
            bc->setPressureConditionOnBlockBoundaries(lattice, boundingBox, outlet, boundary::dirichlet);
            setBoundaryVelocity(lattice, outlet, velocity);
        } else if (param.outflowBcType == 2) {
            bc->setVelocityConditionOnBlockBoundaries(lattice, boundingBox, outlet, boundary::neumann);
            setBoundaryVelocity(lattice, outlet, velocity);
        }

        setBoundaryDensity(lattice, box, param.rho_LB);
    }
}

void initializeSimulation(SimulationParameters const& param, bool continueSimulation, plint& iniIter,
        MultiLevelCoupling3D<T,DESCRIPTOR,RESCALER>& lattices,
        std::vector<MultiBlock3D*>& checkpointBlocks)
{
    if (!continueSimulation) {
        for (plint iLevel = 0; iLevel <= param.finestLevel; iLevel++) {

            Array<T,3> velocity(param.inletVelocity_LB);

            MultiBlockLattice3D<T,DESCRIPTOR>& lattice = lattices.getLevel(iLevel);
            initializeAtEquilibrium(lattice, lattice.getBoundingBox(), param.rho_LB, velocity);
        }

        // We do NOT want to call internal data processors here...
        lattices.initialize();
        lattices.initializeTensorFields();
    } else {
        pcout << std::endl;
        pcout << "Reading state of the simulation from file: " << param.xmlContinueFileName << std::endl;
        loadState(checkpointBlocks, iniIter, param.saveDynamicContent, param.xmlContinueFileName);
        for (plint iLevel = 0; iLevel <= param.finestLevel; iLevel++) {
            MultiBlockLattice3D<T,DESCRIPTOR>& lattice = lattices.getLevel(iLevel);
            plint iniIterAtLevel = iniIter * util::intTwoToThePower(iLevel);
            lattice.resetTime(iniIterAtLevel);
        }

        lattices.initializeTensorFields();
        //lattices.initialize();
    }

    pcout << std::endl;
}

void writeResults(SimulationParameters const& param, 
    MultiLevelCoupling3D<T,DESCRIPTOR,RESCALER>& lattices,
    std::auto_ptr<MultiLevelTensorField3D<T,3> > &avgVel, 
    plint iter)
{
    bool crop = true;
    for (plint iDomain = 0; iDomain < (plint) param.outputDomainNames.size(); iDomain++) {
        std::string fname = createFileName(param.outputDomainNames[iDomain] + "_", iter, param.fileNamePadding);
        SparseVtkImageOutput3D sparseOut(fname);

        std::auto_ptr<MultiLevelTensorFieldForOutput3D<T,3> > velocity =
            computeVelocity(lattices, param.outputDomains.find(param.maxOutputLevel)->second[iDomain],
                    param.maxOutputLevel, crop);

        std::auto_ptr<MultiLevelScalarFieldForOutput3D<T> > density = 
            computeDensity(lattices, param.outputDomains.find(param.maxOutputLevel)->second[iDomain], 
                param.maxOutputLevel, crop);

        std::auto_ptr<MultiLevelTensorFieldForOutput3D<T,3> > outAvgVel;

        if (param.computeAverages) {
            outAvgVel = exportForOutput(
                    *extractSubDomain(*avgVel, param.outputDomains.find(param.maxOutputLevel)->second[iDomain],
                        param.maxOutputLevel), 
                    param.outputDomains.find(param.maxOutputLevel)->second[iDomain],
                    param.maxOutputLevel, crop);

        }

        for (plint iLevel = param.minOutputLevel; iLevel <= param.maxOutputLevel; iLevel++) {
            T dx = param.dxFinest * util::intTwoToThePower(param.finestLevel - iLevel);
            T dt = param.dtFinest * util::intTwoToThePower(param.finestLevel - iLevel);

            T pressureScale = param.rho * (dx * dx) / (dt * dt) * DESCRIPTOR<T>::cs2;
            T pressureOffset = param.ambientPressure - param.rho_LB * pressureScale;
            
            Group3D vtkGroup;
            // You can add scalar- and tensor-fields to the group with the usual "group.add()". The advantage of addTransform
            // is that is also converts the type to float and multiplies by a scale factor (and adds an offset).
            addTransform<T,float,3>(vtkGroup, velocity->getLevel(iLevel), "velocity", dx / dt);
            addTransform<T,float>(vtkGroup, density->getLevel(iLevel), "pressure", pressureScale, pressureOffset);
            if (param.computeAverages) {
                addTransform<T,float,3>(vtkGroup, outAvgVel->getLevel(iLevel), "avgVel", dx / dt);
            }

            // "pointData = true" is the usual VTK output. "pointData = false" inhibits interpolations, and is useful for debugging.
            bool pointData = false;
            sparseOut.writeVtkBlock(vtkGroup, dx, param.physicalLocation, iLevel, pointData);
            // ... Add more data at different levels.

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
        pcerr << "Usage: " << argv[0] << " xml-input-file-name [restart]" << std::endl;
        exit(1);
    }

    std::string xmlInputFileName;
    xmlInputFileName = std::string(argv[1]);
    abortIfCannotOpenFileForReading(xmlInputFileName);

    bool continueSimulation = false;
    if (argc == 3) {
        std::string cmd(argv[2]);
        if (cmd == "restart") {
            continueSimulation = true;
        }
    }

    int nproc = global::mpi().getSize();

    //global::mpi().barrier();
    global::timer("init").start();

    // Set the simulation parameters.

    SimulationParameters param;

    readUserDefinedSimulationParameters(xmlInputFileName, param);

    if (continueSimulation) {
        abortIfCannotOpenFileForReading(param.xmlContinueFileName);
    }

    createOctreeGridStructure(param);
    calculateDerivedSimulationParameters(param);

    global::directories().setOutputDir(param.outDir);
    global::IOpolicy().activateParallelIO(param.useParallelIO);

    if (nproc != param.ogs.getNumProcesses()) {
        pcerr << "The number of processes used is not the same as the one provided in the grid-structure files." << std::endl;
        exit(1);
    }

    // The "order" is about how Palabos rescales the populations. For now, we work at order 0.
    // RESCALER means: we use convective scaling,
    // Palabos scales only the populations and no external scalars, and
    // this works for BGK but not for MRT.
    plint order = 0;

    pcout << std::endl;
    pcout << "Generating the lattices." << std::endl;

    order = 1;
    MultiLevelCoupling3D<T,DESCRIPTOR,RESCALER> lattices(
            param.ogs, new ConsistentSmagorinskyCompleteRegularizedBGKdynamics<T,DESCRIPTOR>(
                    param.omega[0], 0.14 ),
            order);   
    pcout << "CompleteRegularizedBGKdynamics" << std::endl;
    param.incompressibleModel = false;

    for (plint iLevel = 0; iLevel < lattices.getNumLevels(); iLevel++) {
        pcout << "Info for lattice at level: " << iLevel << std::endl;
        pcout << getMultiBlockInfo(lattices.getLevel(iLevel)) << std::endl;
    }

    printSimulationParameters(param);


    // Immersed surfaces only for finest level.

    pcout << "Reading the immersed surface geometries." << std::endl;
    pcout << "Generating fluid blocks at finest level." << std::endl;

    std::vector<MultiBlock3D*> rhoBarJarg;
    // Here we assume that the object is far from the outlet. This is why we can compute
    // the rhoBarJ field to be used with the off-lattice BC, before the FluidPressureOutlet3D
    // is executed.
    plint numScalars = 4;
    plint extendedTensorEnvelopeWidth = 2;  // Extrapolated BCs. guo = 2, generalized = 3
    MultiNTensorField3D<T> *rhoBarJfield = 
        generateMultiNTensorField3D<T>(lattices.getLevel(param.finestLevel), extendedTensorEnvelopeWidth, numScalars);
    rhoBarJfield->toggleInternalStatistics(false);
    rhoBarJarg.push_back(rhoBarJfield);
    integrateProcessingFunctional(new PackedRhoBarJfunctional3D<T,DESCRIPTOR>(),
            lattices.getLevel(param.finestLevel).getBoundingBox(), lattices.getLevel(param.finestLevel), *rhoBarJfield, 0);

    pcout << "Implementing the geometry at level." << std::endl;

    // The next few lines of code are typical. They transform the surface geometry of the
    //   stl given by the user to more efficient data structures that are internally
    //   used by palabos. The TriangleBoundary3D structure will be later used to assign
    //   proper boundary conditions.
    TriangleSet<T> triangleSet(param.staticSurfaceFileName, DBL);
    triangleSet.translate(-param.physicalLocation);
    triangleSet.scale((T)1/param.dxFinest);  

    DEFscaledMesh<T> defMesh(triangleSet, 0, 0, 1, Dot3D(0,0,0));
    defMesh.setDx(param.dxFinest);
    defMesh.setPhysicalLocation(param.physicalLocation);
    TriangleBoundary3D<T> boundary(defMesh);
    boundary.getMesh().inflate(0.05);


    // The aneurysm simulation is an interior (as opposed to exterior) flow problem. For
    // this reason, the lattice nodes that lay inside the computational domain must
    // be identified and distinguished from the ones that lay outside of it. This is
    // handled by the following voxelization process.
    const int flowType = voxelFlag::outside;
    const int borderWidth = 1;
    const int extendedEnvelopeWidth = 2;
    const int blockSize = 0;
    VoxelizedDomain3D<T> voxelizedDomain (
            boundary, flowType, lattices.getLevel(param.finestLevel).getBoundingBox(), borderWidth, extendedEnvelopeWidth, blockSize );
    voxelizedDomain.reparallelize(param.ogs.getMultiBlockManagement(
                param.finestLevel, lattices.getLevel(param.finestLevel).getBoundingBox(),extendedEnvelopeWidth));

    defineDynamics(lattices.getLevel(param.finestLevel), 
               voxelizedDomain.getVoxelMatrix(), 
               lattices.getLevel(param.finestLevel).getBoundingBox(),
               new NoDynamics<T,DESCRIPTOR>((T) 1), voxelFlag::inside); 

    boundary.getMesh().writeAsciiSTL(param.outDir+"flowMesh.stl");

    pcout << getMultiBlockInfo(voxelizedDomain.getVoxelMatrix()) << std::endl;
    // The boundary condition algorithm for the object.
    BoundaryProfiles3D<T,Velocity> profiles;
    profiles.setWallProfile(new NoSlipProfile3D<T>());

    // Filipova BC needs no dynamics inside the voxel object
    defineDynamics(lattices.getLevel(param.finestLevel), 
           voxelizedDomain.getVoxelMatrix(), 
           lattices.getLevel(param.finestLevel).getBoundingBox(),
           new NoDynamics<T,DESCRIPTOR>((T) 1), voxelFlag::innerBorder);     

    // Filippova boundary condition
    FilippovaHaenelModel3D<T,DESCRIPTOR>* model = new FilippovaHaenelModel3D<T,DESCRIPTOR>(
            new TriangleFlowShape3D<T,Array<T,3> >(voxelizedDomain.getBoundary(), profiles),
            flowType);

    OffLatticeBoundaryCondition3D<T,DESCRIPTOR,Velocity>* boundaryCondition = 
        new OffLatticeBoundaryCondition3D<T,DESCRIPTOR,Velocity>(
            model->clone(), voxelizedDomain, lattices.getLevel(param.finestLevel));
    boundaryCondition->insert(rhoBarJarg);
    PLB_ASSERT(boundaryCondition != 0);


    if (!continueSimulation) {
        pcout << "Generating outer domain zones." << std::endl;
        for (plint iLevel = 0; iLevel <= param.finestLevel; iLevel++) {
            createZones(param, lattices.getLevel(iLevel), iLevel);
        }
    }

    pcout << "Generating outer domain boundary conditions." << std::endl;

    OnLatticeBoundaryCondition3D<T,DESCRIPTOR> *bc = 
        createInterpBoundaryCondition3D<T,DESCRIPTOR>();
    applyOuterBoundaryConditions(param, lattices, bc);
    delete bc; bc = 0;

    // Generation of statistics MultiLevel fields.

    std::auto_ptr<MultiLevelTensorField3D<T,3> > vel        ;
    std::auto_ptr<MultiLevelTensorField3D<T,3> > avgVelocity;

    std::vector<MultiLevel3D *> avgVelocityArgs; 

    Box3D coarsestBoundingBox = lattices.getOgs().getClosedCover(0);
    if (param.computeAverages) {
        vel         = generateMultiLevelTensorField3D<T,3>(lattices.getOgs(), coarsestBoundingBox, 0, Array<T,3>(0.0,0.0,0.0));
        avgVelocity = generateMultiLevelTensorField3D<T,3>(lattices.getOgs(), coarsestBoundingBox, 0, Array<T,3>(0.0,0.0,0.0));

        avgVelocityArgs.push_back(vel.get());
        avgVelocityArgs.push_back(avgVelocity.get());
    }

    // Initialization.

    std::vector<MultiBlock3D*> checkpointBlocks;

    for (plint iLevel = 0; iLevel <= param.finestLevel; iLevel++) {
        checkpointBlocks.push_back(&lattices.getLevel(iLevel));
    }

    if (param.computeAverages) {
        for (plint iLevel = 0; iLevel <= param.finestLevel; iLevel++) {
            checkpointBlocks.push_back(&avgVelocity->getLevel(iLevel));
        }
    }

    plint iniIter = 0;

    initializeSimulation(param, continueSimulation, iniIter, 
            lattices, checkpointBlocks);

    // Use "collideAndStream" at all levels except the finest one, at 
    // which "executeInternalProcessors" is used instead.
    std::map<plint, bool> useExecuteInternalProcessors;
    for (plint iLevel = 0; iLevel <= param.finestLevel; iLevel++) {
        useExecuteInternalProcessors[iLevel] = false;
    }

    // Prepare files.
    std::string fileName = param.outDir + "average_energy_finest_level.dat";
    plb_ofstream energyFile(fileName.c_str(), continueSimulation ? std::ofstream::app : std::ofstream::out);

    // integration of point measures and creation of output files if needed
    std::vector<std::vector<plint> > ids;
    std::vector<plb_ofstream *> probesFileNames(param.finestLevel+1);
    std::vector<std::vector<std::vector<T> > > results(param.finestLevel+1);
    const plint statsId = -200;

    bool computeStats = false; 

    
    fileName = param.outDir + param.surfaceName + "_total_force.dat";
    plb_ofstream forces(fileName.c_str(), continueSimulation ? std::ofstream::app : std::ofstream::out);

    // Starting iterations.
    global::timer("init").stop();

    pcout << "The full initialization phase took " << global::timer("init").getTime() << " seconds on "
          << nproc << " processes." << std::endl;
    pcout << std::endl;

    pcout << "Starting simulation." << std::endl;
    pcout << std::endl;
    bool stopExecution = false;
    bool checkForErrors = true;
    plint iter = iniIter;
    bool avgProcIntegrated = false;
    std::vector<plint> extProcFunIds;
    for (; iter < param.maxIter && !stopExecution; iter++) {
        if (iter % param.statIter == 0 && iter != 0) {
            pcout << "At coarsest level iteration: " << iter << ", t = " << iter * param.dtCoarsest << std::endl;
            T energy = boundaryCondition->computeAverageEnergy() *
                param.rho * (param.dxFinest * param.dxFinest) / (param.dtFinest * param.dtFinest);

            pcout << "Average kinetic energy at the finest level: " << energy << std::endl;
            energyFile << (double) (iter * param.dtCoarsest) << " " << (double) energy << std::endl;

            // Forces on immersed surfaces.

            T forceConversion = param.rho * (param.dxFinest * param.dxFinest * param.dxFinest * param.dxFinest) /
                (param.dtFinest * param.dtFinest);
            Array<T,3> force = forceConversion*boundaryCondition->getForceOnObject();
            forces << (double) (iter * param.dtCoarsest) << " " << force[0] << " " << force[1] << " " << force[2] << " " << empirical_sphere_drag(param.Re) << std::endl;
            
            if (iter > 0) {
                T totTime = global::timer("lb-iter").getTime();
                pcout << "Total coarsest level iteration: " << totTime / (T) iter << std::endl;
            }

            pcout << std::endl;
        }   

        // With grid-refinement, the order of the integration of data processors is different.
        // ExternalRhoJcollideAndStream3D is integrated at the end, and not at the beginning
        // of the cycle, as it happens in codes that do not use grid-refinement. This means
        // that in these cases the boundary conditions are imposed "before" collide-and-stream
        // in the cycle. So, the results are exported before the BCs are enforced, and this
        // might cause visualization problems.
        if (iter % param.outIter == 0) {
            pcout << "Output to disk at coarsest level iteration: " << iter << ", t = " << iter * param.dtCoarsest << std::endl;

            writeResults(param, lattices, avgVelocity, iter);

        }

        if ((param.cpIter > 0 && iter % param.cpIter == 0 && iter != iniIter) || iter == param.maxIter - 1) {
            pcout << "Saving the state of the simulation at coarsest level iteration: " << iter << std::endl;
            saveState(checkpointBlocks, iter, param.saveDynamicContent, param.xmlContinueFileName,
                    param.baseFileName, param.fileNamePadding);
            pcout << std::endl;
        }

        if (iter % param.abIter == 0) {
            stopExecution = abortExecution(param.abortFileName, checkpointBlocks, iter,
                    param.saveDynamicContent, param.xmlContinueFileName,
                    param.baseFileName, param.fileNamePadding);

            if (stopExecution) {
                pcout << "Aborting execution at iteration: " << iter << std::endl;
                pcout << std::endl;
            }
        }

        global::timer("lb-iter").start();
        global::timer("solver").start();
        lattices.collideAndStream(0, useExecuteInternalProcessors, extProcFunIds, computeStats, statsId, ids, results);  
        global::timer("solver").stop();
        global::timer("lb-iter").stop();

        if (checkForErrors) {
            abortIfErrorsOccurred();
            checkForErrors = false;
        }

        if (param.computeAverages && (iter == param.avgIter || (iter > param.avgIter && continueSimulation && !avgProcIntegrated))) {
            avgProcIntegrated = true;
            plint compFields = -100;
            plint compAvgs = compFields - 1;
            integrateProcessingFunctional(
                new BoxVelocityFunctional3D<T,DESCRIPTOR>(), coarsestBoundingBox, 0, lattices, *vel, compFields);
            integrateProcessingFunctional(
                new UpdateAveTensorTransientStatistics3D<T,3>(iter-param.avgIter+1), coarsestBoundingBox, 0, lattices, avgVelocityArgs, param.ogs.getNumLevels(), compAvgs);

            extProcFunIds.push_back(compFields);
            extProcFunIds.push_back(compAvgs);
        }
    }

    pcout << "The " << iter - iniIter << " iterations at the coarsest level, took " << global::timer("solver").getTime()
        << " seconds on " << nproc << " processes." << std::endl;

    plb_ofstream summary("execution_summary.txt");
    summary << "Summary of execution of the solver: " << argv[0] << std::endl;
    summary << "Number of processes: " << nproc << std::endl;
    summary << "Total time of the initialization phase           : " << (double) global::timer("init").getTime() << " s" << std::endl;
    summary << "Total time of the pure solution phase (no output): " << (double) global::timer("solver").getTime() << " s" << std::endl;
    summary.close();
    energyFile.close();

    return 0;
}
