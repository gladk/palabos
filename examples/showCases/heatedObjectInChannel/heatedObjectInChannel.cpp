
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
 * Internal flow around a 3D object with convective heat transfer.
 * This example demonstrates many features of Palabos:
 * Efficiently solving an advection-diffusion equation, coupled with the fluid.
 * Configuration with an input XML file.
 * Loading geometries from STL files.
 * Using the voxelizer.
 * Using off-lattice boundary conditions.
 * Imposing sophisticated outflow boundary conditions.
 * Imposing time dependent analytically defined inlet boundary conditions.
 * Imposing an analytically defined initial condition.
 * Checkpointing (saving the state of the simulation and restarting).
 * Parallel I/O.
 * */

/* 
 * Description and usage:
 *
 * This code solves the flow with heat transfer around a 3D object inside a
 * rectangular channel.
 *
 * The code needs to be compiled only once, and then subsequent runs are configured
 * completely through the input XML file. The user can provide his own geometry
 * as an STL file. There are two demo input XML files: "steady.xml" and "unsteady.xml".
 *
 * To run the code with, say 4, processes use:
 *
 * mpirun -np 4 ./heatedObjectInChannel steady.xml
 *
 * The input XML files contain many comments on how to configure the simulation.
 * To stop and save the state of the simulation, just create a file called "abort"
 * in the working directory. To restart the simulation use:
 *
 * mpirun -np 4 ./heatedObjectInChannel steady.xml continue.xml
 *
 * The file "continue.xml" is generated automatically.
 *
 */


#include "palabos3D.h"
#include "palabos3D.hh"

#include <cmath>


using namespace plb;

typedef double T;
typedef Array<T,3> Velocity;


#define NMAX 150
#define PADD 8

#define DESCRIPTOR descriptors::D3Q19Descriptor
#define TEMPERATURE_DESCRIPTOR descriptors::AdvectionDiffusionD3Q7Descriptor
//#define TEMPERATURE_DESCRIPTOR descriptors::AdvectionDiffusionD3Q19Descriptor


static plint xDirection = 0;

static plint borderWidth           = 1;  // Because Guo acts in a one-cell layer.
                                         // Requirement: margin>=borderWidth.
static plint margin                = 1;  // Extra margin of allocated cells around the object.
static plint blockSize             = 0;  // Size of blocks in the sparse/parallel representation.
static plint envelopeWidth         = 1;
static plint extendedEnvelopeWidth = 2;  // Extrapolated BCs.

static T pi = std::acos((T) -1);

static std::string outDir("./tmp/");


// This structure holds all the user-defined parameters, and some
// derived values needed for the simulation.
struct Param {
    std::string file;               // STL file with the object geometry.
    T inflationParameter;           // Parameter to inflate the geometry before voxelization.
    Precision precision;            // Precision for geometric operations.
    T cx, cy, cz;                   // Position of the center of the object.
    T lx, ly, lz;                   // Size of computational domain.

    T uInf;                         // x velocity at the inflow.
    T nu;                           // Kinematic viscosity.

    T Prandtl;                      // Prandtl number.
    plint temperatureTimeFactor;    // Ratio between the temperature and the fluid time steps.
    T T_inf;                        // Ambient temperature.
    T T_wall;                       // Wall temperature.

    T characteristicLength;         // Length to define the resolution.
    plint resolution;               // Number of lattice nodes along the characteristicLength.
    T uLB;                          // Velocity in lattice units.
    T cSmago;                       // Parameter for the Smagorinsky LES model.

    plint maxIter, startIter;       // Time for events in lattice units.
    plint statIter, vtkIter;
    plint cpIter;
    bool useRhoBarJ;                // Use a rhoBar-j formulation or not?
    bool useParallelIO;             // Use MPI I/O or not?

    std::string abortFile;          // File to signal program abortion.
    std::string continueFile;       // File to continue the simulation.
    std::string checkpointFile;     // Base filename for the checkpoint files.
    bool saveDynamicContent;        // Save dynamic content in the checkpoint files or not?
    plint abIter;                   // Frequency to check for user-driven program abortion.

    plint nx, ny, nz;               // Grid resolution of bounding box.
    T nuLB;                         // Kinematic viscosity in lattice units.
    T omega;                        // Relaxation parameter for the fluid.
    T omega_temperature;            // Relaxation parameter for the temperature.
    T dx, dt;                       // Discrete space and time steps.

    Param()
    { }
    
    Param(std::string xmlFname)
    {
        XMLreader document(xmlFname);

        document["geometry"]["file"].read(file);
        document["geometry"]["inflationParameter"].read(inflationParameter);
        std::string precisionStr;
        document["geometry"]["precision"].read(precisionStr);
        PLB_ASSERT(precisionStr == "FLT" || precisionStr == "DBL" || precisionStr == "LDBL" || precisionStr == "INF");
        if (precisionStr == "FLT") {
            precision = FLT;
        } else if (precisionStr == "DBL") {
            precision = DBL;
        } else if (precisionStr == "LDBL") {
            precision = LDBL;
        } else {
            precision = INF;
        }
        document["geometry"]["center"]["x"].read(cx);
        document["geometry"]["center"]["y"].read(cy);
        document["geometry"]["center"]["z"].read(cz);
        document["geometry"]["domain"]["x"].read(lx);
        document["geometry"]["domain"]["y"].read(ly);
        document["geometry"]["domain"]["z"].read(lz);

        document["fluid"]["uInf"].read(uInf);
        document["fluid"]["nu"].read(nu);

        document["temperature"]["Prandtl"].read(Prandtl);
        document["temperature"]["temperatureTimeFactor"].read(temperatureTimeFactor);
        document["temperature"]["T_inf"].read(T_inf);
        document["temperature"]["T_wall"].read(T_wall);

        document["numerics"]["characteristicLength"].read(characteristicLength);
        document["numerics"]["resolution"].read(resolution);
        document["numerics"]["uLB"].read(uLB);
        document["numerics"]["cSmago"].read(cSmago);

        document["simulation"]["maxIter"].read(maxIter);
        document["simulation"]["startIter"].read(startIter);
        document["simulation"]["statIter"].read(statIter);
        document["simulation"]["vtkIter"].read(vtkIter);
        document["simulation"]["cpIter"].read(cpIter);
        document["simulation"]["useRhoBarJ"].read(useRhoBarJ);
        document["simulation"]["useParallelIO"].read(useParallelIO);

        abortFile = "abort";
        continueFile = "continue.xml";
        checkpointFile = "checkpoint_";
        saveDynamicContent = true;
        abIter = 16;

        computeLBparameters();
    }

    void computeLBparameters()
    {
        T uInfNorm = std::fabs(uInf);
        dx = characteristicLength / (T) (resolution - 1);
        dt = uLB / uInfNorm * dx;
        nuLB = nu * dt / (dx * dx);
        omega = (T) 1 / (DESCRIPTOR<T>::invCs2 * nuLB + (T) 0.5);
        nx = util::roundToInt(lx / dx) + 1;
        ny = util::roundToInt(ly / dx) + 1;
        nz = util::roundToInt(lz / dx) + 1;

        T k = nuLB / Prandtl;
        omega_temperature = (T) 1 / (TEMPERATURE_DESCRIPTOR<T>::invCs2 * k + (T) 0.5);
    }

    Box3D bbox()
    {
        return Box3D(0,    nx-1, 0,    ny-1, 0,    nz-1);
    }
    Box3D inlet()
    {
        return Box3D(0,    0,    1,    ny-2, 1,    nz-2);
    }
    Box3D outlet()
    {
        return Box3D(nx-1, nx-1, 1,    ny-2, 1,    nz-2);
    }
    Box3D bottom()
    {
        return Box3D(0,    nx-1, 0,    ny-1, 0,    0   );
    }
    Box3D top()
    {
        return Box3D(0,    nx-1, 0,    ny-1, nz-1, nz-1);
    }
    Box3D lateral1()
    {
        return Box3D(0,    nx-1, 0,    0,    1,    nz-2);
    }
    Box3D lateral2()
    {
        return Box3D(0,    nx-1, ny-1, ny-1, 1,    nz-2);
    }

    Box3D vtkSliceX()
    {
        plint cxLB = util::roundToInt(cx / dx);
        return Box3D(cxLB-1, cxLB+1, 0, ny-1, 0, nz-1);
    }
    Box3D vtkSliceY()
    {
        plint cyLB = util::roundToInt(cy / dx);
        return Box3D(0, nx-1, cyLB-1, cyLB+1, 0, nz-1);
    }
    Box3D vtkSliceZ()
    {
        plint czLB = util::roundToInt(cz / dx);
        return Box3D(0, nx-1, 0, ny-1, czLB-1, czLB+1);
    }
};


static Param param;


static T poiseuillePressure(plint maxN)
{
    T a = param.ny - 1;
    T b = param.nz - 1;

    T nu = param.nuLB;
    T uMax = param.uLB;

    T sum = T();
    for (plint iN = 0; iN < maxN; iN += 2) {
        T twoNplusOne = (T) 2 * (T) iN + (T) 1;
        sum += ((T) 1 / (std::pow(twoNplusOne, (T) 3) * std::cosh(twoNplusOne * pi * b / ((T) 2 * a))));
    }
    for (plint iN = 1; iN < maxN; iN += 2) {
        T twoNplusOne = (T) 2 * (T) iN + (T) 1;
        sum -= ((T) 1 / (std::pow(twoNplusOne, (T) 3) * std::cosh(twoNplusOne * pi * b / ((T) 2 * a))));
    }

    T alpha = - (T) 8 * uMax * pi * pi * pi / (a * a * (pi * pi * pi - (T) 32 * sum)); // alpha = -dp/dz / mu

    T deltaP = - (alpha * nu);

    return deltaP;
}

T poiseuilleVelocity(plint iY, plint iZ, plint maxN)
{
    T a = param.ny - 1;
    T b = param.nz - 1;

    T y = (T) iY - a / (T) 2;
    T z = (T) iZ - b / (T) 2;

    T alpha = - poiseuillePressure(maxN) / param.nuLB;

    T sum = T();
    for (plint iN = 0; iN < maxN; iN += 2) {
        T twoNplusOne = (T) 2 * (T) iN + (T) 1;

        sum += (std::cos(twoNplusOne * pi * y / a) * std::cosh(twoNplusOne * pi * z / a) /
               (std::pow(twoNplusOne, (T) 3)       * std::cosh(twoNplusOne * pi * b / ((T) 2 * a))));
    }
    for (plint iN = 1; iN < maxN; iN += 2) {
        T twoNplusOne = (T) 2 * (T) iN + (T) 1;

        sum -= (std::cos(twoNplusOne * pi * y / a) * std::cosh(twoNplusOne * pi * z / a) /
               (std::pow(twoNplusOne, (T) 3)       * std::cosh(twoNplusOne * pi * b / ((T) 2 * a))));
    }

    sum *= ((T) 4 * alpha * a *a / std::pow(pi, (T) 3));
    sum += (alpha / (T) 2 * (y * y - a * a / (T) 4));

    return sum;
}

template <typename T>
class SquarePoiseuilleDensityAndVelocity {
public:
    SquarePoiseuilleDensityAndVelocity(T C_, plint maxN_)
        : C(C_),
          maxN(maxN_)
    { }
    void operator()(plint iX, plint iY, plint iZ, T &rho, Array<T,3>& u) const
    {
        rho = (T) 1;
        u[0] = C * poiseuilleVelocity(iY, iZ, maxN);
        u[1] = T();
        u[2] = T();
    }
private:
    T C;
    plint maxN;
};

template <typename T>
class SquarePoiseuilleVelocity {
public:
    SquarePoiseuilleVelocity(T C_, plint maxN_)
        : C(C_),
          maxN(maxN_)
    { }
    void operator()(plint iX, plint iY, plint iZ, Array<T,3>& u) const
    {
        u[0] = C * poiseuilleVelocity(iY, iZ, maxN);
        u[1] = T();
        u[2] = T();
    }
private:
    T C;
    plint maxN;
};

// Instantiate the boundary conditions for the outer domain, for the fluid.
void outerDomainBoundaries(MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
        OnLatticeBoundaryCondition3D<T,DESCRIPTOR>& boundaryCondition)
{
    Array<T,3> zero((T) 0, (T) 0, (T) 0);
    T C = util::sinIncreasingFunction<T>(0, param.startIter);

    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, param.inlet(), boundary::dirichlet);
    setBoundaryVelocity(lattice, param.inlet(), SquarePoiseuilleVelocity<T>(C, NMAX));

    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, param.bottom(),   boundary::dirichlet);
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, param.top(),      boundary::dirichlet);
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, param.lateral1(), boundary::dirichlet);
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, param.lateral2(), boundary::dirichlet);
    setBoundaryVelocity(lattice, param.bottom(),   zero);
    setBoundaryVelocity(lattice, param.top(),      zero);
    setBoundaryVelocity(lattice, param.lateral1(), zero);
    setBoundaryVelocity(lattice, param.lateral2(), zero);

    integrateProcessingFunctional(new FluidPressureOutlet3D<T,DESCRIPTOR,0,+1>(), param.outlet(), lattice, 0);
    setBoundaryVelocity(lattice, param.outlet(), SquarePoiseuilleVelocity<T>(C, NMAX));
}

// Instantiate the boundary conditions for the outer domain, for the temperature.
void outerTemperatureBoundaries(MultiBlockLattice3D<T,TEMPERATURE_DESCRIPTOR>& lattice)
{
    instantiateOuterNLDboundary(lattice, lattice.getBoundingBox());

    setAD_NLDboundaryDynamics(lattice, param.inlet(),    boundary::dirichlet);
    setAD_NLDboundaryDynamics(lattice, param.bottom(),   boundary::dirichlet);
    setAD_NLDboundaryDynamics(lattice, param.top(),      boundary::dirichlet);
    setAD_NLDboundaryDynamics(lattice, param.lateral1(), boundary::dirichlet);
    setAD_NLDboundaryDynamics(lattice, param.lateral2(), boundary::dirichlet);
    setBoundaryDensity(lattice, param.inlet(),    param.T_inf);
    setBoundaryDensity(lattice, param.bottom(),   param.T_inf);
    setBoundaryDensity(lattice, param.top(),      param.T_inf);
    setBoundaryDensity(lattice, param.lateral1(), param.T_inf);
    setBoundaryDensity(lattice, param.lateral2(), param.T_inf);

    setAD_NLDboundaryDynamics(lattice, param.outlet(), boundary::neumann);
    setBoundaryDensity(lattice, param.outlet(), param.T_inf);
}

// Write VTK file for the flow around the object, to be viewed with Paraview.
void writeVTK(OffLatticeBoundaryCondition3D<T,DESCRIPTOR,Velocity>& fbc,
        OffLatticeBoundaryCondition3D<T,TEMPERATURE_DESCRIPTOR,Array<T,2> >& tbc,
        plint iT, std::string appStr="")
{
    VtkImageOutput3D<T> vtkOut_x(createFileName(appStr+"x_slice_", iT, PADD), param.dx);
    vtkOut_x.writeData<float>(*fbc.computeVelocityNorm(param.vtkSliceX()), "velocityNorm", param.dx / param.dt);
    vtkOut_x.writeData<3,float>(*fbc.computeVelocity(param.vtkSliceX()), "velocity", param.dx / param.dt);
    vtkOut_x.writeData<float>(*fbc.computePressure(param.vtkSliceX()), "pressure",
            param.dx * param.dx / (param.dt * param.dt));
    vtkOut_x.writeData<float>(*tbc.computeDensity(param.vtkSliceX(), param.T_wall), "temperature", 1.0);

    VtkImageOutput3D<T> vtkOut_y(createFileName(appStr+"y_slice_", iT, PADD), param.dx);
    vtkOut_y.writeData<float>(*fbc.computeVelocityNorm(param.vtkSliceY()), "velocityNorm", param.dx / param.dt);
    vtkOut_y.writeData<3,float>(*fbc.computeVelocity(param.vtkSliceY()), "velocity", param.dx / param.dt);
    vtkOut_y.writeData<float>(*fbc.computePressure(param.vtkSliceY()), "pressure",
            param.dx * param.dx / (param.dt * param.dt));
    vtkOut_y.writeData<float>(*tbc.computeDensity(param.vtkSliceY(), param.T_wall), "temperature", 1.0);

    VtkImageOutput3D<T> vtkOut_z(createFileName(appStr+"z_slice_", iT, PADD), param.dx);
    vtkOut_z.writeData<float>(*fbc.computeVelocityNorm(param.vtkSliceZ()), "velocityNorm", param.dx / param.dt);
    vtkOut_z.writeData<3,float>(*fbc.computeVelocity(param.vtkSliceZ()), "velocity", param.dx / param.dt);
    vtkOut_z.writeData<float>(*fbc.computePressure(param.vtkSliceZ()), "pressure",
            param.dx * param.dx / (param.dt * param.dt));
    vtkOut_z.writeData<float>(*tbc.computeDensity(param.vtkSliceZ(), param.T_wall), "temperature", 1.0);
}

void run(std::string continueFileName)
{
    /*
     * Geometry processing.
     */

    pcout << std::endl << "Processing the object geometry." << std::endl;
    Array<T,3> center(param.cx, param.cy, param.cz);
    Array<T,3> centerLB(center / param.dx);
    TriangleSet<T>* triangleSet = new TriangleSet<T>(param.file, param.precision);

    // Place the object in the correct place in the simulation domain.
    // Here the "geometric center" of the object is computed manually,
    // by computing first its bounding cuboid. In cases that the STL
    // file with the geometry of the object contains its center as
    // the point, say (0, 0, 0), then the following variable
    // "objectCenter" must be set to (0, 0, 0) manually.
    Cuboid<T> bCuboid = triangleSet->getBoundingCuboid();
    Array<T,3> objectCenter = (T) 0.5 * (bCuboid.lowerLeftCorner + bCuboid.upperRightCorner);
    triangleSet->translate(-objectCenter);
    triangleSet->scale((T) 1 / param.dx); // In lattice units from now on...
    triangleSet->translate(centerLB);
    triangleSet->writeBinarySTL(outDir+"object_LB.stl");

    // The DEFscaledMesh, and the triangle-boundary are more sophisticated data
    // structures used internally by Palabos to treat the boundary.
    DEFscaledMesh<T>* defMesh = new DEFscaledMesh<T>(*triangleSet, 0, xDirection, margin, Dot3D(0, 0, 0));
    delete triangleSet;
    triangleSet = 0;
    TriangleBoundary3D<T> triangleBoundary(*defMesh);
    delete defMesh;
    defMesh = 0;
    triangleBoundary.getMesh().inflate(param.inflationParameter);

    // Voxelize the domain means: decide which lattice nodes are inside the solid
    // and which are outside.
    pcout << std::endl << "Voxelizing the simulation domain." << std::endl;
    int flowType = voxelFlag::outside;
    VoxelizedDomain3D<T> voxelizedDomain(triangleBoundary, flowType, param.bbox(), borderWidth,
            extendedEnvelopeWidth, blockSize);
    pcout << getMultiBlockInfo(voxelizedDomain.getVoxelMatrix()) << std::endl;
    {
        VtkImageOutput3D<T> vtkOut(outDir + "voxels_full_domain", param.dx);
        vtkOut.writeData<float>(*copyConvert<int,T>(voxelizedDomain.getVoxelMatrix(), param.bbox()), "voxel", 1.0);
    }

    /*
     * Generating fluid lattice and boundary conditions.
     */

    pcout << "Generating fluid lattice and boundary conditions." << std::endl;

    // Fluid lattice.
    Dynamics<T,DESCRIPTOR>* dynamics = 0;
    bool velIsJ = false;
    if (!util::isZero(param.cSmago)) {
        //dynamics = new SmagorinskyBGKdynamics<T,DESCRIPTOR>(param.omega, param.cSmago);
        //velIsJ = false;
        //dynamics = new SmagorinskyRegularizedDynamics<T,DESCRIPTOR>(param.omega, param.cSmago);
        //velIsJ = false;
        dynamics = new SmagorinskyIncBGKdynamics<T,DESCRIPTOR>(param.omega, param.cSmago);
        velIsJ = true;
    } else {
        //dynamics = new BGKdynamics<T,DESCRIPTOR>(param.omega);
        //velIsJ = false;
        //dynamics = new RegularizedBGKdynamics<T,DESCRIPTOR>(param.omega);
        //velIsJ = false;
        dynamics = new IncBGKdynamics<T,DESCRIPTOR>(param.omega);
        velIsJ = true;
    }

    MultiBlockLattice3D<T,DESCRIPTOR>* lattice = 0;
    if (param.useRhoBarJ) {
        lattice = generateMultiBlockLattice<T,DESCRIPTOR>(voxelizedDomain.getVoxelMatrix(),
                envelopeWidth, new NoDynamics<T,DESCRIPTOR>((T) 1)).release();
    } else {
        lattice = new MultiBlockLattice3D<T,DESCRIPTOR>((MultiBlock3D&) voxelizedDomain.getVoxelMatrix());
    }
    defineDynamics(*lattice, lattice->getBoundingBox(), dynamics->clone());
    delete dynamics;
    dynamics = 0;
    defineDynamics(*lattice, voxelizedDomain.getVoxelMatrix(), lattice->getBoundingBox(),
            new NoDynamics<T,DESCRIPTOR>((T) 1), voxelFlag::inside);
    lattice->toggleInternalStatistics(false);
    lattice->periodicity().toggleAll(false);

    MultiNTensorField3D<T>* rhoBarJfield = 0;
    std::vector<MultiBlock3D*> rhoBarJarg;
    if (param.useRhoBarJ) {
        // Here we assume that the object is far from the outlet. This is why we can compute
        // the rhoBarJ field to be used with the off-lattice BC, before the FluidPressureOutlet3D
        // is executed.
        plint numScalars = 4;
        rhoBarJfield = generateMultiNTensorField3D<T>(*lattice, extendedEnvelopeWidth, numScalars);
        rhoBarJarg.push_back(rhoBarJfield);
        integrateProcessingFunctional(new PackedRhoBarJfunctional3D<T,DESCRIPTOR>(),
                lattice->getBoundingBox(), *lattice, *rhoBarJfield, 0);
    }

    // The boundary condition algorithm for the object.
    BoundaryProfiles3D<T,Velocity> profiles;
    profiles.setWallProfile(new NoSlipProfile3D<T>());
    bool useAllDirections = true;
    GuoOffLatticeModel3D<T,DESCRIPTOR>* offLatticeModel = new GuoOffLatticeModel3D<T,DESCRIPTOR>(
            new TriangleFlowShape3D<T,Array<T,3> >(voxelizedDomain.getBoundary(), profiles),
            flowType, useAllDirections);
    offLatticeModel->selectSecondOrder(true);
    offLatticeModel->selectUseRegularizedModel(true);
    offLatticeModel->selectComputeStat(false);  // No force computation on object required.
    offLatticeModel->setVelIsJ(velIsJ);
    OffLatticeBoundaryCondition3D<T,DESCRIPTOR,Velocity>* boundaryCondition = 
        new OffLatticeBoundaryCondition3D<T,DESCRIPTOR,Velocity>(offLatticeModel, voxelizedDomain, *lattice);
    if (param.useRhoBarJ) {
        boundaryCondition->insert(rhoBarJarg);
    } else {
        boundaryCondition->insert();
    }

    // The boundary condition algorithm or the outer domain.
    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* outerBoundaryCondition =
        createLocalBoundaryCondition3D<T,DESCRIPTOR>();
    outerDomainBoundaries(*lattice, *outerBoundaryCondition);

    /*
     * Generating temperature lattice and boundary conditions.
     */

    pcout << "Generating temperature lattice and boundary conditions." << std::endl;

    // Temperature lattice.
    MultiBlockLattice3D<T,TEMPERATURE_DESCRIPTOR>* temperatureLattice = 
        new MultiBlockLattice3D<T,TEMPERATURE_DESCRIPTOR>((MultiBlock3D&) voxelizedDomain.getVoxelMatrix());
    defineDynamics(*temperatureLattice, temperatureLattice->getBoundingBox(),
            new AdvectionDiffusionRLBdynamics<T,TEMPERATURE_DESCRIPTOR>(param.omega_temperature));
    defineDynamics(*temperatureLattice, voxelizedDomain.getVoxelMatrix(), temperatureLattice->getBoundingBox(),
            new NoDynamics<T,TEMPERATURE_DESCRIPTOR>(param.T_wall), voxelFlag::inside);
    temperatureLattice->toggleInternalStatistics(false);
    temperatureLattice->periodicity().toggleAll(false);

    // The boundary condition algorithm for the object.
    BoundaryProfiles3D<T,Array<T,2> > temperatureProfiles;
    temperatureProfiles.setWallProfile(new ScalarDirichletProfile3D<T>(param.T_wall));
    GuoAdvDiffOffLatticeModel3D<T,TEMPERATURE_DESCRIPTOR>* advDiffOffLatticeModel =
        new GuoAdvDiffOffLatticeModel3D<T,TEMPERATURE_DESCRIPTOR>(
                new TriangleFlowShape3D<T,Array<T,2> >(
                    voxelizedDomain.getBoundary(), temperatureProfiles), flowType);
    advDiffOffLatticeModel->selectSecondOrder(true);
    OffLatticeBoundaryCondition3D<T,TEMPERATURE_DESCRIPTOR,Array<T,2> >* temperatureBoundaryCondition =
        new OffLatticeBoundaryCondition3D<T,TEMPERATURE_DESCRIPTOR,Array<T,2> >(
                advDiffOffLatticeModel, voxelizedDomain, *temperatureLattice);
    temperatureBoundaryCondition->insert();

    // The boundary condition algorithm or the outer domain.
    outerTemperatureBoundaries(*temperatureLattice);

    // Include coupling between fluid and temperature.
    integrateProcessingFunctional(
            new LatticeToPassiveAdvDiff3D<T,DESCRIPTOR,TEMPERATURE_DESCRIPTOR>((T) param.temperatureTimeFactor),
            lattice->getBoundingBox(), *lattice, *temperatureLattice, 1);

    // Minimal output.
    T Re = std::fabs(param.uInf) * param.characteristicLength / param.nu;
    pcout << "Reynolds number = " << Re << std::endl;
    pcout << "fluid tau       = " << (T) 1 / param.omega << std::endl;
    pcout << "temperature tau = " << (T) 1 / param.omega_temperature << std::endl;
    pcout << "dx              = " << param.dx << std::endl;
    pcout << "dt              = " << param.dt << std::endl;
    pcout << "dt/dx^2         = " << param.dt / util::sqr(param.dx) << std::endl;

    /*
     * Initialization.
     */

    std::vector<MultiBlock3D*> checkpointBlocks;
    checkpointBlocks.push_back(lattice);
    checkpointBlocks.push_back(temperatureLattice);

    bool continueSimulation = false;
    if (continueFileName != "") {
        continueSimulation = true;
    }

    plint iter = 0;
    if (continueSimulation) {
        pcout << "Reading state of the simulation from file: " << continueFileName << std::endl;
        loadState(checkpointBlocks, iter, param.saveDynamicContent, continueFileName);
        lattice->resetTime(iter);
        temperatureLattice->resetTime(iter);
    } else {
        T C = util::sinIncreasingFunction<T>(iter, param.startIter);
        initializeAtEquilibrium(*lattice, lattice->getBoundingBox(),
                SquarePoiseuilleDensityAndVelocity<T>(C, NMAX));
        lattice->initialize();
        initializeAtEquilibrium(*temperatureLattice, temperatureLattice->getBoundingBox(),
                param.T_inf, Array<T,3>((T) 0, (T) 0, (T) 0));
        temperatureLattice->initialize();
    }

    /*
     * Lattice-Boltzmann Cycles.
     */

    plb_ofstream energyFile;
    if (continueSimulation) {
        energyFile.open((outDir+"average_energy.dat").c_str(), std::ostream::app);
    } else {
        energyFile.open((outDir+"average_energy.dat").c_str());
    }

    plb_ofstream temperatureFile;
    if (continueSimulation) {
        temperatureFile.open((outDir+"average_temperature.dat").c_str(), std::ostream::app);
    } else {
        temperatureFile.open((outDir+"average_temperature.dat").c_str());
    }

    pcout << "Starting simulation." << std::endl;
    bool checkForErrors = true;
    bool stopProgram = false;
    for (plint i = iter; i < param.maxIter && !stopProgram; ++i) {
        bool output = (i == param.maxIter - 1) || stopProgram;

        if (i <= param.startIter) {
            T C = util::sinIncreasingFunction<T>(i, param.startIter);
            setBoundaryVelocity(*lattice, param.inlet(), SquarePoiseuilleVelocity<T>(C, NMAX));
        }

        if (i % param.statIter == 0 || output) {
            pcout << "At iteration: " << i << ", time: " << i * param.dt << std::endl;
            T avEnergy = boundaryCondition->computeAverageEnergy() * util::sqr(param.dx / param.dt);
            T avTemp = temperatureBoundaryCondition->computeAverageDensity();
            pcout << "Average energy = " << avEnergy << std::endl;
            pcout << "Average temperature = " << avTemp << std::endl;
            energyFile << i * param.dt << " " << avEnergy << std::endl;
            temperatureFile << i * param.dt << " " << avTemp << std::endl;
            if (i != 0) {
                pcout << "Time per lattice Boltzmann cycle: " << global::timer("cycle").getTime() / (T) i << std::endl;
            }
            pcout << std::endl;
        }
        if (i % param.vtkIter == 0 || output) {
            pcout << "Writing VTK at iteration: " << i << std::endl;
            writeVTK(*boundaryCondition, *temperatureBoundaryCondition, i, outDir);
            pcout << std::endl;
        }
        if ((param.cpIter > 0 && i % param.cpIter == 0 && i != iter) || output) {
            pcout << "Saving the state of the simulation at iteration: " << i << std::endl;
            saveState(checkpointBlocks, i, param.saveDynamicContent, param.continueFile,
                    param.checkpointFile, PADD);
            pcout << std::endl;
        }
        if (i % param.abIter == 0) {
            stopProgram = abortExecution(param.abortFile, checkpointBlocks, i, param.saveDynamicContent,
                    param.continueFile, param.checkpointFile, PADD);

            if (stopProgram) {
                pcout << "Aborting execution at iteration: " << i << std::endl;
                pcout << std::endl;
            }
        }

        global::timer("cycle").start();
        if (i >= param.startIter && i % param.temperatureTimeFactor == 0) {
            temperatureLattice->collideAndStream();
        }
        lattice->collideAndStream();
        global::timer("cycle").stop();

        if (checkForErrors) {
            abortIfErrorsOccurred();
            checkForErrors = false;
        }
    }

    energyFile.close();
    temperatureFile.close();

    delete temperatureBoundaryCondition;
    delete temperatureLattice;
    delete outerBoundaryCondition;
    delete boundaryCondition;
    delete rhoBarJfield;
    delete lattice;
}


int main(int argc, char* argv[])
{
    plbInit(&argc, &argv);

    // For better debugging.
    enableCoreDumps();
    unbufferOutputStdStreams();

    // The try-catch blocks catch exceptions in case an error occurs,
    // and terminate the program properly with a nice error message.

    // 1. Read command-line parameter: the input file name.
    std::string xmlFileName;
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

    // Some clusters have fast parallel input/output facilities
    // which Palabos can exploit.
    global::IOpolicy().activateParallelIO(param.useParallelIO);

    std::string continueFileName = "";
    try {
        global::argv(2).read(continueFileName);
    }
    catch (PlbIOException& exception) { }

    // 3. Execute the main program.
    try {
        run(continueFileName);
    }
    catch (PlbIOException& exception) {
        pcout << exception.what() << std::endl;
        return -1;
    }

    return 0;
}
