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
  * This code solves the steady flow inside an aneurysm. It introduces several
  * new concepts like Guo off lattice boundary conditions, reading of
  * surface geometry STL files, smooth grid refinement and voxelization.
  * Make sure to unpack the file aneurysm.stl.tgz before running the
  * simulation.
  **/

#include "palabos3D.h"
#include "palabos3D.hh"

using namespace plb;
using namespace std;

typedef double T;
typedef Array<T,3> Velocity;
#define DESCRIPTOR descriptors::D3Q19Descriptor

plint extraLayer      = 0;  // Make the bounding box larger; for visualization purposes
                            //   only. For the simulation, it is OK to have extraLayer=0.
const plint blockSize = 20; // Zero means: no sparse representation.
const plint envelopeWidth = 1;  // For standard BGK dynamics.
const plint extendedEnvelopeWidth = 2;  // Because the Guo off lattice boundary condition
                                        //   needs 2-cell neighbor access.

bool performOutput = false;
bool doImages = false;
bool useAllDirections = false;
bool useRegularizedWall = false;
bool useIncompressible = false;
bool poiseuilleInlet = false;
bool convectiveScaling = false;

T kinematicViscosity       = 0.;
T averageInletVelocity     = 0.;
plint referenceResolution  = 0.;
T nuLB                     = 0.;
T fluidDensity             = 0.;
T volume                   = 0.;
T userDefinedInletDiameter = 0.;

plint referenceDirection = 0;
plint openingSortDirection = 0;

T simTime = 0;
plint startLevel = 0;
plint maxLevel   = 0;
T epsilon = 0.;

TriangleSet<T>* triangleSet = 0;
T currentTime = 0;

// Structure which defines an ``opening''. The surface geometry of the aneurysm,
//   as given by the user in the form of an STL file, contains holes, which in 
//   the specific simulation represent inlets and outlets.
template<typename T>
struct Opening {
    bool inlet;
    Array<T,3> center;
    T innerRadius;
};

std::vector<Opening<T> > openings;

void iniLattice( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                 VoxelizedDomain3D<T>& voxelizedDomain )
{
    // Switch all remaining outer cells to no-dynamics, except the outer
    //   boundary layer, and keep the rest as BGKdynamics.
    defineDynamics(lattice, voxelizedDomain.getVoxelMatrix(), lattice.getBoundingBox(),
                   new NoDynamics<T,DESCRIPTOR>, voxelFlag::outside);
    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), (T) 1., Array<T,3>((T) 0.,(T) 0.,(T) 0.));
    lattice.initialize();
}

// This function assigns proper boundary conditions to the openings of the surface geometry
//   of the aneurysm. Which opening is inlet and which is outlet is defined by the user in
//   the input XML file. For the inlet, there is a choice between a Poiseuille velocity
//   profile and a simple plug velocity profile. At the outlets a Neumann boundary
//   condition with constant pressure is prescribed.
void setOpenings (
    std::vector<BoundaryProfile3D<T,Velocity>*>& inletOutlets,
    TriangleBoundary3D<T>& boundary, T uLB, T dx, T dt )
{
    for (pluint i=0; i<openings.size(); ++i) {
        Opening<T>& opening = openings[i];
        opening.center = computeBaryCenter (
                boundary.getMesh(),
                boundary.getInletOutlet(openingSortDirection)[i] );
        opening.innerRadius = computeInnerRadius (
                boundary.getMesh(),
                boundary.getInletOutlet(openingSortDirection)[i] );

        if (opening.inlet) {
            if (poiseuilleInlet) {
                inletOutlets.push_back (
                        new PoiseuilleProfile3D<T>(uLB) );
            }
            else {
                inletOutlets.push_back (
                        new VelocityPlugProfile3D<T>(uLB) );
            }
        }
        else {
            inletOutlets.push_back (
                    new DensityNeumannBoundaryProfile3D<T> );
        }
    }
}

// This function outputs velocity, vorticity and pressure data, at selected
//   points of the computational domain, given their coordinates in physical units.
std::vector<T> pointMeasures (
        MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
        Array<T,3> location, T dx, T dt )
{
    std::vector<Array<T,3> > physicalPositions, positions;
    physicalPositions.push_back(Array<T,3>(0.022046, 0.015072, 0.044152));
    physicalPositions.push_back(Array<T,3>(0.027132, 0.049947, 0.095012));
    physicalPositions.push_back(Array<T,3>(0.034398, 0.056487, 0.057957));
    physicalPositions.push_back(Array<T,3>(0.031492, 0.025971, 0.084113));
    physicalPositions.push_back(Array<T,3>(0.025679, 0.025971, 0.091379));
    physicalPositions.push_back(Array<T,3>(0.018413, 0.011439, 0.076848));
    positions.resize(physicalPositions.size());

    for (pluint i=0; i<physicalPositions.size(); ++i) {
        positions[i] = (physicalPositions[i]-location)/dx;
    }

    std::vector<Array<T,3> > velocities = velocitySingleProbes(lattice, positions);
    std::vector<Array<T,3> > vorticities = vorticitySingleProbes(lattice, positions);
    std::vector<T> densities = densitySingleProbes(lattice, positions);

    std::vector<T> data;
    for (pluint i=0; i<physicalPositions.size(); ++i) {
        Array<T,3> pos = physicalPositions[i];
        Array<T,3> vel = velocities[i]*dx/dt;
        Array<T,3> vort = vorticities[i]/dt;
        T pressure = DESCRIPTOR<T>::cs2*(densities[i]-1.)*dx*dx/(dt*dt)*fluidDensity;
        if (performOutput) {
            pcout << "Pos ("
                  << pos[0] << "," << pos[1] << "," << pos[2]
                  << "); Velocity ("
                  << vel[0] << "," << vel[1] << "," << vel[2]
                  << "); Vorticity ("
                  << vort[0] << "," << vort[1] << "," << vort[2]
                  << "); Pressure " << pressure << std::endl;
        }
        data.push_back(norm(vel));
        data.push_back(norm(vort));
        data.push_back(pressure);
    }
    return data;
}

void writeImages (
         MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
         MultiScalarField3D<int>& flagMatrix,
         Box3D const& imageDomain, Box3D const& vtkDomain, std::string fname, Array<T,3> location, T dx, T dt )
{
    ImageWriter<T> imageWriter("leeloo");
    imageWriter.writeScaledPpm(fname, *computeVelocityNorm(lattice, imageDomain));

    VtkImageOutput3D<T> vtkOut(fname, dx, location);
    vtkOut.writeData<float>(*computeDensity(lattice, vtkDomain), "p", util::sqr(dx/dt)*fluidDensity);
    vtkOut.writeData<float>(*computeVelocityNorm(lattice, vtkDomain), "u", dx/dt);
    vtkOut.writeData<float>(*copyConvert<int,T>(*extractSubDomain(flagMatrix, vtkDomain)), "flag", 1.);
}

// This function produces images at predefined yz, xz and xy planes. The coordinates of the planes are given
//   in physical coordinates, and the output variables are velocity, vorticity and pressure.
void writeImages (
         MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
         MultiScalarField3D<int>& flag,
         plint level, Array<T,3> location, T dx, T dt )
{
    plint nx = lattice.getNx();
    plint ny = lattice.getNy();
    plint nz = lattice.getNz();
    Array<T,3> yz_plane(0.016960, 0.032604, 0.057772);
    Array<T,3> xz_plane(0.026725, 0.017978, 0.057772);
    Array<T,3> xy_plane(0.026725, 0.032604, 0.084113);

    Array<T,3> lyz_plane((yz_plane-location)/dx);
    Array<T,3> lxz_plane((xz_plane-location)/dx);
    Array<T,3> lxy_plane((xy_plane-location)/dx);

    Box3D yz_imageDomain (
            util::roundToInt(lyz_plane[0]), util::roundToInt(lyz_plane[0]),
            0, ny-1, 0, nz-1 );
    Box3D xz_imageDomain (
            0, nx-1,
            util::roundToInt(lxz_plane[1]), util::roundToInt(lxz_plane[1]),
            0, nz-1 );
    Box3D xy_imageDomain (
            0, nx-1, 0, ny-1,
            util::roundToInt(lxy_plane[2]), util::roundToInt(lxy_plane[2]) );

    Box3D yz_vtkDomain (
            util::roundToInt(lyz_plane[0])-3, util::roundToInt(lyz_plane[0])+3,
            0, ny-1, 0, nz-1 );
    Box3D xz_vtkDomain (
            0, nx-1,
            util::roundToInt(lxz_plane[1])-3, util::roundToInt(lxz_plane[1])+3,
            0, nz-1 );
    Box3D xy_vtkDomain (
            0, nx-1, 0, ny-1,
            util::roundToInt(lxy_plane[2])-3, util::roundToInt(lxy_plane[2])+3 );

    writeImages(lattice, flag, xy_imageDomain, xy_vtkDomain, "xy_"+util::val2str(level), location, dx, dt);
    writeImages(lattice, flag, xz_imageDomain, xz_vtkDomain, "xz_"+util::val2str(level), location, dx, dt);
    writeImages(lattice, flag, yz_imageDomain, yz_vtkDomain, "yz_"+util::val2str(level), location, dx, dt);
}


// This is the function that prepares and performs the actual simulation.
std::auto_ptr<MultiBlockLattice3D<T,DESCRIPTOR> > run (
        plint level, MultiBlockLattice3D<T,DESCRIPTOR>* iniVal=0 )
{
    plint margin = 3; // Extra margin of allocated cells around the obstacle. 
    plint borderWidth = 1; // Because the Guo boundary condition acts in a one-cell layer.
                           // Requirement: margin>=borderWidth.

    // The resolution is doubled at each coordinate direction with the increase of the
    //   resolution level by one. The parameter ``referenceResolution'' is by definition
    //   the resolution at grid refinement level 0.
    plint resolution = referenceResolution * util::twoToThePower(level);

    // The next few lines of code are typical. They transform the surface geometry of the
    //   aneurysm given by the user to more efficient data structures that are internally
    //   used by palabos. The TriangleBoundary3D structure will be later used to assign
    //   proper boundary conditions.
    DEFscaledMesh<T>* defMesh =
        new DEFscaledMesh<T>(*triangleSet, resolution, referenceDirection, margin, extraLayer);
    TriangleBoundary3D<T> boundary(*defMesh);
    delete defMesh;
    boundary.getMesh().inflate();

    // When convective scaling is used (relationship of dt with respect to dx as the grid is
    //   refined) the value of the kinematic viscosity must be also properly adjusted.
    T nuLB_ = nuLB;
    if (convectiveScaling) {
        nuLB_ = nuLB * util::twoToThePower(level);
    }
    T dx = boundary.getDx();
    T dt = nuLB_ / kinematicViscosity *dx*dx;
    T uAveLB = averageInletVelocity *dt/dx;
    T omega = 1./(3.*nuLB_+0.5);
    Array<T,3> location(boundary.getPhysicalLocation());


    pcout << "uLB=" << uAveLB << std::endl;
    pcout << "nuLB=" << nuLB_ << std::endl;
    pcout << "tau=" << 1./omega << std::endl;
    if (performOutput) {
        pcout << "dx=" << dx << std::endl;
        pcout << "dt=" << dt << std::endl;
    }

    // The aneurysm simulation is an interior (as opposed to exterior) flow problem. For
    //   this reason, the lattice nodes that lay inside the computational domain must
    //   be identified and distinguished from the ones that lay outside of it. This is
    //   handled by the following voxelization process.
    const int flowType = voxelFlag::inside;
    VoxelizedDomain3D<T> voxelizedDomain (
            boundary, flowType, extraLayer, borderWidth, extendedEnvelopeWidth, blockSize );
    if (performOutput) {
        pcout << getMultiBlockInfo(voxelizedDomain.getVoxelMatrix()) << std::endl;
    }

    MultiScalarField3D<int> flagMatrix((MultiBlock3D&)voxelizedDomain.getVoxelMatrix());
    setToConstant(flagMatrix, voxelizedDomain.getVoxelMatrix(),
                  voxelFlag::inside, flagMatrix.getBoundingBox(), 1);
    setToConstant(flagMatrix, voxelizedDomain.getVoxelMatrix(),
                  voxelFlag::innerBorder, flagMatrix.getBoundingBox(), 1);
    pcout << "Number of fluid cells: " << computeSum(flagMatrix) << std::endl;

    Dynamics<T,DESCRIPTOR>* dynamics = 0;
    if (useIncompressible) {
        dynamics = new IncBGKdynamics<T,DESCRIPTOR>(omega); // In this model velocity equals momentum.
    }
    else {
        dynamics = new BGKdynamics<T,DESCRIPTOR>(omega); // In this model velocity equals momentum
                                                         //   divided by density.
    }
    std::auto_ptr<MultiBlockLattice3D<T,DESCRIPTOR> > lattice 
        = generateMultiBlockLattice<T,DESCRIPTOR> (
                voxelizedDomain.getVoxelMatrix(), envelopeWidth, dynamics );
    lattice->toggleInternalStatistics(false);

    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
        = createLocalBoundaryCondition3D<T,DESCRIPTOR>();

    Array<T,3> inletRealPos(0.024, 0.018, 0.002);
    Array<T,3> outlet1RealPos(0.038, 0.06, 0.055);
    Array<T,3> outlet2RealPos(0.007,0.06,0.105);
    T diameterReal = 0.01;

    Array<T,3> inletPos((inletRealPos-location)/dx);
    Array<T,3> outlet1Pos((outlet1RealPos-location)/dx);
    Array<T,3> outlet2Pos((outlet2RealPos-location)/dx);
    plint diameter = util::roundToInt(diameterReal/dx);

    Box3D inletDomain(util::roundToInt(inletPos[0]-diameter), util::roundToInt(inletPos[0]+diameter),
                      util::roundToInt(inletPos[1]-diameter), util::roundToInt(inletPos[1]+diameter),
                      util::roundToInt(inletPos[2]), util::roundToInt(inletPos[2]));
    Box3D behindInlet(inletDomain.x0, inletDomain.x1,
                      inletDomain.y0, inletDomain.y1,
                      inletDomain.z0-diameter, inletDomain.z0-1);

    Box3D outlet1Domain(util::roundToInt(outlet1Pos[0]-diameter), util::roundToInt(outlet1Pos[0]+diameter),
                        util::roundToInt(outlet1Pos[1]-diameter), util::roundToInt(outlet1Pos[1]+diameter),
                        util::roundToInt(outlet1Pos[2]), util::roundToInt(outlet1Pos[2]));
    Box3D behindOutlet1(outlet1Domain.x0, outlet1Domain.x1,
                        outlet1Domain.y0, outlet1Domain.y1,
                        outlet1Domain.z0-diameter, outlet1Domain.z0-1);

    Box3D outlet2Domain(util::roundToInt(outlet2Pos[0]), util::roundToInt(outlet2Pos[0]),
                        util::roundToInt(outlet2Pos[1]-diameter), util::roundToInt(outlet2Pos[1]+diameter),
                        util::roundToInt(outlet2Pos[2]-diameter), util::roundToInt(outlet2Pos[2]+diameter));
    Box3D behindOutlet2(outlet2Domain.x0-diameter, outlet2Domain.x1-1,
                        outlet2Domain.y0, outlet2Domain.y1,
                        outlet2Domain.z0, outlet2Domain.z1);


    boundaryCondition->addVelocityBoundary2N(inletDomain, *lattice);
    setBoundaryVelocity(*lattice, inletDomain, Array<T,3>((T)0.,(T)0.,uAveLB));
    boundaryCondition->addPressureBoundary2N(outlet1Domain, *lattice);
    setBoundaryDensity(*lattice, outlet1Domain, (T)1.);
    boundaryCondition->addPressureBoundary0N(outlet2Domain, *lattice);
    setBoundaryDensity(*lattice, outlet2Domain, (T)1.);

    defineDynamics(*lattice, flagMatrix, lattice->getBoundingBox(), new BounceBack<T,DESCRIPTOR>(1.), 0);
    defineDynamics(*lattice, behindInlet, new BounceBack<T,DESCRIPTOR>(1.));
    defineDynamics(*lattice, behindOutlet1, new BounceBack<T,DESCRIPTOR>(1.));
    defineDynamics(*lattice, behindOutlet2, new BounceBack<T,DESCRIPTOR>(1.));

    iniLattice(*lattice, voxelizedDomain);
    if(iniVal) {
        Box3D toDomain(lattice->getBoundingBox());
        Box3D fromDomain(toDomain.shift(margin,margin,margin)); // During rescaling, the margin doubled in size,
                                                                //   an effect which is cancelled here through a shift.
        copy(*iniVal, fromDomain, *lattice, toDomain, modif::staticVariables);
    }

    // The ValueTracer is needed to check when a chosen quantity (in our case the average energy)
    //   has converged, so to conclude that steady state has been reached for the specific grid
    //   refinement level and stop the simulation.
    plint convergenceIter=20;
    util::ValueTracer<T> velocityTracer(0.05*convergenceIter, resolution, epsilon);
    global::timer("iteration").restart();
    plint i = util::roundToInt(currentTime/dt);
    lattice->resetTime(i);

    pcout << "Saving a " << lattice->getNx() << " by " << lattice->getNy()
          << " by " << lattice->getNz() << " lattice." << std::endl;
    global::timer("io").start();
    parallelIO::save(*lattice, "checkpoint", false);
    pcout << "Total time for i/o: " << global::timer("io").getTime() << std::endl;

    // Collision and streaming iterations.
    pcout << "Starting 100 iterations" << std::endl;
    global::timer("global").start();
    for (plint i=0; i<100; ++i) {

        lattice->collideAndStream();

        ++i;
        currentTime = i*dt;
    }
    pcout << "End of 100 iterations" << std::endl;
    pcout << "Total time of execution: " << global::timer("global").getTime() << std::endl;

    return lattice;
}

// Read the user input XML file provided at the command-line.
void readParameters(XMLreader const& document)
{
    std::string meshFileName;
    std::vector<std::string> openingType;
    document["geometry"]["mesh"].read(meshFileName);
    document["geometry"]["inletDiameter"].read(userDefinedInletDiameter);
    document["geometry"]["averageInletVelocity"].read(averageInletVelocity);
    document["geometry"]["openings"]["sortDirection"].read(openingSortDirection);
    document["geometry"]["openings"]["type"].read(openingType);

    document["fluid"]["kinematicViscosity"].read(kinematicViscosity);
    document["fluid"]["density"].read(fluidDensity);
    document["fluid"]["volume"].read(volume);

    document["numerics"]["referenceDirection"].read(referenceDirection);
    document["numerics"]["referenceResolution"].read(referenceResolution);
    document["numerics"]["nuLB"].read(nuLB);

    document["simulation"]["simTime"].read(simTime);
    document["simulation"]["maxLevel"].read(maxLevel);
    document["simulation"]["epsilon"].read(epsilon);

    document["simulation"]["performOutput"].read(performOutput);
    document["simulation"]["doImages"].read(doImages);
    document["simulation"]["useAllDirections"].read(useAllDirections);
    document["simulation"]["useRegularizedWall"].read(useRegularizedWall);
    document["simulation"]["useIncompressible"].read(useIncompressible);
    document["simulation"]["poiseuilleInlet"].read(poiseuilleInlet);
    document["simulation"]["convectiveScaling"].read(convectiveScaling);

    // At this part, the surface geometry of the aneurysm (as given by the user in
    //   the form of an ASCII or binary STL file) is read into a data structure
    //   comprised by a set of triangles. The DBL constant means that double
    //   precision accuracy will be used (generally the recommended choice).
    triangleSet = new TriangleSet<T>(meshFileName, DBL);
    pcout << "Reynolds number, based on provided inlet diameter: "
          << averageInletVelocity*userDefinedInletDiameter/kinematicViscosity
          << std::endl;
    plbIOError(openingSortDirection<0 || openingSortDirection>2,
               "Sort-direction of opening must be 0 (x), 1 (y), or 2 (z).");
    // The surface geometry, as provided by the STL file, must contain openings,
    //   namely inlets and outlets. On these openings, appropriate boundary conditions
    //   will be imposed by palabos. Which opening is inlet and which is outlet, is
    //   identified by the user in the input XML file.
    openings.resize(openingType.size());
    for (pluint i=0; i<openingType.size(); ++i) {
        std::string next_opening = util::tolower(openingType[i]);
        if (next_opening=="inlet") {
            openings[i].inlet = true;
        }
        else if (next_opening=="outlet") {
            openings[i].inlet = false;
        }
        else {
            plbIOError("Unknown opening type.");
        }
    }
}


int main(int argc, char* argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./");
    global::IOpolicy().activateParallelIO(true);

    string paramXmlFileName;
    try {
        global::argv(1).read(paramXmlFileName);
    }
    catch (PlbIOException& exception) {
        pcout << "Wrong parameters; the syntax is: " 
              << (std::string)global::argv(0) << " parameter-input-file.xml" << std::endl;
        return -1;
    }

    // Read the parameter XML input file. (Lots of comments are included there too).
    try {
        XMLreader document(paramXmlFileName);
        readParameters(paramXmlFileName);
    }
    catch (PlbIOException& exception) {
        pcout << "Error in input file " << paramXmlFileName
              << ": " << exception.what() << std::endl;
        return -1;
    }

    plint iniLevel=0;
    std::auto_ptr<MultiBlockLattice3D<T,DESCRIPTOR> > iniConditionLattice(0);
    // This code incorporates the concept of smooth grid refinement until convergence is
    //   achieved. The word ``smooth'' indicates that as the refinement level increases
    //   by one, the whole grid doubles in each direction. When the grid is refined, both
    //   dx and dt have to change. Whether dt is changed as dx^2 (diffusive behavior)
    //   or as dx (convective behavior), is controlled by the input variable
    //   ``convectiveScaling'' (the recommended choice is not to use convective scaling).
    try {
        for (plint level=iniLevel; level<=maxLevel; ++level) {
            pcout << std::endl << "Running new simulation at level " << level << std::endl;
            std::auto_ptr<MultiBlockLattice3D<T,DESCRIPTOR> > convergedLattice (
                    run(level, iniConditionLattice.get()) );
            if (level != maxLevel) {
                plint dxScale = -1;
                plint dtScale = -2;
                if (convectiveScaling) {
                    dtScale = -1;
                }
                // The converged simulation of the previous grid level is used as the initial condition
                //   for the simulation at the next grid level (after appropriate interpolation has
                //   taken place).
                iniConditionLattice = std::auto_ptr<MultiBlockLattice3D<T,DESCRIPTOR> > (
                        refine(*convergedLattice, dxScale, dtScale, new BGKdynamics<T,DESCRIPTOR>(1.)) );
            }
        }
    }
    catch(PlbException& exception) {
        pcout << exception.what() << std::endl;
        return -1;
    }
}

