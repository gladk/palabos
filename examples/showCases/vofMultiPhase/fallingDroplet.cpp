/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2015 FlowKit Sarl
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

using namespace plb;
using namespace std;

//#define MRT

#ifdef MRT
    #define DESCRIPTOR descriptors::ForcedMRTD3Q19Descriptor
#else
    #define DESCRIPTOR descriptors::ForcedD3Q19Descriptor
#endif

#define PADDING 8
#define CBUFSIZ 256

typedef double T;

// Smagorinsky constant for LES model.
const T cSmago = 0.14;

// Physical dimensions of the system (in meters).
const T lx = 0.02;
const T ly = 0.02;
const T lz = 0.06;

const T radius = 0.002;

const T rhoEmpty = 1.0;
const T rho = 1000.0; // Water density.
    
plint writeImagesIter;
plint statIter;

plint maxIter;
plint N;
plint nx, ny, nz, radiusLB;
T delta_t, delta_x;
Array<T,3> externalForce;
Array<T,3> externalForce2;
T nuPhys, nuLB, tau, omega, surfaceTensionPhys, surfaceTensionLB, contactAngle;
T nuPhys2, nuLB2, tau2, omega2;
T densityRatio, viscosityRatio;

TwoPhaseModel model;

void setupParameters()
{
    delta_x = radius / (N - 1.0);
    nx = util::roundToInt(lx / delta_x);
    ny = util::roundToInt(ly / delta_x);
    nz = util::roundToInt(lz / delta_x) + 1;

    radiusLB = util::roundToInt(radius / delta_x);

    // Gravity in lattice units.
    T gLB = 9.8 * delta_t * delta_t / delta_x;
    externalForce = Array<T,3>(0., 0., -gLB);
    tau            = (nuPhys*DESCRIPTOR<T>::invCs2*delta_t)/(delta_x*delta_x) + 0.5;
    omega          = 1./tau;    
    nuLB           = (tau-0.5)*DESCRIPTOR<T>::cs2; // Viscosity in lattice units.

    nuPhys2 = nuPhys * viscosityRatio;
    tau2           = (nuPhys2*DESCRIPTOR<T>::invCs2*delta_t)/(delta_x*delta_x) + 0.5;
    omega2         = 1./tau2;    
    nuLB2          = (tau2-0.5)*DESCRIPTOR<T>::cs2; // Viscosity in lattice units.
    externalForce2 = densityRatio * externalForce;
    
    surfaceTensionLB = (rhoEmpty / rho) * (delta_t*delta_t) / (delta_x*delta_x*delta_x) * surfaceTensionPhys;

    writeImagesIter = 40;
    statIter = 80;
}

bool insideFluid(T x, T y, T z)
{
    Array<T,3> pos(x, y, z);
    Array<T,3> center(nx/2.0, ny/2.0, 3.*nz/4.0);
    T r = norm(pos-center);
    if (r <= radiusLB) {
        return true;
    }
    return false;
}

void writeTwoPhaseVTK(TwoPhaseFields3D<T,DESCRIPTOR> *fields, plint iT)
{
    std::auto_ptr<MultiScalarField3D<T> > smoothVF(lbmSmoothen<T,DESCRIPTOR>(fields->volumeFraction,
                fields->volumeFraction.getBoundingBox().enlarge(-1)));

    std::vector<T> isoLevels;
    isoLevels.push_back(0.5);

    typedef TriangleSet<T>::Triangle Triangle;
    std::vector<Triangle> triangles;
    isoSurfaceMarchingCube(triangles, *smoothVF, isoLevels, smoothVF->getBoundingBox().enlarge(-4));
    {
        TriangleSet<T> triangleSet(triangles);
        triangleSet.scale(delta_x);
        triangleSet.writeBinarySTL(createFileName("out_interface_", iT, PADDING)+".stl");
    }

    plint centx = util::roundToInt(0.5 * nx);
    plint centy = util::roundToInt(0.5 * ny);
    plint centz = util::roundToInt(0.5 * nz);

    Box3D box_x(centx - 1, centx + 1, 0, ny - 1, 0, nz - 1);
    Box3D box_y(0, nx - 1, centy - 1, centy + 1, 0, nz - 1);
    Box3D box_z(0, nx - 1, 0, ny - 1, centz - 1, centz + 1);

    char fname_x[CBUFSIZ];
    char fname_y[CBUFSIZ];
    char fname_z[CBUFSIZ];
    sprintf(fname_x, "out_phases_x_");
    sprintf(fname_y, "out_phases_y_");
    sprintf(fname_z, "out_phases_z_");

    bool computeFluid1 = true;
    bool computeFluid2 = true;

    {
        VtkImageOutput3D<T> vtkOut_x(createFileName(fname_x, iT, PADDING), delta_x);
        std::auto_ptr<MultiTensorField3D<T,3> > vx = fields->computeVelocity(box_x, computeFluid1, computeFluid2);
        std::auto_ptr<MultiScalarField3D<T> > px = fields->computePressure(box_x, computeFluid1, computeFluid2);
        vtkOut_x.writeData<3,float>(*vx, "v", delta_x / delta_t);
        vtkOut_x.writeData<float>(*px, "p", rho * (delta_x * delta_x) / (delta_t * delta_t));
        vtkOut_x.writeData<float>(*extractSubDomain(fields->volumeFraction, box_x), "vf", 1.0);
    }

    {
        VtkImageOutput3D<T> vtkOut_y(createFileName(fname_y, iT, PADDING), delta_x);
        std::auto_ptr<MultiTensorField3D<T,3> > vy = fields->computeVelocity(box_y, computeFluid1, computeFluid2);
        std::auto_ptr<MultiScalarField3D<T> > py = fields->computePressure(box_y, computeFluid1, computeFluid2);
        vtkOut_y.writeData<3,float>(*vy, "v", delta_x / delta_t);
        vtkOut_y.writeData<float>(*py, "p", rho * (delta_x * delta_x) / (delta_t * delta_t));
        vtkOut_y.writeData<float>(*extractSubDomain(fields->volumeFraction, box_y), "vf", 1.0);
    }

    {
        VtkImageOutput3D<T> vtkOut_z(createFileName(fname_z, iT, PADDING), delta_x);
        std::auto_ptr<MultiTensorField3D<T,3> > vz = fields->computeVelocity(box_z, computeFluid1, computeFluid2);
        std::auto_ptr<MultiScalarField3D<T> > pz = fields->computePressure(box_z, computeFluid1, computeFluid2);
        vtkOut_z.writeData<3,float>(*vz, "v", delta_x / delta_t);
        vtkOut_z.writeData<float>(*pz, "p", rho * (delta_x * delta_x) / (delta_t * delta_t));
        vtkOut_z.writeData<float>(*extractSubDomain(fields->volumeFraction, box_z), "vf", 1.0);
    }
}


int main(int argc, char **argv)
{
    plbInit(&argc, &argv);

    std::string modelName;
    try {
        global::argv(1).read(nuPhys);
        global::argv(2).read(surfaceTensionPhys);
        global::argv(3).read(contactAngle);
        global::argv(4).read(viscosityRatio);
        global::argv(5).read(densityRatio);
        global::argv(6).read(N);
        global::argv(7).read(delta_t);
        global::argv(8).read(modelName);
        global::argv(9).read(maxIter);
    }
    catch(PlbIOException& except) {
        pcout << except.what() << std::endl;
        pcout << "The parameters for this program are :\n";
        pcout << "1. Kinematic viscosity in physical units (m^2/s).\n";
        pcout << "2. Surface tension in physical units (N/m).\n";
        pcout << "3. Contact angle (in degrees).\n";
        pcout << "4. Kinematic viscosity Ratio.\n";
        pcout << "5. Density Ratio.\n";
        pcout << "6. Number of lattice nodes for the sphere radius.\n";
        pcout << "7. Time step in physical units (s).\n"; 
        pcout << "8. Model (kinetic / dynamic / constRho).\n";
        pcout << "9. Maximum number of iterations.\n";
        pcout << "Reasonable parameters on a desktop computer are: "
              << (std::string)global::argv(0)
              << " 1.e-3 0.0728 90.0 2.0 0.01 9 1.e-5 dynamic 500000\n";
        exit (EXIT_FAILURE);
    }

    model = stringToTwoPhaseModel(modelName);
    
    setupParameters();
    
    pcout << "delta_t = " << delta_t << endl;
    pcout << "delta_x = " << delta_x << endl;
    pcout << "externalForce = " << externalForce[2] << endl;
    pcout << "externalForce2 = " << externalForce2[2] << endl;
    pcout << "relaxation time = " << tau << endl;
    pcout << "relaxation time fluid 2 = " << tau2 << endl;
    pcout << "kinematic viscosity physical units = " << nuPhys << endl;
    pcout << "kinematic viscosity fluid 2, physical units = " << nuPhys2 << endl;
    pcout << "kinematic viscosity lattice units = " << nuLB << endl;
    pcout << "kinematic viscosity 2, lattice units = " << nuLB2 << endl;
    pcout << "surface tension, lattice units = " << surfaceTensionLB << endl;
    pcout << "Size of the domain is " << nx << " x " << ny << " x " << nz << endl;
    
    SparseBlockStructure3D blockStructure(createRegularDistribution3D(nx, ny, nz));

#ifdef MRT
    Dynamics<T,DESCRIPTOR>* dynamics
        = new SmagorinskyDynamics<T,DESCRIPTOR>(new VariableOmegaMRTdynamics<T,DESCRIPTOR>(omega), omega, cSmago);

    Dynamics<T,DESCRIPTOR>* dynamics2
        = new SmagorinskyDynamics<T,DESCRIPTOR>(new VariableOmegaMRTdynamics<T,DESCRIPTOR>(omega2), omega2, cSmago);
#else
    Dynamics<T,DESCRIPTOR>* dynamics
        = new SmagorinskyBGKdynamics<T,DESCRIPTOR>(omega, cSmago);

    Dynamics<T,DESCRIPTOR>* dynamics2
        = new SmagorinskyBGKdynamics<T,DESCRIPTOR>(omega2, cSmago);
#endif

    // If surfaceTensionLB is 0, then the surface tension algorithm is deactivated.
    // If contactAngle is less than 0, then the contact angle algorithm is deactivated.
    TwoPhaseFields3D<T,DESCRIPTOR> fields(blockStructure, dynamics->clone(), dynamics2->clone(), rhoEmpty, densityRatio,
                                          surfaceTensionLB, contactAngle, externalForce, externalForce2, model);
    pcout << "Setting initial fluid flags." << std::endl;

    // Initialization

    analyticalIniVolumeFraction(fields.volumeFraction, fields.flag, insideFluid, 32);

    Box3D bottom  (0,    nx-1, 0,    ny-1, 0,    0);
    Box3D top     (0,    nx-1, 0,    ny-1, nz-1, nz-1);

    setToConstant(fields.flag, bottom,   (int) twoPhaseFlag::wall);
    setToConstant(fields.flag, top,      (int) twoPhaseFlag::wall);

    fields.periodicityToggle(0, true);
    fields.periodicityToggle(1, true);
    fields.periodicityToggle(2, false);
    
    fields.partiallyDefaultInitialize();

    for (plint iT = 0; iT <= maxIter; ++iT) {
        if (iT % writeImagesIter == 0) {
            writeTwoPhaseVTK(&fields, iT);
        }

        if (iT % statIter == 0) {
            pcout << "At iteration " << iT << ", t = " << iT * delta_t << std::endl;
            T avE = computeAverageEnergy(fields.lattice);
            avE *= delta_x * delta_x / (delta_t * delta_t);
            pcout << "Average kinetic energy of fluid 1: " << avE << std::endl;
            T avE2 = computeAverageEnergy(*fields.lattice2);
            avE2 *= delta_x * delta_x / (delta_t * delta_t);
            pcout << "Average kinetic energy of fluid 2: " << avE2 << std::endl;
            plint numIntCells = fields.lattice.getInternalStatistics().getIntSum(0);
            pcout << "Number of interface cells: " << numIntCells << std::endl;
            pcout << std::endl;
        }

        // This includes the collision-streaming cycle, plus all free-surface operations.
        // Do collision-streaming for lattice 2.
        fields.lattice2->executeInternalProcessors();

        // Do collision-streaming for lattice 1.
        // Execute all the multi-phase data processors.
        fields.lattice.executeInternalProcessors();
        fields.lattice.evaluateStatistics();
        fields.lattice.incrementTime();

        fields.lattice2->evaluateStatistics();
        fields.lattice2->incrementTime();
    }

    return 0;
}

