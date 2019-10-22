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

/* This code was written with help of a Fortran code kindly provided
 * by Prof. Taehun Lee, and contains important contributions by
 * Andrea Parmigiani.
 *
 *  LITERATURE
 *  ==========
 *  - T. Lee and C.-L. Lin, J Comp Phys 206 (2005), 16-47.
 *  - T. Lee, Comp Math App 58 (2009), 987-994.
 */

#include <cmath>
#include "palabos3D.h"
#include "palabos3D.hh"

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR descriptors::D3Q19Descriptor


/**** Begin User Defined Parameters. *****/

plint max_iter;
plint nx, ny, nz;
plint cx0, cy0, cz0; // Center of the drop.
plint radius; // Drop radius.

T delta;      // Interface thickness.
T sigma;      // Surface tension.
T density_ratio, viscosity_ratio;
T RT;         // Speed of sound.
T Pe;         // M*beta;
T M;          // Mobility;
T ini_vel;    // Velocity to kick off bubble after some time steps.

/**** End User Defined Parameters. *****/

T rho_h, rho_l; // Density of heavy and light phase.
T tau_h, tau_l; // Relaxation time of heavy and light phase.
T beta, kappa;
T concentr_h = rho_h;  // Initial concentration for heavy phase (center of the drop).
T concentr_l = 0.;     // Initial concentration for light phase.

MultiBlockLattice3D<T,DESCRIPTOR> *g;  // "The fluid" --> velocity u and pressure p1.
MultiBlockLattice3D<T,DESCRIPTOR> *f;  // "The concentration" --> Conc. C, density rho and
                                       //    potential mu.
MultiScalarField3D<T> *C;   // Concentration of the heavy fluid.
MultiScalarField3D<T> *rho; // Total density, summed over both components.
MultiScalarField3D<T> *mu, *laplaceMu;   // Chemical potential and its Laplacian.
MultiScalarField3D<T> *p1;               // Flow pressure.
MultiTensorField3D<T,3> *gradC, *gradMu; // Gradients of C and mu.
MultiTensorField3D<T,3> *u;              // Flow velocity.

// The arguments to be delivered to the data-processors.
std::vector<MultiBlock3D*> C_processorArguments;
std::vector<MultiBlock3D*> gradC_rho_mu_processorArguments;
std::vector<MultiBlock3D*> gradMu_processorArguments;
std::vector<MultiBlock3D*> u_p_processorArguments;
std::vector<MultiBlock3D*> gradMu_u_p_processorArguments;
std::vector<MultiBlock3D*> heLeeProcessorArguments;


/******* Input Output parameters *******************/
plint getStatistics, getImages, getVTK ;
std::string fNameOut;
/********* Lattices and ScalarField Parallel subdivision *************/
plint blocks_xDir, blocks_yDir, blocks_zDir;

  
void initializeParameters()
{
    
    blocks_xDir=2;
    blocks_yDir=2; 
    blocks_zDir=2;
    
    max_iter = 20000;
    nx = 251;
    ny = 151;
    nz = 151;

    cx0 = (plint)(nx*0.29);
    cy0 = ny/2;
    cz0 = nz/2;
    //radius = 25.;

    //delta   = 5.;
    //sigma   = 1.e-4;
    density_ratio = 10.;
    viscosity_ratio = 10.;
    RT = 1./3.;
    Pe  = 0.02;
    ini_vel = 0.0125;

    rho_h = 1.;
    rho_l = rho_h/density_ratio;

    tau_l = 0.5;
    T nu_l = rho_l*tau_l/3.;         // Viscosity of light phase.
    T nu_h = nu_l*viscosity_ratio;   // Viscosity of heavy phase.
    tau_h = nu_h/rho_h*3.;

    concentr_h = rho_h;
    concentr_l = 0.;

    // Eq.9 from LeePaper. Not clear why he is using concentr_h-concentr_l here.
    beta = 12.*sigma/std::pow((T)(concentr_h-concentr_l),(T)4.)/delta;   
    kappa = beta*util::sqr(delta)*util::sqr(concentr_h-concentr_l)/8.;
    M = Pe/beta;    
    
    getStatistics =    10;
    getImages     =    120;
    getVTK        =    50;
}

MultiBlockLattice3D<T,DESCRIPTOR>* createLattice()
{
    // Envelope of width 2 for the next-to-nearest-neighbor finite
    //   difference schemes.
    plint envelopeWidth = 2;
    SparseBlockStructure3D blockStructure (
            createRegularDistribution3D(nx, ny, nz) );

    MultiBlockLattice3D<T,DESCRIPTOR>* lattice =
            new MultiBlockLattice3D<T,DESCRIPTOR> (
                    MultiBlockManagement3D (
                        blockStructure,
                        defaultMultiBlockPolicy3D().getThreadAttribution(),
                        envelopeWidth ),
                    defaultMultiBlockPolicy3D().getBlockCommunicator(),
                    defaultMultiBlockPolicy3D().getCombinedStatistics(),
                    defaultMultiBlockPolicy3D().getMultiCellAccess<T,DESCRIPTOR>(),
                    new NoDynamics<T,DESCRIPTOR>() );
    lattice->periodicity().toggleAll(true);
    //lattice->periodicity().toggle(true,2);
    return lattice;
}

MultiScalarField3D<T>* createScalarField()
{
    // Envelope of width 2 for the next-to-nearest-neighbor finite
    //   difference schemes.
    plint envelopeWidth = 2;
    SparseBlockStructure3D blockStructure (
            createRegularDistribution3D(nx, ny, nz) );

    MultiScalarField3D<T>* field =
               new MultiScalarField3D<T> (
                    MultiBlockManagement3D (
                        blockStructure,
                        defaultMultiBlockPolicy3D().getThreadAttribution(),
                        envelopeWidth ),
                    defaultMultiBlockPolicy3D().getBlockCommunicator(),
                    defaultMultiBlockPolicy3D().getCombinedStatistics(),
                    defaultMultiBlockPolicy3D().getMultiScalarAccess<T>() );
    field->periodicity().toggleAll(true);
    return field;
}

MultiTensorField3D<T,3>* createVectorField()
{
    // Envelope of width 2 for the next-to-nearest-neighbor finite
    //   difference schemes.
    plint envelopeWidth = 2;
    SparseBlockStructure3D blockStructure (
            createRegularDistribution3D(nx, ny, nz) );

    MultiTensorField3D<T,3>* field =
            new MultiTensorField3D<T,3> (
                    MultiBlockManagement3D (
                        blockStructure,
                        defaultMultiBlockPolicy3D().getThreadAttribution(),
                        envelopeWidth ),
                    defaultMultiBlockPolicy3D().getBlockCommunicator(),
                    defaultMultiBlockPolicy3D().getCombinedStatistics(),
                    defaultMultiBlockPolicy3D().getMultiTensorAccess<T,3>() );
    field->periodicity().toggleAll(true);
    return field;
}

void createFields()
{
    f = createLattice();
    g = createLattice();
    C = createScalarField();
    rho = createScalarField();
    mu = createScalarField();
    laplaceMu = createScalarField();
    p1 = createScalarField();
    gradC = createVectorField();
    gradMu = createVectorField();
    u = createVectorField();
}

void createProcessorArguments()
{
    // Remember that all processors take f as their
    // first argument, even if they don't need it. In this
    // way, the processors can be added to f for automatic
    // execution.

    C_processorArguments.push_back(f);
    C_processorArguments.push_back(laplaceMu);
    C_processorArguments.push_back(C);

    gradC_rho_mu_processorArguments.push_back(f);
    gradC_rho_mu_processorArguments.push_back(C);
    gradC_rho_mu_processorArguments.push_back(gradC);
    gradC_rho_mu_processorArguments.push_back(rho);
    gradC_rho_mu_processorArguments.push_back(mu);

    gradMu_processorArguments.push_back(f);
    gradMu_processorArguments.push_back(mu);
    gradMu_processorArguments.push_back(gradMu);
    gradMu_processorArguments.push_back(laplaceMu);

    u_p_processorArguments.push_back(f);
    u_p_processorArguments.push_back(g);
    u_p_processorArguments.push_back(C);
    u_p_processorArguments.push_back(rho);
    u_p_processorArguments.push_back(gradC);
    u_p_processorArguments.push_back(gradMu);
    u_p_processorArguments.push_back(u);
    u_p_processorArguments.push_back(p1);

    gradMu_u_p_processorArguments.push_back(f);
    gradMu_u_p_processorArguments.push_back(g);
    gradMu_u_p_processorArguments.push_back(C);
    gradMu_u_p_processorArguments.push_back(rho);
    gradMu_u_p_processorArguments.push_back(gradC);
    gradMu_u_p_processorArguments.push_back(mu);
    gradMu_u_p_processorArguments.push_back(gradMu);
    gradMu_u_p_processorArguments.push_back(laplaceMu);
    gradMu_u_p_processorArguments.push_back(u);
    gradMu_u_p_processorArguments.push_back(p1);

    heLeeProcessorArguments.push_back(f);
    heLeeProcessorArguments.push_back(g);
    heLeeProcessorArguments.push_back(C);
    heLeeProcessorArguments.push_back(rho);
    heLeeProcessorArguments.push_back(gradC);
    heLeeProcessorArguments.push_back(mu);
    heLeeProcessorArguments.push_back(gradMu);
    heLeeProcessorArguments.push_back(laplaceMu);
    heLeeProcessorArguments.push_back(u);
    heLeeProcessorArguments.push_back(p1);
}

void addCouplings()
{
    // The data processors must be added at different levels in order to
    //   interleave a communication step between each processor.
    //   It is necessary to start with the level 1, because the level 0
    //   has a special role (after level 0, the communication step would
    //   be called only for f, and not for the other modified blocks).
    integrateProcessingFunctional(
            new Compute_C_processor<T,DESCRIPTOR>(M),
            C->getBoundingBox(),
            C_processorArguments,
            1 );
    integrateProcessingFunctional(
            new Compute_gradC_rho_mu_processor<T>(beta, kappa, rho_h, rho_l),
            f->getBoundingBox(),
            gradC_rho_mu_processorArguments,
            2 );
    integrateProcessingFunctional(
            new Compute_gradMu_laplaceMu_u_p1_processor<T,DESCRIPTOR>(rho_h, rho_l, RT),
            f->getBoundingBox(),
            gradMu_u_p_processorArguments,
            3 );
    integrateProcessingFunctional(
            new HeLeeCollisionProcessor<T,DESCRIPTOR>(rho_h, rho_l, tau_h, tau_l, M, RT),
            f->getBoundingBox(),
            heLeeProcessorArguments,
            4 );
}

void cleanup()
{
    delete f;
    delete g;
    delete C;
    delete rho;
    delete mu;
    delete laplaceMu;
    delete p1;
    delete gradC;
    delete gradMu;
    delete u;
}

// Definition of a drop as initial condition (taken from Lee's Fortran code).
T C_InitialDrop(plint iX, plint iY, plint iZ)
{
    T concentr_mean = 0.5*(concentr_h+concentr_l);
    T concentr_diff = 0.5*(concentr_h-concentr_l);
    T c = (std::sqrt((T)util::sqr(iX-cx0)+(T)util::sqr(iY-cy0)+(T)util::sqr(iZ-cz0))-radius) / delta * 2.;
    return concentr_mean-concentr_diff*std::tanh(c);
}

// Definition of a drop as initial condition (taken from Lee's Fortran code).
T C_InitialEllipsoid(plint iX, plint iY, plint iZ)
{
    T concentr_mean = 0.5*(concentr_h+concentr_l);
    T concentr_diff = 0.5*(concentr_h-concentr_l);
    T c = (std::sqrt(2*(T)util::sqr(iX-cx0)+(T)util::sqr(iY-cy0)+(T)util::sqr(iZ-cz0))-radius) / delta * 2.;
    return concentr_mean-concentr_diff*std::tanh(c);
}

// Definition of a drop as initial condition (taken from Lee's Fortran code).
T C_InitialTwoEllipsoids(plint iX, plint iY, plint iZ)
{
    T concentr_mean = 0.5*(concentr_h+concentr_l);
    T concentr_diff = 0.5*(concentr_h-concentr_l);
    T c1 = (std::sqrt(2*(T)util::sqr(iX-cx0)+(T)util::sqr(iY-cy0)+(T)util::sqr(iZ-cz0))-radius) / delta * 2.;
    T c2 = (std::sqrt((T)util::sqr(iX-nx+cx0)+2*(T)util::sqr(iY-cy0)+(T)util::sqr(iZ-cz0))-radius) / delta * 2.;
    return 2*concentr_mean-concentr_diff*(std::tanh(c1)+std::tanh(c2));
}

void initialCondition()
{
    // 1. Initialize the three macroscopic variables.
    setToConstant(*p1, p1->getBoundingBox(), (T)0.);
    setToConstant<T,3>(*u, Box3D(0,nx/2,0,ny-1,0,nz-1), Array<T,3>((T)0.02,(T)0.,(T)0.));
    setToConstant<T,3>(*u, Box3D(nx/2+1,nx-1,0,ny-1,0,nz-1), Array<T,3>((T)-0.02,(T)0.,(T)0.));
    //setToFunction(*C, C->getBoundingBox(), C_InitialDrop);
    setToFunction(*C, C->getBoundingBox(), C_InitialTwoEllipsoids);

    // 2. Compute the derived quantities.
    applyProcessingFunctional (
            new Compute_gradC_rho_mu_processor<T>(beta, kappa, rho_h, rho_l),
            f->getBoundingBox(),
            gradC_rho_mu_processorArguments );

    applyProcessingFunctional (
            new Compute_gradMu_laplaceMu_processor<T,DESCRIPTOR>(),
            f->getBoundingBox(),
            gradMu_processorArguments );

    // 3. Initialize populations at equilibrium.
    bool onlySetToEquilibrium = true;
    applyProcessingFunctional(
            new HeLeeCollisionProcessor<T,DESCRIPTOR> (
                rho_h, rho_l, tau_h, tau_l, M, RT, onlySetToEquilibrium ),
            f->getBoundingBox(),
            heLeeProcessorArguments );
}

void writeGifs(plint iter)
{
    const plint imSize = 600;
    Box3D slice(0, nx-1, 0, ny-1, nz/2,nz/2);
    ImageWriter<T> imageWriter("leeloo");

    imageWriter.writePpm( createFileName("rho", iter, 6),
                                *extractSubDomain(*rho, slice),
                                imSize, imSize );
    imageWriter.writePpm( createFileName("C", iter, 6),
                                *extractSubDomain(*C, slice),
                                imSize, imSize );
}

void writeVTKscalarField(MultiScalarField3D<T>& scalarField, std::string name, plint iter)
{
 
    T dx = (T)1/(T)(scalarField.getNx()-1);
    // Write full image
    VtkImageOutput3D<T> vtkOut(createFileName("vtk_"+name+"_", iter, 8), dx);
    vtkOut.writeData<float>(scalarField, name,(T)1 );
}

void readCommandLine(int argc, char* argv[])
{
    try {
        global::argv(1).read(sigma);
        global::argv(2).read(delta);
        global::argv(3).read(radius);
        global::argv(4).read(fNameOut);
    }
    catch(...) {
        pcout << "Error : Wrong parameters specified." << endl;
       
        pcout << "1 : surface tension" << endl;
        pcout << "2 : thickness bubble" << endl;
        pcout << "3 : bubble radius" << endl;
        pcout << "4 : output folder name" << endl;
        pcout << "Possible parameters are: " << argv[0] << " 1.e-4  5  15  tmp" << endl; 
        exit(1);
    }
}

int main(int argc, char* argv[]) {

    
    plbInit(&argc, &argv);
    defaultMultiBlockPolicy3D().toggleBlockingCommunication(true);
    
    readCommandLine(argc, argv);
    global::directories().setOutputDir("./"+fNameOut+"/");
    
    initializeParameters();
    createFields();
    createProcessorArguments();
    initialCondition();
    addCouplings();

    pcout << getMultiBlockInfo(*f) << std::endl;
    
    double elapsed=0;
    for (plint iter=0; iter<max_iter; ++iter) {
        global::timer("iteration").restart();
        
        if (iter%getImages==0) {
            if (iter>0) {
                pcout << "Mega Site updates per second: "
                      << (double)g->getBoundingBox().nCells() / elapsed / 1.e6
                      << std::endl;
            }
            pcout << "Iteration step " << iter << endl;
            writeGifs(iter);      
            writeVTKscalarField(*C, "C", iter);
        }
  
        global::timer("sups").restart();
        g->collideAndStream();
        f->collideAndStream();
        elapsed = global::timer("sups").stop();
    }

    cleanup();
}
