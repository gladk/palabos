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

/* Main author: Orestis Malaspinas
 */

/** \file
 * A fluid constrained between a hot bottom wall (no-slip for the velocity) and a cold
 * top wall (no-slip for the velocity). The lateral walls are periodic. Under the
 * influence of gravity, convection rolls are formed. Thermal effects are modelled
 * by means of a Boussinesq approximation: the fluid is incompressible, and the influence
 * of the temperature is visible only through a body-force term, representing buoyancy
 * effects. The temperature field obeys an advection-diffusion equation.
 *
 * The simulation is first created in a fully symmetric manner. The symmetry is therefore
 * not spontaneously broken; while the temperature drops linearly between the hot and
 * and cold wall, the convection rolls fail to appear at this point. In a second stage, a
 * random noise is added to trigger the instability.
 *
 * This application is technically a bit more advanced than the other ones, because it
 * illustrates the concept of data processors. In the present case, they are used to
 * create the initial condition, and to trigger the instability.
 **/

#include "palabos2D.h"
#include "palabos2D.hh"

#include <cstdlib>
#include <iostream>

using namespace plb;
using namespace std;

typedef double T;

#define NSDESCRIPTOR descriptors::ForcedD2Q9Descriptor
#define ADESCRIPTOR descriptors::AdvectionDiffusionD2Q5Descriptor

#define ADYNAMICS AdvectionDiffusionBGKdynamics
#define NSDYNAMICS GuoExternalForceBGKdynamics


/// Initialization of the temperature field.
template<typename T, template<typename NSU> class nsDescriptor,
                     template<typename ADU> class adDescriptor>
struct IniTemperatureRayleighBenardProcessor2D :
    public BoxProcessingFunctional2D_L<T,adDescriptor> 
{
    IniTemperatureRayleighBenardProcessor2D(RayleighBenardFlowParam<T,nsDescriptor,adDescriptor> parameters_)
        : parameters(parameters_)
    { }
    virtual void process(Box2D domain, BlockLattice2D<T,adDescriptor>& adLattice)
    {
        Dot2D absoluteOffset = adLattice.getLocation();
        
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                plint absoluteY = absoluteOffset.y + iY;
                T temperature = parameters.getHotTemperature() 
                        - parameters.getDeltaTemperature() / (T)(parameters.getNy()-1) * (T)absoluteY;
                
                Array<T,adDescriptor<T>::d> jEq(0.0, 0.0);
                adLattice.get(iX,iY).defineDensity(temperature);
                iniCellAtEquilibrium(adLattice.get(iX,iY), temperature, jEq);
            }
        }
    }
    virtual IniTemperatureRayleighBenardProcessor2D<T,nsDescriptor,adDescriptor>* clone() const
    {
        return new IniTemperatureRayleighBenardProcessor2D<T,nsDescriptor,adDescriptor>(*this);
    }
    
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
    }
    
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulkAndEnvelope;
    }
    
private :
    RayleighBenardFlowParam<T,nsDescriptor,adDescriptor> parameters;
};

/// Perturbation of the temperature field to instantiate the instability.
template<typename T, template<typename NSU> class nsDescriptor,
                     template<typename ADU> class adDescriptor>
struct PerturbTemperatureRayleighBenardProcessor2D : public BoxProcessingFunctional2D_L<T,adDescriptor> 
{
    PerturbTemperatureRayleighBenardProcessor2D (
            RayleighBenardFlowParam<T,nsDescriptor,adDescriptor> parameters_ )
        : parameters(parameters_)
    { }
    virtual void process(Box2D domain,BlockLattice2D<T,adDescriptor>& lattice)
    {
        Dot2D absoluteOffset = lattice.getLocation();
        
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                plint absoluteX = absoluteOffset.x + iX;
                plint absoluteY = absoluteOffset.y + iY;
                
                if ((absoluteX == (parameters.getNx()-1)/2) && (absoluteY == 1))
                {
                    T temperature = T();
                    temperature = parameters.getHotTemperature() * 1.1;
                    
                    Array<T,adDescriptor<T>::d> jEq(0.0,0.0);
                    lattice.get(iX,iY).defineDensity(temperature);
                    iniCellAtEquilibrium(lattice.get(iX,iY), temperature, jEq);
                }
            }
        }
    }
    virtual PerturbTemperatureRayleighBenardProcessor2D<T,nsDescriptor,adDescriptor>* clone() const
    {
        return new PerturbTemperatureRayleighBenardProcessor2D<T,nsDescriptor,adDescriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
    }
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulk;
    }
private :
    RayleighBenardFlowParam<T,nsDescriptor,adDescriptor> parameters;
};

void rayleighBenardSetup (
              MultiBlockLattice2D<T, NSDESCRIPTOR>& nsLattice,
              MultiBlockLattice2D<T, ADESCRIPTOR>& adLattice,
              OnLatticeBoundaryCondition2D<T,NSDESCRIPTOR>& nsBoundaryCondition,
              OnLatticeAdvectionDiffusionBoundaryCondition2D<T,ADESCRIPTOR>& adBoundaryCondition,
              RayleighBenardFlowParam<T,NSDESCRIPTOR,ADESCRIPTOR> &parameters )
{
    plint nx = parameters.getNx();
    plint ny = parameters.getNy();
    
    Box2D top   (0, nx-1, ny-1, ny-1);
    Box2D bottom(0, nx-1, 0,    0);
    
    nsBoundaryCondition.addVelocityBoundary1N(bottom, nsLattice);
    nsBoundaryCondition.addVelocityBoundary1P(top,    nsLattice);
    
    adBoundaryCondition.addTemperatureBoundary1N(bottom, adLattice);
    adBoundaryCondition.addTemperatureBoundary1P(top,    adLattice);

    initializeAtEquilibrium(nsLattice, nsLattice.getBoundingBox(), (T)1., Array<T,2>((T)0.,(T)0.) );
    
    applyProcessingFunctional (
            new IniTemperatureRayleighBenardProcessor2D<T,NSDESCRIPTOR,ADESCRIPTOR>(parameters), 
            adLattice.getBoundingBox(),  adLattice );
    
    nsLattice.initialize();
    adLattice.initialize();
    
}

void writeVTK(MultiBlockLattice2D<T,NSDESCRIPTOR>& nsLattice,
              MultiBlockLattice2D<T,ADESCRIPTOR>& adLattice,
              RayleighBenardFlowParam<T,NSDESCRIPTOR,ADESCRIPTOR> const& parameters, plint iter)
{
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();
    
    VtkImageOutput2D<T> vtkOut(createFileName("vtk", iter, 6), dx);
    // Temperature is the order-0 moment of the advection-diffusion model. It can 
    //    therefore be computed with the function "computeDensity".
    vtkOut.writeData<float>(*computeDensity(adLattice), "temperature", (T)1);
    vtkOut.writeData<2,float>(*computeVelocity(nsLattice), "velocity", dx/dt);
}

void writeGif(MultiBlockLattice2D<T,NSDESCRIPTOR>& nsLattice,
              MultiBlockLattice2D<T,ADESCRIPTOR>& adLattice,int iT)
{
    const plint imSize = 600;
    const plint nx = nsLattice.getNx();
    const plint ny = nsLattice.getNy();
    Box2D slice(0, nx-1, 0, ny-1);
    ImageWriter<T> imageWriter("leeloo.map");
    imageWriter.writeScaledGif(createFileName("u", iT, 6),
                               *computeVelocityNorm(nsLattice, slice),
                               imSize, imSize);
    // Temperature is the order-0 moment of the advection-diffusion model. It can 
    //    therefore be computed with the function "computeDensity".
    imageWriter.writeScaledGif(createFileName("temperature", iT, 6),
                               *computeDensity(adLattice, slice),
                               imSize, imSize);
}

int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);
    
    global::timer("simTime").start();
    
    T Ra=0.;
    try {
        global::argv(1).read(Ra);
    }
    catch(PlbIOException& exception) {
        pcout << exception.what() << endl;
        pcout << "The structure of the input parameters should be : "
              << (string)global::argv(0) << " Ra" << endl;;
        // Exit the program, because wrong input data is a fatal error.
        exit(1);
    }
    
    const T lx  = 2.0;
    const T ly  = 1.0;
    const T uMax  = 0.1;
    const T Pr = 1.0;
    
    const T hotTemperature = 1.0;
    const T coldTemperature = 0.0;
    const plint resolution = 50;

    global::directories().setOutputDir("./tmp/");
    
    RayleighBenardFlowParam<T,NSDESCRIPTOR,ADESCRIPTOR> parameters (
         Ra, 
         Pr, 
         uMax,
         coldTemperature,
         hotTemperature, 
         resolution, 
         lx, 
         ly
   );
                                        
    writeLogFile(parameters,"palabos.log");
    
    const double rayleigh = parameters.getResolution() * parameters.getResolution() * 
            parameters.getResolution() * parameters.getDeltaTemperature() * 
            parameters.getLatticeGravity() / (parameters.getLatticeNu()*parameters.getLatticeKappa());

    const double prandtl = parameters.getLatticeNu() / parameters.getLatticeKappa();

    pcout << "Ra:" << rayleigh << "; Pr:" << prandtl << endl;
        
    plint nx = parameters.getNx();
    plint ny = parameters.getNy();
    
    T nsOmega = parameters.getSolventOmega();
    T adOmega = parameters.getTemperatureOmega();
    
    MultiBlockLattice2D<T, NSDESCRIPTOR> nsLattice (
            nx,ny,new NSDYNAMICS<T, NSDESCRIPTOR>(nsOmega) );
    // Use periodic boundary conditions.
    nsLattice.periodicity().toggleAll(true);
            
    MultiBlockLattice2D<T, ADESCRIPTOR> adLattice (
            nx,ny, new ADYNAMICS<T, ADESCRIPTOR>(adOmega) );
    // Use periodic boundary conditions.
    adLattice.periodicity().toggleAll(true);
            
    OnLatticeBoundaryCondition2D<T,NSDESCRIPTOR>*
        nsBoundaryCondition = createLocalBoundaryCondition2D<T,NSDESCRIPTOR>();
        
    OnLatticeAdvectionDiffusionBoundaryCondition2D<T,ADESCRIPTOR>*
        adBoundaryCondition = createLocalAdvectionDiffusionBoundaryCondition2D<T,ADESCRIPTOR>();
    
    // Turn off internal statistics in order to improve parallel efficiency
    //   (it is not used anyway).
    nsLattice.toggleInternalStatistics(false);
    adLattice.toggleInternalStatistics(false);

    rayleighBenardSetup(nsLattice, adLattice,*nsBoundaryCondition, *adBoundaryCondition, parameters);
    
    Array<T,NSDESCRIPTOR<T>::d> forceOrientation(T(),(T)1);
    plint processorLevel = 1;
    integrateProcessingFunctional (
            new BoussinesqThermalProcessor2D<T,NSDESCRIPTOR,ADESCRIPTOR> (
                parameters.getLatticeGravity(), parameters.getAverageTemperature(),
                parameters.getDeltaTemperature(),forceOrientation ),
            nsLattice.getBoundingBox(),
            nsLattice, adLattice, processorLevel );
    
    T tIni = global::timer("simTime").stop();
    pcout << "time elapsed for rayleighBenardSetup:" << tIni << endl;
    global::timer("simTime").start();
    
    plint evalTime =10000;
    plint iT = 0;
    plint maxT = 1000000000;
    plint statIter = 100;
    plint saveIter = 10000;
    util::ValueTracer<T> converge((T)1,(T)100,1.0e-6);
    bool firstTimeConverged = false;

    // Main loop over time iterations.
    for (iT = 0; iT <= maxT; ++iT) 
    {
        if (iT == (evalTime))
        {
            T tEval = global::timer("simTime").stop();
            T remainTime = (tEval - tIni) / (T)evalTime * (T)maxT/(T)3600;
            global::timer("simTime").start();
            pcout << "Remaining " << (plint)remainTime << " hours, and ";
            pcout << (plint)((T)60*(remainTime - (T)((plint)remainTime))+0.5) << " minutes." << endl;
        }
        if (iT % statIter == 0)
        {
            int yDirection = 1;
            T nusselt = computeNusseltNumber (
                       nsLattice, adLattice, nsLattice.getBoundingBox(),
                       yDirection, parameters.getDeltaX(), 
                       parameters.getLatticeKappa(), parameters.getDeltaTemperature() );
            converge.takeValue(nusselt,true);
        }
        if (converge.hasConverged())
        {
            if (!firstTimeConverged)
            {
                firstTimeConverged = true;
                converge.resetValues();
                converge.setEpsilon(1.0e-14);
                applyProcessingFunctional (
                    new PerturbTemperatureRayleighBenardProcessor2D<T,NSDESCRIPTOR,ADESCRIPTOR>(parameters), 
                    adLattice.getBoundingBox(), adLattice);
                pcout << "Intermediate convergence.\n";
            }
            else
            {
                pcout << "Simulation is over.\n";
                break;
            }
        }
        if (iT % saveIter == 0)
        {
            pcout << "At time " << iT * parameters.getDeltaT() << std::endl;
            pcout << "Writing VTK..." << endl;
            writeVTK(nsLattice, adLattice, parameters, iT);
          
            pcout << "Writing gif..." << endl;
            writeGif(nsLattice,adLattice,iT);
        }
        
        // Lattice Boltzmann iteration step.
        adLattice.collideAndStream();
        nsLattice.collideAndStream();
    }
    
    writeGif(nsLattice,adLattice,iT);
    
    T tEnd = global::timer("simTime").stop();
    
    T totalTime = tEnd-tIni;
    T nx1000 = nsLattice.getNx()/(T)1000;
    T ny1000 = nsLattice.getNy()/(T)1000;
    pcout << "N=" << resolution << endl;
    pcout << "number of processors: " << global::mpi().getSize() << endl;
    pcout << "simulation time: " << totalTime << endl;
    pcout << "total time: " << tEnd << endl;
    pcout << "total iterations: " << iT << endl;
    pcout << "Msus: " << nx1000*ny1000*(T)iT/totalTime << endl;
}
