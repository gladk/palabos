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

#include "palabos2D.h"
#include "palabos2D.hh"

#include <cstdlib>
#include <iostream>
#include <complex>

using namespace plb;
using namespace std;

typedef double T;
const T pi = (T)4.0*std::atan((T)1.0);

#define NSDESCRIPTOR descriptors::ForcedD2Q9Descriptor
#define DYNAMICS GuoExternalForceBGKdynamics<T, NSDESCRIPTOR>(omega)
// #define DYNAMICS HeExternalForceBGKdynamics<T, NSDESCRIPTOR>(omega)
// #define DYNAMICS ShanExternalForceBGKdynamics<T, NSDESCRIPTOR>(omega)

/// Velocity on the parabolic Womersley profile
T womersleyVelocity(plint iY, T t, T A, T omega, T alpha, IncomprFlowParam<T> const& parameters) 
{
    const complex<T> I(0.,1.);
    T y = (T)iY / parameters.getResolution();

    return ( A / (I * omega) * std::exp(I * omega * t) * 
            ( (T)1. - std::cosh(std::sqrt((T)2.)*(y-(T).5)*(alpha+I*alpha)) /
                std::cosh(std::sqrt((T)2.)/(T)2. * (alpha+I*alpha)) ) ).real();
}

/// Time dependent but space-independent Womersley force.
T womersleyForce(T t, T A, T omega, IncomprFlowParam<T> const& parameters) 
{
    return A * std::cos(omega * t);
}

/// A functional, used to initialize the velocity for the boundary conditions
template<typename T>
class WomersleyVelocity {
public:
    WomersleyVelocity(IncomprFlowParam<T> parameters_, T alpha_, T t_)
    :   parameters(parameters_), alpha(alpha_),  
        omega((T)4*util::sqr(alpha)*parameters.getLatticeNu()/(T)(util::sqr(parameters.getResolution()))), 
        A(8. * parameters.getLatticeNu() * parameters.getLatticeU() / ( (T)(util::sqr(parameters.getResolution())) )),
        t(t_)
    { }
    
    void operator()(plint iX, plint iY, Array<T,2>& u) const {
        u[0] = womersleyVelocity(iY, t, A, omega, alpha, parameters);
        u[1] = T();
    }
private:
    IncomprFlowParam<T> parameters;
    T alpha, omega, A;
    T t;
};

void channelSetup( MultiBlockLattice2D<T,NSDESCRIPTOR>& lattice,
                   IncomprFlowParam<T> const& parameters,
                   OnLatticeBoundaryCondition2D<T,NSDESCRIPTOR>& boundaryCondition,
                   T alpha, T frequency, T amplitude)
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    
    Box2D bottom(   0,nx-1,   0,   0);
    Box2D top(   0,nx-1,   ny-1,   ny-1);
    
    boundaryCondition.addVelocityBoundary1N(bottom, lattice);
    boundaryCondition.addPressureBoundary1P(top,    lattice);
    
    Array<T,2> u((T)0.,(T)0.);
    setBoundaryVelocity( lattice, lattice.getBoundingBox(), u );
    initializeAtEquilibrium(lattice,lattice.getBoundingBox(),(T)1.0,u);

    Array<T,NSDESCRIPTOR<T>::d> force(womersleyForce((T)0, amplitude, frequency, parameters),0.);
    setExternalVector(lattice,lattice.getBoundingBox(),NSDESCRIPTOR<T>::ExternalField::forceBeginsAt,force);
    
    lattice.initialize();
}

template<class BlockLatticeT>
void writeVTK(BlockLatticeT& lattice,
            IncomprFlowParam<T> const& parameters, plint iter)
{
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();

    VtkImageOutput2D<T> vtkOut(createFileName("vtk", iter, 6), dx);
    vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", dx/dt);
}

template<class BlockLatticeT>
void writeGif(BlockLatticeT& lattice,plint iT)
{
    const plint imSize = 600;
    ImageWriter<T> imageWriter("leeloo.map");
    imageWriter.writeScaledGif(createFileName("u", iT, 6),
                            *computeVelocityNorm(lattice),
                            imSize, imSize);
}

T computeRMSerror ( MultiBlockLattice2D<T,NSDESCRIPTOR>& lattice,
                    IncomprFlowParam<T> const& parameters,
                    T alpha, plint iT, bool createImage=false)
{
    MultiTensorField2D<T,2> analyticalVelocity(lattice);
    setToFunction( analyticalVelocity, analyticalVelocity.getBoundingBox(),
                   WomersleyVelocity<T>(parameters,alpha,(T)iT) );
    MultiTensorField2D<T,2> numericalVelocity(lattice);
    computeVelocity(lattice, numericalVelocity, lattice.getBoundingBox());

           // Divide by lattice velocity to normalize the error
    return 1./parameters.getLatticeU() *
           // Compute RMS difference between analytical and numerical solution
               std::sqrt( computeAverage( *computeNormSqr (
                              *subtract(analyticalVelocity, numericalVelocity)
                         ) ) );
}

int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);

    if (argc != 2)
    {
        pcout << "Error : Wrong parameters specified." << endl;
        pcout << "1 : N." << endl;
        exit(1);
    }

    const plint N = atoi(argv[1]);
    
    const T Re = 1.0;
    const T alpha = 1.0; // womersley number

    const plint Nref = 10;

    const T uMaxRef = 0.01;

    const T uMax = uMaxRef /(T)N * (T)Nref; // needed to avoid compressibility errors.

    const T lx  = 100.0;
    const T ly  = 1.0;
    pcout << "uMaxRef=" << uMaxRef << std::endl;
    pcout << "uMax=" << uMax << std::endl;

    global::directories().setOutputDir("./tmp/");

    IncomprFlowParam<T> parameters(uMax, Re, N, lx, ly);
    
//     The frequency of the force (lattice units)
    T frequency = (T)4*alpha*alpha*parameters.getLatticeNu() 
                   / (T)(parameters.getResolution()*parameters.getResolution());
                   
//     The amplitude of the forcing term (lattice units)
    T amplitude = 8. * parameters.getLatticeNu() * parameters.getLatticeU() 
                   / ( (T)(parameters.getResolution()*parameters.getResolution()) );
                   
//     Period of the force (lattice units)
    plint tPeriod = (plint)((T)2*pi/frequency + 0.5);

    writeLogFile(parameters,"palabos.log");

    plint nx = parameters.getNx();
    plint ny = parameters.getNy();

    T omega = parameters.getOmega();

    MultiBlockLattice2D<T, NSDESCRIPTOR> lattice (
            nx,ny,new DYNAMICS );
            
    OnLatticeBoundaryCondition2D<T,NSDESCRIPTOR>*
        boundaryCondition = createLocalBoundaryCondition2D<T,NSDESCRIPTOR>();

    lattice.periodicity().toggle(0,true);
    
    channelSetup( lattice, parameters, *boundaryCondition, alpha, frequency, amplitude);

    pcout << "Starting simulation" << endl;

    const plint maxIter = tPeriod * 100;
    //const plint tSave = tPeriod / 24;

    
    T error = T();
    
    lattice.resetTime(1);
    
    pcout << "Omega = " << omega << ", it period = " << tPeriod << endl;

    util::ValueTracer<T> converge(uMax,N,1.0e-3);
    plint iT = 0;
    for (iT = 0; iT < maxIter; ++iT) {
//         Updating the force in the whole domain
        Array<T,NSDESCRIPTOR<T>::d> force(womersleyForce((T)iT, amplitude, frequency, parameters),0.);
        setExternalVector(lattice,lattice.getBoundingBox(),
                          NSDESCRIPTOR<T>::ExternalField::forceBeginsAt,force);
        
        T errorTemp = computeRMSerror( lattice,parameters,alpha,iT);
        error += errorTemp;
        
        //if (iT % tSave == 0) {
        //    pcout << "Writing Gif at time : " << iT << std::endl;
        //    writeGif(lattice,iT);
        //}
        
        if (iT % tPeriod == 0)
        {
//             The error is averaged over one period
            error /= (T)(tPeriod);
            pcout << "For N = " << N << ", Error = " << error << endl;
            converge.takeValue(error,true);
            if (converge.hasConverged())
            {
                cout << "Simulation converged!\n";
                break;
            }
            error = T();
        }
        
        lattice.collideAndStream();
    }
    
    pcout << "For N = " << N << ", Error = " << computeRMSerror ( lattice,parameters,alpha,iT, true) << endl;
}

