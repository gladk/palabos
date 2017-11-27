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
 * Flow around a 2D cylinder inside a channel, with the creation of a von
 * Karman vortex street. This example makes use of bounce-back nodes to
 * describe the shape of the cylinder. The outlet is modeled through a
 * Neumann (zero velocity-gradient) condition.
 */

#include "palabos2D.h"
#include "palabos2D.hh"   // include full template code
#include "complexDataAnalysisWrapper2D.h"
#include "complexDataAnalysisWrapper2D.hh"

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace plb;
using namespace plb::descriptors;
using namespace std;

typedef double T;
typedef plb::Complex<T> PlbT;
#define DESCRIPTOR D2Q9Descriptor

/// Velocity on the parabolic Poiseuille profile
PlbT poiseuilleVelocity(plint iY, IncomprFlowParam<PlbT> const& parameters) {
    PlbT y = (PlbT)iY / (PlbT)parameters.getResolution();
    return 4.*parameters.getLatticeU() * (y-y*y);
}

/// Linearly decreasing pressure profile
PlbT poiseuillePressure(plint iX, IncomprFlowParam<PlbT> const& parameters) {
    PlbT Lx = parameters.getNx()-1;
    PlbT Ly = parameters.getNy()-1;
    return 8.*parameters.getLatticeNu()*parameters.getLatticeU() / (Ly*Ly) * (Lx/(T)2-(T)iX);
}

/// Convert pressure to density according to ideal gas law
PlbT poiseuilleDensity(plint iX, IncomprFlowParam<PlbT> const& parameters) {
    return poiseuillePressure(iX,parameters)*DESCRIPTOR<PlbT>::invCs2 + (PlbT)1;
}

/// A functional, used to initialize the velocity for the boundary conditions
template<typename PlbT>
class PoiseuilleVelocity {
public:
    PoiseuilleVelocity(IncomprFlowParam<PlbT> parameters_)
        : parameters(parameters_)
    { }
    void operator()(plint iX, plint iY, Array<PlbT,2>& u) const {
        u[0] = poiseuilleVelocity(iY, parameters);
        u[1] = PlbT();
    }
private:
    IncomprFlowParam<PlbT> parameters;
};

/// A functional, used to initialize a pressure boundary to constant density
template<typename PlbT>
class ConstantDensity {
public:
    ConstantDensity(PlbT density_)
        : density(density_)
    { }
    PlbT operator()(plint iX, plint iY) const {
        return density;
    }
private:
    PlbT density;
};

/// A functional, used to create an initial condition for the density and velocity
template<typename PlbT>
class PoiseuilleVelocityAndDensity {
public:
    PoiseuilleVelocityAndDensity(IncomprFlowParam<PlbT> parameters_)
        : parameters(parameters_)
    { }
    void operator()(plint iX, plint iY, PlbT& rho, Array<PlbT,2>& u) const {
        rho = poiseuilleDensity(iX,parameters);
        u[0] = poiseuilleVelocity(iY, parameters);
        u[1] = PlbT();
    }
private:
    IncomprFlowParam<PlbT> parameters;
};

template<typename PlbT>
class CylinderShapeDomain2D : public plb::DomainFunctional2D {
public:
    CylinderShapeDomain2D(plb::plint cx_, plb::plint cy_, plb::plint radius)
        : cx(cx_),
          cy(cy_),
          radiusSqr(plb::util::sqr(radius))
    { }
    virtual bool operator() (plb::plint iX, plb::plint iY) const {
        return plb::util::sqr(iX-cx) + plb::util::sqr(iY-cy) <= radiusSqr;
    }
    virtual CylinderShapeDomain2D<PlbT>* clone() const {
        return new CylinderShapeDomain2D<PlbT>(*this);
    }
private:
    plb::plint cx;
    plb::plint cy;
    plb::plint radiusSqr;
};


/// A functional, used to instantiate bounce-back nodes at the locations of the cylinder
void cylinderSetup( MultiBlockLattice2D<PlbT,DESCRIPTOR>& lattice,
                    IncomprFlowParam<PlbT> const& parameters,
                    OnLatticeBoundaryCondition2D<PlbT,DESCRIPTOR>& boundaryCondition )
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    Box2D outlet(nx-1,nx-1, 1, ny-2);

    // Create Velocity boundary conditions everywhere
    boundaryCondition.setVelocityConditionOnBlockBoundaries (
            lattice, Box2D(0, 0, 1, ny-2) );
    boundaryCondition.setVelocityConditionOnBlockBoundaries (
            lattice, Box2D(0, nx-1, 0, 0) );
    boundaryCondition.setVelocityConditionOnBlockBoundaries (
            lattice, Box2D(0, nx-1, ny-1, ny-1) );
    // .. except on right boundary, where we prefer an outflow condition
    //    (zero velocity-gradient).
    boundaryCondition.setVelocityConditionOnBlockBoundaries (
            lattice, Box2D(nx-1, nx-1, 1, ny-2), boundary::outflow );

    setBoundaryVelocity (
            lattice, lattice.getBoundingBox(),
            PoiseuilleVelocity<PlbT>(parameters) );
    setBoundaryDensity (
            lattice, outlet,
            ConstantDensity<PlbT>(1.) );
    initializeAtEquilibrium (
            lattice, lattice.getBoundingBox(),
            PoiseuilleVelocityAndDensity<PlbT>(parameters) );

    plint cx     = nx/4;
    plint cy     = ny/2+2; // cy is slightly offset to avoid full symmetry,
                          //   and to get a Von Karman Vortex street.
    plint radius = cy/4;
    defineDynamics(lattice, lattice.getBoundingBox(),
                   new CylinderShapeDomain2D<T>(cx,cy,radius),
                   new plb::BounceBack<PlbT,DESCRIPTOR>);

    lattice.initialize();
}

/*
void writeGif(MultiBlockLattice2D<PlbT,DESCRIPTOR>& lattice, plint iter)
{
    ImageWriter<PlbT> imageWriter("leeloo");
    imageWriter.writeScaledGif(createFileName("u", iter, 6),
                               *computeVelocityNorm(lattice) );
}
*/
void writeVTK(MultiBlockLattice2D<PlbT,DESCRIPTOR>& lattice,
              IncomprFlowParam<PlbT> const& parameters, plint iter)
{
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();
    VtkImageOutput2D<T> vtkOut(createFileName("vtk", iter, 6), dx);
    vtkOut.writeData<T>(*realPart<PlbT,T>(*computeVelocityNorm(lattice),lattice.getBoundingBox()), "velocityNorm", dx/dt);
    //vtkOut.writeData<2,T>(*computeVelocity(lattice), "velocity", dx/dt);
}

// void writeVTS(MultiBlockLattice2D<PlbT,DESCRIPTOR>& lattice,
//               IncomprFlowParam<PlbT> const& parameters, plint iter)
// {
//     T dx = parameters.getDeltaX();
//     T dt = parameters.getDeltaT();
//     VtkStructuredOutput2D<T> vtkOut(createFileName("vts", iter, 6), dx);
//     vtkOut.writeData<T>(*realPart<PlbT,T>(*computeVelocityNorm(lattice),lattice.getBoundingBox()), "velocityNormRe", dx/dt);
// 	vtkOut.writeData<T>(*imaginaryPart<PlbT,T>(*computeVelocityNorm(lattice),lattice.getBoundingBox()), "velocityNormIm", dx/dt);
//     vtkOut.writeData<2,T>(*realPart<PlbT,T,2>(*computeVelocity(lattice),lattice.getBoundingBox()), "velocityRe", dx/dt);
// 	vtkOut.writeData<2,T>(*imaginaryPart<PlbT,T,2>(*computeVelocity(lattice),lattice.getBoundingBox()), "velocityIm", dx/dt);
// }


int main(int argc, char* argv[]) {
    plbInit(&argc, &argv);

    global::directories().setOutputDir("./tmp/");
    PlbT Re(600.,600*1e-3);
	pcout << Re.real() <<"  "<<Re.imaginary()<<endl;
	IncomprFlowParam<PlbT> parameters(
            (PlbT) 1e-2,  // uMax
            (PlbT) Re,  // Re
            100,       // N
            6.,        // lx
            1.         // ly 
    );
    const T logT     = (T)0.02;
    const T imSave   = (T)0.06;
    const T vtkSave  = (T)0.06;
    const T maxT     = (T)20.1;
	pcout<<"  " << parameters.getDeltaT().imaginary()<<endl;

    //writeLogFile(parameters, "Poiseuille flow");

    MultiBlockLattice2D<PlbT, DESCRIPTOR> lattice (
            parameters.getNx(), parameters.getNy(),
            new BGKdynamics<PlbT,DESCRIPTOR>(parameters.getOmega()) );
	pcout << parameters.getOmega().real()<<"   "<< parameters.getOmega().imaginary()<< endl;

	

    OnLatticeBoundaryCondition2D<PlbT,DESCRIPTOR>*
        boundaryCondition = createLocalBoundaryCondition2D<PlbT,DESCRIPTOR>();

    cylinderSetup(lattice, parameters, *boundaryCondition);
	
	//Array<PlbT,2> velocity;
//    lattice.get(10, 10).computeVelocity(velocity);
//    pcout << "Velocity in the middle of the lattice: ("
//	<< velocity[0].real()<<" "<<velocity[0].imaginary() << "," << velocity[1].real()<<" "<<velocity[1].imaginary()  << ")" << endl;

    // Main loop over time iterations.
    for (plint iT=0; (T)iT*parameters.getDeltaT().real()<maxT; ++iT) {
        // At this point, the state of the lattice corresponds to the
        //   discrete time iT. However, the stored averages (getStoredAverageEnergy
        //   and getStoredAverageDensity) correspond to the previous time iT-1.

       if (iT%parameters.nStep(imSave)==0) {
            //pcout << "Saving Gif ..." << endl;
            //writeGif(lattice, iT);
        }

        //if (iT%parameters.nStep(vtkSave)==0 && iT>=0) {
        //    pcout << "Saving VTK file ..." << endl;
        //    writeVTK(lattice, parameters, iT);
        //}
		if (iT%parameters.nStep(vtkSave)==0 && iT>=0) {
            pcout << "Saving VTK file ..." << endl;
            writeVTK(lattice, parameters, iT);
			Array<PlbT,2> velocity;
		    lattice.get(10, 10).computeVelocity(velocity);
			pcout << "Velocity in the middle of the lattice: ("
			      << velocity[0].real()<<" "<<velocity[0].imaginary() << "," << velocity[1].real()<<" "<<velocity[1].imaginary()  << ")" << endl;
        }

        if (iT%parameters.nStep(logT)==0) {
            pcout << "step " << iT
                  << "; t=" << (T)iT*parameters.getDeltaT().real();
        }

        // Lattice Boltzmann iteration step.
        lattice.collideAndStream();
//		lattice.get(10, 10).computeVelocity(velocity);
//		pcout << "Velocity in the middle of the lattice: ("
//		<< velocity[0].real()<<" "<<velocity[0].imaginary() << "," << velocity[1].real()<<" "<<velocity[1].imaginary()  << ")" << endl;

        // At this point, the state of the lattice corresponds to the
        //   discrete time iT+1, and the stored averages are upgraded to time iT.
        if (iT%parameters.nStep(logT)==0) {
            pcout << "; av energy ="
                  << setprecision(10) << getStoredAverageEnergy<PlbT>(lattice)
                  << "; av rho ="
                  << getStoredAverageDensity<PlbT>(lattice) << endl;
        }
    }
    
    delete boundaryCondition;
}
