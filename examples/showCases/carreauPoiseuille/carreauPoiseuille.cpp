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
#ifndef PLB_PRECOMPILED // Unless precompiled version is used,
#include "palabos2D.hh"   // include full template code
#endif

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "functions.h"
#include "newtonRaphson.h"
#include "trapeziumIntegration.h"

using namespace plb;
using namespace plb::descriptors;
using namespace std;

typedef double T;
#define DESCRIPTOR D2Q9Descriptor

template <typename T>
T computePressureGradient(CarreauFlowParam<T> parameters, T tol, plint maxIter)
{
    T kappa0 = -1.0;
    T kappa = (T)10 * kappa0;
    T uMax = (T)1;

    for (plint iPop = 0; iPop < maxIter; ++iPop)
    {
        TrapeziumIntegration<T> ti0(new NewtonRaphson<T>(
                                           new CarreauFunction_2<T>(parameters.getLatticeNu0()/
                                         (parameters.getDeltaT()/(parameters.getDeltaX()*parameters.getDeltaX())),
                                          parameters.getLatticeLambda()*parameters.getDeltaT(),
                                          parameters.getExponent(),
                                           kappa0,parameters.getLy()),tol,maxIter) , 0.0,
                                          (parameters.getNy()-1)*maxIter);
        T uMaxTmp0 = ti0(parameters.getLy()/(T)2);

        TrapeziumIntegration<T> ti(new NewtonRaphson<T>(
                                           new CarreauFunction_2<T>(parameters.getLatticeNu0()/
                                         (parameters.getDeltaT()/(parameters.getDeltaX()*parameters.getDeltaX())),
                                          parameters.getLatticeLambda()*parameters.getDeltaT(),
                                          parameters.getExponent(),
                                           kappa,parameters.getLy()),tol,maxIter) , 0.0,
                                          (parameters.getNy()-1)*maxIter);
        T uMaxTmp = ti(parameters.getLy()/(T)2);

        if ( ((uMaxTmp - uMax) > 0 && (uMaxTmp0 - uMax) < 0) || ((uMaxTmp - uMax) < 0 && (uMaxTmp0 - uMax) > 0)) {
            break;
        }
        else if ((uMaxTmp - uMax) > 0 && (uMaxTmp0 - uMax) > 0) {
            kappa = kappa / (T)10;
        }
        else if ((uMaxTmp - uMax) < 0 && (uMaxTmp0 - uMax) < 0) {
            kappa = kappa * (T)10;
        }
    }

    for (plint iPop = 0; iPop < maxIter; ++iPop)
    {
        TrapeziumIntegration<T> ti(new NewtonRaphson<T>(
                                       new CarreauFunction_2<T>(parameters.getLatticeNu0()/
                                         (parameters.getDeltaT()/(parameters.getDeltaX()*parameters.getDeltaX())),
                                          parameters.getLatticeLambda()*parameters.getDeltaT(),
                                          parameters.getExponent(),
                                       kappa,parameters.getLy()),tol,maxIter) , 0.0,
                                          (parameters.getNy()-1)*maxIter);

        T uMaxTmp = ti(parameters.getLy()/(T)2);
        if (std::fabs(uMaxTmp - uMax) < tol) {
            return kappa;
        }
        else {
            if (uMaxTmp > uMax) {
                T kappaTmp = kappa;
                kappa -= (kappa - kappa0)/(T)2;
                TrapeziumIntegration<T> ti0(new NewtonRaphson<T>(
                                       new CarreauFunction_2<T>(parameters.getLatticeNu0()/
                                         (parameters.getDeltaT()/(parameters.getDeltaX()*parameters.getDeltaX())),
                                          parameters.getLatticeLambda()*parameters.getDeltaT(),
                                          parameters.getExponent(),
                                          kappa,parameters.getLy()),tol,maxIter) , 0.0,
                                          (parameters.getNy()-1)*maxIter);
                T uMaxTmp0 = ti0(parameters.getLy()/(T)2);
                if (uMaxTmp0 < uMax) {
                    kappa0 = kappaTmp;
                }
            }
            else {
                T kappaTmp = kappa0;
                kappa += (kappa - kappa0)/(T)2;
                TrapeziumIntegration<T> ti0(new NewtonRaphson<T>(
                                       new CarreauFunction_2<T>(parameters.getLatticeNu0()/
                                         (parameters.getDeltaT()/(parameters.getDeltaX()*parameters.getDeltaX())),
                                          parameters.getLatticeLambda()*parameters.getDeltaT(),
                                          parameters.getExponent(),
                                       kappa,parameters.getLy()),tol,maxIter) , 0.0,
                                        (parameters.getNy()-1)*maxIter);
                T uMaxTmp0 = ti0(parameters.getLy()/(T)2);
                if (uMaxTmp0 < uMax) {
                    kappa0 = kappa;
                    kappa = kappaTmp;
                }
            }
        }
    }

    pcout << "Compute pressure did not converge." << std::endl;
    exit(1);

    return kappa;
}

template <typename T, int nDim>
class CarreauVelocity
{
public :
    CarreauVelocity(CarreauFlowParam<T> parameters_, T kappa_, T tol_, plint maxIter_) :
        parameters(parameters_), kappa(kappa_), tol(tol_), maxIter(maxIter_)
    {   }

    void operator()(int iX, int iY, Array<T,nDim> &u)
    {
        TrapeziumIntegration<T> ti(new NewtonRaphson<T>(
                                   new CarreauFunction_2<T>(parameters.getLatticeNu0()/
                                        (parameters.getDeltaT()/(parameters.getDeltaX()*parameters.getDeltaX())),
                                        parameters.getLatticeLambda()*parameters.getDeltaT(),
                                        parameters.getExponent(),
                                        kappa,parameters.getLy()),tol,maxIter) , 0.0, iY *maxIter);

        u[0] = ti((T)iY/(T)(parameters.getNy()-1))*parameters.getLatticeU();
        u[1] = T();
    }
private :
    CarreauFlowParam<T> parameters;
    T kappa, tol;
    plint maxIter;
};

template<typename T, int nDim>
struct IniCarreauVelocityProcessor2D : public BoxProcessingFunctional2D_T<T,nDim>
{
    IniCarreauVelocityProcessor2D(CarreauFlowParam<T> parameters_,T kappa_, T tol_, plint maxIter_);
    
    virtual void process(Box2D domain,TensorField2D<T,nDim>& field);
    virtual IniCarreauVelocityProcessor2D<T,nDim>* clone() const;
    
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    
private :
    CarreauFlowParam<T> parameters;
    
    T kappa, tol;
    plint maxIter;
};

template<typename T, int nDim>
IniCarreauVelocityProcessor2D<T,nDim>::IniCarreauVelocityProcessor2D(
        CarreauFlowParam<T> parameters_,T kappa_, T tol_, plint maxIter_) :
            parameters(parameters_), kappa(kappa_), tol(tol_), maxIter(maxIter_)
{ }

template<typename T, int nDim>
void IniCarreauVelocityProcessor2D<T,nDim>::
        process(Box2D domain, TensorField2D<T,nDim>& field)
{
    enum {x,y}; //spatial coordinates

    Dot2D absoluteOffset = field.getLocation();

    std::vector<T> xvel, yvel;
    for (plint iY=domain.y0; iY<=domain.y1; ++iY)
    {
        Array<T,nDim> velocity;
        plint realY = iY + absoluteOffset.y;

        CarreauVelocity<T,nDim> cv(parameters,kappa,tol,realY*maxIter);

        cv(0,realY,velocity);

        xvel.push_back(velocity[x]);
        yvel.push_back(velocity[y]);
    }


    for (plint iX=domain.x0; iX<=domain.x1; ++iX)
    {
        unsigned iV = 0;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY)
        {
            field.get(iX,iY)[x] = xvel[iV];
            field.get(iX,iY)[y] = yvel[iV];

            ++iV;
        }
    }
}

template<typename T, int nDim>
IniCarreauVelocityProcessor2D<T,nDim>*
    IniCarreauVelocityProcessor2D<T,nDim>::clone() const
{
    return new IniCarreauVelocityProcessor2D<T,nDim>(*this);
}

template<typename T, int nDim>
void IniCarreauVelocityProcessor2D<T,nDim>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor, int nDim>
struct IniCarreauProcessor2D : public BoxProcessingFunctional2D_LT<T,Descriptor,T,nDim>
{
    IniCarreauProcessor2D(CarreauFlowParam<T> parameters_,T kappa_, T tol_, plint maxIter_);

    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
                         TensorField2D<T,nDim>& field);
    virtual IniCarreauProcessor2D<T,Descriptor,nDim>* clone() const;

    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;

private :
    CarreauFlowParam<T> parameters;

    T kappa, tol;
    plint maxIter;
};

template<typename T, template<typename U> class Descriptor, int nDim>
IniCarreauProcessor2D<T,Descriptor,nDim>::IniCarreauProcessor2D(
        CarreauFlowParam<T> parameters_,T kappa_, T tol_, plint maxIter_) :
            parameters(parameters_), kappa(kappa_), tol(tol_), maxIter(maxIter_)
{  }

template<typename T, template<typename U> class Descriptor, int nDim>
void IniCarreauProcessor2D<T,Descriptor,nDim>::
        process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice,TensorField2D<T,nDim>& field)
{
    enum {x,y}; //spatial coordinates

    Dot2D relativeOffset = computeRelativeDisplacement(lattice, field);

    CarreauVelocity<T,nDim> cv(parameters,kappa,tol,maxIter);

    for (plint iX=domain.x0; iX<=domain.x1; ++iX)
    {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY)
        {
            plint relativeX = iX + relativeOffset.x;
            plint relativeY = iY + relativeOffset.y;

            T rho = (T)1;

            lattice.get(iX,iY).defineDensity(rho);
            lattice.get(iX,iY).defineVelocity(field.get(relativeX,relativeY));

            iniCellAtEquilibrium(lattice.get(iX,iY), rho, field.get(relativeX,relativeY));
        }
    }
}

template<typename T, template<typename U> class Descriptor, int nDim>
IniCarreauProcessor2D<T,Descriptor,nDim>*
    IniCarreauProcessor2D<T,Descriptor,nDim>::clone() const
{
    return new IniCarreauProcessor2D<T,Descriptor,nDim>(*this);
}

template<typename T, template<typename U> class Descriptor, int nDim>
void IniCarreauProcessor2D<T,Descriptor,nDim>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

void definePoiseuilleGeometry( MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
                               MultiTensorField2D<T,2> &velField,
                               CarreauFlowParam<T> const& parameters,
                               OnLatticeBoundaryCondition2D<T,DESCRIPTOR>& boundaryCondition,
                               T tol, plint maxIter)
{
    setCompositeDynamics (
            lattice,
            lattice.getBoundingBox(),
            new CarreauDynamics<T,DESCRIPTOR,1>(new NoDynamics<T,DESCRIPTOR>) );

    // Create Velocity boundary conditions
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice);

    T kappa = computePressureGradient(parameters,tol,200);

    pcout << "gradP = " << kappa << endl;

    applyProcessingFunctional(new IniCarreauVelocityProcessor2D<T,2>(parameters,kappa, tol, maxIter),
                            velField.getBoundingBox(),velField);

    applyProcessingFunctional(new IniCarreauProcessor2D<T,DESCRIPTOR,2>(parameters,kappa, tol, maxIter),
                            lattice.getBoundingBox(),lattice,velField);

    lattice.initialize();

    pcout << "lattice totally initialized." << endl;
}

T computeRMSerror ( MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
                    MultiTensorField2D<T,2> &analyticalVelocity,
                    CarreauFlowParam<T> const& parameters )
{
    MultiTensorField2D<T,2> numericalVelocity(lattice);
    computeVelocity(lattice, numericalVelocity, lattice.getBoundingBox());
    
    // Divide by lattice velocity to normalize the error
    return 1./parameters.getLatticeU() *
    // Compute RMS difference between analytical and numerical solution
    std::sqrt( computeAverage( *computeNormSqr(
        *subtract(analyticalVelocity, numericalVelocity)
    ) ) );
}

void writeVTK(MultiBlockLattice2D<T,DESCRIPTOR>& lattice,
            MultiTensorField2D<T,2> &field,
            CarreauFlowParam<T> const& parameters, plint iter)
{
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();
    VtkImageOutput2D<T> vtkOut(createFileName("vtk", iter, 6), dx);

    vtkOut.writeData<2,float>(*computeVelocity(lattice), "vel", dx/dt);
    vtkOut.writeData<2,float>(field, "ana_vel", dx/dt);
    vtkOut.writeData<float>(*computeOmega(lattice), "omega", (T)1);
}

int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");
    global::timer("simTime").start();

    if (argc != 2)
    {
        pcout << "Error N must be specified." << std::endl;
        exit(1);
    }

    const plint Nref = 50;
    const T uMaxRef = 0.01;
    const plint N = atoi(argv[1]);

    const T uMax = uMaxRef * (T)Nref / (T)N;

    CarreauFlowParam<T> parameters(
            (T) uMax,    // uMax
            (T) 1.e0,    // Re
            (T)10.0,     // Cu
            (T)0.0,      // NuInf (Only the case nuInf = 0 has been implemented for the velocity inlet bc)
            (T)0.5,      // n
             N,          // N
             1.,         // lx
             1.          // ly
    );
    
    pcout << "nu0 = " << parameters.getLatticeNu0() << ", nuInf = " << parameters.getLatticeNuInf() << std::endl;
    const plint maxT     = 10000000;

    writeLogFile(parameters, "Carreau Poseuille Flow");

    global::CarreauParameters().setNu0(parameters.getLatticeNu0());
    global::CarreauParameters().setNuInf(parameters.getLatticeNuInf());
    global::CarreauParameters().setLambda(parameters.getLatticeLambda());
    global::CarreauParameters().setExponent(parameters.getExponent());

    MultiBlockLattice2D<T, DESCRIPTOR> lattice (
            parameters.getNx(), parameters.getNy(),
            new BGKdynamics<T,DESCRIPTOR>(parameters.getOmega0()) );
    
    defineDynamics(lattice, lattice.getBoundingBox(), new BGKdynamics<T,DESCRIPTOR>(parameters.getOmega0()) );

    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>*
        boundaryCondition = createLocalBoundaryCondition2D<T,DESCRIPTOR>();

    T tolAll = 1.0e-11;
    plint maxIterAll = 100;

    MultiTensorField2D<T,2> velField(parameters.getNx(), parameters.getNy());

    definePoiseuilleGeometry(lattice, velField, parameters, *boundaryCondition, tolAll, maxIterAll);

    util::ValueTracer<T> converge(uMax,Nref,1.0e-4);

    T tIni = global::timer("simTime").stop();
    pcout << "time elapsed for iniGeometry:" << tIni << endl;
    global::timer("simTime").start();

    plint iT = 0;
    for (iT=0; iT<maxT; ++iT)
    {
        converge.takeValue(getStoredAverageEnergy(lattice),true);
        if (iT % 1000 == 0)
        {
            pcout << iT << " : Writing image." << endl;
            writeVTK(lattice,velField, parameters,iT);
        }
        
        if (converge.hasConverged())
        {
            pcout << "Simulation converged." << endl;
            break;
        }

        lattice.collideAndStream();
    }

    T tEnd = global::timer("simTime").stop();

    T totalTime = tEnd-tIni;
    T N1000 = lattice.getNx()/(T)1000;
    pcout << "N=" << N << endl;
    pcout << "number of processors: " << global::mpi().getSize() << endl;
    pcout << "simulation time: " << totalTime << endl;
    pcout << "total time: " << tEnd << endl;
    pcout << "total iterations: " << iT << endl;
    pcout << "Msus: " << N1000*N1000*(T)iT/totalTime << endl;
    
    plb_ofstream fout("tmp/vel.dat");
    fout << *computeVelocity(lattice);
    fout.close();
    fout.open("tmp/ana_vel.dat");
    fout << velField;
    fout.close();

    pcout << "For N = " << N << ", Error = " << computeRMSerror(lattice, velField,parameters ) << endl;

    delete boundaryCondition;
}
