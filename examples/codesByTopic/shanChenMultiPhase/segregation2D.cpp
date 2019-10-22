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

using namespace plb;
using namespace std;

typedef double T;

#define DESCRIPTOR descriptors::ShanChenD2Q9Descriptor

template<typename T, template<typename U> class Descriptor>
class RandomInitializer : public BoxProcessingFunctional2D_L<T,Descriptor> 
{
public :
    RandomInitializer(T rho0_, T maxRho_) : rho0(rho0_), maxRho(maxRho_)
    { };
    virtual void process(Box2D domain,BlockLattice2D<T,Descriptor>& lattice)
    {
        for (plint iX = domain.x0; iX <= domain.x1; ++iX) {
            for (plint iY = domain.y0; iY <= domain.y1; ++iY)
            {
                T rho = rho0 + ((T)rand()/(T)RAND_MAX)*maxRho;
                Array<T,2> zeroVelocity (0.,0.);
                
                iniCellAtEquilibrium(lattice.get(iX,iY), rho, zeroVelocity);
                
                lattice.get(iX,iY).setExternalField (
                        Descriptor<T>::ExternalField::densityBeginsAt,
                        Descriptor<T>::ExternalField::sizeOfDensity, &rho );
                lattice.get(iX,iY).setExternalField (
                        Descriptor<T>::ExternalField::momentumBeginsAt,
                        Descriptor<T>::ExternalField::sizeOfMomentum, &zeroVelocity[0] );
            }
        }
    };
    virtual RandomInitializer<T,Descriptor>* clone() const
    {
        return new RandomInitializer<T,Descriptor>(*this);
    };
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const
    {
        modified[0] = modif::staticVariables;
    };
    
private :
    T rho0, maxRho;
};


int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");

    srand(global::mpi().getRank() + 3);

    // For the choice of the parameters G, rho0, and psi0, we refer to the book
    //   Michael C. Sukop and Daniel T. Thorne (2006), 
    //   Lattice Boltzmann Modeling; an Introduction for Geoscientists and Engineers.
    //   Springer-Verlag Berlin/Heidelberg.
    
    const T omega = 1.0;
    const int nx   = 400;
    const int ny   = 400;
    const T G      = -120.0;
    const int maxIter  = 100001;
    const int saveIter = 100;
    const int statIter = 100;
    
    const T rho0 = 200.0;
    const T deltaRho = 1.0;
    const T psi0 = 4.0;

    MultiBlockLattice2D<T, DESCRIPTOR> lattice (
            nx,ny, new ExternalMomentBGKdynamics<T, DESCRIPTOR>(omega) );
            
    lattice.periodicity().toggleAll(true);

    // Use a random initial condition, to activate the phase separation.
    applyProcessingFunctional(new RandomInitializer<T,DESCRIPTOR>(rho0,deltaRho), 
                              lattice.getBoundingBox(),lattice);

    // Add the data processor which implements the Shan/Chen interaction potential.
    plint processorLevel = 1;
    integrateProcessingFunctional (
            new ShanChenSingleComponentProcessor2D<T,DESCRIPTOR> (
                G, new interparticlePotential::PsiShanChen94<T>(psi0,rho0) ),
            lattice.getBoundingBox(), lattice, processorLevel );

    lattice.initialize();
    
    pcout << "Starting simulation" << endl;
    for (int iT=0; iT<maxIter; ++iT) {
        if (iT%statIter==0) {
            auto_ptr<MultiScalarField2D<T> > rho( computeDensity(lattice) );
            pcout << iT << ": Average rho fluid one = " << computeAverage(*rho) << endl;
            pcout << "Minimum density: " << computeMin(*rho) << endl;
            pcout << "Maximum density: " << computeMax(*rho) << endl;
        }
        if (iT%saveIter == 0) {
            ImageWriter<T>("leeloo").writeScaledGif (
                    createFileName("rho", iT, 6), *computeDensity(lattice) );
        }

        lattice.collideAndStream();
    }
}

