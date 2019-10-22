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

#ifndef POISEUILLE_HH

#include "poiseuille.h"

/// Velocity on the parabolic Poiseuille profile
template<typename T>
T poiseuilleVelocity(int iY, plb::IncomprFlowParam<T> const& parameters) {
    T y = (T)iY / parameters.getResolution();
    return 4.*parameters.getLatticeU() * (y-y*y);
}

/// Linearly decreasing pressure profile
template<typename T>
T poiseuillePressure(int iX, plb::IncomprFlowParam<T> const& parameters) {
    T Lx = parameters.getNx()-1;
    T Ly = parameters.getNy()-1;
    return 8.*parameters.getLatticeNu()*parameters.getLatticeU() / (Ly*Ly) * (Lx/(T)2-(T)iX);
}

/// Convert pressure to density according to ideal gas law
template<typename T, template<typename U> class Descriptor>
T poiseuilleDensity(int iX, plb::IncomprFlowParam<T> const& parameters) {
    return poiseuillePressure(iX,parameters)*Descriptor<T>::invCs2 + (T)1;
}

/// A functional, used to initialize the velocity for the boundary conditions
template<typename T>
class PoiseuilleVelocity {
public:
    PoiseuilleVelocity(plb::IncomprFlowParam<T> parameters_)
        : parameters(parameters_)
    { }
    void operator()(int iX, int iY, plb::Array<T,2>& u) const {
        u[0] = poiseuilleVelocity(iY, parameters);
        u[1] = T();
    }
private:
    plb::IncomprFlowParam<T> parameters;
};

/// A functional, used to initialize the density for the boundary conditions
template<typename T, template<typename U> class Descriptor>
class PoiseuilleDensity {
public:
    PoiseuilleDensity(plb::IncomprFlowParam<T> parameters_)
        : parameters(parameters_)
    { }
    T operator()(int iX, int iY) const {
        return poiseuilleDensity<T,Descriptor>(iY, parameters);
    }
private:
    plb::IncomprFlowParam<T> parameters;
};

/// A functional, used to create an initial condition for the density and velocity
template<typename T, template<typename U> class Descriptor>
class PoiseuilleVelocityAndDensity {
public:
    PoiseuilleVelocityAndDensity(plb::IncomprFlowParam<T> parameters_)
        : parameters(parameters_)
    { }
    void operator()(int iX, int iY, T& rho, plb::Array<T,2>& u) const {
        rho = poiseuilleDensity<T,Descriptor>(iX,parameters);
        u[0] = poiseuilleVelocity<T>(iY, parameters);
        u[1] = T();
    }
private:
    plb::IncomprFlowParam<T> parameters;
};

template<typename T, template<typename U> class Descriptor>
void createPoiseuilleBoundaries( plb::MultiBlockLattice2D<T,Descriptor>& lattice,
                                 plb::IncomprFlowParam<T> const& parameters,
                                 plb::OnLatticeBoundaryCondition2D<T,Descriptor>& boundaryCondition )
{
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice);

    setBoundaryVelocity (
            lattice, lattice.getBoundingBox(),
            PoiseuilleVelocity<T>(parameters) );
}

template<typename T, template<typename U> class Descriptor>
void createPoiseuilleInitialValues( plb::MultiBlockLattice2D<T,Descriptor>& lattice,
                                    plb::IncomprFlowParam<T> const& parameters )
{
    initializeAtEquilibrium (
            lattice, lattice.getBoundingBox(),
            PoiseuilleVelocityAndDensity<T,Descriptor>(parameters) );
}

#endif  // POISEUILLE_HH
