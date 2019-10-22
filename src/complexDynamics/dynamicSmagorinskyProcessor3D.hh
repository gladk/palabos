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

#ifndef DYNAMIC_SMAGORINSKY_PROCESSOR_3D_HH
#define DYNAMIC_SMAGORINSKY_PROCESSOR_3D_HH

#include "complexDynamics/dynamicSmagorinskyProcessor3D.h"
#include "latticeBoltzmann/momentTemplates.h"

namespace plb {

template<typename T, template<typename U> class Descriptor> 
ComputeSmagoViscosityFunctional3D<T,Descriptor>::ComputeSmagoViscosityFunctional3D(T omega0_, T cSmago_)
    : omega0(omega0_),
      preFactor(computePreFactor(omega0_, cSmago_))
{ }

template<typename T, template<typename U> class Descriptor> 
void ComputeSmagoViscosityFunctional3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Cell<T,Descriptor>& cell = lattice.get(iX,iY,iZ);
                T rhoBar;
                Array<T,Descriptor<T>::d> j;
                Array<T,SymmetricTensor<T,Descriptor>::n> PiNeq;
                momentTemplates<T,Descriptor>::compute_rhoBar_j_PiNeq(cell, rhoBar, j, PiNeq);
                T PiNeqNormSqr = SymmetricTensor<T,Descriptor>::tensorNormSqr(PiNeq);
                T PiNeqNorm    = std::sqrt(PiNeqNormSqr);
                T alpha        = preFactor * Descriptor<T>::invRho(rhoBar);
                T linearTerm   = alpha*PiNeqNorm;
                T squareTerm   = (T)2*alpha*alpha*PiNeqNormSqr;
                // In the following formula, the square-root appearing in the explicit form of omega
                //   is developed to second-order.
                T omega = omega0*(1-linearTerm+squareTerm);
                *cell.getExternal(Descriptor<T>::ExternalField::omegaBeginsAt) = omega;
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor> 
ComputeSmagoViscosityFunctional3D<T,Descriptor>* ComputeSmagoViscosityFunctional3D<T,Descriptor>::clone() const
{
    return new ComputeSmagoViscosityFunctional3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor> 
void ComputeSmagoViscosityFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor> 
BlockDomain::DomainT ComputeSmagoViscosityFunctional3D<T,Descriptor>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

template<typename T, template<typename U> class Descriptor> 
T ComputeSmagoViscosityFunctional3D<T,Descriptor>::computePreFactor(T omega0, T cSmago)
{
    return (T)0.5 * util::sqr(cSmago*omega0*Descriptor<T>::invCs2);
}

} // namespace plb

#endif  // DYNAMIC_SMAGORINSKY_PROCESSOR_3D_HH

