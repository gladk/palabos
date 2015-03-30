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

#ifndef FREE_SURFACE_BOUNDARY_CONDITION_3D_HH
#define FREE_SURFACE_BOUNDARY_CONDITION_3D_HH

#include "multiPhysics/freeSurfaceBoundaryCondition3D.h"
#include <cmath>
#include <iostream>

namespace plb {

template<typename T, template<typename U> class Descriptor>
FreeSurfaceFadingArea3D<T,Descriptor>::FreeSurfaceFadingArea3D(T factor_)
    : factor(factor_)
{ }

template<typename T, template<typename U> class Descriptor>
void FreeSurfaceFadingArea3D<T,Descriptor>::process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice)
{
    std::vector<T> decomposedVariables;

    enum {
        forceOffset          = Descriptor<T>::ExternalField::forceBeginsAt,
        momentumStoredOffset = Descriptor<T>::ExternalField::momentumBeginsAt,
        densityStoredOffset  = Descriptor<T>::ExternalField::densityBeginsAt,
    };
    
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Cell<T,Descriptor>& cell = lattice.get(iX,iY,iZ);
                plint order = 0;
                cell.getDynamics().decompose(cell, decomposedVariables, order);
                
                T density = Descriptor<T>::fullRho(decomposedVariables[0]);
                if (density > T(0)) density *= factor;
                decomposedVariables[0] = Descriptor<T>::rhoBar(density);
                cell.getDynamics().recompose(cell, decomposedVariables, order);

                *cell.getExternal(densityStoredOffset) = density;                    
                
                Array<T,Descriptor<T>::d> j;
                j.resetToZero();
                T rhoBar; 
                momentTemplates<T,Descriptor>::get_rhoBar_j(cell, rhoBar, j);
                  
                // TODO: What about mass, volumeFraction, flagStatus?
                j.to_cArray(cell.getExternal(momentumStoredOffset));
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
FreeSurfaceFadingArea3D<T,Descriptor>* FreeSurfaceFadingArea3D<T,Descriptor>::clone() const {
    return new FreeSurfaceFadingArea3D<T,Descriptor>(*this);
}


template<typename T, template<typename U> class Descriptor>
void RemoveMass3D<T,Descriptor>::processGenericBlocks(Box3D domain,std::vector<AtomicBlock3D*> atomicBlocks)
{
    using namespace twoPhaseFlag;
    FreeSurfaceProcessorParam3D<T,Descriptor> param(atomicBlocks);
            
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                //param.attributeDynamics(iX,iY,iZ, new NoDynamics<T,Descriptor>((T)1.));
                param.setDensity(iX,iY,iZ, (T)1.);
                param.setMomentum(iX,iY,iZ, Array<T,3>((T)0.,(T)0.,(T)0.));
                param.mass(iX,iY,iZ) = (T)0;
                param.volumeFraction(iX,iY,iZ) = (T)0;
                //param.flag(iX,iY,iZ) = empty;
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
RemoveMass3D<T,Descriptor>* RemoveMass3D<T,Descriptor>::clone() const {
    return new RemoveMass3D<T,Descriptor>(*this);
}


template<typename T, template<typename U> class Descriptor>
void PouringLiquid3D<T,Descriptor>::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks)
{
    using namespace twoPhaseFlag;
    FreeSurfaceProcessorParam3D<T,Descriptor> param(atomicBlocks);
    
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                T iniRho = T(1);
                param.attributeDynamics (
                        iX,iY,iZ, dynamicsTemplate->clone() );
                iniCellAtEquilibrium(param.cell(iX,iY,iZ), iniRho, injectionVelocity);
                param.setDensity(iX,iY,iZ, iniRho);
                param.setMomentum(iX,iY,iZ, iniRho*injectionVelocity);
                param.mass(iX,iY,iZ) = iniRho;
                param.volumeFraction(iX,iY,iZ) = (T)1;
                param.flag(iX,iY,iZ) = fluid;
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
PouringLiquid3D<T,Descriptor>* PouringLiquid3D<T,Descriptor>::clone() const {
    return new PouringLiquid3D<T,Descriptor>(*this);
}


template<typename T, template<typename U> class Descriptor>
void ShortenBounceBack3D<T,Descriptor>::processGenericBlocks(Box3D domain,std::vector<AtomicBlock3D*> atomicBlocks)
{
    using namespace twoPhaseFlag;
    typedef Descriptor<T> D;
    FreeSurfaceProcessorParam3D<T,Descriptor> param(atomicBlocks);
    Box3D extDomain = domain.enlarge(1);

    for (plint iX=extDomain.x0; iX<=extDomain.x1; ++iX) {
        bool xBoundary = iX==extDomain.x0 || iX==extDomain.x1;
        for (plint iY=extDomain.y0; iY<=extDomain.y1; ++iY) {
            bool yBoundary = xBoundary || iY==extDomain.y0 || iY==extDomain.y1;
            for (plint iZ=extDomain.z0; iZ<=extDomain.z1; ++iZ) {
                if (param.flag(iX,iY,iZ)==wall) {
                    bool zBoundary = yBoundary || iZ==extDomain.z0 || iZ==extDomain.z1;
                    for (plint iNeighbor=1; iNeighbor<D::q; ++iNeighbor) {
                        plint nextX = iX+D::c[iNeighbor][0];
                        plint nextY = iY+D::c[iNeighbor][1];
                        plint nextZ = iZ+D::c[iNeighbor][2];
                        if(!zBoundary || contained(nextX,nextY,nextZ, domain)) {
                            if (isWet(param.flag(nextX,nextY,nextZ))) {
                                plint opp = indexTemplates::opposite<D>(iNeighbor);
                                param.cell(nextX,nextY,nextZ)[iNeighbor] = param.cell(iX,iY,iZ)[opp];
                            }
                        }
                    }
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
ShortenBounceBack3D<T,Descriptor>* ShortenBounceBack3D<T,Descriptor>::clone() const {
    return new ShortenBounceBack3D<T,Descriptor>(*this);
}

}  // namespace plb

#endif  // FREE_SURFACE_BOUNDARY_CONDITION_3D_HH

