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

#ifndef ADVECTION_DIFFUSION_PROCESSOR_3D_HH
#define ADVECTION_DIFFUSION_PROCESSOR_3D_HH

#include "complexDynamics/advectionDiffusionProcessor3D.h"

namespace plb {

// This data processor uses symmetric finite differences to compute a gradient.
// It cannot be applied on any part of the boundary of the global simulation
// domain, except if this boundary is periodic.
template<typename T, template<typename U> class Descriptor>
void SetEffectiveDiffusivity3D<T,Descriptor>::process (
        Box3D domain,
        BlockLattice3D<T,Descriptor>& lattice,
        ScalarField3D<T>& rhoBarField )
{
    Dot3D offset = computeRelativeDisplacement(lattice, rhoBarField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        plint iX2 = iX+offset.x;
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            plint iY2 = iY+offset.y;
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iZ2 = iZ+offset.z;
                T gradX = 0.5* (rhoBarField.get(iX2+1, iY2, iZ2) - rhoBarField.get(iX2-1, iY2, iZ2));
                T gradY = 0.5* (rhoBarField.get(iX2, iY2+1, iZ2) - rhoBarField.get(iX2, iY2-1, iZ2));
                T gradZ = 0.5* (rhoBarField.get(iX2, iY2, iZ2+1) - rhoBarField.get(iX2, iY2, iZ2-1));
                T normGradT = std::sqrt(gradX*gradX + gradY*gradY + gradZ*gradZ);

                // Model: d = d0 ( 1 + C^2 h |gradT|/T0 )
                T omega = omega0 / (1. + (1.-0.5*omega0)*cSmagoSqr*invT0*normGradT);

                lattice.get(iX,iY,iZ).getDynamics().setOmega(omega);
            }
        }
    }
}

} // namespace plb

#endif  // ADVECTION_DIFFUSION_PROCESSOR_3D_HH

