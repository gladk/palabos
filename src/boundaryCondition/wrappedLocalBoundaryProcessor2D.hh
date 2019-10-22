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

#ifndef WRAPPED_LOCAL_BOUNDARY_PROCESSOR_2D_HH
#define WRAPPED_LOCAL_BOUNDARY_PROCESSOR_2D_HH

#include "boundaryCondition/wrappedLocalBoundaryProcessor2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "core/processorIdentifiers2D.h"
#include <typeinfo>

namespace plb {

///////////  WrappedLocalBoundaryFunctional2D ///////////////////////////////////

template<typename T, template<typename U> class Descriptor>
const int WrappedLocalBoundaryFunctional2D<T,Descriptor>::staticId =
    meta::registerProcessor2D < WrappedLocalBoundaryFunctional2D<T,Descriptor>, T, Descriptor> (std::string("WrappedLocalBoundary2D"));

template<typename T, template<typename U> class Descriptor>
void WrappedLocalBoundaryFunctional2D<T,Descriptor>::process (
        Box2D domain, BlockLattice2D<T,Descriptor>& lattice )
{
    PLB_PRECONDITION(domain.x0==domain.x1 || domain.y0==domain.y1);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            Cell<T,Descriptor>& cell = lattice.get(iX,iY);
            Dynamics<T,Descriptor>& dynamics = cell.getDynamics();
            // Make a safety check to avoid errors in case the user replaces
            //   the original composite dynamics by something else; for
            //   example, when bounce-back is used as an "eraser tool" at
            //   inlets/outlets.
            if (dynamics.isComposite()) {
                dynamic_cast<CompositeDynamics<T,Descriptor>&>(dynamics)
                    .prepareCollision(cell);
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
WrappedLocalBoundaryFunctional2D<T,Descriptor>*
    WrappedLocalBoundaryFunctional2D<T,Descriptor>::clone() const
{
    return new WrappedLocalBoundaryFunctional2D<T,Descriptor>(*this);
}

}  // namespace plb

#endif  // WRAPPED_LOCAL_BOUNDARY_PROCESSOR_2D_HH
