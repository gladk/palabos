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
 * Helper functions for domain initialization -- header file.
 */
#ifndef BC_WRAPPER_2D_H
#define BC_WRAPPER_2D_H

#include "boundaryCondition/boundaryCondition2D.h"

namespace plb {

template<typename T, template<typename U> class Descriptor>
class OuterBoxBC {
public:
    OuterBoxBC(OnLatticeBoundaryCondition2D<T,Descriptor>* bc_);
    OuterBoxBC(OuterBoxBC<T,Descriptor> const& rhs);
    ~OuterBoxBC();
    OuterBoxBC<T,Descriptor>& operator=(OuterBoxBC<T,Descriptor> const& rhs);
    OuterBoxBC<T,Descriptor>* clone() const;
    void setVelocityCondition( MultiBlockLattice2D<T,Descriptor>& lattice,
                               Box2D domain );
    void setPressureCondition( MultiBlockLattice2D<T,Descriptor>& lattice,
                               Box2D domain );
private:
    OnLatticeBoundaryCondition2D<T,Descriptor>* bc;
};

template<typename T, template<typename U> class Descriptor>
OuterBoxBC<T,Descriptor>* generateRegularizedBC();

}  // namespace plb

#endif  // BC_WRAPPER_2D_H
