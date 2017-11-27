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

/* Main author: Orestis Malaspinas.
 */

/** \file
 * Neumann boundary conditions for temperature -- header file.
 */
#ifndef ADIABATIC_BOUNDARY_PROCESSOR_2D_H
#define ADIABATIC_BOUNDARY_PROCESSOR_2D_H

#include "core/globalDefs.h"
#include "core/blockLatticeBase2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataProcessorWrapper2D.h"

namespace plb {

template<typename T, template<typename U> class Descriptor, int direction, int orientation> 
class FlatAdiabaticBoundaryFunctional2D : public BoxProcessingFunctional2D_L<T,Descriptor>
{
public:
    virtual plint extent() const { return 2; }
    virtual plint extent(int whichDirection) const { return 2; }
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice);
    virtual FlatAdiabaticBoundaryFunctional2D<T,Descriptor,direction,orientation>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

}  // namespace plb

#endif  // ADIABATIC_BOUNDARY_PROCESSOR_2D_H
