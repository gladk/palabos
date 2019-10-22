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

#ifndef WRAPPED_LOCAL_BOUNDARY_PROCESSOR_2D_H
#define WRAPPED_LOCAL_BOUNDARY_PROCESSOR_2D_H

#include "core/globalDefs.h"
#include "atomicBlock/dataProcessingFunctional2D.h"
#include "atomicBlock/blockLattice2D.h"

namespace plb {

/**
* This class wraps the dynamics of a local boundary condition to
* present it in terms of a data processor. The advantage is that it
* "repairs" the unknown populations on a wall right after propagation,
* in case they should be used by another data processor.
*/
template<typename T, template<typename U> class Descriptor>
class WrappedLocalBoundaryFunctional2D : public BoxProcessingFunctional2D_L<T,Descriptor>
{
public:
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice);
    virtual WrappedLocalBoundaryFunctional2D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
    }
    virtual int getStaticId() const { return staticId; }
    static const int staticId;
};

}  // namespace plb

#endif  // WRAPPED_LOCAL_BOUNDARY_PROCESSOR_2D_H
