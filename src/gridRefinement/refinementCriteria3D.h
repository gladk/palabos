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

#ifndef REFINEMENT_CRITERIA_3D_H
#define REFINEMENT_CRITERIA_3D_H

#include "core/globalDefs.h"
#include "atomicBlock/dataProcessingFunctional3D.h"

namespace plb {

template<typename T, template<typename U> class Descriptor>
class ComputeRefinementRvalueFunctional3D : public BoxProcessingFunctional3D_LS<T, Descriptor, T> {
public:
    ComputeRefinementRvalueFunctional3D(T knudsen_);
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice, ScalarField3D<T>& scalarField);
    virtual ComputeRefinementRvalueFunctional3D* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    T knudsen;
};

template<typename T, template<typename U> class Descriptor>
std::auto_ptr<MultiScalarField3D<T> > computeRvalues(MultiBlockLattice3D<T,Descriptor>& lattice, T knudsen, Box3D const& domain);

} // namespace plb

#endif  // REFINEMENT_CRITERIA_3D_H
