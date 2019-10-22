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

#ifndef MAKE_SPARSE_2D_H
#define MAKE_SPARSE_2D_H

#include "core/globalDefs.h"
#include "atomicBlock/dataProcessingFunctional2D.h"
#include "multiBlock/multiContainerBlock2D.h"
#include "multiBlock/multiDataField2D.h"

namespace plb {

template<typename T>
class ComputeSparsityFunctional2D :
          public PlainReductiveBoxProcessingFunctional2D
{
public:
    ComputeSparsityFunctional2D();
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> fields);
    virtual ComputeSparsityFunctional2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
    pluint getNumBlocks() const;
private:
    plint numBlocksId;
};
    
template<typename T>
MultiBlockManagement2D computeSparseManagement (
        MultiScalarField2D<T>& field, plint newEnvelopeWidth );

}  // namespace plb

#endif  // MAKE_SPARSE_2D_H

