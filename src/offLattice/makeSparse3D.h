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

#ifndef MAKE_SPARSE_3D_H
#define MAKE_SPARSE_3D_H

#include "core/globalDefs.h"
#include "atomicBlock/dataProcessingFunctional3D.h"
#include "multiBlock/multiContainerBlock3D.h"
#include "multiBlock/multiDataField3D.h"

namespace plb {

template<typename T>
class ComputeSparsityFunctional3D :
          public PlainReductiveBoxProcessingFunctional3D
{
public:
    ComputeSparsityFunctional3D();
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual ComputeSparsityFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
    pluint getNumBlocks() const;
private:
    plint numBlocksId;
};
    
template<typename T>
MultiBlockManagement3D computeSparseManagement (
        MultiScalarField3D<T>& field, plint newEnvelopeWidth );

}  // namespace plb

#endif  // MAKE_SPARSE_3D_H

