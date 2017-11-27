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

/* Main author: Daniel Lagrava
 */

/** \file
 * Dynamics and data processors used to implement 3D grid refinement -- header file.
 */

#ifndef COARSE_GRID_PROCESSORS_3D_H
#define COARSE_GRID_PROCESSORS_3D_H

#include "core/globalDefs.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataProcessorWrapper3D.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiGrid/gridRefinement.h"
#include "multiGrid/gridRefinementDynamics.h"
#include <vector>

namespace plb {

/// Coupling to be added to fine lattice: copies data to coarse lattice after numTimeSteps interations
template<typename T, template<typename U> class Descriptor1, template<typename U> class Descriptor2>
class CopyFineToCoarse3D : public BoxProcessingFunctional3D_LL<T,Descriptor1,T,Descriptor2>
{
public:
    CopyFineToCoarse3D( RescaleEngine<T,Descriptor1>* rescaleEngine_,
                        plint numTimeSteps_, plint executionTime_ );
    virtual ~CopyFineToCoarse3D();
    CopyFineToCoarse3D(CopyFineToCoarse3D<T,Descriptor1,Descriptor2> const& rhs);
    CopyFineToCoarse3D<T,Descriptor1,Descriptor2>& operator=(CopyFineToCoarse3D<T,Descriptor1,Descriptor2> const& rhs);
    virtual void process( Box3D domain,
                          BlockLattice3D<T,Descriptor1>& fineLattice,
                          BlockLattice3D<T,Descriptor2>& coarseLattice );
    virtual CopyFineToCoarse3D<T,Descriptor1,Descriptor2>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulk;
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
        modified[1] = modif::staticVariables;
    }
private:
    RescaleEngine<T,Descriptor1>* rescaleEngine;
    plint numTimeSteps;
    plint executionTime;
};


/// Coupling to be added to fine lattice: copies data from fine to coarse and does a filtering operation
template<typename T, template<typename U> class Descriptor1, template<typename U> class Descriptor2>
class CopyFineToCoarseWithFiltering3D : public BoxProcessingFunctional3D_LL<T,Descriptor1,T,Descriptor2>
{
public:
    CopyFineToCoarseWithFiltering3D( RescaleEngine<T,Descriptor1>* rescaleEngine_,
                        plint numTimeSteps_, plint executionTime_, std::vector<plint> indices_  );
    virtual ~CopyFineToCoarseWithFiltering3D();
    CopyFineToCoarseWithFiltering3D(CopyFineToCoarseWithFiltering3D<T,Descriptor1,Descriptor2> const& rhs);
    CopyFineToCoarseWithFiltering3D<T,Descriptor1,Descriptor2>&
operator=(CopyFineToCoarseWithFiltering3D<T,Descriptor1,Descriptor2> const& rhs);
    virtual void process( Box3D domain,
                          BlockLattice3D<T,Descriptor1>& fineLattice,
                          BlockLattice3D<T,Descriptor2>& coarseLattice );
    virtual CopyFineToCoarseWithFiltering3D<T,Descriptor1,Descriptor2>* clone() const;
    virtual BlockDomain::DomainT appliesTo() const {
        return BlockDomain::bulk;
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
        modified[1] = modif::staticVariables;
    }
private:
    RescaleEngine<T,Descriptor1>* rescaleEngine;
    plint numTimeSteps;
    plint executionTime;
    
    // a set of indices to perform the filtering
    std::vector<plint> indices; 
    
};


}  // namespace plb

#endif  // COARSE_GRID_PROCESSORS_3D_H
