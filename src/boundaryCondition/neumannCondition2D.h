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
 * Neumann boundary conditions -- header file.
 */
#ifndef NEUMANN_CONDITION_2D_H
#define NEUMANN_CONDITION_2D_H

#include "core/globalDefs.h"
#include "core/blockLatticeBase2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataProcessorWrapper2D.h"
#include "boundaryCondition/boundaryDynamics.h"

namespace plb {

template<typename T, template<typename U> class Descriptor, int direction, int orientation>
class CopyUnknownPopulationsFunctional2D : public BoxProcessingFunctional2D_L<T,Descriptor>
{
public:
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice);
    virtual CopyUnknownPopulationsFunctional2D<T,Descriptor,direction,orientation>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor, int normalX, int normalY>
class CopyAllPopulationsFunctional2D : public BoxProcessingFunctional2D_L<T,Descriptor>
{
public:
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice);
    virtual CopyAllPopulationsFunctional2D<T,Descriptor,normalX,normalY>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor, int normalX, int normalY>
class CopyVelocityFunctional2D : public BoxProcessingFunctional2D_L<T,Descriptor>
{
public:
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice);
    virtual CopyVelocityFunctional2D<T,Descriptor,normalX,normalY>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor, int normalX, int normalY>
class CopyTangentialVelocityFunctional2D : public BoxProcessingFunctional2D_L<T,Descriptor>
{
public:
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice);
    virtual CopyTangentialVelocityFunctional2D<T,Descriptor,normalX,normalY>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor, int normalX, int normalY>
class CopyNormalVelocityFunctional2D : public BoxProcessingFunctional2D_L<T,Descriptor>
{
public:
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice);
    virtual CopyNormalVelocityFunctional2D<T,Descriptor,normalX,normalY>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor, int normalX, int normalY>
class CopyDensityFunctional2D : public BoxProcessingFunctional2D_L<T,Descriptor>
{
public:
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice);
    virtual CopyDensityFunctional2D<T,Descriptor,normalX,normalY>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

}  // namespace plb

#endif  // NEUMANN_CONDITION_2D_H
