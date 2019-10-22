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

/* The original version of this file was written by Orestis Malaspinas
 * and Andrea Parmigiani.
 */

#ifndef SHAN_CHEN_PROCESSOR_3D_H
#define SHAN_CHEN_PROCESSOR_3D_H

#include "core/globalDefs.h"
#include "atomicBlock/dataProcessorWrapper3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "multiPhysics/interparticlePotential.h"

namespace plb {

/// Shan-Chen coupling for multi-component flow with or without external force
template<typename T, template<typename U> class Descriptor>
class ShanChenMultiComponentProcessor3D :
    public LatticeBoxProcessingFunctional3D<T,Descriptor>
{
public:
    /// With these constructors, space- and time-dependent values of the
    ///   relaxation parameters omega are accounted for.
    ShanChenMultiComponentProcessor3D(T G_);
    ShanChenMultiComponentProcessor3D(std::vector<std::vector<T> > const& speciesG_);
    /// With these constructors, the values of the relaxation parameters omega are
    ///   taken to be species-dependent, but not space- or time-dependent. Their
    ///   value is imposed in the constructor.
    ShanChenMultiComponentProcessor3D(T G_, std::vector<T> const& imposedOmega_);
    ShanChenMultiComponentProcessor3D(std::vector<std::vector<T> > const& speciesG_, std::vector<T> const& imposedOmega_);
    virtual void process(Box3D domain, std::vector<BlockLattice3D<T,Descriptor>*> lattices );
    virtual ShanChenMultiComponentProcessor3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    T G;
    std::vector<T> speciesG;
    std::vector<T> imposedOmega;
};

/// Shan-Chen coupling for single-component flow with or without external force
template<typename T, template<typename U> class Descriptor>
class ShanChenSingleComponentProcessor3D : public BoxProcessingFunctional3D_L<T,Descriptor> {
public:
    ShanChenSingleComponentProcessor3D(T G_, interparticlePotential::PsiFunction<T>* Psi_);
    virtual ~ShanChenSingleComponentProcessor3D();
    ShanChenSingleComponentProcessor3D(ShanChenSingleComponentProcessor3D<T,Descriptor> const& rhs);
    ShanChenSingleComponentProcessor3D& operator=(ShanChenSingleComponentProcessor3D<T,Descriptor> const& rhs);
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice );
    virtual ShanChenSingleComponentProcessor3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    T G;
    interparticlePotential::PsiFunction<T>* Psi;
};

/// Shan-Chen coupling for multi-component flow with or without external force
/// but with external rhoBar and j.
template<typename T, template<typename U> class Descriptor>
class ShanChenExternalMultiComponentProcessor3D :
    public BoxProcessingFunctional3D
{
public:
    /// With these constructors, space- and time-dependent values of the
    ///   relaxation parameters omega are accounted for.
    ShanChenExternalMultiComponentProcessor3D(T G_);
    ShanChenExternalMultiComponentProcessor3D(std::vector<std::vector<T> > const& speciesG_);
    /// With these constructors, the values of the relaxation parameters omega are
    ///   taken to be species-dependent, but not space- or time-dependent. Their
    ///   value is imposed in the constructor.
    ShanChenExternalMultiComponentProcessor3D(T G_, std::vector<T> const& imposedOmega_);
    ShanChenExternalMultiComponentProcessor3D(std::vector<std::vector<T> > const& speciesG_, std::vector<T> const& imposedOmega_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks);
    virtual ShanChenExternalMultiComponentProcessor3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    T G;
    std::vector<T> speciesG;
    std::vector<T> imposedOmega;
};

}

#endif  // SHAN_CHEN_PROCESSOR_3D_H
