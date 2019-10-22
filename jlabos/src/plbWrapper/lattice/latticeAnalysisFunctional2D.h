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
#ifndef SWIG_LATTICE_ANALYSIS_FUNCTIONAL_2D_H
#define SWIG_LATTICE_ANALYSIS_FUNCTIONAL_2D_H

#include "core/globalDefs.h"
#include "atomicBlock/dataProcessingFunctional2D.h"
#include "atomicBlock/reductiveDataProcessingFunctional2D.h"

namespace plb {


template<typename T, template<typename U> class Descriptor> 
class N_BoxDensityFunctional2D : public BoxProcessingFunctional2D_LN<T,Descriptor,T>
{
public:
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
                                       NTensorField2D<T>& scalarField);
    virtual N_BoxDensityFunctional2D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class Masked_N_BoxDensityFunctional2D : public MaskedBoxProcessingFunctional2D_LN<T,Descriptor,T>
{
public:
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
                                       NTensorField2D<T>& scalarField,
                                       NTensorField2D<int>& mask );
    virtual Masked_N_BoxDensityFunctional2D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class N_BoxKineticEnergyFunctional2D : public BoxProcessingFunctional2D_LN<T,Descriptor,T>
{
public:
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
                                       NTensorField2D<T>& scalarField);
    virtual N_BoxKineticEnergyFunctional2D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class Masked_N_BoxKineticEnergyFunctional2D : public MaskedBoxProcessingFunctional2D_LN<T,Descriptor,T>
{
public:
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
                                       NTensorField2D<T>& scalarField,
                                       NTensorField2D<int>& mask );
    virtual Masked_N_BoxKineticEnergyFunctional2D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class N_BoxVelocityNormFunctional2D : public BoxProcessingFunctional2D_LN<T,Descriptor,T>
{
public:
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
                                       NTensorField2D<T>& scalarField);
    virtual N_BoxVelocityNormFunctional2D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class Masked_N_BoxVelocityNormFunctional2D : public MaskedBoxProcessingFunctional2D_LN<T,Descriptor,T>
{
public:
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
                                       NTensorField2D<T>& scalarField,
                                       NTensorField2D<int>& mask );
    virtual Masked_N_BoxVelocityNormFunctional2D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class N_BoxVelocityComponentFunctional2D :
        public BoxProcessingFunctional2D_LN<T,Descriptor,T>
{
public:
    N_BoxVelocityComponentFunctional2D(int iComponent_);
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
                                       NTensorField2D<T>& scalarField);
    virtual N_BoxVelocityComponentFunctional2D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    int iComponent;
};

template<typename T, template<typename U> class Descriptor> 
class Masked_N_BoxVelocityComponentFunctional2D :
        public MaskedBoxProcessingFunctional2D_LN<T,Descriptor,T>
{
public:
    Masked_N_BoxVelocityComponentFunctional2D(int iComponent_);
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
                                       NTensorField2D<T>& scalarField,
                                       NTensorField2D<int>& mask );
    virtual Masked_N_BoxVelocityComponentFunctional2D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    int iComponent;
};

template<typename T, template<typename U> class Descriptor> 
class N_BoxVelocityFunctional2D :
    public BoxProcessingFunctional2D_LN<T,Descriptor, T>
{
public:
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
                                       NTensorField2D<T>& tensorField);
    virtual N_BoxVelocityFunctional2D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class Masked_N_BoxVelocityFunctional2D :
    public MaskedBoxProcessingFunctional2D_LN<T,Descriptor, T>
{
public:
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
                                       NTensorField2D<T>& tensorField,
                                       NTensorField2D<int>& mask );
    virtual Masked_N_BoxVelocityFunctional2D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class N_BoxPiNeqFunctional2D :
    public BoxProcessingFunctional2D_LN<T,Descriptor, T>
{
public:
    virtual void process(Box2D domain,
                         BlockLattice2D<T,Descriptor>& lattice,
                         NTensorField2D<T>& PiNeq);
    virtual N_BoxPiNeqFunctional2D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class Masked_N_BoxPiNeqFunctional2D :
    public MaskedBoxProcessingFunctional2D_LN<T,Descriptor, T>
{
public:
    virtual void process(Box2D domain,
                         BlockLattice2D<T,Descriptor>& lattice,
                         NTensorField2D<T>& PiNeq,
                         NTensorField2D<int>& mask);
    virtual Masked_N_BoxPiNeqFunctional2D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class N_BoxShearStressFunctional2D :
    public BoxProcessingFunctional2D_LN<T,Descriptor, T>
{
public:
    virtual void process(Box2D domain,
                         BlockLattice2D<T,Descriptor>& lattice,
                         NTensorField2D<T>& ShearStress);
    virtual N_BoxShearStressFunctional2D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class Masked_N_BoxShearStressFunctional2D :
    public MaskedBoxProcessingFunctional2D_LN<T,Descriptor, T>
{
public:
    virtual void process(Box2D domain,
                         BlockLattice2D<T,Descriptor>& lattice,
                         NTensorField2D<T>& ShearStress,
                         NTensorField2D<int>& mask);
    virtual Masked_N_BoxShearStressFunctional2D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class N_BoxStrainRateFromStressFunctional2D :
    public BoxProcessingFunctional2D_LN<T,Descriptor, T>
{
public:
    virtual void process(Box2D domain,
                         BlockLattice2D<T,Descriptor>& lattice,
                         NTensorField2D<T>& S);
    virtual N_BoxStrainRateFromStressFunctional2D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class Masked_N_BoxStrainRateFromStressFunctional2D :
    public MaskedBoxProcessingFunctional2D_LN<T,Descriptor, T>
{
public:
    virtual void process(Box2D domain,
                         BlockLattice2D<T,Descriptor>& lattice,
                         NTensorField2D<T>& S,
                         NTensorField2D<int>& mask );
    virtual Masked_N_BoxStrainRateFromStressFunctional2D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class N_BoxPopulationFunctional2D : public BoxProcessingFunctional2D_LN<T,Descriptor,T>
{
public:
    N_BoxPopulationFunctional2D(plint iComponent_);
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
                                       NTensorField2D<T>& populations);
    virtual N_BoxPopulationFunctional2D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    plint iComponent;
};

template<typename T, template<typename U> class Descriptor> 
class Masked_N_BoxPopulationFunctional2D : public MaskedBoxProcessingFunctional2D_LN<T,Descriptor,T>
{
public:
    Masked_N_BoxPopulationFunctional2D(plint iComponent_);
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
                                       NTensorField2D<T>& populations,
                                       NTensorField2D<int>& mask);
    virtual Masked_N_BoxPopulationFunctional2D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    plint iComponent;
};

template<typename T, template<typename U> class Descriptor> 
class N_BoxPopulationsFunctional2D : public BoxProcessingFunctional2D_LN<T,Descriptor,T>
{
public:
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
                                       NTensorField2D<T>& population);
    virtual N_BoxPopulationsFunctional2D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class Masked_N_BoxPopulationsFunctional2D : public MaskedBoxProcessingFunctional2D_LN<T,Descriptor,T>
{
public:
    virtual void process(Box2D domain, BlockLattice2D<T,Descriptor>& lattice,
                                       NTensorField2D<T>& population,
                                       NTensorField2D<int>& mask);
    virtual Masked_N_BoxPopulationsFunctional2D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

}  // namespace plb

#endif  // SWIG_LATTICE_ANALYSIS_FUNCTIONAL_2D_H
