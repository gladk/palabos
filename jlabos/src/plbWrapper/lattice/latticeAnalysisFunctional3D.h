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
#ifndef SWIG_LATTICE_ANALYSIS_FUNCTIONAL_3D_H
#define SWIG_LATTICE_ANALYSIS_FUNCTIONAL_3D_H

#include "core/globalDefs.h"
#include "atomicBlock/dataProcessingFunctional3D.h"
#include "atomicBlock/reductiveDataProcessingFunctional3D.h"

namespace plb {


template<typename T, template<typename U> class Descriptor> 
class N_BoxDensityFunctional3D : public BoxProcessingFunctional3D_LN<T,Descriptor,T>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       NTensorField3D<T>& scalarField);
    virtual N_BoxDensityFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class Masked_N_BoxDensityFunctional3D : public MaskedBoxProcessingFunctional3D_LN<T,Descriptor,T>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       NTensorField3D<T>& scalarField,
                                       NTensorField3D<int>& mask );
    virtual Masked_N_BoxDensityFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class N_BoxKineticEnergyFunctional3D : public BoxProcessingFunctional3D_LN<T,Descriptor,T>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       NTensorField3D<T>& scalarField);
    virtual N_BoxKineticEnergyFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class Masked_N_BoxKineticEnergyFunctional3D : public MaskedBoxProcessingFunctional3D_LN<T,Descriptor,T>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       NTensorField3D<T>& scalarField,
                                       NTensorField3D<int>& mask );
    virtual Masked_N_BoxKineticEnergyFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class N_BoxVelocityNormFunctional3D : public BoxProcessingFunctional3D_LN<T,Descriptor,T>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       NTensorField3D<T>& scalarField);
    virtual N_BoxVelocityNormFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class Masked_N_BoxVelocityNormFunctional3D : public MaskedBoxProcessingFunctional3D_LN<T,Descriptor,T>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       NTensorField3D<T>& scalarField,
                                       NTensorField3D<int>& mask );
    virtual Masked_N_BoxVelocityNormFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class N_BoxVelocityComponentFunctional3D :
        public BoxProcessingFunctional3D_LN<T,Descriptor,T>
{
public:
    N_BoxVelocityComponentFunctional3D(int iComponent_);
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       NTensorField3D<T>& scalarField);
    virtual N_BoxVelocityComponentFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    int iComponent;
};

template<typename T, template<typename U> class Descriptor> 
class Masked_N_BoxVelocityComponentFunctional3D :
        public MaskedBoxProcessingFunctional3D_LN<T,Descriptor,T>
{
public:
    Masked_N_BoxVelocityComponentFunctional3D(int iComponent_);
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       NTensorField3D<T>& scalarField,
                                       NTensorField3D<int>& mask );
    virtual Masked_N_BoxVelocityComponentFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    int iComponent;
};

template<typename T, template<typename U> class Descriptor> 
class N_BoxVelocityFunctional3D :
    public BoxProcessingFunctional3D_LN<T,Descriptor, T>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       NTensorField3D<T>& tensorField);
    virtual N_BoxVelocityFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class Masked_N_BoxVelocityFunctional3D :
    public MaskedBoxProcessingFunctional3D_LN<T,Descriptor, T>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       NTensorField3D<T>& tensorField,
                                       NTensorField3D<int>& mask );
    virtual Masked_N_BoxVelocityFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class N_BoxPiNeqFunctional3D :
    public BoxProcessingFunctional3D_LN<T,Descriptor, T>
{
public:
    virtual void process(Box3D domain,
                         BlockLattice3D<T,Descriptor>& lattice,
                         NTensorField3D<T>& PiNeq);
    virtual N_BoxPiNeqFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class Masked_N_BoxPiNeqFunctional3D :
    public MaskedBoxProcessingFunctional3D_LN<T,Descriptor, T>
{
public:
    virtual void process(Box3D domain,
                         BlockLattice3D<T,Descriptor>& lattice,
                         NTensorField3D<T>& PiNeq,
                         NTensorField3D<int>& mask);
    virtual Masked_N_BoxPiNeqFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class N_BoxShearStressFunctional3D :
    public BoxProcessingFunctional3D_LN<T,Descriptor, T>
{
public:
    virtual void process(Box3D domain,
                         BlockLattice3D<T,Descriptor>& lattice,
                         NTensorField3D<T>& ShearStress);
    virtual N_BoxShearStressFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class Masked_N_BoxShearStressFunctional3D :
    public MaskedBoxProcessingFunctional3D_LN<T,Descriptor, T>
{
public:
    virtual void process(Box3D domain,
                         BlockLattice3D<T,Descriptor>& lattice,
                         NTensorField3D<T>& ShearStress,
                         NTensorField3D<int>& mask);
    virtual Masked_N_BoxShearStressFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class N_BoxStrainRateFromStressFunctional3D :
    public BoxProcessingFunctional3D_LN<T,Descriptor, T>
{
public:
    virtual void process(Box3D domain,
                         BlockLattice3D<T,Descriptor>& lattice,
                         NTensorField3D<T>& S);
    virtual N_BoxStrainRateFromStressFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class Masked_N_BoxStrainRateFromStressFunctional3D :
    public MaskedBoxProcessingFunctional3D_LN<T,Descriptor, T>
{
public:
    virtual void process(Box3D domain,
                         BlockLattice3D<T,Descriptor>& lattice,
                         NTensorField3D<T>& S,
                         NTensorField3D<int>& mask );
    virtual Masked_N_BoxStrainRateFromStressFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class N_BoxPopulationFunctional3D : public BoxProcessingFunctional3D_LN<T,Descriptor,T>
{
public:
    N_BoxPopulationFunctional3D(plint iComponent_);
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       NTensorField3D<T>& populations);
    virtual N_BoxPopulationFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    plint iComponent;
};

template<typename T, template<typename U> class Descriptor> 
class Masked_N_BoxPopulationFunctional3D : public MaskedBoxProcessingFunctional3D_LN<T,Descriptor,T>
{
public:
    Masked_N_BoxPopulationFunctional3D(plint iComponent_);
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       NTensorField3D<T>& populations,
                                       NTensorField3D<int>& mask);
    virtual Masked_N_BoxPopulationFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    plint iComponent;
};

template<typename T, template<typename U> class Descriptor> 
class N_BoxPopulationsFunctional3D : public BoxProcessingFunctional3D_LN<T,Descriptor,T>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       NTensorField3D<T>& population);
    virtual N_BoxPopulationsFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class Masked_N_BoxPopulationsFunctional3D : public MaskedBoxProcessingFunctional3D_LN<T,Descriptor,T>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       NTensorField3D<T>& population,
                                       NTensorField3D<int>& mask);
    virtual Masked_N_BoxPopulationsFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};


template<typename T, template<typename U> class Descriptor> 
class UPO_Rhs_Functional3D : public BoxProcessingFunctional3D_NN<T,T>
{
public:
    UPO_Rhs_Functional3D(T omega_);
    virtual void process(Box3D domain, NTensorField3D<T>& lattice,
                                       NTensorField3D<T>& result);
    virtual UPO_Rhs_Functional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T omega;
};

template<typename T, template<typename U> class Descriptor> 
class Masked_UPO_Rhs_Functional3D : public MaskedBoxProcessingFunctional3D_NN<T,T>
{
public:
    Masked_UPO_Rhs_Functional3D(T omega_);
    virtual void process(Box3D domain, NTensorField3D<T>& lattice,
                                       NTensorField3D<T>& result,
                                       NTensorField3D<int>& mask);
    virtual Masked_UPO_Rhs_Functional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T omega;
};

template<typename T, template<typename U> class Descriptor> 
class UPO_ApplyJ_Functional3D : public NTensorFieldBoxProcessingFunctional3D<T>
{
public:
    UPO_ApplyJ_Functional3D(T omega_);
    virtual void process(Box3D domain, std::vector<NTensorField3D<T>*> lattices);
    virtual UPO_ApplyJ_Functional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T omega;
};

template<typename T, template<typename U> class Descriptor> 
class Masked_UPO_ApplyJ_Functional3D : public MaskedNTensorFieldBoxProcessingFunctional3D<T>
{
public:
    Masked_UPO_ApplyJ_Functional3D(T omega_);
    virtual void process( Box3D domain,
                          std::vector<NTensorField3D<T>*> lattices,
                          NTensorField3D<int>& mask );
    virtual Masked_UPO_ApplyJ_Functional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T omega;
};

template<typename T, template<typename U> class Descriptor> 
class UPO_EnergyDerivative_Functional3D : public BoxProcessingFunctional3D_NN<T,T>
{
public:
    virtual void process( Box3D domain,
                          NTensorField3D<T>& lattice,
                          NTensorField3D<T>& derivative );
    virtual UPO_EnergyDerivative_Functional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class Masked_UPO_EnergyDerivative_Functional3D : public MaskedBoxProcessingFunctional3D_NN<T,T>
{
public:
    virtual void process( Box3D domain,
                          NTensorField3D<T>& lattice,
                          NTensorField3D<T>& derivative,
                          NTensorField3D<int>& mask );
    virtual Masked_UPO_EnergyDerivative_Functional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

}  // namespace plb

#endif  // SWIG_LATTICE_ANALYSIS_FUNCTIONAL_3D_H
