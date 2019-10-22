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
#ifndef SWIG_DATA_ANALYSIS_FUNCTIONAL_2D_H
#define SWIG_DATA_ANALYSIS_FUNCTIONAL_2D_H

#include "core/globalDefs.h"
#include "atomicBlock/dataProcessingFunctional2D.h"
#include "atomicBlock/reductiveDataProcessingFunctional2D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"

namespace plb {

template<typename T>
class BoxNTensorSumFunctional2D : public ReductiveBoxProcessingFunctional2D_N<T>
{
public:
    BoxNTensorSumFunctional2D(plint ndim);
    virtual void process(Box2D domain, NTensorField2D<T>& tensorField);
    virtual BoxNTensorSumFunctional2D<T>* clone() const;
    std::vector<T> getSumVector() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    std::vector<plint> sumVectorId;
};

template<typename T>
class MaskedBoxNTensorSumFunctional2D : public MaskedReductiveBoxProcessingFunctional2D_N<T>
{
public:
    MaskedBoxNTensorSumFunctional2D(plint ndim);
    virtual void process( Box2D domain,
                          NTensorField2D<T>& tensorField,
                          NTensorField2D<int>& mask);
    virtual MaskedBoxNTensorSumFunctional2D<T>* clone() const;
    std::vector<T> getSumVector() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    std::vector<plint> sumVectorId;
};

template<typename T>
class BoxNTensorMinFunctional2D : public ReductiveBoxProcessingFunctional2D_N<T>
{
public:
    BoxNTensorMinFunctional2D(plint ndim);
    virtual void process(Box2D domain, NTensorField2D<T>& tensorField);
    virtual BoxNTensorMinFunctional2D<T>* clone() const;
    std::vector<T> getMinVector() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    std::vector<plint> maxVectorId;
};

template<typename T>
class MaskedBoxNTensorMinFunctional2D : public MaskedReductiveBoxProcessingFunctional2D_N<T>
{
public:
    MaskedBoxNTensorMinFunctional2D(plint ndim);
    virtual void process( Box2D domain,
                          NTensorField2D<T>& tensorField,
                          NTensorField2D<int>& mask);
    virtual MaskedBoxNTensorMinFunctional2D<T>* clone() const;
    std::vector<T> getMinVector() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    std::vector<plint> maxVectorId;
};

template<typename T>
class BoxNTensorMaxFunctional2D : public ReductiveBoxProcessingFunctional2D_N<T>
{
public:
    BoxNTensorMaxFunctional2D(plint ndim);
    virtual void process(Box2D domain, NTensorField2D<T>& tensorField);
    virtual BoxNTensorMaxFunctional2D<T>* clone() const;
    std::vector<T> getMaxVector() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    std::vector<plint> maxVectorId;
};

template<typename T>
class MaskedBoxNTensorMaxFunctional2D : public MaskedReductiveBoxProcessingFunctional2D_N<T>
{
public:
    MaskedBoxNTensorMaxFunctional2D(plint ndim);
    virtual void process( Box2D domain,
                          NTensorField2D<T>& tensorField,
                          NTensorField2D<int>& mask );
    virtual MaskedBoxNTensorMaxFunctional2D<T>* clone() const;
    std::vector<T> getMaxVector() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    std::vector<plint> maxVectorId;
};

template<typename T>
class BoundedBoxNTensorSumFunctional2D : public BoundedReductiveBoxProcessingFunctional2D_N<T>
{
public:
    BoundedBoxNTensorSumFunctional2D(plint ndim);
    virtual void processBulk(Box2D domain, NTensorField2D<T>& tensorField);
    virtual void processEdge( int direction, int orientation, Box2D domain,
                              NTensorField2D<T>& tensorField );
    virtual void processCorner( int normalX, int normalY, Box2D domain,
                                NTensorField2D<T>& tensorField );
    virtual BoundedBoxNTensorSumFunctional2D<T>* clone() const;
    std::vector<T> getSumVector() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    std::vector<plint> sumVectorId;
};

template<typename T>
class BoundedMaskedBoxNTensorSumFunctional2D : public BoundedMaskedReductiveBoxProcessingFunctional2D_N<T>
{
public:
    BoundedMaskedBoxNTensorSumFunctional2D(plint ndim);
    virtual void processBulk( Box2D domain,
                              NTensorField2D<T>& tensorField,
                              NTensorField2D<int>& mask);
    virtual void processEdge( int direction, int orientation, Box2D domain,
                              NTensorField2D<T>& tensorField,
                              NTensorField2D<int>& mask );
    virtual void processCorner( int normalX, int normalY, Box2D domain,
                                NTensorField2D<T>& tensorField,
                                NTensorField2D<int>& mask );
    virtual BoundedMaskedBoxNTensorSumFunctional2D<T>* clone() const;
    std::vector<T> getSumVector() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    std::vector<plint> sumVectorId;
};


template<typename T1, typename T2>
class MaskedCopyConvertNTensorFunctional2D : public MaskedBoxProcessingFunctional2D_NN<T1,T2>
{
public:
    virtual void process( Box2D domain,
                          NTensorField2D<T1>& field1, NTensorField2D<T2>& field2,
                          NTensorField2D<int>& mask );
    virtual MaskedCopyConvertNTensorFunctional2D<T1,T2>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_plus_alpha_NTensor2D : public BoxProcessingFunctional2D_NN<T,T>
{
public:
    A_plus_alpha_NTensor2D(std::vector<T> const& alpha_);
    virtual void process(Box2D domain, NTensorField2D<T>& A,
                                       NTensorField2D<T>& result);
    virtual A_plus_alpha_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Masked_A_plus_alpha_NTensor2D : public MaskedBoxProcessingFunctional2D_NN<T,T>
{
public:
    Masked_A_plus_alpha_NTensor2D(std::vector<T> const& alpha_);
    virtual void process(Box2D domain, NTensorField2D<T>& A,
                                       NTensorField2D<T>& result,
                                       NTensorField2D<int>& mask);
    virtual Masked_A_plus_alpha_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class A_minus_alpha_NTensor2D : public BoxProcessingFunctional2D_NN<T,T>
{
public:
    A_minus_alpha_NTensor2D(std::vector<T> const& alpha_);
    virtual void process(Box2D domain, NTensorField2D<T>& A,
                                       NTensorField2D<T>& result);
    virtual A_minus_alpha_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Masked_A_minus_alpha_NTensor2D : public MaskedBoxProcessingFunctional2D_NN<T,T>
{
public:
    Masked_A_minus_alpha_NTensor2D(std::vector<T> const& alpha_);
    virtual void process(Box2D domain, NTensorField2D<T>& A,
                                       NTensorField2D<T>& result,
                                       NTensorField2D<int>& mask);
    virtual Masked_A_minus_alpha_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Alpha_minus_A_NTensor2D : public BoxProcessingFunctional2D_NN<T,T>
{
public:
    Alpha_minus_A_NTensor2D(std::vector<T> const& alpha_);
    virtual void process(Box2D domain, NTensorField2D<T>& A,
                                       NTensorField2D<T>& result);
    virtual Alpha_minus_A_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Masked_Alpha_minus_A_NTensor2D : public MaskedBoxProcessingFunctional2D_NN<T,T>
{
public:
    Masked_Alpha_minus_A_NTensor2D(std::vector<T> const& alpha_);
    virtual void process(Box2D domain, NTensorField2D<T>& A,
                                       NTensorField2D<T>& result,
                                       NTensorField2D<int>& mask);
    virtual Masked_Alpha_minus_A_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class A_times_alpha_NTensor2D : public BoxProcessingFunctional2D_NN<T,T>
{
public:
    A_times_alpha_NTensor2D(std::vector<T> const& alpha_);
    virtual void process(Box2D domain, NTensorField2D<T>& A,
                                       NTensorField2D<T>& result);
    virtual A_times_alpha_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Masked_A_times_alpha_NTensor2D : public MaskedBoxProcessingFunctional2D_NN<T,T>
{
public:
    Masked_A_times_alpha_NTensor2D(std::vector<T> const& alpha_);
    virtual void process(Box2D domain, NTensorField2D<T>& A,
                                       NTensorField2D<T>& result,
                                       NTensorField2D<int>& mask);
    virtual Masked_A_times_alpha_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class A_dividedBy_alpha_NTensor2D : public BoxProcessingFunctional2D_NN<T,T>
{
public:
    A_dividedBy_alpha_NTensor2D(std::vector<T> const& alpha_);
    virtual void process(Box2D domain, NTensorField2D<T>& A,
                                       NTensorField2D<T>& result);
    virtual A_dividedBy_alpha_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Masked_A_dividedBy_alpha_NTensor2D : public MaskedBoxProcessingFunctional2D_NN<T,T>
{
public:
    Masked_A_dividedBy_alpha_NTensor2D(std::vector<T> const& alpha_);
    virtual void process(Box2D domain, NTensorField2D<T>& A,
                                       NTensorField2D<T>& result,
                                       NTensorField2D<int>& mask);
    virtual Masked_A_dividedBy_alpha_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Alpha_dividedBy_A_NTensor2D : public BoxProcessingFunctional2D_NN<T,T>
{
public:
    Alpha_dividedBy_A_NTensor2D(std::vector<T> const& alpha_);
    virtual void process(Box2D domain, NTensorField2D<T>& A,
                                       NTensorField2D<T>& result);
    virtual Alpha_dividedBy_A_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Masked_Alpha_dividedBy_A_NTensor2D : public MaskedBoxProcessingFunctional2D_NN<T,T>
{
public:
    Masked_Alpha_dividedBy_A_NTensor2D(std::vector<T> const& alpha_);
    virtual void process(Box2D domain, NTensorField2D<T>& A,
                                       NTensorField2D<T>& result,
                                       NTensorField2D<int>& mask);
    virtual Masked_Alpha_dividedBy_A_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class A_toThePower_alpha_NTensor2D : public BoxProcessingFunctional2D_NN<T,T>
{
public:
    A_toThePower_alpha_NTensor2D(std::vector<T> const& alpha_);
    virtual void process(Box2D domain, NTensorField2D<T>& A,
                                       NTensorField2D<T>& result);
    virtual A_toThePower_alpha_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Masked_A_toThePower_alpha_NTensor2D : public MaskedBoxProcessingFunctional2D_NN<T,T>
{
public:
    Masked_A_toThePower_alpha_NTensor2D(std::vector<T> const& alpha_);
    virtual void process(Box2D domain, NTensorField2D<T>& A,
                                       NTensorField2D<T>& result,
                                       NTensorField2D<int>& mask);
    virtual Masked_A_toThePower_alpha_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Alpha_toThePower_A_NTensor2D : public BoxProcessingFunctional2D_NN<T,T>
{
public:
    Alpha_toThePower_A_NTensor2D(std::vector<T> const& alpha_);
    virtual void process(Box2D domain, NTensorField2D<T>& A,
                                       NTensorField2D<T>& result);
    virtual Alpha_toThePower_A_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Masked_Alpha_toThePower_A_NTensor2D : public MaskedBoxProcessingFunctional2D_NN<T,T>
{
public:
    Masked_Alpha_toThePower_A_NTensor2D(std::vector<T> const& alpha_);
    virtual void process(Box2D domain, NTensorField2D<T>& A,
                                       NTensorField2D<T>& result,
                                       NTensorField2D<int>& mask);
    virtual Masked_Alpha_toThePower_A_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class A_plus_alpha_inplace_NTensor2D : public BoxProcessingFunctional2D_N<T>
{
public:
    A_plus_alpha_inplace_NTensor2D(std::vector<T> const& alpha_);
    virtual void process(Box2D domain, NTensorField2D<T>& A);
    virtual A_plus_alpha_inplace_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Masked_A_plus_alpha_inplace_NTensor2D : public MaskedBoxProcessingFunctional2D_N<T>
{
public:
    Masked_A_plus_alpha_inplace_NTensor2D(std::vector<T> const& alpha_);
    virtual void process( Box2D domain,
                          NTensorField2D<T>& A,
                          NTensorField2D<int>& mask);
    virtual Masked_A_plus_alpha_inplace_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class A_minus_alpha_inplace_NTensor2D : public BoxProcessingFunctional2D_N<T>
{
public:
    A_minus_alpha_inplace_NTensor2D(std::vector<T> const& alpha_);
    virtual void process(Box2D domain, NTensorField2D<T>& A);
    virtual A_minus_alpha_inplace_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Masked_A_minus_alpha_inplace_NTensor2D : public MaskedBoxProcessingFunctional2D_N<T>
{
public:
    Masked_A_minus_alpha_inplace_NTensor2D(std::vector<T> const& alpha_);
    virtual void process( Box2D domain,
                          NTensorField2D<T>& A,
                          NTensorField2D<int>& mask );
    virtual Masked_A_minus_alpha_inplace_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class A_times_alpha_inplace_NTensor2D : public BoxProcessingFunctional2D_N<T>
{
public:
    A_times_alpha_inplace_NTensor2D(std::vector<T> const& alpha_);
    virtual void process(Box2D domain, NTensorField2D<T>& A);
    virtual A_times_alpha_inplace_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Masked_A_times_alpha_inplace_NTensor2D : public MaskedBoxProcessingFunctional2D_N<T>
{
public:
    Masked_A_times_alpha_inplace_NTensor2D(std::vector<T> const& alpha_);
    virtual void process( Box2D domain,
                          NTensorField2D<T>& A,
                          NTensorField2D<int>& masked );
    virtual Masked_A_times_alpha_inplace_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class A_dividedBy_alpha_inplace_NTensor2D : public BoxProcessingFunctional2D_N<T>
{
public:
    A_dividedBy_alpha_inplace_NTensor2D(std::vector<T> const& alpha_);
    virtual void process(Box2D domain, NTensorField2D<T>& A);
    virtual A_dividedBy_alpha_inplace_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Masked_A_dividedBy_alpha_inplace_NTensor2D : public MaskedBoxProcessingFunctional2D_N<T>
{
public:
    Masked_A_dividedBy_alpha_inplace_NTensor2D(std::vector<T> const& alpha_);
    virtual void process( Box2D domain,
                          NTensorField2D<T>& A,
                          NTensorField2D<int>& mask );
    virtual Masked_A_dividedBy_alpha_inplace_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class A_toThePower_alpha_inplace_NTensor2D : public BoxProcessingFunctional2D_N<T>
{
public:
    A_toThePower_alpha_inplace_NTensor2D(std::vector<T> const& alpha_);
    virtual void process(Box2D domain, NTensorField2D<T>& A);
    virtual A_toThePower_alpha_inplace_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Masked_A_toThePower_alpha_inplace_NTensor2D : public MaskedBoxProcessingFunctional2D_N<T>
{
public:
    Masked_A_toThePower_alpha_inplace_NTensor2D(std::vector<T> const& alpha_);
    virtual void process( Box2D domain,
                          NTensorField2D<T>& A,
                          NTensorField2D<int>& mask );
    virtual Masked_A_toThePower_alpha_inplace_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};


template<typename T>
class A_plus_B_NTensor2D : public NTensorFieldBoxProcessingFunctional2D<T>
{
public:
    virtual void process(Box2D domain, std::vector<NTensorField2D<T>*> tensorFields);
    virtual A_plus_B_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class Masked_A_plus_B_NTensor2D : public MaskedNTensorFieldBoxProcessingFunctional2D<T>
{
public:
    virtual void process( Box2D domain,
                          std::vector<NTensorField2D<T>*> tensorFields,
                          NTensorField2D<int>& mask);
    virtual Masked_A_plus_B_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_minus_B_NTensor2D : public NTensorFieldBoxProcessingFunctional2D<T>
{
public:
    virtual void process(Box2D domain, std::vector<NTensorField2D<T>*> tensorFields);
    virtual A_minus_B_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class Masked_A_minus_B_NTensor2D : public MaskedNTensorFieldBoxProcessingFunctional2D<T>
{
public:
    virtual void process( Box2D domain,
                          std::vector<NTensorField2D<T>*> tensorFields,
                          NTensorField2D<int>& mask );
    virtual Masked_A_minus_B_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_times_B_NTensor2D : public NTensorFieldBoxProcessingFunctional2D<T>
{
public:
    virtual void process(Box2D domain, std::vector<NTensorField2D<T>*> tensorFields);
    virtual A_times_B_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class Masked_A_times_B_NTensor2D : public MaskedNTensorFieldBoxProcessingFunctional2D<T>
{
public:
    virtual void process( Box2D domain,
                          std::vector<NTensorField2D<T>*> tensorFields,
                          NTensorField2D<int>& mask );
    virtual Masked_A_times_B_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_dividedBy_B_NTensor2D : public NTensorFieldBoxProcessingFunctional2D<T>
{
public:
    virtual void process(Box2D domain, std::vector<NTensorField2D<T>*> tensorFields);
    virtual A_dividedBy_B_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class Masked_A_dividedBy_B_NTensor2D : public MaskedNTensorFieldBoxProcessingFunctional2D<T>
{
public:
    virtual void process( Box2D domain,
                          std::vector<NTensorField2D<T>*> tensorFields,
                          NTensorField2D<int>& mask );
    virtual Masked_A_dividedBy_B_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_toThePower_B_NTensor2D : public NTensorFieldBoxProcessingFunctional2D<T>
{
public:
    virtual void process(Box2D domain, std::vector<NTensorField2D<T>*> tensorFields);
    virtual A_toThePower_B_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class Masked_A_toThePower_B_NTensor2D : public MaskedNTensorFieldBoxProcessingFunctional2D<T>
{
public:
    virtual void process( Box2D domain,
                          std::vector<NTensorField2D<T>*> tensorFields,
                          NTensorField2D<int>& mask );
    virtual Masked_A_toThePower_B_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_equals_B_NTensor2D : public BoxProcessingFunctional2D
{
public:
    virtual void processGenericBlocks (
            Box2D domain, std::vector<AtomicBlock2D*> blocks);
    virtual A_equals_B_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class Masked_A_equals_B_NTensor2D : public BoxProcessingFunctional2D
{
public:
    virtual void processGenericBlocks (
            Box2D domain, std::vector<AtomicBlock2D*> blocks);
    virtual Masked_A_equals_B_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_equals_alpha_NTensor2D : public BoxProcessingFunctional2D_NN<T,int>
{
public:
    A_equals_alpha_NTensor2D(std::vector<T> const& alpha_);
    virtual void process(Box2D domain, NTensorField2D<T>& A,
                                       NTensorField2D<int>& result);
    virtual A_equals_alpha_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Masked_A_equals_alpha_NTensor2D : public MaskedBoxProcessingFunctional2D_NN<T,int>
{
public:
    Masked_A_equals_alpha_NTensor2D(std::vector<T> const& alpha_);
    virtual void process(Box2D domain, NTensorField2D<T>& A,
                                       NTensorField2D<int>& result,
                                       NTensorField2D<int>& mask );
    virtual Masked_A_equals_alpha_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class A_lessThan_B_NTensor2D : public BoxProcessingFunctional2D
{
public:
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> blocks);
    virtual A_lessThan_B_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class Masked_A_lessThan_B_NTensor2D : public BoxProcessingFunctional2D
{
public:
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> blocks);
    virtual Masked_A_lessThan_B_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_lessThan_alpha_NTensor2D : public BoxProcessingFunctional2D_NN<T,int>
{
public:
    A_lessThan_alpha_NTensor2D(std::vector<T> const& alpha_);
    virtual void process(Box2D domain, NTensorField2D<T>& A,
                                       NTensorField2D<int>& result);
    virtual A_lessThan_alpha_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Masked_A_lessThan_alpha_NTensor2D : public MaskedBoxProcessingFunctional2D_NN<T,int>
{
public:
    Masked_A_lessThan_alpha_NTensor2D(std::vector<T> const& alpha_);
    virtual void process(Box2D domain, NTensorField2D<T>& A,
                                       NTensorField2D<int>& result,
                                       NTensorField2D<int>& mask);
    virtual Masked_A_lessThan_alpha_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class A_lessEqual_B_NTensor2D : public BoxProcessingFunctional2D
{
public:
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> tensorFields);
    virtual A_lessEqual_B_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class Masked_A_lessEqual_B_NTensor2D : public BoxProcessingFunctional2D
{
public:
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> tensorFields);
    virtual Masked_A_lessEqual_B_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_lessEqual_alpha_NTensor2D : public BoxProcessingFunctional2D_NN<T,int>
{
public:
    A_lessEqual_alpha_NTensor2D(std::vector<T> const& alpha_);
    virtual void process(Box2D domain, NTensorField2D<T>& A,
                                       NTensorField2D<int>& result);
    virtual A_lessEqual_alpha_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Masked_A_lessEqual_alpha_NTensor2D : public MaskedBoxProcessingFunctional2D_NN<T,int>
{
public:
    Masked_A_lessEqual_alpha_NTensor2D(std::vector<T> const& alpha_);
    virtual void process(Box2D domain, NTensorField2D<T>& A,
                                       NTensorField2D<int>& result,
                                       NTensorField2D<int>& mask );
    virtual Masked_A_lessEqual_alpha_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class A_greaterThan_B_NTensor2D : public BoxProcessingFunctional2D
{
public:
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> tensorFields);
    virtual A_greaterThan_B_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class Masked_A_greaterThan_B_NTensor2D : public BoxProcessingFunctional2D
{
public:
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> tensorFields);
    virtual Masked_A_greaterThan_B_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_greaterThan_alpha_NTensor2D : public BoxProcessingFunctional2D_NN<T,int>
{
public:
    A_greaterThan_alpha_NTensor2D(std::vector<T> const& alpha_);
    virtual void process(Box2D domain, NTensorField2D<T>& A,
                                       NTensorField2D<int>& result);
    virtual A_greaterThan_alpha_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Masked_A_greaterThan_alpha_NTensor2D : public MaskedBoxProcessingFunctional2D_NN<T,int>
{
public:
    Masked_A_greaterThan_alpha_NTensor2D(std::vector<T> const& alpha_);
    virtual void process(Box2D domain, NTensorField2D<T>& A,
                                       NTensorField2D<int>& result,
                                       NTensorField2D<int>& mask);
    virtual Masked_A_greaterThan_alpha_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class A_greaterEqual_B_NTensor2D : public BoxProcessingFunctional2D
{
public:
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> tensorFields);
    virtual A_greaterEqual_B_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class Masked_A_greaterEqual_B_NTensor2D : public BoxProcessingFunctional2D
{
public:
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> tensorFields);
    virtual Masked_A_greaterEqual_B_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_greaterEqual_alpha_NTensor2D : public BoxProcessingFunctional2D_NN<T,int>
{
public:
    A_greaterEqual_alpha_NTensor2D(std::vector<T> const& alpha_);
    virtual void process(Box2D domain, NTensorField2D<T>& A,
                                       NTensorField2D<int>& result);
    virtual A_greaterEqual_alpha_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Masked_A_greaterEqual_alpha_NTensor2D : public MaskedBoxProcessingFunctional2D_NN<T,int>
{
public:
    Masked_A_greaterEqual_alpha_NTensor2D(std::vector<T> const& alpha_);
    virtual void process(Box2D domain, NTensorField2D<T>& A,
                                       NTensorField2D<int>& result,
                                       NTensorField2D<int>& mask);
    virtual Masked_A_greaterEqual_alpha_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class A_and_B_NTensor2D : public BoxProcessingFunctional2D
{
public:
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> tensorFields);
    virtual A_and_B_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class Masked_A_and_B_NTensor2D : public BoxProcessingFunctional2D
{
public:
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> tensorFields);
    virtual Masked_A_and_B_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_or_B_NTensor2D : public BoxProcessingFunctional2D
{
public:
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> tensorFields);
    virtual A_or_B_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class Masked_A_or_B_NTensor2D : public BoxProcessingFunctional2D
{
public:
    virtual void processGenericBlocks(Box2D domain, std::vector<AtomicBlock2D*> tensorFields);
    virtual Masked_A_or_B_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class Not_A_NTensor2D : public BoxProcessingFunctional2D_NN<T,int>
{
public:
    virtual void process(Box2D domain, NTensorField2D<T>& A,
                                       NTensorField2D<int>& result);
    virtual Not_A_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class Masked_Not_A_NTensor2D : public MaskedBoxProcessingFunctional2D_NN<T,int>
{
public:
    virtual void process(Box2D domain, NTensorField2D<T>& A,
                                       NTensorField2D<int>& result,
                                       NTensorField2D<int>& mask );
    virtual Masked_Not_A_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_plus_B_inplace_NTensor2D : public BoxProcessingFunctional2D_NN<T,T>
{
public:
    virtual void process(Box2D domain, NTensorField2D<T>& A, NTensorField2D<T>& B);
    virtual A_plus_B_inplace_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class Masked_A_plus_B_inplace_NTensor2D : public MaskedBoxProcessingFunctional2D_NN<T,T>
{
public:
    virtual void process(Box2D domain, NTensorField2D<T>& A,
                                       NTensorField2D<T>& B,
                                       NTensorField2D<int>& mask);
    virtual Masked_A_plus_B_inplace_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_minus_B_inplace_NTensor2D : public BoxProcessingFunctional2D_NN<T,T>
{
public:
    virtual void process(Box2D domain, NTensorField2D<T>& A, NTensorField2D<T>& B);
    virtual A_minus_B_inplace_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class Masked_A_minus_B_inplace_NTensor2D : public MaskedBoxProcessingFunctional2D_NN<T,T>
{
public:
    virtual void process(Box2D domain,
                         NTensorField2D<T>& A, NTensorField2D<T>& B,
                         NTensorField2D<int>& mask );
    virtual Masked_A_minus_B_inplace_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_times_B_inplace_NTensor2D : public BoxProcessingFunctional2D_NN<T,T>
{
public:
    virtual void process(Box2D domain, NTensorField2D<T>& A, NTensorField2D<T>& B);
    virtual A_times_B_inplace_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class Masked_A_times_B_inplace_NTensor2D : public MaskedBoxProcessingFunctional2D_NN<T,T>
{
public:
    virtual void process(Box2D domain,
                         NTensorField2D<T>& A, NTensorField2D<T>& B,
                         NTensorField2D<int>& mask );
    virtual Masked_A_times_B_inplace_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_dividedBy_B_inplace_NTensor2D : public BoxProcessingFunctional2D_NN<T,T>
{
public:
    virtual void process(Box2D domain, NTensorField2D<T>& A, NTensorField2D<T>& B);
    virtual A_dividedBy_B_inplace_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class Masked_A_dividedBy_B_inplace_NTensor2D : public MaskedBoxProcessingFunctional2D_NN<T,T>
{
public:
    virtual void process(Box2D domain,
                         NTensorField2D<T>& A, NTensorField2D<T>& B,
                         NTensorField2D<int>& mask);
    virtual Masked_A_dividedBy_B_inplace_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_toThePower_B_inplace_NTensor2D : public BoxProcessingFunctional2D_NN<T,T>
{
public:
    virtual void process(Box2D domain, NTensorField2D<T>& A, NTensorField2D<T>& B);
    virtual A_toThePower_B_inplace_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class Masked_A_toThePower_B_inplace_NTensor2D : public MaskedBoxProcessingFunctional2D_NN<T,T>
{
public:
    virtual void process(Box2D domain,
                         NTensorField2D<T>& A, NTensorField2D<T>& B,
                         NTensorField2D<int>& mask);
    virtual Masked_A_toThePower_B_inplace_NTensor2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class ExtractNTensorComponentFunctional2D :
    public BoxProcessingFunctional2D_NN<T,T>
{
public:
    ExtractNTensorComponentFunctional2D(int iComponent_);
    virtual void process(Box2D domain, NTensorField2D<T>& scalarField,
                                       NTensorField2D<T>& tensorField);
    virtual ExtractNTensorComponentFunctional2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    int iComponent;
};

template<typename T>
class MaskedExtractNTensorComponentFunctional2D :
    public MaskedBoxProcessingFunctional2D_NN<T,T>
{
public:
    MaskedExtractNTensorComponentFunctional2D(int iComponent_);
    virtual void process(Box2D domain, NTensorField2D<T>& scalarField,
                                       NTensorField2D<T>& tensorField,
                                       NTensorField2D<int>& mask);
    virtual MaskedExtractNTensorComponentFunctional2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    int iComponent;
};

template<typename T>
class ComputeNTensorNormFunctional2D : public BoxProcessingFunctional2D_NN<T,T>
{
public:
    virtual void process(Box2D domain, NTensorField2D<T>& scalarField,
                                       NTensorField2D<T>& tensorField);
    virtual ComputeNTensorNormFunctional2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class MaskedComputeNTensorNormFunctional2D : public MaskedBoxProcessingFunctional2D_NN<T,T>
{
public:
    virtual void process(Box2D domain, NTensorField2D<T>& scalarField,
                                       NTensorField2D<T>& tensorField,
                                       NTensorField2D<int>& mask);
    virtual MaskedComputeNTensorNormFunctional2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class ComputeNTensorNormSqrFunctional2D : public BoxProcessingFunctional2D_NN<T,T>
{
public:
    virtual void process(Box2D domain, NTensorField2D<T>& scalarField,
                                       NTensorField2D<T>& tensorField);
    virtual ComputeNTensorNormSqrFunctional2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class MaskedComputeNTensorNormSqrFunctional2D : public MaskedBoxProcessingFunctional2D_NN<T,T>
{
public:
    virtual void process(Box2D domain, NTensorField2D<T>& scalarField,
                                       NTensorField2D<T>& tensorField,
                                       NTensorField2D<int>& mask);
    virtual MaskedComputeNTensorNormSqrFunctional2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class ComputeSymmetricNTensorNormFunctional2D :
    public BoxProcessingFunctional2D_NN<T,T>
{
public:
    virtual void process(Box2D domain, NTensorField2D<T>& scalarField,
                                       NTensorField2D<T>& tensorField);
    virtual ComputeSymmetricNTensorNormFunctional2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class MaskedComputeSymmetricNTensorNormFunctional2D :
    public MaskedBoxProcessingFunctional2D_NN<T,T>
{
public:
    virtual void process(Box2D domain, NTensorField2D<T>& scalarField,
                                       NTensorField2D<T>& tensorField,
                                       NTensorField2D<int>& mask);
    virtual MaskedComputeSymmetricNTensorNormFunctional2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class ComputeSymmetricNTensorNormSqrFunctional2D :
    public BoxProcessingFunctional2D_NN<T,T>
{
public:
    virtual void process(Box2D domain, NTensorField2D<T>& scalarField,
                                       NTensorField2D<T>& tensorField);
    virtual ComputeSymmetricNTensorNormSqrFunctional2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class MaskedComputeSymmetricNTensorNormSqrFunctional2D :
    public MaskedBoxProcessingFunctional2D_NN<T,T>
{
public:
    virtual void process(Box2D domain, NTensorField2D<T>& scalarField,
                                       NTensorField2D<T>& tensorField,
                                       NTensorField2D<int>& mask);
    virtual MaskedComputeSymmetricNTensorNormSqrFunctional2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class ComputeSymmetricNTensorTraceFunctional2D :
    public BoxProcessingFunctional2D_NN<T,T>
{
public:
    virtual void process(Box2D domain, NTensorField2D<T>& scalarField,
                                       NTensorField2D<T>& tensorField);
    virtual ComputeSymmetricNTensorTraceFunctional2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class MaskedComputeSymmetricNTensorTraceFunctional2D :
    public MaskedBoxProcessingFunctional2D_NN<T,T>
{
public:
    virtual void process(Box2D domain, NTensorField2D<T>& scalarField,
                                       NTensorField2D<T>& tensorField,
                                       NTensorField2D<int>& mask);
    virtual MaskedComputeSymmetricNTensorTraceFunctional2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class BoxBulkNTensorVorticityFunctional2D : public BoxProcessingFunctional2D_NN<T,T>
{
public:
    virtual void process(Box2D domain, NTensorField2D<T>& vorticity,
                                       NTensorField2D<T>& velocity);
    virtual BoxBulkNTensorVorticityFunctional2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class MaskedBoxBulkNTensorVorticityFunctional2D : public MaskedBoxProcessingFunctional2D_NN<T,T>
{
public:
    virtual void process(Box2D domain, NTensorField2D<T>& vorticity,
                                       NTensorField2D<T>& velocity,
                                       NTensorField2D<int>& mask);
    virtual MaskedBoxBulkNTensorVorticityFunctional2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class BoxNTensorVorticityFunctional2D : public BoundedBoxProcessingFunctional2D_NN<T,T>
{
public:
    virtual void processBulk( Box2D domain, NTensorField2D<T>& vorticity,
                                            NTensorField2D<T>& velocity );
    virtual void processEdge( int direction, int orientation, Box2D domain,
                              NTensorField2D<T>& vorticity,
                              NTensorField2D<T>& velocity );
    virtual void processCorner( int normalX, int normalY, Box2D domain,
                                NTensorField2D<T>& vorticity,
                                NTensorField2D<T>& velocity );
    virtual BoxNTensorVorticityFunctional2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class MaskedBoxNTensorVorticityFunctional2D : public BoundedMaskedBoxProcessingFunctional2D_NN<T,T>
{
public:
    virtual void processBulk( Box2D domain, NTensorField2D<T>& vorticity,
                                            NTensorField2D<T>& velocity,
                                            NTensorField2D<int>& mask );
    virtual void processEdge( int direction, int orientation, Box2D domain,
                              NTensorField2D<T>& vorticity,
                              NTensorField2D<T>& velocity,
                              NTensorField2D<int>& mask );
    virtual void processCorner( int normalX, int normalY, Box2D domain,
                                NTensorField2D<T>& vorticity,
                                NTensorField2D<T>& velocity,
                                NTensorField2D<int>& mask );
    virtual MaskedBoxNTensorVorticityFunctional2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class BoxBulkNTensorStrainRateFunctional2D : public BoxProcessingFunctional2D_NN<T,T>
{
public:
    virtual void process(Box2D domain, NTensorField2D<T>& velocity,
                                       NTensorField2D<T>& S);
    virtual BoxBulkNTensorStrainRateFunctional2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class MaskedBoxBulkNTensorStrainRateFunctional2D : public MaskedBoxProcessingFunctional2D_NN<T,T>
{
public:
    virtual void process(Box2D domain, NTensorField2D<T>& velocity,
                                       NTensorField2D<T>& S,
                                       NTensorField2D<int>& mask);
    virtual MaskedBoxBulkNTensorStrainRateFunctional2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class BoxNTensorStrainRateFunctional2D : public BoundedBoxProcessingFunctional2D_NN<T,T>
{
public:
    virtual void processBulk( Box2D domain, NTensorField2D<T>& velocity,
                                            NTensorField2D<T>& S );
    virtual void processEdge( int direction, int orientation, Box2D domain,
                              NTensorField2D<T>& velocity,
                              NTensorField2D<T>& S );
    virtual void processCorner( int normalX, int normalY, Box2D domain,
                                NTensorField2D<T>& velocity,
                                NTensorField2D<T>& S );
    virtual BoxNTensorStrainRateFunctional2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class MaskedBoxNTensorStrainRateFunctional2D : public BoundedMaskedBoxProcessingFunctional2D_NN<T,T>
{
public:
    virtual void processBulk( Box2D domain, NTensorField2D<T>& velocity,
                                            NTensorField2D<T>& S,
                                            NTensorField2D<int>& mask );
    virtual void processEdge( int direction, int orientation, Box2D domain,
                              NTensorField2D<T>& velocity,
                              NTensorField2D<T>& S,
                              NTensorField2D<int>& mask );
    virtual void processCorner( int normalX, int normalY, Box2D domain,
                                NTensorField2D<T>& velocity,
                                NTensorField2D<T>& S,
                                NTensorField2D<int>& mask );
    virtual MaskedBoxNTensorStrainRateFunctional2D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

}  // namespace plb

#endif  // SWIG_DATA_ANALYSIS_FUNCTIONAL_2D_H
