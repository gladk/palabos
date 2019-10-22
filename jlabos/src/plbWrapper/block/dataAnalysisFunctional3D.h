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
#ifndef SWIG_DATA_ANALYSIS_FUNCTIONAL_3D_H
#define SWIG_DATA_ANALYSIS_FUNCTIONAL_3D_H

#include "core/globalDefs.h"
#include "atomicBlock/dataProcessingFunctional3D.h"
#include "atomicBlock/reductiveDataProcessingFunctional3D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"

namespace plb {

template<typename T>
class BoxNTensorSumFunctional3D : public ReductiveBoxProcessingFunctional3D_N<T>
{
public:
    BoxNTensorSumFunctional3D(plint ndim);
    virtual void process(Box3D domain, NTensorField3D<T>& tensorField);
    virtual BoxNTensorSumFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    std::vector<T> getSumVector() const;
private:
    std::vector<plint> sumVectorId;
};

template<typename T>
class MaskedBoxNTensorSumFunctional3D : public MaskedReductiveBoxProcessingFunctional3D_N<T>
{
public:
    MaskedBoxNTensorSumFunctional3D(plint ndim);
    virtual void process( Box3D domain,
                          NTensorField3D<T>& tensorField,
                          NTensorField3D<int>& mask);
    virtual MaskedBoxNTensorSumFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    std::vector<T> getSumVector() const;
private:
    std::vector<plint> sumVectorId;
};

template<typename T>
class BoxNTensorMinFunctional3D : public ReductiveBoxProcessingFunctional3D_N<T>
{
public:
    BoxNTensorMinFunctional3D(plint ndim);
    virtual void process(Box3D domain, NTensorField3D<T>& tensorField);
    virtual BoxNTensorMinFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    std::vector<T> getMinVector() const;
private:
    std::vector<plint> maxVectorId;
};

template<typename T>
class MaskedBoxNTensorMinFunctional3D : public MaskedReductiveBoxProcessingFunctional3D_N<T>
{
public:
    MaskedBoxNTensorMinFunctional3D(plint ndim);
    virtual void process( Box3D domain,
                          NTensorField3D<T>& tensorField,
                          NTensorField3D<int>& mask);
    virtual MaskedBoxNTensorMinFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    std::vector<T> getMinVector() const;
private:
    std::vector<plint> maxVectorId;
};

template<typename T>
class BoxNTensorMaxFunctional3D : public ReductiveBoxProcessingFunctional3D_N<T>
{
public:
    BoxNTensorMaxFunctional3D(plint ndim);
    virtual void process(Box3D domain, NTensorField3D<T>& tensorField);
    virtual BoxNTensorMaxFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    std::vector<T> getMaxVector() const;
private:
    std::vector<plint> maxVectorId;
};

template<typename T>
class MaskedBoxNTensorMaxFunctional3D : public MaskedReductiveBoxProcessingFunctional3D_N<T>
{
public:
    MaskedBoxNTensorMaxFunctional3D(plint ndim);
    virtual void process( Box3D domain,
                          NTensorField3D<T>& tensorField,
                          NTensorField3D<int>& mask );
    virtual MaskedBoxNTensorMaxFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    std::vector<T> getMaxVector() const;
private:
    std::vector<plint> maxVectorId;
};

template<typename T>
class BoundedBoxNTensorSumFunctional3D : public BoundedReductiveBoxProcessingFunctional3D_N<T>
{
public:
    BoundedBoxNTensorSumFunctional3D(plint ndim);
    virtual void processBulk(Box3D domain, NTensorField3D<T>& tensorField);
    virtual void processPlane( int direction, int orientation, Box3D domain,
                               NTensorField3D<T>& tensorField );
    virtual void processEdge( int plane, int normal1, int normal2,
                              Box3D domain, NTensorField3D<T>& tensorField );
    virtual void processCorner( int normalX, int normalY, int normalZ,
                                Box3D domain, NTensorField3D<T>& tensorField );
    virtual BoundedBoxNTensorSumFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    std::vector<T> getSumVector() const;
private:
    std::vector<plint> sumVectorId;
};

template<typename T>
class BoundedMaskedBoxNTensorSumFunctional3D : public BoundedMaskedReductiveBoxProcessingFunctional3D_N<T>
{
public:
    BoundedMaskedBoxNTensorSumFunctional3D(plint ndim);
    virtual void processBulk(Box3D domain, NTensorField3D<T>& tensorField, NTensorField3D<int>& mask);
    virtual void processPlane( int direction, int orientation, Box3D domain,
                               NTensorField3D<T>& tensorField, NTensorField3D<int>& mask );
    virtual void processEdge( int plane, int normal1, int normal2,
                              Box3D domain, NTensorField3D<T>& tensorField, NTensorField3D<int>& mask );
    virtual void processCorner( int normalX, int normalY, int normalZ,
                                Box3D domain, NTensorField3D<T>& tensorField, NTensorField3D<int>& mask );
    virtual BoundedMaskedBoxNTensorSumFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    std::vector<T> getSumVector() const;
private:
    std::vector<plint> sumVectorId;
};


/* *************** Data Functionals for scalar-fields **************** */

template<typename T1, typename T2>
class MaskedCopyConvertNTensorFunctional3D : public MaskedBoxProcessingFunctional3D_NN<T1,T2>
{
public:
    virtual void process( Box3D domain,
                          NTensorField3D<T1>& field1, NTensorField3D<T2>& field2,
                          NTensorField3D<int>& mask );
    virtual MaskedCopyConvertNTensorFunctional3D<T1,T2>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_plus_alpha_NTensor3D : public BoxProcessingFunctional3D_NN<T,T>
{
public:
    A_plus_alpha_NTensor3D(std::vector<T> const& alpha_);
    virtual void process(Box3D domain, NTensorField3D<T>& A,
                                       NTensorField3D<T>& result);
    virtual A_plus_alpha_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Masked_A_plus_alpha_NTensor3D : public MaskedBoxProcessingFunctional3D_NN<T,T>
{
public:
    Masked_A_plus_alpha_NTensor3D(std::vector<T> const& alpha_);
    virtual void process(Box3D domain, NTensorField3D<T>& A,
                                       NTensorField3D<T>& result,
                                       NTensorField3D<int>& mask);
    virtual Masked_A_plus_alpha_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class A_minus_alpha_NTensor3D : public BoxProcessingFunctional3D_NN<T,T>
{
public:
    A_minus_alpha_NTensor3D(std::vector<T> const& alpha_);
    virtual void process(Box3D domain, NTensorField3D<T>& A,
                                       NTensorField3D<T>& result);
    virtual A_minus_alpha_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Masked_A_minus_alpha_NTensor3D : public MaskedBoxProcessingFunctional3D_NN<T,T>
{
public:
    Masked_A_minus_alpha_NTensor3D(std::vector<T> const& alpha_);
    virtual void process(Box3D domain, NTensorField3D<T>& A,
                                       NTensorField3D<T>& result,
                                       NTensorField3D<int>& mask);
    virtual Masked_A_minus_alpha_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Alpha_minus_A_NTensor3D : public BoxProcessingFunctional3D_NN<T,T>
{
public:
    Alpha_minus_A_NTensor3D(std::vector<T> const& alpha_);
    virtual void process(Box3D domain, NTensorField3D<T>& A,
                                       NTensorField3D<T>& result);
    virtual Alpha_minus_A_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Masked_Alpha_minus_A_NTensor3D : public MaskedBoxProcessingFunctional3D_NN<T,T>
{
public:
    Masked_Alpha_minus_A_NTensor3D(std::vector<T> const& alpha_);
    virtual void process(Box3D domain, NTensorField3D<T>& A,
                                       NTensorField3D<T>& result,
                                       NTensorField3D<int>& mask);
    virtual Masked_Alpha_minus_A_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class A_times_alpha_NTensor3D : public BoxProcessingFunctional3D_NN<T,T>
{
public:
    A_times_alpha_NTensor3D(std::vector<T> const& alpha_);
    virtual void process(Box3D domain, NTensorField3D<T>& A,
                                       NTensorField3D<T>& result);
    virtual A_times_alpha_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Masked_A_times_alpha_NTensor3D : public MaskedBoxProcessingFunctional3D_NN<T,T>
{
public:
    Masked_A_times_alpha_NTensor3D(std::vector<T> const& alpha_);
    virtual void process(Box3D domain, NTensorField3D<T>& A,
                                       NTensorField3D<T>& result,
                                       NTensorField3D<int>& mask);
    virtual Masked_A_times_alpha_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class A_dividedBy_alpha_NTensor3D : public BoxProcessingFunctional3D_NN<T,T>
{
public:
    A_dividedBy_alpha_NTensor3D(std::vector<T> const& alpha_);
    virtual void process(Box3D domain, NTensorField3D<T>& A,
                                       NTensorField3D<T>& result);
    virtual A_dividedBy_alpha_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Masked_A_dividedBy_alpha_NTensor3D : public MaskedBoxProcessingFunctional3D_NN<T,T>
{
public:
    Masked_A_dividedBy_alpha_NTensor3D(std::vector<T> const& alpha_);
    virtual void process(Box3D domain, NTensorField3D<T>& A,
                                       NTensorField3D<T>& result,
                                       NTensorField3D<int>& mask);
    virtual Masked_A_dividedBy_alpha_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Alpha_dividedBy_A_NTensor3D : public BoxProcessingFunctional3D_NN<T,T>
{
public:
    Alpha_dividedBy_A_NTensor3D(std::vector<T> const& alpha_);
    virtual void process(Box3D domain, NTensorField3D<T>& A,
                                       NTensorField3D<T>& result);
    virtual Alpha_dividedBy_A_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Masked_Alpha_dividedBy_A_NTensor3D : public MaskedBoxProcessingFunctional3D_NN<T,T>
{
public:
    Masked_Alpha_dividedBy_A_NTensor3D(std::vector<T> const& alpha_);
    virtual void process(Box3D domain, NTensorField3D<T>& A,
                                       NTensorField3D<T>& result,
                                       NTensorField3D<int>& mask);
    virtual Masked_Alpha_dividedBy_A_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class A_toThePower_alpha_NTensor3D : public BoxProcessingFunctional3D_NN<T,T>
{
public:
    A_toThePower_alpha_NTensor3D(std::vector<T> const& alpha_);
    virtual void process(Box3D domain, NTensorField3D<T>& A,
                                       NTensorField3D<T>& result);
    virtual A_toThePower_alpha_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Masked_A_toThePower_alpha_NTensor3D : public MaskedBoxProcessingFunctional3D_NN<T,T>
{
public:
    Masked_A_toThePower_alpha_NTensor3D(std::vector<T> const& alpha_);
    virtual void process(Box3D domain, NTensorField3D<T>& A,
                                       NTensorField3D<T>& result,
                                       NTensorField3D<int>& mask);
    virtual Masked_A_toThePower_alpha_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Alpha_toThePower_A_NTensor3D : public BoxProcessingFunctional3D_NN<T,T>
{
public:
    Alpha_toThePower_A_NTensor3D(std::vector<T> const& alpha_);
    virtual void process(Box3D domain, NTensorField3D<T>& A,
                                       NTensorField3D<T>& result);
    virtual Alpha_toThePower_A_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Masked_Alpha_toThePower_A_NTensor3D : public MaskedBoxProcessingFunctional3D_NN<T,T>
{
public:
    Masked_Alpha_toThePower_A_NTensor3D(std::vector<T> const& alpha_);
    virtual void process(Box3D domain, NTensorField3D<T>& A,
                                       NTensorField3D<T>& result,
                                       NTensorField3D<int>& mask);
    virtual Masked_Alpha_toThePower_A_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class A_plus_alpha_inplace_NTensor3D : public BoxProcessingFunctional3D_N<T>
{
public:
    A_plus_alpha_inplace_NTensor3D(std::vector<T> const& alpha_);
    virtual void process(Box3D domain, NTensorField3D<T>& A);
    virtual A_plus_alpha_inplace_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Masked_A_plus_alpha_inplace_NTensor3D : public MaskedBoxProcessingFunctional3D_N<T>
{
public:
    Masked_A_plus_alpha_inplace_NTensor3D(std::vector<T> const& alpha_);
    virtual void process( Box3D domain,
                          NTensorField3D<T>& A,
                          NTensorField3D<int>& mask);
    virtual Masked_A_plus_alpha_inplace_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class A_minus_alpha_inplace_NTensor3D : public BoxProcessingFunctional3D_N<T>
{
public:
    A_minus_alpha_inplace_NTensor3D(std::vector<T> const& alpha_);
    virtual void process(Box3D domain, NTensorField3D<T>& A);
    virtual A_minus_alpha_inplace_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Masked_A_minus_alpha_inplace_NTensor3D : public MaskedBoxProcessingFunctional3D_N<T>
{
public:
    Masked_A_minus_alpha_inplace_NTensor3D(std::vector<T> const& alpha_);
    virtual void process( Box3D domain,
                          NTensorField3D<T>& A,
                          NTensorField3D<int>& mask );
    virtual Masked_A_minus_alpha_inplace_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class A_times_alpha_inplace_NTensor3D : public BoxProcessingFunctional3D_N<T>
{
public:
    A_times_alpha_inplace_NTensor3D(std::vector<T> const& alpha_);
    virtual void process(Box3D domain, NTensorField3D<T>& A);
    virtual A_times_alpha_inplace_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Masked_A_times_alpha_inplace_NTensor3D : public MaskedBoxProcessingFunctional3D_N<T>
{
public:
    Masked_A_times_alpha_inplace_NTensor3D(std::vector<T> const& alpha_);
    virtual void process( Box3D domain,
                          NTensorField3D<T>& A,
                          NTensorField3D<int>& masked );
    virtual Masked_A_times_alpha_inplace_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class A_dividedBy_alpha_inplace_NTensor3D : public BoxProcessingFunctional3D_N<T>
{
public:
    A_dividedBy_alpha_inplace_NTensor3D(std::vector<T> const& alpha_);
    virtual void process(Box3D domain, NTensorField3D<T>& A);
    virtual A_dividedBy_alpha_inplace_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Masked_A_dividedBy_alpha_inplace_NTensor3D : public MaskedBoxProcessingFunctional3D_N<T>
{
public:
    Masked_A_dividedBy_alpha_inplace_NTensor3D(std::vector<T> const& alpha_);
    virtual void process( Box3D domain,
                          NTensorField3D<T>& A,
                          NTensorField3D<int>& mask );
    virtual Masked_A_dividedBy_alpha_inplace_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class A_toThePower_alpha_inplace_NTensor3D : public BoxProcessingFunctional3D_N<T>
{
public:
    A_toThePower_alpha_inplace_NTensor3D(std::vector<T> const& alpha_);
    virtual void process(Box3D domain, NTensorField3D<T>& A);
    virtual A_toThePower_alpha_inplace_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Masked_A_toThePower_alpha_inplace_NTensor3D : public MaskedBoxProcessingFunctional3D_N<T>
{
public:
    Masked_A_toThePower_alpha_inplace_NTensor3D(std::vector<T> const& alpha_);
    virtual void process( Box3D domain,
                          NTensorField3D<T>& A,
                          NTensorField3D<int>& mask );
    virtual Masked_A_toThePower_alpha_inplace_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};


template<typename T>
class A_plus_B_NTensor3D : public NTensorFieldBoxProcessingFunctional3D<T>
{
public:
    virtual void process(Box3D domain, std::vector<NTensorField3D<T>*> tensorFields);
    virtual A_plus_B_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class Masked_A_plus_B_NTensor3D : public MaskedNTensorFieldBoxProcessingFunctional3D<T>
{
public:
    virtual void process( Box3D domain,
                          std::vector<NTensorField3D<T>*> tensorFields,
                          NTensorField3D<int>& mask);
    virtual Masked_A_plus_B_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_minus_B_NTensor3D : public NTensorFieldBoxProcessingFunctional3D<T>
{
public:
    virtual void process(Box3D domain, std::vector<NTensorField3D<T>*> tensorFields);
    virtual A_minus_B_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class Masked_A_minus_B_NTensor3D : public MaskedNTensorFieldBoxProcessingFunctional3D<T>
{
public:
    virtual void process( Box3D domain,
                          std::vector<NTensorField3D<T>*> tensorFields,
                          NTensorField3D<int>& mask );
    virtual Masked_A_minus_B_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_times_B_NTensor3D : public NTensorFieldBoxProcessingFunctional3D<T>
{
public:
    virtual void process(Box3D domain, std::vector<NTensorField3D<T>*> tensorFields);
    virtual A_times_B_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class Masked_A_times_B_NTensor3D : public MaskedNTensorFieldBoxProcessingFunctional3D<T>
{
public:
    virtual void process( Box3D domain,
                          std::vector<NTensorField3D<T>*> tensorFields,
                          NTensorField3D<int>& mask );
    virtual Masked_A_times_B_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_dividedBy_B_NTensor3D : public NTensorFieldBoxProcessingFunctional3D<T>
{
public:
    virtual void process(Box3D domain, std::vector<NTensorField3D<T>*> tensorFields);
    virtual A_dividedBy_B_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class Masked_A_dividedBy_B_NTensor3D : public MaskedNTensorFieldBoxProcessingFunctional3D<T>
{
public:
    virtual void process( Box3D domain,
                          std::vector<NTensorField3D<T>*> tensorFields,
                          NTensorField3D<int>& mask );
    virtual Masked_A_dividedBy_B_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_toThePower_B_NTensor3D : public NTensorFieldBoxProcessingFunctional3D<T>
{
public:
    virtual void process(Box3D domain, std::vector<NTensorField3D<T>*> tensorFields);
    virtual A_toThePower_B_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class Masked_A_toThePower_B_NTensor3D : public MaskedNTensorFieldBoxProcessingFunctional3D<T>
{
public:
    virtual void process( Box3D domain,
                          std::vector<NTensorField3D<T>*> tensorFields,
                          NTensorField3D<int>& mask );
    virtual Masked_A_toThePower_B_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_equals_B_NTensor3D : public BoxProcessingFunctional3D
{
public:
    virtual void processGenericBlocks (
            Box3D domain, std::vector<AtomicBlock3D*> blocks);
    virtual A_equals_B_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class Masked_A_equals_B_NTensor3D : public BoxProcessingFunctional3D
{
public:
    virtual void processGenericBlocks (
            Box3D domain, std::vector<AtomicBlock3D*> blocks);
    virtual Masked_A_equals_B_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_equals_alpha_NTensor3D : public BoxProcessingFunctional3D_NN<T,int>
{
public:
    A_equals_alpha_NTensor3D(std::vector<T> const& alpha_);
    virtual void process(Box3D domain, NTensorField3D<T>& A,
                                       NTensorField3D<int>& result);
    virtual A_equals_alpha_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Masked_A_equals_alpha_NTensor3D : public MaskedBoxProcessingFunctional3D_NN<T,int>
{
public:
    Masked_A_equals_alpha_NTensor3D(std::vector<T> const& alpha_);
    virtual void process(Box3D domain, NTensorField3D<T>& A,
                                       NTensorField3D<int>& result,
                                       NTensorField3D<int>& mask );
    virtual Masked_A_equals_alpha_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class A_lessThan_B_NTensor3D : public BoxProcessingFunctional3D
{
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks);
    virtual A_lessThan_B_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class Masked_A_lessThan_B_NTensor3D : public BoxProcessingFunctional3D
{
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks);
    virtual Masked_A_lessThan_B_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_lessThan_alpha_NTensor3D : public BoxProcessingFunctional3D_NN<T,int>
{
public:
    A_lessThan_alpha_NTensor3D(std::vector<T> const& alpha_);
    virtual void process(Box3D domain, NTensorField3D<T>& A,
                                       NTensorField3D<int>& result);
    virtual A_lessThan_alpha_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Masked_A_lessThan_alpha_NTensor3D : public MaskedBoxProcessingFunctional3D_NN<T,int>
{
public:
    Masked_A_lessThan_alpha_NTensor3D(std::vector<T> const& alpha_);
    virtual void process(Box3D domain, NTensorField3D<T>& A,
                                       NTensorField3D<int>& result,
                                       NTensorField3D<int>& mask);
    virtual Masked_A_lessThan_alpha_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class A_lessEqual_B_NTensor3D : public BoxProcessingFunctional3D
{
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> tensorFields);
    virtual A_lessEqual_B_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class Masked_A_lessEqual_B_NTensor3D : public BoxProcessingFunctional3D
{
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> tensorFields);
    virtual Masked_A_lessEqual_B_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_lessEqual_alpha_NTensor3D : public BoxProcessingFunctional3D_NN<T,int>
{
public:
    A_lessEqual_alpha_NTensor3D(std::vector<T> const& alpha_);
    virtual void process(Box3D domain, NTensorField3D<T>& A,
                                       NTensorField3D<int>& result);
    virtual A_lessEqual_alpha_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Masked_A_lessEqual_alpha_NTensor3D : public MaskedBoxProcessingFunctional3D_NN<T,int>
{
public:
    Masked_A_lessEqual_alpha_NTensor3D(std::vector<T> const& alpha_);
    virtual void process(Box3D domain, NTensorField3D<T>& A,
                                       NTensorField3D<int>& result,
                                       NTensorField3D<int>& mask );
    virtual Masked_A_lessEqual_alpha_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class A_greaterThan_B_NTensor3D : public BoxProcessingFunctional3D
{
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> tensorFields);
    virtual A_greaterThan_B_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class Masked_A_greaterThan_B_NTensor3D : public BoxProcessingFunctional3D
{
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> tensorFields);
    virtual Masked_A_greaterThan_B_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_greaterThan_alpha_NTensor3D : public BoxProcessingFunctional3D_NN<T,int>
{
public:
    A_greaterThan_alpha_NTensor3D(std::vector<T> const& alpha_);
    virtual void process(Box3D domain, NTensorField3D<T>& A,
                                       NTensorField3D<int>& result);
    virtual A_greaterThan_alpha_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Masked_A_greaterThan_alpha_NTensor3D : public MaskedBoxProcessingFunctional3D_NN<T,int>
{
public:
    Masked_A_greaterThan_alpha_NTensor3D(std::vector<T> const& alpha_);
    virtual void process(Box3D domain, NTensorField3D<T>& A,
                                       NTensorField3D<int>& result,
                                       NTensorField3D<int>& mask);
    virtual Masked_A_greaterThan_alpha_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class A_greaterEqual_B_NTensor3D : public BoxProcessingFunctional3D
{
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> tensorFields);
    virtual A_greaterEqual_B_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class Masked_A_greaterEqual_B_NTensor3D : public BoxProcessingFunctional3D
{
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> tensorFields);
    virtual Masked_A_greaterEqual_B_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_greaterEqual_alpha_NTensor3D : public BoxProcessingFunctional3D_NN<T,int>
{
public:
    A_greaterEqual_alpha_NTensor3D(std::vector<T> const& alpha_);
    virtual void process(Box3D domain, NTensorField3D<T>& A,
                                       NTensorField3D<int>& result);
    virtual A_greaterEqual_alpha_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class Masked_A_greaterEqual_alpha_NTensor3D : public MaskedBoxProcessingFunctional3D_NN<T,int>
{
public:
    Masked_A_greaterEqual_alpha_NTensor3D(std::vector<T> const& alpha_);
    virtual void process(Box3D domain, NTensorField3D<T>& A,
                                       NTensorField3D<int>& result,
                                       NTensorField3D<int>& mask);
    virtual Masked_A_greaterEqual_alpha_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    std::vector<T> alpha;
};

template<typename T>
class A_and_B_NTensor3D : public BoxProcessingFunctional3D
{
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> tensorFields);
    virtual A_and_B_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class Masked_A_and_B_NTensor3D : public BoxProcessingFunctional3D
{
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> tensorFields);
    virtual Masked_A_and_B_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_or_B_NTensor3D : public BoxProcessingFunctional3D
{
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> tensorFields);
    virtual A_or_B_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class Masked_A_or_B_NTensor3D : public BoxProcessingFunctional3D
{
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> tensorFields);
    virtual Masked_A_or_B_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class Not_A_NTensor3D : public BoxProcessingFunctional3D_NN<T,int>
{
public:
    virtual void process(Box3D domain, NTensorField3D<T>& A,
                                       NTensorField3D<int>& result);
    virtual Not_A_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class Masked_Not_A_NTensor3D : public MaskedBoxProcessingFunctional3D_NN<T,int>
{
public:
    virtual void process(Box3D domain, NTensorField3D<T>& A,
                                       NTensorField3D<int>& result,
                                       NTensorField3D<int>& mask );
    virtual Masked_Not_A_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_plus_B_inplace_NTensor3D : public BoxProcessingFunctional3D_NN<T,T>
{
public:
    virtual void process(Box3D domain, NTensorField3D<T>& A, NTensorField3D<T>& B);
    virtual A_plus_B_inplace_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class Masked_A_plus_B_inplace_NTensor3D : public MaskedBoxProcessingFunctional3D_NN<T,T>
{
public:
    virtual void process(Box3D domain, NTensorField3D<T>& A,
                                       NTensorField3D<T>& B,
                                       NTensorField3D<int>& mask);
    virtual Masked_A_plus_B_inplace_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_minus_B_inplace_NTensor3D : public BoxProcessingFunctional3D_NN<T,T>
{
public:
    virtual void process(Box3D domain, NTensorField3D<T>& A, NTensorField3D<T>& B);
    virtual A_minus_B_inplace_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class Masked_A_minus_B_inplace_NTensor3D : public MaskedBoxProcessingFunctional3D_NN<T,T>
{
public:
    virtual void process(Box3D domain,
                         NTensorField3D<T>& A, NTensorField3D<T>& B,
                         NTensorField3D<int>& mask );
    virtual Masked_A_minus_B_inplace_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_times_B_inplace_NTensor3D : public BoxProcessingFunctional3D_NN<T,T>
{
public:
    virtual void process(Box3D domain, NTensorField3D<T>& A, NTensorField3D<T>& B);
    virtual A_times_B_inplace_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class Masked_A_times_B_inplace_NTensor3D : public MaskedBoxProcessingFunctional3D_NN<T,T>
{
public:
    virtual void process(Box3D domain,
                         NTensorField3D<T>& A, NTensorField3D<T>& B,
                         NTensorField3D<int>& mask );
    virtual Masked_A_times_B_inplace_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_dividedBy_B_inplace_NTensor3D : public BoxProcessingFunctional3D_NN<T,T>
{
public:
    virtual void process(Box3D domain, NTensorField3D<T>& A, NTensorField3D<T>& B);
    virtual A_dividedBy_B_inplace_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class Masked_A_dividedBy_B_inplace_NTensor3D : public MaskedBoxProcessingFunctional3D_NN<T,T>
{
public:
    virtual void process(Box3D domain,
                         NTensorField3D<T>& A, NTensorField3D<T>& B,
                         NTensorField3D<int>& mask);
    virtual Masked_A_dividedBy_B_inplace_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_toThePower_B_inplace_NTensor3D : public BoxProcessingFunctional3D_NN<T,T>
{
public:
    virtual void process(Box3D domain, NTensorField3D<T>& A, NTensorField3D<T>& B);
    virtual A_toThePower_B_inplace_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class Masked_A_toThePower_B_inplace_NTensor3D : public MaskedBoxProcessingFunctional3D_NN<T,T>
{
public:
    virtual void process(Box3D domain,
                         NTensorField3D<T>& A, NTensorField3D<T>& B,
                         NTensorField3D<int>& mask);
    virtual Masked_A_toThePower_B_inplace_NTensor3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class ExtractNTensorComponentFunctional3D :
    public BoxProcessingFunctional3D_NN<T,T>
{
public:
    ExtractNTensorComponentFunctional3D(int iComponent_);
    virtual void process(Box3D domain, NTensorField3D<T>& scalarField,
                                       NTensorField3D<T>& tensorField);
    virtual ExtractNTensorComponentFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    int iComponent;
};

template<typename T>
class MaskedExtractNTensorComponentFunctional3D :
    public MaskedBoxProcessingFunctional3D_NN<T,T>
{
public:
    MaskedExtractNTensorComponentFunctional3D(int iComponent_);
    virtual void process(Box3D domain, NTensorField3D<T>& scalarField,
                                       NTensorField3D<T>& tensorField,
                                       NTensorField3D<int>& mask);
    virtual MaskedExtractNTensorComponentFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    int iComponent;
};

template<typename T>
class ComputeNTensorNormFunctional3D : public BoxProcessingFunctional3D_NN<T,T>
{
public:
    virtual void process(Box3D domain, NTensorField3D<T>& scalarField,
                                       NTensorField3D<T>& tensorField);
    virtual ComputeNTensorNormFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class MaskedComputeNTensorNormFunctional3D : public MaskedBoxProcessingFunctional3D_NN<T,T>
{
public:
    virtual void process(Box3D domain, NTensorField3D<T>& scalarField,
                                       NTensorField3D<T>& tensorField,
                                       NTensorField3D<int>& mask);
    virtual MaskedComputeNTensorNormFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class ComputeNTensorNormSqrFunctional3D : public BoxProcessingFunctional3D_NN<T,T>
{
public:
    virtual void process(Box3D domain, NTensorField3D<T>& scalarField,
                                       NTensorField3D<T>& tensorField);
    virtual ComputeNTensorNormSqrFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class MaskedComputeNTensorNormSqrFunctional3D : public MaskedBoxProcessingFunctional3D_NN<T,T>
{
public:
    virtual void process(Box3D domain, NTensorField3D<T>& scalarField,
                                       NTensorField3D<T>& tensorField,
                                       NTensorField3D<int>& mask);
    virtual MaskedComputeNTensorNormSqrFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class ComputeSymmetricNTensorNormFunctional3D :
    public BoxProcessingFunctional3D_NN<T,T>
{
public:
    virtual void process(Box3D domain, NTensorField3D<T>& scalarField,
                                       NTensorField3D<T>& tensorField);
    virtual ComputeSymmetricNTensorNormFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class MaskedComputeSymmetricNTensorNormFunctional3D :
    public MaskedBoxProcessingFunctional3D_NN<T,T>
{
public:
    virtual void process(Box3D domain, NTensorField3D<T>& scalarField,
                                       NTensorField3D<T>& tensorField,
                                       NTensorField3D<int>& mask);
    virtual MaskedComputeSymmetricNTensorNormFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class ComputeSymmetricNTensorNormSqrFunctional3D :
    public BoxProcessingFunctional3D_NN<T,T>
{
public:
    virtual void process(Box3D domain, NTensorField3D<T>& scalarField,
                                       NTensorField3D<T>& tensorField);
    virtual ComputeSymmetricNTensorNormSqrFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class MaskedComputeSymmetricNTensorNormSqrFunctional3D :
    public MaskedBoxProcessingFunctional3D_NN<T,T>
{
public:
    virtual void process(Box3D domain, NTensorField3D<T>& scalarField,
                                       NTensorField3D<T>& tensorField,
                                       NTensorField3D<int>& mask);
    virtual MaskedComputeSymmetricNTensorNormSqrFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class ComputeSymmetricNTensorTraceFunctional3D :
    public BoxProcessingFunctional3D_NN<T,T>
{
public:
    virtual void process(Box3D domain, NTensorField3D<T>& scalarField,
                                       NTensorField3D<T>& tensorField);
    virtual ComputeSymmetricNTensorTraceFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class MaskedComputeSymmetricNTensorTraceFunctional3D :
    public MaskedBoxProcessingFunctional3D_NN<T,T>
{
public:
    virtual void process(Box3D domain, NTensorField3D<T>& scalarField,
                                       NTensorField3D<T>& tensorField,
                                       NTensorField3D<int>& mask);
    virtual MaskedComputeSymmetricNTensorTraceFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class BoxBulkNTensorVorticityFunctional3D : public BoxProcessingFunctional3D_NN<T,T>
{
public:
    virtual void process(Box3D domain, NTensorField3D<T>& vorticity,
                                       NTensorField3D<T>& velocity);
    virtual BoxBulkNTensorVorticityFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class MaskedBoxBulkNTensorVorticityFunctional3D : public MaskedBoxProcessingFunctional3D_NN<T,T>
{
public:
    virtual void process(Box3D domain, NTensorField3D<T>& vorticity,
                                       NTensorField3D<T>& velocity,
                                       NTensorField3D<int>& mask);
    virtual MaskedBoxBulkNTensorVorticityFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class BoxNTensorVorticityFunctional3D :
    public BoundedBoxProcessingFunctional3D_NN<T,T>
{
public:
    virtual void processBulk( Box3D domain, NTensorField3D<T>& velocity,
                                            NTensorField3D<T>& vorticity );
    virtual void processPlane( int direction, int orientation, Box3D domain,
                               NTensorField3D<T>& velocity,
                               NTensorField3D<T>& vorticity );
    virtual void processEdge( int plane, int normal1, int normal2, Box3D domain,
                              NTensorField3D<T>& velocity,
                              NTensorField3D<T>& vorticity );
    virtual void processCorner( int normalX, int normalY, int normalZ, Box3D domain,
                                NTensorField3D<T>& velocity,
                                NTensorField3D<T>& vorticity );
    virtual BoxNTensorVorticityFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class MaskedBoxNTensorVorticityFunctional3D :
    public BoundedMaskedBoxProcessingFunctional3D_NN<T,T>
{
public:
    virtual void processBulk( Box3D domain, NTensorField3D<T>& velocity,
                                            NTensorField3D<T>& vorticity,
                                            NTensorField3D<int>& mask);
    virtual void processPlane( int direction, int orientation, Box3D domain,
                               NTensorField3D<T>& velocity,
                               NTensorField3D<T>& vorticity,
                               NTensorField3D<int>& mask);
    virtual void processEdge( int plane, int normal1, int normal2, Box3D domain,
                              NTensorField3D<T>& velocity,
                              NTensorField3D<T>& vorticity,
                              NTensorField3D<int>& mask);
    virtual void processCorner( int normalX, int normalY, int normalZ, Box3D domain,
                                NTensorField3D<T>& velocity,
                                NTensorField3D<T>& vorticity,
                                NTensorField3D<int>& mask);
    virtual MaskedBoxNTensorVorticityFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class BoxBulkNTensorStrainRateFunctional3D : public BoxProcessingFunctional3D_NN<T,T>
{
public:
    virtual void process(Box3D domain, NTensorField3D<T>& velocity,
                                       NTensorField3D<T>& S);
    virtual BoxBulkNTensorStrainRateFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class MaskedBoxBulkNTensorStrainRateFunctional3D : public MaskedBoxProcessingFunctional3D_NN<T,T>
{
public:
    virtual void process(Box3D domain, NTensorField3D<T>& velocity,
                                       NTensorField3D<T>& S,
                                       NTensorField3D<int>& mask);
    virtual MaskedBoxBulkNTensorStrainRateFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class BoxNTensorStrainRateFunctional3D :
    public BoundedBoxProcessingFunctional3D_NN<T,T>
{
public:
    virtual void processBulk( Box3D domain, NTensorField3D<T>& velocity,
                                            NTensorField3D<T>& S );
    virtual void processPlane( int direction, int orientation, Box3D domain,
                               NTensorField3D<T>& velocity,
                               NTensorField3D<T>& S );
    virtual void processEdge( int plane, int normal1, int normal2, Box3D domain,
                              NTensorField3D<T>& velocity,
                              NTensorField3D<T>& S );
    virtual void processCorner( int normalX, int normalY, int normalZ, Box3D domain,
                                NTensorField3D<T>& velocity,
                                NTensorField3D<T>& S );
    virtual BoxNTensorStrainRateFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class MaskedBoxNTensorStrainRateFunctional3D :
    public BoundedMaskedBoxProcessingFunctional3D_NN<T,T>
{
public:
    virtual void processBulk( Box3D domain, NTensorField3D<T>& velocity,
                                            NTensorField3D<T>& S,
                                            NTensorField3D<int>& mask);
    virtual void processPlane( int direction, int orientation, Box3D domain,
                               NTensorField3D<T>& velocity,
                               NTensorField3D<T>& S,
                               NTensorField3D<int>& mask);
    virtual void processEdge( int plane, int normal1, int normal2, Box3D domain,
                              NTensorField3D<T>& velocity,
                              NTensorField3D<T>& S,
                              NTensorField3D<int>& mask);
    virtual void processCorner( int normalX, int normalY, int normalZ, Box3D domain,
                                NTensorField3D<T>& velocity,
                                NTensorField3D<T>& S,
                                NTensorField3D<int>& mask);
    virtual MaskedBoxNTensorStrainRateFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class UPO_ScalarProductFunctional3D : public ReductiveBoxProcessingFunctional3D_NN<T,T>
{
public:
    UPO_ScalarProductFunctional3D();
    virtual void process( Box3D domain,
                          NTensorField3D<T>& a,
                          NTensorField3D<T>& b );
    virtual UPO_ScalarProductFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    T getSum() const;
private:
    plint sumId;
};

template<typename T>
class Masked_UPO_ScalarProductFunctional3D : public MaskedReductiveBoxProcessingFunctional3D_NN<T,T>
{
public:
    Masked_UPO_ScalarProductFunctional3D();
    virtual void process( Box3D domain,
                          NTensorField3D<T>& a,
                          NTensorField3D<T>& b,
                          NTensorField3D<int>& mask );
    virtual Masked_UPO_ScalarProductFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    T getSum() const;
private:
    plint sumId;
};

}  // namespace plb

#endif  // SWIG_DATA_ANALYSIS_FUNCTIONAL_3D_H
