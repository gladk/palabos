/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2015 FlowKit Sarl
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
 * Data processors for data analysis -- header file.
 */

#ifndef DATA_ANALYSIS_FUNCTIONAL_3D_H
#define DATA_ANALYSIS_FUNCTIONAL_3D_H

#include "core/globalDefs.h"
#include "core/array.h"
#include "atomicBlock/dataProcessingFunctional3D.h"
#include "atomicBlock/reductiveDataProcessingFunctional3D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"

namespace plb {

/* *************** PART I ******************************************** */
/* *************** Analysis of the block-lattice ********************* */
/* ******************************************************************* */

/* *************** Reductive Data Functionals for BlockLattice ******* */

template<typename T, template<typename U> class Descriptor> 
class BoxSumRhoBarFunctional3D :
    public ReductiveBoxProcessingFunctional3D_L<T,Descriptor>
{
public:
    BoxSumRhoBarFunctional3D();
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice);
    virtual BoxSumRhoBarFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
    }
    T getSumRhoBar() const;
private:
    plint sumRhoBarId;
};

template<typename T, template<typename U> class Descriptor> 
class BoxSumEnergyFunctional3D :
    public ReductiveBoxProcessingFunctional3D_L<T,Descriptor>
{
public:
    BoxSumEnergyFunctional3D();
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice);
    virtual BoxSumEnergyFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
    }
    T getSumEnergy() const;
    virtual void getDimensionsX(std::vector<int>& dimensions) const;
    virtual void getDimensionsT(std::vector<int>& dimensions) const;
private:
    plint sumEnergyId;
};

template<typename T, template<typename U> class Descriptor, class BoolMask> 
class CountLatticeElementsFunctional3D :
    public ReductiveBoxProcessingFunctional3D_L<T,Descriptor>
{
public:
    CountLatticeElementsFunctional3D(BoolMask boolMask_);
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice);
    virtual CountLatticeElementsFunctional3D<T,Descriptor,BoolMask>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
    }
    plint getCount() const;
private:
    plint countId;
    BoolMask boolMask;
};


template<typename T, template<typename U> class Descriptor> 
class DensitySingleProbe3D :
    public ReductiveBoxProcessingFunctional3D_L<T,Descriptor>
{
public:
    DensitySingleProbe3D(std::vector<Array<T,3> > const& positions_);
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice);
    virtual DensitySingleProbe3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
    }
    std::vector<T> getDensities() const;
private:
    std::vector<Array<T,3> > positions;
    std::vector<plint> densityIds;
};

template<typename T, template<typename U> class Descriptor> 
class VelocitySingleProbe3D :
    public ReductiveBoxProcessingFunctional3D_L<T,Descriptor>
{
public:
    VelocitySingleProbe3D(std::vector<Array<T,3> > const& positions_);
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice);
    virtual VelocitySingleProbe3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
    }
    std::vector<Array<T,3> > getVelocities() const;
private:
    std::vector<Array<T,3> > positions;
    std::vector<Array<plint,3> > velIds;
};

template<typename T, template<typename U> class Descriptor> 
class VorticitySingleProbe3D :
    public ReductiveBoxProcessingFunctional3D_L<T,Descriptor>
{
public:
    VorticitySingleProbe3D(std::vector<Array<T,3> > const& positions_);
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice);
    virtual VorticitySingleProbe3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
    }
    std::vector<Array<T,3> > getVorticities() const;
private:
    std::vector<Array<T,3> > positions;
    std::vector<Array<plint,3> > vorticityIds;
};

/* *************** Data Functionals for BlockLattice ***************** */

template<typename T, template<typename U> class Descriptor> 
class CopyPopulationsFunctional3D :
    public BoxProcessingFunctional3D_LL<T,Descriptor,T,Descriptor>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& latticeFrom,
                                       BlockLattice3D<T,Descriptor>& latticeTo);
    virtual CopyPopulationsFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T1, template<typename U1> class Descriptor1, typename T2, template<typename U2> class Descriptor2> 
class CopyConvertPopulationsFunctional3D : public BoxProcessingFunctional3D_LL<T1,Descriptor1,T2,Descriptor2>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T1,Descriptor1>& latticeFrom,
                         BlockLattice3D<T2,Descriptor2>& latticeTo);
    virtual CopyConvertPopulationsFunctional3D<T1,Descriptor1,T2,Descriptor2>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class LatticeCopyAllFunctional3D :
    public BoxProcessingFunctional3D_LL<T,Descriptor,T,Descriptor>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& latticeFrom,
                                       BlockLattice3D<T,Descriptor>& latticeTo);
    virtual LatticeCopyAllFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class LatticeRegenerateFunctional3D :
    public BoxProcessingFunctional3D_LL<T,Descriptor,T,Descriptor>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& latticeFrom,
                                       BlockLattice3D<T,Descriptor>& latticeTo);
    virtual LatticeRegenerateFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
};

template<typename T, template<typename U> class Descriptor> 
class BoxDensityFunctional3D : public BoxProcessingFunctional3D_LS<T,Descriptor,T>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       ScalarField3D<T>& scalarField);
    virtual BoxDensityFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class BoxRhoBarFunctional3D : public BoxProcessingFunctional3D_LS<T,Descriptor,T>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       ScalarField3D<T>& scalarField);
    virtual BoxRhoBarFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class BoxRhoBarJfunctional3D : public BoxProcessingFunctional3D
{
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual BoxRhoBarJfunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
};

template<typename T, template<typename U> class Descriptor> 
class PackedRhoBarJfunctional3D : public BoxProcessingFunctional3D_LN<T,Descriptor,T>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       NTensorField3D<T>& rhoBarJ);
    virtual PackedRhoBarJfunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
};

template<typename T>
class DensityFromRhoBarJfunctional3D : public BoxProcessingFunctional3D_SN<T,T>
{
public:
    virtual void process(Box3D domain, ScalarField3D<T>& density,
                                       NTensorField3D<T>& rhoBarJ);
    virtual DensityFromRhoBarJfunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
};

template<typename T>
class VelocityFromRhoBarJfunctional3D : public BoxProcessingFunctional3D
{
public:
    VelocityFromRhoBarJfunctional3D(bool velIsJ_ = false);
    // 1: velocity (TensorField), 2: rhoBarJ (NTensorField)
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual VelocityFromRhoBarJfunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    bool velIsJ;
};

template<typename T, template<typename U> class Descriptor> 
class BoxKineticEnergyFunctional3D : public BoxProcessingFunctional3D_LS<T,Descriptor,T>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       ScalarField3D<T>& scalarField);
    virtual BoxKineticEnergyFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class BoxVelocityNormFunctional3D : public BoxProcessingFunctional3D_LS<T,Descriptor,T>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       ScalarField3D<T>& scalarField);
    virtual BoxVelocityNormFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class BoxVelocityComponentFunctional3D :
    public BoxProcessingFunctional3D_LS<T,Descriptor,T>
{
public:
    BoxVelocityComponentFunctional3D(int iComponent);
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       ScalarField3D<T>& scalarField);
    virtual BoxVelocityComponentFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    int iComponent;
};

template<typename T, template<typename U> class Descriptor> 
class BoxVelocityFunctional3D :
    public BoxProcessingFunctional3D_LT<T,Descriptor, T,Descriptor<T>::d>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       TensorField3D<T, Descriptor<T>::d>& tensorField);
    virtual BoxVelocityFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class BoxTemperatureFunctional3D : public BoxProcessingFunctional3D_LS<T,Descriptor,T>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       ScalarField3D<T>& scalarField);
    virtual BoxTemperatureFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class BoxPiNeqFunctional3D :
    public BoxProcessingFunctional3D_LT<T,Descriptor, T,SymmetricTensor<T,Descriptor>::n>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       TensorField3D<T, SymmetricTensor<T,Descriptor>::n>& PiNeq);
    virtual BoxPiNeqFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class BoxShearStressFunctional3D :
public BoxProcessingFunctional3D_LT<T,Descriptor, T,SymmetricTensor<T,Descriptor>::n>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                         TensorField3D<T, SymmetricTensor<T,Descriptor>::n>& PiNeq);
    virtual BoxShearStressFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class BoxStrainRateFromStressFunctional3D :
    public BoxProcessingFunctional3D_LT<T,Descriptor, T,SymmetricTensor<T,Descriptor>::n>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       TensorField3D<T, SymmetricTensor<T,Descriptor>::n>& S);
    virtual BoxStrainRateFromStressFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class BoxPopulationFunctional3D : public BoxProcessingFunctional3D_LS<T,Descriptor,T>
{
public:
    BoxPopulationFunctional3D(plint iComponent_);
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       ScalarField3D<T>& population);
    virtual BoxPopulationFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    plint iComponent;
};

template<typename T, template<typename U> class Descriptor> 
class BoxEquilibriumFunctional3D : public BoxProcessingFunctional3D_LS<T,Descriptor,T>
{
public:
    BoxEquilibriumFunctional3D(plint iComponent_);
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       ScalarField3D<T>& equilibrium);
    virtual BoxEquilibriumFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    plint iComponent;
};

template<typename T, template<typename U> class Descriptor> 
class BoxAllPopulationsFunctional3D : public BoxProcessingFunctional3D_LT<T,Descriptor,T,Descriptor<T>::q>
{
public:
    BoxAllPopulationsFunctional3D();
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       TensorField3D<T,Descriptor<T>::q>& populations);
    virtual BoxAllPopulationsFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class BoxAllEquilibriumFunctional3D : public BoxProcessingFunctional3D_LT<T,Descriptor,T,Descriptor<T>::q>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       TensorField3D<T,Descriptor<T>::q>& equilibrium);
    virtual BoxAllEquilibriumFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class BoxAllNonEquilibriumFunctional3D : public BoxProcessingFunctional3D_LT<T,Descriptor,T,Descriptor<T>::q>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       TensorField3D<T,Descriptor<T>::q>& nonEquilibrium);
    virtual BoxAllNonEquilibriumFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class BoxAllPopulationsToLatticeFunctional3D : public BoxProcessingFunctional3D_LT<T,Descriptor,T,Descriptor<T>::q>
{
public:
    BoxAllPopulationsToLatticeFunctional3D();
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                         TensorField3D<T,Descriptor<T>::q>& populations);
    virtual BoxAllPopulationsToLatticeFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor> 
class BoxOmegaFunctional3D : public BoxProcessingFunctional3D_LS<T,Descriptor,T>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                         ScalarField3D<T>& scalarField);
    virtual BoxOmegaFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor>
class BoxKinematicViscosityFunctional3D : public BoxProcessingFunctional3D_LS<T,Descriptor,T>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                         ScalarField3D<T>& scalarField);
    virtual BoxKinematicViscosityFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, template<typename U> class Descriptor>
class BoxExternalForceFunctional3D :
    public BoxProcessingFunctional3D_LT<T,Descriptor, T,Descriptor<T>::d>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       TensorField3D<T,Descriptor<T>::d>& tensorField);
    virtual BoxExternalForceFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
};

template<typename T, template<typename U> class Descriptor>
class BoxExternalScalarFunctional3D :
    public BoxProcessingFunctional3D_LS<T,Descriptor,T>
{
public:
    BoxExternalScalarFunctional3D(int whichScalar_);
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       ScalarField3D<T>& scalarField);
    virtual BoxExternalScalarFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    int whichScalar;
};

template<typename T, template<typename U> class Descriptor>
class BoxExternalVectorFunctional3D :
    public BoxProcessingFunctional3D_LT<T,Descriptor,T,Descriptor<T>::d>
{
public:
    BoxExternalVectorFunctional3D(int vectorBeginsAt_);
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       TensorField3D<T,Descriptor<T>::d>& tensorField);
    virtual BoxExternalVectorFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    int vectorBeginsAt;
};

template<typename T, template<typename U> class Descriptor>
class BoxDynamicParameterFunctional3D :
    public BoxProcessingFunctional3D_LS<T,Descriptor,T>
{
public:
    BoxDynamicParameterFunctional3D(plint whichScalar_);
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       ScalarField3D<T>& scalarField);
    virtual BoxDynamicParameterFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    plint whichParameter;
};

template<typename T, template<typename U> class Descriptor>
class BoxDynamicViscosityFunctional3D :
    public BoxProcessingFunctional3D_LS<T,Descriptor,T>
{
public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       ScalarField3D<T>& scalarField);
    virtual BoxDynamicViscosityFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
};

/* *************** PART II ******************************************* */
/* *************** Analysis of the scalar-field ********************** */
/* ******************************************************************* */

/* *************** Reductive Data Functionals for ScalarField ******** */

/** Attention: No matter what the type of T is (even if it is an integer type),
 *    the sum is computed in double-precision floating point numbers, and
 *    converted to T at the end (and rounded, if T is an integer).
 **/
template<typename T>
class BoxScalarSumFunctional3D : public ReductiveBoxProcessingFunctional3D_S<T>
{
public:
    BoxScalarSumFunctional3D();
    virtual void process(Box3D domain, ScalarField3D<T>& scalarField);
    virtual BoxScalarSumFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
    }
    T getSumScalar() const;
private:
    plint sumScalarId;
};

/** Attention: No matter what the type of T is (even if it is an integer type),
 *    the average is computed in double-precision floating point numbers, and
 *    converted to T at the end (and rounded, if T is an integer).
 **/
template<typename T>
class MaskedBoxScalarAverageFunctional3D : public ReductiveBoxProcessingFunctional3D_SS<T,int>
{
public:
    MaskedBoxScalarAverageFunctional3D(int flag_);
    virtual void process( Box3D domain,
                          ScalarField3D<T>& scalarField,
                          ScalarField3D<int>& mask );
    virtual MaskedBoxScalarAverageFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
        modified[1] = modif::nothing;
    }
    T getAverageScalar() const;
private:
    plint averageScalarId;
    int flag;
};

/** Attention: No matter what the type of T is (even if it is an integer type),
 *    the min is computed in double-precision floating point numbers, and
 *    converted to T at the end (and rounded, if T is an integer).
 **/
template<typename T>
class BoxScalarMinFunctional3D : public ReductiveBoxProcessingFunctional3D_S<T>
{
public:
    BoxScalarMinFunctional3D();
    virtual void process(Box3D domain, ScalarField3D<T>& scalarField);
    virtual BoxScalarMinFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
    }
    T getMinScalar() const;
private:
    plint maxScalarId;
};

/** Attention: No matter what the type of T is (even if it is an integer type),
 *    the max is computed in double-precision floating point numbers, and
 *    converted to T at the end (and rounded, if T is an integer).
 **/
template<typename T>
class BoxScalarMaxFunctional3D : public ReductiveBoxProcessingFunctional3D_S<T>
{
public:
    BoxScalarMaxFunctional3D();
    virtual void process(Box3D domain, ScalarField3D<T>& scalarField);
    virtual BoxScalarMaxFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
    }
    T getMaxScalar() const;
private:
    plint maxScalarId;
};

/** Attention: No matter what the type of T is (even if it is an integer type),
 *    the sum is computed in double-precision floating point numbers, and
 *    converted to T at the end (and rounded, if T is an integer).
 **/
template<typename T>
class BoundedBoxScalarSumFunctional3D :
    public BoundedReductiveBoxProcessingFunctional3D_S<T>
{
public:
    BoundedBoxScalarSumFunctional3D();
    virtual void processBulk(Box3D domain, ScalarField3D<T>& scalarField);
    virtual void processPlane( int direction, int orientation, Box3D domain,
                               ScalarField3D<T>& scalarField );
    virtual void processEdge( int plane, int normal1, int normal2,
                              Box3D domain, ScalarField3D<T>& scalarField );
    virtual void processCorner( int normalX, int normalY, int normalZ,
                                Box3D domain, ScalarField3D<T>& scalarField );
    virtual BoundedBoxScalarSumFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
    }
    T getSumScalar() const;
private:
    plint sumScalarId;
};


template<typename T, class BoolMask> 
class CountScalarElementsFunctional3D : public ReductiveBoxProcessingFunctional3D_S<T>
{
public:
    CountScalarElementsFunctional3D(BoolMask boolMask_);
    virtual void process(Box3D domain, ScalarField3D<T>& field);
    virtual CountScalarElementsFunctional3D<T,BoolMask>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
    }
    plint getCount() const;
private:
    plint countId;
    BoolMask boolMask;
};


/* *************** Data Functionals for scalar-fields **************** */

template<typename T>
class ComputeScalarSqrtFunctional3D : public BoxProcessingFunctional3D_SS<T,T>
{
public:
    virtual void process(Box3D domain, ScalarField3D<T>& A,
                                       ScalarField3D<T>& B);
    virtual ComputeScalarSqrtFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class ComputeAbsoluteValueFunctional3D : public BoxProcessingFunctional3D_SS<T,T>
{
public:
    virtual void process(Box3D domain, ScalarField3D<T>& A,
                                       ScalarField3D<T>& B);
    virtual ComputeAbsoluteValueFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T1, typename T2>
class CopyConvertScalarFunctional3D : public BoxProcessingFunctional3D_SS<T1,T2>
{
public:
    virtual void process(Box3D domain, ScalarField3D<T1>& field1, ScalarField3D<T2>& field2);
    virtual CopyConvertScalarFunctional3D<T1,T2>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class ExtractScalarSubDomainFunctional3D : public BoxProcessingFunctional3D_SS<T,T>
{
public:
    virtual void process(Box3D domain, ScalarField3D<T>& field1,
                                       ScalarField3D<T>& field2);
    virtual ExtractScalarSubDomainFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, class Function>
class ApplyScalarFunctional3D : public BoxProcessingFunctional3D_S<T> {
public:
    ApplyScalarFunctional3D(Function f_);
    virtual void process(Box3D domain, ScalarField3D<T>& field);
    virtual ApplyScalarFunctional3D<T,Function>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    Function f;
};

template<typename T, class EvalFunction>
class EvaluateScalarFunctional3D : public BoxProcessingFunctional3D_SS<T,T> {
public:
    EvaluateScalarFunctional3D(EvalFunction f_);
    virtual void process(Box3D domain, ScalarField3D<T>& field,
                                       ScalarField3D<T>& result);
    virtual EvaluateScalarFunctional3D<T,EvalFunction>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    EvalFunction f;
};

template<typename T>
class A_plus_alpha_functional3D : public BoxProcessingFunctional3D_SS<T,T>
{
public:
    A_plus_alpha_functional3D(T alpha_);
    virtual void process(Box3D domain, ScalarField3D<T>& A,
                                       ScalarField3D<T>& result);
    virtual A_plus_alpha_functional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T alpha;

};

template<typename T>
class A_lt_alpha_functional3D : public BoxProcessingFunctional3D_SS<T,int>
{
public:
    A_lt_alpha_functional3D(T alpha_);
    virtual void process(Box3D domain, ScalarField3D<T>& A,
                                       ScalarField3D<int>& result);
    virtual A_lt_alpha_functional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T alpha;
};

template<typename T>
class A_gt_alpha_functional3D : public BoxProcessingFunctional3D_SS<T,int>
{
public:
    A_gt_alpha_functional3D(T alpha_);
    virtual void process(Box3D domain, ScalarField3D<T>& A,
                                       ScalarField3D<int>& result);
    virtual A_gt_alpha_functional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T alpha;
};

template<typename T>
class A_minus_alpha_functional3D : public BoxProcessingFunctional3D_SS<T,T>
{
public:
    A_minus_alpha_functional3D(T alpha_);
    virtual void process(Box3D domain, ScalarField3D<T>& A,
                                       ScalarField3D<T>& result);
    virtual A_minus_alpha_functional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T alpha;
};

template<typename T>
class Alpha_minus_A_functional3D : public BoxProcessingFunctional3D_SS<T,T>
{
public:
    Alpha_minus_A_functional3D(T alpha_);
    virtual void process(Box3D domain, ScalarField3D<T>& A,
                                       ScalarField3D<T>& result);
    virtual Alpha_minus_A_functional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T alpha;
};

template<typename T>
class A_times_alpha_functional3D : public BoxProcessingFunctional3D_SS<T,T>
{
public:
    A_times_alpha_functional3D(T alpha_);
    virtual void process(Box3D domain, ScalarField3D<T>& A,
                                       ScalarField3D<T>& result);
    virtual A_times_alpha_functional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T alpha;
};

template<typename T>
class A_dividedBy_alpha_functional3D : public BoxProcessingFunctional3D_SS<T,T>
{
public:
    A_dividedBy_alpha_functional3D(T alpha_);
    virtual void process(Box3D domain, ScalarField3D<T>& A,
                                       ScalarField3D<T>& result);
    virtual A_dividedBy_alpha_functional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T alpha;
};

template<typename T>
class Alpha_dividedBy_A_functional3D : public BoxProcessingFunctional3D_SS<T,T>
{
public:
    Alpha_dividedBy_A_functional3D(T alpha_);
    virtual void process(Box3D domain, ScalarField3D<T>& A,
                                       ScalarField3D<T>& result);
    virtual Alpha_dividedBy_A_functional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T alpha;
};

template<typename T>
class A_plus_alpha_inplace_functional3D : public BoxProcessingFunctional3D_S<T>
{
public:
    A_plus_alpha_inplace_functional3D(T alpha_);
    virtual void process(Box3D domain, ScalarField3D<T>& A);
    virtual A_plus_alpha_inplace_functional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T alpha;
};

template<typename T>
class A_minus_alpha_inplace_functional3D : public BoxProcessingFunctional3D_S<T>
{
public:
    A_minus_alpha_inplace_functional3D(T alpha_);
    virtual void process(Box3D domain, ScalarField3D<T>& A);
    virtual A_minus_alpha_inplace_functional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T alpha;
};

template<typename T>
class A_times_alpha_inplace_functional3D : public BoxProcessingFunctional3D_S<T>
{
public:
    A_times_alpha_inplace_functional3D(T alpha_);
    virtual void process(Box3D domain, ScalarField3D<T>& A);
    virtual A_times_alpha_inplace_functional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T alpha;
};

template<typename T>
class A_dividedBy_alpha_inplace_functional3D : public BoxProcessingFunctional3D_S<T>
{
public:
    A_dividedBy_alpha_inplace_functional3D(T alpha_);
    virtual void process(Box3D domain, ScalarField3D<T>& A);
    virtual A_dividedBy_alpha_inplace_functional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T alpha;
};


template<typename T>
class A_lt_B_functional3D : public BoxProcessingFunctional3D
{
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> scalarFields);
    virtual A_lt_B_functional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_gt_B_functional3D : public BoxProcessingFunctional3D
{
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> scalarFields);
    virtual A_gt_B_functional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_plus_B_functional3D : public ScalarFieldBoxProcessingFunctional3D<T>
{
public:
    virtual void process(Box3D domain, std::vector<ScalarField3D<T>*> scalarFields);
    virtual A_plus_B_functional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_minus_B_functional3D : public ScalarFieldBoxProcessingFunctional3D<T>
{
public:
    virtual void process(Box3D domain, std::vector<ScalarField3D<T>*> scalarFields);
    virtual A_minus_B_functional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_times_B_functional3D : public ScalarFieldBoxProcessingFunctional3D<T>
{
public:
    virtual void process(Box3D domain, std::vector<ScalarField3D<T>*> scalarFields);
    virtual A_times_B_functional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_dividedBy_B_functional3D : public ScalarFieldBoxProcessingFunctional3D<T>
{
public:
    virtual void process(Box3D domain, std::vector<ScalarField3D<T>*> scalarFields);
    virtual A_dividedBy_B_functional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_plus_B_inplace_functional3D : public BoxProcessingFunctional3D_SS<T,T>
{
public:
    virtual void process(Box3D domain, ScalarField3D<T>& A, ScalarField3D<T>& B);
    virtual A_plus_B_inplace_functional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_minus_B_inplace_functional3D : public BoxProcessingFunctional3D_SS<T,T>
{
public:
    virtual void process(Box3D domain, ScalarField3D<T>& A, ScalarField3D<T>& B);
    virtual A_minus_B_inplace_functional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_times_B_inplace_functional3D : public BoxProcessingFunctional3D_SS<T,T>
{
public:
    virtual void process(Box3D domain, ScalarField3D<T>& A, ScalarField3D<T>& B);
    virtual A_times_B_inplace_functional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class A_dividedBy_B_inplace_functional3D : public BoxProcessingFunctional3D_SS<T,T>
{
public:
    virtual void process(Box3D domain, ScalarField3D<T>& A, ScalarField3D<T>& B);
    virtual A_dividedBy_B_inplace_functional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class UniformlyBoundScalarField3D : public BoxProcessingFunctional3D_S<T>
{
public:
    UniformlyBoundScalarField3D(T bound_)
        : bound(bound_)
    { }
    virtual void process(Box3D domain, ScalarField3D<T>& data);
    virtual UniformlyBoundScalarField3D<T>* clone() const {
        return new UniformlyBoundScalarField3D(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
    }
private:
    T bound;
};

template<typename T, template<typename U> class Descriptor>
class LBMsmoothen3D : public BoxProcessingFunctional3D_SS<T,T> {
public:
    virtual void process(Box3D domain, ScalarField3D<T>& data, ScalarField3D<T>& result);
    virtual LBMsmoothen3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
};

template<typename T>
class Smoothen3D : public BoxProcessingFunctional3D_SS<T,T> {
public:
    virtual void process(Box3D domain, ScalarField3D<T>& data, ScalarField3D<T>& result);
    virtual Smoothen3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
};

template<typename T>
class MollifyScalar3D : public BoxProcessingFunctional3D {
public:
    MollifyScalar3D(T l_, plint d_, Box3D globalDomain_, int exclusionFlag_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual MollifyScalar3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    T l;                // Characteristic length (compact support) of the mollifier.
    plint d;            // Mollifying neighborhood (d <= envelope).
    Box3D globalDomain; // The globalDomain must be at most as big as the whole simulation
                        // domain for non-periodic problems, and bigger than the whole simulation
                        // domain plus the envelope (per periodic direction) for periodic problems.
    int exclusionFlag;  // Flag which is excluded from the mollification proceedure.
};

template<typename T,template<typename U> class Descriptor>
class LBMcomputeGradient3D : public BoxProcessingFunctional3D_ST<T,T,3> {
public:
    virtual void process(Box3D domain, ScalarField3D<T>& scalarField, TensorField3D<T,3>& gradient);
    virtual LBMcomputeGradient3D<T,Descriptor>* clone() const {
        return new LBMcomputeGradient3D(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
        modified[1] = modif::staticVariables;
    }
};

template<typename T>
class UpdateMinScalarTransientStatistics3D : public BoxProcessingFunctional3D {
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks);
    virtual UpdateMinScalarTransientStatistics3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class UpdateMaxScalarTransientStatistics3D : public BoxProcessingFunctional3D {
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks);
    virtual UpdateMaxScalarTransientStatistics3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class UpdateAveScalarTransientStatistics3D : public BoxProcessingFunctional3D {
public:
    UpdateAveScalarTransientStatistics3D(plint n_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks);
    virtual UpdateAveScalarTransientStatistics3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    plint n;
};

template<typename T>
class UpdateRmsScalarTransientStatistics3D : public BoxProcessingFunctional3D {
public:
    UpdateRmsScalarTransientStatistics3D(plint n_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks);
    virtual UpdateRmsScalarTransientStatistics3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    plint n;
};

template<typename T>
class UpdateDevScalarTransientStatistics3D : public BoxProcessingFunctional3D {
public:
    UpdateDevScalarTransientStatistics3D(plint n_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks);
    virtual UpdateDevScalarTransientStatistics3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    plint n;
};


/* *************** PART III ****************************************** */
/* *************** Analysis of the tensor-field ********************** */
/* ******************************************************************* */

template<typename T, int nDim>
class ComputeTensorSqrtFunctional3D : public BoxProcessingFunctional3D_TT<T,nDim,T,nDim>
{
public:
    virtual void process(Box3D domain, TensorField3D<T,nDim>& A,
                                       TensorField3D<T,nDim>& B);
    virtual ComputeTensorSqrtFunctional3D<T,nDim>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim, class BoolMask> 
class CountTensorElementsFunctional3D :
    public ReductiveBoxProcessingFunctional3D_T<T,nDim>
{
public:
    CountTensorElementsFunctional3D(BoolMask boolMask_);
    virtual void process(Box3D domain, TensorField3D<T,nDim>& field);
    virtual CountTensorElementsFunctional3D<T,nDim,BoolMask>* clone() const;
    plint getCount() const;
private:
    plint countId;
    BoolMask boolMask;
};

/** Attention: No matter what the type of T is (even if it is an integer type),
 *    the sum is computed in double-precision floating point numbers, and
 *    converted to T at the end (and rounded, if T is an integer).
 **/
template<typename T,int nDim>
class BoxTensorSumFunctional3D : public ReductiveBoxProcessingFunctional3D_T<T,nDim> {
public:
    BoxTensorSumFunctional3D();
    virtual void process(Box3D domain, TensorField3D<T,nDim>& tensorField);
    virtual BoxTensorSumFunctional3D<T,nDim>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const
    {
        modified[0] = modif::nothing;
    }
    Array<T,nDim> getSumTensor() const;
private:
    Array<plint,nDim> sumTensorId;
};

/** Attention: No matter what the type of T is (even if it is an integer type),
 *    the average is computed in double-precision floating point numbers, and
 *    converted to T at the end (and rounded, if T is an integer).
 **/
template<typename T, int nDim>
class MaskedBoxTensorAverageFunctional3D : public ReductiveBoxProcessingFunctional3D_ST<int,T,nDim> {
public:
    MaskedBoxTensorAverageFunctional3D(int flag_);
    virtual void process( Box3D domain,
                          ScalarField3D<int>& mask,
                          TensorField3D<T,nDim>& tensorField );
    virtual MaskedBoxTensorAverageFunctional3D<T,nDim>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const
    {
        modified[0] = modif::nothing;
        modified[1] = modif::nothing;
    }
    Array<T,nDim> getAverageTensor() const;
private:
    Array<plint, nDim> averageTensorId;
    int flag;
};

template<typename T1, typename T2, int nDim>
class CopyConvertTensorFunctional3D : public BoxProcessingFunctional3D_TT<T1,nDim,T2,nDim>
{
public:
    virtual void process(Box3D domain, TensorField3D<T1,nDim>& field1, TensorField3D<T2,nDim>& field2);
    virtual CopyConvertTensorFunctional3D<T1,T2,nDim>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class ExtractTensorSubDomainFunctional3D :
    public BoxProcessingFunctional3D_TT<T,nDim,T,nDim>
{
public:
    virtual void process(Box3D domain, TensorField3D<T,nDim>& field1,
                                       TensorField3D<T,nDim>& field2);
    virtual ExtractTensorSubDomainFunctional3D<T,nDim>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class ExtractTensorComponentFunctional3D :
    public BoxProcessingFunctional3D_ST<T,T,nDim>
{
public:
    ExtractTensorComponentFunctional3D(int iComponent_);
    virtual void process(Box3D domain, ScalarField3D<T>& scalarField,
                                       TensorField3D<T,nDim>& tensorField );
    virtual ExtractTensorComponentFunctional3D<T,nDim>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    int iComponent;
};

template<typename T, int nDim>
class ComputeNormFunctional3D : public BoxProcessingFunctional3D_ST<T,T,nDim>
{
public:
    virtual void process(Box3D domain, ScalarField3D<T>& scalarField,
                                       TensorField3D<T,nDim>& tensorField);
    virtual ComputeNormFunctional3D<T,nDim>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class ComputeNormSqrFunctional3D : public BoxProcessingFunctional3D_ST<T,T,nDim>
{
public:
    virtual void process(Box3D domain, ScalarField3D<T>& scalarField,
                                       TensorField3D<T,nDim>& tensorField);
    virtual ComputeNormSqrFunctional3D<T,nDim>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class ComputeSymmetricTensorNormFunctional3D :
    public BoxProcessingFunctional3D_ST<T,T,6>
{
public:
    virtual void process(Box3D domain, ScalarField3D<T>& scalarField,
                                       TensorField3D<T,6>& tensorField);
    virtual ComputeSymmetricTensorNormFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class ComputeSymmetricTensorNormSqrFunctional3D :
    public BoxProcessingFunctional3D_ST<T,T,6>
{
public:
    virtual void process(Box3D domain, ScalarField3D<T>& scalarField,
                                       TensorField3D<T,6>& tensorField);
    virtual ComputeSymmetricTensorNormSqrFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class ComputeSymmetricTensorTraceFunctional3D :
    public BoxProcessingFunctional3D_ST<T,T,6>
{
public:
    virtual void process(Box3D domain, ScalarField3D<T>& scalarField,
                                       TensorField3D<T,6>& tensorField);
    virtual ComputeSymmetricTensorTraceFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class BoxBulkGradientFunctional3D :
    public BoxProcessingFunctional3D_ST<T,T,3>
{
public:
    virtual void process(Box3D domain, ScalarField3D<T>& phi,
                                       TensorField3D<T,3>& gradient);
    virtual BoxBulkGradientFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class BoxGradientFunctional3D :
    public BoundedBoxProcessingFunctional3D_ST<T,T,3>
{
public:
    virtual void processBulk( Box3D domain, ScalarField3D<T>& phi,
                                            TensorField3D<T,3>& gradient );
    virtual void processPlane( int direction, int orientation, Box3D domain,
                               ScalarField3D<T>& phi,
                               TensorField3D<T,3>& gradient );
    virtual void processEdge( int plane, int normal1, int normal2, Box3D domain,
                              ScalarField3D<T>& phi,
                              TensorField3D<T,3>& gradient );
    virtual void processCorner( int normalX, int normalY, int normalZ, Box3D domain,
                                ScalarField3D<T>& phi,
                                TensorField3D<T,3>& gradient );
    virtual BoxGradientFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class BoxBulkVorticityFunctional3D :
    public BoxProcessingFunctional3D_TT<T,nDim,T,nDim>
{
public:
    virtual void process(Box3D domain, TensorField3D<T,nDim>& velocity,
                                       TensorField3D<T,nDim>& vorticity);
    virtual BoxBulkVorticityFunctional3D<T,nDim>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class BoxVorticityFunctional3D :
    public BoundedBoxProcessingFunctional3D_TT<T,nDim,T,nDim>
{
public:
    virtual void processBulk( Box3D domain, TensorField3D<T,nDim>& velocity,
                                            TensorField3D<T,nDim>& vorticity );
    virtual void processPlane( int direction, int orientation, Box3D domain,
                               TensorField3D<T,nDim>& velocity,
                               TensorField3D<T,nDim>& vorticity );
    virtual void processEdge( int plane, int normal1, int normal2, Box3D domain,
                              TensorField3D<T,nDim>& velocity,
                              TensorField3D<T,nDim>& vorticity );
    virtual void processCorner( int normalX, int normalY, int normalZ, Box3D domain,
                                TensorField3D<T,nDim>& velocity,
                                TensorField3D<T,nDim>& vorticity );
    virtual BoxVorticityFunctional3D<T,nDim>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class BoxBulkStrainRateFunctional3D :
    public BoxProcessingFunctional3D_TT<T,nDim,T,SymmetricTensorImpl<T,nDim>::n>
{
public:
    virtual void process(Box3D domain, TensorField3D<T,nDim>& velocity,
                                       TensorField3D<T,SymmetricTensorImpl<T,nDim>::n>& S);
    virtual BoxBulkStrainRateFunctional3D<T,nDim>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class BoxStrainRateFunctional3D :
    public BoundedBoxProcessingFunctional3D_TT<T,nDim,T,SymmetricTensorImpl<T,nDim>::n>
{
public:
    virtual void processBulk( Box3D domain, TensorField3D<T,nDim>& velocity,
                                            TensorField3D<T,SymmetricTensorImpl<T,nDim>::n>& S );
    virtual void processPlane( int direction, int orientation, Box3D domain,
                               TensorField3D<T,nDim>& velocity,
                               TensorField3D<T,SymmetricTensorImpl<T,nDim>::n>& S );
    virtual void processEdge( int plane, int normal1, int normal2, Box3D domain,
                              TensorField3D<T,nDim>& velocity,
                              TensorField3D<T,SymmetricTensorImpl<T,nDim>::n>& S );
    virtual void processCorner( int normalX, int normalY, int normalZ, Box3D domain,
                                TensorField3D<T,nDim>& velocity,
                                TensorField3D<T,SymmetricTensorImpl<T,nDim>::n>& S );
    virtual BoxStrainRateFunctional3D<T,nDim>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class BoxQcriterionFunctional3D : public BoxProcessingFunctional3D
{
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual BoxQcriterionFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class BoxComputeInstantaneousReynoldsStressFunctional3D : public BoxProcessingFunctional3D
{
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual BoxComputeInstantaneousReynoldsStressFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T>
class BoxLambda2Functional3D : public BoxProcessingFunctional3D
{
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual BoxLambda2Functional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class BoxBulkDivergenceFunctional3D :
    public BoxProcessingFunctional3D_ST<T,T,nDim>
{
public:
    virtual void process(Box3D domain, ScalarField3D<T>& divergence,
                                       TensorField3D<T,nDim>& velocity);
    virtual BoxBulkDivergenceFunctional3D<T,nDim>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class Tensor_A_plus_B_functional3D : public TensorFieldBoxProcessingFunctional3D<T,nDim>
{
public:
    virtual void process(Box3D domain, std::vector<TensorField3D<T,nDim>*> tensorFields);
    virtual Tensor_A_plus_B_functional3D<T,nDim>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class Tensor_A_minus_B_functional3D : public TensorFieldBoxProcessingFunctional3D<T,nDim>
{
public:
    virtual void process(Box3D domain, std::vector<TensorField3D<T,nDim>*> tensorFields);
    virtual Tensor_A_minus_B_functional3D<T,nDim>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class Tensor_A_times_B_functional3D : public TensorFieldBoxProcessingFunctional3D<T,nDim>
{
public:
    virtual void process(Box3D domain, std::vector<TensorField3D<T,nDim>*> tensorFields);
    virtual Tensor_A_times_B_functional3D<T,nDim>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class TensorProduct_A_A_functional3D : public BoxProcessingFunctional3D
{
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> tensorFields);
    virtual TensorProduct_A_A_functional3D<T,nDim>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class IndexContraction_SymmetricTensor_A_SymmetricTensor_B_functional3D : public BoxProcessingFunctional3D
{
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> tensorFields);
    virtual IndexContraction_SymmetricTensor_A_SymmetricTensor_B_functional3D<T,nDim>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};


template<typename T, int nDim>
class Scalar_A_times_Tensor_B_functional3D : public BoxProcessingFunctional3D
{
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> tensorFields);
    virtual Scalar_A_times_Tensor_B_functional3D<T,nDim>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};


template<typename T, int nDim>
class Tensor_A_dividedBy_B_functional3D : public TensorFieldBoxProcessingFunctional3D<T,nDim>
{
public:
    virtual void process(Box3D domain, std::vector<TensorField3D<T,nDim>*> tensorFields);
    virtual Tensor_A_dividedBy_B_functional3D<T,nDim>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class Tensor_A_plus_B_inplace_functional3D :
    public BoxProcessingFunctional3D_TT<T,nDim,T,nDim>
{
public:
    virtual void process(Box3D domain, TensorField3D<T,nDim>& A,
                                       TensorField3D<T,nDim>& B);
    virtual Tensor_A_plus_B_inplace_functional3D<T,nDim>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class Tensor_A_minus_B_inplace_functional3D :
    public BoxProcessingFunctional3D_TT<T,nDim,T,nDim>
{
public:
    virtual void process(Box3D domain, TensorField3D<T,nDim>& A,
                                       TensorField3D<T,nDim>& B);
    virtual Tensor_A_minus_B_inplace_functional3D<T,nDim>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class Tensor_A_times_B_inplace_functional3D :
    public BoxProcessingFunctional3D_TT<T,nDim,T,nDim>
{
public:
    virtual void process(Box3D domain, TensorField3D<T,nDim>& A,
                                       TensorField3D<T,nDim>& B);
    virtual Tensor_A_times_B_inplace_functional3D<T,nDim>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class Tensor_A_times_alpha_functional3D : public BoxProcessingFunctional3D_TT<T,nDim,T,nDim>
{
public:
    Tensor_A_times_alpha_functional3D(T alpha_);
    virtual void process(Box3D domain, TensorField3D<T,nDim>& A,
                                       TensorField3D<T,nDim>& result);
    virtual Tensor_A_times_alpha_functional3D<T,nDim>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T alpha;
};

template<typename T, int nDim>
class Tensor_A_times_alpha_inplace_functional3D :
    public BoxProcessingFunctional3D_T<T,nDim>
{
public:
    Tensor_A_times_alpha_inplace_functional3D(T alpha_);
    virtual void process(Box3D domain, TensorField3D<T,nDim>& A);
    virtual Tensor_A_times_alpha_inplace_functional3D<T,nDim>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    T alpha;
};

template<typename T, int nDim>
class Tensor_A_dividedBy_B_inplace_functional3D :
    public BoxProcessingFunctional3D_TT<T,nDim,T,nDim>
{
public:
    virtual void process(Box3D domain, TensorField3D<T,nDim>& A,
                                       TensorField3D<T,nDim>& B);
    virtual Tensor_A_dividedBy_B_inplace_functional3D<T,nDim>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

template<typename T, int nDim>
class Normalize_Tensor_functional3D : public BoxProcessingFunctional3D_TT<T,nDim,T,nDim>
{
public:
    Normalize_Tensor_functional3D(Precision precision_ = DBL)
        : precision(precision_)
    {
        PLB_ASSERT(precision == FLT || precision == DBL || precision == LDBL || precision == INF);
    }
    virtual void process(Box3D domain, TensorField3D<T,nDim>& data, TensorField3D<T,nDim>& result);
    virtual Normalize_Tensor_functional3D<T,nDim>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    Precision precision;
};

template<typename T, int nDim>
class Normalize_Tensor_inplace_functional3D :
    public BoxProcessingFunctional3D_T<T,nDim>
{
public:
    Normalize_Tensor_inplace_functional3D(Precision precision_ = DBL)
        : precision(precision_)
    {
        PLB_ASSERT(precision == FLT || precision == DBL || precision == LDBL || precision == INF);
    }
    virtual void process(Box3D domain, TensorField3D<T,nDim>& data);
    virtual Normalize_Tensor_inplace_functional3D<T,nDim>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    Precision precision;
};

template<typename T, int nDim, template<typename U> class Descriptor>
class LBMsmoothenTensor3D : public BoxProcessingFunctional3D_TT<T,nDim,T,nDim> {
public:
    virtual void process(Box3D domain, TensorField3D<T,nDim>& data, TensorField3D<T,nDim>& result);
    virtual LBMsmoothenTensor3D<T,nDim,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
};

template<typename T, int nDim>
class SmoothenTensor3D : public BoxProcessingFunctional3D_TT<T,nDim,T,nDim> {
public:
    virtual void process(Box3D domain, TensorField3D<T,nDim>& data, TensorField3D<T,nDim>& result);
    virtual SmoothenTensor3D<T,nDim>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
};

template<typename T, int nDim>
class MollifyTensor3D : public BoxProcessingFunctional3D {
public:
    MollifyTensor3D(T l_, plint d_, Box3D globalDomain_, int exclusionFlag_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> fields);
    virtual MollifyTensor3D<T,nDim>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
private:
    T l;                // Characteristic length (compact support) of the mollifier.
    plint d;            // Mollifying neighborhood (d <= envelope).
    Box3D globalDomain; // The globalDomain must be at most as big as the whole simulation
                        // domain for non-periodic problems, and bigger than the whole simulation
                        // domain plus the envelope (per periodic direction) for periodic problems.
    int exclusionFlag;  // Flag which is excluded from the mollification proceedure.
};

template<typename T,template<typename U> class Descriptor>
class LBMcomputeDivergence3D : public BoxProcessingFunctional3D_ST<T,T,3> {
public:
    virtual void process(Box3D domain, ScalarField3D<T>& divergence, TensorField3D<T,3>& vectorField);
    virtual LBMcomputeDivergence3D<T,Descriptor>* clone() const {
        return new LBMcomputeDivergence3D(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::staticVariables;
        modified[1] = modif::nothing;
    }
};

}  // namespace plb

#endif  // DATA_ANALYSIS_FUNCTIONAL_3D_H

// Explicitly include generic algorithms which are never precompiled (not even in precompiled version)
#include "dataProcessors/dataAnalysisGenerics3D.h"

