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
 * Helper functions for data field initialization -- generic implementation.
 */
#ifndef DATA_FIELD_INITIALIZER_GENERICS_3D_H
#define DATA_FIELD_INITIALIZER_GENERICS_3D_H

#include "core/globalDefs.h"
#include "atomicBlock/dataProcessorWrapper3D.h"
#include "multiBlock/multiDataProcessorWrapper3D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "multiGrid/multiGridUtil.h"

namespace plb {

/* *************** PART I ******************************************** */
/* *************** Initialization of the block-lattice *************** */
/* ******************************************************************* */

/* ************ Class SetCustomBoundaryVelocityFunctional3D ********** */

template<typename T, template<typename U> class Descriptor, class VelocityFunction>
SetCustomBoundaryVelocityFunctional3D<T,Descriptor,VelocityFunction>::
    SetCustomBoundaryVelocityFunctional3D(VelocityFunction f_)
        : f(f_),
          velocityScale( (T)1 )
{ }

template<typename T, template<typename U> class Descriptor, class VelocityFunction>
SetCustomBoundaryVelocityFunctional3D<T,Descriptor,VelocityFunction>*
    SetCustomBoundaryVelocityFunctional3D<T,Descriptor,VelocityFunction>::clone() const
{
    return new SetCustomBoundaryVelocityFunctional3D<T,Descriptor,VelocityFunction>(*this);
}

template<typename T, template<typename U> class Descriptor, class VelocityFunction>
void SetCustomBoundaryVelocityFunctional3D<T,Descriptor,VelocityFunction>::execute (
        plint iX, plint iY, plint iZ, Cell<T,Descriptor>& cell ) const
{
    Array<T,Descriptor<T>::d> u;
    f(iX, iY, iZ, u);
    u[0] *= velocityScale;
    u[1] *= velocityScale;
    u[2] *= velocityScale;
    cell.defineVelocity(u);
}

template<typename T, template<typename U> class Descriptor, class VelocityFunction>
void SetCustomBoundaryVelocityFunctional3D<T,Descriptor,VelocityFunction>::setscale (
        int dxScale, int dtScale )
{
    int dimDx = 1;
    int dimDt = -1;
    velocityScale = scaleFromReference(dxScale, dimDx,
                                       dtScale, dimDt);
}

/* ************ Class SetCustomBoundaryDensityFunctional3D ********** */

template<typename T, template<typename U> class Descriptor, class DensityFunction>
SetCustomBoundaryDensityFunctional3D<T,Descriptor,DensityFunction>::
    SetCustomBoundaryDensityFunctional3D(DensityFunction f_)
        : f(f_)
{ }

template<typename T, template<typename U> class Descriptor, class DensityFunction>
void SetCustomBoundaryDensityFunctional3D<T,Descriptor,DensityFunction>::
    execute(plint iX, plint iY, plint iZ, Cell<T,Descriptor>& cell) const
{
    // No rescaling needed: rho is scale invariant.
    T rho = f(iX, iY, iZ);
    cell.defineDensity(rho);
}

template<typename T, template<typename U> class Descriptor, class DensityFunction>
SetCustomBoundaryDensityFunctional3D<T,Descriptor,DensityFunction>*
    SetCustomBoundaryDensityFunctional3D<T,Descriptor,DensityFunction>::clone() const
{
    return new SetCustomBoundaryDensityFunctional3D<T,Descriptor,DensityFunction>(*this);
}


/* ************ Class IniCustomEquilibriumFunctional3D ********** */

template<typename T, template<typename U> class Descriptor, class RhoUFunction>
IniCustomEquilibriumFunctional3D<T,Descriptor,RhoUFunction>::
    IniCustomEquilibriumFunctional3D(RhoUFunction f_)
        : f(f_),
          velocityScale( (T)1 )
{ }

template<typename T, template<typename U> class Descriptor, class RhoUFunction>
void IniCustomEquilibriumFunctional3D<T,Descriptor,RhoUFunction>::execute (
        plint iX, plint iY, plint iZ, Cell<T,Descriptor>& cell ) const
{
    Array<T,Descriptor<T>::d> j;
    T rho;
    f(iX, iY, iZ, rho, j);
    for (int iD=0; iD<Descriptor<T>::d; ++iD) {
        j[iD] *= rho;
        j[iD] *= velocityScale;
    }
    T rhoBar = Descriptor<T>::rhoBar(rho);
    T jSqr = VectorTemplate<T,Descriptor>::normSqr(j);
    for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
        cell[iPop] = cell.computeEquilibrium(iPop, rhoBar, j, jSqr);
    }
}

template<typename T, template<typename U> class Descriptor, class RhoUFunction>
IniCustomEquilibriumFunctional3D<T,Descriptor,RhoUFunction>*
    IniCustomEquilibriumFunctional3D<T,Descriptor,RhoUFunction>::clone() const
{
    return new IniCustomEquilibriumFunctional3D<T,Descriptor,RhoUFunction>(*this);
}

template<typename T, template<typename U> class Descriptor, class RhoUFunction>
void IniCustomEquilibriumFunctional3D<T,Descriptor,RhoUFunction>::setscale (
        int dxScale, int dtScale )
{
    int dimDx = 1;
    int dimDt = -1;
    velocityScale = scaleFromReference(dxScale, dimDx,
                                       dtScale, dimDt);
}

/* ************ Class IniCustomThermalEquilibriumFunctional3D ********** */

template<typename T, template<typename U> class Descriptor, class RhoVelTempFunction>
IniCustomThermalEquilibriumFunctional3D<T,Descriptor,RhoVelTempFunction>::
    IniCustomThermalEquilibriumFunctional3D(RhoVelTempFunction f_)
        : f(f_),
          velocityScale( (T)1 )
{ }

template<typename T, template<typename U> class Descriptor, class RhoVelTempFunction>
void IniCustomThermalEquilibriumFunctional3D<T,Descriptor,RhoVelTempFunction>::execute (
        plint iX, plint iY, plint iZ, Cell<T,Descriptor>& cell ) const
{
    Array<T,Descriptor<T>::d> j;
    T rho;
    T temperature;
    f(iX, iY, iZ, rho, j, temperature);
    for (int iD=0; iD<Descriptor<T>::d; ++iD) {
        j[iD] *= rho;
        j[iD] *= velocityScale;
    }
    T rhoBar = Descriptor<T>::rhoBar(rho);
    T thetaBar = temperature-(T)1;
    T jSqr = VectorTemplate<T,Descriptor>::normSqr(j);
    for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
        cell[iPop] = cell.computeEquilibrium(iPop, rhoBar, j, jSqr, thetaBar);
    }
}

template<typename T, template<typename U> class Descriptor, class RhoVelTempFunction>
IniCustomThermalEquilibriumFunctional3D<T,Descriptor,RhoVelTempFunction>*
    IniCustomThermalEquilibriumFunctional3D<T,Descriptor,RhoVelTempFunction>::clone() const
{
    return new IniCustomThermalEquilibriumFunctional3D<T,Descriptor,RhoVelTempFunction>(*this);
}

template<typename T, template<typename U> class Descriptor, class RhoVelTempFunction>
void IniCustomThermalEquilibriumFunctional3D<T,Descriptor,RhoVelTempFunction>::setscale (
        int dxScale, int dtScale )
{
    int dimDx = 1;
    int dimDt = -1;
    velocityScale = scaleFromReference(dxScale, dimDx,
                                       dtScale, dimDt);
}


/* ************* Class IniConstEquilibriumOnDomainFunctional3D ******************* */

template<typename T, template<typename U> class Descriptor, class DomainFunctional>
IniConstEquilibriumOnDomainFunctional3D<T,Descriptor,DomainFunctional>::IniConstEquilibriumOnDomainFunctional3D (
        T density_, Array<T,Descriptor<T>::d> velocity, T temperature,
        DomainFunctional const& domain_)
    : rhoBar(Descriptor<T>::rhoBar(density_)),
      j     (density_*velocity[0], density_*velocity[1], density_*velocity[2]),
      jSqr  (VectorTemplate<T,Descriptor>::normSqr(j)),
      thetaBar(temperature-(T)1),
      domain(domain_)
{ }

template<typename T, template<typename U> class Descriptor, class DomainFunctional>
void IniConstEquilibriumOnDomainFunctional3D<T,Descriptor,DomainFunctional>::process (
        Box3D bbox, BlockLattice3D<T,Descriptor>& lattice )
{
    int dimDx = 1;
    int dimDt = -1;
    Dot3D absPos = lattice.getLocation();
    T scaleFactor = scaleFromReference(this->getDxScale(), dimDx,
                                       this->getDtScale(), dimDt);
    Array<T,3> scaledJ = j*scaleFactor;
    T scaledJsqr = jSqr*scaleFactor*scaleFactor;
    for (plint iX=bbox.x0; iX<=bbox.x1; ++iX) {
        for (plint iY=bbox.y0; iY<=bbox.y1; ++iY) {
            for (plint iZ=bbox.z0; iZ<=bbox.z1; ++iZ) {
                if (domain(iX+absPos.x,iY+absPos.y,iZ+absPos.z)) {
                    for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
                        lattice.get(iX,iY,iZ)[iPop] =
                            lattice.get(iX,iY,iZ).computeEquilibrium(iPop, rhoBar, scaledJ, scaledJsqr, thetaBar);
                    }
                }
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor, class DomainFunctional>
IniConstEquilibriumOnDomainFunctional3D<T,Descriptor,DomainFunctional>*
    IniConstEquilibriumOnDomainFunctional3D<T,Descriptor,DomainFunctional>::clone() const
{
    return new IniConstEquilibriumOnDomainFunctional3D<T,Descriptor,DomainFunctional>(*this);
}

template<typename T, template<typename U> class Descriptor, class DomainFunctional>
void IniConstEquilibriumOnDomainFunctional3D<T,Descriptor,DomainFunctional>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
}



/* *************** PART II ******************************************* */
/* *************** Initialization of the scalar- and tensor-field **** */
/* ******************************************************************* */

/* ************ SetToScalarFunctionFunctional3D ********************** */

template<typename T, class Function>
SetToScalarFunctionFunctional3D<T,Function>::SetToScalarFunctionFunctional3D(Function f_)
    : f(f_)
{ }

template<typename T, class Function>
void SetToScalarFunctionFunctional3D<T,Function>::process (
        Box3D domain, ScalarField3D<T>& field )
{
    Dot3D relativeOffset = field.getLocation();
    Array<plint,3> ofs(relativeOffset.x, relativeOffset.y, relativeOffset.z);
    Array<plint,3> pos;
    for ( pos[0]=domain.x0; pos[0]<=domain.x1; ++pos[0] ) {
        for ( pos[1]=domain.y0; pos[1]<=domain.y1; ++pos[1] ) {
            for ( pos[2]=domain.z0; pos[2]<=domain.z1; ++pos[2] ) {
                    field.get(pos[0],pos[1],pos[2]) = 
                            f(pos[0]+ofs[0], pos[1]+ofs[1], pos[2]+ofs[2]);
            }
        }
    }
}

template<typename T, class Function>
SetToScalarFunctionFunctional3D<T,Function>*
    SetToScalarFunctionFunctional3D<T,Function>::clone() const
{
    return new SetToScalarFunctionFunctional3D<T,Function>(*this);
}

template<typename T, class Function>
BlockDomain::DomainT SetToScalarFunctionFunctional3D<T,Function>::appliesTo() const
{
    // Boundary cannot be included, because periodic boundaries would
    //   get the wrong value.
    return BlockDomain::bulk;
}

template<typename T, class Function>
void SetToScalarFunctionFunctional3D<T,Function>::getTypeOfModification (
        std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::staticVariables;
}


/* ************ AnalyticalSetRhoBarJFunctional3D ********************** */

template<typename T, class Function>
AnalyticalSetRhoBarJFunctional3D<T,Function>::AnalyticalSetRhoBarJFunctional3D(Function const& function_)
    : function(function_)
{ }

template<typename T, class Function>
void AnalyticalSetRhoBarJFunctional3D<T,Function>::process (
        Box3D domain, NTensorField3D<T>& rhoBarJField )
{
    PLB_ASSERT(rhoBarJField.getNdim()==4);
    Dot3D relativeOffset = rhoBarJField.getLocation();
    Array<plint,3> ofs(relativeOffset.x, relativeOffset.y, relativeOffset.z);
    Array<plint,3> pos;
    T rhoBar;
    Array<T,3> j;
    for ( pos[0]=domain.x0; pos[0]<=domain.x1; ++pos[0] ) {
        for ( pos[1]=domain.y0; pos[1]<=domain.y1; ++pos[1] ) {
            for ( pos[2]=domain.z0; pos[2]<=domain.z1; ++pos[2] ) {
                function(pos[0]+ofs[0], pos[1]+ofs[1], pos[2]+ofs[2], rhoBar, j);
                T* rhoBarJ = rhoBarJField.get(pos[0],pos[1],pos[2]);
                *rhoBarJ = rhoBar;
                j.to_cArray(rhoBarJ+1);
            }
        }
    }
}

template<typename T, class Function>
AnalyticalSetRhoBarJFunctional3D<T,Function>*
    AnalyticalSetRhoBarJFunctional3D<T,Function>::clone() const
{
    return new AnalyticalSetRhoBarJFunctional3D<T,Function>(*this);
}

template<typename T, class Function>
void AnalyticalSetRhoBarJFunctional3D<T,Function>::getTypeOfModification (
        std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::staticVariables;
}


/* ************ SetToTensorFunctionFunctional3D ********************** */

template<typename T, int nDim, class Function>
SetToTensorFunctionFunctional3D<T,nDim,Function>::
    SetToTensorFunctionFunctional3D(Function f_)
        : f(f_)
{ }

template<typename T, int nDim, class Function>
void SetToTensorFunctionFunctional3D<T,nDim,Function>::process (
        Box3D domain, TensorField3D<T,nDim>& field )
{
    Dot3D relativeOffset = field.getLocation();
    Array<plint,nDim> ofs(relativeOffset.x, relativeOffset.y, relativeOffset.z);
    Array<plint,nDim> pos;
    Array<T,nDim> value;
    for ( pos[0]=domain.x0; pos[0]<=domain.x1; ++pos[0] ) {
        for ( pos[1]=domain.y0; pos[1]<=domain.y1; ++pos[1] ) {
            for ( pos[2]=domain.z0; pos[2]<=domain.z1; ++pos[2] ) {
                f(pos[0]+ofs[0], pos[1]+ofs[1], pos[2]+ofs[2], value);
                field.get(pos[0],pos[1],pos[2]) = value;
            }
        }
    }
}

template<typename T, int nDim, class Function>
SetToTensorFunctionFunctional3D<T,nDim,Function>*
    SetToTensorFunctionFunctional3D<T,nDim,Function>::clone() const
{
    return new SetToTensorFunctionFunctional3D<T,nDim,Function>(*this);
}

template<typename T, int nDim, class Function>
BlockDomain::DomainT SetToTensorFunctionFunctional3D<T,nDim,Function>::appliesTo() const
{
    // Boundary cannot be included, because periodic boundaries
    //   would get the wrong value.
    return BlockDomain::bulk;
}

template<typename T, int nDim, class Function>
void SetToTensorFunctionFunctional3D<T,nDim,Function>::getTypeOfModification (
        std::vector<modif::ModifT>& modified ) const
{
    modified[0] = modif::staticVariables;
}

}  // namespace plb

#endif  // DATA_FIELD_INITIALIZER_GENERICS_3D_H
