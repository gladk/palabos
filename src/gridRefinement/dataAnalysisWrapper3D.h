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
 * Coupling between grids of different refinement level -- header file.
 */

#ifndef ML_DATA_ANALYSIS_WRAPPERS_3D_H
#define ML_DATA_ANALYSIS_WRAPPERS_3D_H

#include "core/globalDefs.h"
#include "gridRefinement/couplingInterfaceGenerator3D.h"
#include "gridRefinement/couplingActionsGenerator3D.h"
#include "gridRefinement/multiLevelTensorField3D.h"
#include "gridRefinement/multiLevelScalarField3D.h"


#include <map>
#include <vector>

namespace plb {

/* ******************************************************************* */
/* *************** PART I. Multi-block wrappers: Block-Lattice ****** */
/* ******************************************************************* */

/* *************** Density ****************************************** */

template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine>
void computeDensity(
    MultiLevelActions3D<T,Descriptor,Engine>& lattices, 
    MultiLevelScalarField3D<T>& densities, 
    Box3D domain, plint levelOfDomain);

template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine>
std::auto_ptr<MultiLevelScalarField3D<T> >
   computeDensity(MultiLevelActions3D<T,Descriptor,Engine>& lattices, Box3D domain, plint levelOfDomain);

 template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine>
std::auto_ptr<MultiLevelScalarFieldForOutput3D<T> >
   computeDensity(MultiLevelActions3D<T,Descriptor,Engine>& lattices, Box3D domain, plint levelOfDomain, bool crop);

/* *************** Velocity ****************************************** */

template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine>
void computeVelocity(MultiLevelActions3D<T,Descriptor,Engine>& lattices, 
	MultiLevelTensorField3D<T,Descriptor<T>::d >& velocities, 
	Box3D domain, plint levelOfDomain);

template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine>
std::auto_ptr<MultiLevelTensorField3D<T,Descriptor<T>::d> >
   computeVelocity(MultiLevelActions3D<T,Descriptor,Engine>& lattices, Box3D domain, plint levelOfDomain);

template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine>
std::auto_ptr<MultiLevelTensorFieldForOutput3D<T,Descriptor<T>::d> >
   computeVelocity(MultiLevelActions3D<T,Descriptor,Engine>& lattices, Box3D domain, plint levelOfDomain, bool crop);
/* *************** Density ****************************************** */

template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine>
void computeDensity(
    MultiLevelCoupling3D<T,Descriptor,Engine>& lattices, 
    MultiLevelScalarField3D<T>& densities, 
    Box3D domain, plint levelOfDomain);

template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine>
std::auto_ptr<MultiLevelScalarField3D<T> >
   computeDensity(MultiLevelCoupling3D<T,Descriptor,Engine>& lattices, Box3D domain, plint levelOfDomain);

 template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine>
std::auto_ptr<MultiLevelScalarFieldForOutput3D<T> >
   computeDensity(MultiLevelCoupling3D<T,Descriptor,Engine>& lattices, Box3D domain, plint levelOfDomain, bool crop);

/* *************** Velocity ****************************************** */

template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine>
void computeVelocity(MultiLevelCoupling3D<T,Descriptor,Engine>& lattices, 
	MultiLevelTensorField3D<T,Descriptor<T>::d >& velocities, 
	Box3D domain, plint levelOfDomain);

template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine>
std::auto_ptr<MultiLevelTensorField3D<T,Descriptor<T>::d> >
   computeVelocity(MultiLevelCoupling3D<T,Descriptor,Engine>& lattices, Box3D domain, plint levelOfDomain);

template<typename T, template<typename U> class Descriptor,
    template<typename T2, template<typename U2> class Descriptor2> class Engine>
std::auto_ptr<MultiLevelTensorFieldForOutput3D<T,Descriptor<T>::d> >
   computeVelocity(MultiLevelCoupling3D<T,Descriptor,Engine>& lattices, Box3D domain, plint levelOfDomain, bool crop);

/* ******************************************************************* */
/* *************** PART II. Multi-block wrappers: Scalar-Field ******** */
/* ******************************************************************* */

/* *************** copy ****************************************** */

template<typename T>
void copy(MultiLevelScalarField3D<T> &from, 
          const Box3D &fromDomain, plint levelOfFromDomain,
          MultiLevelScalarFieldForOutput3D<T> &to, 
          const Box3D &toDomain, plint levelOfToDomain);

/* *************** exportForOutput ****************************************** */

template<typename T>
std::auto_ptr<MultiLevelScalarFieldForOutput3D<T> > 
	exportForOutput(MultiLevelScalarField3D<T> &from, 
          const Box3D &fromDomain, plint levelOfFromDomain, bool crop);

/* *************** Extract Sub-MultiLevelScalarField *************************** */

template<typename T>
void extractSubDomain( MultiLevelScalarField3D<T> &field,
                       MultiLevelScalarField3D<T> &extractedField,
                       const Box3D &domain, plint levelOfDomain);

template<typename T>
std::auto_ptr<MultiLevelScalarField3D<T> > extractSubDomain(
	MultiLevelScalarField3D<T>& field, 
	const Box3D &domain, plint levelOfDomain);

/* *************** add ****************************************** */

template<typename T>
void add(MultiLevelScalarField3D<T> &A, 
         MultiLevelScalarField3D<T> &B, 
         MultiLevelScalarField3D<T> &aPlusB,
         const Box3D &domain, plint levelOfDomain);

template<typename T>
std::auto_ptr<MultiLevelScalarField3D<T> > 
	add(MultiLevelScalarField3D<T> &A, 
     	MultiLevelScalarField3D<T> &B, 
        const Box3D &domain, plint levelOfDomain);

/* *************** addInPlace ****************************************** */

template<typename T>
void addInPlace(MultiLevelScalarField3D<T> &A, 
                MultiLevelScalarField3D<T> &B, 
                const Box3D &domain, plint levelOfDomain);

/* *************** subtract ****************************************** */

template<typename T>
void subtract(MultiLevelScalarField3D<T> &A, 
              MultiLevelScalarField3D<T> &B, 
           	  MultiLevelScalarField3D<T> &aMinusB,
           	  const Box3D &domain, plint levelOfDomain);

template<typename T>
std::auto_ptr<MultiLevelScalarField3D<T> >  
	subtract(MultiLevelScalarField3D<T> &A, 
             MultiLevelScalarField3D<T> &B, 
           	 const Box3D &domain, plint levelOfDomain);

/* *************** subtractInPlace ****************************************** */

template<typename T>
void subtractInPlace(MultiLevelScalarField3D<T> &A, 
              MultiLevelScalarField3D<T> &B, 
           	  const Box3D &domain, plint levelOfDomain);

/* *************** multiply ****************************************** */

template<typename T>
void multiply(MultiLevelScalarField3D<T> &A, 
              MultiLevelScalarField3D<T> &B, 
           	  MultiLevelScalarField3D<T> &aMinusB,
           	  const Box3D &domain, plint levelOfDomain);

template<typename T>
std::auto_ptr<MultiLevelScalarField3D<T> >   
	multiply(MultiLevelScalarField3D<T> &A, 
             MultiLevelScalarField3D<T> &B, 
           	 const Box3D &domain, plint levelOfDomain);

/* *************** multiplyInPlace ****************************************** */

template<typename T>
void multiplyInPlace(MultiLevelScalarField3D<T> &A, 
                     MultiLevelScalarField3D<T> &B, 
           	         const Box3D &domain, plint levelOfDomain);

/* *************** multiply by scalar ****************************************** */

template<typename T>
void multiply(MultiLevelScalarField3D<T> &A, 
              T alpha, MultiLevelScalarField3D<T> &aTimesB,
           	  const Box3D &domain, plint levelOfDomain);

template<typename T>
std::auto_ptr<MultiLevelScalarField3D<T> >   
	multiply(MultiLevelScalarField3D<T> &A, 
             T alpha,
           	 const Box3D &domain, plint levelOfDomain);

/* *************** multiplyInPlace by scalar ****************************************** */

template<typename T>
void multiplyInPlace(MultiLevelScalarField3D<T> &A, 
                     T alpha, 
           	         const Box3D &domain, plint levelOfDomain);

/* *************** divide ****************************************** */

template<typename T>
void divide(MultiLevelScalarField3D<T> &A, 
            MultiLevelScalarField3D<T> &B, 
           	MultiLevelScalarField3D<T> &aMinusB,
           	const Box3D &domain, plint levelOfDomain);

template<typename T>
std::auto_ptr<MultiLevelScalarField3D<T> >  
	divide(MultiLevelScalarField3D<T> &A, 
           MultiLevelScalarField3D<T> &B, 
           const Box3D &domain, plint levelOfDomain);

/* *************** divideInPlace ****************************************** */

template<typename T>
void divideInPlace(MultiLevelScalarField3D<T> &A, 
                   MultiLevelScalarField3D<T> &B, 
           	       const Box3D &domain, plint levelOfDomain);

/* ******************************************************************* */
/* *************** PART III. Multi-block wrappers: Tensor-field ******* */
/* ******************************************************************* */

/* *************** copy ****************************************** */

template<typename T, int nDim>
void copy(MultiLevelTensorField3D<T,nDim> &from, 
          const Box3D &fromDomain, plint levelOfFromDomain,
          MultiLevelTensorFieldForOutput3D<T,nDim> &to, 
          const Box3D &toDomain, plint levelOfToDomain);

/* *************** exportForOutput ****************************************** */

template<typename T, int nDim>
std::auto_ptr<MultiLevelTensorFieldForOutput3D<T,nDim> > 
	exportForOutput(MultiLevelTensorField3D<T,nDim> &from, 
          const Box3D &fromDomain, plint levelOfFromDomain, bool crop);

/* *************** Extract Sub-TensorField *************************** */

template<typename T, int nDim>
void extractSubDomain( MultiLevelTensorField3D<T,nDim> &field,
                       MultiLevelTensorField3D<T,nDim> &extractedField,
                       const Box3D &domain, plint levelOfDomain);

template<typename T, int nDim>
std::auto_ptr<MultiLevelTensorField3D<T,nDim> > extractSubDomain(
	MultiLevelTensorField3D<T,nDim>& field, 
	const Box3D &domain, plint levelOfDomain);

/* ***** Component (multilevel scalar-field) out of a multilevel-tensor-field *** */

template<typename T, int nDim>
void extractComponent(
	MultiLevelTensorField3D<T,nDim>& tensorField, 
	MultiLevelScalarField3D<T>& component, 
	const Box3D &domain, plint levelOfDomain, int iComponent);

template<typename T, int nDim>
std::auto_ptr<MultiLevelScalarField3D<T> > extractComponent(
	MultiLevelTensorField3D<T,nDim>& tensorField, 
	const Box3D &domain, plint levelOfDomain, int iComponent);

/* *************** add ****************************************** */

template<typename T, int nDim>
void add(MultiLevelTensorField3D<T,nDim> &A, 
         MultiLevelTensorField3D<T,nDim> &B, 
         MultiLevelTensorField3D<T,nDim> &aPlusB,
         const Box3D &domain, plint levelOfDomain);

template<typename T, int nDim>
std::auto_ptr<MultiLevelTensorField3D<T,nDim> > 
	add(MultiLevelTensorField3D<T,nDim> &A, 
     	MultiLevelTensorField3D<T,nDim> &B, 
        const Box3D &domain, plint levelOfDomain);

/* *************** addInPlace ****************************************** */

template<typename T, int nDim>
void addInPlace(MultiLevelTensorField3D<T,nDim> &A, 
                MultiLevelTensorField3D<T,nDim> &B, 
                const Box3D &domain, plint levelOfDomain);

/* *************** subtract ****************************************** */

template<typename T, int nDim>
void subtract(MultiLevelTensorField3D<T,nDim> &A, 
              MultiLevelTensorField3D<T,nDim> &B, 
           	  MultiLevelTensorField3D<T,nDim> &aMinusB,
           	  const Box3D &domain, plint levelOfDomain);

template<typename T, int nDim>
std::auto_ptr<MultiLevelTensorField3D<T,nDim> >  
	subtract(MultiLevelTensorField3D<T,nDim> &A, 
             MultiLevelTensorField3D<T,nDim> &B, 
           	 const Box3D &domain, plint levelOfDomain);

/* *************** subtractInPlace ****************************************** */

template<typename T, int nDim>
void subtractInPlace(MultiLevelTensorField3D<T,nDim> &A, 
              MultiLevelTensorField3D<T,nDim> &B, 
           	  const Box3D &domain, plint levelOfDomain);

/* *************** multiply ****************************************** */

template<typename T, int nDim>
void multiply(MultiLevelTensorField3D<T,nDim> &A, 
              MultiLevelTensorField3D<T,nDim> &B, 
           	  MultiLevelTensorField3D<T,nDim> &aMinusB,
           	  const Box3D &domain, plint levelOfDomain);

template<typename T, int nDim>
std::auto_ptr<MultiLevelTensorField3D<T,nDim> >   
	multiply(MultiLevelTensorField3D<T,nDim> &A, 
             MultiLevelTensorField3D<T,nDim> &B, 
           	 const Box3D &domain, plint levelOfDomain);

/* *************** multiplyInPlace ****************************************** */

template<typename T, int nDim>
void multiplyInPlace(MultiLevelTensorField3D<T,nDim> &A, 
                     MultiLevelTensorField3D<T,nDim> &B, 
           	         const Box3D &domain, plint levelOfDomain);

/* *************** multiply by scalar ****************************************** */

template<typename T, int nDim>
void multiply(MultiLevelTensorField3D<T,nDim> &A, 
              T alpha, MultiLevelTensorField3D<T,nDim> &aTimesB,
           	  const Box3D &domain, plint levelOfDomain);

template<typename T, int nDim>
std::auto_ptr<MultiLevelTensorField3D<T,nDim> >   
	multiply(MultiLevelTensorField3D<T,nDim> &A, 
             T alpha,
           	 const Box3D &domain, plint levelOfDomain);

/* *************** multiplyInPlace by scalar ****************************************** */

template<typename T, int nDim>
void multiplyInPlace(MultiLevelTensorField3D<T,nDim> &A, 
                     T alpha, 
           	         const Box3D &domain, plint levelOfDomain);

/* *************** divide ****************************************** */

template<typename T, int nDim>
void divide(MultiLevelTensorField3D<T,nDim> &A, 
            MultiLevelTensorField3D<T,nDim> &B, 
           	MultiLevelTensorField3D<T,nDim> &aMinusB,
           	const Box3D &domain, plint levelOfDomain);

template<typename T, int nDim>
std::auto_ptr<MultiLevelTensorField3D<T,nDim> >  
	divide(MultiLevelTensorField3D<T,nDim> &A, 
           MultiLevelTensorField3D<T,nDim> &B, 
           const Box3D &domain, plint levelOfDomain);

/* *************** divideInPlace ****************************************** */

template<typename T, int nDim>
void divideInPlace(MultiLevelTensorField3D<T,nDim> &A, 
                   MultiLevelTensorField3D<T,nDim> &B, 
           	       const Box3D &domain, plint levelOfDomain);

/* *************** Vorticity ****************************************** */

template<typename T>
std::auto_ptr<MultiLevelTensorField3D<T,3> >
   computeVorticity(const MultiLevelTensorField3D<T,3>& velocities, Box3D domain, plint levelOfDomain);

template<typename T>
std::auto_ptr<MultiLevelTensorFieldForOutput3D<T,3> >
   computeVorticity(const MultiLevelTensorField3D<T,3>& velocities, Box3D domain, plint levelOfDomain, bool crop);




}  // namespace plb

#endif  // DATA_ANALYSIS_WRAPPERS_3D_H
