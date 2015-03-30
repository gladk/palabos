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
 * Helper functions for domain initialization -- header file.
 */
#ifndef FINITE_DIFFERENCE_WRAPPER_3D_H
#define FINITE_DIFFERENCE_WRAPPER_3D_H

#include "core/globalDefs.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataField3D.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiDataField3D.h"
#include "dataProcessors/dataAnalysisFunctional3D.h"
#include <memory>


namespace plb {

template<typename T>
void computeXderivative(MultiScalarField3D<T>& value, MultiScalarField3D<T>& derivative, Box3D const& domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeXderivative(MultiScalarField3D<T>& value, Box3D const& domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeXderivative(MultiScalarField3D<T>& value);

template<typename T>
void computeYderivative(MultiScalarField3D<T>& value, MultiScalarField3D<T>& derivative, Box3D const& domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeYderivative(MultiScalarField3D<T>& value, Box3D const& domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeYderivative(MultiScalarField3D<T>& value);

template<typename T>
void computeZderivative(MultiScalarField3D<T>& value, MultiScalarField3D<T>& derivative, Box3D const& domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeZderivative(MultiScalarField3D<T>& value, Box3D const& domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeZderivative(MultiScalarField3D<T>& value);

template<typename T>
void computeGradientNorm(MultiScalarField3D<T>& value, MultiScalarField3D<T>& derivative, Box3D const& domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeGradientNorm(MultiScalarField3D<T>& value, Box3D const& domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeGradientNorm(MultiScalarField3D<T>& value);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computePoissonRHS(MultiTensorField3D<T,3>& velocity, Box3D const& domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computePoissonRHS(MultiTensorField3D<T,3>& velocity);


template<typename T>
void poissonIterate(MultiScalarField3D<T>& oldPressure, MultiScalarField3D<T>& newPressure,
                    MultiScalarField3D<T>& rhs, T beta, Box3D const& domain, plint boundaryWidth = 1);

template<typename T>
T computePoissonResidue(MultiScalarField3D<T>& pressure, MultiScalarField3D<T>& rhs, Box3D const& domain);

// ========================================================================= //
// PERIODIC VERSIONS OF THE DERIVATIVES AND POISSON SCHEMES //
// ========================================================================= //

template<typename T>
void computeXperiodicDerivative(MultiScalarField3D<T>& value, MultiScalarField3D<T>& derivative, Box3D const& domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeXperiodicDerivative(MultiScalarField3D<T>& value, Box3D const& domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeXperiodicDerivative(MultiScalarField3D<T>& value);

template<typename T>
void computeYperiodicDerivative(MultiScalarField3D<T>& value, MultiScalarField3D<T>& derivative, Box3D const& domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeYperiodicDerivative(MultiScalarField3D<T>& value, Box3D const& domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeYperiodicDerivative(MultiScalarField3D<T>& value);

template<typename T>
void computeZperiodicDerivative(MultiScalarField3D<T>& value, MultiScalarField3D<T>& derivative, Box3D const& domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeZperiodicDerivative(MultiScalarField3D<T>& value, Box3D const& domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computeZperiodicDerivative(MultiScalarField3D<T>& value);

template<typename T>
void computePeriodicGradientNorm(MultiScalarField3D<T>& value, MultiScalarField3D<T>& derivative, Box3D const& domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computePeriodicGradientNorm(MultiScalarField3D<T>& value, Box3D const& domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computePeriodicGradientNorm(MultiScalarField3D<T>& value);

template<typename T>
void computePeriodicGradient(MultiScalarField3D<T>& value, MultiTensorField3D<T,3>& derivative, Box3D const& domain);

template<typename T>
std::auto_ptr<MultiTensorField3D<T,2> > computePeriodicGradient(MultiScalarField3D<T>& value, Box3D const& domain);

template<typename T>
std::auto_ptr<MultiTensorField3D<T,2> > computePeriodicGradient(MultiScalarField3D<T>& value);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computePeriodicPoissonRHS(MultiTensorField3D<T,3>& velocity, Box3D const& domain);

template<typename T>
std::auto_ptr<MultiScalarField3D<T> > computePeriodicPoissonRHS(MultiTensorField3D<T,3>& velocity);


template<typename T>
void periodicPoissonIterate(MultiScalarField3D<T>& oldPressure, MultiScalarField3D<T>& newPressure,
                           MultiScalarField3D<T>& rhs, T beta, Box3D const& domain);

}  // namespace plb

#endif  // FINITE_DIFFERENCE_WRAPPER_3D_H
