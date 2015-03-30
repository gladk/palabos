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

/* Orestis Malaspinas contributed this code.
 */

/** \file
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- header file.
 */
#ifndef ENTROPIC_LB_DYNAMICS_H
#define ENTROPIC_LB_DYNAMICS_H

#include "core/globalDefs.h"
#include "basicDynamics/isoThermalDynamics.h"

namespace plb {

template<typename T, template<typename U> class Descriptor> class Cell;


/// Implementation of the entropic collision step
template<typename T, template<typename U> class Descriptor>
class EntropicDynamics : public IsoThermalBulkDynamics<T,Descriptor> 
{
public:
/* *************** Construction / Destruction ************************ */
    EntropicDynamics(T omega_);

    /// Clone the object on its dynamic type.
    virtual EntropicDynamics<T,Descriptor>* clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

/* *************** Collision and Equilibrium ************************* */

    /// Implementation of the collision step
    virtual void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics& statistics_);

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;
private:
    /// computes the entropy function H(f)=sum_i f_i*ln(f_i/t_i)
    T computeEntropy(Array<T,Descriptor<T>::q> const& f);
    /// computes the entropy growth H(f)-H(f-alpha*fNeq)
    T computeEntropyGrowth(Array<T,Descriptor<T>::q> const& f, Array<T,Descriptor<T>::q> const& fNeq, T alpha);
    /// computes the entropy growth derivative
    /// dH/dalpha=-sum_i fNeq_i*ln((f_i-alpha*fNeq_i)/t_i)
    T computeEntropyGrowthDerivative(Array<T,Descriptor<T>::q> const& f, Array<T,Descriptor<T>::q> const& fNeq, T alpha);
    /// Get the alpha parameter
    bool getAlpha(T &alpha, Array<T,Descriptor<T>::q> const& f, Array<T,Descriptor<T>::q> const& fNeq);
private:
    static int id;
};

/// Implementation of the forced entropic collision step
template<typename T, template<typename U> class Descriptor>
class ForcedEntropicDynamics : public IsoThermalBulkDynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    ForcedEntropicDynamics(T omega_);

    /// Clone the object on its dynamic type.
    virtual ForcedEntropicDynamics<T,Descriptor>* clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

/* *************** Collision and Equilibrium ************************* */

    /// Implementation of the collision step
    virtual void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics& statistics_);

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;
private:
    /// computes the entropy function H(f)=sum_i f_i*ln(f_i/t_i)
    T computeEntropy(Array<T,Descriptor<T>::q> const& f);
    /// computes the entropy growth H(f)-H(f-alpha*fNeq)
    T computeEntropyGrowth(Array<T,Descriptor<T>::q> const& f, Array<T,Descriptor<T>::q> const& fNeq, T alpha);
    /// computes the entropy growth derivative
    /// dH/dalpha=-sum_i fNeq_i*ln((f_i-alpha*fNeq_i)/t_i)
    T computeEntropyGrowthDerivative(Array<T,Descriptor<T>::q> const& f, Array<T,Descriptor<T>::q> const& fNeq, T alpha);
    /// Get the alpha parameter
    bool getAlpha(T &alpha, Array<T,Descriptor<T>::q> const& f, Array<T,Descriptor<T>::q> const& fNeq);
    
    static const int forceBeginsAt = Descriptor<T>::ExternalField::forceBeginsAt;
    static const int sizeOfForce   = Descriptor<T>::ExternalField::sizeOfForce;
    static int id;
};

}  // namespace plb

#endif  // ENTROPIC_LB_DYNAMICS_H
