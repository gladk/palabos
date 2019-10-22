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

/* Orestis Malaspinas contributed this code.
 */

/** \file
 * This object is a MRT LB dynamics as described in D.Yu et al. in
 * Progress in Aerospace Sciences 39 (2003) 329-367.
 * The external force terms are added here.
 */
#ifndef EXTERNAL_FORCE_MRT_DYNAMICS_H
#define EXTERNAL_FORCE_MRT_DYNAMICS_H

#include "core/globalDefs.h"
#include "basicDynamics/isoThermalDynamics.h"
#include "basicDynamics/externalForceDynamics.h"
#include "mrtDynamics.h"

namespace plb {
    
/// Implementation of the MRT collision step
template<typename T, template<typename U> class Descriptor>
class GuoExternalForceMRTdynamics : public ExternalForceDynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    GuoExternalForceMRTdynamics(T omega_);
    GuoExternalForceMRTdynamics(HierarchicUnserializer& unserializer);
    
    /// Clone the object on its dynamic type.
    virtual GuoExternalForceMRTdynamics<T,Descriptor>* clone() const;

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
    static int id;
};

/// Implementation of the MRT collision step with external force
/// and Smagorinsky model
template<typename T, template<typename U> class Descriptor>
class GuoExternalForceSmagorinskyMRTdynamics : public ExternalForceDynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    GuoExternalForceSmagorinskyMRTdynamics(T omega_, T cSmago_);
    GuoExternalForceSmagorinskyMRTdynamics(HierarchicUnserializer& unserializer);

    /// Clone the object on its dynamic type.
    virtual GuoExternalForceSmagorinskyMRTdynamics<T,Descriptor>* clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    /// Serialize the dynamics object.
    virtual void serialize(HierarchicSerializer& serializer) const;
    /// Un-Serialize the dynamics object.
    virtual void unserialize(HierarchicUnserializer& unserializer);

/* *************** Collision and Equilibrium ************************* */

    /// Implementation of the collision step
    virtual void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics& statistics_);

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;
private:
    T cSmago;
    static int id;
};


/// Implementation of the quasi incompressible MRT collision step with external force
/// and Smagorinsky model
template<typename T, template<typename U> class Descriptor>
class GuoExternalForceSmagorinskyIncMRTdynamics : public IncExternalForceDynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    GuoExternalForceSmagorinskyIncMRTdynamics(T omega_, T cSmago_);
    GuoExternalForceSmagorinskyIncMRTdynamics(HierarchicUnserializer& unserializer);

    /// Clone the object on its dynamic type.
    virtual GuoExternalForceSmagorinskyIncMRTdynamics<T,Descriptor>* clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    /// Serialize the dynamics object.
    virtual void serialize(HierarchicSerializer& serializer) const;
    /// Un-Serialize the dynamics object.
    virtual void unserialize(HierarchicUnserializer& unserializer);

    /// Say if velocity in this dynamics is computed as "j" (the order-1 moment
    ///   of the populations) or as "j/rho".
    virtual bool velIsJ() const;

/* *************** Collision and Equilibrium ************************* */

    /// Implementation of the collision step
    virtual void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics& statistics_);

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;
private:
    T cSmago;
    static int id;
};

/// Implementation of the MRT collision step with external force
/// and ConsistentSmagorinsky model
template<typename T, template<typename U> class Descriptor>
class GuoExternalForceConsistentSmagorinskyMRTdynamics : public ExternalForceDynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    GuoExternalForceConsistentSmagorinskyMRTdynamics(T omega_, T cSmago_);
    GuoExternalForceConsistentSmagorinskyMRTdynamics(HierarchicUnserializer& unserializer);

    /// Clone the object on its dynamic type.
    virtual GuoExternalForceConsistentSmagorinskyMRTdynamics<T,Descriptor>* clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    /// Serialize the dynamics object.
    virtual void serialize(HierarchicSerializer& serializer) const;
    /// Un-Serialize the dynamics object.
    virtual void unserialize(HierarchicUnserializer& unserializer);

/* *************** Collision and Equilibrium ************************* */

    /// Implementation of the collision step
    virtual void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics& statistics_);

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;
private:
    T cSmago;
    static int id;
};


/// Implementation of the quasi incompressible MRT collision step with external force
/// and ConsistentSmagorinsky model
template<typename T, template<typename U> class Descriptor>
class GuoExternalForceConsistentSmagorinskyIncMRTdynamics : public IncExternalForceDynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    GuoExternalForceConsistentSmagorinskyIncMRTdynamics(T omega_, T cSmago_);
    GuoExternalForceConsistentSmagorinskyIncMRTdynamics(HierarchicUnserializer& unserializer);

    /// Clone the object on its dynamic type.
    virtual GuoExternalForceConsistentSmagorinskyIncMRTdynamics<T,Descriptor>* clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    /// Serialize the dynamics object.
    virtual void serialize(HierarchicSerializer& serializer) const;
    /// Un-Serialize the dynamics object.
    virtual void unserialize(HierarchicUnserializer& unserializer);

    /// Say if velocity in this dynamics is computed as "j" (the order-1 moment
    ///   of the populations) or as "j/rho".
    virtual bool velIsJ() const;

/* *************** Collision and Equilibrium ************************* */

    /// Implementation of the collision step
    virtual void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics& statistics_);

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;
private:
    T cSmago;
    static int id;
};


/// Implementation of the MRT collision step with External Momenta
template<typename T, template<typename U> class Descriptor>
class GuoExternalForceAndMomentMRTdynamics : public ExternalForceDynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    GuoExternalForceAndMomentMRTdynamics(T omega_);
    GuoExternalForceAndMomentMRTdynamics(HierarchicUnserializer& unserializer);
    
    /// Clone the object on its dynamic type.
    virtual GuoExternalForceAndMomentMRTdynamics<T,Descriptor>* clone() const;

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
    static int id;
};





/// Implementation of the MRT collision step with He External force
template<typename T, template<typename U> class Descriptor>
class HeExternalForceMRTdynamics : public ExternalForceDynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    HeExternalForceMRTdynamics(T omega_);
    HeExternalForceMRTdynamics(HierarchicUnserializer& unserializer);
    
    /// Clone the object on its dynamic type.
    virtual HeExternalForceMRTdynamics<T,Descriptor>* clone() const;

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
    static int id;
};

}  // namespace plb

#endif  // MRT_DYNAMICS_H

