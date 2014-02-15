/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2013 FlowKit Sarl
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
    GuoExternalForceMRTdynamics(plint externalParam_);
    GuoExternalForceMRTdynamics(MRTparam<T,Descriptor>* param_);
    GuoExternalForceMRTdynamics(HierarchicUnserializer& unserializer);
    GuoExternalForceMRTdynamics(GuoExternalForceMRTdynamics<T,Descriptor> const& rhs);
    GuoExternalForceMRTdynamics<T,Descriptor>& operator=(GuoExternalForceMRTdynamics<T,Descriptor> const& rhs);
    virtual ~GuoExternalForceMRTdynamics<T,Descriptor>();
    void swap(GuoExternalForceMRTdynamics<T,Descriptor>& rhs);
    
    /// Clone the object on its dynamic type.
    virtual GuoExternalForceMRTdynamics<T,Descriptor>* clone() const;

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
                                 
    MRTparam<T,Descriptor> const& getMrtParameter() const;
private:
    MRTparam<T,Descriptor>* param;
    plint externalParam;
    static int id;
};

template<typename T, template<typename U> class Descriptor>
class VariableOmegaForcedMRTdynamics : public ExternalForceDynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    VariableOmegaForcedMRTdynamics(T omega);
    
    /// Clone the object on its dynamic type.
    virtual VariableOmegaForcedMRTdynamics<T,Descriptor>* clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

/* *************** Collision and Equilibrium ************************* */

    /// Implementation of the collision step
    virtual void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics& statistics_);

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;
                                 
    MRTparam<T,Descriptor> const& getMrtParameter() const;
private:
    static MRTparam<T,Descriptor> param;
    static int id;
};

/// Implementation of the MRT collision step with external force
/// and Smagorinsky model
template<typename T, template<typename U> class Descriptor>
class GuoExternalForceSmagorinskyMRTdynamics : public ExternalForceDynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    GuoExternalForceSmagorinskyMRTdynamics(plint externalParam_, T cSmago_);
    GuoExternalForceSmagorinskyMRTdynamics(MRTparam<T,Descriptor>* param_, T cSmago_);
    GuoExternalForceSmagorinskyMRTdynamics(HierarchicUnserializer& unserializer);
    GuoExternalForceSmagorinskyMRTdynamics(GuoExternalForceSmagorinskyMRTdynamics<T,Descriptor> const& rhs);
    GuoExternalForceSmagorinskyMRTdynamics<T,Descriptor>& operator=(GuoExternalForceSmagorinskyMRTdynamics<T,Descriptor> const& rhs);
    virtual ~GuoExternalForceSmagorinskyMRTdynamics<T,Descriptor>();
    void swap(GuoExternalForceSmagorinskyMRTdynamics<T,Descriptor>& rhs);

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

    MRTparam<T,Descriptor> const& getMrtParameter() const;
private:
    MRTparam<T,Descriptor>* param;
    plint externalParam;
    T cSmago;
    static int id;
};


/// Implementation of the quasi incompressible MRT collision step with external force
/// and Smagorinsky model
template<typename T, template<typename U> class Descriptor>
class GuoExternalForceSmagorinskyQuasiIncMRTdynamics : public ExternalForceDynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    GuoExternalForceSmagorinskyQuasiIncMRTdynamics(plint externalParam_, T cSmago_);
    GuoExternalForceSmagorinskyQuasiIncMRTdynamics(MRTparam<T,Descriptor>* param_, T cSmago_);
    GuoExternalForceSmagorinskyQuasiIncMRTdynamics(HierarchicUnserializer& unserializer);
    GuoExternalForceSmagorinskyQuasiIncMRTdynamics(GuoExternalForceSmagorinskyQuasiIncMRTdynamics<T,Descriptor> const& rhs);
    GuoExternalForceSmagorinskyQuasiIncMRTdynamics<T,Descriptor>& operator=(GuoExternalForceSmagorinskyQuasiIncMRTdynamics<T,Descriptor> const& rhs);
    virtual ~GuoExternalForceSmagorinskyQuasiIncMRTdynamics<T,Descriptor>();
    void swap(GuoExternalForceSmagorinskyQuasiIncMRTdynamics<T,Descriptor>& rhs);

    /// Clone the object on its dynamic type.
    virtual GuoExternalForceSmagorinskyQuasiIncMRTdynamics<T,Descriptor>* clone() const;

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

    MRTparam<T,Descriptor> const& getMrtParameter() const;
private:
    MRTparam<T,Descriptor>* param;
    plint externalParam;
    T cSmago;
    static int id;
};

/// Implementation of the quasi incompressible MRT collision step with external force
/// and Smagorinsky model
template<typename T, template<typename U> class Descriptor>
class GuoExternalForceConsistentSmagorinskyMRTdynamics : public ExternalForceDynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    GuoExternalForceConsistentSmagorinskyMRTdynamics(plint externalParam_, T cSmago_);
    GuoExternalForceConsistentSmagorinskyMRTdynamics(MRTparam<T,Descriptor>* param_, T cSmago_);
    GuoExternalForceConsistentSmagorinskyMRTdynamics(HierarchicUnserializer& unserializer);
    GuoExternalForceConsistentSmagorinskyMRTdynamics(GuoExternalForceConsistentSmagorinskyMRTdynamics<T,Descriptor> const& rhs);
    GuoExternalForceConsistentSmagorinskyMRTdynamics<T,Descriptor>& operator=(GuoExternalForceConsistentSmagorinskyMRTdynamics<T,Descriptor> const& rhs);
    virtual ~GuoExternalForceConsistentSmagorinskyMRTdynamics<T,Descriptor>();
    void swap(GuoExternalForceConsistentSmagorinskyMRTdynamics<T,Descriptor>& rhs);

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

    MRTparam<T,Descriptor> const& getMrtParameter() const;
private:
    MRTparam<T,Descriptor>* param;
    plint externalParam;
    T cSmago;
    static int id;
};

/// Implementation of the quasi incompressible MRT collision step with external force
/// and Smagorinsky model
template<typename T, template<typename U> class Descriptor>
class GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics : public ExternalForceDynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics(plint externalParam_, T cSmago_);
    GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics(MRTparam<T,Descriptor>* param_, T cSmago_);
    GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics(HierarchicUnserializer& unserializer);
    GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics(GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor> const& rhs);
    GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>& operator=(GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor> const& rhs);
    virtual ~GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>();
    void swap(GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>& rhs);

    /// Clone the object on its dynamic type.
    virtual GuoExternalForceConsistentSmagorinskyQuasiIncMRTdynamics<T,Descriptor>* clone() const;

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

    MRTparam<T,Descriptor> const& getMrtParameter() const;
private:
    MRTparam<T,Descriptor>* param;
    plint externalParam;
    T cSmago;
    static int id;
};

}  // namespace plb

#endif  // MRT_DYNAMICS_H

