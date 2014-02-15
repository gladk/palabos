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
 * Progress in Aerospace Sciences 39 (2003) 329-367
 */
#ifndef MRT_DYNAMICS_H
#define MRT_DYNAMICS_H

#include "core/globalDefs.h"
#include "basicDynamics/isoThermalDynamics.h"

namespace plb {

template<typename T, template<typename U> class Descriptor>
class MRTparam {
public:
    typedef T MatT[Descriptor<T>::q][Descriptor<T>::q];
public:
    MRTparam();
    virtual ~MRTparam() { }
    MRTparam(T omega_);
    MRTparam(T omega_, T lambda_);
    MRTparam(std::vector<T> s_);
    MRTparam(HierarchicUnserializer& unserializer);
    void serialize(HierarchicSerializer& serializer) const;
    void unserialize(HierarchicUnserializer& unserializer);
    MRTparam<T,Descriptor>* clone() const;
    /// Set local value of any generic parameter.
    virtual void setParameter(plint whichParameter, T value);
    /// Get local value of any generic parameter.
    virtual T getParameter(plint whichParameter) const;
    /// Get local relaxation parameter for shear viscosity.
    virtual T getOmega() const;
    /// Set local relaxation parameter for shear viscosity.
    virtual void setOmega(T omega_);
    /// Get local relaxation parameter for bulk viscosity.
    T getLambda() const;
    /// Set local relaxation parameter for bulk viscosity.
    void setLambda(T lambda_);
    /// Get local relaxation parameter for q parameter.
    T getOmegaQ() const;
    /// Set local relaxation parameter for q parameter.
    void setOmegaQ(T q_);
    /// Get local relaxation parameter for epsilon parameter.
    T getOmegaEpsilon() const;
    /// Set local relaxation parameter for epsilon parameter.
    void setOmegaEpsilon(T epsilon_);
    /// Returns the relaxation frequencies vector.
    std::vector<T> const& getS() const;
    /// Returns relaxation time matrix.
    MatT& getInvM();
private:
    /// Initializes the S vector from omega and lambda.
    void iniRelaxationVector();
    /// Precomputes relaxation time matrix for arbitrary relaxation times.
    void precompute_invM_S();
private:
    T invM_S[Descriptor<T>::q][Descriptor<T>::q]; // relaxation time matrix.
    std::vector<T> s;
};

template<typename T, template<typename U> class Descriptor>
class MRTparamList {
public:
    void set(plint id, MRTparam<T,Descriptor> const& param);
    MRTparam<T,Descriptor>& get(plint id);
    MRTparam<T,Descriptor> const& get(plint id) const;
private:
    std::vector<MRTparam<T,Descriptor> > parameters;
};

template<typename T, template<typename U> class Descriptor>
MRTparamList<T,Descriptor>& mrtParam();

/// Implementation of the MRT collision step
template<typename T, template<typename U> class Descriptor>
class MRTdynamics : public IsoThermalBulkDynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    MRTdynamics(plint externalParam_);
    MRTdynamics(MRTparam<T,Descriptor>* param_);
    MRTdynamics(HierarchicUnserializer& unserializer);
    MRTdynamics(MRTdynamics<T,Descriptor> const& rhs);
    MRTdynamics<T,Descriptor>& operator=(MRTdynamics<T,Descriptor> const& rhs);
    virtual ~MRTdynamics<T,Descriptor>();
    void swap(MRTdynamics<T,Descriptor>& rhs);
    
    /// Clone the object on its dynamic type.
    virtual MRTdynamics<T,Descriptor>* clone() const;

    virtual void setOmega(T omega_);

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

    /// Implementation of the collision step, with imposed macroscopic variables
    virtual void collideExternal (
            Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
            T thetaBar, BlockStatistics& stat );

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;
                                 
    MRTparam<T,Descriptor> const& getMrtParameter() const;
private:
    MRTparam<T,Descriptor>* param;
    plint externalParam;
    static int id;
};

/// Implementation of the MRT collision step
template<typename T, template<typename U> class Descriptor>
class VariableOmegaMRTdynamics : public IsoThermalBulkDynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    VariableOmegaMRTdynamics(T omega);

    /// Clone the object on its dynamic type.
    virtual VariableOmegaMRTdynamics<T,Descriptor>* clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

/* *************** Collision and Equilibrium ************************* */

    /// Implementation of the collision step
    virtual void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics& statistics_);

    /// Implementation of the collision step, with imposed macroscopic variables
    virtual void collideExternal (
            Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
            T thetaBar, BlockStatistics& stat );

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;
                                 
    MRTparam<T,Descriptor> const& getMrtParameter() const;
private:
    static MRTparam<T,Descriptor> param;
    static int id;
};

/// Implementation of the MRT collision step
template<typename T, template<typename U> class Descriptor>
class ExternalVelocityMRTdynamics : public IsoThermalBulkDynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    ExternalVelocityMRTdynamics(plint externalParam_);
    ExternalVelocityMRTdynamics(MRTparam<T,Descriptor>* param_);
    ExternalVelocityMRTdynamics(HierarchicUnserializer& unserializer);
    ExternalVelocityMRTdynamics(ExternalVelocityMRTdynamics<T,Descriptor> const& rhs);
    ExternalVelocityMRTdynamics<T,Descriptor>& operator=(ExternalVelocityMRTdynamics<T,Descriptor> const& rhs);
    virtual ~ExternalVelocityMRTdynamics<T,Descriptor>();
    void swap(ExternalVelocityMRTdynamics<T,Descriptor>& rhs);
    
    /// Clone the object on its dynamic type.
    virtual ExternalVelocityMRTdynamics<T,Descriptor>* clone() const;

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
                                 
/* *************** Moments ******************************************* */

    // The function computeDensity() is not overridden, and the default
    //   implementation is kept, for two reasons. First, it is equivalent to access
    //   rho from the external scalar or to recompute it from the bulk (this
    //   is not the case for the velocity, in a Shan/Chen multicomponent model).
    //   Second, the Shan/Chen data-processor needs computeDensity() to be
    //   default implemented, because it uses this function to treat walls with
    //   a virtual-density mechanism.

    /** Accesses velocity from external scalar. **/
    virtual void computeVelocity( Cell<T,Descriptor> const& cell,
                                  Array<T,Descriptor<T>::d>& u ) const;
    /// Compute order-0 moment rho-bar
    /** Accesses rhoBar from external scalar. **/
    virtual T computeRhoBar(Cell<T,Descriptor> const& cell) const;

    /// Compute order-0 moment rho-bar and order-1 moment j
    /** Accesses rhoBar and j from external scalar. **/
    virtual void computeRhoBarJ(Cell<T,Descriptor> const& cell,
                                T& rhoBar, Array<T,Descriptor<T>::d>& j) const;
                                 
public:
    MRTparam<T,Descriptor> const &getMrtParameter() const {
        return *param;
    }
private:
    MRTparam<T,Descriptor>* param;
    plint externalParam;
    static int id;
};

/// Implementation of incompressible MRT dynamics.
/** This is the MRT equivalent of IncBGKdynamics: the "rho" moment of the
 *  populations appears only as a pressure term in the equilibrium, while
 *  the other terms are multiplied by the constant rho0.
 **/
template<typename T, template<typename U> class Descriptor>
class IncMRTdynamics : public IsoThermalBulkDynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    IncMRTdynamics(plint externalParam_);
    IncMRTdynamics(MRTparam<T,Descriptor>* param_);
    IncMRTdynamics(HierarchicUnserializer& unserializer);
    IncMRTdynamics(IncMRTdynamics<T,Descriptor> const& rhs);
    IncMRTdynamics<T,Descriptor>& operator=(IncMRTdynamics<T,Descriptor> const& rhs);
    virtual ~IncMRTdynamics<T,Descriptor>();
    void swap(IncMRTdynamics<T,Descriptor>& rhs);
    
    /// Clone the object on its dynamic type.
    virtual IncMRTdynamics<T,Descriptor>* clone() const;

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
                                 
    virtual bool velIsJ() const;

/* *************** Macroscopic variables ***************************** */
    
    /// Velocity is equal to j, not u.
    virtual void computeVelocity( Cell<T,Descriptor> const& cell,
                                  Array<T,Descriptor<T>::d>& u ) const;
public:
    MRTparam<T,Descriptor> const &getMrtParameter() const {
        return *param;
    }
private:
    MRTparam<T,Descriptor>* param;
    plint externalParam;
    static int id;
};


/// Implementation of incompressible MRT dynamics, with velocity taken from an external scalar.
/** This is the MRT equivalent of IncBGKdynamics: the "rho" moment of the
 *  populations appears only as a pressure term in the equilibrium, while
 *  the other terms are multiplied by the constant rho0.
 **/
template<typename T, template<typename U> class Descriptor>
class ExternalVelocityIncMRTdynamics : public IsoThermalBulkDynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ************************ */
    ExternalVelocityIncMRTdynamics(plint externalParam_);
    ExternalVelocityIncMRTdynamics(MRTparam<T,Descriptor>* param_);
    ExternalVelocityIncMRTdynamics(HierarchicUnserializer& unserializer);
    ExternalVelocityIncMRTdynamics(ExternalVelocityIncMRTdynamics<T,Descriptor> const& rhs);
    ExternalVelocityIncMRTdynamics<T,Descriptor>& operator=(ExternalVelocityIncMRTdynamics<T,Descriptor> const& rhs);
    virtual ~ExternalVelocityIncMRTdynamics<T,Descriptor>();
    void swap(ExternalVelocityIncMRTdynamics<T,Descriptor>& rhs);
    
    /// Clone the object on its dynamic type.
    virtual ExternalVelocityIncMRTdynamics<T,Descriptor>* clone() const;

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
                                 
    virtual bool velIsJ() const;
    
/* *************** Moments ******************************************* */

    // The function computeDensity() is not overridden, and the default
    //   implementation is kept, for two reasons. First, it is equivalent to access
    //   rho from the external scalar or to recompute it from the bulk (this
    //   is not the case for the velocity, in a Shan/Chen multicomponent model).
    //   Second, the Shan/Chen data-processor needs computeDensity() to be
    //   default implemented, because it uses this function to treat walls with
    //   a virtual-density mechanism.

    /** Accesses velocity from external scalar. **/
    virtual void computeVelocity( Cell<T,Descriptor> const& cell,
                                  Array<T,Descriptor<T>::d>& u ) const;
    /// Compute order-0 moment rho-bar
    /** Accesses rhoBar from external scalar. **/
    virtual T computeRhoBar(Cell<T,Descriptor> const& cell) const;

    /// Compute order-0 moment rho-bar and order-1 moment j
    /** Accesses rhoBar and j from external scalar. **/
    virtual void computeRhoBarJ(Cell<T,Descriptor> const& cell,
                                T& rhoBar, Array<T,Descriptor<T>::d>& j) const;
public:
    MRTparam<T,Descriptor> const &getMrtParameter() const {
        return *param;
    }
private:
    MRTparam<T,Descriptor>* param;
    plint externalParam;
    static int id;
};

}  // namespace plb

#endif  // MRT_DYNAMICS_H

