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
 * A collection of dynamics classes (like BGK) with which a Cell
 * object can be instantiated -- header file.
 */
#ifndef DYNAMICS_H
#define DYNAMICS_H

#include "core/globalDefs.h"
#include "core/util.h"
#include "core/blockStatistics.h"
#include "core/array.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"

namespace plb {

namespace dynamicParams {
    // Use 0-99 for relaxation parameters
    const plint omega_shear = 0;
    const plint omega_bulk  = 1;
    
    const plint omega_epsilon  = 2;
    const plint omega_q        = 3;
    
    const plint psi        = 4; // used for complete trt dynamics

    // Use 100-199 for material constants
    const plint sqrSpeedOfSound = 100; // Speed of sound squared

    // Use 1000 and higher for custom user-defined constants
    const plint smagorinskyConstant = 1010;
    const plint dynamicOmega = 1011;
}

template<typename T, template<typename U> class Descriptor> class Cell;
class HierarchicSerializer;
class HierarchicUnserializer;


/// Interface for the dynamics classes
template<typename T, template<typename U> class Descriptor>
struct Dynamics {
/* *************** Construction, Descruction, and unique identifier ********* */

    /// Destructor: virtual to enable inheritance
    virtual ~Dynamics() { }

    /// Clone the object, based on its dynamic type
    virtual Dynamics<T,Descriptor>* clone() const =0;

    /// Get an ID for the class which is unique for a given choice of T and Descriptor.
    virtual int getId() const;

    /// Say if velocity in this dynamics is computed as "j" (the order-1 moment
    ///   of the populations) or as "j/rho".
    virtual bool velIsJ() const;

    /// Say if this dynamics is based on another dynamics class (and is convertible
    ///   to a CompositeDynamics).
    virtual bool isComposite() const;

    /// Say if this class can be used as a base for a composite dynamics.
    virtual bool isComposeable() const;

    /// Say if the dynamics represents a boundary condition.
    virtual bool isBoundary() const;

    // Say if the dynamics has non-local components.
    virtual bool isNonLocal() const;

    // Say if the dynamics implements advection-diffusion.
    virtual bool isAdvectionDiffusion() const;

    // Say if the dynamics implements an entropic model.
    virtual bool isEntropic() const;

    // Say if the dynamics has meaningful moments (for example for BounceBack,
    // SpecularReflection and NoDynamics the moments of the populations have
    // no meaning).
    virtual bool hasMoments() const;

    /// Serialize the dynamics object.
    virtual void serialize(HierarchicSerializer& serializer) const;
    /// Un-Serialize the dynamics object.
    virtual void unserialize(HierarchicUnserializer& unserializer);

/* *************** Collision, Equilibrium, and Non-equilibrium ************** */

    /// Implementation of the collision step
    virtual void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics& statistics_) =0;

    /// Implementation of the collision step, with imposed macroscopic variables
    virtual void collideExternal(Cell<T,Descriptor>& cell, T rhoBar,
                         Array<T,Descriptor<T>::d> const& j, T thetaBar, BlockStatistics& stat);

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const =0;

    virtual void computeEquilibria( Array<T,Descriptor<T>::q>& fEq,  T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                    T jSqr, T thetaBar=T() ) const;

    /// Re-compute particle populations from the leading moments
    virtual void regularize(Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                            T jSqr, Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq, T thetaBar=T() ) const =0;

/* *************** Computation of macroscopic variables ********************* */

    /// Compute the local particle density in lattice units
    virtual T computeDensity(Cell<T,Descriptor> const& cell) const =0;
    /// Compute the local pressure in lattice units
    virtual T computePressure(Cell<T,Descriptor> const& cell) const =0; /// Compute the local fluid velocity in lattice units
    virtual void computeVelocity( Cell<T,Descriptor> const& cell,
                                  Array<T,Descriptor<T>::d>& u ) const =0;
    /// Compute the temperature in lattice units
    virtual T computeTemperature(Cell<T,Descriptor> const& cell) const =0;
    /// Compute the "off-equilibrium part of Pi"
    virtual void computePiNeq (
        Cell<T,Descriptor> const& cell, Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq ) const =0;
    /// Compute the deviatoric stress tensor
    virtual void computeShearStress (
        Cell<T,Descriptor> const& cell, Array<T,SymmetricTensor<T,Descriptor>::n>& stress ) const =0;
    /// Compute the heat flux in lattice units
    virtual void computeHeatFlux( Cell<T,Descriptor> const& cell,
                                  Array<T,Descriptor<T>::d>& q ) const =0;
    /// Compute additional user-defined moments
    virtual void computeMoment( Cell<T,Descriptor> const& cell,
                                plint momentId, T* moment ) const =0;

/* *************** Access to Dynamics variables, e.g. omega ***************** */

    /// Set all relaxation frequencies (in case of an SRT model they should all equal omega).
    virtual void setRelaxationFrequencies(Array<T, Descriptor<T>::q> const& frequencies);

    /// Get all relaxation frequencies (in case of an SRT model they are all equal to omega).
    virtual Array<T, Descriptor<T>::q> getRelaxationFrequencies() const;

    /// Get local relaxation parameter of the dynamics.
    virtual T getOmega() const =0;

    /// Set local relaxation parameter of the dynamics.
    virtual void setOmega(T omega_) =0;

    /// Get local value of any generic parameter.
    virtual T getParameter(plint whichParameter) const;

    /// Get local value of any generic parameter which depends on the cell value.
    virtual T getDynamicParameter(plint whichParameter, Cell<T,Descriptor> const& cell) const;

    /// Set local value of any generic parameter.
    virtual void setParameter(plint whichParameter, T value);

/* *************** Switch between population and moment representation ****** */

    /// Number of variables required to decompose a population representation into moments.
    /** Example: Take the BGK dynamics on a D2Q9 lattice. The dynamics is physically
     *    significant up to order-1 in the Chapman-Enskog expansion. Therefore, a decomposition
     *    means:
     *    - At order 0: Decompose into rho, u, and fNeq (1+2+9=12 variables)
     *    - At order 1: Decompose into rho, u, and PiNeq (1+2+3=6 variables)
     *    - At higher order: Decompose according to order 1.
     */
    virtual plint numDecomposedVariables(plint order) const =0;

    /// Decompose from population representation into moment representation.
    /**   \sa numDecomposedVariables()
     *    \param rawData This vector will be resized automatically inside decomponse(),
     *                   if it doesn't already have the right size. This convention should
     *                   be respected by classes inheriting from Dynamics.
     */
    virtual void decompose(Cell<T,Descriptor> const& cell, std::vector<T>& rawData, plint order) const =0;

    /// Recompose from moment representation to population representation.
    /**   \sa numDecomposedVariables()
     *    This process is also known as "regularization step", and this function is therefore
     *    equivalent to regularize(), although one or the other function may be more useful
     *    in a specific context, due to the form of the parameters.
     */
    virtual void recompose(Cell<T,Descriptor>& cell, std::vector<T> const& rawData, plint order) const =0;

    /// Change the space and time scales of the variables in moment representation.
    /**   \sa numDecomposedVariables()
     *    \param xDxInv Factor by which space scale is multiplied.
     *    \param xDt Factor by which time scale is multiplied.
     */
    virtual void rescale(std::vector<T>& rawData, T xDxInv, T xDt, plint order) const =0;

    /// Adapt internal variables to a new grid level.
    /** In the default implementation, only omega is rescaled.
     **/
    virtual void rescale(int dxScale, int dtScale);

/* *************** Access Cell raw data through Dynamics ******************** */

    /// Access particle populations through the dynamics object.
    /** Default implementation: access cell directly. */
    virtual void getPopulations(Cell<T,Descriptor> const& cell, Array<T,Descriptor<T>::q>& f) const;

    /// Access external fields through the dynamics object.
    /** Default implementation: access cell directly. */
    virtual void getExternalField (
            Cell<T,Descriptor> const& cell, plint pos, plint size, T* ext ) const;

    /// Define particle populations through the dynamics object.
    /** Default implementation: access cell directly. */
    virtual void setPopulations(Cell<T,Descriptor>& cell, Array<T,Descriptor<T>::q> const& f);

    /// Define external fields through the dynamics object.
    /** Default implementation: access cell directly. */
    virtual void setExternalField (
            Cell<T,Descriptor>& cell, plint pos, plint size, const T* ext);


/* *************** Define macroscopic variables, e.g. on boundaries ********* */

    /// Define density, if possible. Does nothing by default.
    virtual void defineDensity(Cell<T,Descriptor>& cell, T density);

    /// Define velocity, if possible. Does nothing by default.
    virtual void defineVelocity(Cell<T,Descriptor>& cell, Array<T,Descriptor<T>::d> const& u);

    /// Define temperature, if possible. Does nothing by default.
    virtual void defineTemperature(Cell<T,Descriptor>& cell, T temperature);

    /// Define heat flux, if possible. Does nothing by default.
    virtual void defineHeatFlux(Cell<T,Descriptor>& cell, Array<T,Descriptor<T>::d> const& q);

    /// Define deviatoric stress, if possible. Does nothing by default.
    virtual void definePiNeq(Cell<T,Descriptor>& cell,
                                        Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq);

    /// Define user define moment, if possible. Does nothing by default.
    virtual void defineMoment(Cell<T,Descriptor>& cell, plint momentId, T const* value);

/* *************** Additional moments, intended for internal use ************ */

    /// Compute order-0 moment rho-bar
    virtual T computeRhoBar(Cell<T,Descriptor> const& cell) const =0;

    /// Compute order-0 moment rho-bar and order-1 moment j
    virtual void computeRhoBarJ(Cell<T,Descriptor> const& cell,
                                T& rhoBar, Array<T,Descriptor<T>::d>& j) const =0;

    /// Compute e-bar, which is related to the internal energy
    virtual T computeEbar(Cell<T,Descriptor> const& cell) const =0;

    /// Compute order-0 moment rho-bar, order-1 moment j, and order-2
    ///   off-equilibrium moment PiNeq.
    virtual void computeRhoBarJPiNeq(Cell<T,Descriptor> const& cell,
                                     T& rhoBar, Array<T,Descriptor<T>::d>& j,
                                     Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq) const =0;

/* ********* Update local variables when the f[iPop] are all known ********** */
 // (needed e.g. when the dynamics of a cell is reset after initialize(),
 // but the collide-step relies on variables, which are usually only calculated
 // within integrate process functionals
    /// by default, nothing happens
    void updateLocalVariables(Cell<T,Descriptor>& cell){ }
};


/// Common base for bulk dynamics
/** In this class, all methods have a default implementation, except
 *  for clone(), collide(), computeEquilibrium(), and regularize(),
 *  which are explicitly dynamics dependent.
 *  In this way, it is easy to inherit from BasicBulkDynamics to define
 *  a new Dynamics class. If you don't care about temperature, for
 *  example, you don't need to define methods such as computeTemperature().
 *  Computation of density, velocity, and pressure default to their BGK,
 *  whereas the higher order moments default to 0.
 *  Also, this class stores the value of the "main" relaxation parameter
 *  omega, and offers access to this parameter.
 */
template<typename T, template<typename U> class Descriptor>
class BasicBulkDynamics : public Dynamics<T,Descriptor> {
public:
/* *************** Construction and Destruction ***************************** */
    BasicBulkDynamics(T omega_);

/* *************** Computation of macroscopic variables ********************* */

    /// Compute the local particle density in lattice units
    virtual T computeDensity(Cell<T,Descriptor> const& cell) const;
    /// Compute the local pressure in lattice units
    virtual T computePressure(Cell<T,Descriptor> const& cell) const;
    /// Compute the local fluid velocity in lattice units
    virtual void computeVelocity( Cell<T,Descriptor> const& cell,
                                  Array<T,Descriptor<T>::d>& u ) const;
    /// Compute the temperature in lattice units
    virtual T computeTemperature(Cell<T,Descriptor> const& cell) const;
    /// Compute the "off-equilibrium part of Pi"
    virtual void computePiNeq (
        Cell<T,Descriptor> const& cell, Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq ) const;
    /// Compute the deviatoric stress tensor
    virtual void computeShearStress (
        Cell<T,Descriptor> const& cell, Array<T,SymmetricTensor<T,Descriptor>::n>& stress ) const;
    /// Compute the heat flux in lattice units
    virtual void computeHeatFlux( Cell<T,Descriptor> const& cell,
                                  Array<T,Descriptor<T>::d>& q ) const;
    /// Compute additional user-defined moments
    virtual void computeMoment( Cell<T,Descriptor> const& cell,
                                plint momentId, T* moment ) const;

/* *************** Access to Dynamics variables, e.g. omega ***************** */
    /// Get local relaxation parameter of the dynamics
    virtual T getOmega() const;

    /// Set local relaxation parameter of the dynamics
    virtual void setOmega(T omega_);

/* *************** Additional moments, intended for internal use ************ */

    /// Compute order-0 moment rho-bar
    virtual T computeRhoBar(Cell<T,Descriptor> const& cell) const;

    /// Compute order-0 moment rho-bar and order-1 moment j
    virtual void computeRhoBarJ(Cell<T,Descriptor> const& cell,
                                T& rhoBar, Array<T,Descriptor<T>::d>& j) const;

    /// Compute order-0 moment rho-bar, order-1 moment j, and order-2
    ///   off-equilibrium moment PiNeq.
    virtual void computeRhoBarJPiNeq(Cell<T,Descriptor> const& cell,
                                     T& rhoBar, Array<T,Descriptor<T>::d>& j,
                                     Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq) const;

    /// Compute e-bar, which is related to the internal energy
    virtual T computeEbar(Cell<T,Descriptor> const& cell) const;

private:
    T omega; ///< Main relaxation parameter, related to shear viscosity
};


/// Base class for composite dynamics
/** During the collision of a composite dynamics, a preparation step is
 *  invoked primary to collision. Then, another dynamics is called to
 *  fullfill collision.
 **/
template<typename T, template<typename U> class Descriptor>
class CompositeDynamics : public Dynamics<T,Descriptor> {
public:
/* *************** Construction and Destruction ***************************** */

    CompositeDynamics(Dynamics<T,Descriptor>* baseDynamics_, bool automaticPrepareCollision_);
    CompositeDynamics(CompositeDynamics<T,Descriptor> const& rhs);
    virtual ~CompositeDynamics();
    CompositeDynamics& operator=(CompositeDynamics<T,Descriptor> const& rhs);
    virtual CompositeDynamics<T,Descriptor>* clone() const =0;
    /// Clone the object with new BaseDynamics
    virtual CompositeDynamics<T,Descriptor>* cloneWithNewBase (
            Dynamics<T,Descriptor>* baseDynamics_ ) const;
    virtual bool isComposite() const;
    /// Serialize the dynamics object.
    virtual void serialize(HierarchicSerializer& serializer) const;
    /// Un-Serialize the dynamics object.
    virtual void unserialize(HierarchicUnserializer& unserializer);

/* *************** Access to base Dynamics ********************************** */

    virtual void replaceBaseDynamics(Dynamics<T,Descriptor>* newBaseDynamics);
    Dynamics<T,Descriptor>& getBaseDynamics();
    Dynamics<T,Descriptor> const& getBaseDynamics() const;

/* *************** Methods to be overloaded to configure behavior  ********** */

    /// Do something before execution of main collision step.
    virtual void prepareCollision(Cell<T,Descriptor>& cell) =0;

/* *************** Collision, Equilibrium, and Non-equilibrium ************** */

    /// Implementation of the collision step
    virtual void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics& statistics);

    /// Implementation of the collision step, with imposed macroscopic variables
    virtual void collideExternal(Cell<T,Descriptor>& cell, T rhoBar,
                         Array<T,Descriptor<T>::d> const& j, T thetaBar, BlockStatistics& stat);

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;

    virtual void computeEquilibria( Array<T,Descriptor<T>::q>& fEq,  T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                    T jSqr, T thetaBar=T() ) const;

    /// Re-compute particle populations from the leading moments
    virtual void regularize(Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                            T jSqr, Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq, T thetaBar=T() ) const;

/* *************** Computation of macroscopic variables ********************* */
    /// Say if velocity in this dynamics is computed as "j" (the order-1 moment
    ///   of the populations) or as "j/rho".
    virtual bool velIsJ() const;

    /// Compute the local particle density in lattice units
    virtual T computeDensity(Cell<T,Descriptor> const& cell) const;

    /// Compute the local pressure in lattice units
    virtual T computePressure(Cell<T,Descriptor> const& cell) const;

    /// Compute the local fluid velocity in lattice units
    virtual void computeVelocity( Cell<T,Descriptor> const& cell,
                                  Array<T,Descriptor<T>::d>& u ) const;

    /// Compute the temperature in lattice units
    virtual T computeTemperature(Cell<T,Descriptor> const& cell) const;

    /// Compute the "off-equilibrium part of Pi"
    virtual void computePiNeq (
        Cell<T,Descriptor> const& cell, Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq ) const;
    /// Compute the deviatoric stress tensor
    virtual void computeShearStress (
        Cell<T,Descriptor> const& cell, Array<T,SymmetricTensor<T,Descriptor>::n>& stress ) const;

    /// Compute the heat flux in lattice units
    virtual void computeHeatFlux( Cell<T,Descriptor> const& cell,
                                  Array<T,Descriptor<T>::d>& q ) const;

    /// Compute additional user-defined moments
    virtual void computeMoment( Cell<T,Descriptor> const& cell,
                                plint momentId, T* moment ) const;

/* *************** Switch between population and moment representation ****** */

    /// Number of variables required to decompose a population representation into moments.
    virtual plint numDecomposedVariables(plint order) const;

    /// Decompose from population representation into moment representation.
    virtual void decompose(Cell<T,Descriptor> const& cell, std::vector<T>& rawData, plint order) const;

    /// Recompose from moment representation to population representation.
    virtual void recompose(Cell<T,Descriptor>& cell, std::vector<T> const& rawData, plint order) const;

    /// Change the space and time scales of the variables in moment representation.
    virtual void rescale(std::vector<T>& rawData, T xDxInv, T xDt, plint order) const;
    virtual void rescale(int dxScale, int dtScale) {
        Dynamics<T,Descriptor>::rescale(dxScale, dtScale);
    }

/* *************** Additional moments, intended for internal use ************ */

    /// Compute order-0 moment rho-bar
    virtual T computeRhoBar(Cell<T,Descriptor> const& cell) const;

    /// Compute order-0 moment rho-bar and order-1 moment j
    virtual void computeRhoBarJ(Cell<T,Descriptor> const& cell,
                                T& rhoBar, Array<T,Descriptor<T>::d>& j) const;

    /// Compute order-0 moment rho-bar, order-1 moment j, and order-2
    ///   off-equilibrium moment PiNeq.
    virtual void computeRhoBarJPiNeq(Cell<T,Descriptor> const& cell,
                                     T& rhoBar, Array<T,Descriptor<T>::d>& j,
                                     Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq) const;

    /// Compute e-bar, which is related to the internal energy
    virtual T computeEbar(Cell<T,Descriptor> const& cell) const;

/* *************** Access to Dynamics variables, e.g. omega ***************** */

    /// Get local relaxation parameter of the dynamics
    virtual T getOmega() const;

    /// Set local relaxation parameter of the dynamics
    virtual void setOmega(T omega_);

    /// Get local value of any generic parameter
    virtual T getParameter(plint whichParameter) const;

    /// Set local value of any generic parameter
    virtual void setParameter(plint whichParameter, T value);

/* *************** Access Cell raw data through Dynamics ********************* */

    /// Access particle populations through the dynamics object.
    virtual void getPopulations(Cell<T,Descriptor> const& cell, Array<T,Descriptor<T>::q>& f) const;

    /// Access external fields through the dynamics object.
    virtual void getExternalField (
            Cell<T,Descriptor> const& cell, plint pos, plint size, T* ext ) const;

    /// Define particle populations through the dynamics object.
    virtual void setPopulations(Cell<T,Descriptor>& cell, Array<T,Descriptor<T>::q> const& f);

    /// Define external fields through the dynamics object.
    virtual void setExternalField (
            Cell<T,Descriptor>& cell, plint pos, plint size, const T* ext);

/* *************** Define macroscopic variables, e.g. on boundaries ********* */

    /// Define density, if possible. Does nothing by default.
    virtual void defineDensity(Cell<T,Descriptor>& cell, T density);

    /// Define velocity, if possible. Does nothing by default.
    virtual void defineVelocity(Cell<T,Descriptor>& cell, Array<T,Descriptor<T>::d> const& u);

    /// Define temperature, if possible. Does nothing by default.
    virtual void defineTemperature(Cell<T,Descriptor>& cell, T temperature);

    /// Define heat flux, if possible. Does nothing by default.
    virtual void defineHeatFlux(Cell<T,Descriptor>& cell, Array<T,Descriptor<T>::d> const& q);

    /// Define deviatoric stress, if possible. Does nothing by default.
    virtual void definePiNeq(Cell<T,Descriptor>& cell,
                                        Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq);

    /// Define user define moment, if possible. Does nothing by default.
    virtual void defineMoment(Cell<T,Descriptor>& cell, plint momentId, T const* value);

    /// Returns flag on whether prepareCollision() is invoked before collision or not.
    bool doesAutomaticPrepareCollision() const;

    /// Toggle flag on whether prepareCollision() is invoked before collision or not.
    void toggleAutomaticPrepareCollision(bool flag);
private:
    Dynamics<T,Descriptor>* baseDynamics;
    bool automaticPrepareCollision;
};

/// Base class for composite dynamics which pre-attributes values to distribution functions.
/** During the collision of a composite dynamics, a completion scheme is
 *  first invoked to assign new values to the particle populations. Then,
 *  another dynamics is called to fullfill collision.
 **/
template<typename T, template<typename U> class Descriptor>
class PreparePopulationsDynamics : public CompositeDynamics<T,Descriptor> {
public:
/* *************** Construction and Destruction ***************************** */
    PreparePopulationsDynamics( Dynamics<T,Descriptor>* baseDynamics_,
                                bool automaticPrepareCollision_=true );
    virtual void prepareCollision(Cell<T,Descriptor>& cell);
    virtual void completePopulations(Cell<T,Descriptor>& cell) const =0;
    virtual PreparePopulationsDynamics<T,Descriptor>* clone() const =0;
};


/// Computation of the macroscopic variables is forwarded to basic dynamics
/** For the computation of macroscopic variables, the behavior is inherited
 *  from CompositeDynamics: forward everything to base dynamics.
 */
template<typename T, template<typename U> class Descriptor>
class BulkCompositeDynamics : public PreparePopulationsDynamics<T,Descriptor> {
public:
/* *************** Construction and Destruction ***************************** */
    BulkCompositeDynamics (Dynamics<T,Descriptor>* baseDynamics_,
                           bool automaticPrepareCollision_=true);
};


/// Implementation of "full-way bounce-back" dynamics
/** This is a very popular way to implement no-slip boundary conditions,
 * because the dynamics are independent of the orientation of the boundary.
 * It is a special case, because it implements no usual LB dynamics.
 * For that reason, it derives directly from the class Dynamics.
 *
 * The code works for both 2D and 3D lattices.
 */
template<typename T, template<typename U> class Descriptor>
class BounceBack : public Dynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ******************************* */

    /// You may fix a fictitious density value on bounce-back nodes via the constructor.
    BounceBack(T rho_=T());
    BounceBack(HierarchicUnserializer& unserializer);

    /// Clone the object on its dynamic type.
    virtual BounceBack<T,Descriptor>* clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    /// Serialize the dynamics object.
    virtual void serialize(HierarchicSerializer& serializer) const;
    /// Un-Serialize the dynamics object.
    virtual void unserialize(HierarchicUnserializer& unserializer);

/* *************** Collision, Equilibrium, and Non-equilibrium ************** */

    /// Implementation of the collision step
    virtual void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics& statistics_);

    /// Implementation of the collision step, with imposed macroscopic variables
    virtual void collideExternal(Cell<T,Descriptor>& cell, T rhoBar,
                         Array<T,Descriptor<T>::d> const& j, T thetaBar, BlockStatistics& stat);

    /// Yields 0
    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;

    /// Does nothing
    virtual void regularize(Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                            T jSqr, Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq, T thetaBar=T() ) const;

/* *************** Computation of macroscopic variables ********************* */

    /// Yields fictitious density
    virtual T computeDensity(Cell<T,Descriptor> const& cell) const;
    /// Yields 0
    virtual T computePressure(Cell<T,Descriptor> const& cell) const;
    /// Yields 0
    virtual void computeVelocity( Cell<T,Descriptor> const& cell,
                                  Array<T,Descriptor<T>::d>& u ) const;
    /// Yields 0
    virtual T computeTemperature(Cell<T,Descriptor> const& cell) const;
    /// Yields 0
    virtual void computePiNeq (
        Cell<T,Descriptor> const& cell, Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq ) const;
    /// Yields 0
    virtual void computeShearStress (
        Cell<T,Descriptor> const& cell, Array<T,SymmetricTensor<T,Descriptor>::n>& stress ) const;
    /// Yields 0
    virtual void computeHeatFlux( Cell<T,Descriptor> const& cell,
                                  Array<T,Descriptor<T>::d>& q ) const;

    /// Does nothing
    virtual void computeMoment( Cell<T,Descriptor> const& cell,
                                plint momentId, T* moment ) const;

/* *************** Access to Dynamics variables, e.g. omega ***************** */

    /// Yields 0
    virtual T getOmega() const;

    /// Does nothing
    virtual void setOmega(T omega_);

/* *************** Switch between population and moment representation ****** */

    /// Yields Descriptor<T>::q + Descriptor<T>::ExternalField::numScalars.
    virtual plint numDecomposedVariables(plint order) const;

    /// Decomposed data is identical with original cell data.
    virtual void decompose(Cell<T,Descriptor> const& cell, std::vector<T>& rawData, plint order) const;

    /// Decomposed data is identical with original cell data.
    virtual void recompose(Cell<T,Descriptor>& cell, std::vector<T> const& rawData, plint order) const;

    /// Nothing happens here.
    virtual void rescale(std::vector<T>& rawData, T xDxInv, T xDt, plint order) const;
    virtual void rescale(int dxScale, int dtScale) {
        Dynamics<T,Descriptor>::rescale(dxScale, dtScale);
    }

    /// BounceBack is a boundary.
    virtual bool isBoundary() const;

    /// For BounceBack the moments of the populations have no meaning.
    virtual bool hasMoments() const;

/* *************** Additional moments, intended for internal use ************ */

    /// Yields fictitious density
    virtual T computeRhoBar(Cell<T,Descriptor> const& cell) const;

    /// Yields fictitious density and 0
    virtual void computeRhoBarJ(Cell<T,Descriptor> const& cell,
                                T& rhoBar, Array<T,Descriptor<T>::d>& j) const;

    /// Compute order-0 moment rho-bar, order-1 moment j, and order-2
    ///   off-equilibrium moment PiNeq.
    virtual void computeRhoBarJPiNeq(Cell<T,Descriptor> const& cell,
                                     T& rhoBar, Array<T,Descriptor<T>::d>& j,
                                     Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq) const;

    /// Yields 0
    virtual T computeEbar(Cell<T,Descriptor> const& cell) const;
private:
    T rho;
private:
    static int id;
};


/// Implementation of specular reflection dynamics
/** This is a very simple way to implement free-slip boundary conditions.
 * It is a special case, because it implements no usual LB dynamics.
 * For that reason, it derives directly from the class Dynamics.
 *
 * The code works for both 2D and 3D lattices.
 */
template<typename T, template<typename U> class Descriptor>
class SpecularReflection : public Dynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ******************************* */

    /// You may fix a fictitious density value on specular reflection nodes via the constructor.
    SpecularReflection(T rho_=T());
    SpecularReflection(Array<bool,Descriptor<T>::d> const& reflectOnPlaneNormalToAxis, T rho_=T());
    SpecularReflection(HierarchicUnserializer& unserializer);

    /// Clone the object on its dynamic type.
    virtual SpecularReflection<T,Descriptor>* clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    /// Serialize the dynamics object.
    virtual void serialize(HierarchicSerializer& serializer) const;
    /// Un-Serialize the dynamics object.
    virtual void unserialize(HierarchicUnserializer& unserializer);

/* *************** Collision, Equilibrium, and Non-equilibrium ************** */

    /// Implementation of the collision step
    virtual void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics& statistics_);

    /// Implementation of the collision step, with imposed macroscopic variables
    virtual void collideExternal(Cell<T,Descriptor>& cell, T rhoBar,
                         Array<T,Descriptor<T>::d> const& j, T thetaBar, BlockStatistics& stat);

    /// Yields 0
    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;

    /// Does nothing
    virtual void regularize(Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                            T jSqr, Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq, T thetaBar=T() ) const;

/* *************** Computation of macroscopic variables ********************* */

    /// Yields fictitious density
    virtual T computeDensity(Cell<T,Descriptor> const& cell) const;
    /// Yields 0
    virtual T computePressure(Cell<T,Descriptor> const& cell) const;
    /// Yields 0
    virtual void computeVelocity( Cell<T,Descriptor> const& cell,
                                  Array<T,Descriptor<T>::d>& u ) const;
    /// Yields 0
    virtual T computeTemperature(Cell<T,Descriptor> const& cell) const;
    /// Yields 0
    virtual void computePiNeq (
        Cell<T,Descriptor> const& cell, Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq ) const;
    /// Yields 0
    virtual void computeShearStress (
        Cell<T,Descriptor> const& cell, Array<T,SymmetricTensor<T,Descriptor>::n>& stress ) const;
    /// Yields 0
    virtual void computeHeatFlux( Cell<T,Descriptor> const& cell,
                                  Array<T,Descriptor<T>::d>& q ) const;

    /// Does nothing
    virtual void computeMoment( Cell<T,Descriptor> const& cell,
                                plint momentId, T* moment ) const;

/* *************** Access to Dynamics variables, e.g. omega ***************** */

    /// Yields 0
    virtual T getOmega() const;

    /// Does nothing
    virtual void setOmega(T omega_);

/* *************** Switch between population and moment representation ****** */

    /// Yields Descriptor<T>::q + Descriptor<T>::ExternalField::numScalars.
    virtual plint numDecomposedVariables(plint order) const;

    /// Decomposed data is identical with original cell data.
    virtual void decompose(Cell<T,Descriptor> const& cell, std::vector<T>& rawData, plint order) const;

    /// Decomposed data is identical with original cell data.
    virtual void recompose(Cell<T,Descriptor>& cell, std::vector<T> const& rawData, plint order) const;

    /// Nothing happens here.
    virtual void rescale(std::vector<T>& rawData, T xDxInv, T xDt, plint order) const;
    virtual void rescale(int dxScale, int dtScale) {
        Dynamics<T,Descriptor>::rescale(dxScale, dtScale);
    }

    /// SpecularReflection is a boundary.
    virtual bool isBoundary() const;

    /// For SpecularReflection the moments of the populations have no meaning.
    virtual bool hasMoments() const;

/* *************** Additional moments, intended for internal use ************ */

    /// Yields fictitious density
    virtual T computeRhoBar(Cell<T,Descriptor> const& cell) const;

    /// Yields fictitious density and 0
    virtual void computeRhoBarJ(Cell<T,Descriptor> const& cell,
                                T& rhoBar, Array<T,Descriptor<T>::d>& j) const;

    /// Compute order-0 moment rho-bar, order-1 moment j, and order-2
    ///   off-equilibrium moment PiNeq.
    virtual void computeRhoBarJPiNeq(Cell<T,Descriptor> const& cell,
                                     T& rhoBar, Array<T,Descriptor<T>::d>& j,
                                     Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq) const;

    /// Yields 0
    virtual T computeEbar(Cell<T,Descriptor> const& cell) const;
private:
    T rho;
    Array<int,Descriptor<T>::d> reflection;
private:
    static int id;
};


/// Implementation of "dead dynamics" which does nothing
template<typename T, template<typename U> class Descriptor>
class NoDynamics : public Dynamics<T,Descriptor> {
public:
/* *************** Construction / Destruction ******************************* */
    NoDynamics(T rho_=T());
    NoDynamics(HierarchicUnserializer& unserializer);

    /// Clone the object on its dynamic type.
    virtual NoDynamics<T,Descriptor>* clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    /// Serialize the dynamics object.
    virtual void serialize(HierarchicSerializer& serializer) const;
    /// Un-Serialize the dynamics object.
    virtual void unserialize(HierarchicUnserializer& unserializer);

/* *************** Collision, Equilibrium, and Non-equilibrium ************** */

    /// Does nothing
    virtual void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics& statistics_);

    /// Does nothing
    virtual void collideExternal(Cell<T,Descriptor>& cell, T rhoBar,
                         Array<T,Descriptor<T>::d> const& j, T thetaBar, BlockStatistics& stat);

    /// Yields 0
    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;

    /// Does nothing
    virtual void regularize(Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                            T jSqr, Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq, T thetaBar=T() ) const;

/* *************** Computation of macroscopic variables ********************* */

    /// Yields 1
    virtual T computeDensity(Cell<T,Descriptor> const& cell) const;
    /// Yields 0
    virtual T computePressure(Cell<T,Descriptor> const& cell) const;
    /// Yields 0
    virtual void computeVelocity( Cell<T,Descriptor> const& cell,
                                  Array<T,Descriptor<T>::d>& u ) const;
    /// Yields 0
    virtual T computeTemperature(Cell<T,Descriptor> const& cell) const;
    /// Yields 0
    virtual void computePiNeq (
        Cell<T,Descriptor> const& cell, Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq ) const;
    /// Yields 0
    virtual void computeShearStress (
        Cell<T,Descriptor> const& cell, Array<T,SymmetricTensor<T,Descriptor>::n>& stress ) const;
    /// Yields 0
    virtual void computeHeatFlux( Cell<T,Descriptor> const& cell,
                                  Array<T,Descriptor<T>::d>& q ) const;

    /// Does nothing
    virtual void computeMoment( Cell<T,Descriptor> const& cell,
                                plint momentId, T* moment ) const;

/* *************** Access to Dynamics variables, e.g. omega ***************** */

    /// Yields 0
    virtual T getOmega() const;

    /// Does nothing
    virtual void setOmega(T omega_);

/* *************** Switch between population and moment representation ****** */

    /// Yields Descriptor<T>::q + Descriptor<T>::ExternalField::numScalars.
    virtual plint numDecomposedVariables(plint order) const;

    /// Decomposed data is identical with original cell data.
    virtual void decompose(Cell<T,Descriptor> const& cell, std::vector<T>& rawData, plint order) const;

    /// Decomposed data is identical with original cell data.
    virtual void recompose(Cell<T,Descriptor>& cell, std::vector<T> const& rawData, plint order) const;

    /// Nothing happens here.
    virtual void rescale(std::vector<T>& rawData, T xDxInv, T xDt, plint order) const;
    virtual void rescale(int dxScale, int dtScale) {
        Dynamics<T,Descriptor>::rescale(dxScale, dtScale);
    }

    /// For NoDynamics the moments of the populations have no meaning.
    virtual bool hasMoments() const;

/* *************** Additional moments, intended for internal use ************ */

    /// Yields rho=1
    virtual T computeRhoBar(Cell<T,Descriptor> const& cell) const;

    /// Yields rho=1, j=0
    virtual void computeRhoBarJ(Cell<T,Descriptor> const& cell,
                                T& rhoBar, Array<T,Descriptor<T>::d>& j) const;

    /// Compute order-0 moment rho-bar, order-1 moment j, and order-2
    ///   off-equilibrium moment PiNeq.
    virtual void computeRhoBarJPiNeq(Cell<T,Descriptor> const& cell,
                                     T& rhoBar, Array<T,Descriptor<T>::d>& j,
                                     Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq) const;

    /// Yields 0
    virtual T computeEbar(Cell<T,Descriptor> const& cell) const;
private:
    T rho;
private:
    static int id;
};

/// Get all the IDs of the dynamics which are part of a composite dynamics construct.
template<typename T, template<typename U> class Descriptor>
void constructIdChain(Dynamics<T,Descriptor> const& dynamics, std::vector<int>& chain);

/// Access the bottom-most dynamics in a hierarchy of composite dynamics.
template<typename T, template<typename U> class Descriptor>
Dynamics<T,Descriptor> const& getBottomMostDynamics(Dynamics<T,Descriptor> const& dynamics);

/// Make a clone of the full composite dynamics, and replace the bottom-most element in the
///   chain (i.e. the one which defines the underlying model, such as BGK).
template<typename T, template<typename U> class Descriptor>
Dynamics<T,Descriptor>* cloneAndReplaceBottomDynamics(Dynamics<T,Descriptor> const& dynamics,
                                                      Dynamics<T,Descriptor>* newBottom);

/// Make a clone of the full composite dynamics, and insert a new dynamics at the top of the
///   chain. If the current top element is not composeable, it is either removed, or the new
///   dynamics is inserted right behind, depending on whether the new dynamics is composeable
///   or not.
template<typename T, template<typename U> class Descriptor>
Dynamics<T,Descriptor>* cloneAndInsertAtTopDynamics(Dynamics<T,Descriptor> const& dynamics,
                                                    CompositeDynamics<T,Descriptor>* newTop);

/// Clone dynamics object, but remove all components which represent boundary (isBoundary()
///   method). Exception is the bottom-most dynamics object which is never removed.
template<typename T, template<typename U> class Descriptor>
Dynamics<T,Descriptor>* removeBoundaryComponents(Dynamics<T,Descriptor> const& dynamics);

/// Serialize the whole dynamics chain into a byte-stream, and append it to
///   the vector data.
template<typename T, template<typename U> class Descriptor>
void serialize(Dynamics<T,Descriptor> const& dynamics, std::vector<char>& data);

/// Serialize all dynamics objects from the list into a single byte-stream.
template<typename T, template<typename U> class Descriptor>
void serialize(std::vector<Dynamics<T,Descriptor>*> const& dynamics, std::vector<char>& data);

/// Unserialize a single dynamics object. Attention: the type of the dynamics
///   must match the serialized content. Otherwise, the result is undefined,
///   except in debug mode where you get an error message. Optionally, you can
///   indicate the position in the byte-stream from which to unserialize.
///   Returns the position from which the next dynamics object can be unserialized,
///   if there is any.
template<typename T, template<typename U> class Descriptor>
pluint unserialize( Dynamics<T,Descriptor>& dynamics,
                    std::vector<char> const& data,
                    pluint serializerPos=0 );

/// Unserialize all data into newly generated dynamics objects.
template<typename T, template<typename U> class Descriptor>
void generateAndUnserializeDynamics (
        std::vector<char> const& data,
        std::vector<Dynamics<T,Descriptor>*>& dynamics);

}  // namespace plb

#endif  // DYNAMICS_H
