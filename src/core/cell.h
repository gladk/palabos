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
 * Definition of a LB cell -- header file.
 */
#ifndef CELL_H
#define CELL_H

#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include "latticeBoltzmann/externalFields.h"
#include "core/dynamics.h"
#include "core/array.h"
#include "core/hierarchicSerializer.h"

namespace plb {

/// Helper class: allocation of memory for external fields in a cell
template<typename T, typename ExternalField>
class ExternalFieldArray {
public:
    T* get(plint index) {
        PLB_PRECONDITION( index < ExternalField::numScalars );
        return data+index;
    }
    T const* get(plint index) const {
        PLB_PRECONDITION( index < ExternalField::numScalars );
        return data+index;
    }
private:
    T data[ExternalField::numScalars];
};

/// Specialization of ExternalFieldArray, when no external field is present
template<typename T>
class ExternalFieldArray<T,descriptors::NoExternalField> {
public:
    T* get(pluint index) {
        PLB_PRECONDITION( false );
        static T data = T();
        return &data;
    }
    T const* get(pluint index) const {
        PLB_PRECONDITION( false );
        static T data = T();
        return &data;
    }
};

/// Helper class to evaluate statically the full memory requirement of a cell.
template<typename T, template<typename U> class Descriptor>
struct CellInfo {
    static const int n = Descriptor<T>::numPop + Descriptor<T>::ExternalField::numScalars;
};

/// A LB lattice cell.
/** A cell contains the q values of the distribution functions f on
 * one lattice point, as well as a pointer to the dynamics of the
 * cell. Thanks to this pointer, one can have a space dependend de-
 * finition of the dynamics. This mechanism is useful e.g. for the
 * implementation of boundary conditions, or an inhomogeneous body
 * force.
 *
 * The dynamics object is not owned by the class, it is not
 * destructed in the Cell destructor.
 *
 * This class is not intended to be derived from.
 */
template<typename T, template<typename U> class Descriptor>
class Cell {
public:
    /// Additional per-cell scalars for external fields, e.g. forces
    typedef ExternalFieldArray<T, typename Descriptor<T>::ExternalField> External;
public:
    /// Default constructor.
    Cell();
    /// Constructor, to be used whenever possible.
    Cell(Dynamics<T,Descriptor>* dynamics_);
public:
    /// Read-write access to distribution functions.
    /** \param iPop index of the accessed distribution function */
    T& operator[](plint iPop) {
        PLB_PRECONDITION( iPop < Descriptor<T>::numPop );
        return f[iPop];
    }
    /// Read-only access to distribution functions.
    /** \param iPop index of the accessed distribution function */
    T const& operator[](plint iPop) const {
        PLB_PRECONDITION( iPop < Descriptor<T>::numPop );
        return f[iPop];
    }
    /// Another way to get direct access to the f's, as in operator[]
    Array<T,Descriptor<T>::numPop>& getRawPopulations() {
        return f;
    }
    /// Another way to get direct, const access to the f's, as in operator[]
    Array<T,Descriptor<T>::numPop> const& getRawPopulations() const {
        return f;
    }
    /// Attribute all f-values from another cell to the present one.
    /** \return a reference to *this
     */
    Cell<T,Descriptor>& attributeF(Cell<T,Descriptor> const& rhs) {
        f = rhs.getRawPopulations();
        return *this;
    };
    /// Attribute all f-values and external scalars from another cell to the present one.
    /** \return a reference to *this
     * This is similar to the assignment operator operator= (which is
     * created by default and copies all the f's, as well as the dynamics
     * object), except that the dynamics object is not copied.
     */
    Cell<T,Descriptor>& attributeValues(Cell<T,Descriptor> const& rhs) {
        attributeF(rhs);
        for (plint iExt=0; iExt < Descriptor<T>::ExternalField::numScalars; ++iExt) {
            *external.get(iExt) = *rhs.external.get(iExt);
        }
        return *this;
    };
    /// Get a pointer to an external field
    T* getExternal(plint offset) {
        PLB_PRECONDITION( offset < Descriptor<T>::ExternalField::numScalars );
        return external.get(offset);
    }
    /// Get a const pointer to an external field
    T const* getExternal(plint offset) const {
        PLB_PRECONDITION( offset < Descriptor<T>::ExternalField::numScalars );
        return external.get(offset);
    }
    /// Get a reference to non-modifiable dynamics
    Dynamics<T,Descriptor> const& getDynamics() const;
    /// Get a a reference to dynamics
    Dynamics<T,Descriptor>& getDynamics();
    /// Request whether this cell does statistics measurements
    bool takesStatistics() const {
        return takesStat;
    }
    /// Specify whether this cell does statistics measurements
    void specifyStatisticsStatus(bool status) {
        takesStat = status;
    }
private:
    /// You can't use this method. Use BlockLattice::attributeDynamics instead.
    /** This is one of the rare cases of a method accepting a pointer but
     *  not managing memory itself. Memory of the dynamics is handled by
     *  the BlockLattice.
     */
    void attributeDynamics(Dynamics<T,Descriptor>* dynamics_);
    // Declare the BlockLatticeXD as a friend, to enable access to attributeDynamics.
    template<typename T_, template<typename U_> class Descriptor_> friend class BlockLattice2D;
    template<typename T_, template<typename U_> class Descriptor_> friend class BlockLattice3D;
#ifdef PLB_MPI_PARALLEL
    template<typename T_, template<typename U_> class Descriptor_> friend class ParallelCellAccess2D;
    template<typename T_, template<typename U_> class Descriptor_> friend class ParallelCellAccess3D;
#endif

// The following helper functions forward the function call
// to the Dynamics object
public:
    /// Apply LB collision to the cell according to local dynamics.
    void collide(BlockStatistics& statistics) {
        PLB_PRECONDITION( dynamics );
        dynamics->collide(*this, statistics);
    }

    /// Compute equilibrium distribution function
    T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                         T jSqr, T thetaBar=T()) const
    {
        PLB_PRECONDITION( dynamics );
        return dynamics->computeEquilibrium(iPop, rhoBar, j, jSqr, thetaBar);
    }

    /// Re-compute particle populations from the leading moments
    void regularize(T rhoBar, Array<T,Descriptor<T>::d> const& j, T jSqr,
                    Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq, T thetaBar=T() )
    {
        PLB_PRECONDITION( dynamics );
        dynamics->regularize(*this, rhoBar, j, jSqr, PiNeq, thetaBar);
    }

    /// Compute particle density on the cell.
    /** \return particle density
     */
    T computeDensity() const {
        PLB_PRECONDITION( dynamics );
        return dynamics->computeDensity(*this);
    }

    /// Compute pressure on the cell.
    /** \return pressure
     */
    T computePressure() const {
        PLB_PRECONDITION( dynamics );
        return dynamics->computePressure(*this);
    }

    /// Compute fluid velocity on the cell.
    /** \param u fluid velocity
     */
    void computeVelocity(Array<T,Descriptor<T>::d>& u) const {
        PLB_PRECONDITION( dynamics );
        dynamics->computeVelocity(*this, u);
    }

    /// Compute Temperature on the cell.
    /** \return temperature
     */
    T computeTemperature() const {
        PLB_PRECONDITION( dynamics );
        return dynamics->computeTemperature(*this);
    }

    /// Compute the "off-equilibrium part of Pi"
    /** \param PiNeq stress tensor */
    void computePiNeq (
            Array<T,SymmetricTensor<T,Descriptor>::n>& PiNeq ) const
    {
        PLB_PRECONDITION( dynamics );
        dynamics->computePiNeq(*this, PiNeq);
    }
    
    /// Compute the deviatoric stress tensor
    /** \param PiNeq stress tensor */
    void computeShearStress (
        Array<T,SymmetricTensor<T,Descriptor>::n>& stress ) const
        {
            PLB_PRECONDITION( dynamics );
            dynamics->computeShearStress(*this, stress);
        }

    /// Compute heat flux on the cell.
    /** \param q heat flux
     */
    void computeHeatFlux(Array<T,Descriptor<T>::d>& q) const {
        PLB_PRECONDITION( dynamics );
        dynamics->computeHeatFlux(*this, q);
    }

    /// Compute user-defined moment on the cell.
    /** \param momentId Identifier for the moment.
     *  \param moment return value: computed moment.
     */
    void computeMoment(plint momentId, T* moment) {
        PLB_PRECONDITION( dynamics );
        dynamics->computeMoment(*this, momentId, moment);
    }

    /// Access particle populations through the dynamics object.
    /** This method is similar to operator[]: it delivers the
     * value of the particle populations. This time, those values
     * are however computed through a virtual call to the dynamics
     * object.
     */
    void getPopulations(Array<T,Descriptor<T>::numPop>& f) const {
        PLB_PRECONDITION( dynamics );
        dynamics->getPopulations(*this, f);
    }

    /// Access external fields through the dynamics object.
    /** This method is similar to getExternal(): it delivers the
     * value of the external fields. This time, those values
     * are however computed through a virtual call to the dynamics
     * object.
     */
    void getExternalField(plint pos, plint size, T* ext) const {
        PLB_PRECONDITION( dynamics );
        dynamics->getExternalField(*this, pos, size, ext);
    }   

    /// Set particle density on the cell.
    /** \param rho particle density
     */
    void defineDensity(T rho) {
        PLB_PRECONDITION( dynamics );
        dynamics->defineDensity(*this, rho);
    }

    /// Set fluid velocity on the cell.
    /** \param u fluid velocity
     */
    void defineVelocity(Array<T,Descriptor<T>::d> const& u) {
        PLB_PRECONDITION( dynamics );
        dynamics->defineVelocity(*this, u);
    }

    /// Set temperature on the cell.
    /** \param temperature Temperature
     */
    void defineTemperature(T temperature) {
        PLB_PRECONDITION( dynamics );
        dynamics->defineTemperature(*this, temperature);
    }

    /// Set heat flux on the cell.
    /** \param q heat flux
     */
    void defineHeatFlux(Array<T,Descriptor<T>::d> const& q) {
        PLB_PRECONDITION( dynamics );
        dynamics->defineHeatFlux(*this, q);
    }

    /// Set components of the deviatoric stress tensor on the cell.
    /** \param PiNeq Deviatoric stress tensor
     * */
    void definePiNeq (
            Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq )
    {
        PLB_PRECONDITION( dynamics );
        dynamics->definePiNeq(*this, PiNeq);
    }

    /// Set generic moment on the cell.
    /** \param momentId Identifier of user-defined moment.
     *  \param value Value of the moment.
     */
    void defineMoment(plint momentId, T const* value)
    {
        PLB_PRECONDITION( dynamics );
        dynamics->defineMoment(*this, momentId, value);
    }

    /// Define particle populations through the dynamics object.
    /** This method is similar to operator[]: it modifies the
     * value of the particle populations. This time, those values
     * are however accessed through a virtual call to the dynamics
     * object.
     */
    void setPopulations(Array<T,Descriptor<T>::numPop> const& f) {
        PLB_PRECONDITION( dynamics );
        dynamics->setPopulations(*this, f);
    }
    /// Define external fields through the dynamics object.
    /** This method is similar to getExternal(): it accesses the
     * value of the external fields. This time, those values
     * are however accessed through a virtual call to the dynamics
     * object.
     */
    void setExternalField(plint pos, plint size, const T* ext) {
        PLB_PRECONDITION( dynamics );
        dynamics->setExternalField(*this, pos, size, ext);
    }
    /// Revert ("bounce-back") the distribution functions.
    void revert();
    void serialize(char* data) const;
    void unSerialize(char const* data);
private:
    void iniPop();
    void iniExternal();
private:
    Array<T,Descriptor<T>::numPop> f;         ///< distribution functions
    External                       external;  ///< external scalars
    bool                           takesStat; ///< is statistics taken?
    Dynamics<T,Descriptor>*        dynamics;  ///< local LB dynamics
};

template<typename T, template<typename U> class Descriptor>
void iniCellAtEquilibrium(Cell<T,Descriptor>& cell, T density, Array<T,Descriptor<T>::d> const& velocity);

template<typename T, template<typename U> class Descriptor>
void iniCellAtEquilibrium(Cell<T,Descriptor>& cell, T density, Array<T,Descriptor<T>::d> const& velocity, T temperature);

}  // namespace plb

#endif  // CELL_H
