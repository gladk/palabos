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

/* Main author: Daniel Lagrava
 **/

/** \file
 * Dynamics classes used to implement grid refinement -- header file.
 */

#ifndef GRID_REFINEMENT_DYNAMICS_H
#define GRID_REFINEMENT_DYNAMICS_H

#include "core/globalDefs.h"
#include "core/blockStatistics.h"
#include "boundaryCondition/boundaryDynamics.h"
#include <vector>

namespace plb {

/// Dynamics attributed to boundary cell of fine lattice with grid refinement
/** This Dynamics
 *  - Stores populations at time t (variable t0)
 *  - Stores populations at time t+1 (variable t1)
 *  - During collision, interpolates between t0 and t1, and then
 *    executes ordinary collision
 */
template<typename T, template<typename U> class Descriptor>
class FineGridBoundaryDynamics : public BoundaryCompositeDynamics<T,Descriptor> {
public:
    /// Constructor
    /** \param referenceLattice_ From the reference lattice, the FineGridBoundaryDynamics
     *                           determines the value of the current iteration.
     *  \param numTimeSteps_ Number of iteration steps leading from time t0 to time t1
     */
    FineGridBoundaryDynamics (
            Dynamics<T,Descriptor>* baseDynamics_,
            TimeCounter const& timeCounter_,
            plint numTimeSteps_,
            plint orderOfDecomposition_ );
    FineGridBoundaryDynamics(HierarchicUnserializer& unserializer);
    virtual FineGridBoundaryDynamics<T,Descriptor>* clone() const;
    bool isComposeable() const;
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Serialize the dynamics object.
    virtual void serialize(HierarchicSerializer& serializer) const;
    /// Un-Serialize the dynamics object.
    virtual void unserialize(HierarchicUnserializer& unserializer);
    /// Execute completion scheme before base collision
    virtual void completePopulations(Cell<T,Descriptor>& cell) const;
    std::vector<T>& getDecomposedValues(plint whichTime);
    std::vector<T> const& getDecomposedValues(plint whichTime) const;
private:
    TimeCounter defaultTimeCounter;
    TimeCounter const& timeCounter;
    plint numTimeSteps;
    plint orderOfDecomposition;
    std::vector<T> decomposedValuesT0;
    std::vector<T> decomposedValuesT1;
private:
    static int id;
};


}  // namespace plb

#endif  // GRID_REFINEMENT_DYNAMICS_H
