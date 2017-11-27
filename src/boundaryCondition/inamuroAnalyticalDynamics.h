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

#ifndef INAMURO_ANALYTICAL_DYNAMICS_H
#define INAMURO_ANALYTICAL_DYNAMICS_H

#include "core/globalDefs.h"
#include "boundaryCondition/boundaryDynamics.h"

namespace plb {

/**
* Implementation of Inamuro velocity boundary condition following
 * the paper
 * "A non-slip boundary condition for lattice Boltzmann simulations",
 * Inamuro, Takaji; Yoshino, Masato; Ogino, Fumimaru, (1995). 
 * This implementation works for the D2Q9 Lattice only.
*/
template<typename T, template<typename U> class Descriptor, int direction, int orientation>
class InamuroAnalyticalVelocityDynamics : public VelocityDirichletBoundaryDynamics<T,Descriptor,direction,orientation>
{
public:
    /// Constructor
    InamuroAnalyticalVelocityDynamics(Dynamics<T,Descriptor>* baseDynamics,
                                      bool automaticPrepareCollision = true);
    /// Clone the object on its dynamic type.
    virtual InamuroAnalyticalVelocityDynamics<T, Descriptor, direction, orientation>* clone() const;
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Execute completion scheme before base collision
    virtual void completePopulations(Cell<T,Descriptor>& cell) const;
private:
    static int id;
};

/**
* Implementation of Inamuro pressure boundary condition following
 * the paper
 * "A non-slip boundary condition for lattice Boltzmann simulations",
 * Inamuro, Takaji; Yoshino, Masato; Ogino, Fumimaru, (1995). 
 * This implementation works for the D2Q9 Lattice only.
*/
template<typename T, template<typename U> class Descriptor, int direction, int orientation>
class InamuroAnalyticalPressureDynamics : public DensityDirichletBoundaryDynamics<T,Descriptor,direction,orientation>
{
public:
    /// Constructor
    InamuroAnalyticalPressureDynamics(Dynamics<T,Descriptor>* baseDynamics,
                                      bool automaticPrepareCollision = true);
    /// Clone the object on its dynamic type.
    virtual InamuroAnalyticalPressureDynamics<T, Descriptor, direction, orientation>* clone() const;
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Execute completion scheme before base collision
    virtual void completePopulations(Cell<T,Descriptor>& cell) const;
private:
    static int id;
};

}  // namespace plb

#endif  // INAMURO_ANALYTICAL_DYNAMICS_H
