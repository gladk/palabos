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

#ifndef ZOU_HE_DYNAMICS_H
#define ZOU_HE_DYNAMICS_H

#include "core/globalDefs.h"
#include "boundaryCondition/boundaryDynamics.h"

namespace plb {

/**
* Implementation of Zou-He boundary condition following
* the paper from Zou and He. This implementation is lattice independent.
* The implementation follow the idea proposed in the paper
* Qisu Zou, Xiaoyi He, 
* "On pressure and velocity boundary conditions for the lattice Boltzmann BGK model",
* Phys. Fluids , (1997), Volume 9, Issue 6, pp. 1591-1598
*/
template<typename T, template<typename U> class Descriptor, int direction, int orientation>
class ZouHeVelocityDynamics : public VelocityDirichletBoundaryDynamics<T,Descriptor,direction,orientation>
{
public:
    /// Constructor
    ZouHeVelocityDynamics(Dynamics<T,Descriptor>* baseDynamics, bool automaticPrepareCollision = true);
    ZouHeVelocityDynamics(HierarchicUnserializer& unserializer);
    /// Clone the object on its dynamic type.
    virtual ZouHeVelocityDynamics<T, Descriptor, direction, orientation>* clone() const;
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Execute completion scheme before base collision
    virtual void completePopulations(Cell<T,Descriptor>& cell) const;
private:
    static int id;
};

/**
* Implementation of Zou-He boundary condition following
* the paper from Zou and He. This implementation is lattice independent.
* The implementation follow the idea proposed in the paper
* Qisu Zou, Xiaoyi He, 
* "On pressure and velocity boundary conditions for the lattice Boltzmann BGK model",
* Phys. Fluids , (1997), Volume 9, Issue 6, pp. 1591-1598
*/
template<typename T, template<typename U> class Descriptor, int direction, int orientation>
class ZouHePressureDynamics : public DensityDirichletBoundaryDynamics<T,Descriptor,direction,orientation>
{
public:
    /// Constructor
    ZouHePressureDynamics(Dynamics<T,Descriptor>* baseDynamics, bool automaticPrepareCollision = true);
    ZouHePressureDynamics(HierarchicUnserializer& unserializer);
    /// Clone the object on its dynamic type.
    virtual ZouHePressureDynamics<T, Descriptor, direction, orientation>* clone() const;
    /// Return a unique ID for this class.
    virtual int getId() const;
    /// Execute completion scheme before base collision
    virtual void completePopulations(Cell<T,Descriptor>& cell) const;
private:
    static int id;
};

}  // namespace plb

#endif  // ZOU_HE_DYNAMICS_H
