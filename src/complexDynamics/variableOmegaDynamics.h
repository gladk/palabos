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

/* Orestis Malaspinas designed some of the classes and concepts contained
 * in this file. */

#ifndef VARIABLE_OMEGA_DYNAMICS_H
#define VARIABLE_OMEGA_DYNAMICS_H

#include "core/globalDefs.h"
#include "core/dynamics.h"

namespace plb {

/// A dynamics for space-dependent relaxation parameter, generic with respect to base dynamics.
/** Important: The philosophy here is that the actual relaxation parameter
 *  used during collision is omega=omega0+deltaOmega. When you refer to an object "o" of
 *  type VariableOmegaDynamics, the "omega" used in "o.setOmega()" and "o.getOmega()" is
 *  omega0, whereas the "omega" in "o.getBaseDynamics().setOmega()" and in
 *  "o.getBaseDynamics().getOmega()" is omega0+deltaOmega.
 **/
template<typename T, template<typename U> class Descriptor>
class VariableOmegaDynamics : public CompositeDynamics<T,Descriptor> {
public:
    VariableOmegaDynamics(Dynamics<T,Descriptor>* baseDynamics_,
                          bool automaticPrepareCollision_=true);
    virtual void prepareCollision(Cell<T,Descriptor>& cell);
    virtual T getOmegaFromCell(Cell<T,Descriptor> const& cell) const =0;
};

/// A dynamics for relaxation parameter dependent on off-equilibrium stress, generic with respect to base dynamics.
/** Read the comments of VariableOmegaDynamics to understand the meaning of omega
 *  in setOmega and getOmega.
 **/
template<typename T, template<typename U> class Descriptor>
class OmegaFromPiDynamics : public VariableOmegaDynamics<T,Descriptor> {
public:
    OmegaFromPiDynamics(Dynamics<T,Descriptor>* baseDynamics_,
                        bool automaticPrepareCollision_=true);
    virtual T getOmegaFromCell(Cell<T,Descriptor> const& cell) const;
    virtual T getOmegaFromPiAndRhoBar(Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq, T rhoBar) const =0;
};

} // namespace plb

#endif  // VARIABLE_OMEGA_DYNAMICS_H

