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

/***
 * This file was added to the Palabos library by Helen Morrison (2014)
 *
 * It implements the dynamics based on
 * Karlin, Bösch & Chikatamarla (2014), Phys. Rev. E 90, 031302(R)
 *
 */

#ifndef KBCDYNAMICS_H
#define KBCDYNAMICS_H

#include "core/globalDefs.h"
#include "basicDynamics/isoThermalDynamics.h"

namespace plb {

template<typename T, template<typename U> class Descriptor> class Cell;

/************* MODEL BASED ON GIBB'S PRINCIPLE ********************************/
/// Implementation of the entropic collision step
template<typename T, template<typename U> class Descriptor>
class KBCDynamics : public IsoThermalBulkDynamics<T,Descriptor>
{
public:
/* *************** Construction / Destruction ************************ */
    KBCDynamics(T omega_);
    KBCDynamics(HierarchicUnserializer& unserializer);

    /// Clone the object on its dynamic type.
    virtual KBCDynamics<T,Descriptor>* clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    /// TODO: check if this has to be true.
    virtual bool isEntropic() const{
        return false;
    }

/* *************** Collision and Equilibrium ************************* */
    /// Implementation of the collision step
    virtual void collide(Cell<T,Descriptor>& cell,
                         BlockStatistics& statistics_);

    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(plint iPop, T rhoBar, Array<T,Descriptor<T>::d> const& j,
                                 T jSqr, T thetaBar=T()) const;
private:
    T computeGamma(Array<T,Descriptor<T>::q> const& deltaS,
                   Array<T,Descriptor<T>::q> const& deltaH,
                   Array<T,Descriptor<T>::q> const& fEq);

    static int id;
};

}  // namespace plb

#endif // KBCDYNAMICS_H
