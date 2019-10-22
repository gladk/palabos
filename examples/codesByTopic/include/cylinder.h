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

#ifndef CYLINDER_H
#define CYLINDER_H

#include "palabos2D.h"

template<typename T>
class CylinderShapeDomain2D : public plb::DomainFunctional2D {
public:
    CylinderShapeDomain2D(plb::plint cx_, plb::plint cy_, plb::plint radius)
        : cx(cx_),
          cy(cy_),
          radiusSqr(plb::util::sqr(radius))
    { }
    virtual bool operator() (plb::plint iX, plb::plint iY) const {
        return plb::util::sqr(iX-cx) + plb::util::sqr(iY-cy) <= radiusSqr;
    }
    virtual CylinderShapeDomain2D<T>* clone() const {
        return new CylinderShapeDomain2D<T>(*this);
    }
private:
    plb::plint cx;
    plb::plint cy;
    plb::plint radiusSqr;
};

/// Convert pressure to density according to ideal gas law
template<typename T, template<typename U> class Descriptor>
void createCylinder( plb::MultiBlockLattice2D<T,Descriptor>& lattice,
                     plb::plint cx, plb::plint cy, plb::plint radius );

#endif  // CYLINDER_H
