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

#ifndef LINEAR_ALGEBRA_HH
#define LINEAR_ALGEBRA_HH

#include "algorithm/linearAlgebra.h"

namespace plb {

template<typename T>
void gramSchmidt(Array<T,3> const& v1Unit, Array<T,3>& v2Unit, Array<T,3>& v3Unit) {
    Array<T,3> e1((T)1,T(),T()), e2(T(),(T)1,T()), e3(T(),T(),(T)1);
    T v11 = dot(v1Unit,e1), av11 = fabs(v11);
    T v12 = dot(v1Unit,e2), av12 = fabs(v12);
    T v13 = dot(v1Unit,e3), av13 = fabs(v13);

    if (av11 > av12) {
        if (av11 > av13) {  // e1 best aligned with v1Unit.
            v2Unit = e2 - v12*v1Unit;
            v3Unit = e3 - v13*v1Unit;
        }
        else { // e3 best aligned with v1Unit.
            v2Unit = e1 - v11*v1Unit;
            v3Unit = e2 - v12*v1Unit;
        }
    }
    else { // av12 >= av11
        if (av12 > av13) {  // e2 best aligned with v1Unit.
            v2Unit = e1 - v11*v1Unit;
            v3Unit = e3 - v13*v1Unit;
        }
        else { // e3 best aligned with v1Unit.
            v2Unit = e1 - v11*v1Unit;
            v3Unit = e2 - v12*v1Unit;
        }
    }

    v2Unit /= norm(v2Unit);
    v3Unit -= dot(v2Unit,v3Unit)*v2Unit;
    v3Unit /= norm(v3Unit);
}

} // namespace plb

#endif  // LINEAR_ALGEBRA_HH

