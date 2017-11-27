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
 * Helper functions for domain initialization -- header file.
 */
#ifndef PLB_MATH_H
#define PLB_MATH_H

#include "core/globalDefs.h"
#include <cmath>

namespace plb {

// Templatize computation of the power, because
//   the integer case must be treated separately.

template<typename T>
inline T customPower(T a, T b) {
    return pow(a,b);
}

template<>
inline int customPower(int a, int b) {
    double result = pow((double)a,(double)b);
    return (int)(result+0.5);
}

template<typename T>
inline void customInPlacePower(T& a, T b) {
    a = pow(a,b);
}

template<>
inline void customInPlacePower(int& a, int b) {
    double result = pow((double)a,(double)b);
    a = (int)(result+0.5);
}

}  // namespace plb

#endif  // PLB_MATH_H
