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

#ifndef INTERPOLATION_HELPER_H
#define INTERPOLATION_HELPER_H

#include "core/globalDefs.h"
#include "core/geometry3D.h"
#include "multiGrid/gridRefinement.h"
#include <vector>

namespace plb {

/// bicuadratic asymetrical interpolation using 9 coarse points
template<typename T>
std::vector<T> cornerInterpolation(T f[3][3]);

template<typename T>
std::vector<T> helperCornerInterpolation(T f[3][3]);

/// copy the decomposed populations to the given cell
template<typename T, template<typename U> class Descriptor>
void copyPopulations(std::vector<T>& decomposedValues, Cell<T,Descriptor>& cell);

template<typename T>
std::vector<T> asymetricCubicInterpolation(T f[4][3]);

template<typename T>
std::vector<T> symetricCubicInterpolation(T f[4][4]);

} // namespace plb

#endif // INTERPOLATION_HELPER_H

