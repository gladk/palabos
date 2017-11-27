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

#ifndef NEXT_NEIGHBORS_3D_HH
#define NEXT_NEIGHBORS_3D_HH

#include "offLattice/nextNeighbors3D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "latticeBoltzmann/externalFieldAccess.h"
#include <cmath>

namespace plb {

template <typename T>
const int NextNeighbor<T>::numNeighbors = 26;

template <typename T>
const int NextNeighbor<T>::c[26][3] =
{             { 0, 0, 1}, { 0, 0,-1},
  { 0, 1, 0}, { 0, 1, 1}, { 0, 1,-1},
  { 0,-1, 0}, { 0,-1, 1}, { 0,-1,-1},

  { 1, 0, 0}, { 1, 0, 1}, { 1, 0,-1},
  { 1, 1, 0}, { 1, 1, 1}, { 1, 1,-1},
  { 1,-1, 0}, { 1,-1, 1}, { 1,-1,-1},

  {-1, 0, 0}, {-1, 0, 1}, {-1, 0,-1},
  {-1, 1, 0}, {-1, 1, 1}, {-1, 1,-1},
  {-1,-1, 0}, {-1,-1, 1}, {-1,-1,-1}
};

template <typename T>
const T NextNeighbor<T>::d1 = (T) 1;
template <typename T>
const T NextNeighbor<T>::d2 = std::sqrt((T) 2);
template <typename T>
const T NextNeighbor<T>::d3 = std::sqrt((T) 3);

template <typename T>
const T NextNeighbor<T>::d[26] =
{                     (T) 1,            (T) 1,
    (T) 1,            std::sqrt((T) 2), std::sqrt((T) 2),
    (T) 1,            std::sqrt((T) 2), std::sqrt((T) 2),
    (T) 1,            std::sqrt((T) 2), std::sqrt((T) 2),
    std::sqrt((T) 2), std::sqrt((T) 3), std::sqrt((T) 3),
    std::sqrt((T) 2), std::sqrt((T) 3), std::sqrt((T) 3),
    (T) 1,            std::sqrt((T) 2), std::sqrt((T) 2),
    std::sqrt((T) 2), std::sqrt((T) 3), std::sqrt((T) 3),
    std::sqrt((T) 2), std::sqrt((T) 3), std::sqrt((T) 3) };

template <typename T>
const T NextNeighbor<T>::id1 = (T) 1;
template <typename T>
const T NextNeighbor<T>::id2 = (T) 1 / std::sqrt((T) 2);
template <typename T>
const T NextNeighbor<T>::id3 = (T) 1 / std::sqrt((T) 3);

template <typename T>
const T NextNeighbor<T>::invD[26] =
{                             (T) 1,                    (T) 1,
    (T) 1,                    (T) 1 / std::sqrt((T) 2), (T) 1 / std::sqrt((T) 2),
    (T) 1,                    (T) 1 / std::sqrt((T) 2), (T) 1 / std::sqrt((T) 2),
    (T) 1,                    (T) 1 / std::sqrt((T) 2), (T) 1 / std::sqrt((T) 2),
    (T) 1 / std::sqrt((T) 2), (T) 1 / std::sqrt((T) 3), (T) 1 / std::sqrt((T) 3),
    (T) 1 / std::sqrt((T) 2), (T) 1 / std::sqrt((T) 3), (T) 1 / std::sqrt((T) 3),
    (T) 1,                    (T) 1 / std::sqrt((T) 2), (T) 1 / std::sqrt((T) 2),
    (T) 1 / std::sqrt((T) 2), (T) 1 / std::sqrt((T) 3), (T) 1 / std::sqrt((T) 3),
    (T) 1 / std::sqrt((T) 2), (T) 1 / std::sqrt((T) 3), (T) 1 / std::sqrt((T) 3) };

template<typename T, template<typename U> class Descriptor>
NextNeighborPop<T,Descriptor>::NextNeighborPop()
{
    typedef Descriptor<T> D;
    for (plint iNeighbor=0; iNeighbor<NextNeighbor<T>::numNeighbors; ++iNeighbor) {
        int const* c = NextNeighbor<T>::c[iNeighbor];
        ids[iNeighbor] = -1;
        for (plint iPop=0; iPop<D::q; ++iPop) {
            if (D::c[iPop][0]==c[0] && D::c[iPop][1]==c[1] && D::c[iPop][2]==c[2]) {
                ids[iNeighbor] = iPop;
                break;
            }
        }
    }
}

}  // namespace plb

#endif  // NEXT_NEIGHBORS_3D_HH
