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

#ifndef SWIG_MULTI_BLOCK_INFO_3D_HH
#define SWIG_MULTI_BLOCK_INFO_3D_HH

#include "plbWrapper/block/multiBlockInfo3D.h"

namespace plb {

template<typename T>
MultiNTensorInfo3D<T>::MultiNTensorInfo3D(MultiNTensorField3D<T> const& multiBlock)
    : nx(),
      ny(),
      nz(),
      numBlocks(),
      numAllocatedCells()
{
    getMultiBlockInfo(multiBlock, nx, ny, nz, numBlocks,
                      smallest, largest, numAllocatedCells);
}

template<typename T>
plint MultiNTensorInfo3D<T>::getNx() const {
    return nx;
}

template<typename T>
plint MultiNTensorInfo3D<T>::getNy() const {
    return ny;
}

template<typename T>
plint MultiNTensorInfo3D<T>::getNz() const {
    return nz;
}

template<typename T>
plint MultiNTensorInfo3D<T>::getNumBlocks() const {
    return numBlocks;
}

template<typename T>
plint MultiNTensorInfo3D<T>::getNumAllocatedCells() const {
    return numAllocatedCells;
}

template<typename T>
Box3D MultiNTensorInfo3D<T>::getSmallestBlock() const {
    return smallest;
}

template<typename T>
Box3D MultiNTensorInfo3D<T>::getLargestBlock() const {
    return largest;
}

}  // namespace plb

#endif  // SWIG_MULTI_BLOCK_INFO_3D_HH
