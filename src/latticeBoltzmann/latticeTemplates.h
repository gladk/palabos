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
 * Helper functions for the implementation of lattice operations. This file is all
 * about efficiency. The generic template code is specialized for commonly
 * used Lattices, so that a maximum performance can be taken out of each
 * case.
 */
#ifndef LATTICE_TEMPLATES_H
#define LATTICE_TEMPLATES_H

#include "core/globalDefs.h"
#include "core/cell.h"
#include "core/util.h"

namespace plb {

/// Helper functions with full-lattice access
template<typename T, template<typename U> class Descriptor>
struct latticeTemplates {

/// Swap ("bounce-back") values of a cell (2D), and apply streaming step
static void swapAndStream2D(Cell<T,Descriptor> **grid, plint iX, plint iY)
{
    const plint half = Descriptor<T>::q/2;
    for (plint iPop=1; iPop<=half; ++iPop) {
        plint nextX = iX + Descriptor<T>::c[iPop][0];
        plint nextY = iY + Descriptor<T>::c[iPop][1];
        T fTmp                   = grid[iX][iY][iPop];
        grid[iX][iY][iPop]       = grid[iX][iY][iPop+half];
        grid[iX][iY][iPop+half]  = grid[nextX][nextY][iPop];
        grid[nextX][nextY][iPop] = fTmp;
     }
}

/// Swap ("bounce-back") values of a cell (3D), and apply streaming step
static void swapAndStream3D(Cell<T,Descriptor> ***grid,
                            plint iX, plint iY, plint iZ)
{
    const plint half = Descriptor<T>::q/2;
    for (plint iPop=1; iPop<=half; ++iPop) {
        plint nextX = iX + Descriptor<T>::c[iPop][0];
        plint nextY = iY + Descriptor<T>::c[iPop][1];
        plint nextZ = iZ + Descriptor<T>::c[iPop][2];
        T fTmp                          = grid[iX][iY][iZ][iPop];
        grid[iX][iY][iZ][iPop]          = grid[iX][iY][iZ][iPop+half];
        grid[iX][iY][iZ][iPop+half]     = grid[nextX][nextY][nextZ][iPop];
        grid[nextX][nextY][nextZ][iPop] = fTmp;
     }
}


};

}  // namespace plb

#include "latticeBoltzmann/latticeTemplates2D.h"
#include "latticeBoltzmann/latticeTemplates3D.h"

#endif
