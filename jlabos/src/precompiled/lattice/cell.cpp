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
 * Definition of a LB cell -- template instantiation.
 */
#include "core/cell.h"
#include "core/cell.hh"
#include "latticeBoltzmann/nearestNeighborLattices2D.h"
#include "latticeBoltzmann/nearestNeighborLattices2D.hh"
#include "latticeBoltzmann/nearestNeighborLattices3D.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.hh"
#include "latticeBoltzmann/mrtLattices.h"
#include "latticeBoltzmann/mrtLattices.hh"
#include "multiPhysics/shanChenLattices2D.h"
#include "multiPhysics/shanChenLattices3D.h"

namespace plb {
#ifdef COMPILE_2D

    template class Cell<FLOAT_T, descriptors::DESCRIPTOR_2D>;

    void iniCellAtEquilibrium(
            Cell<FLOAT_T, descriptors::DESCRIPTOR_2D>& cell, FLOAT_T density, const FLOAT_T velocity[2] );
    
#if NUMBIT_2D == 9
    template class Cell<FLOAT_T, descriptors::ShanChenD2Q9Descriptor>;
    template class Cell<FLOAT_T, descriptors::MRTD2Q9Descriptor>;
    template class Cell<FLOAT_T, descriptors::ForcedShanChenD2Q9Descriptor>;
#endif

#endif  // COMPILE_2D

#ifdef COMPILE_3D

    template class Cell<FLOAT_T, descriptors::DESCRIPTOR_3D>;

    void iniCellAtEquilibrium (
            Cell<FLOAT_T, descriptors::DESCRIPTOR_3D>& cell, FLOAT_T density, const FLOAT_T velocity[3] );

#if NUMBIT_3D == 19
    template class Cell<FLOAT_T, descriptors::MRTD3Q19Descriptor>;
    template class Cell<FLOAT_T, descriptors::ShanChenD3Q19Descriptor>;
    template class Cell<FLOAT_T, descriptors::ForcedShanChenD3Q19Descriptor>;
#endif

#endif  // COMPILE_3D
    
}
