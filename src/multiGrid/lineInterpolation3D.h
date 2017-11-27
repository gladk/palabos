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

/* Main author: Daniel Lagrava
 **/

/** \file
 * Interpolations over a line in coarse and fine coordinates in 3D -- header file.
 */
#ifndef LINE_INTERPOLATION_3D_H
#define LINE_INTERPOLATION_3D_H

#include "core/geometry3D.h"
#include "multiBlock/multiBlockLattice3D.h"


namespace plb {

/* ***************** X coordinate **************** */
template<typename T, template<typename U> class Descriptor>
void cubicInterpolationOverCoarseLineX(Box3D domain,   
         BlockLattice3D<T,Descriptor>& coarseLattice,
         BlockLattice3D<T,Descriptor>& fineLattice,
         RescaleEngine<T,Descriptor>* rescaleEngine);

template<typename T, template<typename U> class Descriptor>
void cubicInterpolationOverFineLineX(Box3D domain,
         BlockLattice3D<T,Descriptor>& coarseLattice,
         BlockLattice3D<T,Descriptor>& fineLattice,
         RescaleEngine<T,Descriptor>* rescaleEngine);
         
template<typename T, template<typename U> class Descriptor>
void cubicInterpolationForEdgeX(Box3D domain,
         BlockLattice3D<T,Descriptor>& coarseLattice,
         BlockLattice3D<T,Descriptor>& fineLattice,
         RescaleEngine<T,Descriptor>* rescaleEngine, plint orientation );
         

/* ***************** Y coordinate **************** */
template<typename T, template<typename U> class Descriptor>
void cubicInterpolationOverCoarseLineY(Box3D domain,
         BlockLattice3D<T,Descriptor>& coarseLattice,
         BlockLattice3D<T,Descriptor>& fineLattice,
         RescaleEngine<T,Descriptor>* rescaleEngine);

template<typename T, template<typename U> class Descriptor>
void cubicInterpolationOverFineLineY(Box3D domain,
         BlockLattice3D<T,Descriptor>& coarseLattice,
         BlockLattice3D<T,Descriptor>& fineLattice,
         RescaleEngine<T,Descriptor>* rescaleEngine);
         
template<typename T, template<typename U> class Descriptor>
void cubicInterpolationForEdgeY(Box3D domain,
         BlockLattice3D<T,Descriptor>& coarseLattice,
         BlockLattice3D<T,Descriptor>& fineLattice,
         RescaleEngine<T,Descriptor>* rescaleEngine, plint orientation );


/* ***************** Z coordinate **************** */
template<typename T, template<typename U> class Descriptor>
void cubicInterpolationOverCoarseLineZ(Box3D domain,
         BlockLattice3D<T,Descriptor>& coarseLattice,
         BlockLattice3D<T,Descriptor>& fineLattice,
         RescaleEngine<T,Descriptor>* rescaleEngine);

template<typename T, template<typename U> class Descriptor>
void cubicInterpolationOverFineLineZ(Box3D domain,
         BlockLattice3D<T,Descriptor>& coarseLattice,
         BlockLattice3D<T,Descriptor>& fineLattice,
         RescaleEngine<T,Descriptor>* rescaleEngine);


} // namespace plb

#endif // LINE_INTERPOLATION_3D_H

