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

/* Main author: Orestis Malaspinas
 */

/** \file
 * Helper functions for domain initialization -- header file.
 */
#ifndef THERMAL_DATA_ANALYSIS_2D_HH
#define THERMAL_DATA_ANALYSIS_2D_HH

#include "core/blockStatistics.h"
#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/reductiveDataProcessorWrapper2D.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "multiPhysics/thermalDataAnalysis2D.h"
#include <cmath>

namespace plb {
    
// The multiBlockLattice2D case

template< typename T,
          template<typename U1> class FluidDescriptor, 
          template<typename U2> class TemperatureDescriptor
        >
T computeNusseltNumber(BlockLattice2D<T,FluidDescriptor>& fluid,
                       BlockLattice2D<T,TemperatureDescriptor>& temperature,
                       Box2D domain,
                       int direction, T deltaX, T kappa, T deltaTemperature)
{
    BoxSumVelocityTemperatureFunctional2D<T,FluidDescriptor,TemperatureDescriptor> functional(direction);
    applyProcessingFunctional(functional, domain, fluid, temperature);
    return (T)1 + functional.getSumVelocityTemperature() / 
            ((T) domain.nCells() * deltaX * deltaTemperature * kappa);
}

template< typename T,
          template<typename U1> class FluidDescriptor, 
          template<typename U2> class TemperatureDescriptor
        >
T computeNusseltNumber(MultiBlockLattice2D<T,FluidDescriptor>& fluid,
                       MultiBlockLattice2D<T,TemperatureDescriptor>& temperature,
                       Box2D domain,
                       int direction, T deltaX, T kappa, T deltaTemperature)
{
    BoxSumVelocityTemperatureFunctional2D<T,FluidDescriptor,TemperatureDescriptor> functional(direction);
    applyProcessingFunctional(functional, domain, fluid, temperature);
    return (T)1 + functional.getSumVelocityTemperature() / 
            ((T) domain.nCells() * deltaX * deltaTemperature * kappa);
}

template< typename T,
          template<typename U1> class FluidDescriptor, 
          template<typename U2> class TemperatureDescriptor
        >
BoxSumVelocityTemperatureFunctional2D<T,FluidDescriptor,TemperatureDescriptor>::BoxSumVelocityTemperatureFunctional2D(int direction_)
    : direction(direction_),
      sumVelocityTemperatureId(this->getStatistics().subscribeSum())
{ }

template< typename T,
          template<typename U1> class FluidDescriptor, 
          template<typename U2> class TemperatureDescriptor
        >
void BoxSumVelocityTemperatureFunctional2D<T,FluidDescriptor,TemperatureDescriptor>::process (
        Box2D domain,
        BlockLattice2D<T,FluidDescriptor>& fluid,
        BlockLattice2D<T,TemperatureDescriptor>& temperature )
{
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) 
    {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) 
        {
            Array<T,FluidDescriptor<T>::d> velocity;
            fluid.get(iX,iY).computeVelocity(velocity);
            T localTemperature = temperature.get(iX,iY).computeDensity();
            
            T velocityTemperature = localTemperature * velocity[direction];
            statistics.gatherSum(sumVelocityTemperatureId, velocityTemperature);
        }
    }
}

template< typename T,
          template<typename U1> class FluidDescriptor, 
          template<typename U2> class TemperatureDescriptor
        >
BoxSumVelocityTemperatureFunctional2D<T,FluidDescriptor,TemperatureDescriptor>*
        BoxSumVelocityTemperatureFunctional2D<T,FluidDescriptor,TemperatureDescriptor>::clone() const
{
    return new BoxSumVelocityTemperatureFunctional2D<T,FluidDescriptor,TemperatureDescriptor>(*this);
}

template< typename T,
          template<typename U1> class FluidDescriptor, 
          template<typename U2> class TemperatureDescriptor
        >
T BoxSumVelocityTemperatureFunctional2D<T,FluidDescriptor,TemperatureDescriptor>::getSumVelocityTemperature() const 
{
    return this->getStatistics().getSum(sumVelocityTemperatureId);
}

}  // namespace plb

#endif  // THERMAL_DATA_ANALYSIS_2D_HH
