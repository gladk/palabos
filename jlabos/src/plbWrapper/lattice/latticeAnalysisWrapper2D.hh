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
#ifndef SWIG_LATTICE_ANALYSIS_WRAPPER_2D_HH
#define SWIG_LATTICE_ANALYSIS_WRAPPER_2D_HH

#include "plbWrapper/lattice/latticeAnalysisWrapper2D.h"
#include "plbWrapper/lattice/latticeAnalysisFunctional2D.h"
#include "multiBlock/multiBlockGenerator2D.h"
#include "dataProcessors/dataAnalysisWrapper2D.h"
#include "multiBlock/reductiveMultiDataProcessorWrapper2D.h"
#include "multiBlock/multiDataProcessorWrapper2D.h"
#include "dataProcessors/dataAnalysisFunctional2D.h"


namespace plb {

template<typename T, template<typename U> class Descriptor>
void computeDensity( MultiBlockLattice2D<T,Descriptor>& lattice,
                     MultiNTensorField2D<T>& density, Box2D domain )
{
    applyProcessingFunctional (
            new N_BoxDensityFunctional2D<T,Descriptor>, domain, lattice, density );
}

template<typename T, template<typename U> class Descriptor>
void maskedComputeDensity( MultiBlockLattice2D<T,Descriptor>& lattice,
                           MultiNTensorField2D<T>& density, MultiNTensorField2D<int>& mask, Box2D domain )
{
    applyProcessingFunctional (
            new Masked_N_BoxDensityFunctional2D<T,Descriptor>, domain, lattice, density, mask );
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* pypal_computeDensity(
        MultiBlockLattice2D<T,Descriptor>& lattice,
        Box2D domain )
{
    MultiNTensorField2D<T>* density =
        generateMultiNTensorField<T>(lattice, domain, 1);
    computeDensity(lattice, *density, domain);
    return density;
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* maskedPypal_computeDensity(
        MultiBlockLattice2D<T,Descriptor>& lattice, MultiNTensorField2D<int>& mask,
        Box2D domain )
{
    MultiNTensorField2D<T>* density =
        generateMultiNTensorField<T>(lattice, domain, 1);
    maskedComputeDensity(lattice, *density, mask, domain);
    return density;
}

template<typename T, template<typename U> class Descriptor>
void computeKineticEnergy (
        MultiBlockLattice2D<T,Descriptor>& lattice, MultiNTensorField2D<T>& energy, Box2D domain )
{
    applyProcessingFunctional (
            new N_BoxKineticEnergyFunctional2D<T,Descriptor>, domain, lattice, energy );
}

template<typename T, template<typename U> class Descriptor>
void maskedComputeKineticEnergy (
        MultiBlockLattice2D<T,Descriptor>& lattice, MultiNTensorField2D<T>& energy, MultiNTensorField2D<int>& mask, Box2D domain )
{
    applyProcessingFunctional (
            new Masked_N_BoxKineticEnergyFunctional2D<T,Descriptor>, domain, lattice, energy, mask );
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* pypal_computeKineticEnergy (
        MultiBlockLattice2D<T,Descriptor>& lattice, Box2D domain )
{
    MultiNTensorField2D<T>* energy =
        generateMultiNTensorField<T>(lattice, domain, 1);
    computeKineticEnergy(lattice, *energy, domain);
    return energy;
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* maskedPypal_computeKineticEnergy (
        MultiBlockLattice2D<T,Descriptor>& lattice, MultiNTensorField2D<int>& mask, Box2D domain )
{
    MultiNTensorField2D<T>* energy =
        generateMultiNTensorField<T>(lattice, domain, 1);
    maskedComputeKineticEnergy(lattice, *energy, mask, domain);
    return energy;
}

template<typename T, template<typename U> class Descriptor>
void computeVelocityNorm (
        MultiBlockLattice2D<T,Descriptor>& lattice,
        MultiNTensorField2D<T>& velocityNorm, Box2D domain)
{
    applyProcessingFunctional (
            new N_BoxVelocityNormFunctional2D<T,Descriptor>, domain, lattice, velocityNorm );
}

template<typename T, template<typename U> class Descriptor>
void maskedComputeVelocityNorm (
        MultiBlockLattice2D<T,Descriptor>& lattice,
        MultiNTensorField2D<T>& velocityNorm, MultiNTensorField2D<int>& mask, Box2D domain)
{
    applyProcessingFunctional (
            new Masked_N_BoxVelocityNormFunctional2D<T,Descriptor>, domain, lattice, velocityNorm, mask );
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* pypal_computeVelocityNorm (
        MultiBlockLattice2D<T,Descriptor>& lattice,
        Box2D domain)
{
    MultiNTensorField2D<T>* velNorm =
        generateMultiNTensorField<T>(lattice, domain, 1);
    computeVelocityNorm(lattice, *velNorm, domain);
    return velNorm;
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* maskedPypal_computeVelocityNorm (
        MultiBlockLattice2D<T,Descriptor>& lattice, MultiNTensorField2D<int>& mask,
        Box2D domain)
{
    MultiNTensorField2D<T>* velNorm =
        generateMultiNTensorField<T>(lattice, domain, 1);
    maskedComputeVelocityNorm(lattice, *velNorm, mask, domain);
    return velNorm;
}

template<typename T, template<typename U> class Descriptor>
void computeVelocityComponent (
        MultiBlockLattice2D<T,Descriptor>& lattice,
        MultiNTensorField2D<T>& velocityComponent,
        Box2D domain, plint iComponent )
{
    applyProcessingFunctional (
            new N_BoxVelocityComponentFunctional2D<T,Descriptor>(iComponent),
            domain, lattice, velocityComponent );
}

template<typename T, template<typename U> class Descriptor>
void maskedComputeVelocityComponent (
        MultiBlockLattice2D<T,Descriptor>& lattice,
        MultiNTensorField2D<T>& velocityComponent, MultiNTensorField2D<int>& mask,
        Box2D domain, plint iComponent )
{
    applyProcessingFunctional (
            new Masked_N_BoxVelocityComponentFunctional2D<T,Descriptor>(iComponent),
            domain, lattice, velocityComponent, mask );
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* pypal_computeVelocityComponent (
        MultiBlockLattice2D<T,Descriptor>& lattice,
        Box2D domain, plint iComponent )
{
    MultiNTensorField2D<T>* velComponent =
        generateMultiNTensorField<T>(lattice, domain, 1);
    computeVelocityComponent(lattice, *velComponent, domain, iComponent);
    return velComponent;
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* maskedPypal_computeVelocityComponent (
        MultiBlockLattice2D<T,Descriptor>& lattice, MultiNTensorField2D<int>& mask,
        Box2D domain, plint iComponent )
{
    MultiNTensorField2D<T>* velComponent =
        generateMultiNTensorField<T>(lattice, domain, 1);
    maskedComputeVelocityComponent(lattice, *velComponent, mask, domain, iComponent);
    return velComponent;
}

template<typename T, template<typename U> class Descriptor>
void computeVelocity(MultiBlockLattice2D<T,Descriptor>& lattice,
                     MultiNTensorField2D<T>& velocity, Box2D domain)
{
    applyProcessingFunctional (
            new N_BoxVelocityFunctional2D<T,Descriptor>, domain, lattice, velocity );
}

template<typename T, template<typename U> class Descriptor>
void maskedComputeVelocity(MultiBlockLattice2D<T,Descriptor>& lattice,
                     MultiNTensorField2D<T>& velocity, MultiNTensorField2D<int>& mask, Box2D domain)
{
    applyProcessingFunctional (
            new Masked_N_BoxVelocityFunctional2D<T,Descriptor>, domain, lattice, velocity, mask );
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* pypal_computeVelocity (
        MultiBlockLattice2D<T,Descriptor>& lattice,
        Box2D domain)
{
    MultiNTensorField2D<T>* velocity =
        generateMultiNTensorField<T>(lattice, domain, 2);
    computeVelocity(lattice, *velocity, domain);
    return velocity;
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* maskedPypal_computeVelocity (
        MultiBlockLattice2D<T,Descriptor>& lattice, MultiNTensorField2D<int>& mask,
        Box2D domain)
{
    MultiNTensorField2D<T>* velocity =
        generateMultiNTensorField<T>(lattice, domain, 2);
    maskedComputeVelocity(lattice, *velocity, mask, domain);
    return velocity;
}

template<typename T, template<typename U> class Descriptor>
void computePiNeq( MultiBlockLattice2D<T,Descriptor>& lattice,
                              MultiNTensorField2D<T>& PiNeq,
                              Box2D domain )
{
    applyProcessingFunctional (
            new N_BoxPiNeqFunctional2D<T,Descriptor>, domain, lattice, PiNeq );
}

template<typename T, template<typename U> class Descriptor>
void maskedComputePiNeq( MultiBlockLattice2D<T,Descriptor>& lattice,
                              MultiNTensorField2D<T>& PiNeq, MultiNTensorField2D<int>& mask,
                              Box2D domain )
{
    applyProcessingFunctional (
            new Masked_N_BoxPiNeqFunctional2D<T,Descriptor>, domain, lattice, PiNeq, mask );
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* pypal_computePiNeq (
        MultiBlockLattice2D<T,Descriptor>& lattice,
        Box2D domain )
{
    MultiNTensorField2D<T>* piNeq =
        generateMultiNTensorField<T>(lattice, domain, 3);
    computePiNeq(lattice, *piNeq, domain);
    return piNeq;
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* maskedPypal_computePiNeq (
        MultiBlockLattice2D<T,Descriptor>& lattice, MultiNTensorField2D<int>& mask,
        Box2D domain )
{
    MultiNTensorField2D<T>* piNeq =
        generateMultiNTensorField<T>(lattice, domain, 3);
    maskedComputePiNeq(lattice, *piNeq, mask, domain);
    return piNeq;
}

template<typename T, template<typename U> class Descriptor>
void computeShearStress( MultiBlockLattice2D<T,Descriptor>& lattice,
                              MultiNTensorField2D<T>& ShearStress,
                              Box2D domain )
{
    applyProcessingFunctional (
            new N_BoxShearStressFunctional2D<T,Descriptor>, domain, lattice, ShearStress );
}

template<typename T, template<typename U> class Descriptor>
void maskedComputeShearStress( MultiBlockLattice2D<T,Descriptor>& lattice,
                              MultiNTensorField2D<T>& ShearStress, MultiNTensorField2D<int>& mask,
                              Box2D domain )
{
    applyProcessingFunctional (
            new Masked_N_BoxShearStressFunctional2D<T,Descriptor>, domain, lattice, ShearStress, mask );
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* pypal_computeShearStress (
        MultiBlockLattice2D<T,Descriptor>& lattice,
        Box2D domain )
{
    MultiNTensorField2D<T>* stress =
        generateMultiNTensorField<T>(lattice, domain, 3);
    computeShearStress(lattice, *stress, domain);
    return stress;
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* maskedPypal_computeShearStress (
        MultiBlockLattice2D<T,Descriptor>& lattice, MultiNTensorField2D<int>& mask,
        Box2D domain )
{
    MultiNTensorField2D<T>* stress =
        generateMultiNTensorField<T>(lattice, domain, 3);
    maskedComputeShearStress(lattice, *stress, mask, domain);
    return stress;
}

template<typename T, template<typename U> class Descriptor>
void computeStrainRateFromStress( MultiBlockLattice2D<T,Descriptor>& lattice,
                                  MultiNTensorField2D<T>& S,
                                  Box2D domain )
{
    applyProcessingFunctional (
            new N_BoxStrainRateFromStressFunctional2D<T,Descriptor>, domain, lattice, S );
}

template<typename T, template<typename U> class Descriptor>
void maskedComputeStrainRateFromStress( MultiBlockLattice2D<T,Descriptor>& lattice,
                                        MultiNTensorField2D<T>& S, MultiNTensorField2D<int>& mask,
                                        Box2D domain )
{
    applyProcessingFunctional (
            new Masked_N_BoxStrainRateFromStressFunctional2D<T,Descriptor>, domain, lattice, S, mask );
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* pypal_computeStrainRateFromStress (
        MultiBlockLattice2D<T,Descriptor>& lattice,
        Box2D domain )
{
    MultiNTensorField2D<T>* strainRate =
        generateMultiNTensorField<T>(lattice, domain, 3);
    computeStrainRateFromStress(lattice, *strainRate, domain);
    return strainRate;
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* maskedPypal_computeStrainRateFromStress (
        MultiBlockLattice2D<T,Descriptor>& lattice, MultiNTensorField2D<int>& mask,
        Box2D domain )
{
    MultiNTensorField2D<T>* strainRate =
        generateMultiNTensorField<T>(lattice, domain, 3);
    maskedComputeStrainRateFromStress(lattice, *strainRate, mask, domain);
    return strainRate;
}

template<typename T, template<typename U> class Descriptor>
void computePopulation( MultiBlockLattice2D<T,Descriptor>& lattice,
                        MultiNTensorField2D<T>& population,
                        Box2D domain, plint iPop )
{
    applyProcessingFunctional (
            new N_BoxPopulationFunctional2D<T,Descriptor>(iPop), domain, lattice, population );
}

template<typename T, template<typename U> class Descriptor>
void maskedComputePopulation( MultiBlockLattice2D<T,Descriptor>& lattice,
                        MultiNTensorField2D<T>& population, MultiNTensorField2D<int>& mask,
                        Box2D domain, plint iPop )
{
    applyProcessingFunctional (
            new Masked_N_BoxPopulationFunctional2D<T,Descriptor>(iPop), domain, lattice, population, mask );
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* pypal_computePopulation (
        MultiBlockLattice2D<T,Descriptor>& lattice,
        Box2D domain, plint iPop )
{
    MultiNTensorField2D<T>* population =
        generateMultiNTensorField<T>(lattice, domain, 1);
    computePopulation(lattice, *population, domain, iPop);
    return population;
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* maskedPypal_computePopulation (
        MultiBlockLattice2D<T,Descriptor>& lattice, MultiNTensorField2D<int>& mask,
        Box2D domain, plint iPop )
{
    MultiNTensorField2D<T>* population =
        generateMultiNTensorField<T>(lattice, domain, 1);
    maskedComputePopulation(lattice, *population, mask, domain, iPop);
    return population;
}

template<typename T, template<typename U> class Descriptor>
void computePopulations( MultiBlockLattice2D<T,Descriptor>& lattice,
                        MultiNTensorField2D<T>& populations,
                        Box2D domain )
{
    applyProcessingFunctional (
            new N_BoxPopulationsFunctional2D<T,Descriptor>, domain, lattice, populations );
}

template<typename T, template<typename U> class Descriptor>
void maskedComputePopulations( MultiBlockLattice2D<T,Descriptor>& lattice,
                        MultiNTensorField2D<T>& populations, MultiNTensorField2D<int>& mask,
                        Box2D domain )
{
    applyProcessingFunctional (
            new Masked_N_BoxPopulationsFunctional2D<T,Descriptor>, domain, lattice, populations, mask );
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* pypal_computePopulations (
        MultiBlockLattice2D<T,Descriptor>& lattice,
        Box2D domain )
{
    MultiNTensorField2D<T>* populations =
        generateMultiNTensorField<T>(lattice, domain, Descriptor<T>::q);
    computePopulations(lattice, *populations, domain);
    return populations;
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField2D<T>* maskedPypal_computePopulations (
        MultiBlockLattice2D<T,Descriptor>& lattice, MultiNTensorField2D<int>& mask,
        Box2D domain )
{
    MultiNTensorField2D<T>* populations =
        generateMultiNTensorField<T>(lattice, domain, Descriptor<T>::q);
    maskedComputePopulations(lattice, *populations, mask, domain);
    return populations;
}

}  // namespace plb

#endif  // SWIG_LATTICE_ANALYSIS_WRAPPER_2D_HH
