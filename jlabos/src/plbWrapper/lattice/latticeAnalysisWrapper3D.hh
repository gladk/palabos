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
#ifndef SWIG_LATTICE_ANALYSIS_WRAPPER_3D_HH
#define SWIG_LATTICE_ANALYSIS_WRAPPER_3D_HH

#include "plbWrapper/lattice/latticeAnalysisWrapper3D.h"
#include "plbWrapper/lattice/latticeAnalysisFunctional3D.h"
#include "multiBlock/multiBlockGenerator3D.h"
#include "dataProcessors/dataAnalysisWrapper3D.h"
#include "multiBlock/reductiveMultiDataProcessorWrapper3D.h"
#include "multiBlock/multiDataProcessorWrapper3D.h"
#include "dataProcessors/dataAnalysisFunctional3D.h"


namespace plb {

template<typename T, template<typename U> class Descriptor>
void computeDensity( MultiBlockLattice3D<T,Descriptor>& lattice,
                     MultiNTensorField3D<T>& density, Box3D domain )
{
    applyProcessingFunctional (
            new N_BoxDensityFunctional3D<T,Descriptor>, domain, lattice, density );
}

template<typename T, template<typename U> class Descriptor>
void maskedComputeDensity( MultiBlockLattice3D<T,Descriptor>& lattice,
                           MultiNTensorField3D<T>& density, MultiNTensorField3D<int>& mask, Box3D domain )
{
    applyProcessingFunctional (
            new Masked_N_BoxDensityFunctional3D<T,Descriptor>, domain, lattice, density, mask );
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* pypal_computeDensity(
        MultiBlockLattice3D<T,Descriptor>& lattice,
        Box3D domain )
{
    MultiNTensorField3D<T>* density =
        generateMultiNTensorField<T>(lattice, domain, 1);
    computeDensity(lattice, *density, domain);
    return density;
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* maskedPypal_computeDensity(
        MultiBlockLattice3D<T,Descriptor>& lattice, MultiNTensorField3D<int>& mask,
        Box3D domain )
{
    MultiNTensorField3D<T>* density =
        generateMultiNTensorField<T>(lattice, domain, 1);
    maskedComputeDensity(lattice, *density, mask, domain);
    return density;
}

template<typename T, template<typename U> class Descriptor>
void computeKineticEnergy (
        MultiBlockLattice3D<T,Descriptor>& lattice, MultiNTensorField3D<T>& energy, Box3D domain )
{
    applyProcessingFunctional (
            new N_BoxKineticEnergyFunctional3D<T,Descriptor>, domain, lattice, energy );
}

template<typename T, template<typename U> class Descriptor>
void maskedComputeKineticEnergy (
        MultiBlockLattice3D<T,Descriptor>& lattice, MultiNTensorField3D<T>& energy, MultiNTensorField3D<int>& mask, Box3D domain )
{
    applyProcessingFunctional (
            new Masked_N_BoxKineticEnergyFunctional3D<T,Descriptor>, domain, lattice, energy, mask );
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* pypal_computeKineticEnergy (
        MultiBlockLattice3D<T,Descriptor>& lattice, Box3D domain )
{
    MultiNTensorField3D<T>* energy =
        generateMultiNTensorField<T>(lattice, domain, 1);
    computeKineticEnergy(lattice, *energy, domain);
    return energy;
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* maskedPypal_computeKineticEnergy (
        MultiBlockLattice3D<T,Descriptor>& lattice, MultiNTensorField3D<int>& mask, Box3D domain )
{
    MultiNTensorField3D<T>* energy =
        generateMultiNTensorField<T>(lattice, domain, 1);
    maskedComputeKineticEnergy(lattice, *energy, mask, domain);
    return energy;
}

template<typename T, template<typename U> class Descriptor>
void computeVelocityNorm (
        MultiBlockLattice3D<T,Descriptor>& lattice,
        MultiNTensorField3D<T>& velocityNorm, Box3D domain)
{
    applyProcessingFunctional (
            new N_BoxVelocityNormFunctional3D<T,Descriptor>, domain, lattice, velocityNorm );
}

template<typename T, template<typename U> class Descriptor>
void maskedComputeVelocityNorm (
        MultiBlockLattice3D<T,Descriptor>& lattice,
        MultiNTensorField3D<T>& velocityNorm, MultiNTensorField3D<int>& mask, Box3D domain)
{
    applyProcessingFunctional (
            new Masked_N_BoxVelocityNormFunctional3D<T,Descriptor>, domain, lattice, velocityNorm, mask );
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* pypal_computeVelocityNorm (
        MultiBlockLattice3D<T,Descriptor>& lattice,
        Box3D domain)
{
    MultiNTensorField3D<T>* velNorm =
        generateMultiNTensorField<T>(lattice, domain, 1);
    computeVelocityNorm(lattice, *velNorm, domain);
    return velNorm;
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* maskedPypal_computeVelocityNorm (
        MultiBlockLattice3D<T,Descriptor>& lattice, MultiNTensorField3D<int>& mask,
        Box3D domain)
{
    MultiNTensorField3D<T>* velNorm =
        generateMultiNTensorField<T>(lattice, domain, 1);
    maskedComputeVelocityNorm(lattice, *velNorm, mask, domain);
    return velNorm;
}

template<typename T, template<typename U> class Descriptor>
void computeVelocityComponent (
        MultiBlockLattice3D<T,Descriptor>& lattice,
        MultiNTensorField3D<T>& velocityComponent,
        Box3D domain, plint iComponent )
{
    applyProcessingFunctional (
            new N_BoxVelocityComponentFunctional3D<T,Descriptor>(iComponent),
            domain, lattice, velocityComponent );
}

template<typename T, template<typename U> class Descriptor>
void maskedComputeVelocityComponent (
        MultiBlockLattice3D<T,Descriptor>& lattice,
        MultiNTensorField3D<T>& velocityComponent, MultiNTensorField3D<int>& mask,
        Box3D domain, plint iComponent )
{
    applyProcessingFunctional (
            new Masked_N_BoxVelocityComponentFunctional3D<T,Descriptor>(iComponent),
            domain, lattice, velocityComponent, mask );
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* pypal_computeVelocityComponent (
        MultiBlockLattice3D<T,Descriptor>& lattice,
        Box3D domain, plint iComponent )
{
    MultiNTensorField3D<T>* velComponent =
        generateMultiNTensorField<T>(lattice, domain, 1);
    computeVelocityComponent(lattice, *velComponent, domain, iComponent);
    return velComponent;
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* maskedPypal_computeVelocityComponent (
        MultiBlockLattice3D<T,Descriptor>& lattice, MultiNTensorField3D<int>& mask,
        Box3D domain, plint iComponent )
{
    MultiNTensorField3D<T>* velComponent =
        generateMultiNTensorField<T>(lattice, domain, 1);
    maskedComputeVelocityComponent(lattice, *velComponent, mask, domain, iComponent);
    return velComponent;
}

template<typename T, template<typename U> class Descriptor>
void computeVelocity(MultiBlockLattice3D<T,Descriptor>& lattice,
                     MultiNTensorField3D<T>& velocity, Box3D domain)
{
    applyProcessingFunctional (
            new N_BoxVelocityFunctional3D<T,Descriptor>, domain, lattice, velocity );
}

template<typename T, template<typename U> class Descriptor>
void maskedComputeVelocity(MultiBlockLattice3D<T,Descriptor>& lattice,
                     MultiNTensorField3D<T>& velocity, MultiNTensorField3D<int>& mask, Box3D domain)
{
    applyProcessingFunctional (
            new Masked_N_BoxVelocityFunctional3D<T,Descriptor>, domain, lattice, velocity, mask );
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* pypal_computeVelocity (
        MultiBlockLattice3D<T,Descriptor>& lattice,
        Box3D domain)
{
    MultiNTensorField3D<T>* velocity =
        generateMultiNTensorField<T>(lattice, domain, 3);
    computeVelocity(lattice, *velocity, domain);
    return velocity;
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* maskedPypal_computeVelocity (
        MultiBlockLattice3D<T,Descriptor>& lattice, MultiNTensorField3D<int>& mask,
        Box3D domain)
{
    MultiNTensorField3D<T>* velocity =
        generateMultiNTensorField<T>(lattice, domain, 3);
    maskedComputeVelocity(lattice, *velocity, mask, domain);
    return velocity;
}

template<typename T, template<typename U> class Descriptor>
void computePiNeq( MultiBlockLattice3D<T,Descriptor>& lattice,
                              MultiNTensorField3D<T>& PiNeq,
                              Box3D domain )
{
    applyProcessingFunctional (
            new N_BoxPiNeqFunctional3D<T,Descriptor>, domain, lattice, PiNeq );
}

template<typename T, template<typename U> class Descriptor>
void maskedComputePiNeq( MultiBlockLattice3D<T,Descriptor>& lattice,
                              MultiNTensorField3D<T>& PiNeq, MultiNTensorField3D<int>& mask,
                              Box3D domain )
{
    applyProcessingFunctional (
            new Masked_N_BoxPiNeqFunctional3D<T,Descriptor>, domain, lattice, PiNeq, mask );
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* pypal_computePiNeq (
        MultiBlockLattice3D<T,Descriptor>& lattice,
        Box3D domain )
{
    MultiNTensorField3D<T>* piNeq =
        generateMultiNTensorField<T>(lattice, domain, 6);
    computePiNeq(lattice, *piNeq, domain);
    return piNeq;
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* maskedPypal_computePiNeq (
        MultiBlockLattice3D<T,Descriptor>& lattice, MultiNTensorField3D<int>& mask,
        Box3D domain )
{
    MultiNTensorField3D<T>* piNeq =
        generateMultiNTensorField<T>(lattice, domain, 6);
    maskedComputePiNeq(lattice, *piNeq, mask, domain);
    return piNeq;
}

template<typename T, template<typename U> class Descriptor>
void computeShearStress( MultiBlockLattice3D<T,Descriptor>& lattice,
                              MultiNTensorField3D<T>& ShearStress,
                              Box3D domain )
{
    applyProcessingFunctional (
            new N_BoxShearStressFunctional3D<T,Descriptor>, domain, lattice, ShearStress );
}

template<typename T, template<typename U> class Descriptor>
void maskedComputeShearStress( MultiBlockLattice3D<T,Descriptor>& lattice,
                              MultiNTensorField3D<T>& ShearStress, MultiNTensorField3D<int>& mask,
                              Box3D domain )
{
    applyProcessingFunctional (
            new Masked_N_BoxShearStressFunctional3D<T,Descriptor>, domain, lattice, ShearStress, mask );
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* pypal_computeShearStress (
        MultiBlockLattice3D<T,Descriptor>& lattice,
        Box3D domain )
{
    MultiNTensorField3D<T>* piNeq =
        generateMultiNTensorField<T>(lattice, domain, 6);
    computeShearStress(lattice, *piNeq, domain);
    return piNeq;
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* maskedPypal_computeShearStress (
        MultiBlockLattice3D<T,Descriptor>& lattice, MultiNTensorField3D<int>& mask,
        Box3D domain )
{
    MultiNTensorField3D<T>* piNeq =
        generateMultiNTensorField<T>(lattice, domain, 6);
    maskedComputeShearStress(lattice, *piNeq, mask, domain);
    return piNeq;
}

template<typename T, template<typename U> class Descriptor>
void computeStrainRateFromStress( MultiBlockLattice3D<T,Descriptor>& lattice,
                                  MultiNTensorField3D<T>& S,
                                  Box3D domain )
{
    applyProcessingFunctional (
            new N_BoxStrainRateFromStressFunctional3D<T,Descriptor>, domain, lattice, S );
}

template<typename T, template<typename U> class Descriptor>
void maskedComputeStrainRateFromStress( MultiBlockLattice3D<T,Descriptor>& lattice,
                                        MultiNTensorField3D<T>& S, MultiNTensorField3D<int>& mask,
                                        Box3D domain )
{
    applyProcessingFunctional (
            new Masked_N_BoxStrainRateFromStressFunctional3D<T,Descriptor>, domain, lattice, S, mask );
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* pypal_computeStrainRateFromStress (
        MultiBlockLattice3D<T,Descriptor>& lattice,
        Box3D domain )
{
    MultiNTensorField3D<T>* strainRate =
        generateMultiNTensorField<T>(lattice, domain, 6);
    computeStrainRateFromStress(lattice, *strainRate, domain);
    return strainRate;
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* maskedPypal_computeStrainRateFromStress (
        MultiBlockLattice3D<T,Descriptor>& lattice, MultiNTensorField3D<int>& mask,
        Box3D domain )
{
    MultiNTensorField3D<T>* strainRate =
        generateMultiNTensorField<T>(lattice, domain, 6);
    maskedComputeStrainRateFromStress(lattice, *strainRate, mask, domain);
    return strainRate;
}

template<typename T, template<typename U> class Descriptor>
void computePopulation( MultiBlockLattice3D<T,Descriptor>& lattice,
                        MultiNTensorField3D<T>& population,
                        Box3D domain, plint iPop )
{
    applyProcessingFunctional (
            new N_BoxPopulationFunctional3D<T,Descriptor>(iPop), domain, lattice, population );
}

template<typename T, template<typename U> class Descriptor>
void maskedComputePopulation( MultiBlockLattice3D<T,Descriptor>& lattice,
                        MultiNTensorField3D<T>& population, MultiNTensorField3D<int>& mask,
                        Box3D domain, plint iPop )
{
    applyProcessingFunctional (
            new Masked_N_BoxPopulationFunctional3D<T,Descriptor>(iPop), domain, lattice, population, mask );
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* pypal_computePopulation (
        MultiBlockLattice3D<T,Descriptor>& lattice,
        Box3D domain, plint iPop )
{
    MultiNTensorField3D<T>* population =
        generateMultiNTensorField<T>(lattice, domain, 1);
    computePopulation(lattice, *population, domain, iPop);
    return population;
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* maskedPypal_computePopulation (
        MultiBlockLattice3D<T,Descriptor>& lattice, MultiNTensorField3D<int>& mask,
        Box3D domain, plint iPop )
{
    MultiNTensorField3D<T>* population =
        generateMultiNTensorField<T>(lattice, domain, 1);
    maskedComputePopulation(lattice, *population, mask, domain, iPop);
    return population;
}

template<typename T, template<typename U> class Descriptor>
void computePopulations( MultiBlockLattice3D<T,Descriptor>& lattice,
                        MultiNTensorField3D<T>& populations,
                        Box3D domain )
{
    applyProcessingFunctional (
            new N_BoxPopulationsFunctional3D<T,Descriptor>, domain, lattice, populations );
}

template<typename T, template<typename U> class Descriptor>
void maskedComputePopulations( MultiBlockLattice3D<T,Descriptor>& lattice,
                        MultiNTensorField3D<T>& populations, MultiNTensorField3D<int>& mask,
                        Box3D domain )
{
    applyProcessingFunctional (
            new Masked_N_BoxPopulationsFunctional3D<T,Descriptor>, domain, lattice, populations, mask );
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* pypal_computePopulations (
        MultiBlockLattice3D<T,Descriptor>& lattice,
        Box3D domain )
{
    MultiNTensorField3D<T>* populations =
        generateMultiNTensorField<T>(lattice, domain, Descriptor<T>::q);
    computePopulations(lattice, *populations, domain);
    return populations;
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* maskedPypal_computePopulations (
        MultiBlockLattice3D<T,Descriptor>& lattice, MultiNTensorField3D<int>& mask,
        Box3D domain )
{
    MultiNTensorField3D<T>* populations =
        generateMultiNTensorField<T>(lattice, domain, Descriptor<T>::q);
    maskedComputePopulations(lattice, *populations, mask, domain);
    return populations;
}


template<typename T, template<typename U> class Descriptor>
void compute_UPO_Rhs(MultiNTensorField3D<T>& lattice,
                     MultiNTensorField3D<T>& result, Box3D domain, T omega)
{
    applyProcessingFunctional (
            new UPO_Rhs_Functional3D<T,Descriptor>(omega), domain, lattice, result );
}

template<typename T, template<typename U> class Descriptor>
void masked_compute_UPO_Rhs(MultiNTensorField3D<T>& lattice,
                            MultiNTensorField3D<T>& result,
                            MultiNTensorField3D<int>& mask, Box3D domain, T omega)
{
    applyProcessingFunctional (
            new Masked_UPO_Rhs_Functional3D<T,Descriptor>(omega), domain, lattice, result, mask );
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* compute_UPO_Rhs (
                           MultiNTensorField3D<T>& lattice, Box3D domain, T omega )
{
    MultiNTensorField3D<T>* result =
        generateMultiNTensorField<T>(lattice, domain, Descriptor<T>::q);
    compute_UPO_Rhs<T,Descriptor>(lattice, *result, domain, omega);
    return result;
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* masked_compute_UPO_Rhs (
                                 MultiNTensorField3D<T>& lattice,
                                 MultiNTensorField3D<int>& mask, Box3D domain, T omega )
{
    MultiNTensorField3D<T>* result =
        generateMultiNTensorField<T>(lattice, domain, Descriptor<T>::q);
    masked_compute_UPO_Rhs<T,Descriptor>(lattice, *result, mask, domain, omega);
    return result;
}



template<typename T, template<typename U> class Descriptor>
void UPO_ApplyJ( MultiNTensorField3D<T>& f,
                 MultiNTensorField3D<T>& g,
                 MultiNTensorField3D<T>& h, Box3D domain, T omega)
{
    std::vector<MultiNTensorField3D<T>*> fields;
    fields.push_back(&f);
    fields.push_back(&g);
    fields.push_back(&h);
    applyProcessingFunctional (
            new UPO_ApplyJ_Functional3D<T,Descriptor>(omega), domain, fields );
}

template<typename T, template<typename U> class Descriptor>
void masked_UPO_ApplyJ( MultiNTensorField3D<T>& f,
                        MultiNTensorField3D<T>& g,
                        MultiNTensorField3D<T>& h,
                        MultiNTensorField3D<int>& mask, Box3D domain, T omega)
{
    std::vector<MultiNTensorField3D<T>*> fields;
    fields.push_back(&f);
    fields.push_back(&g);
    fields.push_back(&h);
    applyProcessingFunctional (
            new Masked_UPO_ApplyJ_Functional3D<T,Descriptor>(omega), domain, fields, mask );
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* UPO_ApplyJ (
                           MultiNTensorField3D<T>& f,
                           MultiNTensorField3D<T>& g,
                           Box3D domain, T omega )
{
    MultiNTensorField3D<T>* result
        = generateIntersectMultiNTensorField<T>(f,g, domain, Descriptor<T>::q);
    UPO_ApplyJ<T,Descriptor>(f, g, *result, domain, omega);
    return result;
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* masked_UPO_ApplyJ (
                             MultiNTensorField3D<T>& f,
                             MultiNTensorField3D<T>& g,
                             MultiNTensorField3D<int>& mask, Box3D domain, T omega )
{
    MultiNTensorField3D<T>* result
        = generateIntersectMultiNTensorField<T>(f,g, domain, Descriptor<T>::q);
    masked_UPO_ApplyJ<T,Descriptor>(f, g, *result, mask, domain, omega);
    return result;
}



template<typename T, template<typename U> class Descriptor>
void compute_UPO_EnergyDerivative (
                 MultiNTensorField3D<T>& lattice,
                 MultiNTensorField3D<T>& result, Box3D domain )
{
    applyProcessingFunctional (
            new UPO_EnergyDerivative_Functional3D<T,Descriptor>, domain, lattice, result );
}


template<typename T, template<typename U> class Descriptor>
void masked_compute_UPO_EnergyDerivative (
                        MultiNTensorField3D<T>& lattice,
                        MultiNTensorField3D<T>& result,
                        MultiNTensorField3D<int>& mask, Box3D domain)
{
    applyProcessingFunctional (
            new Masked_UPO_EnergyDerivative_Functional3D<T,Descriptor>, domain,
            lattice, result, mask );
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* compute_UPO_EnergyDerivative (
                           MultiNTensorField3D<T>& lattice, Box3D domain )
{
    MultiNTensorField3D<T>* result =
        generateMultiNTensorField<T>(lattice, domain, Descriptor<T>::q);
    compute_UPO_EnergyDerivative<T,Descriptor>(lattice, *result, domain);
    return result;
}

template<typename T, template<typename U> class Descriptor>
MultiNTensorField3D<T>* masked_compute_UPO_EnergyDerivative (
                                 MultiNTensorField3D<T>& lattice,
                                 MultiNTensorField3D<int>& mask, Box3D domain )
{
    MultiNTensorField3D<T>* result =
        generateMultiNTensorField<T>(lattice, domain, Descriptor<T>::q);
    masked_compute_UPO_EnergyDerivative<T,Descriptor>(lattice, *result, mask, domain);
    return result;
}

}  // namespace plb

#endif  // SWIG_LATTICE_ANALYSIS_WRAPPER_3D_HH
