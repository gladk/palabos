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
 * Helper functions for domain initialization -- generic implementation.
 */
#ifndef SWIG_LATTICE_INITIALIZER_WRAPPER_2D_HH
#define SWIG_LATTICE_INITIALIZER_WRAPPER_2D_HH

#include "plbWrapper/lattice/latticeInitializerWrapper2D.h"
#include "dataProcessors/dataInitializerWrapper2D.h"
#include "dataProcessors/dataInitializerFunctional2D.h"
#include "core/cell.h"
#include "atomicBlock/dataProcessorWrapper2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataField2D.h"
#include "multiBlock/multiBlockLattice2D.h"
#include "multiBlock/multiDataField2D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"

namespace plb {

template<typename T, template<typename U> class Descriptor>
void pypalDefineDynamics( MultiBlockLattice2D<T,Descriptor>& lattice,
                          Box2D domain, Dynamics<T,Descriptor>* dynamics )
{
    applyProcessingFunctional (
            new InstantiateDynamicsFunctional2D<T,Descriptor> (
                dynamics->clone() ),
            domain, lattice );
}

template<typename T, template<typename U> class Descriptor>
void maskedDefineDynamics( MultiBlockLattice2D<T,Descriptor>& lattice,
                           MultiNTensorField2D<int>& mask,
                           Box2D domain, Dynamics<T,Descriptor>* dynamics )
{
    applyProcessingFunctional (
            new MaskedIniDynamicsFunctional2D<T,Descriptor> (
                dynamics->clone() ),
            domain, lattice, mask );
}

template<typename T, template<class U> class Descriptor>
void setBoundaryVelocity (
        MultiBlockLattice2D<T,Descriptor>& lattice,
        T* velocity, int numDimIs2, Box2D domain )
{
    Array<T,2> velocityVect(velocity[0], velocity[1]);
    applyProcessingFunctional (
            new SetConstBoundaryVelocityFunctional2D<T,Descriptor>(velocityVect), domain, lattice );
}

template<typename T, template<class U> class Descriptor>
void setBoundaryVelocity (
        MultiBlockLattice2D<T,Descriptor>& lattice,
        MultiNTensorField2D<T>& velocity,
        Box2D domain )
{
    applyProcessingFunctional (
            new N_IniBoundaryVelocityFunctional2D<T,Descriptor>, domain, lattice, velocity );
}

template<typename T, template<class U> class Descriptor>
void maskedSetBoundaryVelocity (
        MultiBlockLattice2D<T,Descriptor>& lattice,
        MultiNTensorField2D<T>& velocity,
        MultiNTensorField2D<int>& mask,
        Box2D domain )
{
    applyProcessingFunctional (
            new Masked_N_IniBoundaryVelocityFunctional2D<T,Descriptor>, domain, lattice, velocity, mask );
}

template<typename T, template<class U> class Descriptor>
void initializeAtEquilibrium (
        MultiBlockLattice2D<T,Descriptor>& lattice, T rho,
        T* velocity, int numDimIs2, Box2D domain )
{
    Array<T,2> velocityVect(velocity[0], velocity[1]);
    applyProcessingFunctional (
            new IniConstEquilibriumFunctional2D<T,Descriptor>(rho, velocityVect, (T)1), domain, lattice );
}

template<typename T, template<class U> class Descriptor>
void initializeAtEquilibrium( MultiBlockLattice2D<T,Descriptor>& lattice,
                              MultiNTensorField2D<T>& density,
                              MultiNTensorField2D<T>& velocity,
                              Box2D domain )
{
    std::vector<MultiBlock2D*> fields;
    fields.push_back(&lattice);
    fields.push_back(&density);
    fields.push_back(&velocity);
    applyProcessingFunctional (
            new N_IniEquilibriumFunctional2D<T,Descriptor>, domain, fields );
}

template<typename T, template<class U> class Descriptor>
void maskedInitializeAtEquilibrium( MultiBlockLattice2D<T,Descriptor>& lattice,
                                    MultiNTensorField2D<T>& density,
                                    MultiNTensorField2D<T>& velocity,
                                    MultiNTensorField2D<int>& mask,
                                    Box2D domain )
{
    std::vector<MultiBlock2D*> fields;
    fields.push_back(&lattice);
    fields.push_back(&density);
    fields.push_back(&velocity);
    fields.push_back(&mask);
    applyProcessingFunctional (
            new Masked_N_IniEquilibriumFunctional2D<T,Descriptor>, domain, fields );
}

template<typename T, template<class U> class Descriptor>
void setExternalVector (
        MultiBlockLattice2D<T,Descriptor>& lattice, int vectorStartsAt,
        T* externalVector, int numDimIs2, Box2D domain )
{
    Array<T,2> externalVectorVect(externalVector[0], externalVector[1]);
    applyProcessingFunctional (
            new SetExternalVectorFunctional2D<T,Descriptor>(vectorStartsAt, externalVectorVect),
            domain, lattice );
}


template<typename T, template<class U> class Descriptor>
void setPopulations( MultiBlockLattice2D<T,Descriptor>& lattice,
                     T* populations, int numDimIsQ, Box2D domain )
{
    PLB_PRECONDITION( numDimIsQ == Descriptor<T>::q );
    std::vector<T> populationVect(Descriptor<T>::q);
    for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
        populationVect[iPop] = populations[iPop];
    }
    applyProcessingFunctional (
            new N_IniConstPopulationsFunctional2D<T,Descriptor>(populationVect),
            domain, lattice );
}

template<typename T, template<class U> class Descriptor>
void setPopulations( MultiBlockLattice2D<T,Descriptor>& lattice,
                     MultiNTensorField2D<T>& populations,
                     Box2D domain )
{
    applyProcessingFunctional (
            new N_IniPopulationsFunctional2D<T,Descriptor>, domain,
            lattice, populations );
}

template<typename T, template<class U> class Descriptor>
void maskedSetPopulations( MultiBlockLattice2D<T,Descriptor>& lattice,
                           T* populations, int numDimIsQ,
                           MultiNTensorField2D<int>& mask, Box2D domain )
{
    PLB_PRECONDITION( numDimIsQ == Descriptor<T>::q );
    std::vector<T> populationVect(Descriptor<T>::q);
    for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
        populationVect[iPop] = populations[iPop];
    }
    applyProcessingFunctional (
            new Masked_N_IniConstPopulationsFunctional2D<T,Descriptor>(populationVect),
            domain, lattice, mask );
}

template<typename T, template<class U> class Descriptor>
void maskedSetPopulations( MultiBlockLattice2D<T,Descriptor>& lattice,
                           MultiNTensorField2D<T>& populations,
                           MultiNTensorField2D<int>& mask,
                           Box2D domain )
{
    applyProcessingFunctional (
            new Masked_N_IniPopulationsFunctional2D<T,Descriptor>, domain,
            lattice, populations, mask );
}

}  // namespace plb

#endif  // SWIG_LATTICE_INITIALIZER_WRAPPER_2D_HH
