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
#ifndef SWIG_LATTICE_INITIALIZER_WRAPPER_3D_HH
#define SWIG_LATTICE_INITIALIZER_WRAPPER_3D_HH

#include "plbWrapper/lattice/latticeInitializerWrapper3D.h"
#include "dataProcessors/dataInitializerWrapper3D.h"
#include "dataProcessors/dataInitializerFunctional3D.h"
#include "core/cell.h"
#include "atomicBlock/dataProcessorWrapper3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataField3D.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiDataField3D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"

namespace plb {

template<typename T, template<typename U> class Descriptor>
void pypalDefineDynamics( MultiBlockLattice3D<T,Descriptor>& lattice,
                          Box3D domain, Dynamics<T,Descriptor>* dynamics )
{
    applyProcessingFunctional (
            new InstantiateDynamicsFunctional3D<T,Descriptor> (
                dynamics->clone() ),
            domain, lattice );
}

template<typename T, template<typename U> class Descriptor>
void maskedDefineDynamics( MultiBlockLattice3D<T,Descriptor>& lattice,
                           MultiNTensorField3D<int>& mask,
                           Box3D domain, Dynamics<T,Descriptor>* dynamics )
{
    applyProcessingFunctional (
            new MaskedIniDynamicsFunctional3D<T,Descriptor> (
                dynamics->clone() ),
            domain, lattice, mask );
}

template<typename T, template<class U> class Descriptor>
void setBoundaryVelocity (
        MultiBlockLattice3D<T,Descriptor>& lattice,
        T* velocity, int numDimIs2, Box3D domain )
{
    Array<T,3> velocityVect(velocity[0], velocity[1], velocity[2]);
    applyProcessingFunctional (
            new SetConstBoundaryVelocityFunctional3D<T,Descriptor>(velocityVect), domain, lattice );
}

template<typename T, template<class U> class Descriptor>
void setBoundaryVelocity (
        MultiBlockLattice3D<T,Descriptor>& lattice,
        MultiNTensorField3D<T>& velocity,
        Box3D domain )
{
    applyProcessingFunctional (
            new N_IniBoundaryVelocityFunctional3D<T,Descriptor>, domain, lattice, velocity );
}

template<typename T, template<class U> class Descriptor>
void maskedSetBoundaryVelocity (
        MultiBlockLattice3D<T,Descriptor>& lattice,
        MultiNTensorField3D<T>& velocity,
        MultiNTensorField3D<int>& mask,
        Box3D domain )
{
    applyProcessingFunctional (
            new Masked_N_IniBoundaryVelocityFunctional3D<T,Descriptor>, domain, lattice, velocity, mask );
}

template<typename T, template<class U> class Descriptor>
void initializeAtEquilibrium (
        MultiBlockLattice3D<T,Descriptor>& lattice, T rho,
        T* velocity, int numDimIs2, Box3D domain )
{
    Array<T,3> velocityVect(velocity[0], velocity[1], velocity[2]);
    applyProcessingFunctional (
            new IniConstEquilibriumFunctional3D<T,Descriptor>(rho, velocityVect, (T)1), domain, lattice );
}

template<typename T, template<class U> class Descriptor>
void initializeAtEquilibrium( MultiBlockLattice3D<T,Descriptor>& lattice,
                              MultiNTensorField3D<T>& density,
                              MultiNTensorField3D<T>& velocity,
                              Box3D domain )
{
    std::vector<MultiBlock3D*> fields;
    fields.push_back(&lattice);
    fields.push_back(&density);
    fields.push_back(&velocity);
    applyProcessingFunctional (
            new N_IniEquilibriumFunctional3D<T,Descriptor>, domain, fields );
}

template<typename T, template<class U> class Descriptor>
void maskedInitializeAtEquilibrium( MultiBlockLattice3D<T,Descriptor>& lattice,
                                    MultiNTensorField3D<T>& density,
                                    MultiNTensorField3D<T>& velocity,
                                    MultiNTensorField3D<int>& mask,
                                    Box3D domain )
{
    std::vector<MultiBlock3D*> fields;
    fields.push_back(&lattice);
    fields.push_back(&density);
    fields.push_back(&velocity);
    fields.push_back(&mask);
    applyProcessingFunctional (
            new Masked_N_IniEquilibriumFunctional3D<T,Descriptor>, domain, fields );
}

template<typename T, template<class U> class Descriptor>
void setExternalVector (
        MultiBlockLattice3D<T,Descriptor>& lattice, int vectorStartsAt,
        T* externalVector, int numDimIs2, Box3D domain )
{
    Array<T,3> externalVectorVect(externalVector[0], externalVector[1], externalVector[2]);
    applyProcessingFunctional (
            new SetExternalVectorFunctional3D<T,Descriptor>(vectorStartsAt, externalVectorVect),
            domain, lattice );
}


template<typename T, template<class U> class Descriptor>
void setPopulations( MultiBlockLattice3D<T,Descriptor>& lattice,
                     T* populations, int numDimIsQ, Box3D domain )
{
    PLB_PRECONDITION( numDimIsQ == Descriptor<T>::q );
    std::vector<T> populationVect(Descriptor<T>::q);
    for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
        populationVect[iPop] = populations[iPop];
    }
    applyProcessingFunctional (
            new N_IniConstPopulationsFunctional3D<T,Descriptor>(populationVect),
            domain, lattice );
}

template<typename T, template<class U> class Descriptor>
void setPopulations( MultiBlockLattice3D<T,Descriptor>& lattice,
                     MultiNTensorField3D<T>& populations,
                     Box3D domain )
{
    applyProcessingFunctional (
            new N_IniPopulationsFunctional3D<T,Descriptor>, domain,
            lattice, populations );
}

template<typename T, template<class U> class Descriptor>
void maskedSetPopulations( MultiBlockLattice3D<T,Descriptor>& lattice,
                           T* populations, int numDimIsQ,
                           MultiNTensorField3D<int>& mask, Box3D domain )
{
    PLB_PRECONDITION( numDimIsQ == Descriptor<T>::q );
    std::vector<T> populationVect(Descriptor<T>::q);
    for (plint iPop=0; iPop<Descriptor<T>::q; ++iPop) {
        populationVect[iPop] = populations[iPop];
    }
    applyProcessingFunctional (
            new Masked_N_IniConstPopulationsFunctional3D<T,Descriptor>(populationVect),
            domain, lattice, mask );
}

template<typename T, template<class U> class Descriptor>
void maskedSetPopulations( MultiBlockLattice3D<T,Descriptor>& lattice,
                           MultiNTensorField3D<T>& populations,
                           MultiNTensorField3D<int>& mask,
                           Box3D domain )
{
    applyProcessingFunctional (
            new Masked_N_IniPopulationsFunctional3D<T,Descriptor>, domain,
            lattice, populations, mask );
}

}  // namespace plb

#endif  // SWIG_LATTICE_INITIALIZER_WRAPPER_3D_HH
