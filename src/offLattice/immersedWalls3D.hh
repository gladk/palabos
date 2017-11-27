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

#ifndef IMMERSED_WALLS_3D_HH
#define IMMERSED_WALLS_3D_HH

#include "core/globalDefs.h"
#include "core/array.h"
#include "atomicBlock/dataProcessingFunctional3D.h"
#include "atomicBlock/dataField3D.h"

#include "immersedWalls3D.h"

namespace plb {

/* ******** ReduceAxialTorqueImmersed3D ************************************ */

template<typename T>
ReduceAxialTorqueImmersed3D<T>::ReduceAxialTorqueImmersed3D (
        Array<T,3> const& center_, Array<T,3> const& unitaryAxis_, int reductionFlag_ )
    : center(center_),
      unitaryAxis(unitaryAxis_),
      sum_torque_ids (
            Array<plint,3> (
                this->getStatistics().subscribeSum(),
                this->getStatistics().subscribeSum(),
                this->getStatistics().subscribeSum() ) ),
      reductionFlag(reductionFlag_)
{ }

template<typename T>
void ReduceAxialTorqueImmersed3D<T>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    AtomicContainerBlock3D* container = dynamic_cast<AtomicContainerBlock3D*>(blocks[0]);
    PLB_ASSERT( container );

    ImmersedWallData3D<T>* wallData = 
        dynamic_cast<ImmersedWallData3D<T>*>( container->getData() );
    PLB_ASSERT(wallData);
    std::vector< Array<T,3> > const& vertices = wallData->vertices;
    std::vector< Array<T,3> > const& g = wallData->g;
    std::vector<int> const& flags = wallData->flags;
    Array<T,3> offset = wallData->offset;
    PLB_ASSERT( vertices.size()==g.size() );
    PLB_ASSERT( vertices.size()==flags.size() );

    for (pluint i=0; i<vertices.size(); ++i) {
        Array<T,3> vertex = vertices[i];
        if ( flags[i]==reductionFlag &&
             closedOpenContained(vertex, domain) )
        {
            Array<T,3> physVertex = vertex+offset;
            Array<T,3> r(physVertex-center);
            r -= dot(r,unitaryAxis)*unitaryAxis;
            Array<T,3> torque(crossProduct(r,g[i]));
            this->getStatistics().gatherSum(sum_torque_ids[0], torque[0]);
            this->getStatistics().gatherSum(sum_torque_ids[1], torque[1]);
            this->getStatistics().gatherSum(sum_torque_ids[2], torque[2]);
        }
    }
}

template<typename T>
ReduceAxialTorqueImmersed3D<T>* ReduceAxialTorqueImmersed3D<T>::clone() const {
    return new ReduceAxialTorqueImmersed3D<T>(*this);
}

template<typename T>
void ReduceAxialTorqueImmersed3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing; // Container Block.
}

template<typename T>
BlockDomain::DomainT ReduceAxialTorqueImmersed3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}

template<typename T>
Array<T,3> ReduceAxialTorqueImmersed3D<T>::getSumTorque() const {
    return Array<T,3> (
            this->getStatistics().getSum(sum_torque_ids[0]),
            this->getStatistics().getSum(sum_torque_ids[1]),
            this->getStatistics().getSum(sum_torque_ids[2]) );
}

/* ******** ReduceImmersedForce3D ************************************ */

template<typename T>
ReduceImmersedForce3D<T>::ReduceImmersedForce3D(int reductionFlag_)
    : sum_g_ids(
            Array<plint,3> (
                this->getStatistics().subscribeSum(),
                this->getStatistics().subscribeSum(),
                this->getStatistics().subscribeSum() ) ),
      reductionFlag(reductionFlag_)
{ }

template<typename T>
void ReduceImmersedForce3D<T>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    AtomicContainerBlock3D* container = dynamic_cast<AtomicContainerBlock3D*>(blocks[0]);
    PLB_ASSERT( container );

    ImmersedWallData3D<T>* wallData = 
        dynamic_cast<ImmersedWallData3D<T>*>( container->getData() );
    PLB_ASSERT(wallData);
    std::vector< Array<T,3> > const& vertices = wallData->vertices;
    std::vector< Array<T,3> > const& g = wallData->g;
    std::vector<int> const& flags = wallData->flags;
    PLB_ASSERT( vertices.size()==g.size() );
    PLB_ASSERT( vertices.size()==flags.size() );

    for (pluint i=0; i<vertices.size(); ++i) {
        Array<T,3> vertex = vertices[i];
        if ( flags[i]==reductionFlag &&
             closedOpenContained(vertex, domain) )
        {
            this->getStatistics().gatherSum(sum_g_ids[0], g[i][0]);
            this->getStatistics().gatherSum(sum_g_ids[1], g[i][1]);
            this->getStatistics().gatherSum(sum_g_ids[2], g[i][2]);
        }
    }
}

template<typename T>
ReduceImmersedForce3D<T>* ReduceImmersedForce3D<T>::clone() const {
    return new ReduceImmersedForce3D<T>(*this);
}

template<typename T>
void ReduceImmersedForce3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing; // Container Block.
}

template<typename T>
BlockDomain::DomainT ReduceImmersedForce3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}

template<typename T>
Array<T,3> ReduceImmersedForce3D<T>::getSumG() const {
    return Array<T,3> (
            this->getStatistics().getSum(sum_g_ids[0]),
            this->getStatistics().getSum(sum_g_ids[1]),
            this->getStatistics().getSum(sum_g_ids[2]) );
}

/* ******** ReduceImmersedArea3D ************************************ */

template<typename T>
ReduceImmersedArea3D<T>::ReduceImmersedArea3D(int reductionFlag_)
    : sum_area_id(this->getStatistics().subscribeSum()),
      reductionFlag(reductionFlag_)
{ }

template<typename T>
void ReduceImmersedArea3D<T>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    AtomicContainerBlock3D* container = dynamic_cast<AtomicContainerBlock3D*>(blocks[0]);
    PLB_ASSERT( container );

    ImmersedWallData3D<T>* wallData = 
        dynamic_cast<ImmersedWallData3D<T>*>( container->getData() );
    PLB_ASSERT(wallData);
    std::vector< Array<T,3> > const& vertices = wallData->vertices;
    std::vector<T> const& areas = wallData->areas;
    std::vector<int> const& flags = wallData->flags;
    PLB_ASSERT( vertices.size()==areas.size() );
    PLB_ASSERT( vertices.size()==flags.size() );

    for (pluint i=0; i<vertices.size(); ++i) {
        Array<T,3> vertex = vertices[i];
        if ( flags[i]==reductionFlag &&
             closedOpenContained(vertex, domain) )
        {
            this->getStatistics().gatherSum(sum_area_id, areas[i]);
        }
    }
}

template<typename T>
ReduceImmersedArea3D<T>* ReduceImmersedArea3D<T>::clone() const {
    return new ReduceImmersedArea3D<T>(*this);
}

template<typename T>
void ReduceImmersedArea3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing; // Container Block.
}

template<typename T>
BlockDomain::DomainT ReduceImmersedArea3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}

template<typename T>
T ReduceImmersedArea3D<T>::getSumArea() const {
    return this->getStatistics().getSum(sum_area_id);
}

/* ******** InamuroIteration3D ************************************ */

template<typename T, class VelFunction>
InamuroIteration3D<T,VelFunction>::InamuroIteration3D(VelFunction velFunction_, T tau_, bool incompressibleModel_)
    : velFunction(velFunction_),
      tau(tau_),
      incompressibleModel(incompressibleModel_)
{ }

template<typename T, class VelFunction>
void InamuroIteration3D<T,VelFunction>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==3 );
    ScalarField3D<T>* rhoBar = dynamic_cast<ScalarField3D<T>*>(blocks[0]);
    TensorField3D<T,3>* j = dynamic_cast<TensorField3D<T,3>*>(blocks[1]);
    AtomicContainerBlock3D* container = dynamic_cast<AtomicContainerBlock3D*>(blocks[2]);
    PLB_ASSERT( rhoBar );
    PLB_ASSERT( j );
    PLB_ASSERT( container );

    Dot3D ofsJ = computeRelativeDisplacement(*rhoBar, *j);

    ImmersedWallData3D<T>* wallData = 
        dynamic_cast<ImmersedWallData3D<T>*>( container->getData() );
    PLB_ASSERT(wallData);
    Array<T,3> absOffset = wallData->offset;

    std::vector< Array<T,3> > const& vertices = wallData->vertices;
    std::vector<T> const& areas = wallData->areas;
    PLB_ASSERT( vertices.size()==areas.size() );
    std::vector<Array<T,3> > deltaG(vertices.size());
    std::vector<Array<T,3> >& g = wallData->g;
    PLB_ASSERT( vertices.size()==g.size() );

    // In this iteration, the force is computed for every vertex.
    if (incompressibleModel) {
        for (pluint i=0; i<vertices.size(); ++i) {
            Array<T,3> const& vertex = vertices[i];
            Array<plint,3> intPos (
                    (plint)vertex[0], (plint)vertex[1], (plint)vertex[2] );
            Array<T,3> averageJ; averageJ.resetToZero();
            // Use the weighting function to compute the average momentum
            // and the average density on the surface vertex.
            // x   x . x   x
            for (plint dx=-1; dx<=+2; ++dx) {
                for (plint dy=-1; dy<=+2; ++dy) {
                    for (plint dz=-1; dz<=+2; ++dz) {
                        Array<plint,3> pos(intPos+Array<plint,3>(dx,dy,dz));
                        Array<T,3> nextJ = j->get(pos[0]+ofsJ.x, pos[1]+ofsJ.y, pos[2]+ofsJ.z);
                        Array<T,3> r(pos[0]-vertex[0],pos[1]-vertex[1],pos[2]-vertex[2]);
                        T W = inamuroDeltaFunction<T>().W(r);
                        averageJ += W*nextJ;
                    }
                }
            }
            //averageJ += (T)0.5*g[i];
            Array<T,3> wallVelocity = velFunction(vertex+absOffset);
            deltaG[i] = areas[i]*(wallVelocity-averageJ);
            g[i] += deltaG[i];
        }
    } else { // Compressible model.
        for (pluint i=0; i<vertices.size(); ++i) {
            Array<T,3> const& vertex = vertices[i];
            Array<plint,3> intPos (
                    (plint)vertex[0], (plint)vertex[1], (plint)vertex[2] );
            Array<T,3> averageJ; averageJ.resetToZero();
            T averageRhoBar = T();
            // Use the weighting function to compute the average momentum
            // and the average density on the surface vertex.
            // x   x . x   x
            for (plint dx=-1; dx<=+2; ++dx) {
                for (plint dy=-1; dy<=+2; ++dy) {
                    for (plint dz=-1; dz<=+2; ++dz) {
                        Array<plint,3> pos(intPos+Array<plint,3>(dx,dy,dz));
                        T nextRhoBar = rhoBar->get(pos[0], pos[1], pos[2]);
                        Array<T,3> nextJ = j->get(pos[0]+ofsJ.x, pos[1]+ofsJ.y, pos[2]+ofsJ.z);
                        Array<T,3> r(pos[0]-vertex[0],pos[1]-vertex[1],pos[2]-vertex[2]);
                        T W = inamuroDeltaFunction<T>().W(r);
                        averageJ += W*nextJ;
                        averageRhoBar += W*nextRhoBar;
                    }
                }
            }
            //averageJ += (T)0.5*g[i];
            Array<T,3> wallVelocity = velFunction(vertex+absOffset);
            deltaG[i] = areas[i]*((averageRhoBar+(T)1.)*wallVelocity-averageJ);
            //g[i] += deltaG[i];
            g[i] += deltaG[i]/((T)1.0+averageRhoBar);
        }
    }
    
    // In this iteration, the force is applied from every vertex to the grid nodes.
    for (pluint i=0; i<vertices.size(); ++i) {
        Array<T,3> const& vertex = vertices[i];
        Array<plint,3> intPos (
                (plint)vertex[0], (plint)vertex[1], (plint)vertex[2] );
        for (plint dx=-1; dx<=+2; ++dx) {
            for (plint dy=-1; dy<=+2; ++dy) {
                for (plint dz=-1; dz<=+2; ++dz) {
                    Array<plint,3> pos(intPos+Array<plint,3>(dx,dy,dz));
                    Array<T,3> nextJ = j->get(pos[0]+ofsJ.x, pos[1]+ofsJ.y, pos[2]+ofsJ.z);
                    Array<T,3> r(pos[0]-vertex[0],pos[1]-vertex[1],pos[2]-vertex[2]);
                    T W = inamuroDeltaFunction<T>().W(r);
                    nextJ += tau*W*deltaG[i];
                    j->get(pos[0]+ofsJ.x, pos[1]+ofsJ.y, pos[2]+ofsJ.z) = nextJ;
                }
            }
        }
    }
}

template<typename T, class VelFunction>
InamuroIteration3D<T,VelFunction>* InamuroIteration3D<T,VelFunction>::clone() const {
    return new InamuroIteration3D<T,VelFunction>(*this);
}

template<typename T, class VelFunction>
void InamuroIteration3D<T,VelFunction>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;          // RhoBar
    modified[1] = modif::staticVariables;  // J
    modified[2] = modif::nothing;          // Container Block with triangle data.
}

template<typename T, class VelFunction>
BlockDomain::DomainT InamuroIteration3D<T,VelFunction>::appliesTo() const {
    return BlockDomain::bulk;
}

/* ******** IndexedInamuroIteration3D ************************************ */

template<typename T, class VelFunction>
IndexedInamuroIteration3D<T,VelFunction>::IndexedInamuroIteration3D(VelFunction velFunction_, T tau_, bool incompressibleModel_)
    : velFunction(velFunction_),
      tau(tau_),
      incompressibleModel(incompressibleModel_)
{ }

template<typename T, class VelFunction>
void IndexedInamuroIteration3D<T,VelFunction>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==3 );
    ScalarField3D<T>* rhoBar = dynamic_cast<ScalarField3D<T>*>(blocks[0]);
    TensorField3D<T,3>* j = dynamic_cast<TensorField3D<T,3>*>(blocks[1]);
    AtomicContainerBlock3D* container = dynamic_cast<AtomicContainerBlock3D*>(blocks[2]);
    PLB_ASSERT( rhoBar );
    PLB_ASSERT( j );
    PLB_ASSERT( container );

    Dot3D ofsJ = computeRelativeDisplacement(*rhoBar, *j);

    ImmersedWallData3D<T>* wallData = 
        dynamic_cast<ImmersedWallData3D<T>*>( container->getData() );
    PLB_ASSERT(wallData);

    std::vector< Array<T,3> > const& vertices = wallData->vertices;
    std::vector<T> const& areas = wallData->areas;
    PLB_ASSERT( vertices.size()==areas.size() );
    std::vector<Array<T,3> > deltaG(vertices.size());
    std::vector<Array<T,3> >& g = wallData->g;
    PLB_ASSERT( vertices.size()==g.size() );
    std::vector<pluint> const& globalVertexIds = wallData->globalVertexIds;
    PLB_ASSERT( vertices.size()==globalVertexIds.size() );

    if (incompressibleModel) {
        for (pluint i=0; i<vertices.size(); ++i) {
            Array<T,3> const& vertex = vertices[i];
            Array<plint,3> intPos (
                    (plint)vertex[0], (plint)vertex[1], (plint)vertex[2] );
            Array<T,3> averageJ; averageJ.resetToZero();
            // x   x . x   x
            for (plint dx=-1; dx<=+2; ++dx) {
                for (plint dy=-1; dy<=+2; ++dy) {
                    for (plint dz=-1; dz<=+2; ++dz) {
                        Array<plint,3> pos(intPos+Array<plint,3>(dx,dy,dz));
                        Array<T,3> nextJ = j->get(pos[0]+ofsJ.x, pos[1]+ofsJ.y, pos[2]+ofsJ.z);
                        Array<T,3> r(pos[0]-vertex[0],pos[1]-vertex[1],pos[2]-vertex[2]);
                        T W = inamuroDeltaFunction<T>().W(r);
                        averageJ += W*nextJ;
                    }
                }
            }
            //averageJ += (T)0.5*g[i];
            Array<T,3> wallVelocity = velFunction(globalVertexIds[i]);
            deltaG[i] = areas[i]*(wallVelocity-averageJ);
            g[i] += deltaG[i];
        }
    } else { // Compressible model.
        for (pluint i=0; i<vertices.size(); ++i) {
            Array<T,3> const& vertex = vertices[i];
            Array<plint,3> intPos (
                    (plint)vertex[0], (plint)vertex[1], (plint)vertex[2] );
            Array<T,3> averageJ; averageJ.resetToZero();
            T averageRhoBar = T();
            // x   x . x   x
            for (plint dx=-1; dx<=+2; ++dx) {
                for (plint dy=-1; dy<=+2; ++dy) {
                    for (plint dz=-1; dz<=+2; ++dz) {
                        Array<plint,3> pos(intPos+Array<plint,3>(dx,dy,dz));
                        T nextRhoBar = rhoBar->get(pos[0], pos[1], pos[2]);
                        Array<T,3> nextJ = j->get(pos[0]+ofsJ.x, pos[1]+ofsJ.y, pos[2]+ofsJ.z);
                        Array<T,3> r(pos[0]-vertex[0],pos[1]-vertex[1],pos[2]-vertex[2]);
                        T W = inamuroDeltaFunction<T>().W(r);
                        averageJ += W*nextJ;
                        averageRhoBar += W*nextRhoBar;
                    }
                }
            }
            //averageJ += (T)0.5*g[i];
            Array<T,3> wallVelocity = velFunction(globalVertexIds[i]);
            deltaG[i] = areas[i]*((averageRhoBar+(T)1.)*wallVelocity-averageJ);
            //g[i] += deltaG[i];
            g[i] += deltaG[i]/((T)1.0+averageRhoBar);
        }
    }
    
    for (pluint i=0; i<vertices.size(); ++i) {
        Array<T,3> const& vertex = vertices[i];
        Array<plint,3> intPos (
                (plint)vertex[0], (plint)vertex[1], (plint)vertex[2] );
        for (plint dx=-1; dx<=+2; ++dx) {
            for (plint dy=-1; dy<=+2; ++dy) {
                for (plint dz=-1; dz<=+2; ++dz) {
                    Array<plint,3> pos(intPos+Array<plint,3>(dx,dy,dz));
                    Array<T,3> nextJ = j->get(pos[0]+ofsJ.x, pos[1]+ofsJ.y, pos[2]+ofsJ.z);
                    Array<T,3> r(pos[0]-vertex[0],pos[1]-vertex[1],pos[2]-vertex[2]);
                    T W = inamuroDeltaFunction<T>().W(r);
                    nextJ += tau*W*deltaG[i];
                    j->get(pos[0]+ofsJ.x, pos[1]+ofsJ.y, pos[2]+ofsJ.z) = nextJ;
                }
            }
        }
    }
}

template<typename T, class VelFunction>
IndexedInamuroIteration3D<T,VelFunction>* IndexedInamuroIteration3D<T,VelFunction>::clone() const {
    return new IndexedInamuroIteration3D<T,VelFunction>(*this);
}

template<typename T, class VelFunction>
void IndexedInamuroIteration3D<T,VelFunction>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;          // RhoBar
    modified[1] = modif::staticVariables;  // J
    modified[2] = modif::nothing;          // Container Block with triangle data.
}

template<typename T, class VelFunction>
BlockDomain::DomainT IndexedInamuroIteration3D<T,VelFunction>::appliesTo() const {
    return BlockDomain::bulk;
}

/* ******** ConstVelInamuroIteration3D ************************************ */

template<typename T>
ConstVelInamuroIteration3D<T>::ConstVelInamuroIteration3D(Array<T,3> const& wallVelocity_, T tau_, bool incompressibleModel_)
    : wallVelocity(wallVelocity_),
      tau(tau_),
      incompressibleModel(incompressibleModel_)
{ }

template<typename T>
void ConstVelInamuroIteration3D<T>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==3 );
    ScalarField3D<T>* rhoBar = dynamic_cast<ScalarField3D<T>*>(blocks[0]);
    TensorField3D<T,3>* j = dynamic_cast<TensorField3D<T,3>*>(blocks[1]);
    AtomicContainerBlock3D* container = dynamic_cast<AtomicContainerBlock3D*>(blocks[2]);
    PLB_ASSERT( rhoBar );
    PLB_ASSERT( j );
    PLB_ASSERT( container );

    Dot3D ofsJ = computeRelativeDisplacement(*rhoBar, *j);

    ImmersedWallData3D<T>* wallData = 
        dynamic_cast<ImmersedWallData3D<T>*>( container->getData() );
    PLB_ASSERT(wallData);
    std::vector< Array<T,3> > const& vertices = wallData->vertices;
    std::vector<T> const& areas = wallData->areas;
    PLB_ASSERT( vertices.size()==areas.size() );
    std::vector<Array<T,3> > deltaG(vertices.size());
    std::vector<Array<T,3> >& g = wallData->g;
    PLB_ASSERT( vertices.size()==g.size() );

    if (incompressibleModel) {
        for (pluint i=0; i<vertices.size(); ++i) {
            Array<T,3> const& vertex = vertices[i];
            Array<plint,3> intPos (
                    (plint)vertex[0], (plint)vertex[1], (plint)vertex[2] );
            Array<T,3> averageJ; averageJ.resetToZero();
            // x   x . x   x
            for (plint dx=-1; dx<=+2; ++dx) {
                for (plint dy=-1; dy<=+2; ++dy) {
                    for (plint dz=-1; dz<=+2; ++dz) {
                        Array<plint,3> pos(intPos+Array<plint,3>(dx,dy,dz));
                        Array<T,3> nextJ = j->get(pos[0]+ofsJ.x, pos[1]+ofsJ.y, pos[2]+ofsJ.z);
                        Array<T,3> r(pos[0]-vertex[0],pos[1]-vertex[1],pos[2]-vertex[2]);
                        T W = inamuroDeltaFunction<T>().W(r);
                        averageJ += W*nextJ;
                    }
                }
            }
            //averageJ += (T)0.5*g[i];
            deltaG[i] = areas[i]*(wallVelocity-averageJ);
            g[i] += deltaG[i];
        }
    } else { // Compressible model.
        for (pluint i=0; i<vertices.size(); ++i) {
            Array<T,3> const& vertex = vertices[i];
            Array<plint,3> intPos (
                    (plint)vertex[0], (plint)vertex[1], (plint)vertex[2] );
            Array<T,3> averageJ; averageJ.resetToZero();
            T averageRhoBar = T();
            // x   x . x   x
            for (plint dx=-1; dx<=+2; ++dx) {
                for (plint dy=-1; dy<=+2; ++dy) {
                    for (plint dz=-1; dz<=+2; ++dz) {
                        Array<plint,3> pos(intPos+Array<plint,3>(dx,dy,dz));
                        T nextRhoBar = rhoBar->get(pos[0], pos[1], pos[2]);
                        Array<T,3> nextJ = j->get(pos[0]+ofsJ.x, pos[1]+ofsJ.y, pos[2]+ofsJ.z);
                        Array<T,3> r(pos[0]-vertex[0],pos[1]-vertex[1],pos[2]-vertex[2]);
                        T W = inamuroDeltaFunction<T>().W(r);
                        averageJ += W*nextJ;
                        averageRhoBar += W*nextRhoBar;
                    }
                }
            }
            //averageJ += (T)0.5*g[i];
            deltaG[i] = areas[i]*((averageRhoBar+(T)1.)*wallVelocity-averageJ);
            //g[i] += deltaG[i];
            g[i] += deltaG[i]/((T)1.0+averageRhoBar);
        }
    }
    
    for (pluint i=0; i<vertices.size(); ++i) {
        Array<T,3> const& vertex = vertices[i];
        Array<plint,3> intPos (
                (plint)vertex[0], (plint)vertex[1], (plint)vertex[2] );
        for (plint dx=-1; dx<=+2; ++dx) {
            for (plint dy=-1; dy<=+2; ++dy) {
                for (plint dz=-1; dz<=+2; ++dz) {
                    Array<plint,3> pos(intPos+Array<plint,3>(dx,dy,dz));
                    Array<T,3> nextJ = j->get(pos[0]+ofsJ.x, pos[1]+ofsJ.y, pos[2]+ofsJ.z);
                    Array<T,3> r(pos[0]-vertex[0],pos[1]-vertex[1],pos[2]-vertex[2]);
                    T W = inamuroDeltaFunction<T>().W(r);
                    nextJ += tau*W*deltaG[i];
                    j->get(pos[0]+ofsJ.x, pos[1]+ofsJ.y, pos[2]+ofsJ.z) = nextJ;
                }
            }
        }
    }
}

template<typename T>
ConstVelInamuroIteration3D<T>* ConstVelInamuroIteration3D<T>::clone() const {
    return new ConstVelInamuroIteration3D<T>(*this);
}

template<typename T>
void ConstVelInamuroIteration3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;          // RhoBar
    modified[1] = modif::staticVariables;  // J
    modified[2] = modif::nothing;          // Container Block with triangle data.
}

template<typename T>
BlockDomain::DomainT ConstVelInamuroIteration3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}

/* ******** ComputeImmersedBoundaryForce3D ************************************ */

template<typename T>
void ComputeImmersedBoundaryForce3D<T>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==2 );
    TensorField3D<T,3>* force = dynamic_cast<TensorField3D<T,3>*>(blocks[0]);
    AtomicContainerBlock3D* container = dynamic_cast<AtomicContainerBlock3D*>(blocks[1]);
    PLB_ASSERT(force);
    PLB_ASSERT(container);

    ImmersedWallData3D<T>* wallData = 
        dynamic_cast<ImmersedWallData3D<T>*>( container->getData() );
    PLB_ASSERT(wallData);
    std::vector< Array<T,3> > const& vertices = wallData->vertices;
    std::vector<Array<T,3> >& g = wallData->g;
    PLB_ASSERT( vertices.size()==g.size() );

    for (plint iX = 0; iX < force->getNx(); iX++) {
        for (plint iY = 0; iY < force->getNy(); iY++) {
            for (plint iZ = 0; iZ < force->getNz(); iZ++) {
                force->get(iX, iY, iZ).resetToZero();
            }
        }
    }

    for (pluint i=0; i<vertices.size(); ++i) {
        Array<T,3> const& vertex = vertices[i];
        Array<plint,3> intPos (
                (plint)vertex[0], (plint)vertex[1], (plint)vertex[2] );
        for (plint dx=-1; dx<=+2; ++dx) {
            for (plint dy=-1; dy<=+2; ++dy) {
                for (plint dz=-1; dz<=+2; ++dz) {
                    Array<plint,3> pos(intPos+Array<plint,3>(dx,dy,dz));
                    Array<T,3> r(pos[0]-vertex[0],pos[1]-vertex[1],pos[2]-vertex[2]);
                    T W = inamuroDeltaFunction<T>().W(r);
                    force->get(pos[0], pos[1], pos[2]) += W*g[i];
                }
            }
        }
    }
}

template<typename T>
ComputeImmersedBoundaryForce3D<T>* ComputeImmersedBoundaryForce3D<T>::clone() const {
    return new ComputeImmersedBoundaryForce3D<T>(*this);
}

template<typename T>
void ComputeImmersedBoundaryForce3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;  // Force
    modified[1] = modif::nothing;          // Container Block with triangle data.
}

template<typename T>
BlockDomain::DomainT ComputeImmersedBoundaryForce3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** InstantiateImmersedWallData3D ************************************ */

template<typename T>
InstantiateImmersedWallData3D<T>::InstantiateImmersedWallData3D (
            std::vector< Array<T,3> > const& vertices_,
            std::vector<T> const& areas_,
            std::vector< Array<T,3> > const& normals_)
    : vertices(vertices_),
      areas(areas_),
      normals(normals_)
{
    PLB_ASSERT(vertices.size() == areas.size());
    PLB_ASSERT(normals.size()==0 || normals.size() == areas.size());
}

template<typename T>
void InstantiateImmersedWallData3D<T>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    AtomicContainerBlock3D* container = dynamic_cast<AtomicContainerBlock3D*>(blocks[0]);
    PLB_ASSERT( container );
    bool useNormals = normals.size()>0;
    Dot3D location = container->getLocation();
    Array<T,3> offset(location.x,location.y,location.z);
    ImmersedWallData3D<T>* wallData = new ImmersedWallData3D<T>;
    Box3D extendedEnvelope(domain.enlarge(2));

    for (pluint i=0; i<vertices.size(); ++i) {
        Array<T,3> vertex = vertices[i]-offset;
        // Vertices which are close to the boundaries of the extendedEnvelope
        // are irrelevant, because they will act upon the bulk of the computational
        // domain through an Inamuro kernel, which at this distance is close to zero.
        // It is therefore OK, numerically speaking to exclude an epsilon-margin close
        // to these boundaries. Plus, it is required for technical reasons, because if
        // later on we pass across the boundaries of the extendedEnvelope because
        // of roundoff errors, the code will crash.
        static const T epsilon = 1.e-4;
        if (contained(vertex, extendedEnvelope, epsilon)) {
            wallData->vertices.push_back(vertex);
            wallData->areas.push_back(areas[i]);
            if (useNormals) {
                wallData->normals.push_back(normals[i]);
            }
            wallData->g.push_back(Array<T,3>((T)0.,(T)0.,(T)0.));
            wallData->globalVertexIds.push_back(i);
        }
    }
    wallData->flags = std::vector<int>(wallData->vertices.size(), 0);
    wallData->offset = offset;
    container->setData(wallData);
}

template<typename T>
InstantiateImmersedWallData3D<T>* InstantiateImmersedWallData3D<T>::clone() const {
    return new InstantiateImmersedWallData3D<T>(*this);
}

template<typename T>
void InstantiateImmersedWallData3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;  // Container Block with triangle data.
}

template<typename T>
BlockDomain::DomainT InstantiateImmersedWallData3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** InstantiateImmersedWallDataWithTagging3D ************************************ */

template<typename T>
InstantiateImmersedWallDataWithTagging3D<T>::InstantiateImmersedWallDataWithTagging3D (
            std::vector< Array<T,3> > const& vertices_,
            std::vector<T> const& areas_, int fluidFlag_ )
    : vertices(vertices_),
      areas(areas_),
      fluidFlag(fluidFlag_)
{
    PLB_ASSERT(vertices.size() == areas.size());
}

template<typename T>
void InstantiateImmersedWallDataWithTagging3D<T>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==2 );
    AtomicContainerBlock3D* container = dynamic_cast<AtomicContainerBlock3D*>(blocks[0]);
    PLB_ASSERT( container );
    Dot3D location = container->getLocation();
    Array<T,3> offset(location.x,location.y,location.z);

    ScalarField3D<int>* flagMatrix = dynamic_cast<ScalarField3D<int>*>(blocks[1]);
    PLB_ASSERT(flagMatrix);
    Dot3D ofsFlag = computeRelativeDisplacement(*container, *flagMatrix);
    Array<plint,3> flagDispl(ofsFlag.x,ofsFlag.y,ofsFlag.z);

    ImmersedWallData3D<T>* wallData = new ImmersedWallData3D<T>;
    Box3D extendedEnvelope(domain.enlarge(2));

    for (pluint i=0; i<vertices.size(); ++i) {
        Array<T,3> vertex = vertices[i]-offset;
        // Vertices which are close to the boundaries of the extendedEnvelope
        // are irrelevant, because they will act upon the bulk of the computational
        // domain through an Inamuro kernel, which at this distance is close to zero.
        // It is therefore OK, numerically speaking to exclude an epsilon-margin close
        // to these boundaries. Plus, it is required for technical reasons, because if
        // later on we pass across the boundaries of the extendedEnvelope because
        // of roundoff errors, the code will crash.
        static const T epsilon = 1.e-4;
        if (contained(vertex, extendedEnvelope, epsilon)) {
            wallData->vertices.push_back(vertex);
            wallData->areas.push_back(areas[i]);
            wallData->g.push_back(Array<T,3>((T)0.,(T)0.,(T)0.));
            wallData->globalVertexIds.push_back(i);
            Array<plint,3> intPos (
                    (plint)vertex[0], (plint)vertex[1], (plint)vertex[2] );
            bool hasFluidNeighbor=false;
            for (plint dx=-1; dx<=+2; ++dx) {
                for (plint dy=-1; dy<=+2; ++dy) {
                    for (plint dz=-1; dz<=+2; ++dz) {
                        Array<plint,3> pos(intPos+Array<plint,3>(dx,dy,dz)+flagDispl);
                        if(flagMatrix->get(pos[0],pos[1],pos[2])==fluidFlag) {
                            hasFluidNeighbor=true;
                        }
                    }
                }
            }
            if (hasFluidNeighbor) {
                wallData->flags.push_back(0);
            }
            else {
                wallData->flags.push_back(1);
            }
        }
    }
    wallData->offset = offset;
    container->setData(wallData);
}

template<typename T>
InstantiateImmersedWallDataWithTagging3D<T>* InstantiateImmersedWallDataWithTagging3D<T>::clone() const {
    return new InstantiateImmersedWallDataWithTagging3D<T>(*this);
}

template<typename T>
void InstantiateImmersedWallDataWithTagging3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;  // Container Block with triangle data.
    modified[1] = modif::nothing;  // Flag matrix.
}

template<typename T>
BlockDomain::DomainT InstantiateImmersedWallDataWithTagging3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}

/* ******** InstantiateImmersedWallDataWithIndexedTagging3D ************************************ */

template<typename T>
InstantiateImmersedWallDataWithIndexedTagging3D<T>::InstantiateImmersedWallDataWithIndexedTagging3D (
            std::vector< Array<T,3> > const& vertices_,
            std::vector<T> const& areas_, std::vector<int> const& flags_ )
    : vertices(vertices_),
      areas(areas_),
      flags(flags_)
{
    PLB_ASSERT(vertices.size() == areas.size());
    PLB_ASSERT(vertices.size() == flags.size());
}

template<typename T>
void InstantiateImmersedWallDataWithIndexedTagging3D<T>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    AtomicContainerBlock3D* container = dynamic_cast<AtomicContainerBlock3D*>(blocks[0]);
    PLB_ASSERT( container );
    Dot3D location = container->getLocation();
    Array<T,3> offset(location.x,location.y,location.z);

    ImmersedWallData3D<T>* wallData = new ImmersedWallData3D<T>;
    Box3D extendedEnvelope(domain.enlarge(2));

    for (pluint i=0; i<vertices.size(); ++i) {
        Array<T,3> vertex = vertices[i]-offset;
        // Vertices which are close to the boundaries of the extendedEnvelope
        // are irrelevant, because they will act upon the bulk of the computational
        // domain through an Inamuro kernel, which at this distance is close to zero.
        // It is therefore OK, numerically speaking to exclude an epsilon-margin close
        // to these boundaries. Plus, it is required for technical reasons, because if
        // later on we pass across the boundaries of the extendedEnvelope because
        // of roundoff errors, the code will crash.
        static const T epsilon = 1.e-4;
        if (contained(vertex, extendedEnvelope, epsilon)) {
            wallData->vertices.push_back(vertex);
            wallData->areas.push_back(areas[i]);
            wallData->g.push_back(Array<T,3>((T)0.,(T)0.,(T)0.));
            wallData->flags.push_back(flags[i]);
            wallData->globalVertexIds.push_back(i);
        }
    }
    wallData->offset = offset;
    container->setData(wallData);
}

template<typename T>
InstantiateImmersedWallDataWithIndexedTagging3D<T>* InstantiateImmersedWallDataWithIndexedTagging3D<T>::clone() const {
    return new InstantiateImmersedWallDataWithIndexedTagging3D<T>(*this);
}

template<typename T>
void InstantiateImmersedWallDataWithIndexedTagging3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;  // Container Block with triangle data.
}

template<typename T>
BlockDomain::DomainT InstantiateImmersedWallDataWithIndexedTagging3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** InstantiateSurfaceBlockData3D ************************************ */

template<typename T>
InstantiateSurfaceBlockData3D<T>::InstantiateSurfaceBlockData3D(plint envelopeWidth_, std::vector< Array<T,3> > const& vertices_)
    : envelopeWidth(envelopeWidth_),
      vertices(vertices_)
{
    PLB_ASSERT(envelopeWidth >= 1);
}

template<typename T>
void InstantiateSurfaceBlockData3D<T>::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks)
{
    PLB_PRECONDITION(blocks.size()==1);
    AtomicContainerBlock3D *container = dynamic_cast<AtomicContainerBlock3D*>(blocks[0]);
    PLB_ASSERT(container);

    Dot3D location = container->getLocation();
    Array<T,3> offset(location.x, location.y, location.z);

    SurfaceBlockData3D<T> *surfaceData = new SurfaceBlockData3D<T>;
    Box3D extendedEnvelope(domain.enlarge(envelopeWidth));

    for (pluint i=0; i<vertices.size(); ++i) {
        Array<T,3> vertex = vertices[i]-offset;
        if (containedInclusive(vertex, extendedEnvelope)) {
            surfaceData->vertices.push_back(vertex);
        }
    }
    surfaceData->offset = offset;
    container->setData(surfaceData);
}

template<typename T>
InstantiateSurfaceBlockData3D<T> *InstantiateSurfaceBlockData3D<T>::clone() const
{
    return new InstantiateSurfaceBlockData3D<T>(*this);
}

template<typename T>
void InstantiateSurfaceBlockData3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::staticVariables;  // Container Block with triangle data.
}

template<typename T>
BlockDomain::DomainT InstantiateSurfaceBlockData3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}


/* ******** SurfaceOnLattice3D *************************************** */

template<typename T, typename U>
SurfaceOnLattice3D<T,U>::SurfaceOnLattice3D(U value_, plint envelopeWidth_)
    : value(value_),
      envelopeWidth(envelopeWidth_)
{
    PLB_ASSERT(envelopeWidth >= 1);
}

template<typename T, typename U>
void SurfaceOnLattice3D<T,U>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION(blocks.size() == 2);
    ScalarField3D<U> *surfaceOnLattice = dynamic_cast<ScalarField3D<U>*>(blocks[0]);
    AtomicContainerBlock3D *container = dynamic_cast<AtomicContainerBlock3D*>(blocks[1]);
    PLB_ASSERT(surfaceOnLattice);
    PLB_ASSERT(container);

    SurfaceBlockData3D<T> *surfaceData = dynamic_cast<SurfaceBlockData3D<T>*>(container->getData());
    PLB_ASSERT(surfaceData);

    std::vector< Array<T,3> > const& vertices = surfaceData->vertices;

    plint w = envelopeWidth;
    for (pluint i=0; i<vertices.size(); ++i) {
        Array<T,3> const& vertex = vertices[i];
        Array<plint,3> intPos (
                (plint)vertex[0], (plint)vertex[1], (plint)vertex[2] );
        // x  x  x . x  x  x
        for (plint dx=-(w-1); dx<=+w; ++dx) {
            for (plint dy=-(w-1); dy<=+w; ++dy) {
                for (plint dz=-(w-1); dz<=+w; ++dz) {
                    Array<plint,3> pos(intPos+Array<plint,3>(dx,dy,dz));
                    if (contained(pos[0], pos[1], pos[2], domain)) {
                        surfaceOnLattice->get(pos[0], pos[1], pos[2]) = value;
                    }
                }
            }
        }
    }
}

template<typename T, typename U>
SurfaceOnLattice3D<T,U>* SurfaceOnLattice3D<T,U>::clone() const
{
    return new SurfaceOnLattice3D<T,U>(*this);
}

template<typename T, typename U>
void SurfaceOnLattice3D<T,U>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::staticVariables;  // Surface values (ScalarField3D<U>).
    modified[1] = modif::nothing;          // Container Block with triangle data.
}

template<typename T, typename U>
BlockDomain::DomainT SurfaceOnLattice3D<T,U>::appliesTo() const
{
    return BlockDomain::bulk;
}


/* ******** SurfaceOnLattice3D_N *************************************** */

template<typename T, typename U>
SurfaceOnLattice3D_N<T,U>::SurfaceOnLattice3D_N(U value_, plint envelopeWidth_)
    : value(value_),
      envelopeWidth(envelopeWidth_)
{
    PLB_ASSERT(envelopeWidth >= 1);
}

template<typename T, typename U>
void SurfaceOnLattice3D_N<T,U>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION(blocks.size() == 2);
    NTensorField3D<U> *surfaceOnLattice = dynamic_cast<NTensorField3D<U>*>(blocks[0]);
    AtomicContainerBlock3D *container = dynamic_cast<AtomicContainerBlock3D*>(blocks[1]);
    PLB_ASSERT(surfaceOnLattice);
    PLB_ASSERT(container);

    SurfaceBlockData3D<T> *surfaceData = dynamic_cast<SurfaceBlockData3D<T>*>(container->getData());
    PLB_ASSERT(surfaceData);

    std::vector< Array<T,3> > const& vertices = surfaceData->vertices;

    plint w = envelopeWidth;
    for (pluint i=0; i<vertices.size(); ++i) {
        Array<T,3> const& vertex = vertices[i];
        Array<plint,3> intPos (
                (plint)vertex[0], (plint)vertex[1], (plint)vertex[2] );
        // x  x  x . x  x  x
        for (plint dx=-(w-1); dx<=+w; ++dx) {
            for (plint dy=-(w-1); dy<=+w; ++dy) {
                for (plint dz=-(w-1); dz<=+w; ++dz) {
                    Array<plint,3> pos(intPos+Array<plint,3>(dx,dy,dz));
                    if (contained(pos[0], pos[1], pos[2], domain)) {
                        *surfaceOnLattice->get(pos[0], pos[1], pos[2]) = value;
                    }
                }
            }
        }
    }
}

template<typename T, typename U>
SurfaceOnLattice3D_N<T,U>* SurfaceOnLattice3D_N<T,U>::clone() const
{
    return new SurfaceOnLattice3D_N<T,U>(*this);
}

template<typename T, typename U>
void SurfaceOnLattice3D_N<T,U>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::staticVariables;  // Surface values (NTensorField3D<U>).
    modified[1] = modif::nothing;          // Container Block with triangle data.
}

template<typename T, typename U>
BlockDomain::DomainT SurfaceOnLattice3D_N<T,U>::appliesTo() const
{
    return BlockDomain::bulk;
}


/* ******** ResetForceStatistics3D ************************************ */

template<typename T>
void ResetForceStatistics3D<T>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==1 );
    AtomicContainerBlock3D* container = dynamic_cast<AtomicContainerBlock3D*>(blocks[0]);
    PLB_ASSERT( container );

    ImmersedWallData3D<T>* wallData = 
        dynamic_cast<ImmersedWallData3D<T>*>( container->getData() );
    PLB_ASSERT(wallData);

    std::vector<Array<T,3> >& g = wallData->g;

    for (pluint i = 0; i < g.size(); i++) {
        g[i].resetToZero();
    }
}

template<typename T>
ResetForceStatistics3D<T>* ResetForceStatistics3D<T>::clone() const
{
    return new ResetForceStatistics3D<T>(*this);
}

template<typename T>
void ResetForceStatistics3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;  // Container Block with triangle data.
}

template<typename T>
BlockDomain::DomainT ResetForceStatistics3D<T>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** RecomputeImmersedForce3D ************************************ */

template<typename T, template<typename U> class Descriptor, class NormalFunction>
RecomputeImmersedForce3D<T,Descriptor,NormalFunction>::RecomputeImmersedForce3D(
        NormalFunction normalFunction_, T omega_, T densityOffset_,
        bool incompressibleModel_)
    : normalFunction(normalFunction_),
      omega(omega_),
      rho0(densityOffset_),
      incompressibleModel(incompressibleModel_)
{
    PLB_ASSERT(densityOffset_ > (T) 0);
}

template<typename T, template<typename U> class Descriptor, class NormalFunction>
void RecomputeImmersedForce3D<T,Descriptor,NormalFunction>::processGenericBlocks(
        Box3D domain, std::vector<AtomicBlock3D*> blocks)
{
    PLB_PRECONDITION(blocks.size() == 3);

    ScalarField3D<T>* rhoBar = dynamic_cast<ScalarField3D<T>*>(blocks[0]);
    TensorField3D<T,SymmetricTensorImpl<T,3>::n>* PiNeq = dynamic_cast<TensorField3D<T,SymmetricTensorImpl<T,3>::n>*>(blocks[1]);
    AtomicContainerBlock3D* container = dynamic_cast<AtomicContainerBlock3D*>(blocks[2]);
    PLB_ASSERT(rhoBar);
    PLB_ASSERT(PiNeq);
    PLB_ASSERT(container);

    Dot3D ofsPN = computeRelativeDisplacement(*rhoBar, *PiNeq);

    ImmersedWallData3D<T>* wallData = 
        dynamic_cast<ImmersedWallData3D<T>*>( container->getData() );
    PLB_ASSERT(wallData);

    std::vector< Array<T,3> > const& vertices = wallData->vertices;
    std::vector<T> const& areas = wallData->areas;
    PLB_ASSERT( vertices.size()==areas.size() );
    std::vector<Array<T,3> >& g = wallData->g;
    PLB_ASSERT( vertices.size()==g.size() );
    std::vector<pluint> const& globalVertexIds = wallData->globalVertexIds;
    PLB_ASSERT( vertices.size()==globalVertexIds.size() );

    for (pluint i=0; i<vertices.size(); ++i) {
        Array<T,3> normal = normalFunction(globalVertexIds[i]);

        // Interpolate rhoBar and PiNeq on the vertex position.

        Array<T,3> const& vertex = vertices[i];
        Array<plint,3> intPos (
                (plint)vertex[0], (plint)vertex[1], (plint)vertex[2] );

        T averageRhoBar = 0.0;
        Array<T,SymmetricTensorImpl<T,3>::n> averagePiNeq;
        averagePiNeq.resetToZero();

        // x   x . x   x
        for (plint dx=-1; dx<=+2; ++dx) {
            for (plint dy=-1; dy<=+2; ++dy) {
                for (plint dz=-1; dz<=+2; ++dz) {
                    Array<plint,3> pos(intPos+Array<plint,3>(dx,dy,dz));
                    T nextRhoBar = rhoBar->get(pos[0], pos[1], pos[2]);
                    Array<T,SymmetricTensorImpl<T,3>::n>& nextPiNeq = PiNeq->get(pos[0]+ofsPN.x, pos[1]+ofsPN.y, pos[2]+ofsPN.z);
                    Array<T,3> r(pos[0]-vertex[0],pos[1]-vertex[1],pos[2]-vertex[2]);
                    T W = inamuroDeltaFunction<T>().W(r);
                    averageRhoBar += W * nextRhoBar;
                    averagePiNeq += W * nextPiNeq;
                }
            }
        }

        // Compute the force on the fluid at the vertex position.

        T averageRho = Descriptor<T>::fullRho(averageRhoBar);

        Array<T,3> averagePi_n;
        SymmetricTensorImpl<T,3>::matVectMult(averagePiNeq, normal, averagePi_n);

        // Here we want the force acting on the fluid from the solid. "normal" points towards the fluid,
        // this is why the minus sign in front of the area is needed.
        if (incompressibleModel) {
            g[i] = -areas[i] * (-(averageRho-rho0)*Descriptor<T>::cs2*normal +
                    (omega/(T)2.-(T)1.)*averagePi_n);  // Incompressible vision
        } else {
            g[i] = -areas[i] * (-(averageRho-rho0)*Descriptor<T>::cs2*normal +
                    Descriptor<T>::invRho(averageRhoBar)*(omega/(T)2.-(T)1.)*averagePi_n);  // Compressible vision
        }
    }
}

template<typename T, template<typename U> class Descriptor, class NormalFunction>
RecomputeImmersedForce3D<T,Descriptor,NormalFunction>* RecomputeImmersedForce3D<T,Descriptor,NormalFunction>::clone() const
{
    return new RecomputeImmersedForce3D<T,Descriptor,NormalFunction>(*this);
}

template<typename T, template<typename U> class Descriptor, class NormalFunction>
void RecomputeImmersedForce3D<T,Descriptor,NormalFunction>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;   // RhoBar
    modified[1] = modif::nothing;   // PiNeq
    modified[2] = modif::nothing;   // Container with triangle data
}

template<typename T, template<typename U> class Descriptor, class NormalFunction>
BlockDomain::DomainT RecomputeImmersedForce3D<T,Descriptor,NormalFunction>::appliesTo() const
{
    return BlockDomain::bulk;
}


/* ******** OpenSurfaceImmersedForce3D ************************************ */

template<typename T, template<typename U> class Descriptor, class NormalFunction>
OpenSurfaceImmersedForce3D<T,Descriptor,NormalFunction>::OpenSurfaceImmersedForce3D(
        NormalFunction normalFunction_, T omega_, T densityOffset_,
        bool incompressibleModel_)
    : normalFunction(normalFunction_),
      omega(omega_),
      rho0(densityOffset_),
      incompressibleModel(incompressibleModel_)
{
    PLB_ASSERT(densityOffset_ > (T) 0);
}

template<typename T, template<typename U> class Descriptor, class NormalFunction>
void OpenSurfaceImmersedForce3D<T,Descriptor,NormalFunction>::processGenericBlocks(
        Box3D domain, std::vector<AtomicBlock3D*> blocks)
{
    PLB_PRECONDITION(blocks.size() == 3);

    ScalarField3D<T>* rhoBar = dynamic_cast<ScalarField3D<T>*>(blocks[0]);
    TensorField3D<T,SymmetricTensorImpl<T,3>::n>* PiNeq = dynamic_cast<TensorField3D<T,SymmetricTensorImpl<T,3>::n>*>(blocks[1]);
    AtomicContainerBlock3D* container = dynamic_cast<AtomicContainerBlock3D*>(blocks[2]);
    PLB_ASSERT(rhoBar);
    PLB_ASSERT(PiNeq);
    PLB_ASSERT(container);

    Dot3D ofsPN = computeRelativeDisplacement(*rhoBar, *PiNeq);

    ImmersedWallData3D<T>* wallData = 
        dynamic_cast<ImmersedWallData3D<T>*>( container->getData() );
    PLB_ASSERT(wallData);

    std::vector< Array<T,3> > const& vertices = wallData->vertices;
    std::vector<T> const& areas = wallData->areas;
    PLB_ASSERT( vertices.size()==areas.size() );
    std::vector<Array<T,3> >& g = wallData->g;
    PLB_ASSERT( vertices.size()==g.size() );
    std::vector<pluint> const& globalVertexIds = wallData->globalVertexIds;
    PLB_ASSERT( vertices.size()==globalVertexIds.size() );

    for (pluint i=0; i<vertices.size(); ++i) {
        Array<T,3> normal = normalFunction(globalVertexIds[i]);


        Array<T,3> const& vertex = vertices[i];

        // Interpolate rhoBar and PiNeq on the vertex position.
        Array<plint,3> intPos (
                (plint)vertex[0], (plint)vertex[1], (plint)vertex[2] );

        T averageRhoBar1 = 0.0, averageRhoBar2 = 0.0;
        Array<T,SymmetricTensorImpl<T,3>::n> averagePiNeq1, averagePiNeq2;
        averagePiNeq1.resetToZero();
        averagePiNeq2.resetToZero();
        plint n1=0, n2=0;
        T w1=T(), w2=T();

        // x   x . x   x
        for (plint dx=-1; dx<=+2; ++dx) {
            for (plint dy=-1; dy<=+2; ++dy) {
                for (plint dz=-1; dz<=+2; ++dz) {
                    Array<plint,3> pos(intPos+Array<plint,3>(dx,dy,dz));
                    T nextRhoBar = rhoBar->get(pos[0], pos[1], pos[2]);
                    Array<T,SymmetricTensorImpl<T,3>::n>& nextPiNeq = PiNeq->get(pos[0]+ofsPN.x, pos[1]+ofsPN.y, pos[2]+ofsPN.z);
                    Array<T,3> r(pos[0]-vertex[0],pos[1]-vertex[1],pos[2]-vertex[2]);
                    T W = inamuroDeltaFunction<T>().W(r);
                    if (dot(r, normal)>0) {
                        averageRhoBar1 += W * nextRhoBar;
                        averagePiNeq1 += W * nextPiNeq;
                        w1 += W;
                        ++n1;
                    }
                    else {
                        averageRhoBar2 += W * nextRhoBar;
                        averagePiNeq2 += W * nextPiNeq;
                        w2 += W;
                        ++n2;
                    }
                }
            }
        }

        T averageRho1 = rho0;
        T averageRho2 = rho0;
        if (n1>0) {
            averageRhoBar1 /= w1;
            averageRho1 = Descriptor<T>::fullRho(averageRhoBar1);
            averagePiNeq1 /= w1;
        }
        if (n2>0) {
            averageRhoBar2 /= w2;
            averageRho2 = Descriptor<T>::fullRho(averageRhoBar2);
            averagePiNeq2 /= w2;
        }

        // Compute the force on the fluid at the vertex position.
        Array<T,3> averagePi_n1;
        SymmetricTensorImpl<T,3>::matVectMult(averagePiNeq1, normal, averagePi_n1);
        // Important: on the reverse side, take the negative normal.
        Array<T,3> averagePi_n2;
        SymmetricTensorImpl<T,3>::matVectMult(averagePiNeq2, -normal, averagePi_n2);

        // Here we want the force acting on the fluid from the solid. "normal" points towards the fluid,
        // this is why the minus sign in front of the area is needed.
        if (incompressibleModel) {
            g[i] = -areas[i] * (
                        (averageRho2-averageRho1)*Descriptor<T>::cs2*normal +
                        (omega/(T)2.-(T)1.)*(averagePi_n1+averagePi_n2) );  // Incompressible vision
        } else {
            g[i] = -areas[i] * (
                    (averageRho2-averageRho1)*Descriptor<T>::cs2*normal +
                    Descriptor<T>::invRho(averageRhoBar1)*(omega/(T)2.-(T)1.)*averagePi_n1 +  // Compressible vision
                    Descriptor<T>::invRho(averageRhoBar2)*(omega/(T)2.-(T)1.)*averagePi_n2 );
        }
    }
}

template<typename T, template<typename U> class Descriptor, class NormalFunction>
OpenSurfaceImmersedForce3D<T,Descriptor,NormalFunction>* OpenSurfaceImmersedForce3D<T,Descriptor,NormalFunction>::clone() const
{
    return new OpenSurfaceImmersedForce3D<T,Descriptor,NormalFunction>(*this);
}

template<typename T, template<typename U> class Descriptor, class NormalFunction>
void OpenSurfaceImmersedForce3D<T,Descriptor,NormalFunction>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;   // RhoBar
    modified[1] = modif::nothing;   // PiNeq
    modified[2] = modif::nothing;   // Container with triangle data
}

template<typename T, template<typename U> class Descriptor, class NormalFunction>
BlockDomain::DomainT OpenSurfaceImmersedForce3D<T,Descriptor,NormalFunction>::appliesTo() const
{
    return BlockDomain::bulk;
}

/* ******** TwoPhaseInamuroParam3D ************************************ */

template<typename T>
T TwoPhaseInamuroParam3D<T>::area(plint i) const {
    PLB_ASSERT((pluint)i < numVertices);
    return wallData->areas[i];
}

template<typename T>
Array<T,3>& TwoPhaseInamuroParam3D<T>::g(plint i) {
    PLB_ASSERT((pluint)i < numVertices);
    return wallData->g[i];
}

template<typename T>
Array<T,3> TwoPhaseInamuroParam3D<T>::vertex(plint i) const {
    PLB_ASSERT((pluint)i < numVertices);
    return wallData->vertices[i];
}

template<typename T>
Array<T,3> TwoPhaseInamuroParam3D<T>::absoluteVertex(plint i) const {
    PLB_ASSERT((pluint)i < numVertices);
    return wallData->vertices[i] + absOffset;
}

template<typename T>
Array<plint,3> TwoPhaseInamuroParam3D<T>::intVertex(plint i) const {
    PLB_ASSERT((pluint)i < numVertices);
    Array<T,3> vertex = wallData->vertices[i];
    return Array<plint,3>((plint)vertex[0], (plint)vertex[1], (plint)vertex[2]);
}

template<typename T>
T TwoPhaseInamuroParam3D<T>::rhoBar(plint iX, plint iY, plint iZ) const {
    int flag = getFlag(iX,iY,iZ);
    if (flag==freeSurfaceFlag::empty) {
        return rhoBar2_->get(iX+ofsRhoBar2.x,iY+ofsRhoBar2.y,iZ+ofsRhoBar2.z);
    }
    else {
        return rhoBar_->get(iX,iY,iZ);
    }
}

template<typename T>
Array<T,3> TwoPhaseInamuroParam3D<T>::j(plint iX, plint iY, plint iZ) const {
    int flag = getFlag(iX,iY,iZ);
    if (flag==freeSurfaceFlag::empty) {
        return j2_->get(iX+ofsJ2.x,iY+ofsJ2.y,iZ+ofsJ2.z);
    }
    else {
        return j_->get(iX+ofsJ.x,iY+ofsJ.y,iZ+ofsJ.z);
    }
}

template<typename T>
void TwoPhaseInamuroParam3D<T>::addToJ(plint iX, plint iY, plint iZ, Array<T,3> deltaJ)
{
    int flag = getFlag(iX,iY,iZ);
    if (flag==freeSurfaceFlag::interface) {
        j_->get(iX+ofsJ.x,iY+ofsJ.y,iZ+ofsJ.z) += deltaJ;
        j2_->get(iX+ofsJ2.x,iY+ofsJ2.y,iZ+ofsJ2.z) += deltaJ;
    }
    else if (flag==freeSurfaceFlag::empty) {
        j2_->get(iX+ofsJ2.x,iY+ofsJ2.y,iZ+ofsJ2.z) += deltaJ;
    }
    else {
        j_->get(iX+ofsJ.x,iY+ofsJ.y,iZ+ofsJ.z) += deltaJ;
    }
}

template<typename T>
int TwoPhaseInamuroParam3D<T>::getFlag(plint iX, plint iY, plint iZ) const {
    return flag_->get(iX+ofsFlag.x,iY+ofsFlag.y,iZ+ofsFlag.z);
}

template<typename T>
pluint TwoPhaseInamuroParam3D<T>::getGlobalVertexId(plint i) const {
    PLB_ASSERT((pluint)i < numVertices);
    return wallData->globalVertexIds[i];
}

template<typename T>
T TwoPhaseInamuroParam3D<T>::getTau(plint iX, plint iY, plint iZ) const {
    T vf = volumeFraction_->get(iX+ofsVF.x,iY+ofsVF.y,iZ+ofsVF.z);
    return tau*vf + tau2*(1.-vf);
}

template<typename T>
TwoPhaseInamuroParam3D<T>::TwoPhaseInamuroParam3D(std::vector<AtomicBlock3D*>& blocks, T tau_, T tau2_)
    : tau(tau_), tau2(tau2_)
{
    PLB_PRECONDITION( blocks.size()==7 );
    rhoBar_ = dynamic_cast<ScalarField3D<T>*>(blocks[0]);
    rhoBar2_ = dynamic_cast<ScalarField3D<T>*>(blocks[1]);
    j_ = dynamic_cast<TensorField3D<T,3>*>(blocks[2]);
    j2_ = dynamic_cast<TensorField3D<T,3>*>(blocks[3]);
    flag_ = dynamic_cast<ScalarField3D<int>*>(blocks[4]);
    volumeFraction_ = dynamic_cast<ScalarField3D<T>*>(blocks[5]);
    container = dynamic_cast<AtomicContainerBlock3D*>(blocks[6]);

    PLB_ASSERT( rhoBar_ );
    PLB_ASSERT( rhoBar2_ );
    PLB_ASSERT( j_ );
    PLB_ASSERT( j2_ );
    PLB_ASSERT( flag_ );
    PLB_ASSERT( volumeFraction_ );
    PLB_ASSERT( container );

    ofsRhoBar2 = computeRelativeDisplacement(*rhoBar_, *rhoBar2_);
    ofsJ = computeRelativeDisplacement(*rhoBar_, *j_);
    ofsJ2 = computeRelativeDisplacement(*rhoBar_, *j2_);
    ofsFlag = computeRelativeDisplacement(*rhoBar_, *flag_);
    ofsVF = computeRelativeDisplacement(*rhoBar_, *volumeFraction_);

    wallData = dynamic_cast<ImmersedWallData3D<T>*>( container->getData() );
    PLB_ASSERT(wallData);
    absOffset = wallData->offset;

    numVertices = wallData->vertices.size();
    PLB_ASSERT( numVertices == wallData->areas.size() );
    PLB_ASSERT( numVertices == wallData->g.size() );
}

/* ******** TwoPhaseInamuroIteration3D ************************************ */

template<typename T, class VelFunction>
TwoPhaseInamuroIteration3D<T,VelFunction>::TwoPhaseInamuroIteration3D(VelFunction velFunction_, T tau_, T tau2_)
    : velFunction(velFunction_),
      tau(tau_),
      tau2(tau2_)
{ }

template<typename T, class VelFunction>
void TwoPhaseInamuroIteration3D<T,VelFunction>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    TwoPhaseInamuroParam3D<T> param(blocks, tau, tau2);
    std::vector<Array<T,3> > deltaG(param.getNumVertices());

    for (pluint i=0; i<param.getNumVertices(); ++i) {
        Array<T,3> vertex(param.vertex(i));
        Array<plint,3> intPos(param.intVertex(i));

        Array<T,3> averageJ; averageJ.resetToZero();
        T averageRhoBar = T();
        // x   x . x   x
        for (plint dx=-1; dx<=+2; ++dx) {
            for (plint dy=-1; dy<=+2; ++dy) {
                for (plint dz=-1; dz<=+2; ++dz) {
                    Array<plint,3> pos(intPos+Array<plint,3>(dx,dy,dz));
                    T nextRhoBar = param.rhoBar(pos[0], pos[1], pos[2]);
                    Array<T,3> nextJ = param.j(pos[0], pos[1], pos[2]);
                    Array<T,3> r(pos[0]-vertex[0],pos[1]-vertex[1],pos[2]-vertex[2]);
                    T W = inamuroDeltaFunction<T>().W(r);
                    averageJ += W*nextJ;
                    averageRhoBar += W*nextRhoBar;
                }
            }
        }
        //averageJ += (T)0.5*param.g(i);
        Array<T,3> wallVelocity = velFunction(param.absoluteVertex(i));
        deltaG[i] = param.area(i)*((averageRhoBar+(T)1.)*wallVelocity-averageJ);
        //param.g(i) += deltaG[i];
        param.g(i) += deltaG[i]/((T)1.0+averageRhoBar);
    }
    
    for (pluint i=0; i<param.getNumVertices(); ++i) {
        Array<T,3> vertex(param.vertex(i));
        Array<plint,3> intPos(param.intVertex(i));

        for (plint dx=-1; dx<=+2; ++dx) {
            for (plint dy=-1; dy<=+2; ++dy) {
                for (plint dz=-1; dz<=+2; ++dz) {
                    Array<plint,3> pos(intPos+Array<plint,3>(dx,dy,dz));
                    Array<T,3> r(pos[0]-vertex[0],pos[1]-vertex[1],pos[2]-vertex[2]);
                    T W = inamuroDeltaFunction<T>().W(r);
                    param.addToJ(pos[0],pos[1],pos[2], param.getTau(pos[0],pos[1],pos[2])*W*deltaG[i]);
                }
            }
        }
    }
}

template<typename T, class VelFunction>
TwoPhaseInamuroIteration3D<T,VelFunction>* TwoPhaseInamuroIteration3D<T,VelFunction>::clone() const {
    return new TwoPhaseInamuroIteration3D<T,VelFunction>(*this);
}

template<typename T, class VelFunction>
void TwoPhaseInamuroIteration3D<T,VelFunction>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;          // RhoBar
    modified[1] = modif::nothing;          // RhoBar2
    modified[2] = modif::staticVariables;  // j
    modified[3] = modif::staticVariables;  // j2
    modified[4] = modif::nothing;          // flag
    modified[5] = modif::nothing;          // volume fraction
    modified[6] = modif::nothing;          // Container Block with triangle data.
}

template<typename T, class VelFunction>
BlockDomain::DomainT TwoPhaseInamuroIteration3D<T,VelFunction>::appliesTo() const {
    return BlockDomain::bulk;
}

/* ******** TwoPhaseIndexedInamuroIteration3D ************************************ */

template<typename T, class VelFunction>
TwoPhaseIndexedInamuroIteration3D<T,VelFunction>::TwoPhaseIndexedInamuroIteration3D(VelFunction velFunction_, T tau_, T tau2_)
    : velFunction(velFunction_),
      tau(tau_),
      tau2(tau2_)
{ }

template<typename T, class VelFunction>
void TwoPhaseIndexedInamuroIteration3D<T,VelFunction>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    TwoPhaseInamuroParam3D<T> param(blocks, tau, tau2);
    std::vector<Array<T,3> > deltaG(param.getNumVertices());

    for (pluint i=0; i<param.getNumVertices(); ++i) {
        Array<T,3> vertex(param.vertex(i));
        Array<plint,3> intPos(param.intVertex(i));

        Array<T,3> averageJ; averageJ.resetToZero();
        T averageRhoBar = T();
        // x   x . x   x
        for (plint dx=-1; dx<=+2; ++dx) {
            for (plint dy=-1; dy<=+2; ++dy) {
                for (plint dz=-1; dz<=+2; ++dz) {
                    Array<plint,3> pos(intPos+Array<plint,3>(dx,dy,dz));
                    T nextRhoBar = param.rhoBar(pos[0], pos[1], pos[2]);
                    Array<T,3> nextJ = param.j(pos[0], pos[1], pos[2]);
                    Array<T,3> r(pos[0]-vertex[0],pos[1]-vertex[1],pos[2]-vertex[2]);
                    T W = inamuroDeltaFunction<T>().W(r);
                    averageJ += W*nextJ;
                    averageRhoBar += W*nextRhoBar;
                }
            }
        }
        //averageJ += (T)0.5*param.g(i);
        Array<T,3> wallVelocity = velFunction(param.getGlobalVertexId(i));
        deltaG[i] = param.area(i)*((averageRhoBar+(T)1.)*wallVelocity-averageJ);
        //param.g(i) += deltaG[i];
        param.g(i) += deltaG[i]/((T)1.0+averageRhoBar);
    }
    
    for (pluint i=0; i<param.getNumVertices(); ++i) {
        Array<T,3> vertex(param.vertex(i));
        Array<plint,3> intPos(param.intVertex(i));

        for (plint dx=-1; dx<=+2; ++dx) {
            for (plint dy=-1; dy<=+2; ++dy) {
                for (plint dz=-1; dz<=+2; ++dz) {
                    Array<plint,3> pos(intPos+Array<plint,3>(dx,dy,dz));
                    Array<T,3> r(pos[0]-vertex[0],pos[1]-vertex[1],pos[2]-vertex[2]);
                    T W = inamuroDeltaFunction<T>().W(r);
                    param.addToJ(pos[0],pos[1],pos[2], param.getTau(pos[0],pos[1],pos[2])*W*deltaG[i]);
                }
            }
        }
    }
}

template<typename T, class VelFunction>
TwoPhaseIndexedInamuroIteration3D<T,VelFunction>* TwoPhaseIndexedInamuroIteration3D<T,VelFunction>::clone() const {
    return new TwoPhaseIndexedInamuroIteration3D<T,VelFunction>(*this);
}

template<typename T, class VelFunction>
void TwoPhaseIndexedInamuroIteration3D<T,VelFunction>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;          // RhoBar
    modified[1] = modif::nothing;          // RhoBar2
    modified[2] = modif::staticVariables;  // j
    modified[3] = modif::staticVariables;  // j2
    modified[4] = modif::nothing;          // flag
    modified[5] = modif::nothing;          // volume fraction
    modified[6] = modif::nothing;          // Container Block with triangle data.
}

template<typename T, class VelFunction>
BlockDomain::DomainT TwoPhaseIndexedInamuroIteration3D<T,VelFunction>::appliesTo() const {
    return BlockDomain::bulk;
}

/* ******** TwoPhaseConstVelInamuroIteration3D ************************************ */

template<typename T>
TwoPhaseConstVelInamuroIteration3D<T>::TwoPhaseConstVelInamuroIteration3D(Array<T,3> const& wallVelocity_, T tau_, T tau2_)
    : wallVelocity(wallVelocity_),
      tau(tau_),
      tau2(tau2_)
{ }

template<typename T>
void TwoPhaseConstVelInamuroIteration3D<T>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    TwoPhaseInamuroParam3D<T> param(blocks, tau, tau2);
    std::vector<Array<T,3> > deltaG(param.getNumVertices());

    for (pluint i=0; i<param.getNumVertices(); ++i) {
        Array<T,3> vertex(param.vertex(i));
        Array<plint,3> intPos(param.intVertex(i));

        Array<T,3> averageJ; averageJ.resetToZero();
        T averageRhoBar = T();
        // x   x . x   x
        for (plint dx=-1; dx<=+2; ++dx) {
            for (plint dy=-1; dy<=+2; ++dy) {
                for (plint dz=-1; dz<=+2; ++dz) {
                    Array<plint,3> pos(intPos+Array<plint,3>(dx,dy,dz));
                    T nextRhoBar = param.rhoBar(pos[0], pos[1], pos[2]);
                    Array<T,3> nextJ = param.j(pos[0], pos[1], pos[2]);
                    Array<T,3> r(pos[0]-vertex[0],pos[1]-vertex[1],pos[2]-vertex[2]);
                    T W = inamuroDeltaFunction<T>().W(r);
                    averageJ += W*nextJ;
                    averageRhoBar += W*nextRhoBar;
                }
            }
        }
        //averageJ += (T)0.5*param.g(i);
        deltaG[i] = param.area(i)*((averageRhoBar+(T)1.)*wallVelocity-averageJ);
        //param.g(i) += deltaG[i];
        param.g(i) += deltaG[i]/((T)1.0+averageRhoBar);
    }
    
    for (pluint i=0; i<param.getNumVertices(); ++i) {
        Array<T,3> vertex(param.vertex(i));
        Array<plint,3> intPos(param.intVertex(i));

        for (plint dx=-1; dx<=+2; ++dx) {
            for (plint dy=-1; dy<=+2; ++dy) {
                for (plint dz=-1; dz<=+2; ++dz) {
                    Array<plint,3> pos(intPos+Array<plint,3>(dx,dy,dz));
                    Array<T,3> r(pos[0]-vertex[0],pos[1]-vertex[1],pos[2]-vertex[2]);
                    T W = inamuroDeltaFunction<T>().W(r);
                    param.addToJ(pos[0],pos[1],pos[2], param.getTau(pos[0],pos[1],pos[2])*W*deltaG[i]);
                }
            }
        }
    }
}

template<typename T>
TwoPhaseConstVelInamuroIteration3D<T>* TwoPhaseConstVelInamuroIteration3D<T>::clone() const {
    return new TwoPhaseConstVelInamuroIteration3D<T>(*this);
}

template<typename T>
void TwoPhaseConstVelInamuroIteration3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;          // RhoBar
    modified[1] = modif::nothing;          // RhoBar2
    modified[2] = modif::staticVariables;  // j
    modified[3] = modif::staticVariables;  // j2
    modified[4] = modif::nothing;          // flag
    modified[5] = modif::nothing;          // Volume fraction
    modified[6] = modif::nothing;          // Container Block with triangle data.
}

template<typename T>
BlockDomain::DomainT TwoPhaseConstVelInamuroIteration3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}

}  // namespace plb

#endif  // IMMERSED_WALLS_3D_HH

