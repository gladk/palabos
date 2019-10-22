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

#ifndef BODY_FORCE_3D_HH
#define BODY_FORCE_3D_HH

#include "multiPhysics/freeSurfaceUtil3D.h"
#include "multiPhysics/bodyForce3D.h"

namespace plb {

/* ****************** AddConstForceToMomentum3D *************************************************** */

template<typename T, template<typename U> class Descriptor>
AddConstForceToMomentum3D<T,Descriptor>::AddConstForceToMomentum3D(Array<T,3> const& force_)
    : force(force_)
{ }

template<typename T, template<typename U> class Descriptor>
void AddConstForceToMomentum3D<T,Descriptor>::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks)
{
    PLB_ASSERT(blocks.size() == 3);
    BlockLattice3D<T,Descriptor>* lattice = dynamic_cast<BlockLattice3D<T,Descriptor>*>(blocks[0]);
    PLB_ASSERT(lattice);
    ScalarField3D<T>* rhoBar = dynamic_cast<ScalarField3D<T>*>(blocks[1]);
    PLB_ASSERT(rhoBar);
    TensorField3D<T,3>* j = dynamic_cast<TensorField3D<T,3>*>(blocks[2]);
    PLB_ASSERT(j);

    Dot3D ofsR = computeRelativeDisplacement(*lattice, *rhoBar);
    Dot3D ofsJ = computeRelativeDisplacement(*lattice, *j);

    for (plint iX = domain.x0; iX <= domain.x1; iX++) {
        for (plint iY = domain.y0; iY <= domain.y1; iY++) {
            for (plint iZ = domain.z0; iZ <= domain.z1; iZ++) {
                // Momentum correction.
                Cell<T,Descriptor>& cell = lattice->get(iX, iY, iZ);
                if (cell.getDynamics().hasMoments()) {
                    T dynamicOmega = cell.getDynamics().getDynamicParameter(dynamicParams::dynamicOmega, cell);
                    T tau = 0.0;
                    if (!util::isZero(dynamicOmega)) {
                        tau = (T) 1 / dynamicOmega;
                    } else {
                        tau = (T) 1 / cell.getDynamics().getOmega();
                    }
                    T rho = Descriptor<T>::fullRho(rhoBar->get(iX+ofsR.x, iY+ofsR.y, iZ+ofsR.z));
                    j->get(iX+ofsJ.x, iY+ofsJ.y, iZ+ofsJ.z) += rho * tau * force;
                }
            }
        }
    }   
}

template<typename T, template<typename U> class Descriptor>
AddConstForceToMomentum3D<T,Descriptor>* AddConstForceToMomentum3D<T,Descriptor>::clone() const
{
    return new AddConstForceToMomentum3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void AddConstForceToMomentum3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;           // Lattice
    modified[1] = modif::nothing;           // rhoBar
    modified[2] = modif::staticVariables;   // j
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT AddConstForceToMomentum3D<T,Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template<typename T, template<typename U> class Descriptor>
void addConstForceToMomentum(MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<T>& rhoBar,
        MultiTensorField3D<T,3>& j, Array<T,3> const& force, Box3D const& domain)
{
    std::vector<MultiBlock3D*> args;
    args.push_back(&lattice);
    args.push_back(&rhoBar);
    args.push_back(&j);
    applyProcessingFunctional(new AddConstForceToMomentum3D<T,Descriptor>(force), domain, args);
}

/* ****************** AddCustomForceToMomentum3D *************************************************** */

template<typename T, template<typename U> class Descriptor, class ForceFunction>
AddCustomForceToMomentum3D<T,Descriptor,ForceFunction>::AddCustomForceToMomentum3D(ForceFunction f_)
    : f(f_)
{ }

template<typename T, template<typename U> class Descriptor, class ForceFunction>
void AddCustomForceToMomentum3D<T,Descriptor,ForceFunction>::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks)
{
    PLB_ASSERT(blocks.size() == 3);
    BlockLattice3D<T,Descriptor>* lattice = dynamic_cast<BlockLattice3D<T,Descriptor>*>(blocks[0]);
    PLB_ASSERT(lattice);
    ScalarField3D<T>* rhoBar = dynamic_cast<ScalarField3D<T>*>(blocks[1]);
    PLB_ASSERT(rhoBar);
    TensorField3D<T,3>* j = dynamic_cast<TensorField3D<T,3>*>(blocks[2]);
    PLB_ASSERT(j);

    Dot3D ofsR = computeRelativeDisplacement(*lattice, *rhoBar);
    Dot3D ofsJ = computeRelativeDisplacement(*lattice, *j);

    Dot3D location = lattice->getLocation();

    for (plint iX = domain.x0; iX <= domain.x1; iX++) {
        plint x = iX + location.x;
        for (plint iY = domain.y0; iY <= domain.y1; iY++) {
            plint y = iY + location.y;
            for (plint iZ = domain.z0; iZ <= domain.z1; iZ++) {
                plint z = iZ + location.z;
                // Momentum correction.
                Cell<T,Descriptor>& cell = lattice->get(iX, iY, iZ);
                if (cell.getDynamics().hasMoments()) {
                    T dynamicOmega = cell.getDynamics().getDynamicParameter(dynamicParams::dynamicOmega, cell);
                    T tau = 0.0;
                    if (!util::isZero(dynamicOmega)) {
                        tau = (T) 1 / dynamicOmega;
                    } else {
                        tau = (T) 1 / cell.getDynamics().getOmega();
                    }
                    T rho = Descriptor<T>::fullRho(rhoBar->get(iX+ofsR.x, iY+ofsR.y, iZ+ofsR.z));
                    Array<T,3> force;
                    f(x, y, z, force);
                    j->get(iX+ofsJ.x, iY+ofsJ.y, iZ+ofsJ.z) += rho * tau * force;
                }
            }
        }
    }   
}

template<typename T, template<typename U> class Descriptor, class ForceFunction>
AddCustomForceToMomentum3D<T,Descriptor,ForceFunction>*
    AddCustomForceToMomentum3D<T,Descriptor,ForceFunction>::clone() const
{
    return new AddCustomForceToMomentum3D<T,Descriptor,ForceFunction>(*this);
}

template<typename T, template<typename U> class Descriptor, class ForceFunction>
void AddCustomForceToMomentum3D<T,Descriptor,ForceFunction>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;           // Lattice
    modified[1] = modif::nothing;           // rhoBar
    modified[2] = modif::staticVariables;   // j
}

template<typename T, template<typename U> class Descriptor, class ForceFunction>
BlockDomain::DomainT AddCustomForceToMomentum3D<T,Descriptor,ForceFunction>::appliesTo() const
{
    return BlockDomain::bulk;
}

template<typename T, template<class U> class Descriptor, class ForceFunction>
void addCustomForceToMomentum(MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<T>& rhoBar,
        MultiTensorField3D<T,3>& j, ForceFunction f, Box3D const& domain)
{
    std::vector<MultiBlock3D*> args;
    args.push_back(&lattice);
    args.push_back(&rhoBar);
    args.push_back(&j);
    applyProcessingFunctional(new AddCustomForceToMomentum3D<T,Descriptor,ForceFunction>(f), domain, args);
}

/* ****************** AddForceToMomentum3D *************************************************** */

template<typename T, template<typename U> class Descriptor>
void AddForceToMomentum3D<T,Descriptor>::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks)
{
    PLB_ASSERT(blocks.size() == 4);
    BlockLattice3D<T,Descriptor>* lattice = dynamic_cast<BlockLattice3D<T,Descriptor>*>(blocks[0]);
    PLB_ASSERT(lattice);
    ScalarField3D<T>* rhoBar = dynamic_cast<ScalarField3D<T>*>(blocks[1]);
    PLB_ASSERT(rhoBar);
    TensorField3D<T,3>* j = dynamic_cast<TensorField3D<T,3>*>(blocks[2]);
    PLB_ASSERT(j);
    TensorField3D<T,3>* force = dynamic_cast<TensorField3D<T,3>*>(blocks[3]);
    PLB_ASSERT(force);

    Dot3D ofsR = computeRelativeDisplacement(*lattice, *rhoBar);
    Dot3D ofsJ = computeRelativeDisplacement(*lattice, *j);
    Dot3D ofsF = computeRelativeDisplacement(*lattice, *force);

    for (plint iX = domain.x0; iX <= domain.x1; iX++) {
        for (plint iY = domain.y0; iY <= domain.y1; iY++) {
            for (plint iZ = domain.z0; iZ <= domain.z1; iZ++) {
                // Momentum correction.
                Cell<T,Descriptor>& cell = lattice->get(iX, iY, iZ);
                if (cell.getDynamics().hasMoments()) {
                    T dynamicOmega = cell.getDynamics().getDynamicParameter(dynamicParams::dynamicOmega, cell);
                    T tau = 0.0;
                    if (!util::isZero(dynamicOmega)) {
                        tau = (T) 1 / dynamicOmega;
                    } else {
                        tau = (T) 1 / cell.getDynamics().getOmega();
                    }
                    T rho = Descriptor<T>::fullRho(rhoBar->get(iX+ofsR.x, iY+ofsR.y, iZ+ofsR.z));
                    Array<T,3> const& f = force->get(iX+ofsF.x, iY+ofsF.y, iZ+ofsF.z);
                    j->get(iX+ofsJ.x, iY+ofsJ.y, iZ+ofsJ.z) += rho * tau * f;
                }
            }
        }
    }   
}

template<typename T, template<typename U> class Descriptor>
AddForceToMomentum3D<T,Descriptor>* AddForceToMomentum3D<T,Descriptor>::clone() const
{
    return new AddForceToMomentum3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void AddForceToMomentum3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;           // Lattice
    modified[1] = modif::nothing;           // rhoBar
    modified[2] = modif::staticVariables;   // j
    modified[3] = modif::nothing;           // Force
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT AddForceToMomentum3D<T,Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template<typename T, template<typename U> class Descriptor>
void addForceToMomentum(MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<T>& rhoBar,
        MultiTensorField3D<T,3>& j, MultiTensorField3D<T,3>& force, Box3D const& domain)
{
    std::vector<MultiBlock3D*> args;
    args.push_back(&lattice);
    args.push_back(&rhoBar);
    args.push_back(&j);
    args.push_back(&force);
    applyProcessingFunctional(new AddForceToMomentum3D<T,Descriptor>, domain, args);
}

/* ****************** FreeSurfaceAddForceToMomentum3D *************************************************** */

template<typename T, template<typename U> class Descriptor>
void FreeSurfaceAddForceToMomentum3D<T,Descriptor>::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks)
{
    PLB_ASSERT(blocks.size() == 5);
    BlockLattice3D<T,Descriptor>* lattice = dynamic_cast<BlockLattice3D<T,Descriptor>*>(blocks[0]);
    PLB_ASSERT(lattice);
    ScalarField3D<T>* rhoBar = dynamic_cast<ScalarField3D<T>*>(blocks[1]);
    PLB_ASSERT(rhoBar);
    TensorField3D<T,3>* j = dynamic_cast<TensorField3D<T,3>*>(blocks[2]);
    PLB_ASSERT(j);
    ScalarField3D<int>* flag = dynamic_cast<ScalarField3D<int>*>(blocks[3]);
    PLB_ASSERT(flag);
    TensorField3D<T,3>* force = dynamic_cast<TensorField3D<T,3>*>(blocks[4]);
    PLB_ASSERT(force);

    Dot3D ofsRhoBar = computeRelativeDisplacement(*lattice, *rhoBar);
    Dot3D ofsJ = computeRelativeDisplacement(*lattice, *j);
    Dot3D ofsFlag = computeRelativeDisplacement(*lattice, *flag);
    Dot3D ofsForce = computeRelativeDisplacement(*lattice, *force);

    for (plint iX = domain.x0; iX <= domain.x1; iX++) {
        for (plint iY = domain.y0; iY <= domain.y1; iY++) {
            for (plint iZ = domain.z0; iZ <= domain.z1; iZ++) {
                if (freeSurfaceFlag::isWet(flag->get(iX+ofsFlag.x, iY+ofsFlag.y, iZ+ofsFlag.z))) {
                    // Momentum correction.
                    Cell<T,Descriptor>& cell = lattice->get(iX, iY, iZ);
                    T dynamicOmega = cell.getDynamics().getDynamicParameter(dynamicParams::dynamicOmega, cell);
                    T tau = 0.0;
                    if (!util::isZero(dynamicOmega)) {
                        tau = (T) 1 / dynamicOmega;
                    } else {
                        tau = (T) 1 / cell.getDynamics().getOmega();
                    }
                    T rho = Descriptor<T>::fullRho(rhoBar->get(iX+ofsRhoBar.x, iY+ofsRhoBar.y, iZ+ofsRhoBar.z));
                    Array<T,3> const& f = force->get(iX+ofsForce.x, iY+ofsForce.y, iZ+ofsForce.z);
                    j->get(iX+ofsJ.x, iY+ofsJ.y, iZ+ofsJ.z) += rho * tau * f;
                }
            }
        }
    }   
}

template<typename T, template<typename U> class Descriptor>
FreeSurfaceAddForceToMomentum3D<T,Descriptor>* FreeSurfaceAddForceToMomentum3D<T,Descriptor>::clone() const
{
    return new FreeSurfaceAddForceToMomentum3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void FreeSurfaceAddForceToMomentum3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;           // Lattice
    modified[1] = modif::nothing;           // rhoBar
    modified[2] = modif::staticVariables;   // j
    modified[3] = modif::nothing;           // Flag
    modified[4] = modif::nothing;           // Force
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT FreeSurfaceAddForceToMomentum3D<T,Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template<typename T, template<typename U> class Descriptor>
void freeSurfaceAddForceToMomentum(MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<T>& rhoBar,
        MultiTensorField3D<T,3>& j, MultiScalarField3D<int>& flag, MultiTensorField3D<T,3>& force,
        Box3D const& domain)
{
    std::vector<MultiBlock3D*> args;
    args.push_back(&lattice);
    args.push_back(&rhoBar);
    args.push_back(&j);
    args.push_back(&flag);
    args.push_back(&force);
    applyProcessingFunctional(new FreeSurfaceAddForceToMomentum3D<T,Descriptor>, domain, args);
}

/* ****************** ComputeRotatingFrameForce3D *************************************************** */

template<typename T, template<typename U> class Descriptor>
ComputeRotatingFrameForce3D<T,Descriptor>::ComputeRotatingFrameForce3D(
        Array<T,3> const& constantForce_, Array<T,3> const& angularVelocity_,
        Array<T,3> const& origin_, bool incompressibleModel_)
    : constantForce(constantForce_),
      angularVelocity(angularVelocity_),
      origin(origin_),
      incompressibleModel(incompressibleModel_)
{ }

template<typename T, template<typename U> class Descriptor>
void ComputeRotatingFrameForce3D<T,Descriptor>::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks)
{
    PLB_ASSERT(blocks.size() == 4);
    BlockLattice3D<T,Descriptor>* lattice = dynamic_cast<BlockLattice3D<T,Descriptor>*>(blocks[0]);
    PLB_ASSERT(lattice);
    ScalarField3D<T>* rhoBar = dynamic_cast<ScalarField3D<T>*>(blocks[1]);
    PLB_ASSERT(rhoBar);
    TensorField3D<T,3>* j = dynamic_cast<TensorField3D<T,3>*>(blocks[2]);
    PLB_ASSERT(j);
    TensorField3D<T,3>* force = dynamic_cast<TensorField3D<T,3>*>(blocks[3]);
    PLB_ASSERT(force);

    Dot3D loc = lattice->getLocation();

    Dot3D ofsRhoBar = computeRelativeDisplacement(*lattice, *rhoBar);
    Dot3D ofsJ = computeRelativeDisplacement(*lattice, *j);
    Dot3D ofsForce = computeRelativeDisplacement(*lattice, *force);

    if (!util::isZero(norm(constantForce))) {
        T angularVelocityNorm = norm(angularVelocity);
        if (!util::isZero(angularVelocityNorm)) {
            Array<T,3> normedAxis = angularVelocity / angularVelocityNorm;
            plint t = lattice->getTimeCounter().getTime();
            T theta = -angularVelocityNorm * t;
            constantForce = rotateAtOrigin(constantForce, normedAxis, theta);
        }
    }

    for (plint iX = domain.x0; iX <= domain.x1; iX++) {
        T x = iX + loc.x - origin[0];
        for (plint iY = domain.y0; iY <= domain.y1; iY++) {
            T y = iY + loc.y - origin[1];
            for (plint iZ = domain.z0; iZ <= domain.z1; iZ++) {
                T z = iZ + loc.z - origin[2];

                // Constant force.
                Array<T,3> newForce(constantForce);

                // Coriolis and centripetal forces.
                Array<T,3> r(x, y, z);
                Array<T,3> velocity(j->get(iX+ofsJ.x, iY+ofsJ.y, iZ+ofsJ.z));
                if (!incompressibleModel) {
                    T rho = Descriptor<T>::fullRho(rhoBar->get(iX+ofsRhoBar.x, iY+ofsRhoBar.y, iZ+ofsRhoBar.z));
                    velocity /= rho;
                }
                Array<T,3> const& oldForce = force->get(iX+ofsForce.x, iY+ofsForce.y, iZ+ofsForce.z);
                velocity += (T) 0.5 * oldForce;
                newForce += -((T) 2 * crossProduct(angularVelocity, velocity) +
                                      crossProduct(angularVelocity, crossProduct(angularVelocity, r)));

                force->get(iX+ofsForce.x, iY+ofsForce.y, iZ+ofsForce.z) = newForce;
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
ComputeRotatingFrameForce3D<T,Descriptor>* ComputeRotatingFrameForce3D<T,Descriptor>::clone() const
{
    return new ComputeRotatingFrameForce3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT ComputeRotatingFrameForce3D<T,Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template<typename T, template<typename U> class Descriptor>
void ComputeRotatingFrameForce3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;           // Lattice
    modified[1] = modif::nothing;           // rhoBar
    modified[2] = modif::nothing;           // j
    modified[3] = modif::staticVariables;   // Force
}

template<typename T, template<typename U> class Descriptor>
void computeRotatingFrameForce(MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<T>& rhoBar,
        MultiTensorField3D<T,3>& j, MultiTensorField3D<T,3>& force, Array<T,3> const& constantForce,
        Array<T,3> const& angularVelocity, Array<T,3> const& origin, bool incompressibleModel, Box3D domain)
{
    std::vector<MultiBlock3D*> args;
    args.push_back(&lattice);
    args.push_back(&rhoBar);
    args.push_back(&j);
    args.push_back(&force);
    applyProcessingFunctional(
            new ComputeRotatingFrameForce3D<T,Descriptor>(constantForce, angularVelocity, origin, incompressibleModel),
            domain, args);
}

/* ****************** FreeSurfaceComputeRotatingFrameForce3D *************************************************** */

template<typename T, template<typename U> class Descriptor>
FreeSurfaceComputeRotatingFrameForce3D<T,Descriptor>::FreeSurfaceComputeRotatingFrameForce3D(
        Array<T,3> const& constantForce_, Array<T,3> const& angularVelocity_,
        Array<T,3> const& origin_, bool incompressibleModel_)
    : constantForce(constantForce_),
      angularVelocity(angularVelocity_),
      origin(origin_),
      incompressibleModel(incompressibleModel_)
{ }

template<typename T, template<typename U> class Descriptor>
void FreeSurfaceComputeRotatingFrameForce3D<T,Descriptor>::processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks)
{
    PLB_ASSERT(blocks.size() == 5);
    BlockLattice3D<T,Descriptor>* lattice = dynamic_cast<BlockLattice3D<T,Descriptor>*>(blocks[0]);
    PLB_ASSERT(lattice);
    ScalarField3D<T>* rhoBar = dynamic_cast<ScalarField3D<T>*>(blocks[1]);
    PLB_ASSERT(rhoBar);
    TensorField3D<T,3>* j = dynamic_cast<TensorField3D<T,3>*>(blocks[2]);
    PLB_ASSERT(j);
    ScalarField3D<int>* flag = dynamic_cast<ScalarField3D<int>*>(blocks[3]);
    PLB_ASSERT(flag);
    TensorField3D<T,3>* force = dynamic_cast<TensorField3D<T,3>*>(blocks[4]);
    PLB_ASSERT(force);

    Dot3D loc = lattice->getLocation();

    Dot3D ofsRhoBar = computeRelativeDisplacement(*lattice, *rhoBar);
    Dot3D ofsJ = computeRelativeDisplacement(*lattice, *j);
    Dot3D ofsFlag = computeRelativeDisplacement(*lattice, *flag);
    Dot3D ofsForce = computeRelativeDisplacement(*lattice, *force);

    if (!util::isZero(norm(constantForce))) {
        T angularVelocityNorm = norm(angularVelocity);
        if (!util::isZero(angularVelocityNorm)) {
            Array<T,3> normedAxis = angularVelocity / angularVelocityNorm;
            plint t = lattice->getTimeCounter().getTime();
            T theta = -angularVelocityNorm * t;
            constantForce = rotateAtOrigin(constantForce, normedAxis, theta);
        }
    }

    for (plint iX = domain.x0; iX <= domain.x1; iX++) {
        T x = iX + loc.x - origin[0];
        for (plint iY = domain.y0; iY <= domain.y1; iY++) {
            T y = iY + loc.y - origin[1];
            for (plint iZ = domain.z0; iZ <= domain.z1; iZ++) {
                T z = iZ + loc.z - origin[2];

                Array<T,3> newForce(Array<T,3>::zero());

                if (freeSurfaceFlag::isWet(flag->get(iX+ofsFlag.x, iY+ofsFlag.y, iZ+ofsFlag.z))) {

                    // Constant force.
                    newForce = constantForce;

                    // Coriolis and centripetal forces.
                    Array<T,3> r(x, y, z);
                    Array<T,3> velocity(j->get(iX+ofsJ.x, iY+ofsJ.y, iZ+ofsJ.z));
                    if (!incompressibleModel) {
                        T rho = Descriptor<T>::fullRho(rhoBar->get(iX+ofsRhoBar.x, iY+ofsRhoBar.y, iZ+ofsRhoBar.z));
                        velocity /= rho;
                    }
                    Array<T,3> const& oldForce = force->get(iX+ofsForce.x, iY+ofsForce.y, iZ+ofsForce.z);
                    velocity += (T) 0.5 * oldForce;
                    newForce += -((T) 2 * crossProduct(angularVelocity, velocity) +
                                          crossProduct(angularVelocity, crossProduct(angularVelocity, r)));
                }

                force->get(iX+ofsForce.x, iY+ofsForce.y, iZ+ofsForce.z) = newForce;
            }
        }
    }
}

template<typename T, template<typename U> class Descriptor>
FreeSurfaceComputeRotatingFrameForce3D<T,Descriptor>* FreeSurfaceComputeRotatingFrameForce3D<T,Descriptor>::clone() const
{
    return new FreeSurfaceComputeRotatingFrameForce3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
BlockDomain::DomainT FreeSurfaceComputeRotatingFrameForce3D<T,Descriptor>::appliesTo() const
{
    return BlockDomain::bulk;
}

template<typename T, template<typename U> class Descriptor>
void FreeSurfaceComputeRotatingFrameForce3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;           // Lattice
    modified[1] = modif::nothing;           // rhoBar
    modified[2] = modif::nothing;           // j
    modified[3] = modif::nothing;           // Flag
    modified[4] = modif::staticVariables;   // Force
}

template<typename T, template<typename U> class Descriptor>
void freeSurfaceComputeRotatingFrameForce(MultiBlockLattice3D<T,Descriptor>& lattice, MultiScalarField3D<T>& rhoBar,
        MultiTensorField3D<T,3>& j, MultiScalarField3D<int>& flag, MultiTensorField3D<T,3>& force,
        Array<T,3> const& constantForce, Array<T,3> const& angularVelocity, Array<T,3> const& origin,
        bool incompressibleModel, Box3D domain)
{
    std::vector<MultiBlock3D*> args;
    args.push_back(&lattice);
    args.push_back(&rhoBar);
    args.push_back(&j);
    args.push_back(&flag);
    args.push_back(&force);
    applyProcessingFunctional(
            new FreeSurfaceComputeRotatingFrameForce3D<T,Descriptor>(constantForce, angularVelocity, origin, incompressibleModel),
            domain, args);
}

}  // namespace plb

#endif  // BODY_FORCE_3D_HH
