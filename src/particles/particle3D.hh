/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2013 FlowKit Sarl
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

#ifndef PARTICLE_3D_HH
#define PARTICLE_3D_HH

#include "core/globalDefs.h"
#include "finiteDifference/interpolations3D.h"
#include "particles/particle3D.h"
#include "particles/particleIdentifiers3D.h"
#include <cmath>

namespace plb {

/* *************** class Particle3D ***************************************** */

template<typename T, template<typename U> class Descriptor>
Particle3D<T,Descriptor>::Particle3D()
    : tag(0),
      position(T(),T(),T())
{ }

template<typename T, template<typename U> class Descriptor>
Particle3D<T,Descriptor>::Particle3D(plint tag_, Array<T,3> const& position_)
    : tag(tag_),
      position(position_)
{ }

template<typename T, template<typename U> class Descriptor>
void Particle3D<T,Descriptor>::reset(Array<T,3> const& position_)
{
    position = position_;
}


template<typename T, template<typename U> class Descriptor>
void Particle3D<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    serializer.addValue(tag);
    serializer.addValues<T,3>(position);
}

template<typename T, template<typename U> class Descriptor>
void Particle3D<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
{
    unserializer.readValue(tag);
    unserializer.readValues<T,3>(position);
}

template<typename T, template<typename U> class Descriptor>
plint Particle3D<T,Descriptor>::getTag() const {
    return tag;
}

template<typename T, template<typename U> class Descriptor>
void Particle3D<T,Descriptor>::setTag(plint tag_) {
    tag = tag_;
}


template<typename T, template<typename U> class Descriptor>
bool Particle3D<T,Descriptor>::getVector(plint whichVector, Array<T,3>& vector) const {
    return false;
}

template<typename T, template<typename U> class Descriptor>
bool Particle3D<T,Descriptor>::getScalar(plint whichScalar, T& scalar) const {
    return false;
}

template<typename T, template<typename U> class Descriptor>
bool Particle3D<T,Descriptor>::getTensor(plint whichVector, Array<T,SymmetricTensorImpl<T,3>::n>& tensor) const
{
    return false;
}

template<typename T, template<typename U> class Descriptor>
bool Particle3D<T,Descriptor>::setScalars(std::vector<T> const& scalars)
{
    return false;
}

template<typename T, template<typename U> class Descriptor>
bool Particle3D<T,Descriptor>::setVectors(std::vector<Array<T,3> > const& vectors)
{
    return false;
}

template<typename T, template<typename U> class Descriptor>
bool Particle3D<T,Descriptor>::setTensors (
        std::vector<Array<T,SymmetricTensorImpl<T,3>::n> > const& tensors )
{
    return false;
}

template<typename T, template<typename U> class Descriptor>
void Particle3D<T,Descriptor>::rescale(int dxScale, int dtScale) {
    int dimDx = 1;
    int dimDt = 0;
    T scaleFactor = scaleFromReference(dxScale, dimDx, dtScale, dimDt);
    position *= scaleFactor;
}

/* *************** class PointParticle3D ************************************ */

template<typename T, template<typename U> class Descriptor>
int PointParticle3D<T,Descriptor>::id =
        meta::registerPointParticle3D<T,Descriptor,PointParticle3D<T,Descriptor> >("Point");

template<typename T, template<typename U> class Descriptor>
PointParticle3D<T,Descriptor>::PointParticle3D()
    : Particle3D<T,Descriptor>(),
      velocity(T(),T(),T())
{ }

template<typename T, template<typename U> class Descriptor>
PointParticle3D<T,Descriptor>::PointParticle3D(plint tag_, Array<T,3> const& position_, Array<T,3> const& velocity_)
    : Particle3D<T,Descriptor>(tag_, position_),
      velocity(velocity_)
{ }

template<typename T, template<typename U> class Descriptor>
int PointParticle3D<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void PointParticle3D<T,Descriptor>::reset(Array<T,3> const& position_)
{
    Particle3D<T,Descriptor>::reset(position_);
    velocity.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
void PointParticle3D<T,Descriptor>::velocityToParticle(TensorField3D<T,3>& velocityField, T scaling)
{
    velocity = predictorCorrectorTensorField<T,3>(velocityField, this->getPosition(), scaling);
}

template<typename T, template<typename U> class Descriptor>
void PointParticle3D<T,Descriptor>::rhoBarJtoParticle (
        NTensorField3D<T>& rhoBarJfield, bool velIsJ, T scaling )
{
    T rhoBar;
    Array<T,3> j;
    predictorCorrectorRhoBarJ(rhoBarJfield, this->getPosition(), velIsJ, j, rhoBar);
    if (velIsJ) {
        velocity = j*scaling;
    }
    else {
        velocity = j*scaling*Descriptor<T>::invRho(rhoBar);
    }
}

template<typename T, template<typename U> class Descriptor>
void PointParticle3D<T,Descriptor>::fluidToParticle(BlockLattice3D<T,Descriptor>& fluid, T scaling)
{
    Array<T,3> position1(this->getPosition());
    std::vector<Dot3D> pos(8);
    std::vector<T> weights(8);
    linearInterpolationCoefficients(fluid, position1, pos, weights);
    Array<T,3> tmpVel;
    Array<T,3> velocity1;
    velocity1.resetToZero();
    for (plint iCell=0; iCell<8; ++iCell) {
        // TODO The following condition should never occur, but it is observed
        //      to occur sometimes. In particular, it lead to an assertion being
        //      raised in the two-viscosity code. It is likely that this happens
        //      when the particle is on the edge between two cells, and due to
        //      roundoff errors attributed to the wrong cell. This should be fixed,
        //      and the following condition removed.
        if (contained(pos[iCell].x,pos[iCell].y,pos[iCell].z, fluid.getBoundingBox())) {
            fluid.get(pos[iCell].x,pos[iCell].y,pos[iCell].z).computeVelocity(tmpVel);
        }
        else {
            tmpVel.resetToZero();
        }
        velocity1 += weights[iCell]*tmpVel*scaling;
    }

    Array<T,3> position2(position1+velocity1);
    linearInterpolationCoefficients(fluid, position2, pos, weights);
    Array<T,3> velocity2;
    velocity2.resetToZero();
    for (plint iCell=0; iCell<8; ++iCell) {
        // TODO same as above.
        if (contained(pos[iCell].x,pos[iCell].y,pos[iCell].z, fluid.getBoundingBox())) {
            fluid.get(pos[iCell].x,pos[iCell].y,pos[iCell].z).computeVelocity(tmpVel);
        }
        else {
            tmpVel.resetToZero();
        }
        velocity2 += weights[iCell]*tmpVel*scaling;
    }

    velocity = (velocity1+velocity2)/(T)2;
}

template<typename T, template<typename U> class Descriptor>
void PointParticle3D<T,Descriptor>::advance() {
    PLB_ASSERT( norm(velocity)<1. );
    this->getPosition() += velocity;
}

template<typename T, template<typename U> class Descriptor>
void PointParticle3D<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    Particle3D<T,Descriptor>::serialize(serializer);
    serializer.addValues<T,3>(velocity);
}

template<typename T, template<typename U> class Descriptor>
void PointParticle3D<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
{
    Particle3D<T,Descriptor>::unserialize(unserializer);
    unserializer.readValues<T,3>(velocity);
}

template<typename T, template<typename U> class Descriptor>
PointParticle3D<T,Descriptor>* PointParticle3D<T,Descriptor>::clone() const {
    return new PointParticle3D<T,Descriptor>(*this);
}


template<typename T, template<typename U> class Descriptor>
bool PointParticle3D<T,Descriptor>::getVector(plint whichVector, Array<T,3>& vector) const {
    if (whichVector==0) {
        vector = velocity;
        return true;
    }
    return Particle3D<T,Descriptor>::getVector(whichVector, vector);
}


template<typename T, template<typename U> class Descriptor>
bool PointParticle3D<T,Descriptor>::setVectors(std::vector<Array<T,3> > const& vectors)
{
    if (vectors.size()==1) {
        velocity = vectors[0];
        return true;
    }
    return false;
}

template<typename T, template<typename U> class Descriptor>
void PointParticle3D<T,Descriptor>::rescale(int dxScale, int dtScale) {
    Particle3D<T,Descriptor>::rescale(dxScale, dtScale);
    T dimDx = 1.;
    T dimDt = -1.;
    T scaleFactor = scaleFromReference(dxScale, dimDx, dtScale, dimDt);
    velocity *= scaleFactor;
}


/* *************** class NormedVelocityParticle3D ************************************ */

template<typename T, template<typename U> class Descriptor>
int NormedVelocityParticle3D<T,Descriptor>::id =
        meta::registerGenericParticle3D<T,Descriptor,NormedVelocityParticle3D<T,Descriptor> >("NormedVelocity");

template<typename T, template<typename U> class Descriptor>
NormedVelocityParticle3D<T,Descriptor>::NormedVelocityParticle3D()
    : PointParticle3D<T,Descriptor>(),
      fluidUmax(0.),
      particleUmax(0.),
      exponent(0.)
{ }

template<typename T, template<typename U> class Descriptor>
NormedVelocityParticle3D<T,Descriptor>::NormedVelocityParticle3D (
        plint tag_, Array<T,3> const& position_, Array<T,3> const& velocity_,
        T fluidUmax_, T particleUmax_, T exponent_ )
    : PointParticle3D<T,Descriptor>(tag_, position_, velocity_),
      fluidUmax(fluidUmax_), particleUmax(particleUmax_), exponent(exponent_)
{ }

template<typename T, template<typename U> class Descriptor>
void NormedVelocityParticle3D<T,Descriptor>::advance() {
    T fluidNormu = norm(this->getVelocity());
    PLB_ASSERT( fluidNormu<1. );
    T particleNormU;
    if (fluidNormu >= fluidUmax) {
        particleNormU = particleUmax;
    }
    else {
        particleNormU = pow(fluidNormu/fluidUmax, exponent)*particleUmax;
    }
    PLB_ASSERT( particleNormU<1. );
    this->getPosition() += this->getVelocity()/fluidNormu*particleNormU;
}

template<typename T, template<typename U> class Descriptor>
void NormedVelocityParticle3D<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    PointParticle3D<T,Descriptor>::serialize(serializer);
    serializer.addValue(fluidUmax);
    serializer.addValue(particleUmax);
    serializer.addValue(exponent);
}

template<typename T, template<typename U> class Descriptor>
void NormedVelocityParticle3D<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
{
    PointParticle3D<T,Descriptor>::unserialize(unserializer);
    unserializer.readValue(fluidUmax);
    unserializer.readValue(particleUmax);
    unserializer.readValue(exponent);
}

template<typename T, template<typename U> class Descriptor>
NormedVelocityParticle3D<T,Descriptor>* NormedVelocityParticle3D<T,Descriptor>::clone() const {
    return new NormedVelocityParticle3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void NormedVelocityParticle3D<T,Descriptor>::rescale(int dxScale, int dtScale) {
    PointParticle3D<T,Descriptor>::rescale(dxScale, dtScale);
    T dimDx = 1.;
    T dimDt = -1.;
    T scaleFactor = scaleFromReference(dxScale, dimDx, dtScale, dimDt);
    fluidUmax *= scaleFactor;
    particleUmax *= scaleFactor;
}

template<typename T, template<typename U> class Descriptor>
int NormedVelocityParticle3D<T,Descriptor>::getId() const {
    return id;
}




template<typename T, template<typename U> class Descriptor>
void serialize(Particle3D<T,Descriptor> const& particle, std::vector<char>& data)
{
    HierarchicSerializer serializer(data, particle.getId());
    particle.serialize(serializer);
}

template<typename T, template<typename U> class Descriptor>
void generateAndUnserializeParticles (
        std::vector<char> const& data,
        std::vector<Particle3D<T,Descriptor>*>& particles )
{
    pluint serializerPos = 0;
    while (serializerPos < data.size()) {
        HierarchicUnserializer unserializer(data, serializerPos);
        particles.push_back (
                meta::particleRegistration3D<T,Descriptor>().generate(unserializer) );
        serializerPos = unserializer.getCurrentPos();
    }
}



/* *************** class RestParticle3D ************************************ */

template<typename T, template<typename U> class Descriptor>
int RestParticle3D<T,Descriptor>::id =
        meta::registerGenericParticle3D<T,Descriptor,RestParticle3D<T,Descriptor> >("Rest");

template<typename T, template<typename U> class Descriptor>
RestParticle3D<T,Descriptor>::RestParticle3D()
{ }

template<typename T, template<typename U> class Descriptor>
RestParticle3D<T,Descriptor>::RestParticle3D (
        plint tag_, Array<T,3> const& position )
    : Particle3D<T,Descriptor>(tag_, position)
{ }

template<typename T, template<typename U> class Descriptor>
void RestParticle3D<T,Descriptor>::velocityToParticle(TensorField3D<T,3>& velocityField, T scaling) { }

template<typename T, template<typename U> class Descriptor>
void RestParticle3D<T,Descriptor>::rhoBarJtoParticle (
        NTensorField3D<T>& rhoBarJfield, bool velIsJ, T scaling )
{ }

template<typename T, template<typename U> class Descriptor>
void RestParticle3D<T,Descriptor>::fluidToParticle(BlockLattice3D<T,Descriptor>& fluid, T scaling) { }

template<typename T, template<typename U> class Descriptor>
void RestParticle3D<T,Descriptor>::advance() { }

template<typename T, template<typename U> class Descriptor>
int RestParticle3D<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
RestParticle3D<T,Descriptor>* RestParticle3D<T,Descriptor>::clone() const {
    return new RestParticle3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
bool RestParticle3D<T,Descriptor>::getVector(plint whichVector, Array<T,3>& vector) const {
    return Particle3D<T,Descriptor>::getVector(whichVector, vector);
}


/* *************** class VerletParticle3D ************************************ */

template<typename T, template<typename U> class Descriptor>
int VerletParticle3D<T,Descriptor>::id =
        meta::registerGenericParticle3D<T,Descriptor,VerletParticle3D<T,Descriptor> >("Verlet");

template<typename T, template<typename U> class Descriptor>
VerletParticle3D<T,Descriptor>::VerletParticle3D()
    : v(T(),T(),T()),
      vHalfTime(T(),T(),T()),
      a(T(),T(),T()),
      fluidCompliance((T)1.),
      rho((T)1.),
      invRho((T)1.)
{ }

template<typename T, template<typename U> class Descriptor>
VerletParticle3D<T,Descriptor>::VerletParticle3D (
        plint tag_, Array<T,3> const& position )
    : Particle3D<T,Descriptor>(tag_, position),
      v(T(),T(),T()),
      vHalfTime(T(),T(),T()),
      a(T(),T(),T()),
      fluidCompliance((T)1.),
      rho((T)1.),
      invRho((T)1.)
{ }

template<typename T, template<typename U> class Descriptor>
void VerletParticle3D<T,Descriptor>::velocityToParticle(TensorField3D<T,3>& velocityField, T scaling)
{
    Array<T,3> position(this->getPosition());
    std::vector<Dot3D> pos(8);
    std::vector<T> weights(8);
    linearInterpolationCoefficients(velocityField, position, pos, weights);
    Array<T,3> fluidVelocity;
    fluidVelocity.resetToZero();
    for (plint iCell=0; iCell<8; ++iCell) {
        fluidVelocity += weights[iCell]*velocityField.get(pos[iCell].x,pos[iCell].y,pos[iCell].z)*scaling;
    }

    Array<T,3> force( (fluidVelocity-this->get_v()) *fluidCompliance);
    this->set_a(force / this->get_rho());
}

template<typename T, template<typename U> class Descriptor>
void VerletParticle3D<T,Descriptor>::rhoBarJtoParticle (
        NTensorField3D<T>& rhoBarJfield, bool velIsJ, T scaling )
{
    Array<T,3> position(this->getPosition());
    std::vector<Dot3D> pos(8);
    std::vector<T> weights(8);
    linearInterpolationCoefficients(rhoBarJfield, position, pos, weights);
    Array<T,3> j;
    j.resetToZero();
    T rhoBar = T();
    for (plint iCell=0; iCell<8; ++iCell) {
        T* data = rhoBarJfield.get(pos[iCell].x,pos[iCell].y,pos[iCell].z);
        j.add_from_cArray(data+1, weights[iCell]);
        rhoBar += weights[iCell]*(*data);
    }
    Array<T,3> fluidVelocity;
    if (velIsJ) {
        fluidVelocity = j*scaling;
    }
    else {
        fluidVelocity = j*scaling*Descriptor<T>::invRho(rhoBar);
    }

    Array<T,3> force( (fluidVelocity-this->get_v()) *fluidCompliance);
    this->set_a(force / this->get_rho());
}

template<typename T, template<typename U> class Descriptor>
void VerletParticle3D<T,Descriptor>::fluidToParticle(BlockLattice3D<T,Descriptor>& fluid, T scaling)
{
    Array<T,3> position(this->getPosition());
    std::vector<Dot3D> pos(8);
    std::vector<T> weights(8);
    linearInterpolationCoefficients(fluid, position, pos, weights);
    Array<T,3> tmpVel;
    Array<T,3> fluidVelocity;
    fluidVelocity.resetToZero();
    for (plint iCell=0; iCell<8; ++iCell) {
        if (contained(pos[iCell].x,pos[iCell].y,pos[iCell].z, fluid.getBoundingBox())) {
            fluid.get(pos[iCell].x,pos[iCell].y,pos[iCell].z).computeVelocity(tmpVel);
        }
        else {
            tmpVel.resetToZero();
        }
        fluidVelocity += weights[iCell]*tmpVel*scaling;
    }

    Array<T,3> force( (fluidVelocity-this->get_v()) *fluidCompliance);
    this->set_a(force / this->get_rho());
}

template<typename T, template<typename U> class Descriptor>
void VerletParticle3D<T,Descriptor>::advance() {
    vHalfTime = v + (T)0.5*a;
    this->getPosition() += vHalfTime;
}

template<typename T, template<typename U> class Descriptor>
int VerletParticle3D<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void VerletParticle3D<T,Descriptor>::reset(Array<T,3> const& position_)
{
    Particle3D<T,Descriptor>::reset(position_);
    v.resetToZero();
    vHalfTime.resetToZero();
    a.resetToZero();
}

template<typename T, template<typename U> class Descriptor>
void VerletParticle3D<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    Particle3D<T,Descriptor>::serialize(serializer);
    serializer.addValues<T,3>(v);
    serializer.addValues<T,3>(vHalfTime);
    serializer.addValues<T,3>(a);
    serializer.addValue<T>(fluidCompliance);
    serializer.addValue<T>(rho);
    serializer.addValue<T>(invRho);
}

template<typename T, template<typename U> class Descriptor>
void VerletParticle3D<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
{
    Particle3D<T,Descriptor>::unserialize(unserializer);
    unserializer.readValues<T,3>(v);
    unserializer.readValues<T,3>(vHalfTime);
    unserializer.readValues<T,3>(a);
    unserializer.readValue<T>(fluidCompliance);
    unserializer.readValue<T>(rho);
    unserializer.readValue<T>(invRho);
}

template<typename T, template<typename U> class Descriptor>
VerletParticle3D<T,Descriptor>* VerletParticle3D<T,Descriptor>::clone() const {
    return new VerletParticle3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
bool VerletParticle3D<T,Descriptor>::getVector(plint whichVector, Array<T,3>& vector) const {
    if (whichVector==0) {
        vector = get_v();
        return true;
    }
    else if (whichVector==1) {
        vector = get_a();
        return true;
    }
    return Particle3D<T,Descriptor>::getVector(whichVector, vector);
}

template<typename T, template<typename U> class Descriptor>
bool VerletParticle3D<T,Descriptor>::getScalar(plint whichScalar, T& scalar) const {
    if (whichScalar==0) {
        scalar = rho;
        return true;
    }
    else if (whichScalar==1) {
        scalar = fluidCompliance;
        return true;
    }
    return Particle3D<T,Descriptor>::getScalar(whichScalar, scalar);
}

}  // namespace plb

#endif  // PARTICLE_3D_H
