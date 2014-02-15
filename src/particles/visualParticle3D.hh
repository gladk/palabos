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

#ifndef VISUAL_PARTICLE_3D_HH
#define VISUAL_PARTICLE_3D_HH

#include "core/globalDefs.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "particles/visualParticle3D.h"

#include <limits>
#include <cmath>

namespace plb {

/* *************** class VisualParticle3D ************************************ */

template<typename T, template<typename U> class Descriptor>
int VisualParticle3D<T,Descriptor>::id =
        meta::registerVisualParticle3D<T,Descriptor,VisualParticle3D<T,Descriptor> >("Visual");

template<typename T, template<typename U> class Descriptor>
VisualParticle3D<T,Descriptor>::VisualParticle3D()
{ }

template<typename T, template<typename U> class Descriptor>
VisualParticle3D<T,Descriptor>::VisualParticle3D (
        plint tag_, Array<T,3> const& position )
    : Particle3D<T,Descriptor>(tag_, position)
{ }

template<typename T, template<typename U> class Descriptor>
VisualParticle3D<T,Descriptor>::VisualParticle3D (
        plint tag_, Array<T,3> const& position,
        std::vector<T> const& scalars_,
        std::vector<Array<T,3> > const& vectors_ )
    : Particle3D<T,Descriptor>(tag_, position),
      scalars(scalars_),
      vectors(vectors_)
{ }

template<typename T, template<typename U> class Descriptor>
int VisualParticle3D<T,Descriptor>::getId() const {
    return id;
}

template<typename T, template<typename U> class Descriptor>
void VisualParticle3D<T,Descriptor>::reset(Array<T,3> const& position_)
{
    Particle3D<T,Descriptor>::reset(position_);
    std::fill(scalars.begin(), scalars.end(), T());
    std::fill(vectors.begin(), vectors.end(), Array<T,3>(T(),T(),T()));
}

template<typename T, template<typename U> class Descriptor>
void VisualParticle3D<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    Particle3D<T,Descriptor>::serialize(serializer);
    pluint scalarSize = scalars.size();
    serializer.addValue(scalarSize);
    for (pluint iScalar=0; iScalar<scalars.size(); ++iScalar) {
        serializer.addValue(scalars[iScalar]);
    }
    pluint vectorSize = vectors.size();
    serializer.addValue(vectorSize);
    for (pluint iVector=0; iVector<vectors.size(); ++iVector) {
        serializer.addValues<T,3>(vectors[iVector]);
    }
}

template<typename T, template<typename U> class Descriptor>
void VisualParticle3D<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
{
    Particle3D<T,Descriptor>::unserialize(unserializer);
    pluint scalarSize;
    unserializer.readValue(scalarSize);
    PLB_ASSERT(scalarSize == (pluint)scalars.size());
    for (pluint iScalar=0; iScalar<scalars.size(); ++iScalar) {
        unserializer.readValue(scalars[iScalar]);
    }
    pluint vectorSize;
    unserializer.readValue(vectorSize);
    PLB_ASSERT(vectorSize == (pluint)vectors.size());
    for (pluint iVector=0; iVector<vectors.size(); ++iVector) {
        unserializer.readValues<T,3>(vectors[iVector]);
    }
}

template<typename T, template<typename U> class Descriptor>
VisualParticle3D<T,Descriptor>* VisualParticle3D<T,Descriptor>::clone() const {
    return new VisualParticle3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
bool VisualParticle3D<T,Descriptor>::getVector(plint whichVector, Array<T,3>& vector) const {
    if (whichVector<(plint)vectors.size()) {
        vector = vectors[whichVector];
        return true;
    }
    return Particle3D<T,Descriptor>::getVector(whichVector, vector);
}

template<typename T, template<typename U> class Descriptor>
bool VisualParticle3D<T,Descriptor>::getScalar(plint whichScalar, T& scalar) const {
    if (whichScalar<(plint)scalars.size()) {
        scalar = scalars[whichScalar];
        return true;
    }
    return Particle3D<T,Descriptor>::getScalar(whichScalar, scalar);
}

template<typename T, template<typename U> class Descriptor>
bool VisualParticle3D<T,Descriptor>::setScalars(std::vector<T> const& scalars_)
{
    scalars = scalars_;
    return true;
}

template<typename T, template<typename U> class Descriptor>
bool VisualParticle3D<T,Descriptor>::setVectors(std::vector<Array<T,3> > const& vectors_)
{
    vectors = vectors_;
    return true;
}


/* *************** class MappingParticle3D ************************************ */

template<typename T, template<typename U> class Descriptor>
MappingParticle3D<T,Descriptor>::MappingParticle3D()
    : PointParticle3D<T,Descriptor>(),
      stage(0),
      pos1(T(),T(),T()),
      pos2(T(),T(),T())
{ }

template<typename T, template<typename U> class Descriptor>
MappingParticle3D<T,Descriptor>::MappingParticle3D(plint tag_, Array<T,3> const& position_, Array<T,3> const& velocity_)
    : PointParticle3D<T,Descriptor>(tag_, position_, velocity_),
      stage(0),
      pos1(T(),T(),T()),
      pos2(T(),T(),T())
{ }

template<typename T, template<typename U> class Descriptor>
void MappingParticle3D<T,Descriptor>::advance() {
    // Dynamics of a mapping particle.
    // Stage 0: The particle is behind the first surface and is flowing forward.
    // Stage 1: The particle is behind the second surface and is flowing forward.
    // Stage 2: The particle has crossed the second surface and is
    //          now walking back along a straight line to the first surface.
    // Stage 3: The particle is back at the first surface and does nothing any more.

    // Forward flowing.
    if (stage==0 || stage==1) {
        this->getPosition() += this->getVelocity();
    }
    // Backward walking.
    else if (stage==2) {
        Array<T,3> homeArrow = pos1-this->getPosition();
        T l = norm(homeArrow);
        T step=0.5;
        if (l<step) {
            stage=3;
            this->getPosition() = pos1;
        }
        else {
            this->getPosition() += homeArrow *step/l;
        }
    }
    // Crossed first surface?
    if (stage==0 && crossedSurface1()) {
        stage=1;
        pos1 = this->getPosition();
    }
    // Crossed second surface?
    if (stage==1 && crossedSurface2()) {
        stage=2;
        pos2 = this->getPosition();
    }
}

template<typename T, template<typename U> class Descriptor>
void MappingParticle3D<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    PointParticle3D<T,Descriptor>::serialize(serializer);
    serializer.addValue(stage);
    serializer.addValues<T,3>(pos1);
    serializer.addValues<T,3>(pos2);
}

template<typename T, template<typename U> class Descriptor>
void MappingParticle3D<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
{
    PointParticle3D<T,Descriptor>::unserialize(unserializer);
    unserializer.readValue(stage);
    unserializer.readValues<T,3>(pos1);
    unserializer.readValues<T,3>(pos2);
}

template<typename T, template<typename U> class Descriptor>
void MappingParticle3D<T,Descriptor>::rescale(int dxScale, int dtScale) {
    PointParticle3D<T,Descriptor>::rescale(dxScale, dtScale);
    T posDimDx = 1.;
    T posDimDt = 0.;
    T posScaleFactor = scaleFromReference(dxScale, posDimDx, dtScale, posDimDt);
    pos1 *= posScaleFactor;
    pos2 *= posScaleFactor;
}


/* *************** class MappingParticleZslice3D ************************************ */

template<typename T, template<typename U> class Descriptor>
int MappingParticleZslice3D<T,Descriptor>::id =
        meta::registerGenericParticle3D<T,Descriptor,MappingParticleZslice3D<T,Descriptor> >("MappingZslice");

template<typename T, template<typename U> class Descriptor>
MappingParticleZslice3D<T,Descriptor>::MappingParticleZslice3D()
    : MappingParticle3D<T,Descriptor>(),
      zSlice1(T()), zSlice2(T())
{ }

template<typename T, template<typename U> class Descriptor>
MappingParticleZslice3D<T,Descriptor>::MappingParticleZslice3D (
        plint tag_, Array<T,3> const& position_, Array<T,3> const& velocity_, T zSlice1_, T zSlice2_)
    : MappingParticle3D<T,Descriptor>(tag_, position_, velocity_),
      zSlice1(zSlice1_),
      zSlice2(zSlice2_)
{ }

template<typename T, template<typename U> class Descriptor>
int MappingParticleZslice3D<T,Descriptor>::getId() const 
{
    return id;
}

template<typename T, template<typename U> class Descriptor>
bool MappingParticleZslice3D<T,Descriptor>::crossedSurface1() const {
    if (zSlice2>zSlice1) {
        return this->getPosition()[2] > zSlice1;
    }
    else {
        return this->getPosition()[2] < zSlice1;
    }
}

template<typename T, template<typename U> class Descriptor>
bool MappingParticleZslice3D<T,Descriptor>::crossedSurface2() const {
    if (zSlice2>zSlice1) {
        return this->getPosition()[2] > zSlice2;
    }
    else {
        return this->getPosition()[2] < zSlice2;
    }
}

template<typename T, template<typename U> class Descriptor>
T MappingParticleZslice3D<T,Descriptor>::getSurfaceDistance() const {
    return fabs(zSlice2 - zSlice1);
}

template<typename T, template<typename U> class Descriptor>
MappingParticleZslice3D<T,Descriptor>* MappingParticleZslice3D<T,Descriptor>::clone() const {
    return new MappingParticleZslice3D<T,Descriptor>(*this);
}


template<typename T, template<typename U> class Descriptor>
void MappingParticleZslice3D<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    MappingParticle3D<T,Descriptor>::serialize(serializer);
    serializer.addValue(zSlice1);
    serializer.addValue(zSlice2);
}

template<typename T, template<typename U> class Descriptor>
void MappingParticleZslice3D<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
{
    MappingParticle3D<T,Descriptor>::unserialize(unserializer);
    unserializer.readValue(zSlice1);
    unserializer.readValue(zSlice2);
}


template<typename T, template<typename U> class Descriptor>
void MappingParticleZslice3D<T,Descriptor>::rescale(int dxScale, int dtScale) {
    MappingParticle3D<T,Descriptor>::rescale(dxScale, dtScale);
    T posDimDx = 1.;
    T posDimDt = 0.;
    T posScaleFactor = scaleFromReference(dxScale, posDimDx, dtScale, posDimDt);
    zSlice1 *= posScaleFactor;
    zSlice2 *= posScaleFactor;
}


/* *************** class PlaneMappingParticle3D ************************************ */

template<typename T, template<typename U> class Descriptor>
int PlaneMappingParticle3D<T,Descriptor>::id =
        meta::registerGenericParticle3D<T,Descriptor,PlaneMappingParticle3D<T,Descriptor> >("PlaneMapping");

template<typename T, template<typename U> class Descriptor>
PlaneMappingParticle3D<T,Descriptor>::PlaneMappingParticle3D()
    : PointParticle3D<T,Descriptor>(),
      residenceTime((T) 0),
      advanceBackwardInTime(true),
      timeScaling((T) 1),
      reachedTerminalPlane(false)
{
    terminalPlane.point = Array<T,3>((T) 0, (T) 0, (T) 0);
    terminalPlane.normal = Array<T,3>((T) 1, (T) 0, (T) 0);
    initialPosition = Array<T,3>((T) 0, (T) 0, (T) 0);
    initialVelocity = Array<T,3>((T) 0, (T) 0, (T) 0);
    terminalVelocity = Array<T,3>((T) 0, (T) 0, (T) 0);
}

template<typename T, template<typename U> class Descriptor>
PlaneMappingParticle3D<T,Descriptor>::PlaneMappingParticle3D(plint tag_, Array<T,3> const& position_,
        Array<T,3> const& velocity_, Plane<T> const& terminalPlane_, bool advanceBackwardInTime_, T timeScaling_)
    : PointParticle3D<T,Descriptor>(tag_, position_, velocity_),
      terminalPlane(terminalPlane_),
      initialPosition(position_),
      initialVelocity(velocity_),
      residenceTime((T) 0),
      advanceBackwardInTime(advanceBackwardInTime_),
      timeScaling(timeScaling_),
      reachedTerminalPlane(false)
{
    terminalVelocity = Array<T,3>((T) 0, (T) 0, (T) 0);
    T eps = 100.0 * std::numeric_limits<T>::epsilon();
    T nn = norm(terminalPlane.normal);
    PLB_ASSERT(nn > eps);
    terminalPlane.normal /= nn;
    PLB_ASSERT(timeScaling > (T) 1 || fabs(timeScaling - (T) 1) <= eps);
}

template<typename T, template<typename U> class Descriptor>
PlaneMappingParticle3D<T,Descriptor>* PlaneMappingParticle3D<T,Descriptor>::clone() const
{
    return new PlaneMappingParticle3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
int PlaneMappingParticle3D<T,Descriptor>::getId() const 
{
    return id;
}

template<typename T, template<typename U> class Descriptor>
void PlaneMappingParticle3D<T,Descriptor>::advance()
{
    if (!reachedTerminalPlane) {
        Array<T,3> oldPosition = this->getPosition();
        if (advanceBackwardInTime) {
            this->getPosition() -= this->getVelocity();
        } else {
            this->getPosition() += this->getVelocity();
        }
        Array<T,3> newPosition = this->getPosition();
        residenceTime += timeScaling;
        T distanceFromPlane;
        if (crossedTerminalPlane(distanceFromPlane)) {
            Precision precision;
            if (sizeof(T) == sizeof(float)) {
                precision = FLT;
            } else if (sizeof(T) == sizeof(double)) {
                precision = DBL;
            } else if (sizeof(T) == sizeof(long double)) {
                precision = LDBL;
            } else {
                PLB_ASSERT(false);
            }
            Array<T,3> terminalPosition((T) 0, (T) 0, (T) 0);
            T dt = 0.0;
            if (lineIntersectionWithPlane(terminalPlane, oldPosition, newPosition, precision, terminalPosition) != 1) {
                terminalPosition = newPosition - distanceFromPlane * terminalPlane.normal;
                T normalVelocityMagnitude = fabs(dot<T,3>(this->getVelocity(), terminalPlane.normal));
                dt = distanceFromPlane / normalVelocityMagnitude;
            } else {
                T distance = norm(newPosition - terminalPosition);
                T velocityMagnitude = norm(this->getVelocity());
                dt = distance / velocityMagnitude;
            }
            this->getPosition() = terminalPosition;
            terminalVelocity = this->getVelocity();
            residenceTime -= dt;
            reachedTerminalPlane = true;
            return;
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void PlaneMappingParticle3D<T,Descriptor>::serialize(HierarchicSerializer& serializer) const
{
    PointParticle3D<T,Descriptor>::serialize(serializer);
    serializer.addValues<T,3>(terminalPlane.point);
    serializer.addValues<T,3>(terminalPlane.normal);
    serializer.addValues<T,3>(initialPosition);
    serializer.addValues<T,3>(initialVelocity);
    serializer.addValues<T,3>(terminalVelocity);
    serializer.addValue(residenceTime);
    serializer.addValue(advanceBackwardInTime);
    serializer.addValue(timeScaling);
    serializer.addValue(reachedTerminalPlane);
}

template<typename T, template<typename U> class Descriptor>
void PlaneMappingParticle3D<T,Descriptor>::unserialize(HierarchicUnserializer& unserializer)
{
    PointParticle3D<T,Descriptor>::unserialize(unserializer);
    unserializer.readValues<T,3>(terminalPlane.point);
    unserializer.readValues<T,3>(terminalPlane.normal);
    unserializer.readValues<T,3>(initialPosition);
    unserializer.readValues<T,3>(initialVelocity);
    unserializer.readValues<T,3>(terminalVelocity);
    unserializer.readValue(residenceTime);
    unserializer.readValue(advanceBackwardInTime);
    unserializer.readValue(timeScaling);
    unserializer.readValue(reachedTerminalPlane);
}

template<typename T, template<typename U> class Descriptor>
void PlaneMappingParticle3D<T,Descriptor>::rescale(int dxScale, int dtScale)
{
    PointParticle3D<T,Descriptor>::rescale(dxScale, dtScale);

    T posDimDx = 1.0;
    T posDimDt = 0.0;
    T posScaleFactor = scaleFromReference(dxScale, posDimDx, dtScale, posDimDt);
    initialPosition *= posScaleFactor;
    terminalPlane.point *= posScaleFactor;

    T velDimDx = 1.0;
    T velDimDt = -1.0;
    T velScaleFactor = scaleFromReference(dxScale, velDimDx, dtScale, velDimDt);
    initialVelocity *= velScaleFactor;
    terminalVelocity *= velScaleFactor;

    T timeDimDx = 0.0;
    T timeDimDt = 1.0;
    T timeScaleFactor = scaleFromReference(dxScale, timeDimDx, dtScale, timeDimDt);
    residenceTime *= timeScaleFactor;
    // TODO: Is it correct not to rescale the timeScaling?
}

template<typename T, template<typename U> class Descriptor>
bool PlaneMappingParticle3D<T,Descriptor>::crossedTerminalPlane(T& distanceFromPlane) const
{
    Array<T,3> x = this->getPosition();
    Array<T,3> dx = x - terminalPlane.point;
    distanceFromPlane = dot<T,3>(dx, terminalPlane.normal);

    if (distanceFromPlane >= (T) 0) {
        return true;
    }

    return false;
}

}  // namespace plb

#endif  // VISUAL_PARTICLE_3D_H
