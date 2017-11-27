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


#ifndef VECTOR_FUNCTION_3D_H
#define VECTOR_FUNCTION_3D_H

#include "core/array.h"
#include "core/functions.h"
#include "core/globalDefs.h"
#include "core/util.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiBlock/multiBlockLattice3D.hh"

#include <cmath>

namespace plb {

template<typename T>
struct VectorFunction3D {
    virtual ~VectorFunction3D() { }
    virtual Array<T,3> operator()(Array<T,3> const& position) const = 0;
    virtual VectorFunction3D<T>* clone() const = 0;
};

template<typename T>
class ConstantVectorFunction3D : public VectorFunction3D<T> {
public:
    ConstantVectorFunction3D(Array<T,3> const& constantVector_)
        : constantVector(constantVector_)
    { }

    virtual Array<T,3> operator()(Array<T,3> const& position) const
    {
        return constantVector;
    }

    virtual ConstantVectorFunction3D<T>* clone() const
    {
        return new ConstantVectorFunction3D<T>(*this);
    }
private:
    Array<T,3> constantVector;
};

template<typename T>
class IdentityVectorFunction3D : public VectorFunction3D<T> {
public:
    virtual Array<T,3> operator()(Array<T,3> const& position) const
    {
        return position;
    }

    virtual IdentityVectorFunction3D<T>* clone() const
    {
        return new IdentityVectorFunction3D<T>(*this);
    }
};

template<typename T>
class DiscreteTranslationalPositionFunction3D : public VectorFunction3D<T> {
public:
    DiscreteTranslationalPositionFunction3D(Array<T,3> const& velocity_)
        : velocity(velocity_)
    { }

    Array<T,3> getVelocity() const
    {
        return velocity;
    }

    virtual Array<T,3> operator()(Array<T,3> const& position) const
    {
        return position + velocity;
    }

    virtual DiscreteTranslationalPositionFunction3D<T>* clone() const
    {
        return new DiscreteTranslationalPositionFunction3D<T>(*this);
    }
private:
    Array<T,3> velocity;
};

template<typename T>
class DiscreteRotationalVelocityFunction3D : public VectorFunction3D<T> {
public:
    DiscreteRotationalVelocityFunction3D(Array<T,3> const& angularVelocity_, Array<T,3> const& pointOnRotationAxis_)
        : angularVelocity(angularVelocity_),
          pointOnRotationAxis(pointOnRotationAxis_)
    {
        normAngularVelocity = norm(angularVelocity);
        if (!util::isZero(normAngularVelocity)) {
            rotationAxisUnitVector = angularVelocity / normAngularVelocity;
        } else {
            angularVelocity = Array<T,3>::zero();
            rotationAxisUnitVector = Array<T,3>((T) 1, (T) 0, (T) 0); // Array<T,3>::zero();
            normAngularVelocity = (T) 0;
        }
    }

    DiscreteRotationalVelocityFunction3D(Array<T,3> const& angularVelocity_)
        : angularVelocity(angularVelocity_),
          pointOnRotationAxis(Array<T,3>::zero())
    {
        normAngularVelocity = norm(angularVelocity);
        if (!util::isZero(normAngularVelocity)) {
            rotationAxisUnitVector = angularVelocity / normAngularVelocity;
        } else {
            angularVelocity = Array<T,3>::zero();
            rotationAxisUnitVector = Array<T,3>((T) 1, (T) 0, (T) 0); // Array<T,3>::zero();
            normAngularVelocity = (T) 0;
        }
    }

    Array<T,3> getAngularVelocity() const
    {
        return angularVelocity;
    }

    // The rotation angle between two time instants is defined as the time integral of the
    // norm of the angular velocity.
    T getRotationAngle(T t1, T t2) const
    {
        return normAngularVelocity * (t2 - t1);
    }

    // The current rotation angle is defined as the rotation angle between the current
    // time instant t and and the time instant t + 1.
    T getRotationAngle() const
    {
        return normAngularVelocity;
    }

    virtual Array<T,3> operator()(Array<T,3> const& position) const
    {
        return getDiscreteRotationalVelocity(position, normAngularVelocity, rotationAxisUnitVector, pointOnRotationAxis);
    }

    virtual DiscreteRotationalVelocityFunction3D<T>* clone() const
    {
        return new DiscreteRotationalVelocityFunction3D<T>(*this);
    }
private:
    Array<T,3> angularVelocity, rotationAxisUnitVector, pointOnRotationAxis;
    T normAngularVelocity;
};

template<typename T>
class DiscreteRotationalPositionFunction3D : public VectorFunction3D<T> {
public:
    DiscreteRotationalPositionFunction3D(Array<T,3> const& angularVelocity_, Array<T,3> const& pointOnRotationAxis_)
        : angularVelocity(angularVelocity_),
          pointOnRotationAxis(pointOnRotationAxis_)
    {
        normAngularVelocity = norm(angularVelocity);
        if (!util::isZero(normAngularVelocity)) {
            rotationAxisUnitVector = angularVelocity / normAngularVelocity;
        } else {
            angularVelocity = Array<T,3>::zero();
            rotationAxisUnitVector = Array<T,3>((T) 1, (T) 0, (T) 0); // Array<T,3>::zero();
            normAngularVelocity = (T) 0;
        }
    }

    DiscreteRotationalPositionFunction3D(Array<T,3> const& angularVelocity_)
        : angularVelocity(angularVelocity_),
          pointOnRotationAxis(Array<T,3>::zero())
    {
        normAngularVelocity = norm(angularVelocity);
        if (!util::isZero(normAngularVelocity)) {
            rotationAxisUnitVector = angularVelocity / normAngularVelocity;
        } else {
            angularVelocity = Array<T,3>::zero();
            rotationAxisUnitVector = Array<T,3>((T) 1, (T) 0, (T) 0); // Array<T,3>::zero();
            normAngularVelocity = (T) 0;
        }
    }

    Array<T,3> getAngularVelocity() const
    {
        return angularVelocity;
    }

    // The rotation angle between two time instants is defined as the time integral of the
    // norm of the angular velocity.
    T getRotationAngle(T t1, T t2) const
    {
        return normAngularVelocity * (t2 - t1);
    }

    // The current rotation angle is defined as the rotation angle between the current
    // time instant t and and the time instant t + 1.
    T getRotationAngle() const
    {
        return normAngularVelocity;
    }

    virtual Array<T,3> operator()(Array<T,3> const& position) const
    {
        return getRotatedPosition(position, normAngularVelocity, rotationAxisUnitVector, pointOnRotationAxis);
    }

    virtual DiscreteRotationalPositionFunction3D<T>* clone() const
    {
        return new DiscreteRotationalPositionFunction3D<T>(*this);
    }
private:
    Array<T,3> angularVelocity, rotationAxisUnitVector, pointOnRotationAxis;
    T normAngularVelocity;
};

template<typename T>
class ExactRotationalVelocityFunction3D : public VectorFunction3D<T> {
public:
    ExactRotationalVelocityFunction3D(Array<T,3> const& angularVelocity_, Array<T,3> const& pointOnRotationAxis_)
        : angularVelocity(angularVelocity_),
          pointOnRotationAxis(pointOnRotationAxis_)
    {
        normAngularVelocity = norm(angularVelocity);
        if (!util::isZero(normAngularVelocity)) {
            rotationAxisUnitVector = angularVelocity / normAngularVelocity;
        } else {
            angularVelocity = Array<T,3>::zero();
            rotationAxisUnitVector = Array<T,3>((T) 1, (T) 0, (T) 0); // Array<T,3>::zero();
            normAngularVelocity = (T) 0;
        }
    }

    ExactRotationalVelocityFunction3D(Array<T,3> const& angularVelocity_)
        : angularVelocity(angularVelocity_),
          pointOnRotationAxis(Array<T,3>::zero())
    {
        normAngularVelocity = norm(angularVelocity);
        if (!util::isZero(normAngularVelocity)) {
            rotationAxisUnitVector = angularVelocity / normAngularVelocity;
        } else {
            angularVelocity = Array<T,3>::zero();
            rotationAxisUnitVector = Array<T,3>((T) 1, (T) 0, (T) 0); // Array<T,3>::zero();
            normAngularVelocity = (T) 0;
        }
    }

    Array<T,3> getAngularVelocity() const
    {
        return angularVelocity;
    }

    // The rotation angle between two time instants is defined as the time integral of the
    // norm of the angular velocity.
    T getRotationAngle(T t1, T t2) const
    {
        return normAngularVelocity * (t2 - t1);
    }

    // The current rotation angle is defined as the rotation angle between the current
    // time instant t and and the time instant t + 1.
    T getRotationAngle() const
    {
        return normAngularVelocity;
    }

    virtual Array<T,3> operator()(Array<T,3> const& position) const
    {
        return getExactRotationalVelocity(position, angularVelocity, pointOnRotationAxis);
    }

    virtual ExactRotationalVelocityFunction3D<T>* clone() const
    {
        return new ExactRotationalVelocityFunction3D<T>(*this);
    }
private:
    Array<T,3> angularVelocity, rotationAxisUnitVector, pointOnRotationAxis;
    T normAngularVelocity;
};

template<typename T, template<typename U> class Descriptor> 
class IncreasingDiscreteRotationalVelocityFunction3D : public VectorFunction3D<T> {
public:
    IncreasingDiscreteRotationalVelocityFunction3D(Array<T,3> const& maxAngularVelocity_,
            Array<T,3> const& pointOnRotationAxis_, MultiBlockLattice3D<T,Descriptor> const& lattice_,
            T tOffset_, T maxT_)
        : maxAngularVelocity(maxAngularVelocity_),
          pointOnRotationAxis(pointOnRotationAxis_),
          lattice(lattice_),
          tOffset(tOffset_),
          maxT(maxT_)
    {
        normMaxAngularVelocity = norm(maxAngularVelocity);
        if (!util::isZero(normMaxAngularVelocity)) {
            rotationAxisUnitVector = maxAngularVelocity / normMaxAngularVelocity;
        } else {
            maxAngularVelocity = Array<T,3>::zero();
            rotationAxisUnitVector = Array<T,3>((T) 1, (T) 0, (T) 0); // Array<T,3>::zero();
            normMaxAngularVelocity = (T) 0;
        }
    }

    IncreasingDiscreteRotationalVelocityFunction3D(Array<T,3> const& maxAngularVelocity_,
            MultiBlockLattice3D<T,Descriptor> const& lattice_, T tOffset_, T maxT_)
        : maxAngularVelocity(maxAngularVelocity_),
          pointOnRotationAxis(Array<T,3>::zero()),
          lattice(lattice_),
          tOffset(tOffset_),
          maxT(maxT_)
    {
        normMaxAngularVelocity = norm(maxAngularVelocity);
        if (!util::isZero(normMaxAngularVelocity)) {
            rotationAxisUnitVector = maxAngularVelocity / normMaxAngularVelocity;
        } else {
            maxAngularVelocity = Array<T,3>::zero();
            rotationAxisUnitVector = Array<T,3>((T) 1, (T) 0, (T) 0); // Array<T,3>::zero();
            normMaxAngularVelocity = (T) 0;
        }
    }

    Array<T,3> getAngularVelocity() const
    {
        T t = (T) lattice.getTimeCounter().getTime() + tOffset;
        Array<T,3> angularVelocity = util::sinIncreasingFunction<T>(t, maxT) * maxAngularVelocity;
        return angularVelocity;
    }

    // The rotation angle between two time instants is defined as the time integral of the
    // norm of the angular velocity.
    T getRotationAngle(T t1, T t2) const
    {
        T theta = util::sinIncreasingFunctionIntegral<T>(t1, t2, maxT) * normMaxAngularVelocity;
        return theta;
    }

    // The current rotation angle is defined as the rotation angle between the current
    // time instant t and and the time instant t + 1.
    T getRotationAngle() const
    {
        T t = (T) lattice.getTimeCounter().getTime() + tOffset;
        T theta = getRotationAngle(t, t + (T) 1);
        return theta;
    }

    virtual Array<T,3> operator()(Array<T,3> const& position) const
    {
        return getDiscreteRotationalVelocity(position, getRotationAngle(), rotationAxisUnitVector, pointOnRotationAxis);
    }

    virtual IncreasingDiscreteRotationalVelocityFunction3D<T,Descriptor>* clone() const
    {
        return new IncreasingDiscreteRotationalVelocityFunction3D<T,Descriptor>(*this);
    }
private:
    Array<T,3> maxAngularVelocity, rotationAxisUnitVector, pointOnRotationAxis;
    MultiBlockLattice3D<T,Descriptor> const& lattice;
    T tOffset, maxT, normMaxAngularVelocity;
};

template<typename T, template<typename U> class Descriptor> 
class IncreasingDiscreteRotationalPositionFunction3D : public VectorFunction3D<T> {
public:
    IncreasingDiscreteRotationalPositionFunction3D(Array<T,3> const& maxAngularVelocity_,
            Array<T,3> const& pointOnRotationAxis_, MultiBlockLattice3D<T,Descriptor> const& lattice_,
            T tOffset_, T maxT_)
        : maxAngularVelocity(maxAngularVelocity_),
          pointOnRotationAxis(pointOnRotationAxis_),
          lattice(lattice_),
          tOffset(tOffset_),
          maxT(maxT_)
    {
        normMaxAngularVelocity = norm(maxAngularVelocity);
        if (!util::isZero(normMaxAngularVelocity)) {
            rotationAxisUnitVector = maxAngularVelocity / normMaxAngularVelocity;
        } else {
            maxAngularVelocity = Array<T,3>::zero();
            rotationAxisUnitVector = Array<T,3>((T) 1, (T) 0, (T) 0); // Array<T,3>::zero();
            normMaxAngularVelocity = (T) 0;
        }
    }

    IncreasingDiscreteRotationalPositionFunction3D(Array<T,3> const& maxAngularVelocity_,
            MultiBlockLattice3D<T,Descriptor> const& lattice_, T tOffset_, T maxT_)
        : maxAngularVelocity(maxAngularVelocity_),
          pointOnRotationAxis(Array<T,3>::zero()),
          lattice(lattice_),
          tOffset(tOffset_),
          maxT(maxT_)
    {
        normMaxAngularVelocity = norm(maxAngularVelocity);
        if (!util::isZero(normMaxAngularVelocity)) {
            rotationAxisUnitVector = maxAngularVelocity / normMaxAngularVelocity;
        } else {
            maxAngularVelocity = Array<T,3>::zero();
            rotationAxisUnitVector = Array<T,3>((T) 1, (T) 0, (T) 0); // Array<T,3>::zero();
            normMaxAngularVelocity = (T) 0;
        }
    }

    Array<T,3> getAngularVelocity() const
    {
        T t = (T) lattice.getTimeCounter().getTime() + tOffset;
        Array<T,3> angularVelocity = util::sinIncreasingFunction<T>(t, maxT) * maxAngularVelocity;
        return angularVelocity;
    }

    // The rotation angle between two time instants is defined as the time integral of the
    // norm of the angular velocity.
    T getRotationAngle(T t1, T t2) const
    {
        T theta = util::sinIncreasingFunctionIntegral<T>(t1, t2, maxT) * normMaxAngularVelocity;
        return theta;
    }

    // The current rotation angle is defined as the rotation angle between the current
    // time instant t and and the time instant t + 1.
    T getRotationAngle() const
    {
        T t = (T) lattice.getTimeCounter().getTime() + tOffset;
        T theta = getRotationAngle(t, t + (T) 1);
        return theta;
    }

    virtual Array<T,3> operator()(Array<T,3> const& position) const
    {
        return getRotatedPosition(position, getRotationAngle(), rotationAxisUnitVector, pointOnRotationAxis);
    }

    virtual IncreasingDiscreteRotationalPositionFunction3D<T,Descriptor>* clone() const
    {
        return new IncreasingDiscreteRotationalPositionFunction3D<T,Descriptor>(*this);
    }
private:
    Array<T,3> maxAngularVelocity, rotationAxisUnitVector, pointOnRotationAxis;
    MultiBlockLattice3D<T,Descriptor> const& lattice;
    T tOffset, maxT, normMaxAngularVelocity;
};

template<typename T, template<typename U> class Descriptor> 
class IncreasingExactRotationalVelocityFunction3D : public VectorFunction3D<T> {
public:
    IncreasingExactRotationalVelocityFunction3D(Array<T,3> const& maxAngularVelocity_,
            Array<T,3> const& pointOnRotationAxis_, MultiBlockLattice3D<T,Descriptor> const& lattice_,
            T tOffset_, T maxT_)
        : maxAngularVelocity(maxAngularVelocity_),
          pointOnRotationAxis(pointOnRotationAxis_),
          lattice(lattice_),
          tOffset(tOffset_),
          maxT(maxT_)
    {
        normMaxAngularVelocity = norm(maxAngularVelocity);
        if (!util::isZero(normMaxAngularVelocity)) {
            rotationAxisUnitVector = maxAngularVelocity / normMaxAngularVelocity;
        } else {
            maxAngularVelocity = Array<T,3>::zero();
            rotationAxisUnitVector = Array<T,3>((T) 1, (T) 0, (T) 0); // Array<T,3>::zero();
            normMaxAngularVelocity = (T) 0;
        }
    }

    IncreasingExactRotationalVelocityFunction3D(Array<T,3> const& maxAngularVelocity_,
            MultiBlockLattice3D<T,Descriptor> const& lattice_, T tOffset_, T maxT_)
        : maxAngularVelocity(maxAngularVelocity_),
          pointOnRotationAxis(Array<T,3>::zero()),
          lattice(lattice_),
          tOffset(tOffset_),
          maxT(maxT_)
    {
        normMaxAngularVelocity = norm(maxAngularVelocity);
        if (!util::isZero(normMaxAngularVelocity)) {
            rotationAxisUnitVector = maxAngularVelocity / normMaxAngularVelocity;
        } else {
            maxAngularVelocity = Array<T,3>::zero();
            rotationAxisUnitVector = Array<T,3>((T) 1, (T) 0, (T) 0); // Array<T,3>::zero();
            normMaxAngularVelocity = (T) 0;
        }
    }

    Array<T,3> getAngularVelocity() const
    {
        T t = (T) lattice.getTimeCounter().getTime() + tOffset;
        Array<T,3> angularVelocity = util::sinIncreasingFunction<T>(t, maxT) * maxAngularVelocity;
        return angularVelocity;
    }

    // The rotation angle between two time instants is defined as the time integral of the
    // norm of the angular velocity.
    T getRotationAngle(T t1, T t2) const
    {
        T theta = util::sinIncreasingFunctionIntegral<T>(t1, t2, maxT) * normMaxAngularVelocity;
        return theta;
    }

    // The current rotation angle is defined as the rotation angle between the current
    // time instant t and and the time instant t + 1.
    T getRotationAngle() const
    {
        T t = (T) lattice.getTimeCounter().getTime() + tOffset;
        T theta = getRotationAngle(t, t + (T) 1);
        return theta;
    }

    virtual Array<T,3> operator()(Array<T,3> const& position) const
    {
        return getExactRotationalVelocity(position, getAngularVelocity(), pointOnRotationAxis);
    }

    virtual IncreasingExactRotationalVelocityFunction3D<T,Descriptor>* clone() const
    {
        return new IncreasingExactRotationalVelocityFunction3D<T,Descriptor>(*this);
    }
private:
    Array<T,3> maxAngularVelocity, rotationAxisUnitVector, pointOnRotationAxis;
    MultiBlockLattice3D<T,Descriptor> const& lattice;
    T tOffset, maxT, normMaxAngularVelocity;
};

template<typename T, template<typename U> class Descriptor> 
class HarmonicVectorFunction3D : public VectorFunction3D<T> {
public:
    HarmonicVectorFunction3D(Array<T,3> const& vectorAmplitude_, T angularFrequency_, T phase_,
            MultiBlockLattice3D<T,Descriptor> const& lattice_, T tOffset_)
        : vectorAmplitude(vectorAmplitude_),
          angularFrequency(angularFrequency_),
          phase(phase_),
          lattice(lattice_),
          tOffset(tOffset_)
    { }

    virtual Array<T,3> operator()(Array<T,3> const& position) const
    {
        T t = (T) lattice.getTimeCounter().getTime() + tOffset;
        return std::cos(angularFrequency * t + phase) * vectorAmplitude;
    }

    virtual HarmonicVectorFunction3D<T,Descriptor>* clone() const
    {
        return new HarmonicVectorFunction3D<T,Descriptor>(*this);
    }
private:
    Array<T,3> vectorAmplitude;
    T angularFrequency, phase;
    MultiBlockLattice3D<T,Descriptor> const& lattice;
    T tOffset;
};

} // namespace plb

#endif  // VECTOR_FUNCTION_3D_H
