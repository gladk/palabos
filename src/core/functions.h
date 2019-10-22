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

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "core/array.h"
#include "core/globalDefs.h"

#include <cmath>

namespace plb {

namespace util {

template<typename T, pluint size>
struct TimeDependentFunction {
    virtual ~TimeDependentFunction() { }
    virtual Array<T,size> operator()(T t) const =0;
    virtual TimeDependentFunction<T,size>* clone() const=0;
};

struct SelectInt {
    virtual ~SelectInt() { }
    virtual bool operator() (plint value) const=0;
    virtual SelectInt* clone() const=0;
};

class SelectConstInt : public SelectInt {
public:
    SelectConstInt(plint constVal_) : constVal(constVal_) { }
    virtual bool operator() (plint value) const {
        return value == constVal;
    }
    virtual SelectConstInt* clone() const {
        return new SelectConstInt(*this);
    }
private:
    plint constVal;
};

class SelectLargerEqualInt : public SelectInt {
public:
    SelectLargerEqualInt(plint threshold_) : threshold(threshold_) { }
    virtual bool operator() (plint value) const {
        return value >= threshold;
    }
    virtual SelectLargerEqualInt* clone() const {
        return new SelectLargerEqualInt(*this);
    }
private:
    plint threshold;
};

// Smoothly increasing function from 0.0 to 1.0 using a sine.
template<typename T>
T sinIncreasingFunction(T t, T maxT)
{
    static T pi = std::acos((T) -1);

    PLB_ASSERT(util::greaterEqual(maxT, (T) 0));

    if (t >= maxT) {
        return((T) 1);
    }

    if (t < 0) {
        return((T) 0);
    }

    return(std::sin(pi * t / (2.0 * maxT)));
}

template<typename T>
T sinIncreasingFunctionIntegral(T t1, T t2, T maxT)
{
    static T pi = std::acos((T) -1);

    PLB_ASSERT(util::greaterEqual(maxT, (T) 0));

    if (util::fpequal(t1, t2)) {
        return ((T) 0);
    }

    T sign = (T) 1;
    if (util::lessThan(t2, t1)) {
        std::swap(t1, t2);
        sign = (T) -1;
    }

    T tmin = t1;
    T tmax = t2;

    if (util::lessEqual(t1, (T) 0)) {
        if (util::lessEqual(t2, (T) 0)) {
            return((T) 0);
        } else {
            tmin = (T) 0;
        }
    }

    T rightIntegral = (T) 0;
    if (util::greaterEqual(t2, maxT)) {
        if (util::greaterEqual(t1, maxT)) {
            return sign * (t2 - t1);
        } else {
            tmax = maxT;
            rightIntegral = (t2 - maxT);
        }
    }

    T middleIntegral = ((T) 2 * maxT / pi) * (std::cos(pi * tmin / ((T) 2 * maxT)) - std::cos(pi * tmax / ((T) 2 * maxT)));

    return sign * (middleIntegral + rightIntegral);
}

// Smoothly increasing function from 0.0 to 1.0 using the hyperbolic tangent function.
template<typename T>
T tanhIncreasingFunction(T t, T maxT, T sigma = 6.0)
{
    PLB_ASSERT(util::greaterEqual(maxT, (T) 0));

    if (t >= maxT) {
        return((T) 1);
    }

    if (t < 0) {
        return((T) 0);
    }

    return((T) 0.5 * ((T) 1 + std::tanh(sigma / maxT * (t - (T) 0.5 * maxT))));
}

}  // namespace util

}  // namespace plb

#endif  // FUNCTIONS_H

