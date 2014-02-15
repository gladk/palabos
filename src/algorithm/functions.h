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

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "core/array.h"

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

}  // namespace util

}  // namespace plb

#endif  // FUNCTIONS_H

