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
 * Set of functions commonly used in LB computations
 *  -- header file
 */
#ifndef UTIL_H
#define UTIL_H

#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include "core/runTimeDiagnostics.h"
#include "core/plbTypenames.h"
#include <algorithm>
#include <string>
#include <sstream>
#include <cmath>
#include <cstring>
#include <utility>
#include <vector>
#include <set>
#include <limits>

namespace plb {

namespace util {

/// Compute square of a scalar value
template<typename T>
T sqr(T arg) {
    return arg*arg;
}

/// Kronecker delta
inline plint kronDelta(plint iA, plint iB) {
    if (iA == iB) return 1;
    else return 0;
}



/// Test equality between two boolean values.
/** As we were testing with GCC 4.3.3, we observed that the
 *  equality operator (==) has an unexpected behavior. It
 *  compares the bools as if they were ints, and sometimes
 *  returns false, even if both bools are true. We don't know
 *  if this is a weird definition of the C++ language or a compiler
 *  bug, but in both cases, the present function offers a workaround.
 */
inline bool boolIsEqual(bool val1, bool val2) {
    return ( (val1 && val2) || !(val1 || val2) );
}

/// Round to next plint value
template<typename T>
plint roundToInt(T value) {
    return value>0 ?
               static_cast<plint>(value+(T)0.5) :
               static_cast<plint>(value-(T)0.5);
}

/// Round a signed integer to the next larger value which is divisible by step.
inline plint roundUp(plint value, plint step) {
    plint modulo = value % step;
    plint result = value-modulo;
    if (modulo>0) {
        result += step;
    }
    return result;
}

/// Round a signed integer to the next smaller value which is divisible by step.
inline plint roundDown(plint value, plint step) {
    plint modulo = value % step;
    plint result = value-modulo;
    if (modulo<0) {
        result -= step;
    }
    return result;
}

inline double twoToThePower(int power) {
    switch(power) {
        case 0: return 1.;
        case 1: return 2.;
        case 2: return 4.;
        case 3: return 8.;
        case 4: return 16.;
        case 5: return 32.;
        case 6: return 64.;
        case 7: return 128.;
        case 8: return 256.;
        case -1: return 1./2.;
        case -2: return 1./4.;
        case -3: return 1./8.;
        case -4: return 1./16.;
        case -5: return 1./32.;
        case -6: return 1./64.;
        case -7: return 1./128.;
        case -8: return 1./256.;
        default: return std::pow(2., (double)power);
    }
}

inline plint intTwoToThePower(int power)
{
    PLB_ASSERT( power>=0 );
    switch(power) {
        case 0: return 1;
        case 1: return 2;
        case 2: return 4;
        case 3: return 8;
        case 4: return 16;
        case 5: return 32;
        case 6: return 64;
        case 7: return 128;
        case 8: return 256;
        case 9: return 512;
        case 10: return 1024;
        default: return util::roundToInt(std::pow(2., (double)power));
    }
}

template<typename T>
std::string val2str(T val) {
    std::stringstream valstr;
    valstr << val;
    return valstr.str();
}

template<typename T1, typename T2>
std::string val2str(T1 val1, T2 val2) {
    std::stringstream valstr;
    valstr << val1 << " " << val2;
    return valstr.str();
}

template<typename T1, typename T2, typename T3>
std::string val2str(T1 val1, T2 val2, T3 val3) {
    std::stringstream valstr;
    valstr << val1 << " " << val2 << " " << val3;
    return valstr.str();
}

template<typename T>
void str2val(std::string str, T& val) {
    std::stringstream stream(str);
    if (!(stream >> val)) {
        plbLogicError (
                std::string("Could not convert string \"") + str +
                std::string("\" to type ") +
                std::string(NativeType<T>::getName()) );
    }
}

inline
std::vector<std::string> split(std::string str) {
    std::stringstream stream(str);
    std::string buf;
    std::vector<std::string> result;
    while (stream>>buf) {
        result.push_back(buf);
    }
    return result;
}

template<typename T>
std::string consume(std::string str, T& val) {
    std::vector<std::string> tokens = split(str);
    if (tokens.size()>=1) {
        str2val(tokens[0], val);
    }
    else {
        plbLogicError("Could not read value from empty string.");
    }
    std::string remaining;
    for (plint i=1; i<(plint)tokens.size(); ++i) {
        remaining += tokens[i];
        if (i != (plint)tokens.size()-1) {
            remaining += " ";
        }
    }
    return remaining;
}

template<typename T1, typename T2>
std::string consume(std::string str, T1& val1, T2& val2) {
    std::vector<std::string> tokens = split(str);
    if (tokens.size()>=2) {
        str2val(tokens[0], val1);
        str2val(tokens[1], val2);
    }
    else {
        plbLogicError (
                std::string("Could not read two values from string \"") + str +
                std::string("\"") );
    }
    std::string remaining;
    for (plint i=2; i<(plint)tokens.size(); ++i) {
        remaining += tokens[i];
        if (i != (plint)tokens.size()-1) {
            remaining += " ";
        }
    }
    return remaining;
}

template<typename T1, typename T2, typename T3>
std::string consume(std::string str, T1& val1, T2& val2, T3& val3) {
    std::vector<std::string> tokens = split(str);
    if (tokens.size()>=3) {
        str2val(tokens[0], val1);
        str2val(tokens[1], val2);
        str2val(tokens[2], val3);
    }
    else {
        plbLogicError (
                std::string("Could not read three values from string \"") + str +
                std::string("\"") );
    }
    std::string remaining;
    for (plint i=3; i<(plint)tokens.size(); ++i) {
        remaining += tokens[i];
        if (i != (plint)tokens.size()-1) {
            remaining += " ";
        }
    }
    return remaining;
}

inline std::string tolower(std::string arg) {
    std::string result(arg.size(), ' ');
    std::transform(arg.begin(), arg.end(), result.begin(), ::tolower);
    return result;
}

inline std::string toupper(std::string arg) {
    std::string result(arg.size(), ' ');
    std::transform(arg.begin(), arg.end(), result.begin(), ::toupper);
    return result;
}


/// A simple class for handling buffer memory
/** This class can be seen as a replacement of the std::vector
 *  template. It is less powerful, but at least, unlike std::vector,
 *  it treats the type bool appropriately (this is the only reason
 *  why class Buffer was written). Note that class Buffer is not
 *  STL compatible (it is not an STL container), but it is thread-
 *  safe.
 */
template<typename T>
class Buffer {
public:
    /// Default constructor allcoates no memory.
    Buffer()
        : bufferSize(0),
          data(0)
    { }
    /// Constructor with buffer size; buffer elements are not default-initialized.
    Buffer(pluint bufferSize_)
        : bufferSize(bufferSize_),
          data(new T[bufferSize])
    { }
    Buffer(Buffer<T> const& rhs)
        : bufferSize(rhs.bufferSize),
          data(new T[bufferSize])
    {
        for (pluint iData=0; iData<bufferSize; ++iData) {
            data[iData] = rhs.data[iData];
        }
    }
    ~Buffer() {
        delete [] data;
    }
    /// Assignment-operator is thread-safe.
    Buffer<T>& operator=(Buffer<T> const& rhs) {
        Buffer<T> newBuffer(rhs);
        swap(newBuffer);
        return *this;
    }
    /// Swap with other buffer; this operation cannot throw.
    void swap(Buffer<T>& rhs) {
        std::swap(bufferSize, rhs.bufferSize);
        std::swap(data, rhs.data);
    }
    /// Resize the buffer, and keep values from before.
    void resize(pluint newBufferSize) {
        if (newBufferSize > bufferSize) {
            Buffer<T> newBuffer(newBufferSize);
            for (pluint iData=0; iData<bufferSize; ++iData) {
                newBuffer.data[iData] = data[iData];
            }
            swap(newBuffer);
        }
    }
    /// Resize the buffer, but don't keep values from before.
    void reallocate(pluint newBufferSize) {
        if (newBufferSize > bufferSize) {
            Buffer<T> newBuffer(newBufferSize);
            swap(newBuffer);
        }
    }
    /// Pointer to raw data.
    T* get() {
        return data;
    }
    /// Const-pointer to raw data.
    T const* get() const {
        return data;
    }
    /// Element access.
    T& operator[](pluint index) {
        PLB_PRECONDITION( index<bufferSize );
        return data[index];
    }
    /// Const-element access.
    T const& operator[](pluint index) const {
        PLB_PRECONDITION( index<bufferSize );
        return data[index];
    }
    pluint size() const {
        return bufferSize;
    }
private:
    pluint bufferSize;
    T* data;
};

inline void linearRepartition(plint x0, plint x1, plint nBlocks,
                              std::vector<std::pair<plint,plint> >& ranges)
{
    PLB_PRECONDITION(nBlocks>0);
    plint totalSize = x1-x0+1;
    PLB_PRECONDITION( nBlocks<=totalSize );
    ranges.resize(nBlocks);
    plint basicLength = totalSize/nBlocks;
    plint currentPos=0;
    for (plint iRange=0; iRange<nBlocks; ++iRange) {
        plint currentLength = basicLength;
        if (iRange<totalSize%nBlocks) {
            ++currentLength;
        }
        ranges[iRange] = std::pair<plint,plint>(x0+currentPos, x0+currentPos+currentLength-1);
        currentPos+=currentLength;
    }
}

inline void linearBlockRepartition(plint x0, plint x1,
                                   plint wishedLength,
                                   std::vector<std::pair<plint,plint> >& ranges)
{
    plint totalSize = x1-x0+1;
    plint nBlocks = std::max((plint)1, totalSize/wishedLength);
    util::linearRepartition(x0, x1, nBlocks, ranges);
}

inline void linearEvenBlockRepartition(plint x0, plint x1, plint wishedLength, std::vector<std::pair<plint,plint> >& ranges)
{
    // Repartition of even numbers to blocks that contain even number of elements.
    PLB_PRECONDITION(wishedLength > 0 && wishedLength % 2 == 0);
    plint totalSize = x1 - x0 + 1;
    PLB_PRECONDITION(totalSize > 0 && totalSize % 2 == 0);
    plint nBlocks = std::max((plint) 1, totalSize / wishedLength);
    ranges.resize(nBlocks);
    plint basicLength = nBlocks == 1 ? totalSize : wishedLength;
    plint halfOfRemainder = (totalSize % basicLength) / 2;

    plint currentPos = 0;
    for (plint iRange = 0; iRange < nBlocks; iRange++) {
        plint currentLength = basicLength + 2 * (plint) (halfOfRemainder / nBlocks);
        if (iRange < halfOfRemainder % nBlocks) {
            currentLength += 2;
        }
        ranges[iRange] = std::pair<plint,plint>(x0 + currentPos, x0 + currentPos + currentLength - 1);
        currentPos += currentLength;
    }
    PLB_ASSERT(ranges[nBlocks - 1].second == x1);
}

/// Extend the size of an int-vector, and initialize new values to -1.
inline void extendVectorSize(std::vector<int>& vect, pluint fullSize) {
    pluint currentSize = vect.size();
    if (currentSize<fullSize) {
        vect.resize(fullSize);
        for (pluint i=currentSize; i<fullSize; ++i) {
            vect[i] = -1;
        }
    }
}

/// Test equality of two floating point numbers, with an accuracy
///   that is relative to the value of the arguments.
///   If at least one operand is zero, then an absolute tolerance
///   comparison is performed.
template<typename T>
inline bool fpequal(T x, T y, T eps = std::numeric_limits<T>::epsilon())
{
    // This is an attempt to treat underflow.
    if (std::fabs(x) <= eps) {
        x = (T) 0;
    }

    // This is an attempt to treat underflow.
    if (std::fabs(y) <= eps) {
        y = (T) 0;
    }

    if (x == (T) 0 || y == (T) 0) {
        return (std::fabs(x - y) <= eps);
    } else {
        return (std::fabs(x - y) <= eps * std::fabs(x) && std::fabs(x - y) <= eps * std::fabs(y));
    }
}

template<typename T>
inline bool greaterEqual(T x, T y, T eps = std::numeric_limits<T>::epsilon())
{
    return (x > y) || fpequal(x, y, eps);
}

template<typename T>
inline bool greaterThan(T x, T y, T eps = std::numeric_limits<T>::epsilon())
{
    return (x > y) && !fpequal(x, y, eps);
}

template<typename T>
inline bool lessEqual(T x, T y, T eps = std::numeric_limits<T>::epsilon())
{
    return (x < y) || fpequal(x, y, eps);
}

template<typename T>
inline bool lessThan(T x, T y, T eps = std::numeric_limits<T>::epsilon())
{
    return (x < y) && !fpequal(x, y, eps);
}

/// Test equality of two floating point numbers, with absolute accuracy,
///   independent of the value of the arguments. This is useful in relation
///   with calculations in which all important quantities are of the order 1.
template<typename T>
inline bool fpequal_abs(T x, T y, T eps = std::numeric_limits<T>::epsilon())
{
    return std::fabs(x - y) <= eps;
}

template<typename T>
inline bool greaterEqual_abs(T x, T y, T eps = std::numeric_limits<T>::epsilon())
{
    return (x > y) || fpequal_abs(x, y, eps);
}

template<typename T>
inline bool greaterThan_abs(T x, T y, T eps = std::numeric_limits<T>::epsilon())
{
    return (x > y) && !fpequal_abs(x, y, eps);
}

template<typename T>
inline bool lessEqual_abs(T x, T y, T eps = std::numeric_limits<T>::epsilon())
{
    return (x < y) || fpequal_abs(x, y, eps);
}

template<typename T>
inline bool lessThan_abs(T x, T y, T eps = std::numeric_limits<T>::epsilon())
{
    return (x < y) && !fpequal_abs(x, y, eps);
}

template<typename T>
inline bool isZero(T x, T eps = std::numeric_limits<T>::epsilon())
{
    return fpequal(x, (T) 0, eps);
}

template<typename T>
inline bool isOne(T x, T eps = std::numeric_limits<T>::epsilon())
{
    return fpequal(x, (T) 1, eps);
}

template<typename T>
inline bool isNaN(T x)
{
    // Warning: This check might be hardware dependent.
    return x != x;
}

template<typename T>
inline bool isInfinite(T x)
{
    return !(std::fabs(x) <= std::numeric_limits<T>::max());
}

template<typename T>
inline bool isFiniteNumber(T x)
{
    return (!isNaN(x) && !isInfinite(x));
}


class UniqueId {
public:
    UniqueId() : currentId(0) { }
    id_t getId() {
        if (currentId==std::numeric_limits<id_t>::max()) {
            throw PlbLogicException("Too many unique IDs requested.");
        }
        assignedIds.insert(currentId);
        return currentId++;
    }
    void releaseId(id_t id) {
        std::set<id_t>::iterator it = assignedIds.find(id);
        if (it==assignedIds.end()) {
            throw PlbLogicException("Releasing a non-assigned ID.");
        }
        assignedIds.erase(it);
    }
private:
    UniqueId(UniqueId const& rhs) { }
    UniqueId& operator=(UniqueId const& rhs) { return *this; }
private:
    id_t currentId;
    std::set<id_t> assignedIds;
};

}  // namespace util

}  // namespace plb

#endif  // UTIL_H

