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
 * Class Array, a fixed-width vector type for general use in Palabos.
 */

#ifndef ARRAY_H
#define ARRAY_H

#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include <algorithm>

namespace plb {

/** A simple array class, which is slightly more convenient
 *  to use than the pure C-array. It can be assigned from an
 *  Array of different type, contains bound verifications for
 *  debugging, and is specialized for the 2D and 3D case.
 *  
 *  ATTENTION: Values are not reset to zero in the constructor,
 *  for efficiency reason. If you need them to default to
 *  zero, you need to invoke method resetToZero() explicitly.
 */
template<typename T, pluint size>
class Array {
public:
    Array() { }
    /// Copy-construction is only allowed with same-size Arrays,
    /// to prevent bugs due to cutting off one of the arrays.
    template<typename U>
    Array(Array<U,size> const& rhs) {
        std::copy(&rhs[0], &rhs[0]+size, data);
    }
    /// Assignment is only allowed with same-size Arrays,
    /// to prevent bugs due to cutting off one of the arrays.
    template<typename U>
    Array<T,size>& operator=(Array<U,size> const& rhs) {
        std::copy(&rhs[0], &rhs[0]+size, data);
        return *this;
    }
    T& operator[](pluint index) {
        PLB_PRECONDITION(index<size);
        return data[index];
    }
    T const& operator[](pluint index) const {
        PLB_PRECONDITION(index<size);
        return data[index];
    }
    void from_cArray(T const* cArray) {
        std::copy(cArray, cArray+size, data);
    }
    void add_from_cArray(T const* cArray) {
        for (pluint i=0; i<size; ++i) {
            data[i] += cArray[i];
        }
    }
    void add_from_cArray(T const* cArray, T factor) {
        for (pluint i=0; i<size; ++i) {
            data[i] += factor*cArray[i];
        }
    }
    void to_cArray(T* cArray) const {
        std::copy(data, data+size, cArray);
    }
    void add_to_cArray(T* cArray) const {
        for (pluint i=0; i<size; ++i) {
            cArray[i] += data[i];
        }
    }
    void add_to_cArray(T* cArray, T factor) const {
        for (pluint i=0; i<size; ++i) {
            cArray[i] += factor*data[i];
        }
    }
    void resetToZero() {
        std::fill_n(data, size, T());
    }
    Array<T,size>& operator += (Array<T,size> const& b) {
        for (pluint i=0; i<size; ++i) {
            data[i] += b[i];
        }
        return *this;
    }
    Array<T,size>& operator += (T alpha) {
        for (pluint i=0; i<size; ++i) {
            data[i] += alpha;
        }
        return *this;
    }
    Array<T,size>& operator -= (Array<T,size> const& b) {
        for (pluint i=0; i<size; ++i) {
            data[i] -= b[i];
        }
        return *this;
    }
    Array<T,size>& operator -= (T alpha) {
        for (pluint i=0; i<size; ++i) {
            data[i] -= alpha;
        }
        return *this;
    }
    Array<T,size>& operator *= (Array<T,size> const& b) {
        for (pluint i=0; i<size; ++i) {
            data[i] *= b[i];
        }
        return *this;
    }
    Array<T,size>& operator *= (T alpha) {
        for (pluint i=0; i<size; ++i) {
            data[i] *= alpha;
        }
        return *this;
    }
    Array<T,size>& operator /= (Array<T,size> const& b) {
        for (pluint i=0; i<size; ++i) {
            data[i] /= b[i];
        }
        return *this;
    }
    Array<T,size>& operator /= (T alpha) {
        for (pluint i=0; i<size; ++i) {
            data[i] /= alpha;
        }
        return *this;
    }
public:
    static Array<T,size> zero() {
        Array<T,size> result;
        for (pluint i=0; i<size; ++i) {
            result[i] = T();
        }
        return result;
    }
private:
    T data[size];
};

template <typename T, pluint size>
Array<T,size> operator+(Array<T,size> const& a, Array<T,size> const& b) {
    Array<T,size> result;
    for (pluint i=0; i<size; ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

template<typename T, pluint size>
Array<T,size> operator+(Array<T,size> const& a, T alpha) {
    Array<T,size> result;
    for (pluint i=0; i<size; ++i) {
        result[i] = a[i] + alpha;
    }
    return result;
}

template<typename T, pluint size>
Array<T,size> operator+(T alpha, Array<T,size> const& a) {
    Array<T,size> result;
    for (pluint i=0; i<size; ++i) {
        result[i] = alpha + a[i];
    }
    return result;
}

template<typename T, pluint size>
Array<T,size> operator-(Array<T,size> const& a, Array<T,size> const& b) {
    Array<T,size> result;
    for (pluint i=0; i<size; ++i) {
        result[i] = a[i] - b[i];
    }
    return result;
}

template<typename T, pluint size>
Array<T,size> operator-(Array<T,size> const& a) {
    Array<T,size> result;
    for (pluint i=0; i<size; ++i) {
        result[i] = -a[i];
    }
    return result;
}

template<typename T, pluint size>
Array<T,size> operator-(Array<T,size> const& a, T alpha) {
    Array<T,size> result;
    for (pluint i=0; i<size; ++i) {
        result[i] = a[i] - alpha;
    }
    return result;
}

template<typename T, pluint size>
Array<T,size> operator-(T alpha, Array<T,size> const& a) {
    Array<T,size> result;
    for (pluint i=0; i<size; ++i) {
        result[i] = alpha - a[i];
    }
    return result;
}

template<typename T, pluint size>
Array<T,size> operator*(Array<T,size> const& a, Array<T,size> const& b) {
    Array<T,size> result;
    for (pluint i=0; i<size; ++i) {
        result[i] = a[i] * b[i];
    }
    return result;
}

template<typename T, pluint size>
Array<T,size> operator*(Array<T,size> const& a, T alpha) {
    Array<T,size> result;
    for (pluint i=0; i<size; ++i) {
        result[i] = a[i] * alpha;
    }
    return result;
}

template<typename T, pluint size>
Array<T,size> operator*(T alpha, Array<T,size> const& a) {
    Array<T,size> result;
    for (pluint i=0; i<size; ++i) {
        result[i] = alpha * a[i];
    }
    return result;
}

template<typename T, pluint size>
Array<T,size> operator/(Array<T,size> const& a, Array<T,size> const& b) {
    Array<T,size> result;
    for (pluint i=0; i<size; ++i) {
        result[i] = a[i] / b[i];
    }
    return result;
}

template<typename T, pluint size>
Array<T,size> operator/(Array<T,size> const& a, T alpha) {
    Array<T,size> result;
    for (pluint i=0; i<size; ++i) {
        result[i] = a[i] / alpha;
    }
    return result;
}

template<typename T, pluint size>
inline bool operator==(Array<T,size> const& A, Array<T,size> const& B) {
    for (pluint iA = 0; iA < size; ++iA) {
        if (A[iA] != B[iA]) return false;
    }
    return true;
}

template<typename T, pluint size>
inline bool operator!=(Array<T,size> const& A, Array<T,size> const& B) {
    for (pluint iA = 0; iA < size; ++iA) {
        if (A[iA] != B[iA]) return true;
    }
    return false;
}

template<typename T, pluint size>
inline bool operator<(Array<T,size> const& A, Array<T,size> const& B) {
    for (pluint iA=0; iA<size; ++iA) {
        if (A[iA] < B[iA]) return true;
        if (A[iA] > B[iA]) return false;
    }
    return false;
}

template<typename T, pluint size>
Array<T,size> operator/(T alpha, Array<T,size> const& a) {
    Array<T,size> result;
    for (pluint i=0; i<size; ++i) {
        result[i] = alpha / a[i];
    }
    return result;
}

template<typename T, pluint size>
T maxElement(Array<T,size> const& a) {
    PLB_ASSERT(size > 1);
    T maxEl = std::max(a[0],a[1]);

    for (pluint i=2; i<size; ++i) {
        maxEl = std::max(a[i],maxEl);
    }
    return maxEl;
}


template<typename T>
class Array<T,2> {
public:
    Array() { }
    Array(T x, T y) {
        data[0] = x;
        data[1] = y;
    }
    template<typename U>
    Array(Array<U,2> const& rhs) {
        data[0] = (T) rhs[0];
        data[1] = (T) rhs[1];
    }
    template<typename U>
    Array<T,2>& operator=(Array<U,2> const& rhs) {
        data[0] = (T) rhs[0];
        data[1] = (T) rhs[1];
        return *this;
    }
    T& operator[](pluint index) {
        PLB_PRECONDITION(index<2);
        return data[index];
    }
    T const& operator[](pluint index) const {
        PLB_PRECONDITION(index<2);
        return data[index];
    }
    void from_cArray(T const* cArray) {
        data[0] = cArray[0];
        data[1] = cArray[1];
    }
    void add_from_cArray(T const* cArray) {
        data[0] += cArray[0];
        data[1] += cArray[1];
    }
    void add_from_cArray(T const* cArray, T factor) {
        data[0] += factor*cArray[0];
        data[1] += factor*cArray[1];
    }
    void to_cArray(T* cArray) const {
        cArray[0] = data[0];
        cArray[1] = data[1];
    }
    void add_to_cArray(T* cArray) const {
        cArray[0] += data[0];
        cArray[1] += data[1];
    }
    void add_to_cArray(T* cArray, T factor) const {
        cArray[0] += factor*data[0];
        cArray[1] += factor*data[1];
    }
    void resetToZero() {
        data[0] = T();
        data[1] = T();
    }
    Array<T,2>& operator+=(Array<T,2> const& b) {
        data[0] += b[0];
        data[1] += b[1];
        return *this;
    }
    Array<T,2>& operator+=(T alpha) {
        data[0] += alpha;
        data[1] += alpha;
        return *this;
    }
    Array<T,2>& operator-=(Array<T,2> const& b) {
        data[0] -= b[0];
        data[1] -= b[1];
        return *this;
    }
    Array<T,2>& operator-=(T alpha) {
        data[0] -= alpha;
        data[1] -= alpha;
        return *this;
    }
    Array<T,2>& operator*=(Array<T,2> const& b) {
        data[0] *= b[0];
        data[1] *= b[1];
        return *this;
    }
    Array<T,2>& operator*=(T alpha) {
        data[0] *= alpha;
        data[1] *= alpha;
        return *this;
    }
    Array<T,2>& operator/=(Array<T,2> const& b) {
        data[0] /= b[0];
        data[1] /= b[1];
        return *this;
    }
    Array<T,2>& operator/=(T alpha) {
        data[0] /= alpha;
        data[1] /= alpha;
        return *this;
    }
    static Array<T,2> zero() {
        return Array<T,2>(T(), T());
    }
private:
    T data[2];
};

template<typename T>
Array<T,2> operator+(Array<T,2> const& a, Array<T,2> const& b) {
    return Array<T,2>(a[0]+b[0], a[1]+b[1]);
}

template<typename T>
Array<T,2> operator+(Array<T,2> const& a, T alpha) {
    return Array<T,2>(a[0]+alpha, a[1]+alpha);
}

template<typename T>
Array<T,2> operator+(T alpha, Array<T,2> const& a) {
    return Array<T,2>(alpha+a[0], alpha+a[1]);
}

template<typename T>
Array<T,2> operator-(Array<T,2> const& a, Array<T,2> const& b) {
    return Array<T,2>(a[0]-b[0], a[1]-b[1]);
}

template<typename T>
Array<T,2> operator-(Array<T,2> const& a) {
    return Array<T,2>(-a[0],-a[1]);
}

template<typename T>
Array<T,2> operator-(Array<T,2> const& a, T alpha) {
    return Array<T,2>(a[0]-alpha, a[1]-alpha);
}

template<typename T>
Array<T,2> operator-(T alpha, Array<T,2> const& a) {
    return Array<T,2>(alpha-a[0], alpha-a[1]);
}

template<typename T>
Array<T,2> operator*(Array<T,2> const& a, Array<T,2> const& b) {
    return Array<T,2>(a[0]*b[0], a[1]*b[1]);
}

template<typename T>
Array<T,2> operator*(Array<T,2> const& a, T alpha) {
    return Array<T,2>(a[0]*alpha, a[1]*alpha);
}

template<typename T>
Array<T,2> operator*(T alpha, Array<T,2> const& a) {
    return Array<T,2>(alpha*a[0], alpha*a[1]);
}

template<typename T>
Array<T,2> operator/(Array<T,2> const& a, Array<T,2> const& b) {
    return Array<T,2>(a[0]/b[0], a[1]/b[1]);
}

template<typename T>
Array<T,2> operator/(Array<T,2> const& a, T alpha) {
    return Array<T,2>(a[0]/alpha, a[1]/alpha);
}

template<typename T>
Array<T,2> operator/(T alpha, Array<T,2> const& a) {
    return Array<T,2>(alpha/a[0], alpha/a[1]);
}

template<typename T>
class Array<T,3> {
public:
    Array() { }
    Array(T x, T y, T z) {
        data[0] = x;
        data[1] = y;
        data[2] = z;
    }
    template<typename U>
    Array(Array<U,3> const& rhs) {
        data[0] = (T) rhs[0];
        data[1] = (T) rhs[1];
        data[2] = (T) rhs[2];
    }
    template<typename U>
    Array<T,3>& operator=(Array<U,3> const& rhs) {
        data[0] = (T) rhs[0];
        data[1] = (T) rhs[1];
        data[2] = (T) rhs[2];
        return *this;
    }
    T& operator[](pluint index) {
        PLB_PRECONDITION(index<3);
        return data[index];
    }
    T const& operator[](pluint index) const {
        PLB_PRECONDITION(index<3);
        return data[index];
    }
    void from_cArray(T const* cArray) {
        data[0] = cArray[0];
        data[1] = cArray[1];
        data[2] = cArray[2];
    }
    void add_from_cArray(T const* cArray) {
        data[0] += cArray[0];
        data[1] += cArray[1];
        data[2] += cArray[2];
    }
    void add_from_cArray(T const* cArray, T factor) {
        data[0] += factor*cArray[0];
        data[1] += factor*cArray[1];
        data[2] += factor*cArray[2];
    }
    void to_cArray(T* cArray) const {
        cArray[0] = data[0];
        cArray[1] = data[1];
        cArray[2] = data[2];
    }
    void add_to_cArray(T* cArray) const {
        cArray[0] += data[0];
        cArray[1] += data[1];
        cArray[2] += data[2];
    }
    void add_to_cArray(T* cArray, T factor) const {
        cArray[0] += factor*data[0];
        cArray[1] += factor*data[1];
        cArray[2] += factor*data[2];
    }
    void resetToZero() {
        data[0] = T();
        data[1] = T();
        data[2] = T();
    }
    Array<T,3>& operator += (Array<T,3> const& b) {
        data[0] += b[0];
        data[1] += b[1];
        data[2] += b[2];
        return *this;
    }
    Array<T,3>& operator += (T alpha) {
        data[0] += alpha;
        data[1] += alpha;
        data[2] += alpha;
        return *this;
    }
    Array<T,3>& operator -= (Array<T,3> const& b) {
        data[0] -= b[0];
        data[1] -= b[1];
        data[2] -= b[2];
        return *this;
    }
    Array<T,3>& operator -= (T alpha) {
        data[0] -= alpha;
        data[1] -= alpha;
        data[2] -= alpha;
        return *this;
    }
    Array<T,3>& operator *= (Array<T,3> const& b) {
        data[0] *= b[0];
        data[1] *= b[1];
        data[2] *= b[2];
        return *this;
    }
    Array<T,3>& operator *= (T alpha) {
        data[0] *= alpha;
        data[1] *= alpha;
        data[2] *= alpha;
        return *this;
    }
    Array<T,3>& operator /= (Array<T,3> const& b) {
        data[0] /= b[0];
        data[1] /= b[1];
        data[2] /= b[2];
        return *this;
    }
    Array<T,3>& operator /= (T alpha) {
        data[0] /= alpha;
        data[1] /= alpha;
        data[2] /= alpha;
        return *this;
    }
    static Array<T,3> zero() {
        return Array<T,3>(T(), T(), T());
    }
private:
    T data[3];
};

template<typename T>
Array<T,3> operator+(Array<T,3> const& a, Array<T,3> const& b) {
    return Array<T,3>(a[0]+b[0], a[1]+b[1], a[2]+b[2]);
}

template<typename T>
Array<T,3> operator+(Array<T,3> const& a, T alpha) {
    return Array<T,3>(a[0]+alpha, a[1]+alpha, a[2]+alpha);
}

template<typename T>
Array<T,3> operator+(T alpha, Array<T,3> const& a) {
    return Array<T,3>(alpha+a[0], alpha+a[1], alpha+a[2]);
}

template<typename T>
Array<T,3> operator-(Array<T,3> const& a, Array<T,3> const& b) {
    return Array<T,3>(a[0]-b[0], a[1]-b[1], a[2]-b[2]);
}

template<typename T>
Array<T,3> operator-(Array<T,3> const& a) {
    return Array<T,3>(-a[0], -a[1], -a[2]);
}

template<typename T>
Array<T,3> operator-(Array<T,3> const& a, T alpha) {
    return Array<T,3>(a[0]-alpha, a[1]-alpha, a[2]-alpha);
}

template<typename T>
Array<T,3> operator-(T alpha, Array<T,3> const& a) {
    return Array<T,3>(alpha-a[0], alpha-a[1], alpha-a[2]);
}

template<typename T>
Array<T,3> operator*(Array<T,3> const& a, Array<T,3> const& b) {
    return Array<T,3>(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
}

template<typename T>
Array<T,3> operator*(Array<T,3> const& a, T alpha) {
    return Array<T,3>(a[0]*alpha, a[1]*alpha, a[2]*alpha);
}

template<typename T>
Array<T,3> operator*(T alpha, Array<T,3> const& a) {
    return Array<T,3>(alpha*a[0], alpha*a[1], alpha*a[2]);
}

template<typename T>
Array<T,3> operator/(Array<T,3> const& a, Array<T,3> const& b) {
    return Array<T,3>(a[0]/b[0], a[1]/b[1], a[2]/b[2]);
}

template<typename T>
Array<T,3> operator/(Array<T,3> const& a, T alpha) {
    return Array<T,3>(a[0]/alpha, a[1]/alpha, a[2]/alpha);
}

template<typename T>
Array<T,3> operator/(T alpha, Array<T,3> const& a) {
    return Array<T,3>(alpha/a[0], alpha/a[1], alpha/a[2]);
}


template<typename T>
class Array<T,4> {
public:
    Array() { }
    Array(T a0, T a1, T a2, T a3) {
        data[0] = a0;
        data[1] = a1;
        data[2] = a2;
        data[3] = a3;
    }
    template<typename U>
    Array(Array<U,4> const& rhs) {
        data[0] = (T) rhs[0];
        data[1] = (T) rhs[1];
        data[2] = (T) rhs[2];
        data[3] = (T) rhs[3];
    }
    template<typename U>
    Array<T,4>& operator=(Array<U,4> const& rhs) {
        data[0] = (T) rhs[0];
        data[1] = (T) rhs[1];
        data[2] = (T) rhs[2];
        data[3] = (T) rhs[3];
        return *this;
    }
    T& operator[](pluint index) {
        PLB_PRECONDITION(index<4);
        return data[index];
    }
    T const& operator[](pluint index) const {
        PLB_PRECONDITION(index<4);
        return data[index];
    }
    void from_cArray(T const* cArray) {
        data[0] = cArray[0];
        data[1] = cArray[1];
        data[2] = cArray[2];
        data[3] = cArray[3];
    }
    void add_from_cArray(T const* cArray) {
        data[0] += cArray[0];
        data[1] += cArray[1];
        data[2] += cArray[2];
        data[3] += cArray[3];
    }
    void add_from_cArray(T const* cArray, T factor) {
        data[0] += factor*cArray[0];
        data[1] += factor*cArray[1];
        data[2] += factor*cArray[2];
        data[3] += factor*cArray[3];
    }
    void to_cArray(T* cArray) const {
        cArray[0] = data[0];
        cArray[1] = data[1];
        cArray[2] = data[2];
        cArray[3] = data[3];
    }
    void add_to_cArray(T* cArray) const {
        cArray[0] += data[0];
        cArray[1] += data[1];
        cArray[2] += data[2];
        cArray[3] += data[3];
    }
    void add_to_cArray(T* cArray, T factor) const {
        cArray[0] += factor*data[0];
        cArray[1] += factor*data[1];
        cArray[2] += factor*data[2];
        cArray[3] += factor*data[3];
    }
    void resetToZero() {
        data[0] = T();
        data[1] = T();
        data[2] = T();
        data[3] = T();
    }
    Array<T,4>& operator += (Array<T,4> const& b) {
        data[0] += b[0];
        data[1] += b[1];
        data[2] += b[2];
        data[3] += b[3];
        return *this;
    }
    Array<T,4>& operator += (T alpha) {
        data[0] += alpha;
        data[1] += alpha;
        data[2] += alpha;
        data[3] += alpha;
        return *this;
    }
    Array<T,4>& operator -= (Array<T,4> const& b) {
        data[0] -= b[0];
        data[1] -= b[1];
        data[2] -= b[2];
        data[3] -= b[3];
        return *this;
    }
    Array<T,4>& operator -= (T alpha) {
        data[0] -= alpha;
        data[1] -= alpha;
        data[2] -= alpha;
        data[3] -= alpha;
        return *this;
    }
    Array<T,4>& operator *= (Array<T,4> const& b) {
        data[0] *= b[0];
        data[1] *= b[1];
        data[2] *= b[2];
        data[3] *= b[3];
        return *this;
    }
    Array<T,4>& operator *= (T alpha) {
        data[0] *= alpha;
        data[1] *= alpha;
        data[2] *= alpha;
        data[3] *= alpha;
        return *this;
    }
    Array<T,4>& operator /= (Array<T,4> const& b) {
        data[0] /= b[0];
        data[1] /= b[1];
        data[2] /= b[2];
        data[3] /= b[3];
        return *this;
    }
    Array<T,4>& operator /= (T alpha) {
        data[0] /= alpha;
        data[1] /= alpha;
        data[2] /= alpha;
        data[3] /= alpha;
        return *this;
    }
    static Array<T,4> zero() {
        return Array<T,3>(T(), T(), T(), T());
    }
private:
    T data[4];
};

template<typename T>
Array<T,4> operator+(Array<T,4> const& a, Array<T,4> const& b) {
    return Array<T,4>(a[0]+b[0], a[1]+b[1], a[2]+b[2], a[3]+b[3]);
}

template<typename T>
Array<T,4> operator+(Array<T,4> const& a, T alpha) {
    return Array<T,4>(a[0]+alpha, a[1]+alpha, a[2]+alpha, a[3]+alpha);
}

template<typename T>
Array<T,4> operator+(T alpha, Array<T,4> const& a) {
    return Array<T,4>(alpha+a[0], alpha+a[1], alpha+a[2], alpha+a[3]);
}

template<typename T>
Array<T,4> operator-(Array<T,4> const& a, Array<T,4> const& b) {
    return Array<T,4>(a[0]-b[0], a[1]-b[1], a[2]-b[2], a[3]-b[3]);
}

template<typename T>
Array<T,4> operator-(Array<T,4> const& a) {
    return Array<T,4>(-a[0], -a[1], -a[2], -a[3]);
}

template<typename T>
Array<T,4> operator-(Array<T,4> const& a, T alpha) {
    return Array<T,4>(a[0]-alpha, a[1]-alpha, a[2]-alpha, a[3]-alpha);
}

template<typename T>
Array<T,4> operator-(T alpha, Array<T,4> const& a) {
    return Array<T,4>(alpha-a[0], alpha-a[1], alpha-a[2], alpha-a[3]);
}

template<typename T>
Array<T,4> operator*(Array<T,4> const& a, Array<T,4> const& b) {
    return Array<T,4>(a[0]*b[0], a[1]*b[1], a[2]*b[2], a[3]*b[3]);
}

template<typename T>
Array<T,4> operator*(Array<T,4> const& a, T alpha) {
    return Array<T,4>(a[0]*alpha, a[1]*alpha, a[2]*alpha, a[3]*alpha);
}

template<typename T>
Array<T,4> operator*(T alpha, Array<T,4> const& a) {
    return Array<T,4>(alpha*a[0], alpha*a[1], alpha*a[2], alpha*a[3]);
}

template<typename T>
Array<T,4> operator/(Array<T,4> const& a, Array<T,4> const& b) {
    return Array<T,4>(a[0]/b[0], a[1]/b[1], a[2]/b[2], a[3]/b[3]);
}

template<typename T>
Array<T,4> operator/(Array<T,4> const& a, T alpha) {
    return Array<T,4>(a[0]/alpha, a[1]/alpha, a[2]/alpha, a[3]/alpha);
}

template<typename T>
Array<T,4> operator/(T alpha, Array<T,4> const& a) {
    return Array<T,4>(alpha/a[0], alpha/a[1], alpha/a[2], alpha/a[3]);
}


template<typename T>
class Array<T,6> {
public:
    Array() { }
    Array(T a0, T a1, T a2, T a3, T a4, T a5) {
        data[0] = a0;
        data[1] = a1;
        data[2] = a2;
        data[3] = a3;
        data[4] = a4;
        data[5] = a5;
    }
    template<typename U>
    Array(Array<U,6> const& rhs) {
        data[0] = (T) rhs[0];
        data[1] = (T) rhs[1];
        data[2] = (T) rhs[2];
        data[3] = (T) rhs[3];
        data[4] = (T) rhs[4];
        data[5] = (T) rhs[5];
    }
    template<typename U>
    Array<T,6>& operator=(Array<U,6> const& rhs) {
        data[0] = (T) rhs[0];
        data[1] = (T) rhs[1];
        data[2] = (T) rhs[2];
        data[3] = (T) rhs[3];
        data[4] = (T) rhs[4];
        data[5] = (T) rhs[5];
        return *this;
    }
    T& operator[](pluint index) {
        PLB_PRECONDITION(index<6);
        return data[index];
    }
    T const& operator[](pluint index) const {
        PLB_PRECONDITION(index<6);
        return data[index];
    }
    void from_cArray(T const* cArray) {
        data[0] = cArray[0];
        data[1] = cArray[1];
        data[2] = cArray[2];
        data[3] = cArray[3];
        data[4] = cArray[4];
        data[5] = cArray[5];
    }
    void add_from_cArray(T const* cArray) {
        data[0] += cArray[0];
        data[1] += cArray[1];
        data[2] += cArray[2];
        data[3] += cArray[3];
        data[4] += cArray[4];
        data[5] += cArray[5];
    }
    void add_from_cArray(T const* cArray, T factor) {
        data[0] += factor*cArray[0];
        data[1] += factor*cArray[1];
        data[2] += factor*cArray[2];
        data[3] += factor*cArray[3];
        data[4] += factor*cArray[4];
        data[5] += factor*cArray[5];
    }
    void to_cArray(T* cArray) const {
        cArray[0] = data[0];
        cArray[1] = data[1];
        cArray[2] = data[2];
        cArray[3] = data[3];
        cArray[4] = data[4];
        cArray[5] = data[5];
    }
    void add_to_cArray(T* cArray) const {
        cArray[0] += data[0];
        cArray[1] += data[1];
        cArray[2] += data[2];
        cArray[3] += data[3];
        cArray[4] += data[4];
        cArray[5] += data[5];
    }
    void add_to_cArray(T* cArray, T factor) const {
        cArray[0] += factor*data[0];
        cArray[1] += factor*data[1];
        cArray[2] += factor*data[2];
        cArray[3] += factor*data[3];
        cArray[4] += factor*data[4];
        cArray[5] += factor*data[5];
    }
    void resetToZero() {
        data[0] = T();
        data[1] = T();
        data[2] = T();
        data[3] = T();
        data[4] = T();
        data[5] = T();
    }
    Array<T,6>& operator += (Array<T,6> const& b) {
        data[0] += b[0];
        data[1] += b[1];
        data[2] += b[2];
        data[3] += b[3];
        data[4] += b[4];
        data[5] += b[5];
        return *this;
    }
    Array<T,6>& operator += (T alpha) {
        data[0] += alpha;
        data[1] += alpha;
        data[2] += alpha;
        data[3] += alpha;
        data[4] += alpha;
        data[5] += alpha;
        return *this;
    }
    Array<T,6>& operator -= (Array<T,6> const& b) {
        data[0] -= b[0];
        data[1] -= b[1];
        data[2] -= b[2];
        data[3] -= b[3];
        data[4] -= b[4];
        data[5] -= b[5];
        return *this;
    }
    Array<T,6>& operator -= (T alpha) {
        data[0] -= alpha;
        data[1] -= alpha;
        data[2] -= alpha;
        data[3] -= alpha;
        data[4] -= alpha;
        data[5] -= alpha;
        return *this;
    }
    Array<T,6>& operator *= (Array<T,6> const& b) {
        data[0] *= b[0];
        data[1] *= b[1];
        data[2] *= b[2];
        data[3] *= b[3];
        data[4] *= b[4];
        data[5] *= b[5];
        return *this;
    }
    Array<T,6>& operator *= (T alpha) {
        data[0] *= alpha;
        data[1] *= alpha;
        data[2] *= alpha;
        data[3] *= alpha;
        data[4] *= alpha;
        data[5] *= alpha;
        return *this;
    }
    Array<T,6>& operator /= (Array<T,6> const& b) {
        data[0] /= b[0];
        data[1] /= b[1];
        data[2] /= b[2];
        data[3] /= b[3];
        data[4] /= b[4];
        data[5] /= b[5];
        return *this;
    }
    Array<T,6>& operator /= (T alpha) {
        data[0] /= alpha;
        data[1] /= alpha;
        data[2] /= alpha;
        data[3] /= alpha;
        data[4] /= alpha;
        data[5] /= alpha;
        return *this;
    }
    static Array<T,6> zero() {
        return Array<T,6>(T(), T(), T(), T(), T(), T());
    }
private:
    T data[6];
};

template<typename T>
Array<T,6> operator+(Array<T,6> const& a, Array<T,6> const& b) {
    return Array<T,6>(a[0]+b[0], a[1]+b[1], a[2]+b[2], a[3]+b[3], a[4]+b[4], a[5]+b[5]);
}

template<typename T>
Array<T,6> operator+(Array<T,6> const& a, T alpha) {
    return Array<T,6>(a[0]+alpha, a[1]+alpha, a[2]+alpha, a[3]+alpha, a[4]+alpha, a[5]+alpha);
}

template<typename T>
Array<T,6> operator+(T alpha, Array<T,6> const& a) {
    return Array<T,6>(alpha+a[0], alpha+a[1], alpha+a[2], alpha+a[3], alpha+a[4], alpha+a[5]);
}

template<typename T>
Array<T,6> operator-(Array<T,6> const& a, Array<T,6> const& b) {
    return Array<T,6>(a[0]-b[0], a[1]-b[1], a[2]-b[2], a[3]-b[3], a[4]-b[4], a[5]-b[5]);
}

template<typename T>
Array<T,6> operator-(Array<T,6> const& a) {
    return Array<T,6>(-a[0], -a[1], -a[2], -a[3], -a[4], -a[5]);
}

template<typename T>
Array<T,6> operator-(Array<T,6> const& a, T alpha) {
    return Array<T,6>(a[0]-alpha, a[1]-alpha, a[2]-alpha, a[3]-alpha, a[4]-alpha, a[5]-alpha);
}

template<typename T>
Array<T,6> operator-(T alpha, Array<T,6> const& a) {
    return Array<T,6>(alpha-a[0], alpha-a[1], alpha-a[2], alpha-a[3], alpha-a[4], alpha-a[5]);
}

template<typename T>
Array<T,6> operator*(Array<T,6> const& a, Array<T,6> const& b) {
    return Array<T,6>(a[0]*b[0], a[1]*b[1], a[2]*b[2], a[3]*b[3], a[4]*b[4], a[5]*b[5]);
}

template<typename T>
Array<T,6> operator*(Array<T,6> const& a, T alpha) {
    return Array<T,6>(a[0]*alpha, a[1]*alpha, a[2]*alpha, a[3]*alpha, a[4]*alpha, a[5]*alpha);
}

template<typename T>
Array<T,6> operator*(T alpha, Array<T,6> const& a) {
    return Array<T,6>(alpha*a[0], alpha*a[1], alpha*a[2], alpha*a[3], alpha*a[4], alpha*a[5]);
}

template<typename T>
Array<T,6> operator/(Array<T,6> const& a, Array<T,6> const& b) {
    return Array<T,6>(a[0]/b[0], a[1]/b[1], a[2]/b[2], a[3]/b[3], a[4]/b[4], a[5]/b[5]);
}

template<typename T>
Array<T,6> operator/(Array<T,6> const& a, T alpha) {
    return Array<T,6>(a[0]/alpha, a[1]/alpha, a[2]/alpha, a[3]/alpha, a[4]/alpha, a[5]/alpha);
}

template<typename T>
Array<T,6> operator/(T alpha, Array<T,6> const& a) {
    return Array<T,6>(alpha/a[0], alpha/a[1], alpha/a[2], alpha/a[3], alpha/a[4], alpha/a[5]);
}

} // end namespace plb

#endif

