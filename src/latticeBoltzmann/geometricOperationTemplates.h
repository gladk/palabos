/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2015 FlowKit Sarl
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
 * Templates for common geometric operations (scalar product, vector-
 * matric operations etc.).
 *  -- header file
 */
#ifndef GEOMETRIC_OPERATION_TEMPLATES_H
#define GEOMETRIC_OPERATION_TEMPLATES_H

#include "io/parallelIO.h"
#include "core/globalDefs.h"
#include "core/array.h"
#include "core/util.h"
#include <algorithm>
#include <vector>
#include <cmath>
#include <limits>

namespace plb {

template <typename T, int d> struct VectorTemplateImpl;
template <typename T, template<typename U> class Descriptor> struct SymmetricTensor;
template <typename T, template<typename U> class Descriptor> struct SymmetricRankThreeTensor;
template <typename T, int d> struct SymmetricTensorImpl;
template <typename T, int d> struct SymmetricRankThreeTensorImpl;

template <typename T, template<typename U> class Descriptor>
struct VectorTemplate {
    /// Number of dimensions for current lattice
    static const int d = Descriptor<T>::d;
    /// Compute scalar product between two vectors
    static T scalarProduct(Array<T,d> const& u1, Array<T,d> const& u2) {
        return VectorTemplateImpl<T,d>::scalarProduct(u1,u2);
    }
    /// Compute scalar product between two a c-array and a plb-array
    static T scalarProduct(const T u1[d], Array<T,d> const& u2) {
        return VectorTemplateImpl<T,d>::scalarProduct(u1,u2);
    }
    /// Compute norm-square of a vector
    static T normSqr(Array<T,d> const& u) {
        return VectorTemplateImpl<T,d>::normSqr(u);
    }
    /// Multiply vector elements component-wise by a scalar
    static void multiplyByScalar(Array<T,d>& u, T scalar) {
        VectorTemplateImpl<T,d>::multiplyByScalar(u,scalar);
    }
    /// Multiply vector elements component-wise by a scalar and store in second vector
    static void multiplyByScalar(Array<T,d> const& u, T scalar, Array<T,d>& result) {
        VectorTemplateImpl<T,d>::multiplyByScalar(u,scalar,result);
    }
    /// Symmetric tensor product of two vectors
    static void symTensorProduct( Array<T,d> const &u, Array<T,d> const &v,
                                  Array<T,SymmetricTensor<T,Descriptor>::n> &result)
    {
        VectorTemplateImpl<T,d>::symTensorProduct(u,v,result);
    }
    /// Symmetric tensor product of three vectors
    static void symTensorProduct(Array<T,d> const &u, Array<T,d> const &v, Array<T,d> const &w,
                                 Array<T,SymmetricRankThreeTensor<T,Descriptor>::n> &result)
    {
        VectorTemplateImpl<T,d>::symTensorProduct(u,v,w,result);
    }
};

template <typename T, int d>
struct VectorTemplateImpl {
    static T scalarProduct(Array<T,d> const& u1, Array<T,d> const& u2) {
        T result = T();
        for (int iD=0; iD<d; ++iD) {
            result += u1[iD]*u2[iD];
        }
        return result;
    }
    static T scalarProduct(const T u1[d], Array<T,d> const& u2) {
        T result = T();
        for (int iD=0; iD<d; ++iD) {
            result += u1[iD]*u2[iD];
        }
        return result;
    }
    static T normSqr(Array<T,d> const& u) {
        return scalarProduct(u,u);
    }
    static void multiplyByScalar(Array<T,d>& u, T scalar) {
        for (int iD=0; iD<d; ++iD) {
            u[iD] *= scalar;
        }
    }
    
    static void multiplyByScalar(Array<T,d> const& u, T scalar, Array<T,d>& result) {
        for (int iD=0; iD<d; ++iD) {
            result[iD] = u[iD]*scalar;
        }
    }
    
    static void symTensorProduct(Array<T,d> const &u, Array<T,d> const &v, Array<T,SymmetricTensorImpl<T,d>::n> &result) {
        int iPi = 0;
        for (int iA = 0; iA < d; ++iA) {
            for (int iB = iA; iB < d; ++iB) {
                result[iPi] = u[iA]*v[iB];
                ++iPi;
            }
        }
    }
    
    static void symTensorProduct(Array<T,d> const &u, Array<T,d> const &v, Array<T,d> const &w, 
                                 Array<T,SymmetricRankThreeTensorImpl<T,d>::n> &result) {
        int iPi = 0;
        for (int iA = 0; iA < d; ++iA) {
            for (int iB = iA; iB < d; ++iB) {
                for (int iC = iB; iC < d; ++iC) {
                    result[iPi] = u[iA]*v[iB]*w[iC];
                    ++iPi;
                }
            }
        }
    }
};


template <typename T>
struct VectorTemplateImpl<T,2> {
    static T scalarProduct(Array<T,2> const& u1, Array<T,2> const& u2) {
        return u1[0]*u2[0] + u1[1]*u2[1];
    }
    static T scalarProduct(const T u1[2], Array<T,2> const& u2) {
        return u1[0]*u2[0] + u1[1]*u2[1];
    }
    static T normSqr(Array<T,2> const& u) {
        return u[0]*u[0] + u[1]*u[1];
    }
    static void multiplyByScalar(Array<T,2>& u, T scalar) {
        u[0] *= scalar;
        u[1] *= scalar;
    }
    static void multiplyByScalar(Array<T,2> const& u, T scalar, Array<T,2>& result) {
        result[0] = u[0]*scalar;
        result[1] = u[1]*scalar;
    }
    static void symTensorProduct(Array<T,2> const &u, Array<T,2> const &v, Array<T,3> &result) {
        result[0] = u[0]*v[0];
        result[1] = u[0]*v[1];
        result[2] = u[1]*v[1];
    }
    static void symTensorProduct(Array<T,2> const &u, Array<T,2> const &v, Array<T,2> const &w, 
                                 Array<T,4> &result) {
        result[0] = u[0]*v[0]*w[0];
        result[1] = u[0]*v[0]*w[1];
        result[2] = u[0]*v[1]*w[1];
        result[3] = u[1]*v[1]*w[1];
    }
};


template <typename T>
struct VectorTemplateImpl<T,3> {
    static T scalarProduct(Array<T,3> const& u1, Array<T,3> const& u2) {
        return u1[0]*u2[0] + u1[1]*u2[1] + u1[2]*u2[2];
    }
    static T scalarProduct(const T u1[3], Array<T,3> const& u2) {
        return u1[0]*u2[0] + u1[1]*u2[1] + u1[2]*u2[2];
    }
    static T normSqr(Array<T,3> const& u) {
        return u[0]*u[0] + u[1]*u[1] + u[2]*u[2];
    }
    static void multiplyByScalar(Array<T,3>& u, T scalar) {
        u[0] *= scalar;
        u[1] *= scalar;
        u[2] *= scalar;
    }
    static void multiplyByScalar(Array<T,3> const& u, T scalar, Array<T,3>& result) {
        result[0] = u[0]*scalar;
        result[1] = u[1]*scalar;
        result[2] = u[2]*scalar;
    }
    static void symTensorProduct(Array<T,3> const &u, Array<T,3> const &v, Array<T,6> &result) {
        result[0] = u[0]*v[0];
        result[1] = u[0]*v[1];
        result[2] = u[0]*v[2];
        result[3] = u[1]*v[1];
        result[4] = u[1]*v[2];
        result[5] = u[2]*v[2];
    }
    static void symTensorProduct(Array<T,3> const &u, Array<T,3> const &v, Array<T,3> const &w, 
                                 Array<T,10> &result) {
        result[0] = u[0]*v[0]*w[0];
        result[1] = u[0]*v[0]*w[1];
        result[2] = u[0]*v[0]*w[2];
        result[3] = u[0]*v[1]*w[1];
        result[4] = u[0]*v[1]*w[2];
        result[5] = u[0]*v[2]*w[2];
        result[6] = u[1]*v[1]*w[1];
        result[7] = u[1]*v[1]*w[2];
        result[8] = u[1]*v[2]*w[2];
        result[9] = u[2]*v[2]*w[2];
    }
};


template <typename T, int d> struct SymmetricTensorImpl { };

template<typename T> struct SymmetricTensorImpl<T,2> {
    static const int n = 3;
    enum Indices { xx=0, xy=1, yy=2 };
    static void matVectMult(Array<T,n> const& mat, Array<T,2> const& vect, Array<T,2>& result) {
        result[0] = mat[xx]*vect[0] + mat[xy]*vect[1];
        result[1] = mat[xy]*vect[0] + mat[yy]*vect[1];
    }
    static T tensorNormSqr(Array<T,n> const& mat) {
        return mat[xx]*mat[xx]+mat[yy]*mat[yy] + (T)2*mat[xy]*mat[xy];
    }
    static T contractIndexes(Array<T,n> const& A, Array<T,n> const& B) {
        return A[xx]*B[xx]+A[yy]*B[yy]+ (T)2*A[xy]*B[xy];
    }
    static T trace(Array<T,n> const &A) {
        return A[xx] + A[yy];
    }
    static Array<T,n> id() {
        return Array<T,n>((T)1,T(),(T)1);
    }
};

template<typename T> struct SymmetricTensorImpl<T,3> {
    static const int n = 6;
    enum Indices { xx=0, xy=1, xz=2, yy=3, yz=4, zz=5 };
    static void matVectMult(Array<T,n> const& mat, Array<T,3> const& vect, Array<T,3>& result) {
        result[0] = mat[xx]*vect[0] + mat[xy]*vect[1] + mat[xz]*vect[2];
        result[1] = mat[xy]*vect[0] + mat[yy]*vect[1] + mat[yz]*vect[2];
        result[2] = mat[xz]*vect[0] + mat[yz]*vect[1] + mat[zz]*vect[2];
    }
    static T tensorNormSqr(Array<T,n> const& mat) {
        return mat[xx]*mat[xx]+mat[yy]*mat[yy]+mat[zz]*mat[zz]
               + (T)2*mat[xy]*mat[xy] + (T)2*mat[xz]*mat[xz] + (T)2*mat[yz]*mat[yz];
    }
    static T contractIndexes(Array<T,n> const& A, Array<T,n> const& B) {
        return A[xx]*B[xx]+A[yy]*B[yy]+A[zz]*B[zz]
            + (T)2*(A[xy]*B[xy] + A[xz]*B[xz] + A[yz]*B[yz]);
    }
    static T trace(Array<T,n> const &A) {
        return A[xx] + A[yy] + A[zz];
    }
    static Array<T,n> id() {
        Array<T,n> I;
        
        I[0] = (T)1;
        I[1] = T();
        I[2] = T();
        I[3] = (T)1;
        I[4] = T();
        I[5] = (T)1;
        
        return I;
    }
};

/// Operations on a symmetric tensor which stores only above-or-on-diagonal values
template <typename T, template<typename U> class Descriptor>
struct SymmetricTensor {
    /// Number of dimensions for current lattice
    static const int d = Descriptor<T>::d;
    /// Number of elements (reduced by symmetry)
    static const int n = SymmetricTensorImpl<T,d>::n;
    static void matVectMult(Array<T,n> const& mat, Array<T,d> const& vect, Array<T,d>& result) {
        SymmetricTensorImpl<T,d>::matVectMult(mat, vect, result);
    }
    static T tensorNormSqr(Array<T,n> const& mat) {
        return SymmetricTensorImpl<T,d>::tensorNormSqr(mat);
    }
    static T contractIndexes(Array<T,n> const& A, Array<T,n> const& B) {
        return SymmetricTensorImpl<T,d>::contractIndexes(A,B);
    }
    static T trace(Array<T,n> const &A) {
        return SymmetricTensorImpl<T,d>::trace(A);
    }
    static Array<T,n> id() {
        return SymmetricTensorImpl<T,d>::id();
    }
};

template <typename T, int d> struct SymmetricRankThreeTensorImpl { };

template<typename T> struct SymmetricRankThreeTensorImpl<T,2> {
    static const int n = 4;
    enum Indices { xxx=0, xxy=1, xyy=2, yyy = 3 };
    static void contractLastTwoIndexes(Array<T,n> const& tens, Array<T,2> & res) {
        res[0] = tens[xxx] + tens[xyy];
        res[1] = tens[xxy] + tens[yyy];
    }
    
    static T contractIndexes(Array<T,n> const& A, Array<T,n> const& B) {
        return A[xxx]*B[xxx] + A[yyy]*B[yyy] + (T)3*(A[xxy]*B[xxy] + A[xyy]*B[xyy]);
    }
    
    static void multWithRankTwoSymTensor(Array<T,n> const&A, Array<T,SymmetricTensorImpl<T,2>::n> const &B,
                                         Array<T,2> &x) {
        x[0] = A[xxx]*B[SymmetricTensorImpl<T,2>::xx] + (T)2*A[xxy]*B[SymmetricTensorImpl<T,2>::xy] + A[xyy]*B[SymmetricTensorImpl<T,2>::yy];
        x[1] = A[xxy]*B[SymmetricTensorImpl<T,2>::xx] + (T)2*A[xyy]*B[SymmetricTensorImpl<T,2>::xy] + A[yyy]*B[SymmetricTensorImpl<T,2>::yy];
    }
};

template<typename T> struct SymmetricRankThreeTensorImpl<T,3> {
    static const int n = 10;
    enum Indices { xxx=0, xxy=1, xxz=2, xyy=3, xyz=4, xzz=5, yyy=6, yyz=7, yzz=8, zzz=9 };
    static void contractLastTwoIndexes(Array<T,n> const& tens, Array<T,3> & res) {
        res[0] = tens[xxx] + tens[xyy] + tens[xzz];
        res[1] = tens[xxy] + tens[yyy] + tens[yzz];
        res[2] = tens[xxz] + tens[yyz] + tens[zzz];
    }
    
    static T contractIndexes(Array<T,n> const& A, Array<T,n> const& B) {
        return A[xxx]*B[xxx] + A[yyy]*B[yyy] + A[zzz]*B[zzz] 
            + (T)3*(A[xxy]*B[xxy] + A[xxz]*B[xxz] + A[xyy]*B[xyy] + A[xzz]*B[xzz] + A[yzz]*B[yzz] + A[yyz]*B[yyz])
            + (T)6*A[xyz]*B[xyz];
    }
    
    static void multWithRankTwoSymTensor(Array<T,n> const&A, Array<T,SymmetricTensorImpl<T,3>::n> const &B,
                                         Array<T,3> &x) {
        x[0] = A[xxx]*B[SymmetricTensorImpl<T,3>::xx] + A[xyy]*B[SymmetricTensorImpl<T,3>::yy] + A[xzz]*B[SymmetricTensorImpl<T,3>::zz]
            + (T)2*(A[xxy]*B[SymmetricTensorImpl<T,3>::xy] + A[xxz]*B[SymmetricTensorImpl<T,3>::xz] + A[xyz]*B[SymmetricTensorImpl<T,3>::yz]);
        x[1] = A[yyy]*B[SymmetricTensorImpl<T,3>::yy] + A[xxy]*B[SymmetricTensorImpl<T,3>::xx] + A[yzz]*B[SymmetricTensorImpl<T,3>::zz]
            + (T)2*(A[xyy]*B[SymmetricTensorImpl<T,3>::xy] + A[yyz]*B[SymmetricTensorImpl<T,3>::yz] + A[xyz]*B[SymmetricTensorImpl<T,3>::xz]);
        x[2] = A[zzz]*B[SymmetricTensorImpl<T,3>::zz] + A[xxz]*B[SymmetricTensorImpl<T,3>::xx] + A[yyz]*B[SymmetricTensorImpl<T,3>::yy]
            + (T)2*(A[xzz]*B[SymmetricTensorImpl<T,3>::xz] + A[yzz]*B[SymmetricTensorImpl<T,3>::yz] + A[xyz]*B[SymmetricTensorImpl<T,3>::xy]);
    }
};

/// Operations on a symmetric tensor which stores only above-or-on-diagonal values
template <typename T, template<typename U> class Descriptor>
struct SymmetricRankThreeTensor {
    /// Number of dimensions for current lattice
    static const int d = Descriptor<T>::d;
    /// Number of elements (reduced by symmetry)
    static const int n = SymmetricRankThreeTensorImpl<T,d>::n;
    /// computes the contraction of the last two indexes of the tensor
    static void contractLastTwoIndexes(Array<T,n> const& tens, Array<T,d> &res) {
        SymmetricRankThreeTensorImpl<T,d>::contractLastTwoIndexes(tens, res);
    }
    /// computes the contraction of the last two indexes of the tensor
    static T contractIndexes(Array<T,n> const& A, Array<T,n> const& B) {
        return SymmetricRankThreeTensorImpl<T,d>::contractIndexes(A, B);
    }
    static void multWithRankTwoSymTensor(Array<T,n> const&A, Array<T,SymmetricTensor<T,Descriptor>::n> const &B,
                                         Array<T,d> &x) {
        SymmetricRankThreeTensorImpl<T,d>::multWithRankTwoSymTensor(A, B, x);
    }
};

template <typename T, int d> struct SymmetricRankFourTensorImpl { };

template<typename T> struct SymmetricRankFourTensorImpl<T,2> {
    static const int n = 5;
    enum Indices { xxxx=0, xxxy=1, xxyy=2, xyyy = 3, yyyy = 4 };
    
    static T contractIndexes(Array<T,n> const& A, Array<T,n> const& B) {
        return A[xxxx]*B[xxxx] + A[yyyy]*B[yyyy] + (T)4*(A[xxxy]*B[xxxy] + A[xyyy]*B[xyyy]) + (T)6*A[xxyy]*B[xxyy];
    }
    /// computes the contraction of the last two indexes of the tensor
    static void contractLastTwoIndexes(Array<T,n> const& tens, Array<T,3> &res) {
        res[SymmetricTensorImpl<T,2>::xx] = tens[xxxx]+tens[xxyy];
        res[SymmetricTensorImpl<T,2>::xy] = tens[xxxy]+tens[xyyy];
        res[SymmetricTensorImpl<T,2>::yy] = tens[xxyy]+tens[yyyy];
    }
    /// computes the interior product of symmetric rank 4 with rank three tensor
    static void multWithRankThreeSymTensor(Array<T,n> const&A, 
                                           Array<T,SymmetricRankThreeTensorImpl<T,2>::n> const &B,
                                           Array<T,2> &x) {
        typedef SymmetricRankThreeTensorImpl<T,2> SRT;
        
        x[0] = A[xxxx]*B[SRT::xxx]+(T)3*(A[xxxy]*B[SRT::xxy]+A[xxyy]*B[SRT::xyy])+A[xyyy]*B[SRT::yyy];
        x[1] = A[xxxy]*B[SRT::xxx]+(T)3*(A[xxyy]*B[SRT::xxy]+A[xyyy]*B[SRT::xyy])+A[yyyy]*B[SRT::yyy];
    }
};

template<typename T> struct SymmetricRankFourTensorImpl<T,3> {
    static const int n = 15;
    enum Indices { xxxx=0, xxxy=1, xxxz=2, xxyy=3, xxyz=4, xxzz=5, xyyy=6, xyyz=7, xyzz=8, xzzz=9, 
                   yyyy=10, yyyz=11, yyzz=12, yzzz=13, zzzz=14 };
                   
    static T contractIndexes(Array<T,n> const& A, Array<T,n> const& B) {
        return A[xxxx]*B[xxxx] + A[yyyy]*B[yyyy] + A[zzzz]*B[zzzz] 
            + (T)4*(A[xxxy]*B[xxxy] + A[xxxz]*B[xxxz] + A[xyyy]*B[xyyy] + A[xzzz]*B[xzzz] + A[yyyz]*B[yyyz] + A[yzzz]*B[yzzz])
            + (T)6*(A[xxyy]*B[xxyy] + A[xxzz]*B[xxzz] + A[yyzz]*B[yyzz])
            + (T)12*(A[xxyz]*B[xxyz] + A[xyyz]*B[xyyz] + A[xyzz]*B[xyzz]);
    }
    /// computes the contraction of the last two indexes of the tensor
    static void contractLastTwoIndexes(Array<T,n> const& tens, Array<T,6> &res) {
        res[SymmetricTensorImpl<T,3>::xx] = tens[xxxx]+tens[xxyy]+tens[xxzz];
        res[SymmetricTensorImpl<T,3>::xy] = tens[xxxy]+tens[xyyy]+tens[xyzz];
        res[SymmetricTensorImpl<T,3>::xz] = tens[xxxz]+tens[xyyz]+tens[xzzz];
        res[SymmetricTensorImpl<T,3>::yy] = tens[xxyy]+tens[yyyy]+tens[yyzz];
        res[SymmetricTensorImpl<T,3>::yz] = tens[xxyz]+tens[yyyz]+tens[yzzz];
        res[SymmetricTensorImpl<T,3>::zz] = tens[xxzz]+tens[yyzz]+tens[zzzz];
    }
    static void multWithRankThreeSymTensor(Array<T,n> const&A, 
                                           Array<T,SymmetricRankThreeTensorImpl<T,3>::n> const &B,
                                           Array<T,3> &x) {
        typedef SymmetricRankThreeTensorImpl<T,3> SRT;
        
        x[0] = A[xxxx]*B[SRT::xxx]+A[xyyy]*B[SRT::yyy]+A[xzzz]*B[SRT::zzz]
            +(T)3*(A[xxxy]*B[SRT::xxy]+A[xxyy]*B[SRT::xyy]+A[xxxz]*B[SRT::xxz]+A[xyyz]*B[SRT::yyz]+A[xxzz]*B[SRT::xzz]+A[xyzz]*B[SRT::yzz])
            +(T)6*A  [xxyz]*B[SRT::xyz];
        
        x[1] = A[yyyy]*B[SRT::yyy]+A[xxxy]*B[SRT::xxx]+A[yzzz]*B[SRT::zzz]
            +(T)3*(A[xxyy]*B[SRT::xxy]+A[xyyy]*B[SRT::xyy]+A[yyyz]*B[SRT::yyz]+A[yyzz]*B[SRT::yzz]+A[xyzz]*B[SRT::xzz]+A[xxyz]*B[SRT::xxz])
            +(T)6*A[xyyz]*B[SRT::xyz];
        
        x[2] = A[zzzz]*B[SRT::zzz]+A[xxxz]*B[SRT::xxx]+A[yyyz]*B[SRT::yyy]
            +(T)3*(A[xxyz]*B[SRT::xxy]+A[xyyz]*B[SRT::xyy]+A[xxzz]*B[SRT::xxz]+A[yyzz]*B[SRT::yyz]+A[xzzz]*B[SRT::xzz]+A[yzzz]*B[SRT::yzz])
            +(T)6*A[xyzz]*B[SRT::xyz];
    }
};

/// Operations on a symmetric tensor which stores only above-or-on-diagonal values
template <typename T, template<typename U> class Descriptor>
struct SymmetricRankFourTensor {
    /// Number of dimensions for current lattice
    static const int d = Descriptor<T>::d;
    /// Number of elements (reduced by symmetry)
    static const int n = SymmetricRankFourTensorImpl<T,d>::n;
    
    /// computes the contraction of the last two indexes of the tensor
    static T contractIndexes(Array<T,n> const& A, Array<T,n> const& B) {
        return SymmetricRankFourTensorImpl<T,d>::contractIndexes(A, B);
    }
    /// computes the contraction of the last two indexes of the tensor
    static void contractLastTwoIndexes(Array<T,n> const& tens, Array<T,SymmetricTensor<T,Descriptor>::n> &res) {
        SymmetricRankFourTensorImpl<T,d>::contractLastTwoIndexes(tens, res);
    }
    static void multWithRankThreeSymTensor(Array<T,n> const&A, 
                                           Array<T,SymmetricRankThreeTensor<T,Descriptor>::n> const &B,
                                           Array<T,d> &x) {
        SymmetricRankFourTensorImpl<T,d>::multWithRankThreeSymTensor(A, B, x);
    }
};

template <typename T>
void crossProduct(Array<T,3> const& u1, Array<T,3> const& u2, Array<T,3>& result)
{
    result[0] = u1[1]*u2[2] - u1[2]*u2[1];
    result[1] = u1[2]*u2[0] - u1[0]*u2[2];
    result[2] = u1[0]*u2[1] - u1[1]*u2[0];
}

template <typename T>
Array<T,3> crossProduct(Array<T,3> const& u1, Array<T,3> const& u2)
{
    Array<T,3> result;
    crossProduct(u1, u2, result);

    return result;
}

/// Scalar product between two vectors.
template<typename T, pluint n>
T dot(Array<T,n> const& v1, Array<T,n> const& v2) {
    return VectorTemplateImpl<T,n>::scalarProduct(v1,v2);
}

/// Scalar product between a c-array and a vector.
template<typename T, pluint n>
T dot(T v1[n], Array<T,n> const& v2) {
    return VectorTemplateImpl<T,n>::scalarProduct(v1,v2);
}

template<typename T, pluint n>
T normSqr(Array<T,n> const& v) {
    return VectorTemplateImpl<T,n>::normSqr(v);
}

template<typename T, pluint n>
T norm(Array<T,n> const& v) {
    return std::sqrt(normSqr<T,n>(v));
}

/// Compute the normal vector for a given triangle. If "isAreaWeighted" is false,
///   then the normal has length equal to one. If "isAreaWeighted" is true, then
///   the normal has length equal to twice the area of the triangle.
template<typename T>
Array<T,3> computeTriangleNormal(Array<T,3> const& v0, Array<T,3> const& v1, Array<T,3> const& v2, bool isAreaWeighted)
{
#ifdef PLB_DEBUG
    static T eps = getEpsilon<T>(floatingPointPrecision<T>());
#endif

    Array<T,3> e01 = v1 - v0;
    Array<T,3> e02 = v2 - v0;

    Array<T,3> n;
    crossProduct<T>(e01, e02, n);
    T normN = norm<T,3>(n);
    PLB_ASSERT(normN > eps);
    if (!isAreaWeighted)
        n /= normN;

    return n;
}

template<typename T>
T computeTriangleArea(Array<T,3> const& v0, Array<T,3> const& v1, Array<T,3> const& v2)
{
    Array<T,3> e01 = v1 - v0;
    Array<T,3> e02 = v2 - v0;
    Array<T,3> cross;
    crossProduct<T>(e01, e02, cross);

    return 0.5 * norm<T,3>(cross);
}

template<typename T>
T computeTriangleArea(T *v0, T *v1, T *v2)
{
    int i;
    T e01[3], e02[3], cross[3];

    for (i = 0; i < 3; e01[i] = v1[i] - v0[i], i++)
        ;
    for (i = 0; i < 3; e02[i] = v2[i] - v0[i], i++)
        ;

    cross[0] = e01[1]*e02[2] - e01[2]*e02[1];
    cross[1] = e01[2]*e02[0] - e01[0]*e02[2];
    cross[2] = e01[0]*e02[1] - e01[1]*e02[0];

    return 0.5 * std::sqrt(cross[0]*cross[0] + cross[1]*cross[1] + cross[2]*cross[2]);
}

template<typename T>
T computeTetrahedronSignedVolume(Array<T,3> const& v0, Array<T,3> const& v1, Array<T,3> const& v2, Array<T,3> const& v3)
{
    // The position arrays specify two oriented triangles, namely 
    // (0, 1, 2) and (3, 2, 1). These triangles share
    // the common edge 1-2.

    Array<T,3> e10 = v0 - v1;
    Array<T,3> e12 = v2 - v1;
    Array<T,3> e13 = v3 - v1;
    Array<T,3> cross;
    crossProduct<T>(e12, e13, cross);

    return dot(e10, cross) / 6.0;
}

/// Given two vectors this function computes the angle between them
///   assuming that they originate at the same point in space
template<typename T>
T angleBetweenVectors(Array<T,3> const& v1, Array<T,3> const& v2)
{
    Array<T,3> cross;
    crossProduct<T>(v1, v2, cross); 
    return std::atan2(norm(cross), dot(v1,v2));
}

template<typename T>
Array<T,3> rotateAtOrigin(Array<T,3> const& p, Array<T,3> const& normedAxis, T theta) {
    Array<T,3> const& u = normedAxis;
    T d = std::sqrt(u[1]*u[1] + u[2]*u[2]);

    Array<T,3> q1 = p;
    Array<T,3> q2 = q1;

    T eps = 1.e-2;
    //T eps = (T)100.*std::numeric_limits<T>::epsilon();
    // Rotate about the x-axis to be in xz-plane.
    if (!util::fpequal(d, (T)0., eps)) {
        q2[0] = q1[0];
        q2[1] = q1[1] * u[2] / d - q1[2] * u[1] / d;
        q2[2] = q1[1] * u[1] / d + q1[2] * u[2] / d;
    }

    // Rotate about the y-axis to fall on z-axis.
    q1[0] = q2[0] * d - q2[2] * u[0];
    q1[1] = q2[1];
    q1[2] = q2[0] * u[0] + q2[2] * d;

    // Perform desired rotation.
    T ct = std::cos(theta);
    T st = std::sin(theta);
    q2[0] = q1[0] * ct - q1[1] * st;
    q2[1] = q1[0] * st + q1[1] * ct;
    q2[2] = q1[2];

    // Rotate backward around y-axis.
    q1[0] =   q2[0] * d + q2[2] * u[0];
    q1[1] =   q2[1];
    q1[2] = - q2[0] * u[0] + q2[2] * d;

    q2 = q1;
    // Rotate backward around x-axis.
    if (!util::fpequal(d, (T)0., eps)) {
        q2[0] =   q1[0];
        q2[1] =   q1[1] * u[2] / d + q1[2] * u[1] / d;
        q2[2] = - q1[1] * u[1] / d + q1[2] * u[2] / d;
    }

    return q2;
}

template<typename T>
Array<T,3> rotateWithEulerAngles(Array<T,3> const& p, T phi, T theta, T psi)
{
    T eps = getEpsilon<T>(floatingPointPrecision<T>());
    T pi = std::acos((T) -1.0);

    PLB_ASSERT((theta > (T) 0.0 || util::fpequal(theta, (T) 0.0, eps)) &&
               (theta < pi  || util::fpequal(theta, pi, eps)));

    T a[3][3];
    a[0][0] =  (T) 1.0;
    a[0][1] =  (T) 0.0;
    a[0][2] =  (T) 0.0;
    a[1][0] =  (T) 0.0;
    a[1][1] =  std::cos(theta);
    a[1][2] = -std::sin(theta);
    a[2][0] =  (T) 0.0;
    a[2][1] =  std::sin(theta);
    a[2][2] =  std::cos(theta);

    T b[3][3];
    b[0][0] =  std::cos(phi);
    b[0][1] = -std::sin(phi);
    b[0][2] =  (T) 0.0;
    b[1][0] =  std::sin(phi);
    b[1][1] =  std::cos(phi);
    b[1][2] =  (T) 0.0;
    b[2][0] =  (T) 0.0;
    b[2][1] =  (T) 0.0;
    b[2][2] =  (T) 1.0;

    T c[3][3];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            c[i][j] = (T) 0.0;
            for (int k = 0; k < 3; k++) {
                c[i][j] += a[i][k]*b[k][j];
            }
        }
    }

    b[0][0] =  std::cos(psi);
    b[0][1] = -std::sin(psi);
    b[0][2] =  (T) 0.0;
    b[1][0] =  std::sin(psi);
    b[1][1] =  std::cos(psi);
    b[1][2] =  (T) 0.0;
    b[2][0] =  (T) 0.0;
    b[2][1] =  (T) 0.0;
    b[2][2] =  (T) 1.0;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            a[i][j] = (T) 0.0;
            for (int k = 0; k < 3; k++) {
                a[i][j] += b[i][k]*c[k][j];
            }
        }
    }

    Array<T,3> x;
    for (int i = 0; i < 3; i++) {
        x[i] = (T) 0.0;
        for (int j = 0; j < 3; j++) {
            x[i] += a[i][j]*p[j];
        }
    }

    return x;
}

// This function returns the new position of the point "oldPosition", after a discrete rotation
// with time step equal to 1.0. It needs the angular velocity of the rotation, a normalized vector
// parallel to the axis of rotation, and a point on the axis of rotation.
template<typename T>
Array<T,3> getRotatedPosition(Array<T,3> const& oldPosition, Array<T,3> const& angularVelocity,
        Array<T,3> const& rotationAxisUnitVector, Array<T,3> const& pointOnRotationAxis)
{
    // Standard implementation
    /*
    Array<T,3> dp = oldPosition - pointOnRotationAxis;
    T theta = dot(angularVelocity, rotationAxisUnitVector);
    Array<T,3> rotatedPosition = rotateAtOrigin(dp, rotationAxisUnitVector, theta);
    Array<T,3> newPosition = rotatedPosition + pointOnRotationAxis;
    return(newPosition);
    */

    // Optimized implementation

    static T eps = getEpsilon<T>();
    static T epsSqr = eps*eps;

    Array<T,3> dp = oldPosition - pointOnRotationAxis;
    Array<T,3> parallelDp = dot(dp, rotationAxisUnitVector) * rotationAxisUnitVector;
    Array<T,3> oldR = dp - parallelDp;
    T normSqrOldR = normSqr(oldR);
    if (normSqrOldR <= epsSqr) {
        return(oldPosition);
    }
    Array<T,3> velocity = crossProduct(angularVelocity, oldR);
    Array<T,3> newR = oldR + velocity;
    T normSqrNewR = normSqr(newR);
    Array<T,3> correctedR = std::sqrt(normSqrOldR / normSqrNewR) * newR;
    Array<T,3> newPosition = correctedR + parallelDp + pointOnRotationAxis;
    return(newPosition);
}

// This function returns the discrete rotational velocity of the point "oldPosition", considering
// a time step equal to 1.0. It needs the angular velocity of the rotation, a normalized vector
// parallel to the axis of rotation, and a point on the axis of rotation.
template<typename T>
Array<T,3> getRotationalVelocity(Array<T,3> const& oldPosition, Array<T,3> const& angularVelocity,
        Array<T,3> const& rotationAxisUnitVector, Array<T,3> const& pointOnRotationAxis)
{
    Array<T,3> newPosition = getRotatedPosition(oldPosition, angularVelocity, rotationAxisUnitVector, pointOnRotationAxis);
    Array<T,3> velocity = newPosition - oldPosition;
    return(velocity);
}

}  // namespace plb

#endif
