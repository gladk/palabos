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
#ifndef HERMITE_POLYNOMIALS_TEMPLATES_H
#define HERMITE_POLYNOMIALS_TEMPLATES_H

#include "core/globalDefs.h"
#include "core/array.h"
#include "geometricOperationTemplates.h"
#include "momentTemplates.h"
#include <algorithm>
#include <vector>

namespace plb {

template <typename T, class Descriptor, int d> struct HermiteTemplateImpl;

template <typename T, template<typename U> class Descriptor>
struct HermiteTemplate {
    /// Number of dimensions for current lattice
    static const int d = Descriptor<T>::d;
    /// Compute order0 hermite polynomial
    static T order0(plint iPop) {
        PLB_ASSERT(iPop < Descriptor<T>::q && iPop >= 0);
        return (T)1;
    }
    /// Compute order1 hermite polynomial
    static Array<T,d> order1(plint iPop) {
        PLB_ASSERT(iPop < Descriptor<T>::q && iPop >= 0);
        return HermiteTemplateImpl<T,typename Descriptor<T>::BaseDescriptor,Descriptor<T>::d>::
            order1(iPop);
    }
    
    /// Compute order2 hermite polynomial
    static Array<T,SymmetricTensorImpl<T,d>::n > order2(plint iPop) {
        PLB_ASSERT(iPop < Descriptor<T>::q && iPop >= 0);
        return HermiteTemplateImpl<T,typename Descriptor<T>::BaseDescriptor,Descriptor<T>::d>::
            order2(iPop);
    }
    
    /// Compute order2 hermite polynomial contracted with symmetric tensor (only the hermite part)
    /// Symmetric parts are multiplied bz the appropriate amount
    static Array<T,SymmetricTensorImpl<T,d>::n > contractedOrder2(plint iPop) {
        PLB_ASSERT(iPop < Descriptor<T>::q && iPop >= 0);
        return HermiteTemplateImpl<T,typename Descriptor<T>::BaseDescriptor,Descriptor<T>::d>::
            contractedOrder2(iPop);
    }
    
    /// Compute order3 hermite polynomial
    static Array<T,SymmetricRankThreeTensorImpl<T,d>::n > order3(plint iPop) {
        PLB_ASSERT(iPop < Descriptor<T>::q && iPop >= 0);
        return HermiteTemplateImpl<T,typename Descriptor<T>::BaseDescriptor,Descriptor<T>::d>::
            order3(iPop);
    }
    
    /// Compute order3 hermite polynomial contracted with symmetric tensor (only the hermite part)
    /// Symmetric parts are multiplied bz the appropriate amount
    static Array<T,SymmetricRankThreeTensorImpl<T,d>::n > contractedOrder3(plint iPop) {
        PLB_ASSERT(iPop < Descriptor<T>::q && iPop >= 0);
        return HermiteTemplateImpl<T,typename Descriptor<T>::BaseDescriptor,Descriptor<T>::d>::
            contractedOrder3(iPop);
    }
    
    /// Compute order4 hermite polynomial
    static Array<T,SymmetricRankFourTensorImpl<T,d>::n > order4(plint iPop) {
    PLB_ASSERT(iPop < Descriptor<T>::q && iPop >= 0);
        return HermiteTemplateImpl<T,typename Descriptor<T>::BaseDescriptor,Descriptor<T>::d>::
            order4(iPop);
    }
    
    /// Compute order2 hermite polynomial contracted with symmetric tensor (only the hermite part)
    /// Symmetric parts are multiplied bz the appropriate amount
    static Array<T,SymmetricRankFourTensorImpl<T,d>::n > contractedOrder4(plint iPop) {
        PLB_ASSERT(iPop < Descriptor<T>::q && iPop >= 0);
        return HermiteTemplateImpl<T,typename Descriptor<T>::BaseDescriptor,Descriptor<T>::d>::
            contractedOrder4(iPop);
    }
    
    static plint numDecomposedVariables(plint order) {
        switch (order) {
            case 0 : {
                return 1;
                break;
            }
            case 1: {
                return 1+Descriptor<T>::d;
                break;
            }
            case 2: {
                return 1+Descriptor<T>::d+SymmetricTensor<T,Descriptor>::n;
                break;
            }
            case 3: {
                return 1+Descriptor<T>::d+SymmetricTensor<T,Descriptor>::n+SymmetricRankThreeTensor<T,Descriptor>::n;
                break;
            }
            case 4: {
                return 1+Descriptor<T>::d+SymmetricTensor<T,Descriptor>::n
                                         +SymmetricRankThreeTensor<T,Descriptor>::n
                                         +SymmetricRankFourTensor<T,Descriptor>::n;
                break;
            }
            default :
                return 0;
        }
    }
    
    // specialization in 2D of the decomposition in hermite polynomials
    static std::vector<T> decompose(const Cell<T,Descriptor> &cell, plint order) {
        PLB_ASSERT(order >= 0 && order <= 4 );
        std::vector<T> coeffs;
        coeffs.resize(numDecomposedVariables(order));
        
        if (order >= 0) {
            coeffs[0] = momentTemplates<T,Descriptor>::get_rhoBar(cell);
            if (order >= 1) {
                Array<T,d> a1;
                momentTemplates<T,Descriptor>::get_j(cell,a1);
                a1.to_cArray(&coeffs[1]);
            }
            if (order >= 2) {
                Array<T,SymmetricTensorImpl<T,d>::n> a2;
                a2.resetToZero();
                for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                    Array<T,SymmetricTensorImpl<T,d>::n> H2 = order2(iPop);
                    for (plint iPi = 0; iPi < SymmetricTensorImpl<T,d>::n; ++iPi) {
                        a2[iPi] += H2[iPi]*cell[iPop];
                    }
                }
                a2.to_cArray(&coeffs[1+Descriptor<T>::d]);
                if (order >= 3) {
                    Array<T,SymmetricRankThreeTensorImpl<T,d>::n> a3;
                    a3.resetToZero();
                    for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                        Array<T,SymmetricRankThreeTensorImpl<T,d>::n> H3 = order3(iPop);
                        for (plint iPi = 0; iPi < SymmetricRankThreeTensorImpl<T,d>::n; ++iPi) {
                            a3[iPi] += H3[iPi]*cell[iPop];
                        }
                    }
                    a3.to_cArray(&coeffs[1+Descriptor<T>::d]+SymmetricTensorImpl<T,d>::n);
                    if (order >= 4) {
                        Array<T,SymmetricRankFourTensorImpl<T,d>::n> a4;
                        a4.resetToZero();
                        for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
                            Array<T,SymmetricRankFourTensorImpl<T,d>::n> H4 = order4(iPop);
                            for (plint iPi = 0; iPi < SymmetricRankFourTensorImpl<T,d>::n; ++iPi) {
                                a4[iPi] += H4[iPi]*cell[iPop];
                            }
                        }
                        a4.to_cArray(&coeffs[1+Descriptor<T>::d]+SymmetricTensorImpl<T,d>::n+SymmetricRankThreeTensorImpl<T,d>::n);
                    }
                }
            }
        }
        
        return coeffs;
    }
    
    // specialization in 2D of the decomposition in hermite polynomials
    static void recompose(const std::vector<T> &coeffs, plint order, Cell<T,Descriptor> &cell) {
        PLB_ASSERT(order >= 0 && order <= 4 );
        
        static T invCs4 = Descriptor<T>::invCs2* Descriptor<T>::invCs2 / (T)2;
        static T invCs6 = invCs4 * Descriptor<T>::invCs2 / (T)3;
        static T invCs8 = invCs6 * Descriptor<T>::invCs2 / (T)4;
        
        for (plint iPop = 0; iPop < Descriptor<T>::q; ++iPop) {
            if (order >= 0) {
                cell[iPop] = coeffs[0];
                if (order >= 1) {
                    Array<T,d> a1;
                    a1.from_cArray(&coeffs[1]);
                    T c_u = Descriptor<T>::c[iPop][0] * a1[0];
                    for (plint iD = 1; iD < d; ++iD) c_u += Descriptor<T>::c[iPop][iD] * a1[iD];
                    cell[iPop] += Descriptor<T>::invCs2 * c_u;
                }
                if (order >= 2) {
                    Array<T,SymmetricTensor<T,Descriptor>::n> a2;
                    a2.from_cArray(&coeffs[1+d]);
                    Array<T,SymmetricTensor<T,Descriptor>::n> H2 = order2(iPop);
                    cell[iPop] += invCs4 * SymmetricTensor<T,Descriptor>::contractIndexes(H2,a2);
                    if (order >= 3) {
                        Array<T,SymmetricRankThreeTensor<T,Descriptor>::n> a3;
                        Array<T,SymmetricRankThreeTensor<T,Descriptor>::n> H3 = order3(iPop);
                        a3.from_cArray(&coeffs[1+d+SymmetricTensor<T,Descriptor>::n]);
                        cell[iPop] += invCs6 * SymmetricRankThreeTensor<T,Descriptor>::contractIndexes(H3,a3);
                        if (order >= 4) {
                            Array<T,SymmetricRankFourTensor<T,Descriptor>::n> a4;
                            Array<T,SymmetricRankFourTensor<T,Descriptor>::n> H4 = order4(iPop);
                            a4.from_cArray(&coeffs[1+d+SymmetricTensor<T,Descriptor>::n
                                                      +SymmetricRankThreeTensor<T,Descriptor>::n]);
                            cell[iPop] += invCs8 * SymmetricRankFourTensor<T,Descriptor>::contractIndexes(H4,a4);
                        }
                    }
                }
            }
            cell[iPop] *= Descriptor<T>::t[iPop];
        }
    }
};

template <typename T, class Descriptor>
struct HermiteTemplateImpl<T,Descriptor,2> {
    
    typedef Descriptor D;
    
    static T order0(plint iPop) {
        return (T)1;
    }
    
    static Array<T,2> order1(plint iPop) {
        Array<T,2> H1(D::c[iPop][0],D::c[iPop][1]);
        return H1;
    }
    
    /// Compute order2 hermite polynomial
    static Array<T,SymmetricTensorImpl<T,2>::n > order2(plint iPop) {
        Array<T,SymmetricTensorImpl<T,2>::n> H2;
        H2[SymmetricTensorImpl<T,2>::xx] = D::c[iPop][0]*D::c[iPop][0]-D::cs2;
        H2[SymmetricTensorImpl<T,2>::xy] = D::c[iPop][0]*D::c[iPop][1];
        H2[SymmetricTensorImpl<T,2>::yy] = D::c[iPop][1]*D::c[iPop][1]-D::cs2;
        
        return H2;
    }
    
    /// Compute order2 hermite polynomial contracted with symmetric tensor (only the hermite part)
    /// Symmetric parts are multiplied bz the appropriate amount
    static Array<T,SymmetricTensorImpl<T,2>::n > contractedOrder2(plint iPop) {
        Array<T,SymmetricTensorImpl<T,2>::n> H2;
        H2[SymmetricTensorImpl<T,2>::xx] = D::c[iPop][0]*D::c[iPop][0]-D::cs2;
        H2[SymmetricTensorImpl<T,2>::xy] = (T)2*D::c[iPop][0]*D::c[iPop][1];
        H2[SymmetricTensorImpl<T,2>::yy] = D::c[iPop][1]*D::c[iPop][1]-D::cs2;
        
        return H2;
    }
    
    /// Compute order3 hermite polynomial
    static Array<T,SymmetricRankThreeTensorImpl<T,2>::n > order3(plint iPop) {
        plint cx2 = D::c[iPop][0]*D::c[iPop][0];
        plint cy2 = D::c[iPop][1]*D::c[iPop][1];
        
        Array<T,SymmetricRankThreeTensorImpl<T,2>::n> H3;
        H3[SymmetricRankThreeTensorImpl<T,2>::xxx] = D::c[iPop][0]*(cx2-(T)3*D::cs2);
        H3[SymmetricRankThreeTensorImpl<T,2>::xxy] = D::c[iPop][1]*(cx2-D::cs2);
        H3[SymmetricRankThreeTensorImpl<T,2>::xyy] = D::c[iPop][0]*(cy2-D::cs2);
        H3[SymmetricRankThreeTensorImpl<T,2>::yyy] = D::c[iPop][1]*(cy2-(T)3*D::cs2);
        
        return H3;
    }
    
    /// Compute order3 hermite polynomial contracted with symmetric tensor (only the hermite part)
    /// Symmetric parts are multiplied bz the appropriate amount
    static Array<T,SymmetricRankThreeTensorImpl<T,2>::n > contractedOrder3(plint iPop) {
        plint cx2 = D::c[iPop][0]*D::c[iPop][0];
        plint cy2 = D::c[iPop][1]*D::c[iPop][1];
        
        Array<T,SymmetricRankThreeTensorImpl<T,2>::n> H3;
        H3[SymmetricRankThreeTensorImpl<T,2>::xxx] = D::c[iPop][0]*(cx2-(T)3*D::cs2);
        H3[SymmetricRankThreeTensorImpl<T,2>::xxy] = (T)3*D::c[iPop][1]*(cx2-D::cs2);
        H3[SymmetricRankThreeTensorImpl<T,2>::xyy] = (T)3*D::c[iPop][0]*(cy2-D::cs2);
        H3[SymmetricRankThreeTensorImpl<T,2>::yyy] = D::c[iPop][1]*(cy2-(T)3*D::cs2);
        
        return H3;
    }
    
    /// Compute order4 hermite polynomial
    static Array<T,SymmetricRankFourTensorImpl<T,2>::n > order4(plint iPop) {
        plint cx2 = D::c[iPop][0]*D::c[iPop][0];
        plint cy2 = D::c[iPop][1]*D::c[iPop][1];
        
        Array<T,SymmetricRankFourTensorImpl<T,2>::n> H4;
        H4[SymmetricRankFourTensorImpl<T,2>::xxxx] = 
            cx2 * (cx2 - (T)6*D::cs2) + (T)3*D::cs2*D::cs2;
        H4[SymmetricRankFourTensorImpl<T,2>::xxxy] = 
            D::c[iPop][0] * D::c[iPop][1] * (cx2 - (T)3*D::cs2);
        H4[SymmetricRankFourTensorImpl<T,2>::xxyy] = 
            cx2*cy2 - D::cs2*cx2 - D::cs2*cy2 + D::cs2*D::cs2;;
        H4[SymmetricRankFourTensorImpl<T,2>::xyyy] = 
            D::c[iPop][0] * D::c[iPop][1] * (cy2 - (T)3*D::cs2);
        H4[SymmetricRankFourTensorImpl<T,2>::yyyy] = 
            cy2 * (cy2 - (T)6*D::cs2) + (T)3*D::cs2*D::cs2;
        
        return H4;
    }
    
    /// Compute order4 hermite polynomial contracted with symmetric tensor (only the hermite part)
    /// Symmetric parts are multiplied bz the appropriate amount
    static Array<T,SymmetricRankFourTensorImpl<T,2>::n > contractedOrder4(plint iPop) {
        plint cx2 = D::c[iPop][0]*D::c[iPop][0];
        plint cy2 = D::c[iPop][1]*D::c[iPop][1];
        
        Array<T,SymmetricRankFourTensorImpl<T,2>::n> H4;
        H4[SymmetricRankFourTensorImpl<T,2>::xxxx] = 
            (T)cx2 * (cx2 - (T)6*D::cs2) + (T)3*D::cs2*D::cs2;
        H4[SymmetricRankFourTensorImpl<T,2>::xxxy] = 
            (T)4*D::c[iPop][0] * D::c[iPop][1] * (cx2 - (T)3*D::cs2);
        H4[SymmetricRankFourTensorImpl<T,2>::xxyy] = 
            (T)6*(cx2*cy2 - D::cs2*cx2 - D::cs2*cy2 + D::cs2*D::cs2);
        H4[SymmetricRankFourTensorImpl<T,2>::xyyy] = 
            (T)4*D::c[iPop][0] * D::c[iPop][1] * (cy2 - (T)3*D::cs2);
        H4[SymmetricRankFourTensorImpl<T,2>::yyyy] = 
            (T)cy2 * (cy2 - (T)6*D::cs2) + (T)3*D::cs2*D::cs2;
        
        return H4;
    }
};


template <typename T, class Descriptor>
struct HermiteTemplateImpl<T,Descriptor,3> {
    
    typedef Descriptor D;
    
    static T order0(plint iPop) {
        return (T)1;
    }
    
    static Array<T,3> order1(plint iPop) {
        Array<T,3> H1(D::c[iPop][0],D::c[iPop][1],D::c[iPop][2]);
        return H1;
    }
    
    /// Compute order2 hermite polynomial
    static Array<T,SymmetricTensorImpl<T,3>::n > order2(plint iPop) {
        Array<T,SymmetricTensorImpl<T,3>::n> H2;
        H2[SymmetricTensorImpl<T,3>::xx] = D::c[iPop][0]*D::c[iPop][0]-D::cs2;
        H2[SymmetricTensorImpl<T,3>::xy] = D::c[iPop][0]*D::c[iPop][1];
        H2[SymmetricTensorImpl<T,3>::xz] = D::c[iPop][0]*D::c[iPop][2];
        H2[SymmetricTensorImpl<T,3>::yy] = D::c[iPop][1]*D::c[iPop][1]-D::cs2;
        H2[SymmetricTensorImpl<T,3>::yz] = D::c[iPop][1]*D::c[iPop][2];
        H2[SymmetricTensorImpl<T,3>::zz] = D::c[iPop][2]*D::c[iPop][2]-D::cs2;
        
        return H2;
    }
    
    /// Compute order2 hermite polynomial contracted with symmetric tensor (only the hermite part)
    /// Symmetric parts are multiplied bz the appropriate amount
    static Array<T,SymmetricTensorImpl<T,3>::n > contractedOrder2(plint iPop) {
        Array<T,SymmetricTensorImpl<T,3>::n> H2;
        H2[SymmetricTensorImpl<T,3>::xx] =      D::c[iPop][0]*D::c[iPop][0]-D::cs2;
        H2[SymmetricTensorImpl<T,3>::xy] = (T)2*D::c[iPop][0]*D::c[iPop][1];
        H2[SymmetricTensorImpl<T,3>::xz] = (T)2*D::c[iPop][0]*D::c[iPop][2];
        H2[SymmetricTensorImpl<T,3>::yy] =      D::c[iPop][1]*D::c[iPop][1]-D::cs2;
        H2[SymmetricTensorImpl<T,3>::yz] = (T)2*D::c[iPop][1]*D::c[iPop][2];
        H2[SymmetricTensorImpl<T,3>::zz] =      D::c[iPop][2]*D::c[iPop][2]-D::cs2;
        
        return H2;
    }
    
    /// Compute order3 hermite polynomial
    static Array<T,SymmetricRankThreeTensorImpl<T,3>::n > order3(plint iPop) {
        plint cx2 = D::c[iPop][0]*D::c[iPop][0];
        plint cy2 = D::c[iPop][1]*D::c[iPop][1];
        plint cz2 = D::c[iPop][2]*D::c[iPop][2];
        
        Array<T,SymmetricRankThreeTensorImpl<T,3>::n> H3;
        H3[SymmetricRankThreeTensorImpl<T,3>::xxx] = D::c[iPop][0]*(cx2-(T)3*D::cs2);
        H3[SymmetricRankThreeTensorImpl<T,3>::xxy] = D::c[iPop][1]*(cx2-D::cs2);
        H3[SymmetricRankThreeTensorImpl<T,3>::xxz] = D::c[iPop][2]*(cx2-D::cs2);
        
        H3[SymmetricRankThreeTensorImpl<T,3>::xyy] = D::c[iPop][0]*(cy2-D::cs2);
        H3[SymmetricRankThreeTensorImpl<T,3>::xzz] = D::c[iPop][0]*(cz2-D::cs2);
        
        H3[SymmetricRankThreeTensorImpl<T,3>::xyz] = D::c[iPop][0]*D::c[iPop][1]*D::c[iPop][2];
        
        H3[SymmetricRankThreeTensorImpl<T,3>::yyz] = D::c[iPop][2]*(cy2-D::cs2);
        H3[SymmetricRankThreeTensorImpl<T,3>::yzz] = D::c[iPop][1]*(cz2-D::cs2);

        H3[SymmetricRankThreeTensorImpl<T,3>::yyy] = D::c[iPop][1]*(cy2-(T)3*D::cs2);
        
        H3[SymmetricRankThreeTensorImpl<T,3>::zzz] = D::c[iPop][2]*(cz2-(T)3*D::cs2);
        
        return H3;
    }
    
    /// Compute order3 hermite polynomial contracted with symmetric tensor (only the hermite part)
    /// Symmetric parts are multiplied bz the appropriate amount
    static Array<T,SymmetricRankThreeTensorImpl<T,3>::n > contractedOrder3(plint iPop) {
        plint cx2 = D::c[iPop][0]*D::c[iPop][0];
        plint cy2 = D::c[iPop][1]*D::c[iPop][1];
        plint cz2 = D::c[iPop][2]*D::c[iPop][2];
        
        Array<T,SymmetricRankThreeTensorImpl<T,3>::n> H3;
        H3[SymmetricRankThreeTensorImpl<T,3>::xxx] =      D::c[iPop][0]*(cx2-(T)3*D::cs2);
        H3[SymmetricRankThreeTensorImpl<T,3>::xxy] = (T)3*D::c[iPop][1]*(cx2-D::cs2);
        H3[SymmetricRankThreeTensorImpl<T,3>::xxz] = (T)3*D::c[iPop][2]*(cx2-D::cs2);
        
        H3[SymmetricRankThreeTensorImpl<T,3>::xyy] = (T)3*D::c[iPop][0]*(cy2-D::cs2);
        H3[SymmetricRankThreeTensorImpl<T,3>::xzz] = (T)3*D::c[iPop][0]*(cz2-D::cs2);
        
        H3[SymmetricRankThreeTensorImpl<T,3>::xyz] = (T)6*D::c[iPop][0]*D::c[iPop][1]*D::c[iPop][2];
        
        H3[SymmetricRankThreeTensorImpl<T,3>::yyz] = (T)3*D::c[iPop][2]*(cy2-D::cs2);
        H3[SymmetricRankThreeTensorImpl<T,3>::yzz] = (T)3*D::c[iPop][1]*(cz2-D::cs2);
        
        H3[SymmetricRankThreeTensorImpl<T,3>::yyy] =      D::c[iPop][1]*(cy2-(T)3*D::cs2);
        
        H3[SymmetricRankThreeTensorImpl<T,3>::zzz] =      D::c[iPop][2]*(cz2-(T)3*D::cs2);
        
        return H3;
    }
    
    /// Compute order4 hermite polynomial
    static Array<T,SymmetricRankFourTensorImpl<T,3>::n > order4(plint iPop) {
        plint cx2 = D::c[iPop][0]*D::c[iPop][0];
        plint cy2 = D::c[iPop][1]*D::c[iPop][1];
        plint cz2 = D::c[iPop][2]*D::c[iPop][2];
        
        //plint cxy = D::c[iPop][0]*D::c[iPop][1];
        //plint cxz = D::c[iPop][0]*D::c[iPop][2];
        //plint cyz = D::c[iPop][1]*D::c[iPop][2];
        
        Array<T,SymmetricRankFourTensorImpl<T,3>::n> H4;
        H4[SymmetricRankFourTensorImpl<T,3>::xxxx] = cx2 * (cx2 - (T)6*D::cs2) + (T)3*D::cs2*D::cs2;
        H4[SymmetricRankFourTensorImpl<T,3>::yyyy] = cy2 * (cy2 - (T)6*D::cs2) + (T)3*D::cs2*D::cs2;
        H4[SymmetricRankFourTensorImpl<T,3>::zzzz] = cz2 * (cz2 - (T)6*D::cs2) + (T)3*D::cs2*D::cs2;
        
        H4[SymmetricRankFourTensorImpl<T,3>::xxxy] = D::c[iPop][0] * D::c[iPop][1] * (cx2 - (T)3*D::cs2);
        H4[SymmetricRankFourTensorImpl<T,3>::xxxz] = D::c[iPop][0] * D::c[iPop][2] * (cx2 - (T)3*D::cs2);
        H4[SymmetricRankFourTensorImpl<T,3>::xyyy] = D::c[iPop][0] * D::c[iPop][1] * (cy2 - (T)3*D::cs2);
        H4[SymmetricRankFourTensorImpl<T,3>::xzzz] = D::c[iPop][0] * D::c[iPop][2] * (cz2 - (T)3*D::cs2);
        H4[SymmetricRankFourTensorImpl<T,3>::yyyz] = D::c[iPop][2] * D::c[iPop][1] * (cy2 - (T)3*D::cs2);
        H4[SymmetricRankFourTensorImpl<T,3>::yzzz] = D::c[iPop][2] * D::c[iPop][1] * (cz2 - (T)3*D::cs2);
        
        H4[SymmetricRankFourTensorImpl<T,3>::xxyy] = cx2*cy2 - D::cs2*(cx2 - cy2) + D::cs2*D::cs2;
        H4[SymmetricRankFourTensorImpl<T,3>::xxzz] = cx2*cz2 - D::cs2*(cx2 - cz2) + D::cs2*D::cs2;
        H4[SymmetricRankFourTensorImpl<T,3>::yyzz] = cy2*cz2 - D::cs2*(cy2 - cz2) + D::cs2*D::cs2;
        
        H4[SymmetricRankFourTensorImpl<T,3>::xxyz] = D::c[iPop][1] * D::c[iPop][2] * (cx2 - D::cs2);
        H4[SymmetricRankFourTensorImpl<T,3>::xyyz] = D::c[iPop][0] * D::c[iPop][2] * (cy2 - D::cs2);
        H4[SymmetricRankFourTensorImpl<T,3>::xyzz] = D::c[iPop][0] * D::c[iPop][1] * (cz2 - D::cs2);
        
        return H4;
    }
    
    /// Compute order4 hermite polynomial contracted with symmetric tensor (only the hermite part)
    /// Symmetric parts are multiplied bz the appropriate amount
    static Array<T,SymmetricRankFourTensorImpl<T,3>::n > contractedOrder4(plint iPop) {
        plint cx2 = D::c[iPop][0]*D::c[iPop][0];
        plint cy2 = D::c[iPop][1]*D::c[iPop][1];
        plint cz2 = D::c[iPop][2]*D::c[iPop][2];
        
        //plint cxy = D::c[iPop][0]*D::c[iPop][1];
        //plint cxz = D::c[iPop][0]*D::c[iPop][2];
        //plint cyz = D::c[iPop][1]*D::c[iPop][2];
        
        Array<T,SymmetricRankFourTensorImpl<T,3>::n> H4;
        H4[SymmetricRankFourTensorImpl<T,3>::xxxx] = cx2 * (cx2 - (T)6*D::cs2) + (T)3*D::cs2*D::cs2;
        H4[SymmetricRankFourTensorImpl<T,3>::yyyy] = cy2 * (cy2 - (T)6*D::cs2) + (T)3*D::cs2*D::cs2;
        H4[SymmetricRankFourTensorImpl<T,3>::zzzz] = cz2 * (cz2 - (T)6*D::cs2) + (T)3*D::cs2*D::cs2;
        
        H4[SymmetricRankFourTensorImpl<T,3>::xxxy] = (T)4*D::c[iPop][0] * D::c[iPop][1] * (cx2 - (T)3*D::cs2);
        H4[SymmetricRankFourTensorImpl<T,3>::xxxz] = (T)4*D::c[iPop][0] * D::c[iPop][2] * (cx2 - (T)3*D::cs2);
        H4[SymmetricRankFourTensorImpl<T,3>::xyyy] = (T)4*D::c[iPop][0] * D::c[iPop][1] * (cy2 - (T)3*D::cs2);
        H4[SymmetricRankFourTensorImpl<T,3>::xzzz] = (T)4*D::c[iPop][0] * D::c[iPop][2] * (cz2 - (T)3*D::cs2);
        H4[SymmetricRankFourTensorImpl<T,3>::yyyz] = (T)4*D::c[iPop][2] * D::c[iPop][1] * (cy2 - (T)3*D::cs2);
        H4[SymmetricRankFourTensorImpl<T,3>::yzzz] = (T)4*D::c[iPop][2] * D::c[iPop][1] * (cz2 - (T)3*D::cs2);
        
        H4[SymmetricRankFourTensorImpl<T,3>::xxyy] = (T)6*cx2*cy2 - D::cs2*(cx2 - cy2) + D::cs2*D::cs2;
        H4[SymmetricRankFourTensorImpl<T,3>::xxzz] = (T)6*cx2*cz2 - D::cs2*(cx2 - cz2) + D::cs2*D::cs2;
        H4[SymmetricRankFourTensorImpl<T,3>::yyzz] = (T)6*cy2*cz2 - D::cs2*(cy2 - cz2) + D::cs2*D::cs2;
        
        H4[SymmetricRankFourTensorImpl<T,3>::xxyz] = (T)12*D::c[iPop][1] * D::c[iPop][2] * (cx2 - D::cs2);
        H4[SymmetricRankFourTensorImpl<T,3>::xyyz] = (T)12*D::c[iPop][0] * D::c[iPop][2] * (cy2 - D::cs2);
        H4[SymmetricRankFourTensorImpl<T,3>::xyzz] = (T)12*D::c[iPop][0] * D::c[iPop][1] * (cz2 - D::cs2);
        
        return H4;
    }
};

}  // namespace plb

#endif
