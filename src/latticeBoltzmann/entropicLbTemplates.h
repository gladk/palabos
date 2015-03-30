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

/* Orestis Malaspinas contributed this code.
 */

/** \file
 * A collection of dynamics classes (e.g. BGK) with which a Cell object
 * can be instantiated -- header file.
 */
#ifndef ENTROPIC_LB_HELPERS_H
#define ENTROPIC_LB_HELPERS_H

#include "core/globalDefs.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include <cmath>

namespace plb {

template<typename T, template<typename U> class Descriptor>
struct entropicLbTemplates 
{
    /// Computation of equilibrium distribution
    static T equilibrium( plint iPop, T rho, Array<T,Descriptor<T>::d> const& u)
    {
        typedef Descriptor<T> L;
        const T invCs = std::sqrt(L::invCs2);
        const T sqt3 = std::sqrt((T)3.0);
        T prod = (T)1;
        for (int iD=0; iD < Descriptor<T>::d; ++iD)
        {
            T uc = u[iD] * invCs; // u[iD] / c_s

            prod *= ((T)2 - std::sqrt((T)1.0+uc*uc)) * 
                    std::pow(((T)2.0 / sqt3 * uc + 
                    std::sqrt((T)1.0+uc*uc))/((T)1.0-uc/sqt3),
                        L::c[iPop][iD]/sqt3*invCs);
        }
        return rho*L::t[iPop]*prod-L::SkordosFactor()*L::t[iPop];
    }
    
    /// Computation of equilibrium distribution
    static T equilibriumApprox( plint iPop, T rho, Array<T,Descriptor<T>::d> const& u)
    {
        typedef Descriptor<T> L;

        T uSqr = VectorTemplate<T,Descriptor>::normSqr(u);
        T cu = T();
        for (int iD=0; iD < Descriptor<T>::d; ++iD)
        {
            cu += L::c[iPop][iD]*u[iD];
        }
        
        return rho * L::t[iPop] * ((T)1.0 +
                cu*L::invCs2 - (T)0.5 * uSqr*L::invCs2 + (T)0.5*std::pow(L::invCs2,(T)2)*cu*cu
                - (T)0.5*std::pow(L::invCs2,(T)2)*cu*uSqr + std::pow(cu,(T)3)*std::pow(L::invCs2,(T)3)/(T)6.0
                + (T)0.125*uSqr*uSqr*std::pow(L::invCs2,(T)2) - (T)0.25*cu*cu*uSqr*std::pow(L::invCs2,(T)3)
                + std::pow(cu,(T)4)*std::pow(L::invCs2,(T)4)/(T)24.0)-L::SkordosFactor()*L::t[iPop];
    }
};

}

#include "latticeBoltzmann/entropicLbTemplates2D.h"
#include "latticeBoltzmann/entropicLbTemplates3D.h"

#endif
