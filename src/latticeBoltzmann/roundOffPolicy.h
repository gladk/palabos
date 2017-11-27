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
 * Definition of policies for minimization of round-off errors -- header file
 */
#ifndef ROUND_OFF_POLICY_H
#define ROUND_OFF_POLICY_H

#include "core/globalDefs.h"

namespace plb {


/// Convert a Skordos-optimized representation of f (i.e. fBar) into the conventional one.
/** This function cannot be part of the RoundOffPolicy class, because it needs
 *  awareness of specific lattice information (the t_i).
 */
template<typename T, template<typename U> class Descriptor>
inline T fullF(T fBar, plint iPop) {
    return fBar + Descriptor<T>::SkordosFactor()*Descriptor<T>::t[iPop];
}

/// Convert a conventional representation of the f into a Skordos-optimized one (i.e. fBar).
/** This function cannot be part of the RoundOffPolicy class, because it needs
 *  awareness of specific lattice information (the t_i).
 */
template<typename T, template<typename U> class Descriptor>
inline T fBar(T fullF, plint iPop) {
    return fullF - Descriptor<T>::SkordosFactor()*Descriptor<T>::t[iPop];
}

template<typename T>
struct NoOptimizationRoundOffPolicy {
    static int SkordosFactor() {
        return 0;
    }
    static T rhoBar(T rho) {
        return rho;
    }
    static T fullRho(T rhoBar) {
        return rhoBar;
    }
    static T invRho(T rhoBar) {
        return (T)1 / rhoBar;
    }
    static T rhoMinus1(T rhoBar) {
        return rhoBar - (T)1;
    }
};

template<typename T>
struct DefaultRoundOffPolicy {
    static int SkordosFactor() {
        return 1;
    }
    static T rhoBar(T rho) {
        return rho - (T)1;
    }
    static T fullRho(T rhoBar) {
        return rhoBar + (T)1;
    }
    static T invRho(T rhoBar) {
        return (T)1 / (rhoBar + (T)1);
    }
    static T rhoMinus1(T rhoBar) {
        return rhoBar;
    }
};

template<typename T>
struct LinearRoundOffPolicy {
    static int SkordosFactor() {
        return 1;
    }
    static T rhoBar(T rho) {
        return rho - (T)1;
    }
    static T fullRho(T rhoBar) {
        return rhoBar + (T)1;
    }
    static T invRho(T rhoBar) {
        return (T)1 - rhoBar;
    }
    static T rhoMinus1(T rhoBar) {
        return rhoBar;
    }
};

template<typename T>
struct Order2RoundOffPolicy {
    static int SkordosFactor() {
        return 1;
    }
    static T rhoBar(T rho) {
        return rho - (T)1;
    }
    static T fullRho(T rhoBar) {
        return rhoBar + (T)1;
    }
    static T invRho(T rhoBar) {
        return (T)1 - rhoBar + rhoBar*rhoBar;
    }
    static T rhoMinus1(T rhoBar) {
        return rhoBar;
    }
};

template<typename T>
struct Order3RoundOffPolicy {
    static int SkordosFactor() {
        return 1;
    }
    static T rhoBar(T rho) {
        return rho - (T)1;
    }
    static T fullRho(T rhoBar) {
        return rhoBar + (T)1;
    }
    static T invRho(T rhoBar) {
        return (T)1 - rhoBar + rhoBar*rhoBar - rhoBar*rhoBar*rhoBar;
    }
    static T rhoMinus1(T rhoBar) {
        return rhoBar;
    }
};

}  // namespace plb

#endif
