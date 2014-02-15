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

/** \file
 * 2D specialization of dynamicsTemplates functions.
 */

#ifndef DYNAMICS_TEMPLATES_2D_H
#define DYNAMICS_TEMPLATES_2D_H

#include "core/globalDefs.h"
#include "latticeBoltzmann/nearestNeighborLattices2D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"

namespace plb {

/// Compute Pi tensor efficiently on D2Q9 lattice
template<typename T>
struct neqPiD2Q9 {
      
    typedef SymmetricTensorImpl<T,2> S;

    static T fromPiToFneq0(Array<T,3> const& pi) {
        return (T)2 * (-(T)1/(T)3*pi[S::xx] - (T)1/(T)3*pi[S::yy]);
    }

    static T fromPiToFneq1(Array<T,3> const& pi) {
        return (T)1/(T)4 * ((T)1/(T)3*pi[S::xx] + (T)1/(T)3*pi[S::yy] - pi[S::xy]);
    }

    static T fromPiToFneq2(Array<T,3> const& pi) {
        return (T)1/(T)2 * ((T)2/(T)3*pi[S::xx] - (T)1/(T)3*pi[S::yy]);
    }

    static T fromPiToFneq3(Array<T,3> const& pi) {
        return (T)1/(T)4 * ((T)1/(T)3*pi[S::xx] + (T)1/(T)3*pi[S::yy] + pi[S::xy]);
    }

    static T fromPiToFneq4(Array<T,3> const& pi) {
        return (T)1/(T)2 * (-(T)1/(T)3*pi[S::xx] + (T)2/(T)3*pi[S::yy]);
    }
};  //struct neqPiD2Q9

// Efficient specialization for D2Q9 base lattice
template<typename T> struct dynamicsTemplatesImpl<T, descriptors::D2Q9DescriptorBase<T> > {

typedef descriptors::D2Q9DescriptorBase<T> Descriptor;

static T bgk_ma2_equilibrium(plint iPop, T rhoBar, T invRho, Array<T,2> const& j, T jSqr ) {
    typedef descriptors::D2Q9DescriptorBase<T> L;
    T c_j = L::c[iPop][0]*j[0] + L::c[iPop][1]*j[1];
    return L::t[iPop] * ( rhoBar +
               3.*c_j + invRho*(4.5*c_j*c_j - 1.5*jSqr) );
}

static void bgk_ma2_equilibria( T rhoBar, T invRho, Array<T,Descriptor::d> const& j,
                                T jSqr, Array<T,Descriptor::q>& eqPop )
{
    T t0 = Descriptor::t[0];
    T t1 = Descriptor::t[1];
    T t2 = Descriptor::t[2];

    T kx     = (T)3 * j[0];
    T ky     = (T)3 * j[1];
    T kxSqr_ = invRho / (T)2 * kx*kx;
    T kySqr_ = invRho / (T)2 * ky*ky;
    T kxky_  = invRho * kx*ky;

    T C1 = rhoBar + invRho*(T)3*jSqr;
    T C2, C3;

    // i=0
    C3 = -kxSqr_ - kySqr_;
    eqPop[0] = t0 * (C1+C3);

    // i=1 and i=5
    C2 = -kx + ky;
    C3 = -kxky_;
    eqPop[1] = t1 * (C1+C2+C3);
    eqPop[5] = t1 * (C1-C2+C3);

    // i=2 and i=6
    C2 = -kx;
    C3 = -kySqr_;
    eqPop[2] = t2 * (C1+C2+C3);
    eqPop[6] = t2 * (C1-C2+C3);

    // i=3 and i=7
    C2 = -kx - ky;
    C3 = kxky_;
    eqPop[3] = t1 * (C1+C2+C3);
    eqPop[7] = t1 * (C1-C2+C3);

    // i=4 and i=8
    C2 = -ky;
    C3 = -kxSqr_;
    eqPop[4] = t2 * (C1+C2+C3);
    eqPop[8] = t2 * (C1-C2+C3);
}

static T bgk_ma2_collision_base(Array<T,Descriptor::q>& f, T rhoBar, Array<T,2> const& j, T omega, T invRho ) {
    T one_m_omega = (T)1 - omega;
    T t0_omega = Descriptor::t[0] * omega;
    T t1_omega = Descriptor::t[1] * omega;
    T t2_omega = Descriptor::t[2] * omega;

    T jSqr   = j[0]*j[0] + j[1]*j[1];
    T kx     = (T)3 * j[0];
    T ky     = (T)3 * j[1];
    T kxSqr_ = invRho / (T)2 * kx*kx;
    T kySqr_ = invRho / (T)2 * ky*ky;
    T kxky_  = invRho * kx*ky;

    T C1 = rhoBar + invRho*(T)3*jSqr;
    T C2, C3;

    // i=0
    C3 = -kxSqr_ - kySqr_;
    f[0] *= one_m_omega; f[0] += t0_omega * (C1+C3);

    // i=1 and i=5
    C2 = -kx + ky;
    C3 = -kxky_;
    f[1] *= one_m_omega; f[1] += t1_omega * (C1+C2+C3);
    f[5] *= one_m_omega; f[5] += t1_omega * (C1-C2+C3);

    // i=2 and i=6
    C2 = -kx;
    C3 = -kySqr_;
    f[2] *= one_m_omega; f[2] += t2_omega * (C1+C2+C3);
    f[6] *= one_m_omega; f[6] += t2_omega * (C1-C2+C3);

    // i=3 and i=7
    C2 = -kx - ky;
    C3 = kxky_;
    f[3] *= one_m_omega; f[3] += t1_omega * (C1+C2+C3);
    f[7] *= one_m_omega; f[7] += t1_omega * (C1-C2+C3);

    // i=4 and i=8
    C2 = -ky;
    C3 = -kxSqr_;
    f[4] *= one_m_omega; f[4] += t2_omega * (C1+C2+C3);
    f[8] *= one_m_omega; f[8] += t2_omega * (C1-C2+C3);

    return invRho*invRho*jSqr;
}

static T bgk_ma2_collision(Array<T,Descriptor::q>& f, T rhoBar, Array<T,2> const& j, T omega) {
    return bgk_ma2_collision_base(f, rhoBar, j, omega, Descriptor::invRho(rhoBar));
}

static T bgk_inc_collision(Array<T,Descriptor::q>& f, T rhoBar, Array<T,2> const& j, T omega, T invRho0=(T)1 ) {
    return bgk_ma2_collision_base(f, rhoBar, j, omega, invRho0);
}

static T rlb_collision(Array<T,Descriptor::q>& f, T rhoBar, T invRho, Array<T,2> const& j, Array<T,3> const& PiNeq, T omega ) {
    typedef dynamicsTemplatesImpl<T, descriptors::D2Q9DescriptorBase<T> > DH;
    const T jSqr = j[0]*j[0] + j[1]*j[1];

    T piNeq0 = neqPiD2Q9<T>::fromPiToFneq0(PiNeq);
    T piNeq1 = neqPiD2Q9<T>::fromPiToFneq1(PiNeq);
    T piNeq2 = neqPiD2Q9<T>::fromPiToFneq2(PiNeq);
    T piNeq3 = neqPiD2Q9<T>::fromPiToFneq3(PiNeq);
    T piNeq4 = neqPiD2Q9<T>::fromPiToFneq4(PiNeq);

    f[0] = DH::bgk_ma2_equilibrium(0, rhoBar, invRho, j, jSqr)
                  + ((T)1-omega) * piNeq0;
    f[1] = DH::bgk_ma2_equilibrium(1, rhoBar, invRho, j, jSqr)
                  + ((T)1-omega) * piNeq1;
    f[2] = DH::bgk_ma2_equilibrium(2, rhoBar, invRho, j, jSqr)
                  + ((T)1-omega) * piNeq2;
    f[3] = DH::bgk_ma2_equilibrium(3, rhoBar, invRho, j, jSqr)
                  + ((T)1-omega) * piNeq3;
    f[4] = DH::bgk_ma2_equilibrium(4, rhoBar, invRho, j, jSqr)
                  + ((T)1-omega) * piNeq4;
    f[5] = DH::bgk_ma2_equilibrium(5, rhoBar, invRho, j, jSqr)
                  + ((T)1-omega) * piNeq1;
    f[6] = DH::bgk_ma2_equilibrium(6, rhoBar, invRho, j, jSqr)
                  + ((T)1-omega) * piNeq2;
    f[7] = DH::bgk_ma2_equilibrium(7, rhoBar, invRho, j, jSqr)
                  + ((T)1-omega) * piNeq3;
    f[8] = DH::bgk_ma2_equilibrium(8, rhoBar, invRho, j, jSqr)
                  + ((T)1-omega) * piNeq4;
    return jSqr*invRho*invRho;
}

static T bgk_ma2_constRho_collision(Array<T,Descriptor::q>& f, T rhoBar, Array<T,2> const& j, T ratioRho, T omega) {
    typedef descriptors::D2Q9DescriptorBase<T> L;
    T invRho = L::invRho(rhoBar);
    const T jSqr = j[0]*j[0] + j[1]*j[1];
    for (plint iPop=0; iPop < L::q; ++iPop) {
        T feq = bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr );
        f[iPop] =
          ratioRho*feq + L::t[iPop]*(ratioRho-(T)1) +
          ((T)1-omega)*(f[iPop]-feq);
    }
    return jSqr*invRho*invRho;
}

static T precond_bgk_ma2_equilibrium(plint iPop, T rhoBar, T invRho, Array<T,2> const& j, T jSqr, T invGamma ) {
    typedef descriptors::D2Q9DescriptorBase<T> L;
    T c_j = L::c[iPop][0]*j[0] + L::c[iPop][1]*j[1];
    return L::t[iPop] * ( rhoBar +
               3.*c_j + invGamma*invRho*(4.5*c_j*c_j - 1.5*jSqr) );
}

static T precond_bgk_ma2_collision_base( Array<T,Descriptor::q>& f, T rhoBar, Array<T,2> const& j, T omega,
                                         T invGamma, bool incompressible )
{
    T invRho = incompressible ? (T)1 : Descriptor::invRho(rhoBar);
    T one_m_omega = (T)1 - omega;
    T t0_omega = Descriptor::t[0] * omega;
    T t1_omega = Descriptor::t[1] * omega;
    T t2_omega = Descriptor::t[2] * omega;

    T jSqr   = j[0]*j[0] + j[1]*j[1];
    T kx     = (T)3 * j[0];
    T ky     = (T)3 * j[1];
    T kxSqr_ = invGamma* invRho / (T)2 * kx*kx;
    T kySqr_ = invGamma* invRho / (T)2 * ky*ky;
    T kxky_  = invGamma* invRho * kx*ky;

    T C1 = rhoBar + invGamma*invRho*(T)3*jSqr;
    T C2, C3;

    // i=0
    C3 = -kxSqr_ - kySqr_;
    f[0] *= one_m_omega; f[0] += t0_omega * (C1+C3);

    // i=1 and i=5
    C2 = -kx + ky;
    C3 = -kxky_;
    f[1] *= one_m_omega; f[1] += t1_omega * (C1+C2+C3);
    f[5] *= one_m_omega; f[5] += t1_omega * (C1-C2+C3);

    // i=2 and i=6
    C2 = -kx;
    C3 = -kySqr_;
    f[2] *= one_m_omega; f[2] += t2_omega * (C1+C2+C3);
    f[6] *= one_m_omega; f[6] += t2_omega * (C1-C2+C3);

    // i=3 and i=7
    C2 = -kx - ky;
    C3 = kxky_;
    f[3] *= one_m_omega; f[3] += t1_omega * (C1+C2+C3);
    f[7] *= one_m_omega; f[7] += t1_omega * (C1-C2+C3);

    // i=4 and i=8
    C2 = -ky;
    C3 = -kxSqr_;
    f[4] *= one_m_omega; f[4] += t2_omega * (C1+C2+C3);
    f[8] *= one_m_omega; f[8] += t2_omega * (C1-C2+C3);

    return invRho*invRho*jSqr;
}

static T precond_bgk_ma2_collision( Array<T,Descriptor::q>& f, T rhoBar, Array<T,2> const& j,
                                    T omega, T invGamma )
{
    return precond_bgk_ma2_collision_base(f, rhoBar, j, omega, invGamma, false);
}

};  //struct dynamicsTemplatesImpl<D2Q9DescriptorBase>

}  // namespace plb

#endif  // DYNAMICS_TEMPLATES_2D_H
