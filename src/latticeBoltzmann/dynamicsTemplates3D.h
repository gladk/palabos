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
 * 3D specialization of dynamicsTemplates functions.
 */

#ifndef DYNAMICS_TEMPLATES_3D_H
#define DYNAMICS_TEMPLATES_3D_H

#include "core/globalDefs.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"

namespace plb {

template<typename T>
struct neqPiD3Q27 {

    typedef SymmetricTensorImpl<T,3> S;

    static T fromPiToFneq0(Array<T,6> const& pi) {
        return (T)4/(T)3 * (
                 -(T)1/(T)3*pi[S::xx] - (T)1/(T)3*pi[S::yy] - (T)1/(T)3*pi[S::zz]
               );
    }

    static T fromPiToFneq1(Array<T,6> const& pi) {
        return (T)1/(T)3 * (
                  (T)2/(T)3*pi[S::xx] - (T)1/(T)3*pi[S::yy] - (T)1/(T)3*pi[S::zz]
               );
    }

    static T fromPiToFneq2(Array<T,6> const& pi) {
        return (T)1/(T)3 * (
                  -(T)1/(T)3*pi[S::xx] + (T)2/(T)3*pi[S::yy] - (T)1/(T)3*pi[S::zz]
               );
    }

    static T fromPiToFneq3(Array<T,6> const& pi) {
        return (T)1/(T)3 * (
                  -(T)1/(T)3*pi[S::xx] - (T)1/(T)3*pi[S::yy] + (T)2/(T)3*pi[S::zz]
               );
    }

    static T fromPiToFneq4(Array<T,6> const& pi) {
        return (T)1/(T)12 * (
                  (T)2 * pi[S::xy]
                  + (T)2/(T)3*pi[S::xx] + (T)2/(T)3*pi[S::yy] - (T)1/(T)3*pi[S::zz]
               );
    }

    static T fromPiToFneq5(Array<T,6> const& pi) {
        return (T)1/(T)12 * (
                  - (T)2 * pi[S::xy]
                  + (T)2/(T)3*pi[S::xx] + (T)2/(T)3*pi[S::yy] - (T)1/(T)3*pi[S::zz]
               );
    }

    static T fromPiToFneq6(Array<T,6> const& pi) {
        return (T)1/(T)12 * (
                  (T)2 * pi[S::xz]
                  + (T)2/(T)3*pi[S::xx] - (T)1/(T)3*pi[S::yy] + (T)2/(T)3*pi[S::zz]
               );
    }

    static T fromPiToFneq7(Array<T,6> const& pi) {
        return (T)1/(T)12 * (
                  - (T)2 * pi[S::xz]
                  + (T)2/(T)3*pi[S::xx] - (T)1/(T)3*pi[S::yy] + (T)2/(T)3*pi[S::zz]
               );
    }

    static T fromPiToFneq8(Array<T,6> const& pi) {
        return (T)1/(T)12 * (
                  (T)2 * pi[S::yz]
                  - (T)1/(T)3*pi[S::xx] + (T)2/(T)3*pi[S::yy] + (T)2/(T)3*pi[S::zz]
               );
    }

    static T fromPiToFneq9(Array<T,6> const& pi) {
        return (T)1/(T)12 * (
                  - (T)2 * pi[S::yz]
                  - (T)1/(T)3*pi[S::xx] + (T)2/(T)3*pi[S::yy] + (T)2/(T)3*pi[S::zz]
               );
    }

    static T fromPiToFneq10(Array<T,6> const& pi) {
        return (T)1/(T)48 * (
                  + (T)2*pi[S::xy] + (T)2*pi[S::xz] + (T)2*pi[S::yz]
                  + (T)2/(T)3*pi[S::xx] + (T)2/(T)3*pi[S::yy] + (T)2/(T)3*pi[S::zz]
               );
    }

    static T fromPiToFneq11(Array<T,6> const& pi) {
        return (T)1/(T)48 * (
                  + (T)2*pi[S::xy] - (T)2*pi[S::xz] - (T)2*pi[S::yz]
                  + (T)2/(T)3*pi[S::xx] + (T)2/(T)3*pi[S::yy] + (T)2/(T)3*pi[S::zz]
               );
    }

    static T fromPiToFneq12(Array<T,6> const& pi) {
        return (T)1/(T)48 * (
                  - (T)2*pi[S::xy] + (T)2*pi[S::xz] - (T)2*pi[S::yz]
                  + (T)2/(T)3*pi[S::xx] + (T)2/(T)3*pi[S::yy] + (T)2/(T)3*pi[S::zz]
               );
    }

    static T fromPiToFneq13(Array<T,6> const& pi) {
        return (T)1/(T)48 * (
                  - (T)2*pi[S::xy] - (T)2*pi[S::xz] + (T)2*pi[S::yz]
                  + (T)2/(T)3*pi[S::xx] + (T)2/(T)3*pi[S::yy] + (T)2/(T)3*pi[S::zz]
               );
    }

};

// Efficient specialization for D3Q27 lattice
template<typename T>
struct dynamicsTemplatesImpl<T, descriptors::D3Q27DescriptorBase<T> > {

typedef descriptors::D3Q27DescriptorBase<T> Descriptor;

static T bgk_ma2_equilibrium(plint iPop, T rhoBar, T invRho, Array<T,3> const& j, T jSqr ) {
    T c_j = Descriptor::c[iPop][0]*j[0] + Descriptor::c[iPop][1]*j[1] + Descriptor::c[iPop][2]*j[2];
    return Descriptor::t[iPop] * ( rhoBar + 3.*c_j + invRho*(4.5*c_j*c_j - 1.5*jSqr) );
}

static void bgk_ma2_equilibria( T rhoBar, T invRho, Array<T,Descriptor::d> const& j,
                                T jSqr, Array<T,Descriptor::q>& eqPop )
{
    T t0 = Descriptor::t[0];
    T t1 = Descriptor::t[1];
    T t4 = Descriptor::t[4];
    T t10 = Descriptor::t[10];

    T kx     = (T)3 * j[0];
    T ky     = (T)3 * j[1];
    T kz     = (T)3 * j[2];
    T kxSqr_ = invRho / (T)2 * kx*kx;
    T kySqr_ = invRho / (T)2 * ky*ky;
    T kzSqr_ = invRho / (T)2 * kz*kz;
    T kxky_  = invRho * kx*ky;
    T kxkz_  = invRho * kx*kz;
    T kykz_  = invRho * ky*kz;

    T C1 = rhoBar + invRho*(T)3*jSqr;
    T C2, C3;

    // i=0
    C3 = -kxSqr_ - kySqr_ - kzSqr_;
    eqPop[0] = t0 * (C1+C3);

    // i=1 and i=14
    C2 = -kx;
    C3 = -kySqr_ - kzSqr_;
    eqPop[1]  = t1 * (C1+C2+C3);
    eqPop[14] = t1 * (C1-C2+C3);

    // i=2 and i=15
    C2 = -ky;
    C3 = -kxSqr_ - kzSqr_;
    eqPop[2]  = t1 * (C1+C2+C3);
    eqPop[15] = t1 * (C1-C2+C3);

    // i=3 and i=16
    C2 = -kz;
    C3 = -kxSqr_ - kySqr_;
    eqPop[3]  = t1 * (C1+C2+C3);
    eqPop[16] = t1 * (C1-C2+C3);

    // i=4 and i=17
    C2 = -kx - ky;
    C3 = kxky_ - kzSqr_;
    eqPop[4]  = t4 * (C1+C2+C3);
    eqPop[17] = t4 * (C1-C2+C3);

    // i=5 and i=18
    C2 = -kx + ky;
    C3 = -kxky_ - kzSqr_;
    eqPop[5]  = t4 * (C1+C2+C3);
    eqPop[18] = t4 * (C1-C2+C3);

    // i=6 and i=19
    C2 = -kx - kz;
    C3 = kxkz_ - kySqr_;
    eqPop[6]  = t4 * (C1+C2+C3);
    eqPop[19] = t4 * (C1-C2+C3);

    // i=7 and i=20
    C2 = -kx + kz;
    C3 = -kxkz_ - kySqr_;
    eqPop[7]  = t4 * (C1+C2+C3);
    eqPop[20] = t4 * (C1-C2+C3);

    // i=8 and i=21
    C2 = -ky - kz;
    C3 = kykz_ - kxSqr_;
    eqPop[8]  = t4 * (C1+C2+C3);
    eqPop[21] = t4 * (C1-C2+C3);

    // i=9 and i=22
    C2 = -ky + kz;
    C3 = -kykz_ - kxSqr_;
    eqPop[9]  = t4 * (C1+C2+C3);
    eqPop[22] = t4 * (C1-C2+C3);

    // i=10 and i=23
    C2 = -kx -ky -kz;
    C3 = kxky_ + kxkz_ + kykz_;
    eqPop[10] = t10 * (C1+C2+C3);
    eqPop[23] = t10 * (C1-C2+C3);

    // i=11 and i=24
    C2 = -kx -ky +kz;
    C3 = kxky_ - kxkz_ - kykz_;
    eqPop[11] = t10 * (C1+C2+C3);
    eqPop[24] = t10 * (C1-C2+C3);

    // i=12 and i=25
    C2 = -kx +ky -kz;
    C3 = -kxky_ + kxkz_ - kykz_;
    eqPop[12] = t10 * (C1+C2+C3);
    eqPop[25] = t10 * (C1-C2+C3);

    // i=13 and i=26
    C2 = -kx +ky +kz;
    C3 = -kxky_ - kxkz_ + kykz_;
    eqPop[13] = t10 * (C1+C2+C3);
    eqPop[26] = t10 * (C1-C2+C3);
}

static T bgk_ma2_collision_base(Array<T,Descriptor::q>& f, T rhoBar, Array<T,3> const& j, T omega, T invRho ) {
    T one_m_omega = (T)1 - omega;
    T t0_omega = Descriptor::t[0] * omega;
    T t1_omega = Descriptor::t[1] * omega;
    T t4_omega = Descriptor::t[4] * omega;
    T t10_omega = Descriptor::t[10] * omega;

    T jSqr   = j[0]*j[0] + j[1]*j[1] + j[2]*j[2];
    T kx     = (T)3 * j[0];
    T ky     = (T)3 * j[1];
    T kz     = (T)3 * j[2];
    T kxSqr_ = invRho / (T)2 * kx*kx;
    T kySqr_ = invRho / (T)2 * ky*ky;
    T kzSqr_ = invRho / (T)2 * kz*kz;
    T kxky_  = invRho * kx*ky;
    T kxkz_  = invRho * kx*kz;
    T kykz_  = invRho * ky*kz;

    T C1 = rhoBar + invRho*(T)3*jSqr;
    T C2, C3;

    // i=0
    C3 = -kxSqr_ - kySqr_ - kzSqr_;
    f[0] *= one_m_omega; f[0] += t0_omega * (C1+C3);

    // i=1 and i=14
    C2 = -kx;
    C3 = -kySqr_ - kzSqr_;
    f[1]  *= one_m_omega; f[1]  += t1_omega * (C1+C2+C3);
    f[14] *= one_m_omega; f[14] += t1_omega * (C1-C2+C3);

    // i=2 and i=15
    C2 = -ky;
    C3 = -kxSqr_ - kzSqr_;
    f[2]  *= one_m_omega; f[2]  += t1_omega * (C1+C2+C3);
    f[15] *= one_m_omega; f[15] += t1_omega * (C1-C2+C3);

    // i=3 and i=16
    C2 = -kz;
    C3 = -kxSqr_ - kySqr_;
    f[3]  *= one_m_omega; f[3]  += t1_omega * (C1+C2+C3);
    f[16] *= one_m_omega; f[16] += t1_omega * (C1-C2+C3);

    // i=4 and i=17
    C2 = -kx - ky;
    C3 = kxky_ - kzSqr_;
    f[4]  *= one_m_omega; f[4]  += t4_omega * (C1+C2+C3);
    f[17] *= one_m_omega; f[17] += t4_omega * (C1-C2+C3);

    // i=5 and i=18
    C2 = -kx + ky;
    C3 = -kxky_ - kzSqr_;
    f[5]  *= one_m_omega; f[5]  += t4_omega * (C1+C2+C3);
    f[18] *= one_m_omega; f[18] += t4_omega * (C1-C2+C3);

    // i=6 and i=19
    C2 = -kx - kz;
    C3 = kxkz_ - kySqr_;
    f[6]  *= one_m_omega; f[6]  += t4_omega * (C1+C2+C3);
    f[19] *= one_m_omega; f[19] += t4_omega * (C1-C2+C3);

    // i=7 and i=20
    C2 = -kx + kz;
    C3 = -kxkz_ - kySqr_;
    f[7]  *= one_m_omega; f[7]  += t4_omega * (C1+C2+C3);
    f[20] *= one_m_omega; f[20] += t4_omega * (C1-C2+C3);

    // i=8 and i=21
    C2 = -ky - kz;
    C3 = kykz_ - kxSqr_;
    f[8]  *= one_m_omega; f[8]  += t4_omega * (C1+C2+C3);
    f[21] *= one_m_omega; f[21] += t4_omega * (C1-C2+C3);

    // i=9 and i=22
    C2 = -ky + kz;
    C3 = -kykz_ - kxSqr_;
    f[9]  *= one_m_omega; f[9]  += t4_omega * (C1+C2+C3);
    f[22] *= one_m_omega; f[22] += t4_omega * (C1-C2+C3);

    // i=10 and i=23
    C2 = -kx -ky -kz;
    C3 = kxky_ + kxkz_ + kykz_;
    f[10] *= one_m_omega; f[10] += t10_omega * (C1+C2+C3);
    f[23] *= one_m_omega; f[23] += t10_omega * (C1-C2+C3);

    // i=11 and i=24
    C2 = -kx -ky +kz;
    C3 = kxky_ - kxkz_ - kykz_;
    f[11] *= one_m_omega; f[11] += t10_omega * (C1+C2+C3);
    f[24] *= one_m_omega; f[24] += t10_omega * (C1-C2+C3);

    // i=12 and i=25
    C2 = -kx +ky -kz;
    C3 = -kxky_ + kxkz_ - kykz_;
    f[12] *= one_m_omega; f[12] += t10_omega * (C1+C2+C3);
    f[25] *= one_m_omega; f[25] += t10_omega * (C1-C2+C3);

    // i=13 and i=26
    C2 = -kx +ky +kz;
    C3 = -kxky_ - kxkz_ + kykz_;
    f[13] *= one_m_omega; f[13] += t10_omega * (C1+C2+C3);
    f[26] *= one_m_omega; f[26] += t10_omega * (C1-C2+C3);

    return invRho*invRho*jSqr;
}

static T bgk_ma2_collision(Array<T,Descriptor::q>& f, T rhoBar, Array<T,3> const& j, T omega) {
    T invRho = Descriptor::invRho(rhoBar);
    const T jSqr = VectorTemplateImpl<T,Descriptor::d>::normSqr(j);
    for (plint iPop=0; iPop < Descriptor::q; ++iPop) {
        f[iPop] *= (T)1-omega;
        f[iPop] += omega * dynamicsTemplatesImpl<T,Descriptor>::bgk_ma2_equilibrium (
                                iPop, rhoBar, invRho, j, jSqr );
    }
    return jSqr*invRho*invRho;
    //return bgk_ma2_collision_base(f, rhoBar, j, omega, Descriptor::invRho(rhoBar));
}

static T bgk_inc_collision(Array<T,Descriptor::q>& f, T rhoBar, Array<T,3> const& j, T omega, T invRho0=(T)1 ) {
    return bgk_ma2_collision_base(f, rhoBar, j, omega, invRho0);
}

static T rlb_collision(Array<T,Descriptor::q>& f, T rhoBar, T invRho, Array<T,3> const& j, Array<T,6> const& PiNeq, T omega ) {
    typedef dynamicsTemplatesImpl<T, descriptors::D3Q27DescriptorBase<T> > DH;
    const T jSqr = j[0]*j[0] + j[1]*j[1] + j[2]*j[2];

    T piNeq0 = neqPiD3Q27<T>::fromPiToFneq0(PiNeq);
    T piNeq1 = neqPiD3Q27<T>::fromPiToFneq1(PiNeq);
    T piNeq2 = neqPiD3Q27<T>::fromPiToFneq2(PiNeq);
    T piNeq3 = neqPiD3Q27<T>::fromPiToFneq3(PiNeq);
    T piNeq4 = neqPiD3Q27<T>::fromPiToFneq4(PiNeq);
    T piNeq5 = neqPiD3Q27<T>::fromPiToFneq5(PiNeq);
    T piNeq6 = neqPiD3Q27<T>::fromPiToFneq6(PiNeq);
    T piNeq7 = neqPiD3Q27<T>::fromPiToFneq7(PiNeq);
    T piNeq8 = neqPiD3Q27<T>::fromPiToFneq8(PiNeq);
    T piNeq9 = neqPiD3Q27<T>::fromPiToFneq9(PiNeq);
    T piNeq10 = neqPiD3Q27<T>::fromPiToFneq10(PiNeq);
    T piNeq11 = neqPiD3Q27<T>::fromPiToFneq11(PiNeq);
    T piNeq12 = neqPiD3Q27<T>::fromPiToFneq12(PiNeq);
    T piNeq13 = neqPiD3Q27<T>::fromPiToFneq13(PiNeq);

    f[0]  = DH::bgk_ma2_equilibrium(0, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq0;
    f[1]  = DH::bgk_ma2_equilibrium(1, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq1;
    f[2]  = DH::bgk_ma2_equilibrium(2, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq2;
    f[3]  = DH::bgk_ma2_equilibrium(3, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq3;
    f[4]  = DH::bgk_ma2_equilibrium(4, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq4;
    f[5]  = DH::bgk_ma2_equilibrium(5, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq5;
    f[6]  = DH::bgk_ma2_equilibrium(6, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq6;
    f[7]  = DH::bgk_ma2_equilibrium(7, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq7;
    f[8]  = DH::bgk_ma2_equilibrium(8, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq8;
    f[9]  = DH::bgk_ma2_equilibrium(9, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq9;
    f[10]  = DH::bgk_ma2_equilibrium(10, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq10;
    f[11]  = DH::bgk_ma2_equilibrium(11, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq11;
    f[12]  = DH::bgk_ma2_equilibrium(12, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq12;
    f[13]  = DH::bgk_ma2_equilibrium(13, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq13;
    f[14]  = DH::bgk_ma2_equilibrium(14, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq1;
    f[15]  = DH::bgk_ma2_equilibrium(15, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq2;
    f[16]  = DH::bgk_ma2_equilibrium(16, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq3;
    f[17]  = DH::bgk_ma2_equilibrium(17, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq4;
    f[18]  = DH::bgk_ma2_equilibrium(18, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq5;
    f[19]  = DH::bgk_ma2_equilibrium(19, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq6;
    f[20]  = DH::bgk_ma2_equilibrium(20, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq7;
    f[21]  = DH::bgk_ma2_equilibrium(21, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq8;
    f[22]  = DH::bgk_ma2_equilibrium(22, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq9;
    f[23]  = DH::bgk_ma2_equilibrium(23, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq10;
    f[24]  = DH::bgk_ma2_equilibrium(24, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq11;
    f[25]  = DH::bgk_ma2_equilibrium(25, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq12;
    f[26]  = DH::bgk_ma2_equilibrium(26, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq13;

    return jSqr*invRho*invRho;
}

static T bgk_ma2_constRho_collision(Array<T,Descriptor::q>& f, T rhoBar, Array<T,Descriptor::d> const& j, T ratioRho, T omega) {
    T invRho = Descriptor::invRho(rhoBar);
    const T jSqr = j[0]*j[0] + j[1]*j[1] + j[2]*j[2];
    for (plint iPop=0; iPop < Descriptor::q; ++iPop) {
        T feq = bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr );
        f[iPop] =
          ratioRho*feq + Descriptor::t[iPop]*(ratioRho-(T)1) +
          ((T)1-omega)*(f[iPop]-feq);
    }
    return jSqr*invRho*invRho;
}

static T precond_bgk_ma2_equilibrium(plint iPop, T rhoBar, T invRho, Array<T,3> const& j, T jSqr, T invGamma ) {
    T c_j = Descriptor::c[iPop][0]*j[0] + Descriptor::c[iPop][1]*j[1] + Descriptor::c[iPop][2]*j[2];
    return Descriptor::t[iPop] * ( rhoBar + 3.*c_j + invGamma*invRho*(4.5*c_j*c_j - 1.5*jSqr) );
}

static T precond_bgk_ma2_collision_base (
        Array<T,Descriptor::q>& f, T rhoBar, Array<T,3> const& j, T omega, T invGamma, bool incompressible )
{
    T invRho = incompressible ? (T)1 : Descriptor::invRho(rhoBar);
    T one_m_omega = (T)1 - omega;
    T t0_omega = Descriptor::t[0] * omega;
    T t1_omega = Descriptor::t[1] * omega;
    T t4_omega = Descriptor::t[4] * omega;
    T t10_omega = Descriptor::t[10] * omega;

    T jSqr   = j[0]*j[0] + j[1]*j[1] + j[2]*j[2];
    T kx     = (T)3 * j[0];
    T ky     = (T)3 * j[1];
    T kz     = (T)3 * j[2];
    T kxSqr_ = invGamma* invRho / (T)2 * kx*kx;
    T kySqr_ = invGamma* invRho / (T)2 * ky*ky;
    T kzSqr_ = invGamma* invRho / (T)2 * kz*kz;
    T kxky_  = invGamma* invRho * kx*ky;
    T kxkz_  = invGamma* invRho * kx*kz;
    T kykz_  = invGamma* invRho * ky*kz;

    T C1 = rhoBar + invGamma*invRho*(T)3*jSqr;
    T C2, C3;

    // i=0
    C3 = -kxSqr_ - kySqr_ - kzSqr_;
    f[0] *= one_m_omega; f[0] += t0_omega * (C1+C3);

    // i=1 and i=14
    C2 = -kx;
    C3 = -kySqr_ - kzSqr_;
    f[1]  *= one_m_omega; f[1]  += t1_omega * (C1+C2+C3);
    f[14] *= one_m_omega; f[14] += t1_omega * (C1-C2+C3);

    // i=2 and i=15
    C2 = -ky;
    C3 = -kxSqr_ - kzSqr_;
    f[2]  *= one_m_omega; f[2]  += t1_omega * (C1+C2+C3);
    f[15] *= one_m_omega; f[15] += t1_omega * (C1-C2+C3);

    // i=3 and i=16
    C2 = -kz;
    C3 = -kxSqr_ - kySqr_;
    f[3]  *= one_m_omega; f[3]  += t1_omega * (C1+C2+C3);
    f[16] *= one_m_omega; f[16] += t1_omega * (C1-C2+C3);

    // i=4 and i=17
    C2 = -kx - ky;
    C3 = kxky_ - kzSqr_;
    f[4]  *= one_m_omega; f[4]  += t4_omega * (C1+C2+C3);
    f[17] *= one_m_omega; f[17] += t4_omega * (C1-C2+C3);

    // i=5 and i=18
    C2 = -kx + ky;
    C3 = -kxky_ - kzSqr_;
    f[5]  *= one_m_omega; f[5]  += t4_omega * (C1+C2+C3);
    f[18] *= one_m_omega; f[18] += t4_omega * (C1-C2+C3);

    // i=6 and i=19
    C2 = -kx - kz;
    C3 = kxkz_ - kySqr_;
    f[6]  *= one_m_omega; f[6]  += t4_omega * (C1+C2+C3);
    f[19] *= one_m_omega; f[19] += t4_omega * (C1-C2+C3);

    // i=7 and i=20
    C2 = -kx + kz;
    C3 = -kxkz_ - kySqr_;
    f[7]  *= one_m_omega; f[7]  += t4_omega * (C1+C2+C3);
    f[20] *= one_m_omega; f[20] += t4_omega * (C1-C2+C3);

    // i=8 and i=21
    C2 = -ky - kz;
    C3 = kykz_ - kxSqr_;
    f[8]  *= one_m_omega; f[8]  += t4_omega * (C1+C2+C3);
    f[21] *= one_m_omega; f[21] += t4_omega * (C1-C2+C3);

    // i=9 and i=22
    C2 = -ky + kz;
    C3 = -kykz_ - kxSqr_;
    f[9]  *= one_m_omega; f[9]  += t4_omega * (C1+C2+C3);
    f[22] *= one_m_omega; f[22] += t4_omega * (C1-C2+C3);

    // i=10 and i=23
    C2 = -kx -ky -kz;
    C3 = kxky_ + kxkz_ + kykz_;
    f[10] *= one_m_omega; f[10] += t10_omega * (C1+C2+C3);
    f[23] *= one_m_omega; f[23] += t10_omega * (C1-C2+C3);

    // i=11 and i=24
    C2 = -kx -ky +kz;
    C3 = kxky_ - kxkz_ - kykz_;
    f[11] *= one_m_omega; f[11] += t10_omega * (C1+C2+C3);
    f[24] *= one_m_omega; f[24] += t10_omega * (C1-C2+C3);

    // i=12 and i=25
    C2 = -kx +ky -kz;
    C3 = -kxky_ + kxkz_ - kykz_;
    f[12] *= one_m_omega; f[12] += t10_omega * (C1+C2+C3);
    f[25] *= one_m_omega; f[25] += t10_omega * (C1-C2+C3);

    // i=13 and i=26
    C2 = -kx +ky +kz;
    C3 = -kxky_ - kxkz_ + kykz_;
    f[13] *= one_m_omega; f[13] += t10_omega * (C1+C2+C3);
    f[26] *= one_m_omega; f[26] += t10_omega * (C1-C2+C3);

    return invRho*invRho*jSqr;
}

static T precond_bgk_ma2_collision(Array<T,Descriptor::q>& f, T rhoBar, Array<T,3> const& j, T omega, T invGamma) {
    T invRho = Descriptor::invRho(rhoBar);
    const T jSqr = VectorTemplateImpl<T,Descriptor::d>::normSqr(j);
    for (plint iPop=0; iPop < Descriptor::q; ++iPop) {
        f[iPop] *= (T)1-omega;
        f[iPop] += omega * dynamicsTemplatesImpl<T,Descriptor>::precond_bgk_ma2_equilibrium (
                                iPop, rhoBar, invRho, j, jSqr, invGamma );
    }
    return jSqr*invRho*invRho;
    //return bgk_ma2_collision_base(f, rhoBar, j, omega, false);
}

};  //struct dynamicsTemplatesImpl<D3Q27DescriptorBase>


template<typename T>
struct neqPiD3Q19 {

    typedef SymmetricTensorImpl<T,3> S;

    static T fromPiToFneq0(Array<T,6> const& pi) {
        return (T)3/(T)2 * (
                 -(T)1/(T)3*pi[S::xx] - (T)1/(T)3*pi[S::yy] - (T)1/(T)3*pi[S::zz]
               );
    }

    static T fromPiToFneq1(Array<T,6> const& pi) {
        return (T)1/(T)4 * (
                  (T)2/(T)3*pi[S::xx] - (T)1/(T)3*pi[S::yy] - (T)1/(T)3*pi[S::zz]
               );
    }

    static T fromPiToFneq2(Array<T,6> const& pi) {
        return (T)1/(T)4 * (
                 -(T)1/(T)3*pi[S::xx] + (T)2/(T)3*pi[S::yy] - (T)1/(T)3*pi[S::zz]
               );
    }

    static T fromPiToFneq3(Array<T,6> const& pi) {
        return (T)1/(T)4 * (
                 -(T)1/(T)3*pi[S::xx] - (T)1/(T)3*pi[S::yy] + (T)2/(T)3*pi[S::zz]
               );
    }

    static T fromPiToFneq4(Array<T,6> const& pi) {
        return (T)1/(T)8 * (
                  (T)2/(T)3*pi[S::xx] + (T)2/(T)3*pi[S::yy] - (T)1/(T)3*pi[S::zz]
                  + (T)2*pi[S::xy]
               );
    }

    static T fromPiToFneq5(Array<T,6> const& pi) {
        return (T)1/(T)8 * (
                  (T)2/(T)3*pi[S::xx] + (T)2/(T)3*pi[S::yy] - (T)1/(T)3*pi[S::zz]
                  - (T)2*pi[S::xy]
               );
    }

    static T fromPiToFneq6(Array<T,6> const& pi) {
        return (T)1/(T)8 * (
                  (T)2/(T)3*pi[S::xx] - (T)1/(T)3*pi[S::yy] + (T)2/(T)3*pi[S::zz]
                  + (T)2*pi[S::xz]
               );
    }

    static T fromPiToFneq7(Array<T,6> const& pi) {
        return (T)1/(T)8 * (
                  (T)2/(T)3*pi[S::xx] - (T)1/(T)3*pi[S::yy] + (T)2/(T)3*pi[S::zz]
                  - (T)2*pi[S::xz]
               );
    }

    static T fromPiToFneq8(Array<T,6> const& pi) {
        return (T)1/(T)8 * (
                 -(T)1/(T)3*pi[S::xx] + (T)2/(T)3*pi[S::yy] + (T)2/(T)3*pi[S::zz]
                  + (T)2*pi[S::yz]
               );
    }

    static T fromPiToFneq9(Array<T,6> const& pi) {
        return (T)1/(T)8 * (
                 -(T)1/(T)3*pi[S::xx] + (T)2/(T)3*pi[S::yy] + (T)2/(T)3*pi[S::zz]
                  - (T)2*pi[S::yz]
               );
    }

};  // struct neqPiD3Q19



// Efficient specialization for D3Q19 lattice
template<typename T>
struct dynamicsTemplatesImpl<T, descriptors::D3Q19DescriptorBase<T> > {

typedef descriptors::D3Q19DescriptorBase<T> Descriptor;

static T bgk_ma2_equilibrium(plint iPop, T rhoBar, T invRho, Array<T,3> const& j, T jSqr ) {
    T c_j = Descriptor::c[iPop][0]*j[0] + Descriptor::c[iPop][1]*j[1] + Descriptor::c[iPop][2]*j[2];
    return Descriptor::t[iPop] * ( rhoBar + 3.*c_j + invRho*(4.5*c_j*c_j - 1.5*jSqr) );
}

static void bgk_ma2_equilibria( T rhoBar, T invRho, Array<T,Descriptor::d> const& j,
                                T jSqr, Array<T,Descriptor::q>& eqPop )
{
    T t0 = Descriptor::t[0];
    T t1 = Descriptor::t[1];
    T t4 = Descriptor::t[4];

    T kx     = (T)3 * j[0];
    T ky     = (T)3 * j[1];
    T kz     = (T)3 * j[2];
    T kxSqr_ = invRho / (T)2 * kx*kx;
    T kySqr_ = invRho / (T)2 * ky*ky;
    T kzSqr_ = invRho / (T)2 * kz*kz;
    T kxky_  = invRho * kx*ky;
    T kxkz_  = invRho * kx*kz;
    T kykz_  = invRho * ky*kz;

    T C1 = rhoBar + invRho*(T)3*jSqr;
    T C2, C3;

    // i=0
    C3 = -kxSqr_ - kySqr_ - kzSqr_;
    eqPop[0] = t0 * (C1+C3);

    // i=1 and i=10
    C2 = -kx;
    C3 = -kySqr_ - kzSqr_;
    eqPop[1]  = t1 * (C1+C2+C3);
    eqPop[10] = t1 * (C1-C2+C3);

    // i=2 and i=11
    C2 = -ky;
    C3 = -kxSqr_ - kzSqr_;
    eqPop[2]  = t1 * (C1+C2+C3);
    eqPop[11] = t1 * (C1-C2+C3);

    // i=3 and i=12
    C2 = -kz;
    C3 = -kxSqr_ - kySqr_;
    eqPop[3]  = t1 * (C1+C2+C3);
    eqPop[12] = t1 * (C1-C2+C3);

    // i=4 and i=13
    C2 = -kx - ky;
    C3 = kxky_ - kzSqr_;
    eqPop[4]  = t4 * (C1+C2+C3);
    eqPop[13] = t4 * (C1-C2+C3);

    // i=5 and i=14
    C2 = -kx + ky;
    C3 = -kxky_ - kzSqr_;
    eqPop[5]  = t4 * (C1+C2+C3);
    eqPop[14] = t4 * (C1-C2+C3);

    // i=6 and i=15
    C2 = -kx - kz;
    C3 = kxkz_ - kySqr_;
    eqPop[6]  = t4 * (C1+C2+C3);
    eqPop[15] = t4 * (C1-C2+C3);

    // i=7 and i=16
    C2 = -kx + kz;
    C3 = -kxkz_ - kySqr_;
    eqPop[7]  = t4 * (C1+C2+C3);
    eqPop[16] = t4 * (C1-C2+C3);

    // i=8 and i=17
    C2 = -ky - kz;
    C3 = kykz_ - kxSqr_;
    eqPop[8]  = t4 * (C1+C2+C3);
    eqPop[17] = t4 * (C1-C2+C3);

    // i=9 and i=18
    C2 = -ky + kz;
    C3 = -kykz_ - kxSqr_;
    eqPop[9]  = t4 * (C1+C2+C3);
    eqPop[18] = t4 * (C1-C2+C3);
}

static T bgk_ma2_collision_base(Array<T,Descriptor::q>& f, T rhoBar, Array<T,3> const& j, T omega, T invRho) {
    T one_m_omega = (T)1 - omega;
    T t0_omega = Descriptor::t[0] * omega;
    T t1_omega = Descriptor::t[1] * omega;
    T t4_omega = Descriptor::t[4] * omega;

    T jSqr   = j[0]*j[0] + j[1]*j[1] + j[2]*j[2];
    T kx     = (T)3 * j[0];
    T ky     = (T)3 * j[1];
    T kz     = (T)3 * j[2];
    T kxSqr_ = invRho / (T)2 * kx*kx;
    T kySqr_ = invRho / (T)2 * ky*ky;
    T kzSqr_ = invRho / (T)2 * kz*kz;
    T kxky_  = invRho * kx*ky;
    T kxkz_  = invRho * kx*kz;
    T kykz_  = invRho * ky*kz;

    T C1 = rhoBar + invRho*(T)3*jSqr;
    T C2, C3;

    // i=0
    C3 = -kxSqr_ - kySqr_ - kzSqr_;
    f[0] *= one_m_omega; f[0] += t0_omega * (C1+C3);

    // i=1 and i=10
    C2 = -kx;
    C3 = -kySqr_ - kzSqr_;
    f[1]  *= one_m_omega; f[1]  += t1_omega * (C1+C2+C3);
    f[10] *= one_m_omega; f[10] += t1_omega * (C1-C2+C3);

    // i=2 and i=11
    C2 = -ky;
    C3 = -kxSqr_ - kzSqr_;
    f[2]  *= one_m_omega; f[2]  += t1_omega * (C1+C2+C3);
    f[11] *= one_m_omega; f[11] += t1_omega * (C1-C2+C3);

    // i=3 and i=12
    C2 = -kz;
    C3 = -kxSqr_ - kySqr_;
    f[3]  *= one_m_omega; f[3]  += t1_omega * (C1+C2+C3);
    f[12] *= one_m_omega; f[12] += t1_omega * (C1-C2+C3);

    // i=4 and i=13
    C2 = -kx - ky;
    C3 = kxky_ - kzSqr_;
    f[4]  *= one_m_omega; f[4]  += t4_omega * (C1+C2+C3);
    f[13] *= one_m_omega; f[13] += t4_omega * (C1-C2+C3);

    // i=5 and i=14
    C2 = -kx + ky;
    C3 = -kxky_ - kzSqr_;
    f[5]  *= one_m_omega; f[5]  += t4_omega * (C1+C2+C3);
    f[14] *= one_m_omega; f[14] += t4_omega * (C1-C2+C3);

    // i=6 and i=15
    C2 = -kx - kz;
    C3 = kxkz_ - kySqr_;
    f[6]  *= one_m_omega; f[6]  += t4_omega * (C1+C2+C3);
    f[15] *= one_m_omega; f[15] += t4_omega * (C1-C2+C3);

    // i=7 and i=16
    C2 = -kx + kz;
    C3 = -kxkz_ - kySqr_;
    f[7]  *= one_m_omega; f[7]  += t4_omega * (C1+C2+C3);
    f[16] *= one_m_omega; f[16] += t4_omega * (C1-C2+C3);

    // i=8 and i=17
    C2 = -ky - kz;
    C3 = kykz_ - kxSqr_;
    f[8]  *= one_m_omega; f[8]  += t4_omega * (C1+C2+C3);
    f[17] *= one_m_omega; f[17] += t4_omega * (C1-C2+C3);

    // i=9 and i=18
    C2 = -ky + kz;
    C3 = -kykz_ - kxSqr_;
    f[9]  *= one_m_omega; f[9]  += t4_omega * (C1+C2+C3);
    f[18] *= one_m_omega; f[18] += t4_omega * (C1-C2+C3);

    return invRho*invRho*jSqr;
}

static T bgk_ma2_collision(Array<T,Descriptor::q>& f, T rhoBar, Array<T,3> const& j, T omega) {
    return bgk_ma2_collision_base(f, rhoBar, j, omega, Descriptor::invRho(rhoBar));
}

static T bgk_inc_collision(Array<T,Descriptor::q>& f, T rhoBar, Array<T,3> const& j, T omega, T invRho0=(T)1 ) {
    return bgk_ma2_collision_base(f, rhoBar, j, omega, invRho0);
}

static T rlb_collision(Array<T,Descriptor::q>& f, T rhoBar, T invRho, Array<T,3> const& j, Array<T,6> const& PiNeq, T omega ) {
    typedef dynamicsTemplatesImpl<T, descriptors::D3Q19DescriptorBase<T> > DH;
    const T jSqr = j[0]*j[0] + j[1]*j[1] + j[2]*j[2];

    T piNeq0 = neqPiD3Q19<T>::fromPiToFneq0(PiNeq);
    T piNeq1 = neqPiD3Q19<T>::fromPiToFneq1(PiNeq);
    T piNeq2 = neqPiD3Q19<T>::fromPiToFneq2(PiNeq);
    T piNeq3 = neqPiD3Q19<T>::fromPiToFneq3(PiNeq);
    T piNeq4 = neqPiD3Q19<T>::fromPiToFneq4(PiNeq);
    T piNeq5 = neqPiD3Q19<T>::fromPiToFneq5(PiNeq);
    T piNeq6 = neqPiD3Q19<T>::fromPiToFneq6(PiNeq);
    T piNeq7 = neqPiD3Q19<T>::fromPiToFneq7(PiNeq);
    T piNeq8 = neqPiD3Q19<T>::fromPiToFneq8(PiNeq);
    T piNeq9 = neqPiD3Q19<T>::fromPiToFneq9(PiNeq);

    f[0]  = DH::bgk_ma2_equilibrium(0, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq0;
    f[1]  = DH::bgk_ma2_equilibrium(1, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq1;
    f[2]  = DH::bgk_ma2_equilibrium(2, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq2;
    f[3]  = DH::bgk_ma2_equilibrium(3, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq3;
    f[4]  = DH::bgk_ma2_equilibrium(4, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq4;
    f[5]  = DH::bgk_ma2_equilibrium(5, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq5;
    f[6]  = DH::bgk_ma2_equilibrium(6, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq6;
    f[7]  = DH::bgk_ma2_equilibrium(7, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq7;
    f[8]  = DH::bgk_ma2_equilibrium(8, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq8;
    f[9]  = DH::bgk_ma2_equilibrium(9, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq9;
    f[10]  = DH::bgk_ma2_equilibrium(10, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq1;
    f[11]  = DH::bgk_ma2_equilibrium(11, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq2;
    f[12]  = DH::bgk_ma2_equilibrium(12, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq3;
    f[13]  = DH::bgk_ma2_equilibrium(13, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq4;
    f[14]  = DH::bgk_ma2_equilibrium(14, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq5;
    f[15]  = DH::bgk_ma2_equilibrium(15, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq6;
    f[16]  = DH::bgk_ma2_equilibrium(16, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq7;
    f[17]  = DH::bgk_ma2_equilibrium(17, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq8;
    f[18]  = DH::bgk_ma2_equilibrium(18, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq9;
    return jSqr*invRho*invRho;
}

static T bgk_ma2_constRho_collision (
        Array<T,Descriptor::q>& f, T rhoBar, Array<T,Descriptor::d> const& j, T ratioRho, T omega )
{
    T invRho = Descriptor::invRho(rhoBar);
    const T jSqr = j[0]*j[0] + j[1]*j[1] + j[2]*j[2];
    for (plint iPop=0; iPop < Descriptor::q; ++iPop) {
        T feq = bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr );
        f[iPop] =
          ratioRho*feq + Descriptor::t[iPop]*(ratioRho-(T)1) +
          ((T)1-omega)*(f[iPop]-feq);
    }
    return jSqr*invRho*invRho;
}


static T precond_bgk_ma2_equilibrium(plint iPop, T rhoBar, T invRho, Array<T,3> const& j, T jSqr, T invGamma ) {
    T c_j = Descriptor::c[iPop][0]*j[0] + Descriptor::c[iPop][1]*j[1] + Descriptor::c[iPop][2]*j[2];
    return Descriptor::t[iPop] * ( rhoBar + 3.*c_j + invGamma*invRho*(4.5*c_j*c_j - 1.5*jSqr) );
}

static T precond_bgk_ma2_collision_base( Array<T,Descriptor::q>& f, T rhoBar, Array<T,3> const& j, T omega,
                                         T invGamma, bool incompressible )
{
    T invRho = incompressible ? (T)1 : Descriptor::invRho(rhoBar);
    T one_m_omega = (T)1 - omega;
    T t0_omega = Descriptor::t[0] * omega;
    T t1_omega = Descriptor::t[1] * omega;
    T t4_omega = Descriptor::t[4] * omega;

    T jSqr   = j[0]*j[0] + j[1]*j[1] + j[2]*j[2];
    T kx     = (T)3 * j[0];
    T ky     = (T)3 * j[1];
    T kz     = (T)3 * j[2];
    T kxSqr_ = invGamma* invRho / (T)2 * kx*kx;
    T kySqr_ = invGamma* invRho / (T)2 * ky*ky;
    T kzSqr_ = invGamma* invRho / (T)2 * kz*kz;
    T kxky_  = invGamma* invRho * kx*ky;
    T kxkz_  = invGamma* invRho * kx*kz;
    T kykz_  = invGamma* invRho * ky*kz;

    T C1 = rhoBar + invGamma*invRho*(T)3*jSqr;
    T C2, C3;

    // i=0
    C3 = -kxSqr_ - kySqr_ - kzSqr_;
    f[0] *= one_m_omega; f[0] += t0_omega * (C1+C3);

    // i=1 and i=10
    C2 = -kx;
    C3 = -kySqr_ - kzSqr_;
    f[1]  *= one_m_omega; f[1]  += t1_omega * (C1+C2+C3);
    f[10] *= one_m_omega; f[10] += t1_omega * (C1-C2+C3);

    // i=2 and i=11
    C2 = -ky;
    C3 = -kxSqr_ - kzSqr_;
    f[2]  *= one_m_omega; f[2]  += t1_omega * (C1+C2+C3);
    f[11] *= one_m_omega; f[11] += t1_omega * (C1-C2+C3);

    // i=3 and i=12
    C2 = -kz;
    C3 = -kxSqr_ - kySqr_;
    f[3]  *= one_m_omega; f[3]  += t1_omega * (C1+C2+C3);
    f[12] *= one_m_omega; f[12] += t1_omega * (C1-C2+C3);

    // i=4 and i=13
    C2 = -kx - ky;
    C3 = kxky_ - kzSqr_;
    f[4]  *= one_m_omega; f[4]  += t4_omega * (C1+C2+C3);
    f[13] *= one_m_omega; f[13] += t4_omega * (C1-C2+C3);

    // i=5 and i=14
    C2 = -kx + ky;
    C3 = -kxky_ - kzSqr_;
    f[5]  *= one_m_omega; f[5]  += t4_omega * (C1+C2+C3);
    f[14] *= one_m_omega; f[14] += t4_omega * (C1-C2+C3);

    // i=6 and i=15
    C2 = -kx - kz;
    C3 = kxkz_ - kySqr_;
    f[6]  *= one_m_omega; f[6]  += t4_omega * (C1+C2+C3);
    f[15] *= one_m_omega; f[15] += t4_omega * (C1-C2+C3);

    // i=7 and i=16
    C2 = -kx + kz;
    C3 = -kxkz_ - kySqr_;
    f[7]  *= one_m_omega; f[7]  += t4_omega * (C1+C2+C3);
    f[16] *= one_m_omega; f[16] += t4_omega * (C1-C2+C3);

    // i=8 and i=17
    C2 = -ky - kz;
    C3 = kykz_ - kxSqr_;
    f[8]  *= one_m_omega; f[8]  += t4_omega * (C1+C2+C3);
    f[17] *= one_m_omega; f[17] += t4_omega * (C1-C2+C3);

    // i=9 and i=18
    C2 = -ky + kz;
    C3 = -kykz_ - kxSqr_;
    f[9]  *= one_m_omega; f[9]  += t4_omega * (C1+C2+C3);
    f[18] *= one_m_omega; f[18] += t4_omega * (C1-C2+C3);

    return invRho*invRho*jSqr;
}

static T precond_bgk_ma2_collision(Array<T,Descriptor::q>& f, T rhoBar, Array<T,3> const& j, T omega, T invGamma)
{
    return precond_bgk_ma2_collision_base(f, rhoBar, j, omega, invGamma, false);
}

};  //struct dynamicsTemplatesImpl<D3Q19DescriptorBase>


/// Compute Pi tensor efficiently on D3Q15 lattice
template<typename T>
struct neqPiD3Q15 {

    typedef SymmetricTensorImpl<T,3> S;

    static T fromPiToFneq0(Array<T,6> const& pi) {
        return -(T)1/(T)3*pi[S::xx] - (T)1/(T)3*pi[S::yy] - (T)1/(T)3*pi[S::zz];
    }

    static T fromPiToFneq1(Array<T,6> const& pi) {
        return (T)1/(T)2 * (
                  (T)2/(T)3*pi[S::xx] - (T)1/(T)3*pi[S::yy] - (T)1/(T)3*pi[S::zz]
               );
    }

    static T fromPiToFneq2(Array<T,6> const& pi) {
        return (T)1/(T)2 * (
                 -(T)1/(T)3*pi[S::xx] + (T)2/(T)3*pi[S::yy] - (T)1/(T)3*pi[S::zz]
               );
    }

    static T fromPiToFneq3(Array<T,6> const& pi) {
        return (T)1/(T)2 * (
                 -(T)1/(T)3*pi[S::xx] - (T)1/(T)3*pi[S::yy] + (T)2/(T)3*pi[S::zz]
               );
    }

    static T fromPiToFneq4(Array<T,6> const& pi) {
        return (T)1/(T)16 * (
                  (T)2/(T)3*pi[S::xx] + (T)2/(T)3*pi[S::yy] + (T)2/(T)3*pi[S::zz]
                  + pi[S::xy] + pi[S::xz] + pi[S::yz]
               );
    }

    static T fromPiToFneq5(Array<T,6> const& pi) {
        return (T)1/(T)16 * (
                  (T)2/(T)3*pi[S::xx] + (T)2/(T)3*pi[S::yy] + (T)2/(T)3*pi[S::zz]
                  + pi[S::xy] - pi[S::xz] - pi[S::yz]
               );
    }

    static T fromPiToFneq6(Array<T,6> const& pi) {
        return (T)1/(T)16 * (
                  (T)2/(T)3*pi[S::xx] + (T)2/(T)3*pi[S::yy] + (T)2/(T)3*pi[S::zz]
                  - pi[S::xy] + pi[S::xz] - pi[S::yz]
               );
    }

    static T fromPiToFneq7(Array<T,6> const& pi) {
        return (T)1/(T)16 * (
                  (T)2/(T)3*pi[S::xx] + (T)2/(T)3*pi[S::yy] + (T)2/(T)3*pi[S::zz]
                  - pi[S::xy] - pi[S::xz] + pi[S::yz]
               );
    }

};  // struct neqPiD3Q15


// Efficient specialization for D3Q15 lattice
template<typename T>
struct dynamicsTemplatesImpl<T, descriptors::D3Q15DescriptorBase<T> > {

typedef descriptors::D3Q15DescriptorBase<T> Descriptor;

static T bgk_ma2_equilibrium(plint iPop, T rhoBar, T invRho, Array<T,3> const& j, T jSqr ) {
    T c_j = Descriptor::c[iPop][0]*j[0] + Descriptor::c[iPop][1]*j[1] + Descriptor::c[iPop][2]*j[2];
    return Descriptor::t[iPop] * ( rhoBar + 3.*c_j + invRho*(4.5*c_j*c_j - 1.5*jSqr) );
}

static void bgk_ma2_equilibria( T rhoBar, T invRho, Array<T,Descriptor::d> const& j,
                                T jSqr, Array<T,Descriptor::q>& eqPop )
{
    T t0 = Descriptor::t[0];
    T t1 = Descriptor::t[1];
    T t4 = Descriptor::t[4];

    T kx     = (T)3 * j[0];
    T ky     = (T)3 * j[1];
    T kz     = (T)3 * j[2];
    T kxSqr_ = invRho / (T)2 * kx*kx;
    T kySqr_ = invRho / (T)2 * ky*ky;
    T kzSqr_ = invRho / (T)2 * kz*kz;
    T kxky_  = invRho * kx*ky;
    T kxkz_  = invRho * kx*kz;
    T kykz_  = invRho * ky*kz;

    T C1 = rhoBar + invRho*(T)3*jSqr;
    T C2, C3;

    // i=0
    C3 = -kxSqr_ - kySqr_ - kzSqr_;
    eqPop[0] = t0 * (C1+C3);

    // i=1 and i=8
    C2 = -kx;
    C3 = -kySqr_ - kzSqr_;
    eqPop[1]  = t1 * (C1+C2+C3);
    eqPop[8]  = t1 * (C1-C2+C3);

    // i=2 and i=9
    C2 = -ky;
    C3 = -kxSqr_ - kzSqr_;
    eqPop[2]  = t1 * (C1+C2+C3);
    eqPop[9]  = t1 * (C1-C2+C3);

    // i=3 and i=10
    C2 = -kz;
    C3 = -kxSqr_ - kySqr_;
    eqPop[3]  = t1 * (C1+C2+C3);
    eqPop[10] = t1 * (C1-C2+C3);

    // i=4 and i=11
    C2 = -kx -ky -kz;
    C3 = kxky_ + kxkz_ + kykz_;
    eqPop[4]  = t4 * (C1+C2+C3);
    eqPop[11] = t4 * (C1-C2+C3);

    // i=5 and i=12
    C2 = -kx -ky +kz;
    C3 = kxky_ - kxkz_ - kykz_;
    eqPop[5]  = t4 * (C1+C2+C3);
    eqPop[12] = t4 * (C1-C2+C3);

    // i=6 and i=13
    C2 = -kx +ky -kz;
    C3 = -kxky_ + kxkz_ - kykz_;
    eqPop[6]  = t4 * (C1+C2+C3);
    eqPop[13] = t4 * (C1-C2+C3);

    // i=7 and i=14
    C2 = -kx +ky +kz;
    C3 = -kxky_ - kxkz_ + kykz_;
    eqPop[7]  = t4 * (C1+C2+C3);
    eqPop[14] = t4 * (C1-C2+C3);
}

static T bgk_ma2_collision_base(Array<T,Descriptor::q>& f, T rhoBar, Array<T,3> const& j, T omega, T invRho ) {
    T one_m_omega = (T)1 - omega;
    T t0_omega = Descriptor::t[0] * omega;
    T t1_omega = Descriptor::t[1] * omega;
    T t4_omega = Descriptor::t[4] * omega;

    T jSqr   = j[0]*j[0] + j[1]*j[1] + j[2]*j[2];
    T kx     = (T)3 * j[0];
    T ky     = (T)3 * j[1];
    T kz     = (T)3 * j[2];
    T kxSqr_ = invRho / (T)2 * kx*kx;
    T kySqr_ = invRho / (T)2 * ky*ky;
    T kzSqr_ = invRho / (T)2 * kz*kz;
    T kxky_  = invRho * kx*ky;
    T kxkz_  = invRho * kx*kz;
    T kykz_  = invRho * ky*kz;

    T C1 = rhoBar + invRho*(T)3*jSqr;
    T C2, C3;

    // i=0
    C3 = -kxSqr_ - kySqr_ - kzSqr_;
    f[0] *= one_m_omega; f[0] += t0_omega * (C1+C3);

    // i=1 and i=8
    C2 = -kx;
    C3 = -kySqr_ - kzSqr_;
    f[1]  *= one_m_omega; f[1]  += t1_omega * (C1+C2+C3);
    f[8]  *= one_m_omega; f[8]  += t1_omega * (C1-C2+C3);

    // i=2 and i=9
    C2 = -ky;
    C3 = -kxSqr_ - kzSqr_;
    f[2]  *= one_m_omega; f[2]  += t1_omega * (C1+C2+C3);
    f[9]  *= one_m_omega; f[9]  += t1_omega * (C1-C2+C3);

    // i=3 and i=10
    C2 = -kz;
    C3 = -kxSqr_ - kySqr_;
    f[3]  *= one_m_omega; f[3]  += t1_omega * (C1+C2+C3);
    f[10] *= one_m_omega; f[10] += t1_omega * (C1-C2+C3);

    // i=4 and i=11
    C2 = -kx -ky -kz;
    C3 = kxky_ + kxkz_ + kykz_;
    f[4]  *= one_m_omega; f[4]  += t4_omega * (C1+C2+C3);
    f[11] *= one_m_omega; f[11] += t4_omega * (C1-C2+C3);

    // i=5 and i=12
    C2 = -kx -ky +kz;
    C3 = kxky_ - kxkz_ - kykz_;
    f[5]  *= one_m_omega; f[5]  += t4_omega * (C1+C2+C3);
    f[12] *= one_m_omega; f[12] += t4_omega * (C1-C2+C3);

    // i=6 and i=13
    C2 = -kx +ky -kz;
    C3 = -kxky_ + kxkz_ - kykz_;
    f[6]  *= one_m_omega; f[6]  += t4_omega * (C1+C2+C3);
    f[13] *= one_m_omega; f[13] += t4_omega * (C1-C2+C3);

    // i=7 and i=14
    C2 = -kx +ky +kz;
    C3 = -kxky_ - kxkz_ + kykz_;
    f[7]  *= one_m_omega; f[7]  += t4_omega * (C1+C2+C3);
    f[14] *= one_m_omega; f[14] += t4_omega * (C1-C2+C3);

    return invRho*invRho*jSqr;
}



static T bgk_ma2_collision(Array<T,Descriptor::q>& f, T rhoBar, Array<T,3> const& j, T omega) {
    return bgk_ma2_collision_base(f, rhoBar, j, omega, Descriptor::invRho(rhoBar));
}

static T bgk_inc_collision(Array<T,Descriptor::q>& f, T rhoBar, Array<T,3> const& j, T omega, T invRho0=(T)1 ) {
    return bgk_ma2_collision_base(f, rhoBar, j, omega, invRho0);
}

static T rlbCollision(Array<T,Descriptor::q>& f, T rhoBar, T invRho, Array<T,3> const& j, Array<T,6> const& PiNeq, T omega) {
    typedef dynamicsTemplatesImpl<T, descriptors::D3Q15DescriptorBase<T> > DH;
    const T jSqr = j[0]*j[0] + j[1]*j[1] + j[2]*j[2];

    T piNeq0 = neqPiD3Q15<T>::fromPiToFneq0(PiNeq);
    T piNeq1 = neqPiD3Q15<T>::fromPiToFneq1(PiNeq);
    T piNeq2 = neqPiD3Q15<T>::fromPiToFneq2(PiNeq);
    T piNeq3 = neqPiD3Q15<T>::fromPiToFneq3(PiNeq);
    T piNeq4 = neqPiD3Q15<T>::fromPiToFneq4(PiNeq);
    T piNeq5 = neqPiD3Q15<T>::fromPiToFneq5(PiNeq);
    T piNeq6 = neqPiD3Q15<T>::fromPiToFneq6(PiNeq);
    T piNeq7 = neqPiD3Q15<T>::fromPiToFneq7(PiNeq);

    f[0]  = DH::bgk_ma2_equilibrium(0, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq0;
    f[1]  = DH::bgk_ma2_equilibrium(1, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq1;
    f[2]  = DH::bgk_ma2_equilibrium(2, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq2;
    f[3]  = DH::bgk_ma2_equilibrium(3, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq3;
    f[4]  = DH::bgk_ma2_equilibrium(4, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq4;
    f[5]  = DH::bgk_ma2_equilibrium(5, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq5;
    f[6]  = DH::bgk_ma2_equilibrium(6, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq6;
    f[7]  = DH::bgk_ma2_equilibrium(7, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq7;
    f[8]  = DH::bgk_ma2_equilibrium(8, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq1;
    f[9]  = DH::bgk_ma2_equilibrium(9, rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq2;
    f[10]  = DH::bgk_ma2_equilibrium(10,rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq3;
    f[11]  = DH::bgk_ma2_equilibrium(11,rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq4;
    f[12]  = DH::bgk_ma2_equilibrium(12,rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq5;
    f[13]  = DH::bgk_ma2_equilibrium(13,rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq6;
    f[14]  = DH::bgk_ma2_equilibrium(14,rhoBar, invRho, j, jSqr)
                   + ((T)1-omega) * piNeq7;
    return jSqr*invRho*invRho;
}

static T bgk_ma2_constRho_collision(Array<T,Descriptor::q>& f, T rhoBar, Array<T,Descriptor::d> const& j, T ratioRho, T omega) {
    T invRho = Descriptor::invRho(rhoBar);
    const T jSqr = j[0]*j[0] + j[1]*j[1] + j[2]*j[2];
    for (plint iPop=0; iPop < Descriptor::q; ++iPop) {
        T feq = bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr );
        f[iPop] =
          ratioRho*feq + Descriptor::t[iPop]*(ratioRho-(T)1) +
          ((T)1-omega)*(f[iPop]-feq);
    }
    return jSqr*invRho*invRho;
}


static T precond_bgk_ma2_equilibrium(plint iPop, T rhoBar, T invRho, Array<T,3> const& j, T jSqr, T invGamma )
{
    T c_j = Descriptor::c[iPop][0]*j[0] + Descriptor::c[iPop][1]*j[1] + Descriptor::c[iPop][2]*j[2];
    return Descriptor::t[iPop] * ( rhoBar + 3.*c_j + invGamma*invRho*(4.5*c_j*c_j - 1.5*jSqr) );
}

static T precond_bgk_ma2_collision_base (
        Array<T,Descriptor::q>& f, T rhoBar, Array<T,3> const& j, T omega, T invGamma, bool incompressible )
{
    T invRho = incompressible ? (T)1 : Descriptor::invRho(rhoBar);
    T one_m_omega = (T)1 - omega;
    T t0_omega = Descriptor::t[0] * omega;
    T t1_omega = Descriptor::t[1] * omega;
    T t4_omega = Descriptor::t[4] * omega;

    T jSqr   = j[0]*j[0] + j[1]*j[1] + j[2]*j[2];
    T kx     = (T)3 * j[0];
    T ky     = (T)3 * j[1];
    T kz     = (T)3 * j[2];
    T kxSqr_ = invGamma* invRho / (T)2 * kx*kx;
    T kySqr_ = invGamma* invRho / (T)2 * ky*ky;
    T kzSqr_ = invGamma* invRho / (T)2 * kz*kz;
    T kxky_  = invGamma* invRho * kx*ky;
    T kxkz_  = invGamma* invRho * kx*kz;
    T kykz_  = invGamma* invRho * ky*kz;

    T C1 = rhoBar + invGamma*invRho*(T)3*jSqr;
    T C2, C3;

    // i=0
    C3 = -kxSqr_ - kySqr_ - kzSqr_;
    f[0] *= one_m_omega; f[0] += t0_omega * (C1+C3);

    // i=1 and i=8
    C2 = -kx;
    C3 = -kySqr_ - kzSqr_;
    f[1]  *= one_m_omega; f[1]  += t1_omega * (C1+C2+C3);
    f[8]  *= one_m_omega; f[8]  += t1_omega * (C1-C2+C3);

    // i=2 and i=9
    C2 = -ky;
    C3 = -kxSqr_ - kzSqr_;
    f[2]  *= one_m_omega; f[2]  += t1_omega * (C1+C2+C3);
    f[9]  *= one_m_omega; f[9]  += t1_omega * (C1-C2+C3);

    // i=3 and i=10
    C2 = -kz;
    C3 = -kxSqr_ - kySqr_;
    f[3]  *= one_m_omega; f[3]  += t1_omega * (C1+C2+C3);
    f[10] *= one_m_omega; f[10] += t1_omega * (C1-C2+C3);

    // i=4 and i=11
    C2 = -kx -ky -kz;
    C3 = kxky_ + kxkz_ + kykz_;
    f[4]  *= one_m_omega; f[4]  += t4_omega * (C1+C2+C3);
    f[11] *= one_m_omega; f[11] += t4_omega * (C1-C2+C3);

    // i=5 and i=12
    C2 = -kx -ky +kz;
    C3 = kxky_ - kxkz_ - kykz_;
    f[5]  *= one_m_omega; f[5]  += t4_omega * (C1+C2+C3);
    f[12] *= one_m_omega; f[12] += t4_omega * (C1-C2+C3);

    // i=6 and i=13
    C2 = -kx +ky -kz;
    C3 = -kxky_ + kxkz_ - kykz_;
    f[6]  *= one_m_omega; f[6]  += t4_omega * (C1+C2+C3);
    f[13] *= one_m_omega; f[13] += t4_omega * (C1-C2+C3);

    // i=7 and i=14
    C2 = -kx +ky +kz;
    C3 = -kxky_ - kxkz_ + kykz_;
    f[7]  *= one_m_omega; f[7]  += t4_omega * (C1+C2+C3);
    f[14] *= one_m_omega; f[14] += t4_omega * (C1-C2+C3);

    return invRho*invRho*jSqr;
}

static T precond_bgk_ma2_collision(Array<T,Descriptor::q>& f, T rhoBar, Array<T,3> const& j, T omega, T invGamma)
{
    return precond_bgk_ma2_collision_base(f, rhoBar, j, omega, invGamma, false);
}


};  //struct dynamicsTemplatesImpl<D3Q15DescriptorBase>

}  // namespace plb

#endif  // DYNAMICS_TEMPLATES_3D_H
