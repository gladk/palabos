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

/* Main author: Orestis Malaspinas
 */

/** \file
 * Helper functions for the implementation of LB dynamics. This file is all
 * about efficiency. The generic template code is specialized for commonly
 * used Lattices, so that a maximum performance can be taken out of each
 * case.
 */
#ifndef ADVECTION_DIFFUSION_DYNAMICS_TEMPLATES_3D_H
#define ADVECTION_DIFFUSION_DYNAMICS_TEMPLATES_3D_H

namespace plb {

/// All helper functions are inside this structure
template<typename T> struct advectionDiffusionDynamicsTemplatesImpl<T,descriptors::D3Q7DescriptorBase<T> >
{
    
typedef descriptors::D3Q7DescriptorBase<T> Descriptor;

static T bgk_ma1_equilibrium(plint iPop, T rhoBar, Array<T,Descriptor::d> const& jEq) 
{
    return Descriptor::t[iPop] * (rhoBar + Descriptor::invCs2 * 
            (Descriptor::c[iPop][0]*jEq[0]+Descriptor::c[iPop][1]*jEq[1]+Descriptor::c[iPop][2]*jEq[2]));
}

/// Regularization
static void regularize( Array<T,Descriptor::q>& f, T rhoBar,
                        Array<T,Descriptor::d> const& jAdvDiff,
                        Array<T,Descriptor::d> const& jEq )
{
    f[0] = Descriptor::t[0] * rhoBar;
    
    f[1] = Descriptor::t[1] * (rhoBar - Descriptor::invCs2*jAdvDiff[0]);
    f[2] = Descriptor::t[2] * (rhoBar - Descriptor::invCs2*jAdvDiff[1]);
    f[3] = Descriptor::t[3] * (rhoBar - Descriptor::invCs2*jAdvDiff[2]);
    f[4] = Descriptor::t[4] * (rhoBar + Descriptor::invCs2*jAdvDiff[0]);
    f[5] = Descriptor::t[5] * (rhoBar + Descriptor::invCs2*jAdvDiff[1]);
    f[6] = Descriptor::t[6] * (rhoBar + Descriptor::invCs2*jAdvDiff[2]);
}

static T no_corr_bgk_collision(
        Array<T,Descriptor::q>& f, T rhoBar, Array<T,Descriptor::d> const& jEq, 
        T omega) 
{
    T invRho = Descriptor::invRho(rhoBar);
    const T jSqr = jEq[0]*jEq[0] + jEq[1]*jEq[1] + jEq[2]*jEq[2];
    
    const T oneMinusOmega = (T)1 - omega;
    const T halfOmega = (T)0.5 * omega;
    const T cs2RhoBar = Descriptor::cs2*rhoBar;
    
    f[0] = oneMinusOmega*f[0]+omega*((T)1-(T)3*Descriptor::cs2)*rhoBar;
    
    f[1] = oneMinusOmega*f[1]+halfOmega*(cs2RhoBar-jEq[0]);
    f[2] = oneMinusOmega*f[2]+halfOmega*(cs2RhoBar-jEq[1]);
    f[3] = oneMinusOmega*f[3]+halfOmega*(cs2RhoBar-jEq[2]);
    
    f[4] = oneMinusOmega*f[4]+halfOmega*(cs2RhoBar+jEq[0]);
    f[5] = oneMinusOmega*f[5]+halfOmega*(cs2RhoBar+jEq[1]);
    f[6] = oneMinusOmega*f[6]+halfOmega*(cs2RhoBar+jEq[2]);
    
    return jSqr*invRho*invRho;
}

static T no_corr_bgk_collision(
        Array<T,Descriptor::q>& f, T rhoBar, Array<T,Descriptor::d> const& jEq, 
        T omega, T source) 
{
    T invRho = Descriptor::invRho(rhoBar);
    const T jSqr = jEq[0]*jEq[0] + jEq[1]*jEq[1] + jEq[2]*jEq[2];
    
    const T oneMinusOmega = (T)1 - omega;
    const T halfOmega = (T)0.5 * omega;
    const T cs2RhoBar = Descriptor::cs2*rhoBar;
    const T halfSourceCs2 = (T)0.5*source * Descriptor::cs2;
    
    f[0] = oneMinusOmega*f[0]+((T)1-(T)3*Descriptor::cs2)*(omega*rhoBar+source);
    
    f[1] = oneMinusOmega*f[1]+halfOmega*(cs2RhoBar-jEq[0]) + halfSourceCs2;
    f[2] = oneMinusOmega*f[2]+halfOmega*(cs2RhoBar-jEq[1]) + halfSourceCs2;
    f[3] = oneMinusOmega*f[3]+halfOmega*(cs2RhoBar-jEq[2]) + halfSourceCs2;
    
    f[4] = oneMinusOmega*f[4]+halfOmega*(cs2RhoBar+jEq[0]) + halfSourceCs2;
    f[5] = oneMinusOmega*f[5]+halfOmega*(cs2RhoBar+jEq[1]) + halfSourceCs2;
    f[6] = oneMinusOmega*f[6]+halfOmega*(cs2RhoBar+jEq[2]) + halfSourceCs2;
    
    return jSqr*invRho*invRho;
}

static T no_corr_rlb_collision (
    Array<T,Descriptor::q>& f, T rhoBar, Array<T,Descriptor::d> const& jEq,
    Array<T,Descriptor::d> const& jNeq,T omega )
{
    T invRho = Descriptor::invRho(rhoBar);
    const T jSqr = jEq[0]*jEq[0] + jEq[1]*jEq[1] + jEq[2]*jEq[2];
    
    const T oneHalfMinusHalfOmega = (T)0.5-(T)0.5*omega;
    const T cs2RhoBar = Descriptor::cs2 * rhoBar;
    
    const T jNeqTerm_0 = oneHalfMinusHalfOmega*jNeq[0];
    const T jNeqTerm_1 = oneHalfMinusHalfOmega*jNeq[1];
    const T jNeqTerm_2 = oneHalfMinusHalfOmega*jNeq[2];
    
    f[0] = ((T)1-(T)3*Descriptor::cs2)*rhoBar;
    
    f[1] = -jNeqTerm_0 + (T)0.5*(cs2RhoBar-jEq[0]);
    f[2] = -jNeqTerm_1 + (T)0.5*(cs2RhoBar-jEq[1]);
    f[3] = -jNeqTerm_2 + (T)0.5*(cs2RhoBar-jEq[2]);
    
    f[4] = -f[1] + cs2RhoBar;
    f[5] = -f[2] + cs2RhoBar;
    f[6] = -f[3] + cs2RhoBar;
            
    return jSqr*invRho*invRho;
}

static T no_corr_rlb_collision (
    Array<T,Descriptor::q>& f, T rhoBar, Array<T,Descriptor::d> const& jEq,
    Array<T,Descriptor::d> const& jNeq, T omega, T source )
{
    T invRho = Descriptor::invRho(rhoBar);
    const T jSqr = jEq[0]*jEq[0] + jEq[1]*jEq[1] + jEq[2]*jEq[2];
    
    const T oneHalfMinusHalfOmega = (T)0.5-(T)0.5*omega;
    const T cs2RhoBar = Descriptor::cs2 * rhoBar;
    const T halfSourceCs2 = (T)0.5*source * Descriptor::cs2;
    
    const T jNeqTerm_0 = oneHalfMinusHalfOmega*jNeq[0];
    const T jNeqTerm_1 = oneHalfMinusHalfOmega*jNeq[1];
    const T jNeqTerm_2 = oneHalfMinusHalfOmega*jNeq[2];
    
    f[0] = ((T)1-(T)3*Descriptor::cs2)*(rhoBar+source);
    
    f[1] = -jNeqTerm_0 + (T)0.5*(cs2RhoBar-jEq[0]);
    f[2] = -jNeqTerm_1 + (T)0.5*(cs2RhoBar-jEq[1]);
    f[3] = -jNeqTerm_2 + (T)0.5*(cs2RhoBar-jEq[2]);
    
    f[4] = -f[1] + cs2RhoBar;
    f[5] = -f[2] + cs2RhoBar;
    f[6] = -f[3] + cs2RhoBar;

    f[1] += halfSourceCs2;
    f[2] += halfSourceCs2;
    f[3] += halfSourceCs2;
    f[4] += halfSourceCs2;
    f[5] += halfSourceCs2;
    f[6] += halfSourceCs2;
            
    return jSqr*invRho*invRho;
}

};  // struct advectionDiffusionDynamicsTemplatesImpl

}  // namespace plb

#endif
