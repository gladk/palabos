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

/* Orestis Malaspinas contributed this code.
 */

/** \file
 * Template specializations for some computationally intensive LB 
 * functions of the header file mrtTemplates.h, for the D2Q9 grid.
 */

#ifndef MRT_TEMPLATES_2D_H
#define MRT_TEMPLATES_2D_H

#include "core/globalDefs.h"

namespace plb {

// Efficient specialization for D2Q9 lattice
template<typename T> struct mrtTemplatesImpl<T, descriptors::MRTD2Q9DescriptorBase<T> > {

    typedef descriptors::D2Q9DescriptorBase<T> Descriptor; 

    /// Computation of all equilibrium distribution (in moments space)
    static void computeEquilibrium( Array<T,Descriptor::q>& momentsEq,
                                    T rhoBar, Array<T,2> const& j, T jSqr )
    {
        T invRho = Descriptor::invRho(rhoBar);
        momentsEq[0] = rhoBar;
        momentsEq[1] = (T)3*jSqr*invRho-2*rhoBar;
        momentsEq[2] = -(T)3*jSqr*invRho+rhoBar;
        momentsEq[3] = j[0];
        momentsEq[4] = -j[0];
        momentsEq[5] = j[1];
        momentsEq[6] = -j[1];
        momentsEq[7] = (j[0]*j[0]-j[1]*j[1])*invRho;
        momentsEq[8] = j[1]*j[0]*invRho;
    }
    
    /// Computation of all moments (specialized for d2q9)
    static void computeMoments(Array<T,Descriptor::q>& moments, const Array<T,Descriptor::q>& f)
    {
        T f1_f3 = f[1]-f[3];
        T f1pf3 = f[1]+f[3];
        
        T f2pf4 = f[2]+f[4];
        T f5_f7 = f[5]-f[7];
        T f5pf7 = f[5]+f[7];
        T f6pf8 = f[6]+f[8];
        
        moments[0] = f[0] + f1pf3 + f2pf4 + f5pf7 + f6pf8;
        moments[1] = -(T)4*f[0] + (T)2*(f1pf3 + f5pf7) - f2pf4 - f6pf8;
        moments[2] =  (T)4*f[0] + f1pf3 - (T)2*(f2pf4 + f6pf8) + f5pf7;
        moments[3] = -f1pf3 - f[2] + f5pf7 + f[6];
        moments[4] = -f1pf3 + (T)2*(f[2] - f[6]) + f5pf7;
        moments[5] =  f1_f3 - f[4] - f5_f7 + f[8];
        moments[6] =  f1_f3 + (T)2*(f[4] - f[8]) - f5_f7;
        moments[7] =  f[2] - f[4] + f[6] - f[8];
        moments[8] = -f1_f3 - f5_f7;
    }
    
    /// MRT collision step
    static T mrtCollision( Array<T,Descriptor::q>& f,
                           const T &rhoBar, const Array<T,2> & j,
                           T invM_S[9][9] )
    {
        Array<T,9> moments, momentsEq;

        computeMoments(moments,f);
        T jSqr = VectorTemplateImpl<T,2>::normSqr(j);
        
        computeEquilibrium(momentsEq,rhoBar,j,jSqr);
        
        T mom1 = moments[1] - momentsEq[1];
        T mom2 = moments[2] - momentsEq[2];
        T mom3 = moments[3] - momentsEq[3];
        T mom4 = moments[4] - momentsEq[4];
        T mom5 = moments[5] - momentsEq[5];
        T mom6 = moments[6] - momentsEq[6];
        T mom7 = moments[7] - momentsEq[7];
        T mom8 = moments[8] - momentsEq[8];
        
        f[0] -= invM_S[0][1]*mom1 + 
                   invM_S[0][2]*mom2;
        
        f[1] -= invM_S[1][1]*mom1 +
                   invM_S[1][2]*mom2 +
                   invM_S[1][3]*mom3 +
                   invM_S[1][4]*mom4 + 
                   invM_S[1][5]*mom5 + 
                   invM_S[1][6]*mom6 +
                   invM_S[1][8]*mom8;
        
        f[2] -= invM_S[2][1]*mom1 +
                   invM_S[2][2]*mom2 +
                   invM_S[2][3]*mom3 +
                   invM_S[2][4]*mom4 +
                   invM_S[2][7]*mom7;
        
        f[3] -= invM_S[3][1]*mom1 +
                   invM_S[3][2]*mom2 +
                   invM_S[3][3]*mom3 +
                   invM_S[3][4]*mom4 +
                   invM_S[3][5]*mom5 +
                   invM_S[3][6]*mom6 +
                   invM_S[3][8]*mom8;
        
        f[4] -= invM_S[4][1]*mom1 +
                   invM_S[4][2]*mom2 +
                   invM_S[4][5]*mom5 +
                   invM_S[4][6]*mom6 +
                   invM_S[4][7]*mom7;
        
        f[5] -= invM_S[5][1]*mom1 +
                   invM_S[5][2]*mom2 +
                   invM_S[5][3]*mom3 +
                   invM_S[5][4]*mom4 +
                   invM_S[5][5]*mom5 +
                   invM_S[5][6]*mom6 +
                   invM_S[5][8]*mom8;
        
        f[6] -= invM_S[6][1]*mom1 +
                   invM_S[6][2]*mom2 +
                   invM_S[6][3]*mom3 +
                   invM_S[6][4]*mom4 +
                   invM_S[6][7]*mom7;
        
        f[7] -= invM_S[7][1]*mom1 +
                   invM_S[7][2]*mom2 +
                   invM_S[7][3]*mom3 +
                   invM_S[7][4]*mom4 +
                   invM_S[7][5]*mom5 +
                   invM_S[7][6]*mom6 +
                   invM_S[7][8]*mom8;
        
        f[8] -= invM_S[8][1]*mom1 +
                   invM_S[8][2]*mom2 +
                   invM_S[8][5]*mom5 +
                   invM_S[8][6]*mom6 +
                   invM_S[8][7]*mom7;

        return jSqr;
    }
    
    static T variableOmegaMrtCollision(
            Array<T,Descriptor::q>& f, const T &rhoBar,
            const Array<T,2> & j, T invM_S[9][9], T omega )
    {
        Array<T,9> moments, momentsEq;

        computeMoments(moments,f);
        T jSqr = VectorTemplateImpl<T,2>::normSqr(j);
        
        computeEquilibrium(momentsEq,rhoBar,j,jSqr);
        
        T mom1 = moments[1] - momentsEq[1];
        T mom2 = moments[2] - momentsEq[2];
        T mom3 = moments[3] - momentsEq[3];
        T mom4 = moments[4] - momentsEq[4];
        T mom5 = moments[5] - momentsEq[5];
        T mom6 = moments[6] - momentsEq[6];
        T mom7 = moments[7] - momentsEq[7];
        T mom8 = moments[8] - momentsEq[8];

        // Multiply shear-viscosity indices by omega.
        mom7 *= omega;
        mom8 *= omega;
        
        f[0] -= invM_S[0][1]*mom1 + 
                   invM_S[0][2]*mom2;
        
        f[1] -= invM_S[1][1]*mom1 +
                invM_S[1][2]*mom2 +
                invM_S[1][3]*mom3 +
                invM_S[1][4]*mom4 + 
                invM_S[1][5]*mom5 + 
                invM_S[1][6]*mom6 +
                invM_S[1][8]*mom8;
        
        f[2] -= invM_S[2][1]*mom1 +
                invM_S[2][2]*mom2 +
                invM_S[2][3]*mom3 +
                invM_S[2][4]*mom4 +
                invM_S[2][7]*mom7;
        
        f[3] -= invM_S[3][1]*mom1 +
                invM_S[3][2]*mom2 +
                invM_S[3][3]*mom3 +
                invM_S[3][4]*mom4 +
                invM_S[3][5]*mom5 +
                invM_S[3][6]*mom6 +
                invM_S[3][8]*mom8;
        
        f[4] -= invM_S[4][1]*mom1 +
                invM_S[4][2]*mom2 +
                invM_S[4][5]*mom5 +
                invM_S[4][6]*mom6 +
                invM_S[4][7]*mom7;
        
        f[5] -= invM_S[5][1]*mom1 +
                invM_S[5][2]*mom2 +
                invM_S[5][3]*mom3 +
                invM_S[5][4]*mom4 +
                invM_S[5][5]*mom5 +
                invM_S[5][6]*mom6 +
                invM_S[5][8]*mom8;
        
        f[6] -= invM_S[6][1]*mom1 +
                invM_S[6][2]*mom2 +
                invM_S[6][3]*mom3 +
                invM_S[6][4]*mom4 +
                invM_S[6][7]*mom7;
        
        f[7] -= invM_S[7][1]*mom1 +
                invM_S[7][2]*mom2 +
                invM_S[7][3]*mom3 +
                invM_S[7][4]*mom4 +
                invM_S[7][5]*mom5 +
                invM_S[7][6]*mom6 +
                invM_S[7][8]*mom8;
        
        f[8] -= invM_S[8][1]*mom1 +
                invM_S[8][2]*mom2 +
                invM_S[8][5]*mom5 +
                invM_S[8][6]*mom6 +
                invM_S[8][7]*mom7;

        return jSqr;
    }
    
    static void addGuoForce( Array<T,Descriptor::q>& f, const Array<T,Descriptor::d>& force,
                             Array<T,Descriptor::d> const& u,
                             T invM_S[Descriptor::q][Descriptor::q], T amplitude )
    {
        Array<T,Descriptor::q> forcing;
        T g_u   = force[0] * u[0] + force[1] * u[1];
        T g_u_a = force[0] * u[0] - force[1] * u[1];            
        
        T mom1 = -(T)3 * g_u;
        T mom2 =  (T)3 * g_u;
        T mom3 = -(T)0.5 * force[0];
        T mom4 = -mom3;
        T mom5 = -(T)0.5 * force[1];
        T mom6 = -mom5;
        T mom7 = -g_u_a;
        T mom8 = -(T)0.5 * (force[0]*u[1] + force[1]*u[0]);

        f[0] += invM_S[0][1]*mom1 + 
                   invM_S[0][2]*mom2;
        
        f[1] += invM_S[1][1]*mom1 +
                   invM_S[1][2]*mom2 +
                   invM_S[1][3]*mom3 +
                   invM_S[1][4]*mom4 + 
                   invM_S[1][5]*mom5 + 
                   invM_S[1][6]*mom6 +
                   invM_S[1][8]*mom8;
        
        f[2] += invM_S[2][1]*mom1 +
                   invM_S[2][2]*mom2 +
                   invM_S[2][3]*mom3 +
                   invM_S[2][4]*mom4 +
                   invM_S[2][7]*mom7;
        
        f[3] += invM_S[3][1]*mom1 +
                   invM_S[3][2]*mom2 +
                   invM_S[3][3]*mom3 +
                   invM_S[3][4]*mom4 +
                   invM_S[3][5]*mom5 +
                   invM_S[3][6]*mom6 +
                   invM_S[3][8]*mom8;
        
        f[4] += invM_S[4][1]*mom1 +
                   invM_S[4][2]*mom2 +
                   invM_S[4][5]*mom5 +
                   invM_S[4][6]*mom6 +
                   invM_S[4][7]*mom7;
        
        f[5] += invM_S[5][1]*mom1 +
                   invM_S[5][2]*mom2 +
                   invM_S[5][3]*mom3 +
                   invM_S[5][4]*mom4 +
                   invM_S[5][5]*mom5 +
                   invM_S[5][6]*mom6 +
                   invM_S[5][8]*mom8;
        
        f[6] += invM_S[6][1]*mom1 +
                   invM_S[6][2]*mom2 +
                   invM_S[6][3]*mom3 +
                   invM_S[6][4]*mom4 +
                   invM_S[6][7]*mom7;
        
        f[7] += invM_S[7][1]*mom1 +
                   invM_S[7][2]*mom2 +
                   invM_S[7][3]*mom3 +
                   invM_S[7][4]*mom4 +
                   invM_S[7][5]*mom5 +
                   invM_S[7][6]*mom6 +
                   invM_S[7][8]*mom8;
        
        f[8] += invM_S[8][1]*mom1 +
                   invM_S[8][2]*mom2 +
                   invM_S[8][5]*mom5 +
                   invM_S[8][6]*mom6 +
                   invM_S[8][7]*mom7;
                   
        f[0] += -(T)2/(T)3*g_u;
        f[1] += (-force[0]+force[1]+(T)2*g_u-(T)3*g_u)/(T)12;
        f[2] += (-force[0]+2*force[0]*u[0]-force[1]*u[1])/(T)3;
        f[3] += (-force[0]-force[1]+(T)2*g_u+(T)3*g_u) / (T)12;
        f[4] += -(force[1]+force[0]*u[0] - (T)2*force[1]*u[1]) / (T)3;
        f[5] += (force[0]-force[1]+(T)2*g_u-(T)3*g_u) / (T)12;
        f[6] += (force[0]+(T)2*force[0]*u[0]-force[1]*u[1]) / (T)3;
        f[7] += (force[0]+force[1]+(T)2*g_u+(T)3*g_u) / (T)12;
        f[8] += -(-force[1]+force[0]*u[0]-2*force[1]*u[1])/(T)3;
    }
    
    static void variableOmegaAddGuoForce(
            Array<T,Descriptor::q>& f, const Array<T,Descriptor::d>& force,
            Array<T,Descriptor::d> const& u,
            T invM_S[Descriptor::q][Descriptor::q], T amplitude, T omega )
    {
        Array<T,Descriptor::q> forcing;
        T g_u   = force[0] * u[0] + force[1] * u[1];
        T g_u_a = force[0] * u[0] - force[1] * u[1];            
        
        T mom1 = -(T)3 * g_u;
        T mom2 =  (T)3 * g_u;
        T mom3 = -(T)0.5 * force[0];
        T mom4 = -mom3;
        T mom5 = -(T)0.5 * force[1];
        T mom6 = -mom5;
        T mom7 = -g_u_a;
        T mom8 = -(T)0.5 * (force[0]*u[1] + force[1]*u[0]);

        // Multiply shear-viscosity indices by omega.
        mom7 *= omega;
        mom8 *= omega;
        
        f[0] += invM_S[0][1]*mom1 + 
                   invM_S[0][2]*mom2;
        
        f[1] += invM_S[1][1]*mom1 +
                   invM_S[1][2]*mom2 +
                   invM_S[1][3]*mom3 +
                   invM_S[1][4]*mom4 + 
                   invM_S[1][5]*mom5 + 
                   invM_S[1][6]*mom6 +
                   invM_S[1][8]*mom8;
        
        f[2] += invM_S[2][1]*mom1 +
                   invM_S[2][2]*mom2 +
                   invM_S[2][3]*mom3 +
                   invM_S[2][4]*mom4 +
                   invM_S[2][7]*mom7;
        
        f[3] += invM_S[3][1]*mom1 +
                   invM_S[3][2]*mom2 +
                   invM_S[3][3]*mom3 +
                   invM_S[3][4]*mom4 +
                   invM_S[3][5]*mom5 +
                   invM_S[3][6]*mom6 +
                   invM_S[3][8]*mom8;
        
        f[4] += invM_S[4][1]*mom1 +
                   invM_S[4][2]*mom2 +
                   invM_S[4][5]*mom5 +
                   invM_S[4][6]*mom6 +
                   invM_S[4][7]*mom7;
        
        f[5] += invM_S[5][1]*mom1 +
                   invM_S[5][2]*mom2 +
                   invM_S[5][3]*mom3 +
                   invM_S[5][4]*mom4 +
                   invM_S[5][5]*mom5 +
                   invM_S[5][6]*mom6 +
                   invM_S[5][8]*mom8;
        
        f[6] += invM_S[6][1]*mom1 +
                   invM_S[6][2]*mom2 +
                   invM_S[6][3]*mom3 +
                   invM_S[6][4]*mom4 +
                   invM_S[6][7]*mom7;
        
        f[7] += invM_S[7][1]*mom1 +
                   invM_S[7][2]*mom2 +
                   invM_S[7][3]*mom3 +
                   invM_S[7][4]*mom4 +
                   invM_S[7][5]*mom5 +
                   invM_S[7][6]*mom6 +
                   invM_S[7][8]*mom8;
        
        f[8] += invM_S[8][1]*mom1 +
                   invM_S[8][2]*mom2 +
                   invM_S[8][5]*mom5 +
                   invM_S[8][6]*mom6 +
                   invM_S[8][7]*mom7;
                   
        f[0] += -(T)2/(T)3*g_u;
        f[1] += (-force[0]+force[1]+(T)2*g_u-(T)3*g_u)/(T)12;
        f[2] += (-force[0]+2*force[0]*u[0]-force[1]*u[1])/(T)3;
        f[3] += (-force[0]-force[1]+(T)2*g_u+(T)3*g_u) / (T)12;
        f[4] += -(force[1]+force[0]*u[0] - (T)2*force[1]*u[1]) / (T)3;
        f[5] += (force[0]-force[1]+(T)2*g_u-(T)3*g_u) / (T)12;
        f[6] += (force[0]+(T)2*force[0]*u[0]-force[1]*u[1]) / (T)3;
        f[7] += (force[0]+force[1]+(T)2*g_u+(T)3*g_u) / (T)12;
        f[8] += -(-force[1]+force[0]*u[0]-2*force[1]*u[1])/(T)3;
    }
    
    /// MRT collision step
    static T mrtCollisionWithForce( Array<T,Descriptor::q>& f,
                                    const T &rhoBar, const Array<T,Descriptor::d> & u,
                                    T invM_S[Descriptor::q][Descriptor::q], 
                                    const Array<T,Descriptor::d> &force, T amplitude) 
    {
        Array<T,Descriptor::d> j = Descriptor::fullRho(rhoBar)*u;
        T jSqr = mrtCollision( f, rhoBar, j, invM_S );
        addGuoForce( f, force, u, invM_S, amplitude );
        
        return jSqr;
    }
    
    /// MRT collision step
    static T quasiIncMrtCollisionWithForce( Array<T,Descriptor::q>& f,
                                    const T &rhoBar, const Array<T,Descriptor::d> & u,
                                    T invM_S[Descriptor::q][Descriptor::q], 
                                    const Array<T,Descriptor::d> &force, T amplitude) 
    {
        Array<T,Descriptor::d> j = Descriptor::fullRho(rhoBar)*u;
        T jSqr = quasiIncMrtCollision( f, rhoBar, j, invM_S );
        addGuoForce( f, force, u, invM_S, amplitude );
        
        return jSqr;
    }
    
    static T variableOmegaMrtCollisionWithForce (
            Array<T,Descriptor::q>& f, const T &rhoBar,
            const Array<T,Descriptor::d>&  u, T invM_S[Descriptor::q][Descriptor::q], 
            const Array<T,Descriptor::d>& force, T amplitude, T omega )
    {
        Array<T,Descriptor::d> j = Descriptor::fullRho(rhoBar)*u;
        T jSqr = variableOmegaMrtCollision( f, rhoBar, j, invM_S, omega );
        variableOmegaAddGuoForce( f, force, u, invM_S, amplitude, omega );
        
        return jSqr;
    }
    
    /// Computation of all equilibrium distribution (in moments space)
    static void computeQuasiIncEquilibrium( Array<T,Descriptor::q>& momentsEq,
                                    T rhoBar, Array<T,2> const& j, T jSqr )
    {
        momentsEq[0] = rhoBar;
        momentsEq[1] = (T)3*jSqr-2*rhoBar;
        momentsEq[2] = -(T)3*jSqr+rhoBar;
        momentsEq[3] = j[0];
        momentsEq[4] = -j[0];
        momentsEq[5] = j[1];
        momentsEq[6] = -j[1];
        momentsEq[7] = (j[0]*j[0]-j[1]*j[1]);
        momentsEq[8] = j[1]*j[0];
    }
    
    /// MRT collision step
    static T quasiIncMrtCollision( Array<T,Descriptor::q>& f,
                                   const T &rhoBar, const Array<T,2> & j,
                                   T invM_S[9][9] )
    {
        Array<T,9> moments, momentsEq;

        computeMoments(moments,f);
        T jSqr = VectorTemplateImpl<T,2>::normSqr(j);
        
        computeQuasiIncEquilibrium(momentsEq,rhoBar,j,jSqr);
        
        T mom1 = moments[1] - momentsEq[1];
        T mom2 = moments[2] - momentsEq[2];
        T mom3 = moments[3] - momentsEq[3];
        T mom4 = moments[4] - momentsEq[4];
        T mom5 = moments[5] - momentsEq[5];
        T mom6 = moments[6] - momentsEq[6];
        T mom7 = moments[7] - momentsEq[7];
        T mom8 = moments[8] - momentsEq[8];
        
        f[0] -= invM_S[0][1]*mom1 + 
                   invM_S[0][2]*mom2;
        
        f[1] -= invM_S[1][1]*mom1 +
                   invM_S[1][2]*mom2 +
                   invM_S[1][4]*mom4 + 
                   invM_S[1][6]*mom6 +
                   invM_S[1][8]*mom8;
        
        f[2] -= invM_S[2][1]*mom1 +
                   invM_S[2][2]*mom2 +
                   invM_S[2][4]*mom4 +
                   invM_S[2][7]*mom7;
        
        f[3] -= invM_S[3][1]*mom1 +
                   invM_S[3][2]*mom2 +
                   invM_S[3][4]*mom4 +
                   invM_S[3][6]*mom6 +
                   invM_S[3][8]*mom8;
        
        f[4] -= invM_S[4][1]*mom1 +
                   invM_S[4][2]*mom2 +
                   invM_S[4][6]*mom6 +
                   invM_S[4][7]*mom7;
        
        f[5] -= invM_S[5][1]*mom1 +
                   invM_S[5][2]*mom2 +
                   invM_S[5][4]*mom4 +
                   invM_S[5][6]*mom6 +
                   invM_S[5][8]*mom8;
        
        f[6] -= invM_S[6][1]*mom1 +
                   invM_S[6][2]*mom2 +
                   invM_S[6][4]*mom4 +
                   invM_S[6][7]*mom7;
        
        f[7] -= invM_S[7][1]*mom1 +
                   invM_S[7][2]*mom2 +
                   invM_S[7][4]*mom4 +
                   invM_S[7][6]*mom6 +
                   invM_S[7][8]*mom8;
        
        f[8] -= invM_S[8][1]*mom1 +
                   invM_S[8][2]*mom2 +
                   invM_S[8][6]*mom6 +
                   invM_S[8][7]*mom7;

        return jSqr;
    }

};


}  // namespace plb

#endif
