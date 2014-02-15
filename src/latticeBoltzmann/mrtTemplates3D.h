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
 * functions of the header file mrtTemplates.h, for the D3Q19 grid.
 */

#ifndef MRT_TEMPLATES_3D_H
#define MRT_TEMPLATES_3D_H

#include "core/globalDefs.h"

namespace plb {

// Efficient specialization for D3Q19 lattice
template<typename T> struct mrtTemplatesImpl<T, descriptors::MRTD3Q19DescriptorBase<T> > {
    
    typedef descriptors::D3Q19DescriptorBase<T> Descriptor; 

    /// Computation of all equilibrium distribution (in moments space)
    static void computeEquilibrium( Array<T,Descriptor::q>& momentsEq,
                                    T rhoBar, Array<T,3> const& j, T jSqr )
    {
        T invRho = Descriptor::invRho(rhoBar);
        
        momentsEq[0] = rhoBar;
        momentsEq[1] = (T)19*jSqr*invRho-(T)11*rhoBar;
        momentsEq[2] = -(T)5.5*jSqr*invRho+(T)3*rhoBar;
        momentsEq[3] = j[0];
        momentsEq[4] = -((T)2/(T)3)*j[0];
        momentsEq[5] = j[1];
        momentsEq[6] = -((T)2/3)*j[1];
        momentsEq[7] = j[2];
        momentsEq[8] = -((T)2/(T)3)*j[2];
        momentsEq[9] = ((T)2*j[0]*j[0]-j[1]*j[1]-j[2]*j[2])*invRho;
        momentsEq[10] = (-j[0]*j[0]+(T)0.5*j[1]*j[1]+(T)0.5*j[2]*j[2])*invRho;
        momentsEq[11] = (j[1]*j[1]-j[2]*j[2])*invRho;
        momentsEq[12] = (-(T)0.5*j[1]*j[1]+(T)0.5*j[2]*j[2])*invRho;
        momentsEq[13] = j[1]*j[0]*invRho;
        momentsEq[14] = j[2]*j[1]*invRho;
        momentsEq[15] = j[2]*j[0]*invRho;
        momentsEq[16] = T();
        momentsEq[17] = T();
        momentsEq[18] = T();
        
    }
    
    /// Computation of all equilibrium distribution (in moments space)
    static void computeSmagorinskyEquilibrium( Array<T,Descriptor::q>& momentsEq,
                                                       T rhoBar, Array<T,3> const& j, T jSqr, const Array<T,6> &strain, T cSmago )
    {
        typedef SymmetricTensorImpl<T,3> S;
        
        T invRho = Descriptor::invRho(rhoBar);
        T rho = Descriptor::fullRho(rhoBar);
        T rho2 = rho*rho;
        
        T sNorm = sqrt((T)2*SymmetricTensorImpl<T,3>::tensorNormSqr(strain));
        T smagoFactor = (T)2*cSmago * cSmago * sNorm;
        
        T ux2 = rho2 * smagoFactor * strain[S::xx];
        T uy2 = rho2 * smagoFactor * strain[S::yy];
        T uz2 = rho2 * smagoFactor * strain[S::zz];
        
        T uxuy = rho2 * smagoFactor * strain[S::xy];
        T uyuz = rho2 * smagoFactor * strain[S::yz];
        T uxuz = rho2 * smagoFactor * strain[S::xz];
        
        momentsEq[0] = rhoBar;
        momentsEq[1] = (T)19*(jSqr + ux2+uy2+uz2)*invRho-(T)11*rhoBar;
        momentsEq[2] = -(T)5.5*(jSqr + ux2+uy2+uz2)*invRho + (T)3*rhoBar;
        momentsEq[3] = j[0];
        momentsEq[4] = -((T)2/(T)3)*j[0];
        momentsEq[5] = j[1];
        momentsEq[6] = -((T)2/3)*j[1];
        momentsEq[7] = j[2];
        momentsEq[8] = -((T)2/(T)3)*j[2];
        momentsEq[9] = ((T)2*(j[0]*j[0]+ux2)-(j[1]*j[1]+uy2)-(j[2]*j[2]+uz2))*invRho;
        momentsEq[10] = (-(j[0]*j[0]+ux2)+(T)0.5*(j[1]*j[1]+uy2)+(T)0.5*(j[2]*j[2]+uz2))*invRho;
        momentsEq[11] = (j[1]*j[1]+uy2-(j[2]*j[2]+uz2))*invRho;
        momentsEq[12] = (-(T)0.5*(j[1]*j[1]+uy2)+(T)0.5*(j[2]*j[2]+uz2))*invRho;
        momentsEq[13] = (j[1]*j[0]+uxuy)*invRho;
        momentsEq[14] = (j[2]*j[1]+uyuz)*invRho;
        momentsEq[15] = (j[2]*j[0]+uxuz)*invRho;
        momentsEq[16] = T();
        momentsEq[17] = T();
        momentsEq[18] = T();
    }
    
    /// Computation of all moments (specialized for d3q19)
    static void computeMoments(Array<T,Descriptor::q>& moments, const Array<T,Descriptor::q>& f)
    {
        moments[0] = f[0]+f[1]+f[2]+f[3]+f[4]+f[5]+f[6]+f[7]+f[8]+f[9]+f[10]+f[11]+f[12]+f[13]+f[14]+f[15]+f[16]+f[17]+f[18];
        moments[1] =
            -(T)30*f[0]-(T)11*(f[1]+f[2]+f[3]+f[10]+f[11]+f[12])
            +(T)8*(f[4]+f[5]+f[6]+f[7]+f[8]+f[9]+f[13]+f[14]+f[15]+f[16]+f[17]+f[18]);
        moments[2] = 12*f[0]-4*f[1]-4*f[2]-4*f[3]+f[4]+f[5]+f[6]+f[7]+f[8]+f[9]-4*f[10]-4*f[11]-4*f[12]+f[13]+f[14]+f[15]+f[16]+f[17]+f[18];
        moments[3] = -f[1]-f[4]-f[5]-f[6]-f[7]+f[10]+f[13]+f[14]+f[15]+f[16];
        moments[4] = 4*f[1]-f[4]-f[5]-f[6]-f[7]-4*f[10]+f[13]+f[14]+f[15]+f[16];
        moments[5] = -f[2]-f[4]+f[5]-f[8]-f[9]+f[11]+f[13]-f[14]+f[17]+f[18];
        moments[6] = 4*f[2]-f[4]+f[5]-f[8]-f[9]-4*f[11]+f[13]-f[14]+f[17]+f[18];
        moments[7] = -f[3]-f[6]+f[7]-f[8]+f[9]+f[12]+f[15]-f[16]+f[17]-f[18];
        moments[8] = 4*f[3]-f[6]+f[7]-f[8]+f[9]-4*f[12]+f[15]-f[16]+f[17]-f[18];
        moments[9] = 2*f[1]-f[2]-f[3]+f[4]+f[5]+f[6]+f[7]-2*f[8]-2*f[9]+2*f[10]-f[11]-f[12]+f[13]+f[14]+f[15]+f[16]-2*f[17]-2*f[18];
        moments[10] = -4*f[1]+2*f[2]+2*f[3]+f[4]+f[5]+f[6]+f[7]-2*f[8]-2*f[9]-4*f[10]+2*f[11]+2*f[12]+f[13]+f[14]+f[15]+f[16]-2*f[17]-2*f[18];
        moments[11] = f[2]-f[3]+f[4]+f[5]-f[6]-f[7]+f[11]-f[12]+f[13]+f[14]-f[15]-f[16];
        moments[12] = -2*f[2]+2*f[3]+f[4]+f[5]-f[6]-f[7]-2*f[11]+2*f[12]+f[13]+f[14]-f[15]-f[16];
        moments[13] = f[4]-f[5]+f[13]-f[14];
        moments[14] = f[8]-f[9]+f[17]-f[18];
        moments[15] = f[6]-f[7]+f[15]-f[16];
        moments[16] = -f[4]-f[5]+f[6]+f[7]+f[13]+f[14]-f[15]-f[16];
        moments[17] = f[4]-f[5]-f[8]-f[9]-f[13]+f[14]+f[17]+f[18];
        moments[18] = -f[6]+f[7]+f[8]-f[9]+f[15]-f[16]-f[17]+f[18];
        
    }
    
    /// MRT collision step
    static T mrtCollision( Array<T,Descriptor::q>& f,
                           const T &rhoBar, const Array<T,3> & j,
                           const T invM_S[19][19] )
    {
        
        Array<T,19> moments, momentsEq;
        
        computeMoments(moments,f);
        
        T jSqr = VectorTemplateImpl<T,3>::normSqr(j);
        
        computeEquilibrium(momentsEq,rhoBar,j,jSqr);
        
        T mom1 = moments[1] - momentsEq[1];
        T mom2 = moments[2] - momentsEq[2];
        T mom3 = moments[3] - momentsEq[3];
        T mom4 = moments[4] - momentsEq[4];
        T mom5 = moments[5] - momentsEq[5];
        T mom6 = moments[6] - momentsEq[6];
        T mom7 = moments[7] - momentsEq[7];
        T mom8 = moments[8] - momentsEq[8];
        T mom9 = moments[9] - momentsEq[9];
        T mom10 = moments[10] - momentsEq[10];
        T mom11 = moments[11] - momentsEq[11];
        T mom12 = moments[12] - momentsEq[12];
        T mom13 = moments[13] - momentsEq[13];
        T mom14 = moments[14] - momentsEq[14];
        T mom15 = moments[15] - momentsEq[15];
        T mom16 = moments[16];
        T mom17 = moments[17];
        T mom18 = moments[18];
        
        
        f[0] -= invM_S[0][1]*mom1
        +invM_S[0][2]*mom2;
        
        f[1] -= invM_S[1][1]*mom1
        +invM_S[1][2]*mom2
        +invM_S[1][3]*mom3
        +invM_S[1][4]*mom4
        +invM_S[1][9]*mom9
        +invM_S[1][10]*mom10;
        
        f[2] -= invM_S[2][1]*mom1
        +invM_S[2][2]*mom2
        +invM_S[2][5]*mom5
        +invM_S[2][6]*mom6
        +invM_S[2][9]*mom9
        +invM_S[2][10]*mom10
        +invM_S[2][11]*mom11
        +invM_S[2][12]*mom12;
        
        f[3] -= invM_S[3][1]*mom1
        +invM_S[3][2]*mom2
        +invM_S[3][7]*mom7
        +invM_S[3][8]*mom8
        +invM_S[3][9]*mom9
        +invM_S[3][10]*mom10
        +invM_S[3][11]*mom11
        +invM_S[3][12]*mom12;
        
        f[4] -= invM_S[4][1]*mom1
        +invM_S[4][2]*mom2
        +invM_S[4][3]*mom3
        +invM_S[4][4]*mom4
        +invM_S[4][5]*mom5
        +invM_S[4][6]*mom6
        +invM_S[4][9]*mom9
        +invM_S[4][10]*mom10
        +invM_S[4][11]*mom11
        +invM_S[4][12]*mom12
        +invM_S[4][13]*mom13
        +invM_S[4][16]*mom16
        +invM_S[4][17]*mom17;
        
        f[5] -= invM_S[5][1]*mom1
        +invM_S[5][2]*mom2
        +invM_S[5][3]*mom3
        +invM_S[5][4]*mom4
        +invM_S[5][5]*mom5
        +invM_S[5][6]*mom6
        +invM_S[5][9]*mom9
        +invM_S[5][10]*mom10
        +invM_S[5][11]*mom11
        +invM_S[5][12]*mom12
        +invM_S[5][13]*mom13
        +invM_S[5][16]*mom16
        +invM_S[5][17]*mom17;
        
        f[6] -= invM_S[6][1]*mom1
        +invM_S[6][2]*mom2
        +invM_S[6][3]*mom3
        +invM_S[6][4]*mom4
        +invM_S[6][7]*mom7
        +invM_S[6][8]*mom8
        +invM_S[6][9]*mom9
        +invM_S[6][10]*mom10
        +invM_S[6][11]*mom11
        +invM_S[6][12]*mom12
        +invM_S[6][15]*mom15
        +invM_S[6][16]*mom16
        +invM_S[6][18]*mom18;
        
        f[7] -= invM_S[7][1]*mom1
        +invM_S[7][2]*mom2
        +invM_S[7][3]*mom3
        +invM_S[7][4]*mom4
        +invM_S[7][7]*mom7
        +invM_S[7][8]*mom8
        +invM_S[7][9]*mom9
        +invM_S[7][10]*mom10
        +invM_S[7][11]*mom11
        +invM_S[7][12]*mom12
        +invM_S[7][15]*mom15
        +invM_S[7][16]*mom16
        +invM_S[7][18]*mom18;
        
        f[8] -= invM_S[8][1]*mom1
        +invM_S[8][2]*mom2
        +invM_S[8][5]*mom5
        +invM_S[8][6]*mom6
        +invM_S[8][7]*mom7
        +invM_S[8][8]*mom8
        +invM_S[8][9]*mom9
        +invM_S[8][10]*mom10
        +invM_S[8][14]*mom14
        +invM_S[8][17]*mom17
        +invM_S[8][18]*mom18;
        
        f[9] -= invM_S[9][1]*mom1
        +invM_S[9][2]*mom2
        +invM_S[9][5]*mom5
        +invM_S[9][6]*mom6
        +invM_S[9][7]*mom7
        +invM_S[9][8]*mom8
        +invM_S[9][9]*mom9
        +invM_S[9][10]*mom10
        +invM_S[9][14]*mom14
        +invM_S[9][17]*mom17
        +invM_S[9][18]*mom18;
        
        f[10] -= invM_S[10][1]*mom1
        +invM_S[10][2]*mom2
        +invM_S[10][3]*mom3
        +invM_S[10][4]*mom4
        +invM_S[10][9]*mom9
        +invM_S[10][10]*mom10;
        
        f[11] -= invM_S[11][1]*mom1
        +invM_S[11][2]*mom2
        +invM_S[11][5]*mom5
        +invM_S[11][6]*mom6
        +invM_S[11][9]*mom9
        +invM_S[11][10]*mom10
        +invM_S[11][11]*mom11
        +invM_S[11][12]*mom12;
        
        f[12] -= invM_S[12][1]*mom1
        +invM_S[12][2]*mom2
        +invM_S[12][7]*mom7
        +invM_S[12][8]*mom8
        +invM_S[12][9]*mom9
        +invM_S[12][10]*mom10
        +invM_S[12][11]*mom11
        +invM_S[12][12]*mom12;
        
        f[13] -= invM_S[13][1]*mom1
        +invM_S[13][2]*mom2
        +invM_S[13][3]*mom3
        +invM_S[13][4]*mom4
        +invM_S[13][5]*mom5
        +invM_S[13][6]*mom6
        +invM_S[13][9]*mom9
        +invM_S[13][10]*mom10
        +invM_S[13][11]*mom11
        +invM_S[13][12]*mom12
        +invM_S[13][13]*mom13
        +invM_S[13][16]*mom16
        +invM_S[13][17]*mom17;
        
        f[14] -= invM_S[14][1]*mom1
        +invM_S[14][2]*mom2
        +invM_S[14][3]*mom3
        +invM_S[14][4]*mom4
        +invM_S[14][5]*mom5
        +invM_S[14][6]*mom6
        +invM_S[14][9]*mom9
        +invM_S[14][10]*mom10
        +invM_S[14][11]*mom11
        +invM_S[14][12]*mom12
        +invM_S[14][13]*mom13
        +invM_S[14][16]*mom16
        +invM_S[14][17]*mom17;
        
        f[15] -= invM_S[15][1]*mom1
        +invM_S[15][2]*mom2
        +invM_S[15][3]*mom3
        +invM_S[15][4]*mom4
        +invM_S[15][7]*mom7
        +invM_S[15][8]*mom8
        +invM_S[15][9]*mom9
        +invM_S[15][10]*mom10
        +invM_S[15][11]*mom11
        +invM_S[15][12]*mom12
        +invM_S[15][15]*mom15
        +invM_S[15][16]*mom16
        +invM_S[15][18]*mom18;
        
        f[16] -= invM_S[16][1]*mom1
        +invM_S[16][2]*mom2
        +invM_S[16][3]*mom3
        +invM_S[16][4]*mom4
        +invM_S[16][7]*mom7
        +invM_S[16][8]*mom8
        +invM_S[16][9]*mom9
        +invM_S[16][10]*mom10
        +invM_S[16][11]*mom11
        +invM_S[16][12]*mom12
        +invM_S[16][15]*mom15
        +invM_S[16][16]*mom16
        +invM_S[16][18]*mom18;
        
        f[17] -= invM_S[17][1]*mom1
        +invM_S[17][2]*mom2
        +invM_S[17][5]*mom5
        +invM_S[17][6]*mom6
        +invM_S[17][7]*mom7
        +invM_S[17][8]*mom8
        +invM_S[17][9]*mom9
        +invM_S[17][10]*mom10
        +invM_S[17][14]*mom14
        +invM_S[17][17]*mom17
        +invM_S[17][18]*mom18;
        
        f[18] -= invM_S[18][1]*mom1
        +invM_S[18][2]*mom2
        +invM_S[18][5]*mom5
        +invM_S[18][6]*mom6
        +invM_S[18][7]*mom7
        +invM_S[18][8]*mom8
        +invM_S[18][9]*mom9
        +invM_S[18][10]*mom10
        +invM_S[18][14]*mom14
        +invM_S[18][17]*mom17
        +invM_S[18][18]*mom18;
        
        return jSqr;
    }
    
    /// MRT collision step
    static T smagorinskyMrtCollision( Array<T,Descriptor::q>& f,
                           const T &rhoBar, const Array<T,3> & j,
                           const T invM_S[19][19], const Array<T,6> &strain, T cSmago )
    {
        
        Array<T,19> moments, momentsEq;

        computeMoments(moments,f);
        
        T jSqr = VectorTemplateImpl<T,3>::normSqr(j);
        
        computeSmagorinskyEquilibrium( momentsEq, rhoBar, j, jSqr, strain, cSmago );
        
        T mom1 = moments[1] - momentsEq[1];
        T mom2 = moments[2] - momentsEq[2];
        T mom3 = moments[3] - momentsEq[3];
        T mom4 = moments[4] - momentsEq[4];
        T mom5 = moments[5] - momentsEq[5];
        T mom6 = moments[6] - momentsEq[6];
        T mom7 = moments[7] - momentsEq[7];
        T mom8 = moments[8] - momentsEq[8];
        T mom9 = moments[9] - momentsEq[9];
        T mom10 = moments[10] - momentsEq[10];
        T mom11 = moments[11] - momentsEq[11];
        T mom12 = moments[12] - momentsEq[12];
        T mom13 = moments[13] - momentsEq[13];
        T mom14 = moments[14] - momentsEq[14];
        T mom15 = moments[15] - momentsEq[15];
        T mom16 = moments[16];
        T mom17 = moments[17];
        T mom18 = moments[18];

        
        f[0] -= invM_S[0][1]*mom1
                +invM_S[0][2]*mom2;
        
        f[1] -= invM_S[1][1]*mom1
                +invM_S[1][2]*mom2
                +invM_S[1][3]*mom3
                +invM_S[1][4]*mom4
                +invM_S[1][9]*mom9
                +invM_S[1][10]*mom10;
        
        f[2] -= invM_S[2][1]*mom1
                +invM_S[2][2]*mom2
                +invM_S[2][5]*mom5
                +invM_S[2][6]*mom6
                +invM_S[2][9]*mom9
                +invM_S[2][10]*mom10
                +invM_S[2][11]*mom11
                +invM_S[2][12]*mom12;
        
        f[3] -= invM_S[3][1]*mom1
                +invM_S[3][2]*mom2
                +invM_S[3][7]*mom7
                +invM_S[3][8]*mom8
                +invM_S[3][9]*mom9
                +invM_S[3][10]*mom10
                +invM_S[3][11]*mom11
                +invM_S[3][12]*mom12;
        
        f[4] -= invM_S[4][1]*mom1
                +invM_S[4][2]*mom2
                +invM_S[4][3]*mom3
                +invM_S[4][4]*mom4
                +invM_S[4][5]*mom5
                +invM_S[4][6]*mom6
                +invM_S[4][9]*mom9
                +invM_S[4][10]*mom10
                +invM_S[4][11]*mom11
                +invM_S[4][12]*mom12
                +invM_S[4][13]*mom13
                +invM_S[4][16]*mom16
                +invM_S[4][17]*mom17;
        
        f[5] -= invM_S[5][1]*mom1
                +invM_S[5][2]*mom2
                +invM_S[5][3]*mom3
                +invM_S[5][4]*mom4
                +invM_S[5][5]*mom5
                +invM_S[5][6]*mom6
                +invM_S[5][9]*mom9
                +invM_S[5][10]*mom10
                +invM_S[5][11]*mom11
                +invM_S[5][12]*mom12
                +invM_S[5][13]*mom13
                +invM_S[5][16]*mom16
                +invM_S[5][17]*mom17;
        
        f[6] -= invM_S[6][1]*mom1
                +invM_S[6][2]*mom2
                +invM_S[6][3]*mom3
                +invM_S[6][4]*mom4
                +invM_S[6][7]*mom7
                +invM_S[6][8]*mom8
                +invM_S[6][9]*mom9
                +invM_S[6][10]*mom10
                +invM_S[6][11]*mom11
                +invM_S[6][12]*mom12
                +invM_S[6][15]*mom15
                +invM_S[6][16]*mom16
                +invM_S[6][18]*mom18;
        
        f[7] -= invM_S[7][1]*mom1
                +invM_S[7][2]*mom2
                +invM_S[7][3]*mom3
                +invM_S[7][4]*mom4
                +invM_S[7][7]*mom7
                +invM_S[7][8]*mom8
                +invM_S[7][9]*mom9
                +invM_S[7][10]*mom10
                +invM_S[7][11]*mom11
                +invM_S[7][12]*mom12
                +invM_S[7][15]*mom15
                +invM_S[7][16]*mom16
                +invM_S[7][18]*mom18;
        
        f[8] -= invM_S[8][1]*mom1
                +invM_S[8][2]*mom2
                +invM_S[8][5]*mom5
                +invM_S[8][6]*mom6
                +invM_S[8][7]*mom7
                +invM_S[8][8]*mom8
                +invM_S[8][9]*mom9
                +invM_S[8][10]*mom10
                +invM_S[8][14]*mom14
                +invM_S[8][17]*mom17
                +invM_S[8][18]*mom18;
        
        f[9] -= invM_S[9][1]*mom1
                +invM_S[9][2]*mom2
                +invM_S[9][5]*mom5
                +invM_S[9][6]*mom6
                +invM_S[9][7]*mom7
                +invM_S[9][8]*mom8
                +invM_S[9][9]*mom9
                +invM_S[9][10]*mom10
                +invM_S[9][14]*mom14
                +invM_S[9][17]*mom17
                +invM_S[9][18]*mom18;
        
        f[10] -= invM_S[10][1]*mom1
                +invM_S[10][2]*mom2
                +invM_S[10][3]*mom3
                +invM_S[10][4]*mom4
                +invM_S[10][9]*mom9
                +invM_S[10][10]*mom10;
        
        f[11] -= invM_S[11][1]*mom1
                +invM_S[11][2]*mom2
                +invM_S[11][5]*mom5
                +invM_S[11][6]*mom6
                +invM_S[11][9]*mom9
                +invM_S[11][10]*mom10
                +invM_S[11][11]*mom11
                +invM_S[11][12]*mom12;
        
        f[12] -= invM_S[12][1]*mom1
                +invM_S[12][2]*mom2
                +invM_S[12][7]*mom7
                +invM_S[12][8]*mom8
                +invM_S[12][9]*mom9
                +invM_S[12][10]*mom10
                +invM_S[12][11]*mom11
                +invM_S[12][12]*mom12;
        
        f[13] -= invM_S[13][1]*mom1
                +invM_S[13][2]*mom2
                +invM_S[13][3]*mom3
                +invM_S[13][4]*mom4
                +invM_S[13][5]*mom5
                +invM_S[13][6]*mom6
                +invM_S[13][9]*mom9
                +invM_S[13][10]*mom10
                +invM_S[13][11]*mom11
                +invM_S[13][12]*mom12
                +invM_S[13][13]*mom13
                +invM_S[13][16]*mom16
                +invM_S[13][17]*mom17;
        
        f[14] -= invM_S[14][1]*mom1
                +invM_S[14][2]*mom2
                +invM_S[14][3]*mom3
                +invM_S[14][4]*mom4
                +invM_S[14][5]*mom5
                +invM_S[14][6]*mom6
                +invM_S[14][9]*mom9
                +invM_S[14][10]*mom10
                +invM_S[14][11]*mom11
                +invM_S[14][12]*mom12
                +invM_S[14][13]*mom13
                +invM_S[14][16]*mom16
                +invM_S[14][17]*mom17;
        
        f[15] -= invM_S[15][1]*mom1
                +invM_S[15][2]*mom2
                +invM_S[15][3]*mom3
                +invM_S[15][4]*mom4
                +invM_S[15][7]*mom7
                +invM_S[15][8]*mom8
                +invM_S[15][9]*mom9
                +invM_S[15][10]*mom10
                +invM_S[15][11]*mom11
                +invM_S[15][12]*mom12
                +invM_S[15][15]*mom15
                +invM_S[15][16]*mom16
                +invM_S[15][18]*mom18;
        
        f[16] -= invM_S[16][1]*mom1
                +invM_S[16][2]*mom2
                +invM_S[16][3]*mom3
                +invM_S[16][4]*mom4
                +invM_S[16][7]*mom7
                +invM_S[16][8]*mom8
                +invM_S[16][9]*mom9
                +invM_S[16][10]*mom10
                +invM_S[16][11]*mom11
                +invM_S[16][12]*mom12
                +invM_S[16][15]*mom15
                +invM_S[16][16]*mom16
                +invM_S[16][18]*mom18;
        
        f[17] -= invM_S[17][1]*mom1
                +invM_S[17][2]*mom2
                +invM_S[17][5]*mom5
                +invM_S[17][6]*mom6
                +invM_S[17][7]*mom7
                +invM_S[17][8]*mom8
                +invM_S[17][9]*mom9
                +invM_S[17][10]*mom10
                +invM_S[17][14]*mom14
                +invM_S[17][17]*mom17
                +invM_S[17][18]*mom18;
        
        f[18] -= invM_S[18][1]*mom1
                +invM_S[18][2]*mom2
                +invM_S[18][5]*mom5
                +invM_S[18][6]*mom6
                +invM_S[18][7]*mom7
                +invM_S[18][8]*mom8
                +invM_S[18][9]*mom9
                +invM_S[18][10]*mom10
                +invM_S[18][14]*mom14
                +invM_S[18][17]*mom17
                +invM_S[18][18]*mom18;
        
        return jSqr;
    }
    
    static T variableOmegaMrtCollision( Array<T,Descriptor::q>& f,
                           const T &rhoBar, const Array<T,3> & j,
                           const T invM_S[19][19], T omega )
    {
        
        Array<T,19> moments, momentsEq;

        computeMoments(moments,f);
        
        T jSqr = VectorTemplateImpl<T,3>::normSqr(j);
        
        computeEquilibrium(momentsEq,rhoBar,j,jSqr);
        
        T mom1 = moments[1] - momentsEq[1];
        T mom2 = moments[2] - momentsEq[2];
        T mom3 = moments[3] - momentsEq[3];
        T mom4 = moments[4] - momentsEq[4];
        T mom5 = moments[5] - momentsEq[5];
        T mom6 = moments[6] - momentsEq[6];
        T mom7 = moments[7] - momentsEq[7];
        T mom8 = moments[8] - momentsEq[8];
        T mom9 = moments[9] - momentsEq[9];
        T mom10 = moments[10] - momentsEq[10];
        T mom11 = moments[11] - momentsEq[11];
        T mom12 = moments[12] - momentsEq[12];
        T mom13 = moments[13] - momentsEq[13];
        T mom14 = moments[14] - momentsEq[14];
        T mom15 = moments[15] - momentsEq[15];
        T mom16 = moments[16];
        T mom17 = moments[17];
        T mom18 = moments[18];

        // Multiply shear-viscosity indices by omega.
        mom9  *= omega;
        mom11 *= omega;
        mom13 *= omega;
        mom14 *= omega;
        mom15 *= omega;
        
        f[0] -= invM_S[0][1]*mom1
                +invM_S[0][2]*mom2;
        
        f[1] -= invM_S[1][1]*mom1
                +invM_S[1][2]*mom2
                +invM_S[1][3]*mom3
                +invM_S[1][4]*mom4
                +invM_S[1][9]*mom9
                +invM_S[1][10]*mom10;
        
        f[2] -= invM_S[2][1]*mom1
                +invM_S[2][2]*mom2
                +invM_S[2][5]*mom5
                +invM_S[2][6]*mom6
                +invM_S[2][9]*mom9
                +invM_S[2][10]*mom10
                +invM_S[2][11]*mom11
                +invM_S[2][12]*mom12;
        
        f[3] -= invM_S[3][1]*mom1
                +invM_S[3][2]*mom2
                +invM_S[3][7]*mom7
                +invM_S[3][8]*mom8
                +invM_S[3][9]*mom9
                +invM_S[3][10]*mom10
                +invM_S[3][11]*mom11
                +invM_S[3][12]*mom12;
        
        f[4] -= invM_S[4][1]*mom1
                +invM_S[4][2]*mom2
                +invM_S[4][3]*mom3
                +invM_S[4][4]*mom4
                +invM_S[4][5]*mom5
                +invM_S[4][6]*mom6
                +invM_S[4][9]*mom9
                +invM_S[4][10]*mom10
                +invM_S[4][11]*mom11
                +invM_S[4][12]*mom12
                +invM_S[4][13]*mom13
                +invM_S[4][16]*mom16
                +invM_S[4][17]*mom17;
        
        f[5] -= invM_S[5][1]*mom1
                +invM_S[5][2]*mom2
                +invM_S[5][3]*mom3
                +invM_S[5][4]*mom4
                +invM_S[5][5]*mom5
                +invM_S[5][6]*mom6
                +invM_S[5][9]*mom9
                +invM_S[5][10]*mom10
                +invM_S[5][11]*mom11
                +invM_S[5][12]*mom12
                +invM_S[5][13]*mom13
                +invM_S[5][16]*mom16
                +invM_S[5][17]*mom17;
        
        f[6] -= invM_S[6][1]*mom1
                +invM_S[6][2]*mom2
                +invM_S[6][3]*mom3
                +invM_S[6][4]*mom4
                +invM_S[6][7]*mom7
                +invM_S[6][8]*mom8
                +invM_S[6][9]*mom9
                +invM_S[6][10]*mom10
                +invM_S[6][11]*mom11
                +invM_S[6][12]*mom12
                +invM_S[6][15]*mom15
                +invM_S[6][16]*mom16
                +invM_S[6][18]*mom18;
        
        f[7] -= invM_S[7][1]*mom1
                +invM_S[7][2]*mom2
                +invM_S[7][3]*mom3
                +invM_S[7][4]*mom4
                +invM_S[7][7]*mom7
                +invM_S[7][8]*mom8
                +invM_S[7][9]*mom9
                +invM_S[7][10]*mom10
                +invM_S[7][11]*mom11
                +invM_S[7][12]*mom12
                +invM_S[7][15]*mom15
                +invM_S[7][16]*mom16
                +invM_S[7][18]*mom18;
        
        f[8] -= invM_S[8][1]*mom1
                +invM_S[8][2]*mom2
                +invM_S[8][5]*mom5
                +invM_S[8][6]*mom6
                +invM_S[8][7]*mom7
                +invM_S[8][8]*mom8
                +invM_S[8][9]*mom9
                +invM_S[8][10]*mom10
                +invM_S[8][14]*mom14
                +invM_S[8][17]*mom17
                +invM_S[8][18]*mom18;
        
        f[9] -= invM_S[9][1]*mom1
                +invM_S[9][2]*mom2
                +invM_S[9][5]*mom5
                +invM_S[9][6]*mom6
                +invM_S[9][7]*mom7
                +invM_S[9][8]*mom8
                +invM_S[9][9]*mom9
                +invM_S[9][10]*mom10
                +invM_S[9][14]*mom14
                +invM_S[9][17]*mom17
                +invM_S[9][18]*mom18;
        
        f[10] -= invM_S[10][1]*mom1
                +invM_S[10][2]*mom2
                +invM_S[10][3]*mom3
                +invM_S[10][4]*mom4
                +invM_S[10][9]*mom9
                +invM_S[10][10]*mom10;
        
        f[11] -= invM_S[11][1]*mom1
                +invM_S[11][2]*mom2
                +invM_S[11][5]*mom5
                +invM_S[11][6]*mom6
                +invM_S[11][9]*mom9
                +invM_S[11][10]*mom10
                +invM_S[11][11]*mom11
                +invM_S[11][12]*mom12;
        
        f[12] -= invM_S[12][1]*mom1
                +invM_S[12][2]*mom2
                +invM_S[12][7]*mom7
                +invM_S[12][8]*mom8
                +invM_S[12][9]*mom9
                +invM_S[12][10]*mom10
                +invM_S[12][11]*mom11
                +invM_S[12][12]*mom12;
        
        f[13] -= invM_S[13][1]*mom1
                +invM_S[13][2]*mom2
                +invM_S[13][3]*mom3
                +invM_S[13][4]*mom4
                +invM_S[13][5]*mom5
                +invM_S[13][6]*mom6
                +invM_S[13][9]*mom9
                +invM_S[13][10]*mom10
                +invM_S[13][11]*mom11
                +invM_S[13][12]*mom12
                +invM_S[13][13]*mom13
                +invM_S[13][16]*mom16
                +invM_S[13][17]*mom17;
        
        f[14] -= invM_S[14][1]*mom1
                +invM_S[14][2]*mom2
                +invM_S[14][3]*mom3
                +invM_S[14][4]*mom4
                +invM_S[14][5]*mom5
                +invM_S[14][6]*mom6
                +invM_S[14][9]*mom9
                +invM_S[14][10]*mom10
                +invM_S[14][11]*mom11
                +invM_S[14][12]*mom12
                +invM_S[14][13]*mom13
                +invM_S[14][16]*mom16
                +invM_S[14][17]*mom17;
        
        f[15] -= invM_S[15][1]*mom1
                +invM_S[15][2]*mom2
                +invM_S[15][3]*mom3
                +invM_S[15][4]*mom4
                +invM_S[15][7]*mom7
                +invM_S[15][8]*mom8
                +invM_S[15][9]*mom9
                +invM_S[15][10]*mom10
                +invM_S[15][11]*mom11
                +invM_S[15][12]*mom12
                +invM_S[15][15]*mom15
                +invM_S[15][16]*mom16
                +invM_S[15][18]*mom18;
        
        f[16] -= invM_S[16][1]*mom1
                +invM_S[16][2]*mom2
                +invM_S[16][3]*mom3
                +invM_S[16][4]*mom4
                +invM_S[16][7]*mom7
                +invM_S[16][8]*mom8
                +invM_S[16][9]*mom9
                +invM_S[16][10]*mom10
                +invM_S[16][11]*mom11
                +invM_S[16][12]*mom12
                +invM_S[16][15]*mom15
                +invM_S[16][16]*mom16
                +invM_S[16][18]*mom18;
        
        f[17] -= invM_S[17][1]*mom1
                +invM_S[17][2]*mom2
                +invM_S[17][5]*mom5
                +invM_S[17][6]*mom6
                +invM_S[17][7]*mom7
                +invM_S[17][8]*mom8
                +invM_S[17][9]*mom9
                +invM_S[17][10]*mom10
                +invM_S[17][14]*mom14
                +invM_S[17][17]*mom17
                +invM_S[17][18]*mom18;
        
        f[18] -= invM_S[18][1]*mom1
                +invM_S[18][2]*mom2
                +invM_S[18][5]*mom5
                +invM_S[18][6]*mom6
                +invM_S[18][7]*mom7
                +invM_S[18][8]*mom8
                +invM_S[18][9]*mom9
                +invM_S[18][10]*mom10
                +invM_S[18][14]*mom14
                +invM_S[18][17]*mom17
                +invM_S[18][18]*mom18;
        
        return jSqr;
    }
    
    static void addGuoForce( Array<T,Descriptor::q>& f, const Array<T,Descriptor::d>& force,
                             Array<T,Descriptor::d> const& u,
                             T invM_S[Descriptor::q][Descriptor::q], T amplitude )
    {
        Array<T,Descriptor::q> forcing;
        T g_u   = force[0] * u[0] + force[1] * u[1] + force[2] * u[2];
        
        T mom1 = -(T)19*g_u;
        T mom2 = (T)5.5*g_u;
        T mom3 = -(T)0.5*force[0];
        T mom4 = force[0] / (T)3;
        T mom5 = -(T)0.5*force[1];
        T mom6 = force[1] / (T)3;
        T mom7 = -(T)0.5*force[2];
        T mom8 = force[2] / (T)3;
        T mom9 = -(2*force[0]*u[0]-force[1]*u[1]-force[2]*u[2]);
        T mom10 = (T)0.5*(2*force[0]*u[0]-force[1]*u[1]-force[2]*u[2]);
        T mom11 = -(force[1]*u[1]-force[2]*u[2]);
        T mom12 = (T)0.5*(force[1]*u[1]-force[2]*u[2]);
        T mom13 = -(T)0.5*(force[0]*u[1]+force[1]*u[0]);
        T mom14 = -(T)0.5*(force[2]*u[1]+force[1]*u[2]);
        T mom15 = -(T)0.5*(force[0]*u[2]+force[2]*u[0]);
        
        
        f[0] += invM_S[0][1]*mom1
                +invM_S[0][2]*mom2;
        
        f[1] += invM_S[1][1]*mom1
                +invM_S[1][2]*mom2
                +invM_S[1][3]*mom3
                +invM_S[1][4]*mom4
                +invM_S[1][9]*mom9
                +invM_S[1][10]*mom10;
        
        f[2] += invM_S[2][1]*mom1
                +invM_S[2][2]*mom2
                +invM_S[2][5]*mom5
                +invM_S[2][6]*mom6
                +invM_S[2][9]*mom9
                +invM_S[2][10]*mom10
                +invM_S[2][11]*mom11
                +invM_S[2][12]*mom12;
        
        f[3] += invM_S[3][1]*mom1
                +invM_S[3][2]*mom2
                +invM_S[3][7]*mom7
                +invM_S[3][8]*mom8
                +invM_S[3][9]*mom9
                +invM_S[3][10]*mom10
                +invM_S[3][11]*mom11
                +invM_S[3][12]*mom12;
        
        f[4] += invM_S[4][1]*mom1
                +invM_S[4][2]*mom2
                +invM_S[4][3]*mom3
                +invM_S[4][4]*mom4
                +invM_S[4][5]*mom5
                +invM_S[4][6]*mom6
                +invM_S[4][9]*mom9
                +invM_S[4][10]*mom10
                +invM_S[4][11]*mom11
                +invM_S[4][12]*mom12
                +invM_S[4][13]*mom13;
        
        f[5] += invM_S[5][1]*mom1
                +invM_S[5][2]*mom2
                +invM_S[5][3]*mom3
                +invM_S[5][4]*mom4
                +invM_S[5][5]*mom5
                +invM_S[5][6]*mom6
                +invM_S[5][9]*mom9
                +invM_S[5][10]*mom10
                +invM_S[5][11]*mom11
                +invM_S[5][12]*mom12
                +invM_S[5][13]*mom13;
        
        f[6] += invM_S[6][1]*mom1
                +invM_S[6][2]*mom2
                +invM_S[6][3]*mom3
                +invM_S[6][4]*mom4
                +invM_S[6][7]*mom7
                +invM_S[6][8]*mom8
                +invM_S[6][9]*mom9
                +invM_S[6][10]*mom10
                +invM_S[6][11]*mom11
                +invM_S[6][12]*mom12
                +invM_S[6][15]*mom15;
        
        f[7] += invM_S[7][1]*mom1
                +invM_S[7][2]*mom2
                +invM_S[7][3]*mom3
                +invM_S[7][4]*mom4
                +invM_S[7][7]*mom7
                +invM_S[7][8]*mom8
                +invM_S[7][9]*mom9
                +invM_S[7][10]*mom10
                +invM_S[7][11]*mom11
                +invM_S[7][12]*mom12
                +invM_S[7][15]*mom15;
        
        f[8] += invM_S[8][1]*mom1
                +invM_S[8][2]*mom2
                +invM_S[8][5]*mom5
                +invM_S[8][6]*mom6
                +invM_S[8][7]*mom7
                +invM_S[8][8]*mom8
                +invM_S[8][9]*mom9
                +invM_S[8][10]*mom10
                +invM_S[8][14]*mom14;
        
        f[9] += invM_S[9][1]*mom1
                +invM_S[9][2]*mom2
                +invM_S[9][5]*mom5
                +invM_S[9][6]*mom6
                +invM_S[9][7]*mom7
                +invM_S[9][8]*mom8
                +invM_S[9][9]*mom9
                +invM_S[9][10]*mom10
                +invM_S[9][14]*mom14;
        
        f[10] += invM_S[10][1]*mom1
                +invM_S[10][2]*mom2
                +invM_S[10][3]*mom3
                +invM_S[10][4]*mom4
                +invM_S[10][9]*mom9
                +invM_S[10][10]*mom10;
        
        f[11] += invM_S[11][1]*mom1
                +invM_S[11][2]*mom2
                +invM_S[11][5]*mom5
                +invM_S[11][6]*mom6
                +invM_S[11][9]*mom9
                +invM_S[11][10]*mom10
                +invM_S[11][11]*mom11
                +invM_S[11][12]*mom12;
        
        f[12] += invM_S[12][1]*mom1
                +invM_S[12][2]*mom2
                +invM_S[12][7]*mom7
                +invM_S[12][8]*mom8
                +invM_S[12][9]*mom9
                +invM_S[12][10]*mom10
                +invM_S[12][11]*mom11
                +invM_S[12][12]*mom12;
        
        f[13] += invM_S[13][1]*mom1
                +invM_S[13][2]*mom2
                +invM_S[13][3]*mom3
                +invM_S[13][4]*mom4
                +invM_S[13][5]*mom5
                +invM_S[13][6]*mom6
                +invM_S[13][9]*mom9
                +invM_S[13][10]*mom10
                +invM_S[13][11]*mom11
                +invM_S[13][12]*mom12
                +invM_S[13][13]*mom13;
        
        f[14] += invM_S[14][1]*mom1
                +invM_S[14][2]*mom2
                +invM_S[14][3]*mom3
                +invM_S[14][4]*mom4
                +invM_S[14][5]*mom5
                +invM_S[14][6]*mom6
                +invM_S[14][9]*mom9
                +invM_S[14][10]*mom10
                +invM_S[14][11]*mom11
                +invM_S[14][12]*mom12
                +invM_S[14][13]*mom13;
        
        f[15] += invM_S[15][1]*mom1
                +invM_S[15][2]*mom2
                +invM_S[15][3]*mom3
                +invM_S[15][4]*mom4
                +invM_S[15][7]*mom7
                +invM_S[15][8]*mom8
                +invM_S[15][9]*mom9
                +invM_S[15][10]*mom10
                +invM_S[15][11]*mom11
                +invM_S[15][12]*mom12
                +invM_S[15][15]*mom15;
        
        f[16] += invM_S[16][1]*mom1
                +invM_S[16][2]*mom2
                +invM_S[16][3]*mom3
                +invM_S[16][4]*mom4
                +invM_S[16][7]*mom7
                +invM_S[16][8]*mom8
                +invM_S[16][9]*mom9
                +invM_S[16][10]*mom10
                +invM_S[16][11]*mom11
                +invM_S[16][12]*mom12
                +invM_S[16][15]*mom15;
        
        f[17] += invM_S[17][1]*mom1
                +invM_S[17][2]*mom2
                +invM_S[17][5]*mom5
                +invM_S[17][6]*mom6
                +invM_S[17][7]*mom7
                +invM_S[17][8]*mom8
                +invM_S[17][9]*mom9
                +invM_S[17][10]*mom10
                +invM_S[17][14]*mom14;
        
        f[18] += invM_S[18][1]*mom1
                +invM_S[18][2]*mom2
                +invM_S[18][5]*mom5
                +invM_S[18][6]*mom6
                +invM_S[18][7]*mom7
                +invM_S[18][8]*mom8
                +invM_S[18][9]*mom9
                +invM_S[18][10]*mom10
                +invM_S[18][14]*mom14;
                
                
        static const T oneOver6 = (T)1/(T)6;
        static const T oneOver12 = (T)1/(T)12;
        
        f[0] += -g_u;
        
        f[1] += oneOver6*(force[0]*(-(T)1+2*u[0])-force[1]*u[1]-force[2]*u[2]);
        
        f[2] += -oneOver6*(force[0]*u[0]+force[1]*((T)1-2*u[1])+force[2]*u[2]);
        
        f[3] += -oneOver6*(force[0]*u[0]
                            + force[1]*u[1]
                            + force[2]*((T)1-2*u[2]));
        
        f[4] += oneOver12*(force[0]*(-(T)1+2*u[0]+3*u[1])
                            + force[1]*(-(T)1+2*u[1]+3*u[0])
                            - force[2]*u[2]);
        
        f[5] += oneOver12*( force[0]*(-(T)1+2*u[0]-3*u[1])
                                + force[1]*((T)1+2*u[1]-3*u[0])
                                - force[2]*u[2]);
        
        f[6] += oneOver12*(force[0]*(-(T)1+2*u[0]+3*u[2])
                            - force[1]*u[1]
                            + force[2]*(-(T)1+2*u[2]+3*u[0]));
        
        f[7] += oneOver12*(force[0]*(-(T)1+2*u[0]-3*u[2])
                            - force[1]*u[1]
                            + force[2]*((T)1+2*u[2]-3*u[0]));
        
        f[8] += -oneOver12*(force[0]*u[0]
                                + force[1]*((T)1-2*u[1]-3*u[2])
                                + force[2]*((T)1-2*u[2]-3*u[1]));
        
        f[9] += -oneOver12*(force[0]*u[0]
                                + force[1]*((T)1-2*u[1]+3*u[2])
                                + force[2]*(-(T)1-2*u[2]+3*u[1]));
        
        f[10] += oneOver6*(force[0]*((T)1+2*u[0])
                                -force[1]*u[1]
                                -force[2]*u[2]);
        
        f[11] += -oneOver6*(force[0]*u[0]
                                +force[1]*(-(T)1-2*u[1])
                                +force[2]*u[2]);
        
        f[12] += -oneOver6*(force[0]*u[0]
                                +force[1]*u[1]
                                +force[2]*(-(T)1-2*u[2]));
        
        f[13] += oneOver12*(force[0]*((T)1+2*u[0]+3*u[1])
                                +force[1]*((T)1+2*u[1]+3*u[0])
                                -force[2]*u[2]);
        
        f[14] += oneOver12*(force[0]*((T)1+2*u[0]-3*u[1])
                                +force[1]*(-(T)1+2*u[1]-3*u[0])
                                -force[2]*u[2]);
        
        f[15] += oneOver12*(force[0]*((T)1+2*u[0]+3*u[2])
                                -force[1]*u[1]
                                +force[2]*((T)1+2*u[2]+3*u[0]));
        
        f[16] += oneOver12*(force[0]*((T)1+2*u[0]-3*u[2])
                                -force[1]*u[1]
                                +force[2]*(-(T)1+2*u[2]-3*u[0]));
        
        f[17] += -oneOver12*(force[0]*u[0]
                                +force[1]*(-(T)1-2*u[1]-3*u[2])
                                +force[2]*(-(T)1-2*u[2]-3*u[1]));
        
        f[18] += -oneOver12*(force[0]*u[0]
                                +force[1]*(-(T)1-2*u[1]+3*u[2])
                                +force[2]*((T)1-2*u[2]+3*u[1]));
    }
    
    static void variableOmegaAddGuoForce (
            Array<T,Descriptor::q>& f, const Array<T,Descriptor::d>& force,
            Array<T,Descriptor::d> const& u, T invM_S[Descriptor::q][Descriptor::q],
            T amplitude, T omega )
    {
        Array<T,Descriptor::q> forcing;
        T g_u   = force[0] * u[0] + force[1] * u[1] + force[2] * u[2];
        
        T mom1 = -(T)19*g_u;
        T mom2 = (T)5.5*g_u;
        T mom3 = -(T)0.5*force[0];
        T mom4 = force[0] / (T)3;
        T mom5 = -(T)0.5*force[1];
        T mom6 = force[1] / (T)3;
        T mom7 = -(T)0.5*force[2];
        T mom8 = force[2] / (T)3;
        T mom9 = -(2*force[0]*u[0]-force[1]*u[1]-force[2]*u[2]);
        T mom10 = (T)0.5*(2*force[0]*u[0]-force[1]*u[1]-force[2]*u[2]);
        T mom11 = -(force[1]*u[1]-force[2]*u[2]);
        T mom12 = (T)0.5*(force[1]*u[1]-force[2]*u[2]);
        T mom13 = -(T)0.5*(force[0]*u[1]+force[1]*u[0]);
        T mom14 = -(T)0.5*(force[2]*u[1]+force[1]*u[2]);
        T mom15 = -(T)0.5*(force[0]*u[2]+force[2]*u[0]);
        
        // Multiply shear-viscosity indices by omega.
        mom9  *= omega;
        mom11 *= omega;
        mom13 *= omega;
        mom14 *= omega;
        mom15 *= omega;
        
        f[0] += invM_S[0][1]*mom1
                +invM_S[0][2]*mom2;
        
        f[1] += invM_S[1][1]*mom1
                +invM_S[1][2]*mom2
                +invM_S[1][3]*mom3
                +invM_S[1][4]*mom4
                +invM_S[1][9]*mom9
                +invM_S[1][10]*mom10;
        
        f[2] += invM_S[2][1]*mom1
                +invM_S[2][2]*mom2
                +invM_S[2][5]*mom5
                +invM_S[2][6]*mom6
                +invM_S[2][9]*mom9
                +invM_S[2][10]*mom10
                +invM_S[2][11]*mom11
                +invM_S[2][12]*mom12;
        
        f[3] += invM_S[3][1]*mom1
                +invM_S[3][2]*mom2
                +invM_S[3][7]*mom7
                +invM_S[3][8]*mom8
                +invM_S[3][9]*mom9
                +invM_S[3][10]*mom10
                +invM_S[3][11]*mom11
                +invM_S[3][12]*mom12;
        
        f[4] += invM_S[4][1]*mom1
                +invM_S[4][2]*mom2
                +invM_S[4][3]*mom3
                +invM_S[4][4]*mom4
                +invM_S[4][5]*mom5
                +invM_S[4][6]*mom6
                +invM_S[4][9]*mom9
                +invM_S[4][10]*mom10
                +invM_S[4][11]*mom11
                +invM_S[4][12]*mom12
                +invM_S[4][13]*mom13;
        
        f[5] += invM_S[5][1]*mom1
                +invM_S[5][2]*mom2
                +invM_S[5][3]*mom3
                +invM_S[5][4]*mom4
                +invM_S[5][5]*mom5
                +invM_S[5][6]*mom6
                +invM_S[5][9]*mom9
                +invM_S[5][10]*mom10
                +invM_S[5][11]*mom11
                +invM_S[5][12]*mom12
                +invM_S[5][13]*mom13;
        
        f[6] += invM_S[6][1]*mom1
                +invM_S[6][2]*mom2
                +invM_S[6][3]*mom3
                +invM_S[6][4]*mom4
                +invM_S[6][7]*mom7
                +invM_S[6][8]*mom8
                +invM_S[6][9]*mom9
                +invM_S[6][10]*mom10
                +invM_S[6][11]*mom11
                +invM_S[6][12]*mom12
                +invM_S[6][15]*mom15;
        
        f[7] += invM_S[7][1]*mom1
                +invM_S[7][2]*mom2
                +invM_S[7][3]*mom3
                +invM_S[7][4]*mom4
                +invM_S[7][7]*mom7
                +invM_S[7][8]*mom8
                +invM_S[7][9]*mom9
                +invM_S[7][10]*mom10
                +invM_S[7][11]*mom11
                +invM_S[7][12]*mom12
                +invM_S[7][15]*mom15;
        
        f[8] += invM_S[8][1]*mom1
                +invM_S[8][2]*mom2
                +invM_S[8][5]*mom5
                +invM_S[8][6]*mom6
                +invM_S[8][7]*mom7
                +invM_S[8][8]*mom8
                +invM_S[8][9]*mom9
                +invM_S[8][10]*mom10
                +invM_S[8][14]*mom14;
        
        f[9] += invM_S[9][1]*mom1
                +invM_S[9][2]*mom2
                +invM_S[9][5]*mom5
                +invM_S[9][6]*mom6
                +invM_S[9][7]*mom7
                +invM_S[9][8]*mom8
                +invM_S[9][9]*mom9
                +invM_S[9][10]*mom10
                +invM_S[9][14]*mom14;
        
        f[10] += invM_S[10][1]*mom1
                +invM_S[10][2]*mom2
                +invM_S[10][3]*mom3
                +invM_S[10][4]*mom4
                +invM_S[10][9]*mom9
                +invM_S[10][10]*mom10;
        
        f[11] += invM_S[11][1]*mom1
                +invM_S[11][2]*mom2
                +invM_S[11][5]*mom5
                +invM_S[11][6]*mom6
                +invM_S[11][9]*mom9
                +invM_S[11][10]*mom10
                +invM_S[11][11]*mom11
                +invM_S[11][12]*mom12;
        
        f[12] += invM_S[12][1]*mom1
                +invM_S[12][2]*mom2
                +invM_S[12][7]*mom7
                +invM_S[12][8]*mom8
                +invM_S[12][9]*mom9
                +invM_S[12][10]*mom10
                +invM_S[12][11]*mom11
                +invM_S[12][12]*mom12;
        
        f[13] += invM_S[13][1]*mom1
                +invM_S[13][2]*mom2
                +invM_S[13][3]*mom3
                +invM_S[13][4]*mom4
                +invM_S[13][5]*mom5
                +invM_S[13][6]*mom6
                +invM_S[13][9]*mom9
                +invM_S[13][10]*mom10
                +invM_S[13][11]*mom11
                +invM_S[13][12]*mom12
                +invM_S[13][13]*mom13;
        
        f[14] += invM_S[14][1]*mom1
                +invM_S[14][2]*mom2
                +invM_S[14][3]*mom3
                +invM_S[14][4]*mom4
                +invM_S[14][5]*mom5
                +invM_S[14][6]*mom6
                +invM_S[14][9]*mom9
                +invM_S[14][10]*mom10
                +invM_S[14][11]*mom11
                +invM_S[14][12]*mom12
                +invM_S[14][13]*mom13;
        
        f[15] += invM_S[15][1]*mom1
                +invM_S[15][2]*mom2
                +invM_S[15][3]*mom3
                +invM_S[15][4]*mom4
                +invM_S[15][7]*mom7
                +invM_S[15][8]*mom8
                +invM_S[15][9]*mom9
                +invM_S[15][10]*mom10
                +invM_S[15][11]*mom11
                +invM_S[15][12]*mom12
                +invM_S[15][15]*mom15;
        
        f[16] += invM_S[16][1]*mom1
                +invM_S[16][2]*mom2
                +invM_S[16][3]*mom3
                +invM_S[16][4]*mom4
                +invM_S[16][7]*mom7
                +invM_S[16][8]*mom8
                +invM_S[16][9]*mom9
                +invM_S[16][10]*mom10
                +invM_S[16][11]*mom11
                +invM_S[16][12]*mom12
                +invM_S[16][15]*mom15;
        
        f[17] += invM_S[17][1]*mom1
                +invM_S[17][2]*mom2
                +invM_S[17][5]*mom5
                +invM_S[17][6]*mom6
                +invM_S[17][7]*mom7
                +invM_S[17][8]*mom8
                +invM_S[17][9]*mom9
                +invM_S[17][10]*mom10
                +invM_S[17][14]*mom14;
        
        f[18] += invM_S[18][1]*mom1
                +invM_S[18][2]*mom2
                +invM_S[18][5]*mom5
                +invM_S[18][6]*mom6
                +invM_S[18][7]*mom7
                +invM_S[18][8]*mom8
                +invM_S[18][9]*mom9
                +invM_S[18][10]*mom10
                +invM_S[18][14]*mom14;
                
                
        static const T oneOver6 = (T)1/(T)6;
        static const T oneOver12 = (T)1/(T)12;
        
        f[0] += -g_u;
        
        f[1] += oneOver6*(force[0]*(-(T)1+2*u[0])-force[1]*u[1]-force[2]*u[2]);
        
        f[2] += -oneOver6*(force[0]*u[0]+force[1]*((T)1-2*u[1])+force[2]*u[2]);
        
        f[3] += -oneOver6*(force[0]*u[0]
                            + force[1]*u[1]
                            + force[2]*((T)1-2*u[2]));
        
        f[4] += oneOver12*(force[0]*(-(T)1+2*u[0]+3*u[1])
                            + force[1]*(-(T)1+2*u[1]+3*u[0])
                            - force[2]*u[2]);
        
        f[5] += oneOver12*( force[0]*(-(T)1+2*u[0]-3*u[1])
                                + force[1]*((T)1+2*u[1]-3*u[0])
                                - force[2]*u[2]);
        
        f[6] += oneOver12*(force[0]*(-(T)1+2*u[0]+3*u[2])
                            - force[1]*u[1]
                            + force[2]*(-(T)1+2*u[2]+3*u[0]));
        
        f[7] += oneOver12*(force[0]*(-(T)1+2*u[0]-3*u[2])
                            - force[1]*u[1]
                            + force[2]*((T)1+2*u[2]-3*u[0]));
        
        f[8] += -oneOver12*(force[0]*u[0]
                                + force[1]*((T)1-2*u[1]-3*u[2])
                                + force[2]*((T)1-2*u[2]-3*u[1]));
        
        f[9] += -oneOver12*(force[0]*u[0]
                                + force[1]*((T)1-2*u[1]+3*u[2])
                                + force[2]*(-(T)1-2*u[2]+3*u[1]));
        
        f[10] += oneOver6*(force[0]*((T)1+2*u[0])
                                -force[1]*u[1]
                                -force[2]*u[2]);
        
        f[11] += -oneOver6*(force[0]*u[0]
                                +force[1]*(-(T)1-2*u[1])
                                +force[2]*u[2]);
        
        f[12] += -oneOver6*(force[0]*u[0]
                                +force[1]*u[1]
                                +force[2]*(-(T)1-2*u[2]));
        
        f[13] += oneOver12*(force[0]*((T)1+2*u[0]+3*u[1])
                                +force[1]*((T)1+2*u[1]+3*u[0])
                                -force[2]*u[2]);
        
        f[14] += oneOver12*(force[0]*((T)1+2*u[0]-3*u[1])
                                +force[1]*(-(T)1+2*u[1]-3*u[0])
                                -force[2]*u[2]);
        
        f[15] += oneOver12*(force[0]*((T)1+2*u[0]+3*u[2])
                                -force[1]*u[1]
                                +force[2]*((T)1+2*u[2]+3*u[0]));
        
        f[16] += oneOver12*(force[0]*((T)1+2*u[0]-3*u[2])
                                -force[1]*u[1]
                                +force[2]*(-(T)1+2*u[2]-3*u[0]));
        
        f[17] += -oneOver12*(force[0]*u[0]
                                +force[1]*(-(T)1-2*u[1]-3*u[2])
                                +force[2]*(-(T)1-2*u[2]-3*u[1]));
        
        f[18] += -oneOver12*(force[0]*u[0]
                                +force[1]*(-(T)1-2*u[1]+3*u[2])
                                +force[2]*((T)1-2*u[2]+3*u[1]));
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
    
    /// Smagorinsky MRT collision step
    static T smagorinskyMrtCollisionWithForce( Array<T,Descriptor::q>& f,
                                                       const T &rhoBar, const Array<T,Descriptor::d> & u,
                                                       T invM_S[Descriptor::q][Descriptor::q], 
                                                       const Array<T,SymmetricTensorImpl<T,Descriptor::d>::n > &strain, T cSmago, 
                                                       const Array<T,Descriptor::d> &force, T amplitude) 
    {
        Array<T,Descriptor::d> j = Descriptor::fullRho(rhoBar)*u;
        T jSqr = smagorinskyMrtCollision( f, rhoBar, j, invM_S, strain, cSmago );
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
    
    /// Smagorinsky MRT collision step
    static T quasiIncSmagorinskyMrtCollisionWithForce( Array<T,Descriptor::q>& f,
                                                       const T &rhoBar, const Array<T,Descriptor::d> & u,
                                                       T invM_S[Descriptor::q][Descriptor::q], 
                                                       const Array<T,SymmetricTensorImpl<T,Descriptor::d>::n > &strain, T cSmago, 
                                                       const Array<T,Descriptor::d> &force, T amplitude) 
    {
        Array<T,Descriptor::d> j = Descriptor::fullRho(rhoBar)*u;
        T jSqr = quasiIncSmagorinskyMrtCollision( f, rhoBar, j, invM_S, strain, cSmago );
        addGuoForce( f, force, u, invM_S, amplitude );
        
        return jSqr;
    }
    
    static T variableOmegaMrtCollisionWithForce (
            Array<T,Descriptor::q>& f, const T &rhoBar,
            const Array<T,Descriptor::d> & u, T invM_S[Descriptor::q][Descriptor::q], 
            const Array<T,Descriptor::d> &force, T amplitude, T omega ) 
    {
        Array<T,Descriptor::d> j = Descriptor::fullRho(rhoBar)*u;
        T jSqr = variableOmegaMrtCollision( f, rhoBar, j, invM_S, omega );
        variableOmegaAddGuoForce( f, force, u, invM_S, amplitude, omega );
        
        return jSqr;
    }
    
    /// Computation of all equilibrium distribution (in moments space)
    static void computeQuasiIncEquilibrium( Array<T,Descriptor::q>& momentsEq,
                                            T rhoBar, Array<T,3> const& j, T jSqr )
    {
        momentsEq[0] = rhoBar;
        momentsEq[1] = (T)19*jSqr-(T)11*rhoBar;
        momentsEq[2] = -(T)5.5*jSqr+(T)3*rhoBar;
        momentsEq[3] = j[0];
        momentsEq[4] = -((T)2/(T)3)*j[0];
        momentsEq[5] = j[1];
        momentsEq[6] = -((T)2/3)*j[1];
        momentsEq[7] = j[2];
        momentsEq[8] = -((T)2/(T)3)*j[2];
        momentsEq[9] = ((T)2*j[0]*j[0]-j[1]*j[1]-j[2]*j[2]);
        momentsEq[10] = (-j[0]*j[0]+(T)0.5*j[1]*j[1]+(T)0.5*j[2]*j[2]);
        momentsEq[11] = (j[1]*j[1]-j[2]*j[2]);
        momentsEq[12] = (-(T)0.5*j[1]*j[1]+(T)0.5*j[2]*j[2]);
        momentsEq[13] = j[1]*j[0];
        momentsEq[14] = j[2]*j[1];
        momentsEq[15] = j[2]*j[0];
        momentsEq[16] = T();
        momentsEq[17] = T();
        momentsEq[18] = T();
    }
    
    /// Computation of all equilibrium distribution (in moments space)
    static void computeQuasiIncSmagorinskyEquilibrium( Array<T,Descriptor::q>& momentsEq,
                                                       T rhoBar, Array<T,3> const& j, T jSqr, const Array<T,6> &strain, T cSmago )
    {
        typedef SymmetricTensorImpl<T,3> S;
        T sNorm = sqrt((T)2*SymmetricTensorImpl<T,3>::tensorNormSqr(strain));
        T smagoFactor = (T)2*cSmago * cSmago * sNorm;
        
        T ux2 = smagoFactor * strain[S::xx];
        T uy2 = smagoFactor * strain[S::yy];
        T uz2 = smagoFactor * strain[S::zz];
        
        T uxuy = smagoFactor * strain[S::xy];
        T uyuz = smagoFactor * strain[S::yz];
        T uxuz = smagoFactor * strain[S::xz];
        
        momentsEq[0] = rhoBar;
        momentsEq[1] = (T)19*(jSqr + ux2+uy2+uz2) -(T)11*rhoBar;
        momentsEq[2] = -(T)5.5*(jSqr + ux2+uy2+uz2)+(T)3*rhoBar;
        momentsEq[3] = j[0];
        momentsEq[4] = -((T)2/(T)3)*j[0];
        momentsEq[5] = j[1];
        momentsEq[6] = -((T)2/3)*j[1];
        momentsEq[7] = j[2];
        momentsEq[8] = -((T)2/(T)3)*j[2];
        momentsEq[9] = ((T)2*(j[0]*j[0]+ux2)-(j[1]*j[1]+uy2)-(j[2]*j[2]+uz2));
        momentsEq[10] = (-(j[0]*j[0]+ux2)+(T)0.5*(j[1]*j[1]+uy2)+(T)0.5*(j[2]*j[2]+uz2));
        momentsEq[11] = (j[1]*j[1]+uy2-(j[2]*j[2]+uz2));
        momentsEq[12] = (-(T)0.5*(j[1]*j[1]+uy2)+(T)0.5*(j[2]*j[2]+uz2));
        momentsEq[13] = j[1]*j[0]+uxuy;
        momentsEq[14] = j[2]*j[1]+uyuz;
        momentsEq[15] = j[2]*j[0]+uxuz;
        momentsEq[16] = T();
        momentsEq[17] = T();
        momentsEq[18] = T();
    }
    
    /// MRT collision step
    static T quasiIncMrtCollision( Array<T,Descriptor::q>& f,
                                   const T &rhoBar, const Array<T,3> & j,
                                   const T invM_S[19][19] )
    {
        
        Array<T,19> moments, momentsEq;

        computeMoments(moments,f);
        moments[0] = rhoBar;
        moments[3] = j[0];
        moments[5] = j[1];
        moments[7] = j[2];
        
        T jSqr = VectorTemplateImpl<T,3>::normSqr(j);
        
        computeQuasiIncEquilibrium(momentsEq,rhoBar,j,jSqr);
        
        T mom1 = moments[1] - momentsEq[1];
        T mom2 = moments[2] - momentsEq[2];
        T mom4 = moments[4] - momentsEq[4];
        T mom6 = moments[6] - momentsEq[6];
        T mom8 = moments[8] - momentsEq[8];
        T mom9 = moments[9] - momentsEq[9];
        T mom10 = moments[10] - momentsEq[10];
        T mom11 = moments[11] - momentsEq[11];
        T mom12 = moments[12] - momentsEq[12];
        T mom13 = moments[13] - momentsEq[13];
        T mom14 = moments[14] - momentsEq[14];
        T mom15 = moments[15] - momentsEq[15];
        T mom16 = moments[16];
        T mom17 = moments[17];
        T mom18 = moments[18];

        
        f[0] -= invM_S[0][1]*mom1
                +invM_S[0][2]*mom2;
        
        f[1] -= invM_S[1][1]*mom1
                +invM_S[1][2]*mom2
                +invM_S[1][4]*mom4
                +invM_S[1][9]*mom9
                +invM_S[1][10]*mom10;
        
        f[2] -= invM_S[2][1]*mom1
                +invM_S[2][2]*mom2
                +invM_S[2][6]*mom6
                +invM_S[2][9]*mom9
                +invM_S[2][10]*mom10
                +invM_S[2][11]*mom11
                +invM_S[2][12]*mom12;
        
        f[3] -= invM_S[3][1]*mom1
                +invM_S[3][2]*mom2
                +invM_S[3][8]*mom8
                +invM_S[3][9]*mom9
                +invM_S[3][10]*mom10
                +invM_S[3][11]*mom11
                +invM_S[3][12]*mom12;
        
        f[4] -= invM_S[4][1]*mom1
                +invM_S[4][2]*mom2
                +invM_S[4][4]*mom4
                +invM_S[4][6]*mom6
                +invM_S[4][9]*mom9
                +invM_S[4][10]*mom10
                +invM_S[4][11]*mom11
                +invM_S[4][12]*mom12
                +invM_S[4][13]*mom13
                +invM_S[4][16]*mom16
                +invM_S[4][17]*mom17;
        
        f[5] -= invM_S[5][1]*mom1
                +invM_S[5][2]*mom2
                +invM_S[5][4]*mom4
                +invM_S[5][6]*mom6
                +invM_S[5][9]*mom9
                +invM_S[5][10]*mom10
                +invM_S[5][11]*mom11
                +invM_S[5][12]*mom12
                +invM_S[5][13]*mom13
                +invM_S[5][16]*mom16
                +invM_S[5][17]*mom17;
        
        f[6] -= invM_S[6][1]*mom1
                +invM_S[6][2]*mom2
                +invM_S[6][4]*mom4
                +invM_S[6][8]*mom8
                +invM_S[6][9]*mom9
                +invM_S[6][10]*mom10
                +invM_S[6][11]*mom11
                +invM_S[6][12]*mom12
                +invM_S[6][15]*mom15
                +invM_S[6][16]*mom16
                +invM_S[6][18]*mom18;
        
        f[7] -= invM_S[7][1]*mom1
                +invM_S[7][2]*mom2
                +invM_S[7][4]*mom4
                +invM_S[7][8]*mom8
                +invM_S[7][9]*mom9
                +invM_S[7][10]*mom10
                +invM_S[7][11]*mom11
                +invM_S[7][12]*mom12
                +invM_S[7][15]*mom15
                +invM_S[7][16]*mom16
                +invM_S[7][18]*mom18;
        
        f[8] -= invM_S[8][1]*mom1
                +invM_S[8][2]*mom2
                +invM_S[8][6]*mom6
                +invM_S[8][8]*mom8
                +invM_S[8][9]*mom9
                +invM_S[8][10]*mom10
                +invM_S[8][14]*mom14
                +invM_S[8][17]*mom17
                +invM_S[8][18]*mom18;
        
        f[9] -= invM_S[9][1]*mom1
                +invM_S[9][2]*mom2
                +invM_S[9][6]*mom6
                +invM_S[9][8]*mom8
                +invM_S[9][9]*mom9
                +invM_S[9][10]*mom10
                +invM_S[9][14]*mom14
                +invM_S[9][17]*mom17
                +invM_S[9][18]*mom18;
        
        f[10] -= invM_S[10][1]*mom1
                +invM_S[10][2]*mom2
                +invM_S[10][4]*mom4
                +invM_S[10][9]*mom9
                +invM_S[10][10]*mom10;
        
        f[11] -= invM_S[11][1]*mom1
                +invM_S[11][2]*mom2
                +invM_S[11][6]*mom6
                +invM_S[11][9]*mom9
                +invM_S[11][10]*mom10
                +invM_S[11][11]*mom11
                +invM_S[11][12]*mom12;
        
        f[12] -= invM_S[12][1]*mom1
                +invM_S[12][2]*mom2
                +invM_S[12][8]*mom8
                +invM_S[12][9]*mom9
                +invM_S[12][10]*mom10
                +invM_S[12][11]*mom11
                +invM_S[12][12]*mom12;
        
        f[13] -= invM_S[13][1]*mom1
                +invM_S[13][2]*mom2
                +invM_S[13][4]*mom4
                +invM_S[13][6]*mom6
                +invM_S[13][9]*mom9
                +invM_S[13][10]*mom10
                +invM_S[13][11]*mom11
                +invM_S[13][12]*mom12
                +invM_S[13][13]*mom13
                +invM_S[13][16]*mom16
                +invM_S[13][17]*mom17;
        
        f[14] -= invM_S[14][1]*mom1
                +invM_S[14][2]*mom2
                +invM_S[14][4]*mom4
                +invM_S[14][6]*mom6
                +invM_S[14][9]*mom9
                +invM_S[14][10]*mom10
                +invM_S[14][11]*mom11
                +invM_S[14][12]*mom12
                +invM_S[14][13]*mom13
                +invM_S[14][16]*mom16
                +invM_S[14][17]*mom17;
        
        f[15] -= invM_S[15][1]*mom1
                +invM_S[15][2]*mom2
                +invM_S[15][4]*mom4
                +invM_S[15][8]*mom8
                +invM_S[15][9]*mom9
                +invM_S[15][10]*mom10
                +invM_S[15][11]*mom11
                +invM_S[15][12]*mom12
                +invM_S[15][15]*mom15
                +invM_S[15][16]*mom16
                +invM_S[15][18]*mom18;
        
        f[16] -= invM_S[16][1]*mom1
                +invM_S[16][2]*mom2
                +invM_S[16][4]*mom4
                +invM_S[16][8]*mom8
                +invM_S[16][9]*mom9
                +invM_S[16][10]*mom10
                +invM_S[16][11]*mom11
                +invM_S[16][12]*mom12
                +invM_S[16][15]*mom15
                +invM_S[16][16]*mom16
                +invM_S[16][18]*mom18;
        
        f[17] -= invM_S[17][1]*mom1
                +invM_S[17][2]*mom2
                +invM_S[17][6]*mom6
                +invM_S[17][8]*mom8
                +invM_S[17][9]*mom9
                +invM_S[17][10]*mom10
                +invM_S[17][14]*mom14
                +invM_S[17][17]*mom17
                +invM_S[17][18]*mom18;
        
        f[18] -= invM_S[18][1]*mom1
                +invM_S[18][2]*mom2
                +invM_S[18][6]*mom6
                +invM_S[18][8]*mom8
                +invM_S[18][9]*mom9
                +invM_S[18][10]*mom10
                +invM_S[18][14]*mom14
                +invM_S[18][17]*mom17
                +invM_S[18][18]*mom18;
        
        return jSqr;
    }
    
    
    /// MRT collision step
    static T quasiIncSmagorinskyMrtCollision( Array<T,Descriptor::q>& f,
                                              const T &rhoBar, const Array<T,3> & j,
                                              const T invM_S[19][19], const Array<T,6> &strain, T cSmago )
    {
        
        Array<T,19> moments, momentsEq;

        computeMoments(moments,f);
        moments[0] = rhoBar;
        moments[3] = j[0];
        moments[5] = j[1];
        moments[7] = j[2];
        
        T jSqr = VectorTemplateImpl<T,3>::normSqr(j);
        
        computeQuasiIncSmagorinskyEquilibrium(momentsEq,rhoBar,j,jSqr,strain,cSmago);
        
        T mom1 = moments[1] - momentsEq[1];
        T mom2 = moments[2] - momentsEq[2];
        T mom4 = moments[4] - momentsEq[4];
        T mom6 = moments[6] - momentsEq[6];
        T mom8 = moments[8] - momentsEq[8];
        T mom9 = moments[9] - momentsEq[9];
        T mom10 = moments[10] - momentsEq[10];
        T mom11 = moments[11] - momentsEq[11];
        T mom12 = moments[12] - momentsEq[12];
        T mom13 = moments[13] - momentsEq[13];
        T mom14 = moments[14] - momentsEq[14];
        T mom15 = moments[15] - momentsEq[15];
        T mom16 = moments[16];
        T mom17 = moments[17];
        T mom18 = moments[18];

        
        f[0] -= invM_S[0][1]*mom1
                +invM_S[0][2]*mom2;
        
        f[1] -= invM_S[1][1]*mom1
                +invM_S[1][2]*mom2
                +invM_S[1][4]*mom4
                +invM_S[1][9]*mom9
                +invM_S[1][10]*mom10;
        
        f[2] -= invM_S[2][1]*mom1
                +invM_S[2][2]*mom2
                +invM_S[2][6]*mom6
                +invM_S[2][9]*mom9
                +invM_S[2][10]*mom10
                +invM_S[2][11]*mom11
                +invM_S[2][12]*mom12;
        
        f[3] -= invM_S[3][1]*mom1
                +invM_S[3][2]*mom2
                +invM_S[3][8]*mom8
                +invM_S[3][9]*mom9
                +invM_S[3][10]*mom10
                +invM_S[3][11]*mom11
                +invM_S[3][12]*mom12;
        
        f[4] -= invM_S[4][1]*mom1
                +invM_S[4][2]*mom2
                +invM_S[4][4]*mom4
                +invM_S[4][6]*mom6
                +invM_S[4][9]*mom9
                +invM_S[4][10]*mom10
                +invM_S[4][11]*mom11
                +invM_S[4][12]*mom12
                +invM_S[4][13]*mom13
                +invM_S[4][16]*mom16
                +invM_S[4][17]*mom17;
        
        f[5] -= invM_S[5][1]*mom1
                +invM_S[5][2]*mom2
                +invM_S[5][4]*mom4
                +invM_S[5][6]*mom6
                +invM_S[5][9]*mom9
                +invM_S[5][10]*mom10
                +invM_S[5][11]*mom11
                +invM_S[5][12]*mom12
                +invM_S[5][13]*mom13
                +invM_S[5][16]*mom16
                +invM_S[5][17]*mom17;
        
        f[6] -= invM_S[6][1]*mom1
                +invM_S[6][2]*mom2
                +invM_S[6][4]*mom4
                +invM_S[6][8]*mom8
                +invM_S[6][9]*mom9
                +invM_S[6][10]*mom10
                +invM_S[6][11]*mom11
                +invM_S[6][12]*mom12
                +invM_S[6][15]*mom15
                +invM_S[6][16]*mom16
                +invM_S[6][18]*mom18;
        
        f[7] -= invM_S[7][1]*mom1
                +invM_S[7][2]*mom2
                +invM_S[7][4]*mom4
                +invM_S[7][8]*mom8
                +invM_S[7][9]*mom9
                +invM_S[7][10]*mom10
                +invM_S[7][11]*mom11
                +invM_S[7][12]*mom12
                +invM_S[7][15]*mom15
                +invM_S[7][16]*mom16
                +invM_S[7][18]*mom18;
        
        f[8] -= invM_S[8][1]*mom1
                +invM_S[8][2]*mom2
                +invM_S[8][6]*mom6
                +invM_S[8][8]*mom8
                +invM_S[8][9]*mom9
                +invM_S[8][10]*mom10
                +invM_S[8][14]*mom14
                +invM_S[8][17]*mom17
                +invM_S[8][18]*mom18;
        
        f[9] -= invM_S[9][1]*mom1
                +invM_S[9][2]*mom2
                +invM_S[9][6]*mom6
                +invM_S[9][8]*mom8
                +invM_S[9][9]*mom9
                +invM_S[9][10]*mom10
                +invM_S[9][14]*mom14
                +invM_S[9][17]*mom17
                +invM_S[9][18]*mom18;
        
        f[10] -= invM_S[10][1]*mom1
                +invM_S[10][2]*mom2
                +invM_S[10][4]*mom4
                +invM_S[10][9]*mom9
                +invM_S[10][10]*mom10;
        
        f[11] -= invM_S[11][1]*mom1
                +invM_S[11][2]*mom2
                +invM_S[11][6]*mom6
                +invM_S[11][9]*mom9
                +invM_S[11][10]*mom10
                +invM_S[11][11]*mom11
                +invM_S[11][12]*mom12;
        
        f[12] -= invM_S[12][1]*mom1
                +invM_S[12][2]*mom2
                +invM_S[12][8]*mom8
                +invM_S[12][9]*mom9
                +invM_S[12][10]*mom10
                +invM_S[12][11]*mom11
                +invM_S[12][12]*mom12;
        
        f[13] -= invM_S[13][1]*mom1
                +invM_S[13][2]*mom2
                +invM_S[13][4]*mom4
                +invM_S[13][6]*mom6
                +invM_S[13][9]*mom9
                +invM_S[13][10]*mom10
                +invM_S[13][11]*mom11
                +invM_S[13][12]*mom12
                +invM_S[13][13]*mom13
                +invM_S[13][16]*mom16
                +invM_S[13][17]*mom17;
        
        f[14] -= invM_S[14][1]*mom1
                +invM_S[14][2]*mom2
                +invM_S[14][4]*mom4
                +invM_S[14][6]*mom6
                +invM_S[14][9]*mom9
                +invM_S[14][10]*mom10
                +invM_S[14][11]*mom11
                +invM_S[14][12]*mom12
                +invM_S[14][13]*mom13
                +invM_S[14][16]*mom16
                +invM_S[14][17]*mom17;
        
        f[15] -= invM_S[15][1]*mom1
                +invM_S[15][2]*mom2
                +invM_S[15][4]*mom4
                +invM_S[15][8]*mom8
                +invM_S[15][9]*mom9
                +invM_S[15][10]*mom10
                +invM_S[15][11]*mom11
                +invM_S[15][12]*mom12
                +invM_S[15][15]*mom15
                +invM_S[15][16]*mom16
                +invM_S[15][18]*mom18;
        
        f[16] -= invM_S[16][1]*mom1
                +invM_S[16][2]*mom2
                +invM_S[16][4]*mom4
                +invM_S[16][8]*mom8
                +invM_S[16][9]*mom9
                +invM_S[16][10]*mom10
                +invM_S[16][11]*mom11
                +invM_S[16][12]*mom12
                +invM_S[16][15]*mom15
                +invM_S[16][16]*mom16
                +invM_S[16][18]*mom18;
        
        f[17] -= invM_S[17][1]*mom1
                +invM_S[17][2]*mom2
                +invM_S[17][6]*mom6
                +invM_S[17][8]*mom8
                +invM_S[17][9]*mom9
                +invM_S[17][10]*mom10
                +invM_S[17][14]*mom14
                +invM_S[17][17]*mom17
                +invM_S[17][18]*mom18;
        
        f[18] -= invM_S[18][1]*mom1
                +invM_S[18][2]*mom2
                +invM_S[18][6]*mom6
                +invM_S[18][8]*mom8
                +invM_S[18][9]*mom9
                +invM_S[18][10]*mom10
                +invM_S[18][14]*mom14
                +invM_S[18][17]*mom17
                +invM_S[18][18]*mom18;
        
        return jSqr;
    }

};


}  // namespace plb

#endif
