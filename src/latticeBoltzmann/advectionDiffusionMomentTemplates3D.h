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

/* Main author: Orestis Malaspinas
 */

/** \file
 * Helper functions for the computation of velocity moments for the f's.
 * This file is all about efficiency. The generic template code is specialized
 * for commonly used Lattices, so that a maximum performance can be taken out
 * of each case.
 */
#ifndef ADVECTION_DIFFUSION_MOMENTS_TEMPLATES_3D_H
#define ADVECTION_DIFFUSION_MOMENTS_TEMPLATES_3D_H

namespace plb {

// This structure forwards the calls to the appropriate helper class
template<typename T>
struct advectionDiffusionMomentTemplatesImpl<T,descriptors::D3Q7DescriptorBase<T> >
{
    
typedef descriptors::D3Q7DescriptorBase<T> Descriptor;
    
static void get_rhoBar_jEq(Array<T,Descriptor::q> const& f, T& rhoBar, 
                           Array<T,Descriptor::d>& jEq, const T u[Descriptor::d] )
{
    rhoBar = momentTemplatesImpl<T,Descriptor>::get_rhoBar(f);
    T rho = Descriptor::fullRho(rhoBar);
    jEq[0] = rho * u[0];
    jEq[1] = rho * u[1];
    jEq[2] = rho * u[2];
}
    
static void get_jEq(const T& rhoBar, Array<T,Descriptor::d>& jEq, const T *u)
{
    T rho = Descriptor::fullRho(rhoBar);
    jEq[0] = rho * u[0];
    jEq[1] = rho * u[1];
    jEq[2] = rho * u[2];
}
    
static void get_rhoBar_jEq_jNeq_linear (
        Array<T,Descriptor::q> const& f, T& rhoBar, 
        Array<T,Descriptor::d>& jEq, Array<T,Descriptor::d>& jNeq, const T *u )
{
    rhoBar = momentTemplatesImpl<T,Descriptor>::get_rhoBar(f);
    jEq[0] = u[0];
    jEq[1] = u[1];
    jEq[2] = u[2];
    
    jNeq[0] = - f[1] + f[4] - jEq[0];
    jNeq[1] = - f[2] + f[5] - jEq[1];
    jNeq[2] = - f[3] + f[6] - jEq[2];
}
    
static void get_rhoBar_jEq_jNeq(Array<T,Descriptor::q> const& f, T& rhoBar, 
                                Array<T,Descriptor::d>& jEq, Array<T,Descriptor::d>& jNeq,
                                const T *u)
{
    rhoBar = momentTemplatesImpl<T,Descriptor>::get_rhoBar(f);
    T rho = Descriptor::fullRho(rhoBar);
    jEq[0] = rho * u[0];
    jEq[1] = rho * u[1];
    jEq[2] = rho * u[2];
    
    jNeq[0] = - f[1] + f[4] - jEq[0];
    jNeq[1] = - f[2] + f[5] - jEq[1];
    jNeq[2] = - f[3] + f[6] - jEq[2];
}

};  // struct advectionDiffusionMomentTemplatesImpl

}  // namespace plb

#endif
