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

#ifndef ADVECTION_DIFFUSION_PROCESSOR_3D_H
#define ADVECTION_DIFFUSION_PROCESSOR_3D_H

#include "core/globalDefs.h"
#include "core/dynamics.h"
#include "atomicBlock/dataProcessingFunctional3D.h"


namespace plb {

// This data processor uses symmetric finite differences to compute a gradient.
// It cannot be applied on any part of the boundary of the global simulation
// domain, except if this boundary is periodic.
template<typename T, template<typename U> class Descriptor>
class SetEffectiveDiffusivity3D : public BoxProcessingFunctional3D_LS<T,Descriptor,T>
{
public:
    SetEffectiveDiffusivity3D(T omega0_, T T0_, T cSmago_)
        : omega0(omega0_),
          invT0((T) 1 / T0_),
          cSmagoSqr(cSmago_*cSmago_)
    { }
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice,
                                       ScalarField3D<T>& rhoBar);
    virtual SetEffectiveDiffusivity3D<T,Descriptor>* clone() const
    {
        return new SetEffectiveDiffusivity3D<T,Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const
    {
        modified[0] = modif::dynamicVariables;   // lattice
        modified[1] = modif::nothing;            // rhoBar
    }
private:
    T omega0;
    T invT0;
    T cSmagoSqr;
};
    
} // namespace plb

#endif  // ADVECTION_DIFFUSION_PROCESSOR_3D_H

