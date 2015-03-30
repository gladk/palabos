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
 * Sponge (absorbing) zones, to be mainly used in addition to outflow boundary conditions -- header file.
 */

#ifndef SPONGE_ZONES_3D_H
#define SPONGE_ZONES_3D_H

#include "core/globalDefs.h"
#include "core/dynamics.h"
#include "atomicBlock/dataProcessingFunctional3D.h"

namespace plb {

// Data processor to implement a viscosity sponge zone:
// The relaxation parameter (kinematic viscosity) is progressively increased
// inside the zone.
// The dynamics object of every cell is changed by this data processor, so
// the user must make sure in the construction of the MultiBlockLattice3D that
// each node has its own dynamics object.
template<typename T, template<typename U> class Descriptor>
class ViscositySpongeZone : public BoxProcessingFunctional3D
{
public:
    // Constructor for the tanh sponge function.
    //   Nice value for the translation parameters is 0.5.
    //   Nice value for the scale parameters is 0.12.
    ViscositySpongeZone(plint nx_, plint ny_, plint nz_, T bulkOmega_,
            Array<plint,6> const& numSpongeCells_, Array<T,6> const& translationParameters_,
            Array<T,6> const& scaleParameters_);

    // Constructor for the cos sponge function.
    ViscositySpongeZone(plint nx_, plint ny_, plint nz_, T bulkOmega_,
                        Array<plint,6> const& numSpongeCells_);

    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks);

    virtual ViscositySpongeZone<T,Descriptor>* clone() const
    {
        return new ViscositySpongeZone<T,Descriptor>(*this);
    }

    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const
    {
        modified[0] = modif::dynamicVariables; // Block lattice.
    }
private:
    plint nx, ny, nz;                  // Lattice dimensions.
    T bulkOmega;                       // Value of the relaxation parameter outside of the sponge zone.
    Array<plint,6> numSpongeCells;     // Width of the sponge zones.
    Array<T,6> translationParameters;  // Translation parameters of the tanh sponge functions.
    Array<T,6> scaleParameters;        // Scaling parameters of the tanh sponge functions.
    bool useTanhSpongeFunction;        // Use a tanh sponge function, or a cos sponge function.
};


template<typename T, template<typename U> class Descriptor>
class MaskedViscositySpongeZone : public BoxProcessingFunctional3D
{
public:
    // Constructor for the tanh sponge function.
    //   Nice value for the translation parameters is 0.5.
    //   Nice value for the scale parameters is 0.12.
    MaskedViscositySpongeZone(plint nx_, plint ny_, plint nz_, T bulkOmega_, int flag_,
            Array<plint,6> const& numSpongeCells_, Array<T,6> const& translationParameters_,
            Array<T,6> const& scaleParameters_);

    // Constructor for the cos sponge function.
    MaskedViscositySpongeZone(plint nx_, plint ny_, plint nz_, T bulkOmega_, int flag_,
                        Array<plint,6> const& numSpongeCells_);

    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks);

    virtual MaskedViscositySpongeZone<T,Descriptor>* clone() const
    {
        return new MaskedViscositySpongeZone<T,Descriptor>(*this);
    }

    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const
    {
        modified[0] = modif::dynamicVariables; // Block lattice.
        modified[1] = modif::nothing;          // Flag matrix.
    }
private:
    plint nx, ny, nz;                  // Lattice dimensions.
    T bulkOmega;                       // Value of the relaxation parameter outside of the sponge zone.
    Array<plint,6> numSpongeCells;     // Width of the sponge zones.
    Array<T,6> translationParameters;  // Translation parameters of the tanh sponge functions.
    Array<T,6> scaleParameters;        // Scaling parameters of the tanh sponge functions.
    bool useTanhSpongeFunction;        // Use a tanh sponge function, or a cos sponge function.
    int flag;
};


// Data processor to implement a Smagorinsky sponge zone:
// The Smagorinsky parameter is progressively increased
// inside the zone.
// The dynamics object of every cell is changed by this data processor, so
// the user must make sure in the construction of the MultiBlockLattice3D that
// each node has its own dynamics object.
template<typename T, template<typename U> class Descriptor>
class SmagorinskySpongeZone : public BoxProcessingFunctional3D
{
public:
    // Constructor for the tanh sponge function.
    //   Nice value for the translation parameters is 0.5.
    //   Nice value for the scale parameters is 0.12.
    SmagorinskySpongeZone(plint nx_, plint ny_, plint nz_, T bulkCSmago_, T targetCSmago_,
            Array<plint,6> const& numSpongeCells_, Array<T,6> const& translationParameters_,
            Array<T,6> const& scaleParameters_);

    // Constructor for the cos sponge function.
    SmagorinskySpongeZone(plint nx_, plint ny_, plint nz_, T bulkCSmago_, T targetCSmago_,
            Array<plint,6> const& numSpongeCells_);

    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks);

    virtual SmagorinskySpongeZone<T,Descriptor>* clone() const
    {
        return new SmagorinskySpongeZone<T,Descriptor>(*this);
    }

    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const
    {
        modified[0] = modif::dynamicVariables; // Block lattice.
    }
private:
    plint nx, ny, nz;                  // Lattice dimensions.
    T bulkCSmago, targetCSmago;        // Varying parameter: bulk and target values.
    Array<plint,6> numSpongeCells;     // Width of the sponge zones.
    Array<T,6> translationParameters;  // Translation parameters of the tanh sponge functions.
    Array<T,6> scaleParameters;        // Scaling parameters of the tanh sponge functions.
    bool useTanhSpongeFunction;        // Use a tanh sponge function, or a cos sponge function.
};


template<typename T, template<typename U> class Descriptor>
class MaskedSmagorinskySpongeZone : public BoxProcessingFunctional3D
{
public:
    // Constructor for the tanh sponge function.
    //   Nice value for the translation parameters is 0.5.
    //   Nice value for the scale parameters is 0.12.
    MaskedSmagorinskySpongeZone(plint nx_, plint ny_, plint nz_, T bulkCSmago_, T targetCSmago_, int flag_,
            Array<plint,6> const& numSpongeCells_, Array<T,6> const& translationParameters_,
            Array<T,6> const& scaleParameters_);

    // Constructor for the cos sponge function.
    MaskedSmagorinskySpongeZone(plint nx_, plint ny_, plint nz_, T bulkCSmago_, T targetCSmago_, int flag_,
            Array<plint,6> const& numSpongeCells_);

    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> blocks);

    virtual MaskedSmagorinskySpongeZone<T,Descriptor>* clone() const
    {
        return new MaskedSmagorinskySpongeZone<T,Descriptor>(*this);
    }

    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const
    {
        modified[0] = modif::dynamicVariables; // Block lattice.
        modified[1] = modif::nothing;          // Flag matrix.
    }
private:
    plint nx, ny, nz;                  // Lattice dimensions.
    T bulkCSmago, targetCSmago;        // Varying parameter: bulk and target values.
    Array<plint,6> numSpongeCells;     // Width of the sponge zones.
    Array<T,6> translationParameters;  // Translation parameters of the tanh sponge functions.
    Array<T,6> scaleParameters;        // Scaling parameters of the tanh sponge functions.
    bool useTanhSpongeFunction;        // Use a tanh sponge function, or a cos sponge function.
    int flag;
};

}

#endif  // SPONGE_ZONES_3D_H

