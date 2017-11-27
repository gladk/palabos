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

/* Main author: Daniel Lagrava
 */

/** \file
 * Dynamics and data processors used to implement 3D grid refinement -- header file.
 */

#ifndef HELPER_FINE_GRID_PROCESSORS_3D_H
#define HELPER_FINE_GRID_PROCESSORS_3D_H

#include "core/globalDefs.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataProcessorWrapper3D.h"
#include "multiBlock/multiBlockLattice3D.h"
#include "multiGrid/gridRefinement.h"
#include "multiGrid/gridRefinementDynamics.h"
#include "finiteDifference/fdStencils1D.h"
#include <vector>

namespace plb {

template<typename T, template<typename U> class Descriptor>
void copyPopulations(std::vector<T>& decomposedValues, Cell<T,Descriptor>& cell);

//////////////////////////////////////////////////////////////////////////////
/////////////////////          YZ plane           ////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/// The edges
template<typename T, template<typename U> class Descriptor>
class CubicInterpolationYZLineYHelper3D : public BoxProcessingFunctional3D_LL<T,Descriptor,T,Descriptor>
{
public:
    CubicInterpolationYZLineYHelper3D( plint delta_,RescaleEngine<T,Descriptor>* rescaleEngine_ );
    
    virtual void process( Box3D coarseDomain,
                          BlockLattice3D<T,Descriptor>& coarseLattice,
                          BlockLattice3D<T,Descriptor>& fineLattice );

    virtual CubicInterpolationYZLineYHelper3D<T,Descriptor>* clone() const;
    
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
        modified[1] = modif::allVariables;
    }
    
private:
    plint delta;
    RescaleEngine<T,Descriptor>* rescaleEngine;
    
    Array<Array<T,2>,4> points; // where the interpolated values will be located
};

template<typename T, template<typename U> class Descriptor>
class CubicInterpolationYZLineZHelper3D : public BoxProcessingFunctional3D_LL<T,Descriptor,T,Descriptor>
{
public:
    CubicInterpolationYZLineZHelper3D( plint delta_, RescaleEngine<T,Descriptor>* rescaleEngine_ );
    
    virtual void process( Box3D coarseDomain,
                          BlockLattice3D<T,Descriptor>& coarseLattice,
                          BlockLattice3D<T,Descriptor>& fineLattice );

    virtual CubicInterpolationYZLineZHelper3D<T,Descriptor>* clone() const;
    
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
        modified[1] = modif::allVariables;
    }
    
private:
    
    
private:
    plint delta;
    RescaleEngine<T,Descriptor>* rescaleEngine;
    
    Array<Array<T,2>,4> points; // where the interpolated values will be located
};

/// The corners
template<typename T, template<typename U> class Descriptor>
class CubicInterpolationYZCornerHelper3D : public BoxProcessingFunctional3D_LL<T,Descriptor,T,Descriptor>
{
public:
    CubicInterpolationYZCornerHelper3D(plint deltaY_, plint deltaZ_, RescaleEngine<T,Descriptor>* rescaleEngine_);
    
    virtual ~CubicInterpolationYZCornerHelper3D(){
        delete rescaleEngine;
    }
    
    virtual void process( Box3D coarseDomain,
                          BlockLattice3D<T,Descriptor>& coarseLattice,
                          BlockLattice3D<T,Descriptor>& fineLattice
                          );

    virtual CubicInterpolationYZCornerHelper3D<T,Descriptor>* clone() const;
    
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
        modified[1] = modif::allVariables;
    }
    
private:
    plint deltaY;
    plint deltaZ;
    RescaleEngine<T,Descriptor>* rescaleEngine;
    
    Array<Array<T,2>,3> points; // where the interpolated values will be located
};

//////////////////////////////////////////////////////////////////////////////
/////////////////////          XY plane           ////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/// The edges
template<typename T, template<typename U> class Descriptor>
class CubicInterpolationXYLineXHelper3D : public BoxProcessingFunctional3D_LL<T,Descriptor,T,Descriptor>
{
public:
    CubicInterpolationXYLineXHelper3D( plint delta_,RescaleEngine<T,Descriptor>* rescaleEngine_ );
    
    virtual void process( Box3D coarseDomain,
                          BlockLattice3D<T,Descriptor>& coarseLattice,
                          BlockLattice3D<T,Descriptor>& fineLattice );

    virtual CubicInterpolationXYLineXHelper3D<T,Descriptor>* clone() const;
    
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
        modified[1] = modif::allVariables;
    }
    
private:
    plint delta;
    RescaleEngine<T,Descriptor>* rescaleEngine;
    
    Array<Array<T,2>,4> points; // where the interpolated values will be located
};

template<typename T, template<typename U> class Descriptor>
class CubicInterpolationXYLineYHelper3D : public BoxProcessingFunctional3D_LL<T,Descriptor,T,Descriptor>
{
public:
    CubicInterpolationXYLineYHelper3D( plint delta_, RescaleEngine<T,Descriptor>* rescaleEngine_ );
    
    virtual void process( Box3D coarseDomain,
                          BlockLattice3D<T,Descriptor>& coarseLattice,
                          BlockLattice3D<T,Descriptor>& fineLattice );

    virtual CubicInterpolationXYLineYHelper3D<T,Descriptor>* clone() const;
    
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
        modified[1] = modif::allVariables;
    }
    
private:
    
    
private:
    plint delta;
    RescaleEngine<T,Descriptor>* rescaleEngine;
    
    Array<Array<T,2>,4> points; // where the interpolated values will be located
};



/// The corners
template<typename T, template<typename U> class Descriptor>
class CubicInterpolationXYCornerHelper3D : public BoxProcessingFunctional3D_LL<T,Descriptor,T,Descriptor>
{
public:
    CubicInterpolationXYCornerHelper3D(plint deltaX_, plint deltaY_, RescaleEngine<T,Descriptor>* rescaleEngine_);
    
    virtual ~CubicInterpolationXYCornerHelper3D(){
        delete rescaleEngine;
    }
    
    virtual void process( Box3D coarseDomain,
                          BlockLattice3D<T,Descriptor>& coarseLattice,
                          BlockLattice3D<T,Descriptor>& fineLattice
                          );

    virtual CubicInterpolationXYCornerHelper3D<T,Descriptor>* clone() const;
    
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
        modified[1] = modif::allVariables;
    }
    
private:
    plint deltaX;
    plint deltaY;
    RescaleEngine<T,Descriptor>* rescaleEngine;
    
    Array<Array<T,2>,3> points; // where the interpolated values will be located
};


//////////////////////////////////////////////////////////////////////////////
/////////////////////          XZ plane           ////////////////////////////
//////////////////////////////////////////////////////////////////////////////

/// The edges
template<typename T, template<typename U> class Descriptor>
class CubicInterpolationXZLineXHelper3D : public BoxProcessingFunctional3D_LL<T,Descriptor,T,Descriptor>
{
public:
    CubicInterpolationXZLineXHelper3D( plint delta_,RescaleEngine<T,Descriptor>* rescaleEngine_ );
    
    virtual void process( Box3D coarseDomain,
                          BlockLattice3D<T,Descriptor>& coarseLattice,
                          BlockLattice3D<T,Descriptor>& fineLattice );

    virtual CubicInterpolationXZLineXHelper3D<T,Descriptor>* clone() const;
    
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
        modified[1] = modif::allVariables;
    }
    
private:
    plint delta;
    RescaleEngine<T,Descriptor>* rescaleEngine;
    
    Array<Array<T,2>,4> points; // where the interpolated values will be located
};

template<typename T, template<typename U> class Descriptor>
class CubicInterpolationXZLineZHelper3D : public BoxProcessingFunctional3D_LL<T,Descriptor,T,Descriptor>
{
public:
    CubicInterpolationXZLineZHelper3D( plint delta_, RescaleEngine<T,Descriptor>* rescaleEngine_ );
    
    virtual void process( Box3D coarseDomain,
                          BlockLattice3D<T,Descriptor>& coarseLattice,
                          BlockLattice3D<T,Descriptor>& fineLattice );

    virtual CubicInterpolationXZLineZHelper3D<T,Descriptor>* clone() const;
    
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
        modified[1] = modif::allVariables;
    }
    
private:
    
    
private:
    plint delta;
    RescaleEngine<T,Descriptor>* rescaleEngine;
    
    Array<Array<T,2>,4> points; // where the interpolated values will be located
};



/// The corners
template<typename T, template<typename U> class Descriptor>
class CubicInterpolationXZCornerHelper3D : public BoxProcessingFunctional3D_LL<T,Descriptor,T,Descriptor>
{
public:
    CubicInterpolationXZCornerHelper3D(plint deltaX_, plint deltaZ_, RescaleEngine<T,Descriptor>* rescaleEngine_);
    
    virtual ~CubicInterpolationXZCornerHelper3D(){
        delete rescaleEngine;
    }
    
    virtual void process( Box3D coarseDomain,
                          BlockLattice3D<T,Descriptor>& coarseLattice,
                          BlockLattice3D<T,Descriptor>& fineLattice
                          );

    virtual CubicInterpolationXZCornerHelper3D<T,Descriptor>* clone() const;
    
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        modified[0] = modif::nothing;
        modified[1] = modif::allVariables;
    }
    
private:
    plint deltaX;
    plint deltaZ;
    RescaleEngine<T,Descriptor>* rescaleEngine;
    
    Array<Array<T,2>,3> points; // where the interpolated values will be located
};


}  // namespace plb

#endif  // HELPER_FINE_GRID_PROCESSORS_3D_H
