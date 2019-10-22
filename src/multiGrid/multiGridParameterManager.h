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
 **/

#ifndef MULTI_GRID_PARAMETER_MANAGER_H
#define MULTI_GRID_PARAMETER_MANAGER_H

#include "core/globalDefs.h"
#include "core/units.h"

namespace plb {

/// Interface for a wrapper of refinement parameters (each level posseses specific parameters)
template <typename T>
class RefinementParameters {
  public:
    RefinementParameters( plint levelNumber_, plint referenceLevel_, IncomprFlowParam<T> parameters_);
    RefinementParameters(RefinementParameters<T> const& rhs);
    
    /// This function shall define the scaling between the grids by defining all the parameters
    virtual void createParameters()=0;
    virtual ~RefinementParameters(){}
    
    IncomprFlowParam<T> const& getParameters(plint lvl) const;
      
    plint getReferenceLevel();
    plint getNumLevels();
    
    IncomprFlowParam<T> const& operator[](plint lvl) const;
    
    void putParameter(IncomprFlowParam<T> &newParam);
    
    IncomprFlowParam<T>& getOriginalParameters();
    
  private:
    plint levelNumber;
    plint referenceLevel;
    IncomprFlowParam<T> originalParameters;
    std::vector< IncomprFlowParam<T> > parameters;
    
};

/// In the convective scaling (dx=dt) 
template <typename T>
class ConvectiveRefinementParameters : public RefinementParameters<T> {
  
  public:
    ConvectiveRefinementParameters( plint levels, plint reference, IncomprFlowParam<T> params);
    
    virtual void createParameters();    
};

/// In the convective scaling (dx=dt^2) so viscosity is constant
template <typename T>
class DiffusiveRefinementParameters : public RefinementParameters<T> {
  
  public:
    DiffusiveRefinementParameters(plint levels, plint reference, IncomprFlowParam<T> params);
   
    virtual void createParameters();
};


} //namespace plb

#endif  // MULTI_GRID_PARAMETER_MANAGER_H

