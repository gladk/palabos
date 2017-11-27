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

#ifndef DEFAULT_MULTI_GRID_POLICY_2D_H
#define DEFAULT_MULTI_GRID_POLICY_2D_H

#include "multiBlock/defaultMultiBlockPolicy2D.h"
#include "multiGrid/multiGrid2D.h"
#include "multiGrid/multiGridManagement2D.h"

#include <cmath>

namespace plb {

class DefaultMultiGridPolicy2D {
  public:
    /// return a default communicator for a multiGrid with a number of levels
    template <typename T>
    std::vector< BlockCommunicator2D* > getBlockCommunicator(plint levels){
      PLB_ASSERT(levels>0);
      std::vector< BlockCommunicator2D* > res;
      for (int i = 0; i < levels; ++i){
        res.push_back(defaultMultiBlockPolicy2D().getBlockCommunicator());
      }
      return res;
    }
    
    /// return a default MultiGridManagement2D
    MultiGridManagement2D getMultiGridManagement(plint nx, plint ny, plint numLevels, plint referenceLevel) {
        return MultiGridManagement2D(nx, ny, numLevels, referenceLevel );
    }
    
    /// return a vector of CombinedStatistics
    std::vector<CombinedStatistics*> getCombinedStatistics(plint levels) {
        PLB_ASSERT(levels>0);
        std::vector<CombinedStatistics*> stats(levels);
        for (int i = 0; i < levels; ++i){
            // it will be defaultMultiBlockPolicy2D who will give us serial or parallel
            stats[i] = defaultMultiBlockPolicy2D().getCombinedStatistics();
        }
        return stats;
    }

    
    friend DefaultMultiGridPolicy2D& defaultMultiGridPolicy2D();
    
  private:
    
    DefaultMultiGridPolicy2D()
        :numProcesses(global::mpi().getSize())
    {}
    
  private:
    int               numProcesses;
};

inline DefaultMultiGridPolicy2D& defaultMultiGridPolicy2D(){
  static DefaultMultiGridPolicy2D singleton;
  return singleton;
}


} // namespace plb

#endif  // DEFAULT_MULTI_GRID_POLICY_2D_H

