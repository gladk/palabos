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

/** \file
 * Coupling between grids of different refinement level -- header file.
 */

#ifndef MULTI_LEVEL_NTENSOR_FIELD_3D_H
#define MULTI_LEVEL_NTENSOR_FIELD_3D_H

#include "core/globalDefs.h"

#include "gridRefinement/octreeGridStructure.h"
#include "gridRefinement/multiLevel3D.h"
#include "gridRefinement/boxLogic3D.h"

#include <map>
#include <vector>

namespace plb {

template<typename T>
class MultiLevelNTensorField3D : public MultiLevel3D {
public:
    // The dynamics dyn must be right for the coarsest instantiated level.
    MultiLevelNTensorField3D(OctreeGridStructure& ogs_);
    MultiLevelNTensorField3D(OctreeGridStructure& ogs_, Box3D domain, plint domainLevel);
    MultiLevelNTensorField3D(MultiLevelNTensorField3D<T> const& rhs);
    MultiLevelNTensorField3D<T>& operator=(MultiLevelNTensorField3D<T> const& rhs);
    void swap(MultiLevelNTensorField3D<T>& rhs);
    ~MultiLevelNTensorField3D();

    virtual MultiNTensorField3D<T> const &getLevel(plint iL) const {
        PLB_ASSERT(iL <= getNumLevels() && iL >= 0);
        return *gridLevels[iL];
    }
    
    virtual MultiNTensorField3D<T> &getLevel(plint iL) {
        PLB_ASSERT(iL <= getNumLevels() && iL >= 0);
        return *gridLevels[iL];
    }
    
    plint const getNumLevels() const { 
        return ogs.getNumLevels();
    }
    
    OctreeGridStructure const& getOgs() const {
        return ogs;
    }
private:
    OctreeGridStructure ogs;
private:
    std::vector<MultiNTensorField3D<T> *> gridLevels;
};

template<typename T>
class MultiLevelNTensorFieldForOutput3D : public MultiLevel3D {
public:
    // The dynamics dyn must be right for the coarsest instantiated level.
    MultiLevelNTensorFieldForOutput3D(OctreeGridStructure& ogs_, bool crop);
    MultiLevelNTensorFieldForOutput3D(OctreeGridStructure& ogs_, Box3D domain, plint domainLevel, bool crop);
    MultiLevelNTensorFieldForOutput3D(MultiLevelNTensorFieldForOutput3D<T> const& rhs);
    MultiLevelNTensorFieldForOutput3D<T>& operator=(MultiLevelNTensorFieldForOutput3D<T> const& rhs);
    void swap(MultiLevelNTensorFieldForOutput3D<T>& rhs);
    ~MultiLevelNTensorFieldForOutput3D();

    virtual MultiNTensorField3D<T> const &geLevel(plint iL) const {
        PLB_ASSERT(iL <= getNumLevels() && iL >= 0);
        return *gridLevels[iL];
    }
    
    virtual MultiNTensorField3D<T> &getLevel(plint iL) {
        PLB_ASSERT(iL <= getNumLevels() && iL >= 0);
        return *gridLevels[iL];
    }
    
    plint const getNumLevels() const { 
        return ogs.getNumLevels();
    }
    
    OctreeGridStructure const& getOgs() const {
        return ogs;
    }
private:
    OctreeGridStructure ogs;
private:
    std::vector<MultiNTensorField3D<T> *> gridLevels;
};

}  // namespace plb

#endif  // MULTI_LEVEL_NTENSOR_FIELD_3D_H
