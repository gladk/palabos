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


#ifndef MULTI_BLOCK_IDENTIFIERS_3D_H
#define MULTI_BLOCK_IDENTIFIERS_3D_H

#include "core/globalDefs.h"
#include "multiBlock/multiBlock3D.h"
#include <string>
#include <map>

namespace plb {

namespace meta {

struct MultiBlockGenerator3D {
    virtual ~MultiBlockGenerator3D() { }
    virtual MultiBlock3D* generate(MultiBlockManagement3D const& manager, plint nDim) const =0;
};

template<class MultiBlockGenT>
class SpecificMultiBlockGenerator3D : public MultiBlockGenerator3D
{
public:
    SpecificMultiBlockGenerator3D(MultiBlockGenT multiBlockGen_)
        : multiBlockGen(multiBlockGen_)
    { }
    virtual MultiBlock3D* generate(MultiBlockManagement3D const& manager, plint nDim) const {
        return multiBlockGen(manager, nDim).release();
    }
private:
    MultiBlockGenT multiBlockGen;
};

class MultiBlockRegistration3D {
public:
    typedef std::map<std::string, MultiBlockGenerator3D*> Str_gen_map;
    typedef std::map<std::string, Str_gen_map> Str_str_gen_map;
    typedef std::map<std::string, Str_str_gen_map> Str_str_str_gen_map;
public:
    ~MultiBlockRegistration3D();
    int announce( std::string T_name, std::string Descriptor_name,
                  std::string nameOfMultiBlock, MultiBlockGenerator3D* generator);
    MultiBlock3D* generate (
                      std::string T_name, std::string Descriptor_name,
                      std::string nameOfMultiBlock,
                      MultiBlockManagement3D const& manager, plint nDim=0 );
public:
    /// This default constructor should actually be private, but it is public
    ///  for now to fix a parse error in older GCCs.
    MultiBlockRegistration3D()
        : numRegistered(0)
    { }
private:
    MultiBlockRegistration3D(MultiBlockRegistration3D const& rhs)
        : numRegistered(rhs.numRegistered)
    { }
    MultiBlockRegistration3D& operator=(MultiBlockRegistration3D const& rhs)
    {
        numRegistered = rhs.numRegistered;
        return *this;
    }
private:
    int numRegistered;
    Str_str_str_gen_map generators;
};

MultiBlockRegistration3D& multiBlockRegistration3D();

template<class MultiBlockGenT>
int registerMultiBlock3D (
             std::string T_name, std::string Descriptor_name, std::string nameOfMultiBlock,
             MultiBlockGenT multiBlockGen )
{
    return multiBlockRegistration3D().announce (
               T_name, Descriptor_name, nameOfMultiBlock,
               new SpecificMultiBlockGenerator3D<MultiBlockGenT>(multiBlockGen) );
}

}  // namespace meta

}  // namespace plb

#endif  // MULTI_BLOCK_IDENTIFIERS_3D_H
