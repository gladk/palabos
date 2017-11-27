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


#include "core/multiBlockIdentifiers2D.h"
#include "core/util.h"

namespace plb {

namespace meta {

MultiBlockRegistration2D::~MultiBlockRegistration2D() {
    Str_str_str_gen_map::iterator it1 = generators.begin();
    for (; it1 != generators.end(); ++it1) {
        Str_str_gen_map::iterator it2 = it1->second.begin();
        for (; it2 != it1->second.end(); ++it2) {
            Str_gen_map::iterator it3 = it2->second.begin();
            for (; it3 != it2->second.end(); ++it3) {
                delete it3->second;
            }
        }
    }
}

int MultiBlockRegistration2D::announce (
        std::string T_name, std::string Descriptor_name,
        std::string nameOfMultiBlock, MultiBlockGenerator2D* generator )
{
    if (Descriptor_name=="")  {
        Descriptor_name="NA";
    }
    Str_gen_map::iterator it = 
        generators[T_name][Descriptor_name].find(nameOfMultiBlock);
    if (it != generators[T_name][Descriptor_name].end()) {
        delete it->second;
    }
    generators[T_name][Descriptor_name][nameOfMultiBlock] = generator;
    return numRegistered++;
}

MultiBlock2D* MultiBlockRegistration2D::generate (
                  std::string T_name, std::string Descriptor_name,
                  std::string nameOfMultiBlock, MultiBlockManagement2D const& manager, plint nDim )
{
    if (Descriptor_name=="") {
        Descriptor_name="NA";
    }
    Str_str_str_gen_map::iterator it1 = generators.find(T_name);
    if (it1 == generators.end()) {
        return 0;
    }
    Str_str_gen_map::iterator it2 = it1->second.find(Descriptor_name);
    if (it2 == it1->second.end()) {
        return 0;
    }
    Str_gen_map::iterator it3 = it2->second.find(nameOfMultiBlock);
    if (it3 == it2->second.end()) {
        return 0;
    }
    else {
        return it3->second->generate(manager, nDim);
    }
}

MultiBlockRegistration2D& multiBlockRegistration2D() {
    static MultiBlockRegistration2D instance;
    return instance;
}

}  // namespace meta

}  // namespace plb
