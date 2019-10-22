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


#ifndef DYNAMICS_IDENTIFIERS_HH
#define DYNAMICS_IDENTIFIERS_HH

#include "core/dynamicsIdentifiers.h"
#include "core/runTimeDiagnostics.h"
#include <sstream>

namespace plb {

namespace meta {

template<typename T, template<typename U> class Descriptor>
void createIdIndirection (
        std::map<int,std::string> const& foreignIdToName,
        std::map<int,int>& idIndirect )
{
    idIndirect.clear();
    std::map<int,std::string>::const_iterator it = foreignIdToName.begin();
    for (; it != foreignIdToName.end(); ++it) {
        int foreignId = it->first;
        std::string dynamicsName = it->second;
        int id = dynamicsRegistration<T,Descriptor>().getId(dynamicsName);
        if (id==0) {
            // You are trying to unserialize a dynamics object which has not been
            // compiled into your current object code. Try to explicitly instantiate
            // a dynamics object for all kinds which are going to be unserialized.
            PLB_ASSERT( false );
        }
        idIndirect.insert(std::pair<int,int>(foreignId,id));
    }
}

template<typename T, template<typename U> class Descriptor>
DynamicsRegistration<T,Descriptor>::~DynamicsRegistration()
{
    for (pluint iEntry=0; iEntry<dynamicsByNumber.size(); ++iEntry) {
        delete dynamicsByNumber[iEntry].generator;
    }
}

template<typename T, template<typename U> class Descriptor>
int DynamicsRegistration<T,Descriptor>::announce (
        std::string nameOfDynamics,
        DynamicsGenerator<T,Descriptor>* generator )
{
    Entry entry(nameOfDynamics, generator);
    typename EntryMap::iterator it = dynamicsByName.find(entry);
    if (it != dynamicsByName.end()) {
        plbLogicError( std::string("The dynamics class ") + nameOfDynamics +
                       std::string(" was registered twice") );
    }
    dynamicsByNumber.push_back(entry);
    int nextId = dynamicsByNumber.size();
    dynamicsByName[entry] = nextId;
    return nextId;
}

template<typename T, template<typename U> class Descriptor>
int DynamicsRegistration<T,Descriptor>::getId(std::string name) const
{
    Entry entry(name, 0);
    typename EntryMap::const_iterator it = dynamicsByName.find(entry);
    if (it == dynamicsByName.end()) {
        return 0;
    }
    else {
        return it->second;
    }
}

template<typename T, template<typename U> class Descriptor>
int DynamicsRegistration<T,Descriptor>::getNumId() const
{
    return (int)(dynamicsByNumber.size());
}

template<typename T, template<typename U> class Descriptor>
std::string DynamicsRegistration<T,Descriptor>::getName(int id) const
{
    if (id==0) {
        return std::string("Undefined");
    }
    if (id < 0 || id > (int)dynamicsByNumber.size()) {
        std::stringstream message;
        message << "A dynamics class with ID " << id << " doesn't exist.";
        plbLogicError(message.str());
    }
    return dynamicsByNumber[id-1].name;
}

template<typename T, template<typename U> class Descriptor>
Dynamics<T,Descriptor>* DynamicsRegistration<T,Descriptor>::generate
                             ( HierarchicUnserializer& unserializer )
{
    plint id = unserializer.getId();
    PLB_ASSERT (id>0 && (pluint)id <= dynamicsByNumber.size());
    return dynamicsByNumber[id-1].generator->generate(unserializer);
}

template<typename T, template<typename U> class Descriptor>
typename DynamicsRegistration<T,Descriptor>::EntryMap::const_iterator
    DynamicsRegistration<T,Descriptor>::begin() const
{
    return dynamicsByName.begin();
}

template<typename T, template<typename U> class Descriptor>
typename DynamicsRegistration<T,Descriptor>::EntryMap::const_iterator
    DynamicsRegistration<T,Descriptor>::end() const
{
    return dynamicsByName.end();
}


template<typename T, template<typename U> class Descriptor>
DynamicsRegistration<T,Descriptor>& dynamicsRegistration() {
    static DynamicsRegistration<T,Descriptor> instance;
    return instance;
}

template<typename T, template<typename U> class Descriptor>
std::string constructIdNameChain(std::vector<int> const& ids, std::string separator)
{
    std::string nameChain;
    if (ids.size()>0 && ids[0]>=0) {
        nameChain += dynamicsRegistration<T,Descriptor>().
                         getName(ids[0]);
    }
    for (pluint i=1; i<ids.size(); ++i) {
        if (ids[i]>=0) {
            nameChain += separator;
            nameChain += dynamicsRegistration<T,Descriptor>().
                             getName(ids[i]);
        }
    }
    return nameChain;
}

}  // namespace meta

}  // namespace plb

#endif  // DYNAMICS_IDENTIFIERS_HH
