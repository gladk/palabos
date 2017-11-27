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


#ifndef DYNAMICS_IDENTIFIERS_H
#define DYNAMICS_IDENTIFIERS_H

#include "core/globalDefs.h"
#include "core/dynamics.h"
#include "core/hierarchicSerializer.h"
#include <vector>
#include <string>
#include <map>

namespace plb {

namespace meta {

template<typename T, template<typename U> class Descriptor>
void createIdIndirection (
        std::map<int,std::string> const& foreignIdToName,
        std::map<int,int>& idIndirect );

template<typename T, template<typename U> class Descriptor>
struct DynamicsGenerator {
    virtual ~DynamicsGenerator() { }
    virtual Dynamics<T,Descriptor>*
        generate(HierarchicUnserializer& unserializer) const =0;
};

template<typename T, template<typename U> class Descriptor>
class DynamicsRegistration {
public:
    struct Entry {
        Entry(std::string name_, DynamicsGenerator<T,Descriptor>* generator_)
            : name(name_),
              generator(generator_)
        { }
        std::string name;
        DynamicsGenerator<T,Descriptor>* generator;
    };
    struct EntryLessThan {
        bool operator()(Entry const& entry1, Entry const& entry2) const {
            return entry1.name < entry2.name;
        }
    };
    typedef std::map<Entry,int,EntryLessThan> EntryMap;
public:
    ~DynamicsRegistration();
    int announce(std::string nameOfDynamics,
                 DynamicsGenerator<T,Descriptor>* generator_=0);
    int getId(std::string name) const;
    int getNumId() const;
    std::string getName(int id) const;
    Dynamics<T,Descriptor>* generate(HierarchicUnserializer& unserializer);
    typename EntryMap::const_iterator begin() const;
    typename EntryMap::const_iterator end() const;
public:
    /// This default constructor should actually be private, but it is public
    ///  for now to fix a parse error in older GCCs.
    DynamicsRegistration() { }
private:
    DynamicsRegistration(DynamicsRegistration<T,Descriptor> const& rhs) { }
    DynamicsRegistration<T,Descriptor>& operator= (
            DynamicsRegistration<T,Descriptor> const& rhs )
    {
        return *this;
    }
private:
    EntryMap dynamicsByName;
    std::vector<Entry> dynamicsByNumber;

// TODO: This friend declaration is not properly parsed in older (but not-so-old) GCC
//   compilers. Therefore, it is commented for now, and the default constructor of
//   DynamicsRegistration is public, although it should be private, because
//   DynamicsRegistration is a singleton.
//
// template<typename T_, template<typename U_> class Descriptor_>
// friend DynamicsRegistration<T_,Descriptor_>& dynamicsRegistration();
};

template<typename T, template<typename U> class Descriptor>
DynamicsRegistration<T,Descriptor>& dynamicsRegistration();

template<typename T, template<typename U> class Descriptor>
std::string constructIdNameChain(std::vector<int> const& ids, std::string separator=".");


template< typename T,
          template<typename U> class Descriptor,
          class NoParamDynamics >
class NoParamDynamicsGenerator : public DynamicsGenerator<T,Descriptor>
{
    virtual Dynamics<T,Descriptor>* generate(HierarchicUnserializer& unserializer) const {
        return new NoParamDynamics();
    }
};

template< typename T,
          template<typename U> class Descriptor,
          class OneParamDynamics >
class OneParamDynamicsGenerator : public DynamicsGenerator<T,Descriptor>
{
    virtual Dynamics<T,Descriptor>* generate(HierarchicUnserializer& unserializer) const {
        return new OneParamDynamics(unserializer.readValue<T>());
    }
};

template< typename T,
          template<typename U> class Descriptor,
          class TwoParamDynamics >
class TwoParamDynamicsGenerator : public DynamicsGenerator<T,Descriptor>
{
    virtual Dynamics<T,Descriptor>* generate(HierarchicUnserializer& unserializer) const {
        return new TwoParamDynamics (
                       unserializer.readValue<T>(),
                       unserializer.readValue<T>() );
    }
};

template< typename T,
          template<typename U> class Descriptor,
          class GeneralDynamics >
class GeneralDynamicsGenerator : public DynamicsGenerator<T,Descriptor>
{
    virtual Dynamics<T,Descriptor>* generate(HierarchicUnserializer& unserializer) const {
        return new GeneralDynamics(unserializer);
    }
};

template< typename T,
          template<typename U> class Descriptor,
          class CompDynamics >
class CompositeDynamicsGenerator : public DynamicsGenerator<T,Descriptor>
{
    virtual Dynamics<T,Descriptor>* generate(HierarchicUnserializer& unserializer) const {
        bool automaticPrepareCollision;
        unserializer.readValue(automaticPrepareCollision);
        Dynamics<T,Descriptor>* compositeDynamics = dynamicsRegistration<T,Descriptor>().generate(unserializer);
        return new CompDynamics(compositeDynamics, automaticPrepareCollision);
    }
};

template< typename T,
          template<typename U> class Descriptor,
          class GeneralDynamics >
int registerGeneralDynamics(std::string name) {
    return dynamicsRegistration<T,Descriptor>().announce (
               name, new GeneralDynamicsGenerator<T,Descriptor,GeneralDynamics> );
}

template< typename T,
          template<typename U> class Descriptor,
          class CompDynamics >
int registerCompositeDynamics(std::string name) {
    return dynamicsRegistration<T,Descriptor>().announce (
               name, new CompositeDynamicsGenerator<T,Descriptor,CompDynamics> );
}

}  // namespace meta

}  // namespace plb

#endif  // DYNAMICS_IDENTIFIERS_H

