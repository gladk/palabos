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


#ifndef PARTICLE_IDENTIFIERS_2D_H
#define PARTICLE_IDENTIFIERS_2D_H

#include "core/globalDefs.h"
#include "core/hierarchicSerializer.h"
#include <vector>
#include <string>
#include <map>

namespace plb {

namespace meta {

template<typename T, template<typename U> class Descriptor>
struct ParticleGenerator2D {
    virtual ~ParticleGenerator2D() { }
    virtual Particle2D<T,Descriptor>*
        generate(HierarchicUnserializer& unserializer) const =0;
};

template<typename T, template<typename U> class Descriptor>
class ParticleRegistration2D {
public:
    struct Entry {
        Entry(std::string name_, ParticleGenerator2D<T,Descriptor>* generator_)
            : name(name_),
              generator(generator_)
        { }
        std::string name;
        ParticleGenerator2D<T,Descriptor>* generator;
    };
    struct EntryLessThan {
        bool operator()(Entry const& entry1, Entry const& entry2) const {
            return entry1.name < entry2.name;
        }
    };
    typedef std::map<Entry,int,EntryLessThan> EntryMap;
public:
    ~ParticleRegistration2D();
    int announce(std::string nameOfParticle,
                 ParticleGenerator2D<T,Descriptor>* generator_=0);
    int getId(std::string name) const;
    int getNumId() const;
    std::string getName(int id) const;
    Particle2D<T,Descriptor>* generate(HierarchicUnserializer& unserializer);
    typename EntryMap::const_iterator begin() const;
    typename EntryMap::const_iterator end() const;
public:
    /// This default constructor should actually be private, but it is public
    ///  for now to fix a parse error in older GCCs.
    ParticleRegistration2D() { }
private:
    ParticleRegistration2D(ParticleRegistration2D<T,Descriptor> const& rhs) { }
    ParticleRegistration2D<T,Descriptor>& operator= (
            ParticleRegistration2D<T,Descriptor> const& rhs )
    {
        return *this;
    }
private:
    EntryMap particleByName;
    std::vector<Entry> particleByNumber;

// TODO: This friend declaration is not properly parsed in older (but not-so-old) GCC
//   compilers. Therefore, it is commented for now, and the default constructor of
//   ParticleRegistration is public, although it should be private, because
//   ParticleRegistration is a singleton.
//
// template<typename T_, template<typename U_> class Descriptor_>
// friend ParticleRegistration2D<T_,Descriptor_>& particleRegistration();
};

template<typename T, template<typename U> class Descriptor>
ParticleRegistration2D<T,Descriptor>& particleRegistration2D();


template< typename T,
          template<typename U> class Descriptor,
          class ParticleT >
class GenericParticleGenerator2D : public ParticleGenerator2D<T,Descriptor>
{
    virtual Particle2D<T,Descriptor>* generate (
            HierarchicUnserializer& unserializer ) const
    {
        Particle2D<T,Descriptor>* particle = new ParticleT();
        particle->unserialize(unserializer);
        return particle;
    }
};



template< typename T,
          template<typename U> class Descriptor,
          class PointParticle >
class PointParticleGenerator2D : public ParticleGenerator2D<T,Descriptor>
{
    virtual Particle2D<T,Descriptor>* generate (
            HierarchicUnserializer& unserializer ) const
    {
        plint tag;
        unserializer.readValue(tag);
        Array<T,2> position;
        unserializer.readValues<T,2>(position);
        Array<T,2> velocity;
        unserializer.readValues<T,2>(velocity);
        return new PointParticle(tag, position, velocity);
    }
};


template< typename T,
          template<typename U> class Descriptor,
          class ParticleT >
int registerGenericParticle2D(std::string name) {
    return particleRegistration2D<T,Descriptor>().announce (
               name, new GenericParticleGenerator2D<T,Descriptor,ParticleT> );
}


template< typename T,
          template<typename U> class Descriptor,
          class PointParticle >
int registerPointParticle2D(std::string name) {
    return particleRegistration2D<T,Descriptor>().announce (
               name, new PointParticleGenerator2D<T,Descriptor,PointParticle> );
}

}  // namespace meta

}  // namespace plb

#endif  // PARTICLE_IDENTIFIERS_2D_H
