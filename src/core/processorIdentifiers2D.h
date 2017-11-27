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


#ifndef PROCESSOR_IDENTIFIERS_2D_H
#define PROCESSOR_IDENTIFIERS_2D_H

#include "core/globalDefs.h"
#include "core/plbTypenames.h"
#include "core/util.h"
#include "atomicBlock/dataProcessingFunctional2D.h"
#include <vector>
#include <string>
#include <map>

namespace plb {

namespace meta {

struct ProcessorFactory2D {
    virtual ~ProcessorFactory2D() { }
    virtual BoxProcessingFunctional2D* create(std::string data) const =0;
};

class ProcessorRegistration2D {
public:
    struct Entry {
        Entry(std::string name_, ProcessorFactory2D* factory_)
            : name(name_),
              factory(factory_)
        { }
        std::string name;
        ProcessorFactory2D* factory;
    };
    struct EntryLessThan {
        bool operator()(Entry const& entry1, Entry const& entry2) const {
            return entry1.name < entry2.name;
        }
    };
    typedef std::map<Entry,int,EntryLessThan> EntryMap;
public:
    ProcessorRegistration2D() { }
    ~ProcessorRegistration2D();
    int announce(std::string nameOfProcessor, ProcessorFactory2D* factory_=0);
    int getId(std::string name) const;
    int getNumId() const;
    std::string getName(int id) const;
    BoxProcessingFunctional2D* create(std::string procName, std::string data);
    EntryMap::const_iterator begin() const;
    EntryMap::const_iterator end() const;
private:
    ProcessorRegistration2D(ProcessorRegistration2D const& rhs) { }
    ProcessorRegistration2D& operator= (
            ProcessorRegistration2D const& rhs )
    { return *this; }
private:
    EntryMap processorByName;
    std::vector<Entry> processorByNumber;
};

ProcessorRegistration2D& processorRegistration2D();

template<class FunctionalT>
class FunctionalFactory2D : public ProcessorFactory2D
{
    virtual BoxProcessingFunctional2D* create(std::string data) const {
        BoxProcessingFunctional2D* functional = new FunctionalT;
        functional->unserialize(data);
        return functional;
    }
};

template<class FunctionalT>
int registerProcessor2D(std::string name) {
    return processorRegistration2D().announce(name, new FunctionalFactory2D<FunctionalT>);
}

template <
    class FunctionalT, typename T >
int registerProcessor2D(std::string name) {
    std::string fullName =
        name + "_" + std::string(NativeType<T>::getName());
    return processorRegistration2D().announce(fullName, new FunctionalFactory2D<FunctionalT>);
}

template <
    class FunctionalT,
    typename T1, typename T2 >
int registerProcessor2D(std::string name) {
    std::string fullName =
        name + "_" + std::string(NativeType<T1>::getName())+ "_" + std::string(NativeType<T2>::getName());
    return processorRegistration2D().announce(fullName, new FunctionalFactory2D<FunctionalT>);
}

template <
    class FunctionalT,
    typename T, int val >
int registerProcessor2D(std::string name) {
    std::string fullName =
        name + "_" + std::string(NativeType<T>::getName()) + "_" + util::val2str(val);
    return processorRegistration2D().announce(fullName, new FunctionalFactory2D<FunctionalT>);
}

template <
    class FunctionalT, typename T, template<typename U> class Descriptor >
int registerProcessor2D(std::string name) {
    std::string fullName =
        name + "_" + std::string(NativeType<T>::getName()) + "_" + std::string(Descriptor<T>::name);
    return processorRegistration2D().announce(fullName, new FunctionalFactory2D<FunctionalT>);
}

template <
    class FunctionalT,
    typename T, template<typename U> class Descriptor, int val1, int val2 >
int registerProcessor2D(std::string name) {
    std::string fullName =
        name + "_" + std::string(NativeType<T>::getName()) + "_" +
        std::string(Descriptor<T>::name) + "_" + util::val2str(val1) + "_" + util::val2str(val2);
    return processorRegistration2D().announce(fullName, new FunctionalFactory2D<FunctionalT>);
}

template <
    class FunctionalT,
    typename T, template<typename U> class Descriptor, int val1, int val2, int val3 >
int registerProcessor2D(std::string name) {
    std::string fullName =
        name + "_" + std::string(NativeType<T>::getName()) + "_" +
        std::string(Descriptor<T>::name) + "_" + util::val2str(val1) + "_" + util::val2str(val2) + "_" + util::val2str(val3);
    return processorRegistration2D().announce(fullName, new FunctionalFactory2D<FunctionalT>);
}

}  // namespace meta

}  // namespace plb

#endif  // PROCESSOR_IDENTIFIERS_2D_H
