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


#include "core/processorIdentifiers3D.h"
#include "core/runTimeDiagnostics.h"
#include <sstream>

namespace plb {

namespace meta {

ProcessorRegistration3D::~ProcessorRegistration3D()
{
    for (pluint iEntry=0; iEntry<processorByNumber.size(); ++iEntry) {
        delete processorByNumber[iEntry].factory;
    }
}

int ProcessorRegistration3D::announce (
        std::string nameOfProcessor,
        ProcessorFactory3D* factory_ )
{
    Entry entry(nameOfProcessor, factory_);
    EntryMap::iterator it = processorByName.find(entry);
    if (it != processorByName.end()) {
        plbLogicError( std::string("The processor ") + nameOfProcessor +
                       std::string(" was registered twice") );
    }
    processorByNumber.push_back(entry);
    int nextId = processorByNumber.size();
    processorByName[entry] = nextId;
    return nextId;
}

int ProcessorRegistration3D::getId(std::string name) const
{
    Entry entry(name, 0);
    EntryMap::const_iterator it = processorByName.find(entry);
    if (it == processorByName.end()) {
        return 0;
    }
    else {
        return it->second;
    }
}

int ProcessorRegistration3D::getNumId() const
{
    return (int)(processorByNumber.size());
}

std::string ProcessorRegistration3D::getName(int id) const
{
    if (id==0) {
        return std::string("Undefined");
    }
    if (id < 0 || id > (int)processorByNumber.size()) {
        std::stringstream message;
        message << "A processor with ID " << id << " doesn't exist.";
        plbLogicError(message.str());
    }
    return processorByNumber[id-1].name;
}

BoxProcessingFunctional3D* ProcessorRegistration3D::create(std::string procName, std::string data)
{
    int id = getId(procName);
    if (id==0) {
        plbLogicError(std::string("A processor with name ")+procName+" does not exits.");
    }
    return processorByNumber[id-1].factory->create(data);
}

ProcessorRegistration3D::EntryMap::const_iterator
    ProcessorRegistration3D::begin() const
{
    return processorByName.begin();
}

ProcessorRegistration3D::EntryMap::const_iterator
    ProcessorRegistration3D::end() const
{
    return processorByName.end();
}


ProcessorRegistration3D& processorRegistration3D() {
    static ProcessorRegistration3D instance;
    return instance;
}

}  // namespace meta

}  // namespace plb
