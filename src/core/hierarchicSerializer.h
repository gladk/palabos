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
 * Serialization and unserialization of generic data.
 */
#ifndef HIERARCHIC_SERIALIZER_H
#define HIERARCHIC_SERIALIZER_H

#include "core/dynamics.h"
#include <map>
#include <cstring>

namespace plb {

class HierarchicUnserializer {
public:
    HierarchicUnserializer(std::vector<char> const& data_, pluint currentPos_)
        : data(data_),
          currentPos(currentPos_),
          idIndirect(0)
    {
        readHeader();
    }
    HierarchicUnserializer(std::vector<char> const& data_, pluint currentPos_,
                           std::map<int,int> const* idIndirect_)
        : data(data_),
          currentPos(currentPos_),
          idIndirect(idIndirect_)
    {
        readHeader();
    }
    int getNumObjects() const {
        return numObjects;
    }
    int getId() const {
        if (!idIndirect) {
            return id;
        }
        else {
            std::map<int,int>::const_iterator it = idIndirect->find(id);
            PLB_ASSERT( it != idIndirect->end() );
            return it->second;
        }
    }
    int getNumVal() const {
        return numValInObject;
    }
    pluint getCurrentPos() const {
        return currentPos;
    }
    template<typename T> bool readValue(T& value) {
        read(value);
        --numRepetitions;
        --numValInObject;
        return advance();
    }
    template<typename T> T readValue() {
        T value;
        readValue<T>(value);
        return value;
    }
    template<typename T> bool readValues(std::vector<T>& values) {
        PLB_PRECONDITION((plint)values.size()<=numRepetitions);
        for (pluint iVal=0; iVal<values.size(); ++iVal) {
            read(values[iVal]);
        }
        numRepetitions -= values.size();
        numValInObject -= values.size();
        return advance();
    }
    template<typename T, int ndim> bool readValues(Array<T,ndim>& values) {
        for (plint iVal=0; iVal<ndim; ++iVal) {
            read(values[iVal]);
        }
        numRepetitions -= ndim;
        numValInObject -= ndim;
        return advance();
    }
    void remapId(std::map<int,int> const& indirect) {
        std::map<int,int>::const_iterator it = indirect.find(id);
        PLB_ASSERT( it != indirect.end() );
        id = it->second;
    }
private:
    void readHeader() {
        read(numObjects);
        read(id);
        read(numValInObject);
        read(currentTypeSize);
        read(numRepetitions);
    }
    template<typename T> void read(T& value) {
        PLB_PRECONDITION(currentPos+sizeof(value)<=data.size());
        if (sizeof(value)>0) {
            //value = *( (T const*) (&data[currentPos]) );
            memcpy((void*)(&value), (void*)(&data[currentPos]), sizeof(T));
            currentPos += sizeof(value);
        }
    }
    bool advance() {
        if (numValInObject==0) {
            --numObjects;
            if (numObjects>0) {
                read(id);
                read(numValInObject);
            }
            else {
                // Return if there's no object left, to avoid that the
                //   next repetition is mistakingly parsed.
                return false;
            }
        }
        if (numRepetitions==0) {
            read(currentTypeSize);
            read(numRepetitions);
        }
        return true;
    }
private:
    std::vector<char> const& data;
    pluint currentPos;
    int    numObjects;
    int    id;
    int    numValInObject;
    int    numRepetitions;
    int    currentTypeSize;
    std::map<int,int> const* idIndirect;
};

class HierarchicSerializer {
public:
    HierarchicSerializer(std::vector<char>& data_, int topMostObjectId)
        : data(data_),
          currentPos(data.size())
    {
        numDynamicsPos = assignPointer();
        nextDynamics(topMostObjectId);
    }
    void nextDynamics(int id) {
        write(id);
        numValInDynamicsPos = assignPointer();
        currentTypeSize = 0;
        incrementInt(numDynamicsPos);
    }
    template<typename T> void addValue(T value)
    {
        PLB_PRECONDITION( currentPos > sizeof(int) );
        if (sizeof(value) != currentTypeSize) {
            currentTypeSize = sizeof(value);
            write(currentTypeSize);
            numRepetitionsPos = assignPointer();
        }
        write(value);
        incrementInt(numRepetitionsPos);
        incrementInt(numValInDynamicsPos);
    }
    template<typename T> void addValues(std::vector<T> const& values)
    {
        PLB_PRECONDITION( currentPos > sizeof(int) );
        if (values.empty()) return;
        addValue(values[0]);
        for (pluint iVal=1; iVal<values.size(); ++iVal) {
            write(values[iVal]);
            incrementInt(numRepetitionsPos);
            incrementInt(numValInDynamicsPos);
        }
    }
    template<typename T, int ndim> void addValues(Array<T,ndim> const& values)
    {
        PLB_PRECONDITION( ndim>0 );
        addValue(values[0]);
        for (pluint iVal=1; iVal<ndim; ++iVal) {
            write(values[iVal]);
            incrementInt(numRepetitionsPos);
            incrementInt(numValInDynamicsPos);
        }
    }
private:
    template<typename T> void write(T value) {
        if (sizeof(value)==0) return;
        data.resize(currentPos+sizeof(value));
        //*( (T*) (&data[currentPos]) ) = value;
        memcpy((void*)(&data[currentPos]), (void*)(&value), sizeof(value));
        currentPos += sizeof(T);
    }
    pluint assignPointer() {
        data.resize(currentPos+sizeof(int));
        pluint ptr = currentPos;
        int zero=0;
        memcpy((void*)(&data[ptr]), (void*)(&zero), sizeof(int));
        //*( (int*) (&data[ptr]) ) = 0;
        currentPos += sizeof(int);
        return ptr;
    }
    void incrementInt(pluint pos) {
        PLB_PRECONDITION(pos+sizeof(int)<=data.size());
        //++(*( (int*) (&data[pos]) ) );
        int value;
        memcpy((void*)(&value), (void*)(&data[pos]), sizeof(int));
        ++value;
        memcpy((void*)(&data[pos]), (void*)(&value), sizeof(int));
    }
private:
    std::vector<char>& data;
    pluint currentPos;
    pluint numDynamicsPos;
    pluint numValInDynamicsPos;
    pluint numRepetitionsPos;
    int    currentTypeSize;
};

}  // namespace plb

#endif  // HIERARCHIC_SERIALIZER_H
