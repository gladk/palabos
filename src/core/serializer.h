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


#ifndef SERIALIZER_H
#define SERIALIZER_H

#include "core/globalDefs.h"
#include "core/util.h"

namespace plb {

struct DataSerializer {
    virtual ~DataSerializer() { }
    virtual DataSerializer* clone() const =0;
    virtual pluint getSize() const =0;
    virtual const char* getNextDataBuffer(pluint& bufferSize) const =0;
    virtual bool isEmpty() const =0;
};

struct DataUnSerializer {
    virtual ~DataUnSerializer() { }
    virtual DataUnSerializer* clone() const =0;
    virtual pluint getSize() const =0;
    virtual char* getNextDataBuffer(pluint& bufferSize) =0;
    virtual void commitData() =0;
    virtual bool isFull() const =0;
};

struct SerializedWriter {
    virtual ~SerializedWriter() { }
    virtual SerializedWriter* clone() const =0;
    virtual void writeHeader(pluint dataSize) =0;
    virtual void writeData(char const* dataBuffer, pluint bufferSize) =0;
};


struct SerializedReader {
    virtual ~SerializedReader() { }
    virtual SerializedReader* clone() const =0;
    virtual void readHeader(pluint dataSize) const =0;
    virtual void readData(char* dataBuffer, pluint bufferSize) const =0;
};

/// A serializer-sink for flushing data from a Block into a non-parallel
///   C-array.
template<typename T>
class WriteToSerialArray : public SerializedWriter {
public:
    WriteToSerialArray(T* array_, plint typedSize_);
    virtual WriteToSerialArray* clone() const;
    virtual void writeHeader(pluint dataSize);
    virtual void writeData(char const* dataBuffer, pluint bufferSize);
private:
    char* array;
    plint typedSize;
    plint pos;
};

/// A serializer-source for reading data from a non-parallel
///   C-array into a Block.
template<typename T>
class ReadFromSerialArray : public SerializedReader {
public:
    ReadFromSerialArray(T const* array_, plint typedSize_);
    virtual ReadFromSerialArray<T>* clone() const;
    virtual void readHeader(pluint dataSize) const;
    virtual void readData(char* dataBuffer, pluint bufferSize) const;
private:
    char const* array;
    plint typedSize;
    mutable plint pos;
};


void serializerToUnSerializer(DataSerializer const* serializer, DataUnSerializer* unSerializer, bool mainProcOnly=true);

void serializerToSink(DataSerializer const* serializer, SerializedWriter* sink, bool mainProcOnly=true);

void sourceToUnSerializer(SerializedReader const* source, DataUnSerializer* unSerializer, bool mainProcOnly=true);


} // namespace plb

#endif  // SERIALIZER_H
