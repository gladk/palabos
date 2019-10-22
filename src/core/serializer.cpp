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


#include "parallelism/mpiManager.h"
#include "core/serializer.h"
#include "core/plbDebug.h"
#include <algorithm>

namespace plb {

////////// Free functions ////////////////////////////

void serializerToUnSerializer (
        DataSerializer const* serializer, DataUnSerializer* unSerializer, bool mainProcOnly )
{
    PLB_PRECONDITION( serializer->getSize() == unSerializer->getSize() );
    pluint writePos = 0, readPos = 0;
    pluint serializerBufferSize =0, unSerializerBufferSize =0;
    const char* serializerBuffer =0;
    char* unSerializerBuffer =0;
    while (!unSerializer->isFull()) {
        if (readPos==serializerBufferSize) {
            serializerBuffer = serializer->getNextDataBuffer(serializerBufferSize);
            readPos = 0;
        }
        if (writePos==unSerializerBufferSize) {
            unSerializerBuffer = unSerializer->getNextDataBuffer(unSerializerBufferSize);
            writePos = 0;
        }

        pluint remainToRead = (plint)serializerBufferSize - (plint)readPos;
        pluint remainToWrite = (plint)unSerializerBufferSize - (plint)writePos;
        pluint nextChunk = std::min(remainToRead, remainToWrite);
        for (pluint iChunk=0; iChunk<nextChunk; ++iChunk, ++readPos, ++writePos) {
            if (!mainProcOnly || global::mpi().isMainProcessor()) {
                unSerializerBuffer[writePos] = serializerBuffer[readPos];
            }
        }
        if (writePos==unSerializerBufferSize) {
            unSerializer->commitData();
        }
    }
    delete serializer;
    delete unSerializer;
}

void serializerToSink(DataSerializer const* serializer, SerializedWriter* sink, bool mainProcOnly) {
    if (!mainProcOnly || global::mpi().isMainProcessor()) {
        pluint dataSize = serializer->getSize();
        sink->writeHeader(dataSize);
    }
    while (!serializer->isEmpty()) {
        pluint bufferSize;
        const char* dataBuffer = serializer->getNextDataBuffer(bufferSize);
        if (!mainProcOnly || global::mpi().isMainProcessor()) {
            sink->writeData(dataBuffer, bufferSize);
        }
    }
    delete sink;
    delete serializer;
}

void sourceToUnSerializer(SerializedReader const* source, DataUnSerializer* unSerializer, bool mainProcOnly)
{
    if (!mainProcOnly || global::mpi().isMainProcessor()) {
        pluint dataSize = unSerializer->getSize();
        source->readHeader(dataSize);
    }
    while (!unSerializer->isFull()) {
        pluint bufferSize = 0;
        char* dataBuffer = unSerializer->getNextDataBuffer(bufferSize);
        if (!mainProcOnly || global::mpi().isMainProcessor()) {
            source->readData(dataBuffer, bufferSize);
        }
        unSerializer->commitData();
    }
    delete source;
    delete unSerializer;
}

}  // namespace plb
