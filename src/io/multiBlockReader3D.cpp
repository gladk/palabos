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
 * I/O routines for 3D multiblock -- implementation.
 */

#include "core/globalDefs.h"
#include "io/multiBlockReader3D.h"
#include "io/mpiParallelIO.h"
#include "parallelism/mpiManager.h"
#include "libraryInterfaces/TINYXML_xmlIO.h"
#include "libraryInterfaces/TINYXML_xmlIO.hh"
#include "core/util.h"
#include "core/plbTypenames.h"
#include "core/multiBlockIdentifiers3D.h"
#include "core/processorIdentifiers3D.h"
#include "multiBlock/nonLocalTransfer3D.h"
#include "multiBlock/multiBlockOperations3D.h"
#include "io/plbFiles.h"
#include <numeric>
#include <algorithm>
#include <memory>

namespace plb {

namespace parallelIO {

/***** 1. Multi-Block Reader **************************************************/

void dumpRestoreData( MultiBlock3D& multiBlock, bool dynamicContent,
                      std::vector<plint> const& myBlockIds, std::vector<std::vector<char> > const& data,
                      std::map<int,std::string> const& foreignIds )
{
    modif::ModifT typeOfVariables = dynamicContent ? modif::dataStructure : modif::staticVariables;
    for (pluint iBlock=0; iBlock<myBlockIds.size(); ++iBlock) {
        plint blockId = myBlockIds[iBlock];
        SmartBulk3D bulk(multiBlock.getMultiBlockManagement(), blockId);
        Box3D localBulk(bulk.toLocal(bulk.getBulk()));
        AtomicBlock3D& block = multiBlock.getComponent(blockId);
        block.getDataTransfer().receive(localBulk, data[iBlock], typeOfVariables, foreignIds);
    }
    multiBlock.getBlockCommunicator().duplicateOverlaps(multiBlock, typeOfVariables);
}

void createDynamicsForeignIds3D(FileName fName, std::map<int,std::string>& foreignIds)
{
    foreignIds.clear();
    std::vector<MultiBlock3D::ProcessorStorage3D> processors;
    fName.defaultPath(global::directories().getInputDir());
    fName.setExt("plb");
    XMLreader reader(fName);
    XMLreaderProxy dynReader(0);
    try {
        dynReader = reader["Block3D"]["Data"]["DynamicsDict"];
    }
    catch(PlbIOException const&) {
        return;
    }
    std::vector<XMLreader*> const& dynItems = dynReader.getChildren();
    std::map<int,std::string> foreignIdToName;
    for (pluint i=0; i<dynItems.size(); ++i) {
        std::string dynamicsName = dynItems[i]->getName();
        int dynamicsId;
        XMLreaderProxy(dynItems[i]).read(dynamicsId);
        foreignIds.insert(std::pair<int,std::string>(dynamicsId,dynamicsName));
    }
}

void readXmlProcessors(FileName fName, MultiBlock3D& block) {
    std::vector<MultiBlock3D::ProcessorStorage3D> processors;
    fName.defaultPath(global::directories().getInputDir());
    fName.setExt("plb");
    XMLreader reader(fName);
    XMLreaderProxy procReader(0);
    try {
        procReader = reader["Block3D"]["Data"]["Processor"];
    }
    catch(PlbIOException const&) {
        return;
    }
    std::vector<MultiBlock3D*> processorPartnerList;
    processorPartnerList.push_back(&block);
    for (; procReader.isValid(); procReader = procReader.iterId()) {
        std::string processorName, data;
        Box3D domain; Array<plint,6> domain_arr;
        plint level;
        std::vector<id_t> blocks;
        procReader["Name"].read(processorName);
        procReader["Data"].read(data);
        procReader["Domain"].read<plint,6>(domain_arr);
        domain.from_plbArray(domain_arr);
        procReader["Level"].read(level);
        procReader["Blocks"].read(blocks);

        BoxProcessingFunctional3D* functional = meta::processorRegistration3D().create(processorName, data);
        MultiBlock3D::ProcessorStorage3D newStorage (
                         BoxProcessorGenerator3D(functional, domain),
                         processorPartnerList, level );

        processors.push_back(newStorage);
    }
    for (pluint iProcessor=0; iProcessor<processors.size(); ++iProcessor) {
        DataProcessorGenerator3D* newGenerator = processors[iProcessor].getGenerator().clone();
        if (newGenerator->extract(block.getBoundingBox())) {
            addInternalProcessor( *newGenerator,
                                  processors[iProcessor].getMultiBlocks(),
                                  processors[iProcessor].getLevel() );
        }
        delete newGenerator;
    }
}

void readXmlSpec (
    FileName fName, Box3D& boundingBox, std::vector<plint>& offsets,
    plint& envelopeWidth, plint& gridLevel, plint& cellDim,
    std::string& dataType, std::string& descriptor, std::string& family,
    std::vector<Box3D>& components, bool& dynamicContent, FileName& data_fName )
{
    fName.defaultPath(global::directories().getInputDir());
    fName.setExt("plb");
    XMLreader reader(fName);
    std::string ordering;
    bool forwardOrdering;
    plint numComponents;
    Array<plint,6> boundingBox_array;
    std::string data_fName_str;

    reader["Block3D"]["General"]["Family"].read(family);
    reader["Block3D"]["General"]["Datatype"].read(dataType);
    try {
        reader["Block3D"]["General"]["Descriptor"].read(descriptor);
    }
    catch(PlbIOException const&) {
        descriptor = "NA";
    }
    reader["Block3D"]["General"]["cellDim"].read(cellDim);
    reader["Block3D"]["General"]["dynamicContent"].read(dynamicContent);
    reader["Block3D"]["Structure"]["BoundingBox"].read<plint,6>(boundingBox_array);
    boundingBox.from_plbArray(boundingBox_array);
    reader["Block3D"]["Structure"]["NumComponents"].read(numComponents);
    reader["Block3D"]["Structure"]["EnvelopeWidth"].read(envelopeWidth);
    reader["Block3D"]["Structure"]["GridLevel"].read(gridLevel);
    reader["Block3D"]["Data"]["Offsets"].read(offsets);
    reader["Block3D"]["Data"]["File"].read(data_fName_str);
    // If the filename for the data has no path specification, it is taken
    //   to be in the same directory as the xml file.
    data_fName = FileName(data_fName_str).defaultPath(fName.getPath());
    data_fName.defaultExt("dat");
    try {
        reader["Block3D"]["Data"]["IndexOrdering"].read(ordering);
    }
    catch(PlbIOException const&) {
        ordering = "";
    }


    if( (plint)offsets.size() != numComponents ) {
        plbIOError(std::string("Number of offsets does not match number of components in XML file."));
    }

    forwardOrdering = true;
    if (ordering != "") {
        if( ordering=="zIsFastest") {
            forwardOrdering = true;
        }
        else if( ordering=="xIsFastest") {
            forwardOrdering = false;
        }
        else {
            plbIOError(std::string("Invalid ordering: ")+ordering);
        }
        if (!forwardOrdering) {
            plbIOError(std::string("Backward ordering not accepted while reading data."));
        }
    }

    SparseBlockStructure3D blockStructure(boundingBox);

    XMLreaderProxy comp = reader["Block3D"]["Data"]["Component"];
    for (; comp.isValid(); comp = comp.iterId()) {
        Array<plint,6> component_array;
        comp.read<plint,6>(component_array);
        Box3D component;
        component.from_plbArray(component_array);
        components.push_back(component);
    }
    if( numComponents != (plint)components.size() ) {
        plbIOError(std::string("Actual number of components does not match the claimed number in XML file."));
    }
}

MultiBlock3D* load3D(FileName fName)
{
    Box3D boundingBox;
    std::vector<plint> offsets;
    plint envelopeWidth, gridLevel;
    std::string dataType, descriptor, family;
    FileName data_fName;
    std::vector<Box3D> components;
    bool dynamicContent;
    plint cellDim;
    readXmlSpec( fName, boundingBox, offsets, envelopeWidth, gridLevel, cellDim, dataType,
                 descriptor, family, components, dynamicContent, data_fName );

    SparseBlockStructure3D blockStructure(boundingBox);
    for( plint iComponent=0; iComponent<(plint)components.size(); ++iComponent) {
        blockStructure.addBlock(components[iComponent], iComponent);
    }

    ExplicitThreadAttribution* threadAttribution = new ExplicitThreadAttribution;
    std::vector<std::pair<plint,plint> > blockRanges;
    plint numBlocks = offsets.size();
    plint numRanges = std::min(numBlocks, (plint)global::mpi().getSize());
    util::linearRepartition(0, numBlocks-1, numRanges, blockRanges);
    std::vector<plint> myBlockIds;
    for (plint iThread=0; iThread<(plint)blockRanges.size(); ++iThread) {
        for (plint iBlock=blockRanges[iThread].first; iBlock<=blockRanges[iThread].second; ++iBlock) {
            threadAttribution->addBlock(iBlock, iThread);
            if (iThread==global::mpi().getRank()) {
                myBlockIds.push_back(iBlock);
            }
        }
    }

    MultiBlockManagement3D management(blockStructure, threadAttribution, envelopeWidth, gridLevel);

    MultiBlock3D* newBlock =
                meta::multiBlockRegistration3D().generate (
                    dataType, descriptor, family, management, cellDim );
    PLB_ASSERT( newBlock );
    std::vector<std::vector<char> > data(myBlockIds.size());
    loadRawData( data_fName, myBlockIds, offsets, data);
    std::map<int,std::string> foreignIds;
    createDynamicsForeignIds3D(fName, foreignIds);
    dumpRestoreData(*newBlock, dynamicContent, myBlockIds, data, foreignIds);
    readXmlProcessors(fName, *newBlock);
    return newBlock;
}

void load(FileName fName, MultiBlock3D& intoBlock, bool dynamicContent )
{
    std::auto_ptr<MultiBlock3D> loadedBlock ( load3D(fName) );
    modif::ModifT typeOfVariables = dynamicContent ?
            modif::dataStructure : modif::staticVariables;
    copy_generic( *loadedBlock, loadedBlock->getBoundingBox(),
                  intoBlock, intoBlock.getBoundingBox(), typeOfVariables );
}


SavedFullMultiBlockSerializer3D::SavedFullMultiBlockSerializer3D(FileName fName)
{
    fName.defaultPath(global::directories().getInputDir());
    fName.setExt("plb");
    XMLreader reader(fName);
    std::string family, ordering;
    plint numComponents, numberOfBytes;
    Array<plint,6> boundingBox_array;
    std::string data_fName_str;

    reader["Block3D"]["General"]["Family"].read(family);
    reader["Block3D"]["General"]["Datatype"].read(str_dataType);
    reader["Block3D"]["General"]["cellDim"].read(cellDim);
    reader["Block3D"]["Structure"]["BoundingBox"].read<plint,6>(boundingBox_array);
    boundingBox.from_plbArray(boundingBox_array);
    reader["Block3D"]["Structure"]["NumComponents"].read(numComponents);
    reader["Block3D"]["Data"]["Offsets"].read(numberOfBytes);
    reader["Block3D"]["Data"]["File"].read(data_fName_str);
    reader["Block3D"]["Data"]["IndexOrdering"].read(ordering);

    data_fName = FileName(data_fName_str);
    data_fName.defaultPath(global::directories().getInputDir());
    data_fName.defaultExt("dat");

    if( ordering=="zIsFastest") {
        forwardOrdering = true;
    }
    else if( ordering=="xIsFastest") {
        forwardOrdering = false;
    }
    else {
        plbIOError(std::string("Invalid ordering: ")+ordering);
    }

    if (!( (family=="ScalarField3D") || (family=="TensorField3D") ||
           (family=="NTensorField3D") ) )
    {
        plbIOError ( std::string (
            "Was expecting to read a Palabos data file of type ScalarField3D, "\
            "TensorField3D or NTensorField3D, but got \"" ) + family + "\"" );
    }
    if (!( (str_dataType=="float") || (str_dataType=="double") || (family=="int") ) )
    {
        plbIOError ( std::string (
            "Was expecting to read a Palabos data file of type float, "\
            "double or int, but got \"" ) + family + "\"" );
    }
    if (!( cellDim==1 || cellDim==3 || cellDim==6)) {
        std::stringstream compstr;
        compstr << cellDim;
        plbIOError ( std::string (
            "Was expecting to read a Palabos data file with a structure of 1, "\
            "3, or 6 components, but got " ) + compstr.str() + " components" );
    }
    if (numComponents != 1) {
        plbIOError( "Cannot convert a sparse-block object. Save the structure "\
                    "as a full block." );
    }

    typeSize=NativeTypeConstructor(str_dataType).getTypeSize();
    if (typeSize*cellDim*boundingBox.nCells() != numberOfBytes) {
        plbIOError("Number of bytes does not match up with block dimensions in file.");
    }
    pos = 0;
    sizeOfChunk = 1000000; // Treat 1 MB at a time.
    if (global::mpi().isMainProcessor()) {
        fp = fopen(data_fName.get().c_str(), "rb");
    }
    else {
        fp = 0;
    }
    plbMainProcIOError(!fp, "Could not open file "+data_fName.get());
}

SavedFullMultiBlockSerializer3D::SavedFullMultiBlockSerializer3D (
        SavedFullMultiBlockSerializer3D const& rhs )
    : boundingBox(rhs.boundingBox),
      cellDim(rhs.cellDim), typeSize(rhs.typeSize),
      sizeOfChunk(rhs.sizeOfChunk),
      data_fName(rhs.data_fName),
      str_dataType(rhs.str_dataType),
      forwardOrdering(rhs.forwardOrdering),
      pos(rhs.pos)
{  
    if (global::mpi().isMainProcessor()) {
        fp = fopen(data_fName.get().c_str(), "rb");
#if defined PLB_MAC_OS_X || defined PLB_BSD
        fseek(fp, (long int)pos, SEEK_SET);
#else
        fseeko64(fp, pos, SEEK_SET);
#endif
    }
    else {
        fp = 0;
    }
    plbMainProcIOError(!fp, "Could not open file "+data_fName.get());
}

SavedFullMultiBlockSerializer3D::~SavedFullMultiBlockSerializer3D()
{
    if (global::mpi().isMainProcessor()) {
        fclose(fp);
    }
}

SavedFullMultiBlockSerializer3D* SavedFullMultiBlockSerializer3D::clone() const
{
    return new SavedFullMultiBlockSerializer3D(*this);
}

pluint SavedFullMultiBlockSerializer3D::getSize() const {
    return typeSize*cellDim*boundingBox.nCells();
}

const char* SavedFullMultiBlockSerializer3D::getNextDataBuffer(pluint& bufferSize) const
{
    bufferSize = std::min(sizeOfChunk, (plint)(getSize()-pos));
    buffer.resize(bufferSize);
    plint numRead=0;
    if (global::mpi().isMainProcessor()) {
        numRead = (plint) fread( &buffer[0], 1, bufferSize, fp );
    }
    pos += bufferSize;
    plbMainProcIOError(numRead!=(plint)bufferSize, "Error while reading from file "+data_fName.get());
    return &buffer[0];
}

bool SavedFullMultiBlockSerializer3D::isEmpty() const {
    return pos >= (plint)getSize();
}

plint SavedFullMultiBlockSerializer3D::getCellDim() const {
    return cellDim;
}

std::string SavedFullMultiBlockSerializer3D::dataType() const {
    return str_dataType;
}

Box3D SavedFullMultiBlockSerializer3D::getBoundingBox() const {
    return boundingBox;
}

bool SavedFullMultiBlockSerializer3D::orderingIsForward() const {
    return forwardOrdering;
}


}  // namespace parallelIO

}  // namespace plb
