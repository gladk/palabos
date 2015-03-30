/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2015 FlowKit Sarl
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

/* Main author: Daniel Lagrava
 **/

/** \file
 * Management of multi grid domains -- Implementation
 */

#include "multiGrid/multiGridManagement3D.h"
#include "io/parallelIO.h"

namespace plb {

enum {left=0,right,bottom,top,front,back};

MultiGridManagement3D::MultiGridManagement3D(plint nx, plint ny, plint nz, plint numLevels, plint referenceLevel_)
                                  : referenceLevel(referenceLevel_),
                                    boundingBoxes(numLevels),
                                    coarseGridInterfaces(numLevels),
                                    fineGridInterfaces(numLevels),
                                    bulks(numLevels)                                    
{
    PLB_ASSERT( numLevels > 0 );
    PLB_ASSERT( referenceLevel < numLevels );
    initialize(Box3D(0, nx-1, 0, ny-1, 0, nz-1));
}

MultiGridManagement3D::MultiGridManagement3D(Box3D coarseBoundingBox, plint numLevels, plint referenceLevel_)
                                  : referenceLevel(referenceLevel_),
                                    boundingBoxes(numLevels),
                                    coarseGridInterfaces(numLevels),
                                    fineGridInterfaces(numLevels),
                                    bulks(numLevels)
{
    PLB_ASSERT( numLevels > 0 );
    PLB_ASSERT( referenceLevel < numLevels );
    initialize(coarseBoundingBox);
}

MultiGridManagement3D::~MultiGridManagement3D() {
    delete scaleManager;
}

void MultiGridManagement3D::initialize(Box3D const& level0_box)
{
    scaleManager = global::getDefaultMultiScaleManager().clone();

    // The parameters nx and ny refer to the coarsest lattice (lattice 0).
    boundingBoxes[0] = level0_box;
    // Multiply the values maxX and maxY by two to get, iteratively,
    //   the boundingBoxes of the next finer lattice respectively.
    for (pluint iDim=1; iDim<boundingBoxes.size(); ++iDim) {
        boundingBoxes[iDim] = scaleManager->scaleBox(boundingBoxes[iDim-1], 1);
    }
    
    // Add a block for the full domain of the first instantiated multi-block.
    bulks[referenceLevel].push_back(boundingBoxes[referenceLevel]);
}

MultiGridManagement3D::MultiGridManagement3D(MultiGridManagement3D const& rhs)
  : referenceLevel(rhs.referenceLevel),
    boundingBoxes(rhs.boundingBoxes),
    mpiProcess(rhs.mpiProcess),
    coarseGridInterfaces(rhs.coarseGridInterfaces),
    fineGridInterfaces(rhs.fineGridInterfaces),
    bulks(rhs.bulks),
    scaleManager(rhs.scaleManager->clone())
{}

void MultiGridManagement3D::swap(MultiGridManagement3D& rhs)
{
    std::swap(referenceLevel, rhs.referenceLevel);
    boundingBoxes.swap(rhs.boundingBoxes);
    mpiProcess.swap(rhs.mpiProcess);
    coarseGridInterfaces.swap(rhs.coarseGridInterfaces);
    fineGridInterfaces.swap(rhs.fineGridInterfaces);
    bulks.swap(rhs.bulks);
    std::swap(scaleManager, rhs.scaleManager);
}

Box3D MultiGridManagement3D::getBoundingBox(plint level) const {
    PLB_PRECONDITION(level>=0 && level<(plint)bulks.size());
    return boundingBoxes[level];
}

std::vector<Box3D> const& MultiGridManagement3D::getBulks(plint iLevel) const{
    PLB_PRECONDITION( iLevel>=0 && iLevel<(plint)bulks.size() );
    return bulks[iLevel];
}

std::vector<std::vector<Box3D> > const& MultiGridManagement3D::getBulks() const{
  return bulks;
}

void MultiGridManagement3D::refine(plint coarseLevel, Box3D coarseDomain){
    
    // The finest multi-block, at level blocks.size()-1, cannot be further refined.
    PLB_PRECONDITION( coarseLevel>=0 && coarseLevel<(plint)bulks.size()-1 );
    plint fineLevel = coarseLevel+1;

    // container for booleans that indicate if coarseDomain touches the six walls of the
    // domain. The convention is the following: 0=left,1=right,2=bottom,3=top,4=back,5=front
    bool touches[6] = {false,false,false,false,false,false};
    trimDomain(coarseLevel, coarseDomain, touches);

    // The reduced coarse domain is the one which is going to be excluded from
    //   the original coarse lattice.
    Box3D reducedCoarseDomain(coarseDomain.enlarge(-1));
    // The extended coarse domain it the one which is going to be added
    //   to the original fine lattice.
    Box3D extendedCoarseDomain(coarseDomain.enlarge(1));
    
    // If the domain in question touches a boundary of the multi-block,
    //   both the reduced and the extended coarse domain are
    //   identified with the boundary location.
    if (touches[left]) { // left
        reducedCoarseDomain.x0    -= 1;
        extendedCoarseDomain.x0   += 1;
    }
    if (touches[right]) { // right
        reducedCoarseDomain.x1    += 1;
        extendedCoarseDomain.x1   -= 1;
    }
    if (touches[bottom]) { // bottom
        reducedCoarseDomain.y0    -= 1;
        extendedCoarseDomain.y0   += 1;
    }
    if (touches[top]) { // top
        reducedCoarseDomain.y1    += 1;
        extendedCoarseDomain.y1   -= 1;
    }
    if (touches[front]) { // front
        reducedCoarseDomain.z0    -= 1;
        extendedCoarseDomain.z0   += 1;
    }
    if (touches[back]) { // back
        reducedCoarseDomain.z1    += 1;
        extendedCoarseDomain.z1   -= 1;
    }

    // Extract reduced coarse domain from the original coarse multi-block, for the bulks.
    std::vector<Box3D> exceptedBulks;
    for (pluint iBlock=0; iBlock<bulks[coarseLevel].size(); ++iBlock) {
        except(bulks[coarseLevel][iBlock], reducedCoarseDomain, exceptedBulks);
    }
    exceptedBulks.swap(bulks[coarseLevel]);
    
    // Convert the extended coarse domain to fine units,
    //   and add to the original fine multi-block.
    Box3D extendedFineDomain(extendedCoarseDomain.multiply(2));
    bulks[fineLevel].push_back(extendedFineDomain);
    
    // Define coupling interfaces for all six sides of the coarsened domain, unless they
    //   touch a boundary of the multi-block.
    if (!touches[left]) { ///LEFT
        coarseGridInterfaces[coarseLevel].push_back( Box3D( coarseDomain.x0, coarseDomain.x0,
                                                            coarseDomain.y0, coarseDomain.y1,
                                                            coarseDomain.z0, coarseDomain.z1) );
        fineGridInterfaces[fineLevel].push_back( Box3D( extendedCoarseDomain.x0, extendedCoarseDomain.x0, 
                                                        extendedCoarseDomain.y0, extendedCoarseDomain.y1,
                                                        extendedCoarseDomain.z0, extendedCoarseDomain.z1) );
    }
    if (!touches[right]) { ///RIGHT
        coarseGridInterfaces[coarseLevel].push_back( Box3D( coarseDomain.x1, coarseDomain.x1,
                                                            coarseDomain.y0, coarseDomain.y1,
                                                            coarseDomain.z0, coarseDomain.z1) );
        fineGridInterfaces[fineLevel].push_back( Box3D( extendedCoarseDomain.x1, extendedCoarseDomain.x1, 
                                                        extendedCoarseDomain.y0, extendedCoarseDomain.y1,
                                                        extendedCoarseDomain.z0, extendedCoarseDomain.z1) );
    }
    if (!touches[bottom]) { ///BOTTOM
        coarseGridInterfaces[coarseLevel].push_back( Box3D( coarseDomain.x0, coarseDomain.x1,
                                                            coarseDomain.y0, coarseDomain.y0,
                                                            coarseDomain.z0, coarseDomain.z1) );
        fineGridInterfaces[fineLevel].push_back( Box3D( extendedCoarseDomain.x0, extendedCoarseDomain.x1, 
                                                        extendedCoarseDomain.y0, extendedCoarseDomain.y0,
                                                        extendedCoarseDomain.z0, extendedCoarseDomain.z1) );
    }
    if (!touches[top]) { ///TOP
        coarseGridInterfaces[coarseLevel].push_back( Box3D( coarseDomain.x0, coarseDomain.x1,
                                                            coarseDomain.y1, coarseDomain.y1,
                                                            coarseDomain.z0, coarseDomain.z1) );
        fineGridInterfaces[fineLevel].push_back( Box3D( extendedCoarseDomain.x0, extendedCoarseDomain.x1, 
                                                        extendedCoarseDomain.y1, extendedCoarseDomain.y1,
                                                        extendedCoarseDomain.z0, extendedCoarseDomain.z1) );
    }
    if (!touches[front]) { ///FRONT
        coarseGridInterfaces[coarseLevel].push_back( Box3D( coarseDomain.x0, coarseDomain.x1,
                                                            coarseDomain.y0, coarseDomain.y1,
                                                            coarseDomain.z0, coarseDomain.z0) );
        fineGridInterfaces[fineLevel].push_back( Box3D( extendedCoarseDomain.x0, extendedCoarseDomain.x1, 
                                                        extendedCoarseDomain.y0, extendedCoarseDomain.y1,
                                                        extendedCoarseDomain.z0, extendedCoarseDomain.z0) );
    }
    if (!touches[back]) { ///BACK
        coarseGridInterfaces[coarseLevel].push_back( Box3D( coarseDomain.x0, coarseDomain.x1,
                                                            coarseDomain.y0, coarseDomain.y1,
                                                            coarseDomain.z1, coarseDomain.z1) );
        fineGridInterfaces[fineLevel].push_back( Box3D( extendedCoarseDomain.x0, extendedCoarseDomain.x1, 
                                                        extendedCoarseDomain.y0, extendedCoarseDomain.y1,
                                                        extendedCoarseDomain.z1, extendedCoarseDomain.z1) );
    }
}

void MultiGridManagement3D::coarsen(plint fineLevel, Box3D coarseDomain){
    PLB_PRECONDITION( fineLevel>=1 && fineLevel<(plint)bulks.size() );
    plint coarseLevel = fineLevel-1;

    // First, trim the domain coarseDomain in case it exceeds the extent of the
    //   multi-block, and determine whether coarseDomain touches one of the boundaries
    //   of the multi-block. This information is needed, because the coarse domain
    //   fully replaces the fine domain on boundaries of the multi-block, and there
    //   is therefore no need to create a coarse-fine coupling.
    bool touches[6] = {false,false,false,false,false,false};
    trimDomain(coarseLevel, coarseDomain, touches);

    // Convert the coarse domain to fine units.
    Box3D fineDomain(coarseDomain.multiply(2));
    // The reduced fine domain is the one which is going to be excluded from
    //   the original fine lattice.
    Box3D reducedFineDomain(fineDomain.enlarge(-1));
    
    // The extended coarse domain is the one which is going to be added
    //   to the original coarse lattice.
    Box3D extendedCoarseDomain(coarseDomain.enlarge(1));
    
    // If the domain in question touches a boundary of the multi-block,
    //   both the reduced fine domain and the extended coarse domain are
    //   identified with the boundary location.
    if (touches[left]) { // left
        extendedCoarseDomain.x0 += 1;
        reducedFineDomain.x0 -= 1;
    }
    if (touches[right]) { // right
        extendedCoarseDomain.x1 -= 1;
        reducedFineDomain.x1 += 1;
    }
    if (touches[bottom]) { // bottom
        extendedCoarseDomain.y0 += 1;
        reducedFineDomain.y0 -= 1;
    }
    if (touches[top]) { // top
        extendedCoarseDomain.y1 -= 1;
        reducedFineDomain.y1 += 1;
    }
    if (touches[front]) { // front
        extendedCoarseDomain.z0 += 1;
        reducedFineDomain.z0 -= 1;
    }
    if (touches[back]) { // back
        extendedCoarseDomain.z1 -= 1;
        reducedFineDomain.z1 += 1;
    }
    
    // Extract reduced fine domain from the original fine multi-block.
    std::vector<Box3D> exceptedBlocks;
    for (pluint iBlock=0; iBlock<bulks[fineLevel].size(); ++iBlock) {
        except(bulks[fineLevel][iBlock], reducedFineDomain, exceptedBlocks);
    }
    exceptedBlocks.swap(bulks[fineLevel]);
    
    // Add extended coarse domain to the original coarse multi-block.
    bulks[coarseLevel].push_back(extendedCoarseDomain);
    
    pcout << "new fine grid : " << std::endl;
    for (pluint iComp=0; iComp<bulks[fineLevel].size(); ++iComp){
        Box3D box(bulks[fineLevel][iComp]);
        pcout << " >>>>>> " << box.x0 << " " << box.x1 << " " << box.y0 << " " 
              << box.y1 << " " << box.z0 << " " << box.z1 << std::endl;
    }
    
    
    pcout << "Added to the coarse grid : " << std::endl;
    Box3D box(extendedCoarseDomain);
    pcout << " >>> " << box.x0 << " " << box.x1 << " " << box.y0 << " " 
              << box.y1 << " " << box.z0 << " " << box.z1 << std::endl;
    
    // Define coupling interfaces for all four sides of the refined domain, unless they
    //   touch a boundary of the multi-block.
    if (!touches[left]) { ///LEFT
        fineGridInterfaces[fineLevel].push_back( Box3D( coarseDomain.x0, coarseDomain.x0,
                                                        coarseDomain.y0, coarseDomain.y1,
                                                        coarseDomain.z0, coarseDomain.z1) );
        coarseGridInterfaces[coarseLevel].push_back( Box3D( extendedCoarseDomain.x0, extendedCoarseDomain.x0, 
                                                          extendedCoarseDomain.y0, extendedCoarseDomain.y1,
                                                          extendedCoarseDomain.z0, extendedCoarseDomain.z1) );
    }
    if (!touches[right]) { ///RIGHT
        fineGridInterfaces[fineLevel].push_back( Box3D( coarseDomain.x1, coarseDomain.x1,
                                                          coarseDomain.y0, coarseDomain.y1,
                                                          coarseDomain.z0, coarseDomain.z1) );
        coarseGridInterfaces[coarseLevel].push_back( Box3D( extendedCoarseDomain.x1, extendedCoarseDomain.x1, 
                                                        extendedCoarseDomain.y0, extendedCoarseDomain.y1,
                                                        extendedCoarseDomain.z0, extendedCoarseDomain.z1) );
    }
    if (!touches[bottom]) { ///BOTTOM
        fineGridInterfaces[fineLevel].push_back( Box3D( coarseDomain.x0, coarseDomain.x1,
                                                          coarseDomain.y0, coarseDomain.y0,
                                                          coarseDomain.z0, coarseDomain.z1) );
        coarseGridInterfaces[coarseLevel].push_back( Box3D( extendedCoarseDomain.x0, extendedCoarseDomain.x1, 
                                                        extendedCoarseDomain.y0, extendedCoarseDomain.y0,
                                                        extendedCoarseDomain.z0, extendedCoarseDomain.z1) );
    }
    if (!touches[top]) { ///TOP
        fineGridInterfaces[fineLevel].push_back( Box3D( coarseDomain.x0, coarseDomain.x1,
                                                          coarseDomain.y1, coarseDomain.y1,
                                                          coarseDomain.z0, coarseDomain.z1) );
        coarseGridInterfaces[coarseLevel].push_back( Box3D( extendedCoarseDomain.x0, extendedCoarseDomain.x1, 
                                                          extendedCoarseDomain.y1, extendedCoarseDomain.y1,
                                                          extendedCoarseDomain.z0, extendedCoarseDomain.z1) );
    }
    if (!touches[front]) { ///FRONT
        fineGridInterfaces[fineLevel].push_back( Box3D( coarseDomain.x0, coarseDomain.x1,
                                                          coarseDomain.y0, coarseDomain.y1,
                                                          coarseDomain.z0, coarseDomain.z0) );
        coarseGridInterfaces[coarseLevel].push_back( Box3D( extendedCoarseDomain.x0, extendedCoarseDomain.x1, 
                                                          extendedCoarseDomain.y0, extendedCoarseDomain.y1,
                                                          extendedCoarseDomain.z0, extendedCoarseDomain.z0) );
    }
    if (!touches[back]) { ///BACK
        fineGridInterfaces[fineLevel].push_back( Box3D( coarseDomain.x0, coarseDomain.x1,
                                                          coarseDomain.y0, coarseDomain.y1,
                                                          coarseDomain.z1, coarseDomain.z1) );
        coarseGridInterfaces[coarseLevel].push_back( Box3D( extendedCoarseDomain.x0, extendedCoarseDomain.x1, 
                                                            extendedCoarseDomain.y0, extendedCoarseDomain.y1,
                                                            extendedCoarseDomain.z1, extendedCoarseDomain.z1) );
    }
    
    pcout << "Fine Interfaces : " << std::endl;
    for (pluint iComp=0; iComp < fineGridInterfaces[fineLevel].size(); ++iComp){
        Box3D box(fineGridInterfaces[fineLevel][iComp]);
        pcout << " >>> " << box.x0 << " " << box.x1 << " " << box.y0 << " " 
              << box.y1 << " " << box.z0 << " " << box.z1 << std::endl;
    }
    pcout << "....\n";
    pcout << "Coarse Interfaces : " << std::endl;
    for (pluint iComp=0; iComp < coarseGridInterfaces[coarseLevel].size(); ++iComp){
        Box3D box(coarseGridInterfaces[coarseLevel][iComp]);
        pcout << " >>> " << box.x0 << " " << box.x1 << " " << box.y0 << " " 
              << box.y1 << " " << box.z0 << " " << box.z1 << std::endl;
    }
}

std::vector<std::vector<Box3D> > const& MultiGridManagement3D::getCoarseInterface() const
{
    return coarseGridInterfaces;
}

std::vector<std::vector<Box3D> > const& MultiGridManagement3D::getFineInterface() const
{
    return fineGridInterfaces;
}


void MultiGridManagement3D::trimDomain (plint whichLevel, Box3D& domain,
                                        bool* touches) const
{
    Box3D bbox = boundingBoxes[whichLevel];
    /// LEFT
    if (domain.x0 <= bbox.x0) {
        domain.x0 = bbox.x0; // Trim in case multi-block extent is exceeded.
        touches[left] = true;
    }
    /// BOTTOM
    if (domain.y0 <= bbox.y0) {
        domain.y0 = bbox.y0; // Trim in case multi-block extent is exceeded.
        touches[bottom] = true;
    }
    /// FRONT
    if (domain.z0 <= bbox.z0) {
        domain.z0 = bbox.z0; // Trim in case multi-block extent is exceeded.
        touches[front] = true;
    }
    /// RIGHT
    if (domain.x1 >= bbox.x1) {
        domain.x1 = bbox.x1; // Trim in case multi-block extent is exceeded.
        touches[right] = true;
    }
    /// TOP
    if (domain.y1 >= bbox.y1) {
        domain.y1 = bbox.y1; // Trim in case multi-block extent is exceeded.
        touches[top] = true;
    }
    /// BACK
    if (domain.z1 >= bbox.z1) {
        domain.z1 = bbox.z1; // Trim in case multi-block extent is exceeded.
        touches[back] = true;
    }
}

/**
 * Use one of the Parallelizer2D objects to recompute the parallelization of the
 *   domain.
 */
void MultiGridManagement3D::parallelize(Parallelizer3D* parallelizer){
    parallelizer->parallelize();
    bulks = parallelizer->getRecomputedBlocks();
    mpiProcess = parallelizer->getMpiDistribution();
    delete parallelizer;
}


plint MultiGridManagement3D::getNumLevels() const {
    return (plint) bulks.size();
}

plint MultiGridManagement3D::getReferenceLevel() const {
    return referenceLevel;
}

std::vector<std::vector<plint> > const& MultiGridManagement3D::getMpiProcesses() const{
  return mpiProcess;
}

/// Extract a domain (in coarse coordinates) from a MultiGridManagement3D
/** This is achieved by taking all the geometric information of the MultiGridManagement3D
 *  and intersecting it with the so called coarsestDomain. All the intersecting parts
 *  are kept on the new MultiGridManagement3D.
 */
MultiGridManagement3D extractManagement(MultiGridManagement3D management, Box3D coarsestDomain, bool crop){
    //TODO: implement
    return management;
}


} // namespace plb

