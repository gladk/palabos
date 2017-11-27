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

/* Main author: Orestis Malaspinas */

#ifndef VOXELIZER_H
#define VOXELIZER_H

#include "core/globalDefs.h"
#include "atomicBlock/dataProcessingFunctional3D.h"
#include "multiBlock/multiBlockManagement3D.h"
#include "offLattice/triangleHash.h"
#include <memory>

namespace plb {

namespace voxelFlag {

    /// It is a requirement that "undetermined" equals zero, because the
    ///   default initialization value for scalar-fields is zero.
    static const int undetermined = 0;
    static const int outside      = 1;
    /// Cells which are outside, but which have
    ///   inside neighbors.
    static const int outerBorder  = 2;
    static const int inside       = 3;
    /// Cells which are inside, but which have
    ///   outside neighbors.
    static const int innerBorder  = 4;
    static const int toBeInside   = 5;

    int invert(int arg);
    int bulkFlag(int arg);
    int borderFlag(int arg);
    bool insideFlag(int arg);
    bool outsideFlag(int arg);

}

template<typename T>
std::auto_ptr<MultiScalarField3D<int> > voxelize (
        TriangularSurfaceMesh<T> const& mesh,
        plint symmetricLayer, plint borderWidth );

template<typename T>
std::auto_ptr<MultiScalarField3D<int> > voxelize (
        TriangularSurfaceMesh<T> const& mesh,
        Box3D const& domain, plint borderWidth );

// For faster results, the "seed" should contain at least one of the domain corners.
template<typename T>
std::auto_ptr<MultiScalarField3D<int> > voxelize (
        TriangularSurfaceMesh<T> const& mesh,
        Box3D const& domain, plint borderWidth, Box3D seed );

// For faster results, the "seed" should contain at least one of the domain corners.
template<typename T>
std::auto_ptr<MultiScalarField3D<int> > voxelize (
        TriangularSurfaceMesh<T> const& mesh,
        MultiBlockManagement3D const& management,
        plint borderWidth, Box3D seed );


template<typename T>
std::auto_ptr<MultiScalarField3D<int> > revoxelize (
        TriangularSurfaceMesh<T> const& mesh,
        MultiScalarField3D<int>& oldVoxelMatrix,
        MultiContainerBlock3D& hashContainer, plint borderWidth );

template<typename T>
class VoxelizeMeshFunctional3D : public BoxProcessingFunctional3D {
public:
    VoxelizeMeshFunctional3D (
            TriangularSurfaceMesh<T> const& mesh_, bool useFullVoxelizationRange_ = false );
    virtual void processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks );
    virtual VoxelizeMeshFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    bool checkIfFacetsCrossed (
            AtomicContainerBlock3D& hashContainer,
            Array<T,3> const& point1, Array<T,3> const& point2,
            T& distance, plint& whichTriangle );
    bool distanceToSurface (
            AtomicContainerBlock3D& hashContainer,
            Array<T,3> const& point, T& distance, bool& isBehind ) const;
    bool createVoxelizationRange (
            Box3D const& domain, ScalarField3D<int>& voxels,
            Array<plint,2>& xRange, Array<plint,2>& yRange, Array<plint,2>& zRange );
    bool voxelizeFromNeighbor (
        ScalarField3D<int> const& voxels, AtomicContainerBlock3D& hashContainer,
        Dot3D pos, Dot3D neighbor, int& voxelType );
    void printOffender (
            ScalarField3D<int> const& voxels,
            AtomicContainerBlock3D& hashContainer, Dot3D pos );
private:
    TriangularSurfaceMesh<T> const& mesh;
    bool useFullVoxelizationRange;
};

class UndeterminedToFlagFunctional3D : public BoxProcessingFunctional3D_S<int> {
public:
    UndeterminedToFlagFunctional3D(int flag_);
    virtual UndeterminedToFlagFunctional3D* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual void process(Box3D domain, ScalarField3D<int>& voxels);
private:
    int flag;
};

void convertUndeterminedToFlag(MultiScalarField3D<int>& voxels, int flag);

/// Convert inside flags to innerBoundary, and outside flags to outerBoundary,
///   within a layer of width "borderWidth".
template<typename T>
void detectBorderLine( MultiScalarField3D<T>& voxelMatrix,
                       Box3D const& domain, plint borderWidth );

template<typename T>
class DetectBorderLineFunctional3D : public BoxProcessingFunctional3D_S<T> {
public:
    DetectBorderLineFunctional3D(plint borderWidth_);
    virtual void process(Box3D domain, ScalarField3D<T>& voxels);
    virtual DetectBorderLineFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
private:
    plint borderWidth;
};

///// 'Undo' detectBorderLine (i.e. afterwards there are only inside and outside flags!
template<typename T>
void resetBorderFlags( MultiScalarField3D<T>& voxelMatrix,
                       Box3D const& domain );

template<typename T>
class ResetBorderLineFunctional3D : public BoxProcessingFunctional3D_S<T> {
public:
    ResetBorderLineFunctional3D();
    virtual void process(Box3D domain, ScalarField3D<T>& voxels);
    virtual ResetBorderLineFunctional3D<T>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

} // namespace plb

#endif  // VOXELIZER_H
