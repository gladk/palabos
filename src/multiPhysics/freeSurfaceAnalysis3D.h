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

#ifndef FREE_SURFACE_ANALYSIS_3D_H
#define FREE_SURFACE_ANALYSIS_3D_H

#include "multiBlock/multiBlockLattice3D.h"

namespace plb {

template<typename T, template<typename U> class Descriptor>
T freeSurfaceAverageMass(std::vector<MultiBlock3D*> twoPhaseArgs, Box3D domain);
             
template<typename T, template<typename U> class Descriptor>
class FS_AverageMassFunctional3D : public PlainReductiveBoxProcessingFunctional3D
{
public:
    FS_AverageMassFunctional3D();
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks);
    virtual FS_AverageMassFunctional3D<T,Descriptor>* clone() const;
    T getAverageMass() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        for (pluint i=0; i<modified.size(); ++i) {
            modified[i] = modif::nothing;
        }
    }
private:
    plint averageMassId;
};


template<typename T, template<typename U> class Descriptor>
T freeSurfaceTotalMass(std::vector<MultiBlock3D*> twoPhaseArgs, Box3D domain);
             
template<typename T, template<typename U> class Descriptor>
class FS_TotalMassFunctional3D : public PlainReductiveBoxProcessingFunctional3D
{
public:
    FS_TotalMassFunctional3D();
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks);
    virtual FS_TotalMassFunctional3D<T,Descriptor>* clone() const;
    T getTotalMass() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        for (pluint i=0; i<modified.size(); ++i) {
            modified[i] = modif::nothing;
        }
    }
private:
    plint totalMassId;
};


template<typename T, template<typename U> class Descriptor>
T freeSurfaceAverageDensity(std::vector<MultiBlock3D*> twoPhaseArgs, Box3D domain);
             
template<typename T, template<typename U> class Descriptor>
class FS_AverageDensityFunctional3D : public PlainReductiveBoxProcessingFunctional3D
{
public:
    FS_AverageDensityFunctional3D();
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks);
    virtual FS_AverageDensityFunctional3D<T,Descriptor>* clone() const;
    T getAverageDensity() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        for (pluint i=0; i<modified.size(); ++i) {
            modified[i] = modif::nothing;
        }
    }
private:
    plint averageDensityId;
};


template<typename T, template<typename U> class Descriptor>
T freeSurfaceAverageVolumeFraction(std::vector<MultiBlock3D*> twoPhaseArgs, Box3D domain);
             
template<typename T, template<typename U> class Descriptor>
class FS_AverageVolumeFractionFunctional3D : public PlainReductiveBoxProcessingFunctional3D
{
public:
    FS_AverageVolumeFractionFunctional3D();
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks);
    virtual FS_AverageVolumeFractionFunctional3D<T,Descriptor>* clone() const;
    T getAverageVolumeFraction() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        for (pluint i=0; i<modified.size(); ++i) {
            modified[i] = modif::nothing;
        }
    }
private:
    plint averageVfId;
};


template<typename T, template<typename U> class Descriptor>
plint countFreeSurfaceElements(std::vector<MultiBlock3D*> twoPhaseArgs, plint flagToLookFor, Box3D domain);

template<typename T, template<typename U> class Descriptor>
class CountFreeSurfaceElementsFunctional3D : public PlainReductiveBoxProcessingFunctional3D
{
public:
    CountFreeSurfaceElementsFunctional3D(plint flagToLookFor_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks);
    virtual CountFreeSurfaceElementsFunctional3D<T,Descriptor>* clone() const;
    plint getNumInterfaceCells() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        for (pluint i=0; i<modified.size(); ++i) {
            modified[i] = modif::nothing;
        }
    }
private:
    plint numCellsId;
    plint flagToLookFor;
};
    

template<typename T, template<typename U> class Descriptor>
Array<T,3> freeSurfaceAverageMomentum(std::vector<MultiBlock3D*> twoPhaseArgs, Box3D domain);
             
template<typename T, template<typename U> class Descriptor>
class FS_AverageMomentumFunctional3D : public PlainReductiveBoxProcessingFunctional3D
{
public:
    FS_AverageMomentumFunctional3D();
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks);
    virtual FS_AverageMomentumFunctional3D<T,Descriptor>* clone() const;
    Array<T,3> getAverageMomentum() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        for (pluint i=0; i<modified.size(); ++i) {
            modified[i] = modif::nothing;
        }
    }
private:
    Array<plint,3> averageMomentumId;
};


template<typename T, template<typename U> class Descriptor>
T freeSurfaceAverageHeight(std::vector<MultiBlock3D*> twoPhaseArgs, Box3D domain);
             
template<typename T, template<typename U> class Descriptor>
class FS_AverageHeightFunctional3D : public PlainReductiveBoxProcessingFunctional3D
{
public:
    FS_AverageHeightFunctional3D();
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks);
    virtual FS_AverageHeightFunctional3D<T,Descriptor>* clone() const;
    T getAverageHeight() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        for (pluint i=0; i<modified.size(); ++i) {
            modified[i] = modif::nothing;
        }
    }
private:
    plint averageHeightId;
};


template<typename T, template<typename U> class Descriptor>
T getAverageHeightAtXY(std::vector<MultiBlock3D*> twoPhaseArgs, plint N, Box3D domain);

template<typename T, template<typename U> class Descriptor>
class GetWaterLevelAtxyFunctional3D : public PlainReductiveBoxProcessingFunctional3D
{
public:
    GetWaterLevelAtxyFunctional3D();
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks);
    virtual GetWaterLevelAtxyFunctional3D<T,Descriptor>* clone() const;
    plint getNumFluidCellsAtXY() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        for (pluint i=0; i<modified.size(); ++i) {
            modified[i] = modif::nothing;
        }
    }
private:
    plint numFluidOccupiedCellId;
};

}  // namespace plb

#endif  // FREE_SURFACE_ANALYSIS_3D_H

