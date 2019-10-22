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
 * Helper functions for domain initialization -- header file.
 */
#ifndef DATA_ANALYSIS_FUNCTIONAL_2D_HH
#define DATA_ANALYSIS_FUNCTIONAL_2D_HH

#include "plbWrapper/block/dataAnalysisFunctional2D.h"
#include "plbWrapper/block/plbMath.h"
#include "core/blockStatistics.h"
#include "core/plbDebug.h"
#include "core/util.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "finiteDifference/fdStencils1D.h"
#include "atomicBlock/atomicBlock2D.h"
#include "atomicBlock/blockLattice2D.h"
#include "atomicBlock/dataField2D.h"
#include <cmath>
#include <limits>

namespace plb {

/// Finite Difference operations on data-fields
namespace fdNTensorField {

template<typename T>
inline T bulkXderiv (
        NTensorField2D<T> const& velocity, plint iX, plint iY, int iD )
{
    T dxu = fd::ctl_diff( velocity.get(iX+1,iY)[iD],
                          velocity.get(iX-1,iY)[iD] );
    return dxu;
}

template<typename T>
inline T bulkYderiv (
        NTensorField2D<T> const& velocity, plint iX, plint iY, int iD )
{
    T dyu = fd::ctl_diff( velocity.get(iX,iY+1)[iD],
                          velocity.get(iX,iY-1)[iD] );
    return dyu;
}

template<typename T>
inline T edgeXderiv (
        NTensorField2D<T> const& velocity, int direction, int orientation,
        plint iX, plint iY, int iD )
{
    if (direction==0) {
        return -orientation *
            fd::o1_fwd_diff( velocity.get(iX              ,iY)[iD],
                             velocity.get(iX-1*orientation,iY)[iD] );
    }
    else {
        return bulkXderiv(velocity, iX,iY, iD);
    }
}

template<typename T>
inline T edgeYderiv (
        NTensorField2D<T> const& velocity, int direction, int orientation,
        plint iX, plint iY, int iD )
{
    if (direction==1) {
        return -orientation *
            fd::o1_fwd_diff( velocity.get(iX,iY              )[iD],
                             velocity.get(iX,iY-1*orientation)[iD] );
    }
    else {
        return bulkYderiv(velocity, iX,iY, iD);
    }
}

template<typename T>
inline T cornerXderiv (
        NTensorField2D<T> const& velocity,
        int normalX, int normalY,
        plint iX, plint iY, int iD )
{
    int orientation = normalX;
    return -orientation *
        fd::o1_fwd_diff( velocity.get(iX              ,iY)[iD],
                         velocity.get(iX-1*orientation,iY)[iD] );
}

template<typename T>
inline T cornerYderiv (
        NTensorField2D<T> const& velocity,
        int normalX, int normalY,
        plint iX, plint iY, int iD )
{
    int orientation = normalY;
    return -orientation *
        fd::o1_fwd_diff( velocity.get(iX,iY              )[iD],
                         velocity.get(iX,iY-1*orientation)[iD] );
}

template<typename T>
inline T bulkVorticity(NTensorField2D<T> const& velocity, plint iX, plint iY)
{
    T dxuy = bulkXderiv(velocity, iX,iY, 1);
    T dyux = bulkYderiv(velocity, iX,iY, 0);
    return dxuy - dyux;
}

template<typename T>
inline T edgeVorticity( NTensorField2D<T> const& velocity, int direction, int orientation,
                        plint iX, plint iY )
{
    T dxuy = edgeXderiv(velocity, direction, orientation, iX,iY, 1);
    T dyux = edgeYderiv(velocity, direction, orientation, iX,iY, 0);
    return dxuy - dyux;
}

template<typename T>
inline T cornerVorticity( NTensorField2D<T> const& velocity, int normalX, int normalY,
                          plint iX, plint iY )
{
    T dxuy = cornerXderiv(velocity, normalX, normalY, iX,iY, 1);
    T dyux = cornerYderiv(velocity, normalX, normalY, iX,iY, 0);
    return dxuy - dyux;
}

}  // fdNTensorField


/* *************** Reductive Data Functionals for NTensorField ******** */

template<typename T>
BoxNTensorSumFunctional2D<T>::BoxNTensorSumFunctional2D(plint ndim)
    : sumVectorId(ndim)
{ 
    for (plint iDim=0; iDim<ndim; ++iDim) {
        sumVectorId[iDim] = this->getStatistics().subscribeSum();
    }
}

template<typename T>
void BoxNTensorSumFunctional2D<T>::process (
        Box2D domain, NTensorField2D<T>& vectorField )
{
    PLB_PRECONDITION( vectorField.getNdim() == (plint)sumVectorId.size() );
    plint ndim = vectorField.getNdim();
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iDim=0; iDim<ndim; ++iDim) {
                statistics.gatherSum(sumVectorId[iDim], (double)vectorField.get(iX,iY)[iDim]);
            }
        }
    }
}

template<typename T>
BoxNTensorSumFunctional2D<T>* BoxNTensorSumFunctional2D<T>::clone() const
{
    return new BoxNTensorSumFunctional2D<T>(*this);
}

template<typename T>
void BoxNTensorSumFunctional2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
}

template<typename T>
std::vector<T> BoxNTensorSumFunctional2D<T>::getSumVector() const {
    std::vector<T> sumVector(sumVectorId.size());
    for (pluint iDim=0; iDim<sumVector.size(); ++iDim) {
        double doubleSum = this->getStatistics().getSum(sumVectorId[iDim]);
        // The sum is internally computed on floating-point values. If T is
        //   integer, the value must be rounded at the end.
        if (std::numeric_limits<T>::is_integer) {
            sumVector[iDim] = (T) util::roundToInt(doubleSum);
        }
        else {
            sumVector[iDim] = (T) doubleSum;
        }
    }
    return sumVector;
}


template<typename T>
MaskedBoxNTensorSumFunctional2D<T>::MaskedBoxNTensorSumFunctional2D(plint ndim)
    : sumVectorId(ndim)
{ 
    for (plint iDim=0; iDim<ndim; ++iDim) {
        sumVectorId[iDim] = this->getStatistics().subscribeSum();
    }
}

template<typename T>
void MaskedBoxNTensorSumFunctional2D<T>::process (
        Box2D domain, NTensorField2D<T>& vectorField,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( vectorField.getNdim() == (plint)sumVectorId.size() );
    Dot2D maskOfs = computeRelativeDisplacement(vectorField, mask);
    plint ndim = vectorField.getNdim();
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    statistics.gatherSum(sumVectorId[iDim], (double)vectorField.get(iX,iY)[iDim]);
                }
            }
        }
    }
}

template<typename T>
MaskedBoxNTensorSumFunctional2D<T>* MaskedBoxNTensorSumFunctional2D<T>::clone() const
{
    return new MaskedBoxNTensorSumFunctional2D<T>(*this);
}

template<typename T>
void MaskedBoxNTensorSumFunctional2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
}

template<typename T>
std::vector<T> MaskedBoxNTensorSumFunctional2D<T>::getSumVector() const {
    std::vector<T> sumVector(sumVectorId.size());
    for (pluint iDim=0; iDim<sumVector.size(); ++iDim) {
        double doubleSum = this->getStatistics().getSum(sumVectorId[iDim]);
        // The sum is internally computed on floating-point values. If T is
        //   integer, the value must be rounded at the end.
        if (std::numeric_limits<T>::is_integer) {
            sumVector[iDim] = (T) util::roundToInt(doubleSum);
        }
        else {
            sumVector[iDim] = (T) doubleSum;
        }
    }
    return sumVector;
}


template<typename T>
BoxNTensorMinFunctional2D<T>::BoxNTensorMinFunctional2D(plint ndim)
    : maxVectorId(ndim)
{
    for (plint iDim=0; iDim<ndim; ++iDim) {
        maxVectorId[iDim] = this->getStatistics().subscribeMax();
    }
}

template<typename T>
void BoxNTensorMinFunctional2D<T>::process (
        Box2D domain, NTensorField2D<T>& vectorField )
{
    PLB_PRECONDITION( vectorField.getNdim() == (plint)maxVectorId.size() );
    plint ndim = vectorField.getNdim();
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            // BlockStatistics computes only maximum, no minimum. Therefore,
            //   the relation min(x) = -max(-x) is used.
            for (plint iDim=0; iDim<ndim; ++iDim) {
                statistics.gatherMax(maxVectorId[iDim], -(double)vectorField.get(iX,iY)[iDim]);
            }
        }
    }
}

template<typename T>
BoxNTensorMinFunctional2D<T>* BoxNTensorMinFunctional2D<T>::clone() const
{
    return new BoxNTensorMinFunctional2D<T>(*this);
}

template<typename T>
void BoxNTensorMinFunctional2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;
}

template<typename T>
std::vector<T> BoxNTensorMinFunctional2D<T>::getMinVector() const {
    std::vector<T> minVector(maxVectorId.size());
    for (pluint iDim=0; iDim<minVector.size(); ++iDim) {
        // The minus sign accounts for the relation min(x) = -max(-x).
        double doubleMin = - this->getStatistics().getMax(maxVectorId[iDim]);
        // The min is internally computed on floating-point values. If T is
        //   integer, the value must be rounded at the end.
        if (std::numeric_limits<T>::is_integer) {
            minVector[iDim] = (T) util::roundToInt(doubleMin);
        }
        else {
            minVector[iDim] = (T) doubleMin;
        }
    }
    return minVector;
}


template<typename T>
MaskedBoxNTensorMinFunctional2D<T>::MaskedBoxNTensorMinFunctional2D(plint ndim)
    : maxVectorId(ndim)
{
    for (plint iDim=0; iDim<ndim; ++iDim) {
        maxVectorId[iDim] = this->getStatistics().subscribeMax();
    }
}

template<typename T>
void MaskedBoxNTensorMinFunctional2D<T>::process (
        Box2D domain, NTensorField2D<T>& vectorField,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( vectorField.getNdim() == (plint)maxVectorId.size() );
    Dot2D maskOfs = computeRelativeDisplacement(vectorField, mask);
    plint ndim = vectorField.getNdim();
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                // BlockStatistics computes only maximum, no minimum. Therefore,
                //   the relation min(x) = -max(-x) is used.
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    statistics.gatherMax(maxVectorId[iDim], -(double)vectorField.get(iX,iY)[iDim]);
                }
            }
        }
    }
}

template<typename T>
MaskedBoxNTensorMinFunctional2D<T>* MaskedBoxNTensorMinFunctional2D<T>::clone() const
{
    return new MaskedBoxNTensorMinFunctional2D<T>(*this);
}

template<typename T>
void MaskedBoxNTensorMinFunctional2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
}

template<typename T>
std::vector<T> MaskedBoxNTensorMinFunctional2D<T>::getMinVector() const {
    std::vector<T> minVector(maxVectorId.size());
    for (pluint iDim=0; iDim<minVector.size(); ++iDim) {
        // The minus sign accounts for the relation min(x) = -max(-x).
        double doubleMin = - this->getStatistics().getMax(maxVectorId[iDim]);
        // The min is internally computed on floating-point values. If T is
        //   integer, the value must be rounded at the end.
        if (std::numeric_limits<T>::is_integer) {
            minVector[iDim] = (T) util::roundToInt(doubleMin);
        }
        else {
            minVector[iDim] = (T) doubleMin;
        }
    }
    return minVector;
}


template<typename T>
BoxNTensorMaxFunctional2D<T>::BoxNTensorMaxFunctional2D(plint ndim)
    : maxVectorId(ndim)
{
    for (plint iDim=0; iDim<ndim; ++iDim) {
        maxVectorId[iDim] = this->getStatistics().subscribeMax();
    }
}

template<typename T>
void BoxNTensorMaxFunctional2D<T>::process (
        Box2D domain, NTensorField2D<T>& vectorField )
{
    PLB_PRECONDITION( vectorField.getNdim() == (plint)maxVectorId.size() );
    plint ndim = vectorField.getNdim();
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iDim=0; iDim<ndim; ++iDim) {
                statistics.gatherMax(maxVectorId[iDim], (double)vectorField.get(iX,iY)[iDim]);
            }
        }
    }
}

template<typename T>
BoxNTensorMaxFunctional2D<T>* BoxNTensorMaxFunctional2D<T>::clone() const
{
    return new BoxNTensorMaxFunctional2D<T>(*this);
}

template<typename T>
void BoxNTensorMaxFunctional2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;
}

template<typename T>
std::vector<T> BoxNTensorMaxFunctional2D<T>::getMaxVector() const {
    std::vector<T> maxVector(maxVectorId.size());
    for (pluint iDim=0; iDim<maxVector.size(); ++iDim) {
        double doubleMax = this->getStatistics().getMax(maxVectorId[iDim]);
        // The max is internally computed on floating-point values. If T is
        //   integer, the value must be rounded at the end.
        if (std::numeric_limits<T>::is_integer) {
            maxVector[iDim] = (T) util::roundToInt(doubleMax);
        }
        else {
            maxVector[iDim] = (T) doubleMax;
        }
    }
    return maxVector;
}


template<typename T>
MaskedBoxNTensorMaxFunctional2D<T>::MaskedBoxNTensorMaxFunctional2D(plint ndim)
    : maxVectorId(ndim)
{
    for (plint iDim=0; iDim<ndim; ++iDim) {
        maxVectorId[iDim] = this->getStatistics().subscribeMax();
    }
}

template<typename T>
void MaskedBoxNTensorMaxFunctional2D<T>::process (
        Box2D domain, NTensorField2D<T>& vectorField,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( vectorField.getNdim() == (plint)maxVectorId.size() );
    Dot2D maskOfs = computeRelativeDisplacement(vectorField, mask);
    plint ndim = vectorField.getNdim();
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    statistics.gatherMax(maxVectorId[iDim], (double)vectorField.get(iX,iY)[iDim]);
                }
            }
        }
    }
}

template<typename T>
MaskedBoxNTensorMaxFunctional2D<T>* MaskedBoxNTensorMaxFunctional2D<T>::clone() const
{
    return new MaskedBoxNTensorMaxFunctional2D<T>(*this);
}

template<typename T>
void MaskedBoxNTensorMaxFunctional2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
}

template<typename T>
std::vector<T> MaskedBoxNTensorMaxFunctional2D<T>::getMaxVector() const {
    std::vector<T> maxVector(maxVectorId.size());
    for (pluint iDim=0; iDim<maxVector.size(); ++iDim) {
        double doubleMax = this->getStatistics().getMax(maxVectorId[iDim]);
        // The max is internally computed on floating-point values. If T is
        //   integer, the value must be rounded at the end.
        if (std::numeric_limits<T>::is_integer) {
            maxVector[iDim] = (T) util::roundToInt(doubleMax);
        }
        else {
            maxVector[iDim] = (T) doubleMax;
        }
    }
    return maxVector;
}


template<typename T>
BoundedBoxNTensorSumFunctional2D<T>::BoundedBoxNTensorSumFunctional2D(plint ndim)
    : sumVectorId(ndim)
{
    for (plint iDim=0; iDim<ndim; ++iDim) {
        sumVectorId[iDim] = this->getStatistics().subscribeSum();
    }
}

template<typename T>
void BoundedBoxNTensorSumFunctional2D<T>::processBulk (
        Box2D domain, NTensorField2D<T>& vectorField )
{
    PLB_PRECONDITION( vectorField.getNdim() == (plint)sumVectorId.size() );
    plint ndim = vectorField.getNdim();
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iDim=0; iDim<ndim; ++iDim) {
                statistics.gatherSum(sumVectorId[iDim], (double)vectorField.get(iX,iY)[iDim]);
            }
        }
    }
}

template<typename T>
void BoundedBoxNTensorSumFunctional2D<T>::processEdge (
        int direction, int orientation,
        Box2D domain, NTensorField2D<T>& vectorField )
{
    PLB_PRECONDITION( vectorField.getNdim() == (plint)sumVectorId.size() );
    plint ndim = vectorField.getNdim();
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            // Edge nodes have a weight of 0.5, because only 50% of the
            //   cell centered at the node is inside the computational domain.
            for (plint iDim=0; iDim<ndim; ++iDim) {
                statistics.gatherSum( sumVectorId[iDim],
                                      (double)vectorField.get(iX,iY)[iDim] / 2.);
            }
        }
    }
}

template<typename T>
void BoundedBoxNTensorSumFunctional2D<T>::processCorner (
        int normalX, int normalY,
        Box2D domain, NTensorField2D<T>& vectorField )
{
    PLB_PRECONDITION( vectorField.getNdim() == (plint)sumVectorId.size() );
    plint ndim = vectorField.getNdim();
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            // Corner nodes have a weight of 0.25, because only 25% of the
            //   cell centered at the node is inside the computational domain.
            for (plint iDim=0; iDim<ndim; ++iDim) {
                statistics.gatherSum( sumVectorId[iDim],
                                      (double)vectorField.get(iX,iY)[iDim] / 4.);
            }
        }
    }
}

template<typename T>
BoundedBoxNTensorSumFunctional2D<T>* BoundedBoxNTensorSumFunctional2D<T>::clone() const
{
    return new BoundedBoxNTensorSumFunctional2D<T>(*this);
}

template<typename T>
void BoundedBoxNTensorSumFunctional2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;
}

template<typename T>
std::vector<T> BoundedBoxNTensorSumFunctional2D<T>::getSumVector() const {
    std::vector<T> sumVector(sumVectorId.size());
    for (pluint iDim=0; iDim<sumVector.size(); ++iDim) {
        double doubleSum = this->getStatistics().getSum(sumVectorId[iDim]);
        // The sum is internally computed on floating-point values. If T is
        //   integer, the value must be rounded at the end.
        if (std::numeric_limits<T>::is_integer) {
            sumVector[iDim] = (T) util::roundToInt(doubleSum);
        }
        else {
            sumVector[iDim] = (T) doubleSum;
        }
    }
    return sumVector;
}


template<typename T>
BoundedMaskedBoxNTensorSumFunctional2D<T>::BoundedMaskedBoxNTensorSumFunctional2D(plint ndim)
    : sumVectorId(ndim)
{
    for (plint iDim=0; iDim<ndim; ++iDim) {
        sumVectorId[iDim] = this->getStatistics().subscribeSum();
    }
}

template<typename T>
void BoundedMaskedBoxNTensorSumFunctional2D<T>::processBulk (
        Box2D domain, NTensorField2D<T>& vectorField,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( vectorField.getNdim() == (plint)sumVectorId.size() );
    Dot2D maskOfs = computeRelativeDisplacement(vectorField, mask);
    plint ndim = vectorField.getNdim();
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    statistics.gatherSum(sumVectorId[iDim], (double)vectorField.get(iX,iY)[iDim]);
                }
            }
        }
    }
}

template<typename T>
void BoundedMaskedBoxNTensorSumFunctional2D<T>::processEdge (
        int direction, int orientation,
        Box2D domain, NTensorField2D<T>& vectorField,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( vectorField.getNdim() == (plint)sumVectorId.size() );
    Dot2D maskOfs = computeRelativeDisplacement(vectorField, mask);
    plint ndim = vectorField.getNdim();
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                // Edge nodes have a weight of 0.5, because only 50% of the
                //   cell centered at the node is inside the computational domain.
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    statistics.gatherSum( sumVectorId[iDim],
                                          (double)vectorField.get(iX,iY)[iDim] / 2.);
                }
            }
        }
    }
}

template<typename T>
void BoundedMaskedBoxNTensorSumFunctional2D<T>::processCorner (
        int normalX, int normalY,
        Box2D domain, NTensorField2D<T>& vectorField,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( vectorField.getNdim() == (plint)sumVectorId.size() );
    Dot2D maskOfs = computeRelativeDisplacement(vectorField, mask);
    plint ndim = vectorField.getNdim();
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                // Corner nodes have a weight of 0.25, because only 25% of the
                //   cell centered at the node is inside the computational domain.
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    statistics.gatherSum( sumVectorId[iDim],
                                          (double)vectorField.get(iX,iY)[iDim] / 4.);
                }
            }
        }
    }
}

template<typename T>
BoundedMaskedBoxNTensorSumFunctional2D<T>* BoundedMaskedBoxNTensorSumFunctional2D<T>::clone() const
{
    return new BoundedMaskedBoxNTensorSumFunctional2D<T>(*this);
}

template<typename T>
void BoundedMaskedBoxNTensorSumFunctional2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
}

template<typename T>
std::vector<T> BoundedMaskedBoxNTensorSumFunctional2D<T>::getSumVector() const {
    std::vector<T> sumVector(sumVectorId.size());
    for (pluint iDim=0; iDim<sumVector.size(); ++iDim) {
        double doubleSum = this->getStatistics().getSum(sumVectorId[iDim]);
        // The sum is internally computed on floating-point values. If T is
        //   integer, the value must be rounded at the end.
        if (std::numeric_limits<T>::is_integer) {
            sumVector[iDim] = (T) util::roundToInt(doubleSum);
        }
        else {
            sumVector[iDim] = (T) doubleSum;
        }
    }
    return sumVector;
}


template<typename T1, typename T2>
void MaskedCopyConvertNTensorFunctional2D<T1,T2>::process (
        Box2D domain, NTensorField2D<T1>& field1,
                      NTensorField2D<T2>& field2,
                      NTensorField2D<int>& mask )
{

    PLB_PRECONDITION( field1.getNdim() == field2.getNdim());
    Dot2D maskOfs = computeRelativeDisplacement(field1, mask);
    plint ndim = field1.getNdim();
    Dot2D offset = computeRelativeDisplacement(field1, field2);

    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    for (int iDim=0; iDim<ndim; ++iDim) {
                        *field2.get(iX+offset.x,iY+offset.y) =
                            (T2) *field1.get(iX,iY);
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    for (int iDim=0; iDim<ndim; ++iDim) {
                        field2.get(iX+offset.x,iY+offset.y)[iDim] =
                            (T2) field1.get(iX,iY)[iDim];
                    }
                }
            }
        }
    }
}

template<typename T1, typename T2>
MaskedCopyConvertNTensorFunctional2D<T1,T2>* MaskedCopyConvertNTensorFunctional2D<T1,T2>::clone() const
{
    return new MaskedCopyConvertNTensorFunctional2D<T1,T2>(*this);
}

template<typename T1, typename T2>
void MaskedCopyConvertNTensorFunctional2D<T1,T2>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T1, typename T2>
BlockDomain::DomainT MaskedCopyConvertNTensorFunctional2D<T1,T2>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T>
ExtractNTensorComponentFunctional2D<T>::ExtractNTensorComponentFunctional2D(int iComponent_)
    : iComponent(iComponent_)
{ }

template<typename T>
void ExtractNTensorComponentFunctional2D<T>::process (
        Box2D domain, NTensorField2D<T>& scalarField,
                      NTensorField2D<T>& tensorField )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    PLB_PRECONDITION( tensorField.getNdim()>iComponent );
    Dot2D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            *scalarField.get(iX,iY) = tensorField.get(iX+offset.x,iY+offset.y)[iComponent];
        }
    }
}

template<typename T>
ExtractNTensorComponentFunctional2D<T>* ExtractNTensorComponentFunctional2D<T>::clone() const
{
    return new ExtractNTensorComponentFunctional2D<T>(*this);
}

template<typename T>
void ExtractNTensorComponentFunctional2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT ExtractNTensorComponentFunctional2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T>
MaskedExtractNTensorComponentFunctional2D<T>::MaskedExtractNTensorComponentFunctional2D(int iComponent_)
    : iComponent(iComponent_)
{ }

template<typename T>
void MaskedExtractNTensorComponentFunctional2D<T>::process (
        Box2D domain, NTensorField2D<T>& scalarField,
                      NTensorField2D<T>& tensorField,
                      NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    PLB_PRECONDITION( tensorField.getNdim()>iComponent );
    Dot2D offset = computeRelativeDisplacement(scalarField, tensorField);
    Dot2D maskOfs = computeRelativeDisplacement(scalarField, mask);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                *scalarField.get(iX,iY) = tensorField.get(iX+offset.x,iY+offset.y)[iComponent];
            }
        }
    }
}

template<typename T>
MaskedExtractNTensorComponentFunctional2D<T>* MaskedExtractNTensorComponentFunctional2D<T>::clone() const
{
    return new MaskedExtractNTensorComponentFunctional2D<T>(*this);
}

template<typename T>
void MaskedExtractNTensorComponentFunctional2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT MaskedExtractNTensorComponentFunctional2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T>
void ComputeNTensorNormFunctional2D<T>::process (
        Box2D domain, NTensorField2D<T>& scalarField,
                      NTensorField2D<T>& tensorField )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    plint ndim = tensorField.getNdim();
    Dot2D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iDim=0; iDim<ndim; ++iDim) {
                T normSqr = T();
                T* vector = tensorField.get(iX+offset.x,iY+offset.y);
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    normSqr += vector[iDim]*vector[iDim];
                }
                *scalarField.get(iX,iY) = std::sqrt(normSqr);
            }
        }
    }
}

template<typename T>
ComputeNTensorNormFunctional2D<T>* ComputeNTensorNormFunctional2D<T>::clone() const
{
    return new ComputeNTensorNormFunctional2D<T>(*this);
}

template<typename T>
void ComputeNTensorNormFunctional2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT ComputeNTensorNormFunctional2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T>
void MaskedComputeNTensorNormFunctional2D<T>::process (
        Box2D domain, NTensorField2D<T>& scalarField,
                      NTensorField2D<T>& tensorField,
                      NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    Dot2D maskOfs = computeRelativeDisplacement(scalarField, mask);
    plint ndim = tensorField.getNdim();
    Dot2D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    T normSqr = T();
                    T* vector = tensorField.get(iX+offset.x,iY+offset.y);
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        normSqr += vector[iDim]*vector[iDim];
                    }
                    *scalarField.get(iX,iY) = std::sqrt(normSqr);
                }
            }
        }
    }
}

template<typename T>
MaskedComputeNTensorNormFunctional2D<T>* MaskedComputeNTensorNormFunctional2D<T>::clone() const
{
    return new MaskedComputeNTensorNormFunctional2D<T>(*this);
}

template<typename T>
void MaskedComputeNTensorNormFunctional2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT MaskedComputeNTensorNormFunctional2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T>
void ComputeNTensorNormSqrFunctional2D<T>::process (
        Box2D domain, NTensorField2D<T>& scalarField,
                      NTensorField2D<T>& tensorField )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    plint ndim = tensorField.getNdim();
    Dot2D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iDim=0; iDim<ndim; ++iDim) {
                T normSqr = T();
                T* vector = tensorField.get(iX+offset.x,iY+offset.y);
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    normSqr += vector[iDim]*vector[iDim];
                }
                *scalarField.get(iX,iY) = normSqr;
            }
        }
    }
}

template<typename T>
ComputeNTensorNormSqrFunctional2D<T>* ComputeNTensorNormSqrFunctional2D<T>::clone() const
{
    return new ComputeNTensorNormSqrFunctional2D<T>(*this);
}

template<typename T>
void ComputeNTensorNormSqrFunctional2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT ComputeNTensorNormSqrFunctional2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T>
void MaskedComputeNTensorNormSqrFunctional2D<T>::process (
        Box2D domain, NTensorField2D<T>& scalarField,
                      NTensorField2D<T>& tensorField,
                      NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    Dot2D maskOfs = computeRelativeDisplacement(scalarField, mask);
    plint ndim = tensorField.getNdim();
    Dot2D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    T normSqr = T();
                    T* vector = tensorField.get(iX+offset.x,iY+offset.y);
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        normSqr += vector[iDim]*vector[iDim];
                    }
                    *scalarField.get(iX,iY) = normSqr;
                }
            }
        }
    }
}

template<typename T>
MaskedComputeNTensorNormSqrFunctional2D<T>* MaskedComputeNTensorNormSqrFunctional2D<T>::clone() const
{
    return new MaskedComputeNTensorNormSqrFunctional2D<T>(*this);
}

template<typename T>
void MaskedComputeNTensorNormSqrFunctional2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT MaskedComputeNTensorNormSqrFunctional2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T>
void ComputeSymmetricNTensorNormFunctional2D<T>::process (
        Box2D domain, NTensorField2D<T>& scalarField,
                      NTensorField2D<T>& tensorField )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    typedef SymmetricTensorImpl<T,2> tensor;
    Dot2D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            T* el = tensorField.get(iX+offset.x,iY+offset.y);
            *scalarField.get(iX,iY) = std::sqrt ( 
                    // Count diagonal components once ...
                            util::sqr(el[tensor::xx]) + util::sqr(el[tensor::yy]) +
                    // .. and off-diagonal component twice, due to symmetry.
                    (T)2 * util::sqr(el[tensor::xy]) );
        }
    }
}

template<typename T>
ComputeSymmetricNTensorNormFunctional2D<T>* ComputeSymmetricNTensorNormFunctional2D<T>::clone() const
{
    return new ComputeSymmetricNTensorNormFunctional2D<T>(*this);
}

template<typename T>
void ComputeSymmetricNTensorNormFunctional2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT ComputeSymmetricNTensorNormFunctional2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T>
void MaskedComputeSymmetricNTensorNormFunctional2D<T>::process (
        Box2D domain, NTensorField2D<T>& scalarField,
                      NTensorField2D<T>& tensorField,
                      NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    Dot2D maskOfs = computeRelativeDisplacement(scalarField, mask);
    typedef SymmetricTensorImpl<T,2> tensor;
    Dot2D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                T* el = tensorField.get(iX+offset.x,iY+offset.y);
                *scalarField.get(iX,iY) = std::sqrt ( 
                        // Count diagonal components once ...
                                util::sqr(el[tensor::xx]) + util::sqr(el[tensor::yy]) +
                        // .. and off-diagonal component twice, due to symmetry.
                        (T)2 * util::sqr(el[tensor::xy]) );
            }
        }
    }
}

template<typename T>
MaskedComputeSymmetricNTensorNormFunctional2D<T>* MaskedComputeSymmetricNTensorNormFunctional2D<T>::clone() const
{
    return new MaskedComputeSymmetricNTensorNormFunctional2D<T>(*this);
}

template<typename T>
void MaskedComputeSymmetricNTensorNormFunctional2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT MaskedComputeSymmetricNTensorNormFunctional2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T>
void ComputeSymmetricNTensorNormSqrFunctional2D<T>::process (
        Box2D domain, NTensorField2D<T>& scalarField,
                      NTensorField2D<T>& tensorField )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    typedef SymmetricTensorImpl<T,2> tensor;
    Dot2D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            T* el = tensorField.get(iX+offset.x,iY+offset.y);
            *scalarField.get(iX,iY) = 
                    // Count diagonal components once ...
                            util::sqr(el[tensor::xx]) + util::sqr(el[tensor::yy]) +
                    // .. and off-diagonal components twice, due to symmetry.
                    (T)2 * util::sqr(el[tensor::xy]);
        }
    }
}

template<typename T>
ComputeSymmetricNTensorNormSqrFunctional2D<T>* ComputeSymmetricNTensorNormSqrFunctional2D<T>::clone() const
{
    return new ComputeSymmetricNTensorNormSqrFunctional2D<T>(*this);
}

template<typename T>
void ComputeSymmetricNTensorNormSqrFunctional2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT ComputeSymmetricNTensorNormSqrFunctional2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

template<typename T>
void MaskedComputeSymmetricNTensorNormSqrFunctional2D<T>::process (
        Box2D domain, NTensorField2D<T>& scalarField,
                      NTensorField2D<T>& tensorField,
                      NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    Dot2D maskOfs = computeRelativeDisplacement(scalarField, mask);
    typedef SymmetricTensorImpl<T,2> tensor;
    Dot2D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                T* el = tensorField.get(iX+offset.x,iY+offset.y);
                *scalarField.get(iX,iY) = 
                        // Count diagonal components once ...
                                util::sqr(el[tensor::xx]) + util::sqr(el[tensor::yy]) +
                        // .. and off-diagonal components twice, due to symmetry.
                        (T)2 * util::sqr(el[tensor::xy]);
            }
        }
    }
}

template<typename T>
MaskedComputeSymmetricNTensorNormSqrFunctional2D<T>* MaskedComputeSymmetricNTensorNormSqrFunctional2D<T>::clone() const
{
    return new MaskedComputeSymmetricNTensorNormSqrFunctional2D<T>(*this);
}

template<typename T>
void MaskedComputeSymmetricNTensorNormSqrFunctional2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT MaskedComputeSymmetricNTensorNormSqrFunctional2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

template<typename T>
void ComputeSymmetricNTensorTraceFunctional2D<T>::process (
        Box2D domain, NTensorField2D<T>& scalarField,
                      NTensorField2D<T>& tensorField )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    typedef SymmetricTensorImpl<T,2> tensor;
    Dot2D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            T* el = tensorField.get(iX+offset.x,iY+offset.y);
            *scalarField.get(iX,iY) = el[tensor::xx] + el[tensor::yy];
        }
    }
}

template<typename T>
ComputeSymmetricNTensorTraceFunctional2D<T>* ComputeSymmetricNTensorTraceFunctional2D<T>::clone() const
{
    return new ComputeSymmetricNTensorTraceFunctional2D<T>(*this);
}

template<typename T>
void ComputeSymmetricNTensorTraceFunctional2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT ComputeSymmetricNTensorTraceFunctional2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T>
void MaskedComputeSymmetricNTensorTraceFunctional2D<T>::process (
        Box2D domain, NTensorField2D<T>& scalarField,
                      NTensorField2D<T>& tensorField,
                      NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    Dot2D maskOfs = computeRelativeDisplacement(scalarField, mask);
    typedef SymmetricTensorImpl<T,2> tensor;
    Dot2D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                T* el = tensorField.get(iX+offset.x,iY+offset.y);
                *scalarField.get(iX,iY) = el[tensor::xx] + el[tensor::yy];
            }
        }
    }
}

template<typename T>
MaskedComputeSymmetricNTensorTraceFunctional2D<T>* MaskedComputeSymmetricNTensorTraceFunctional2D<T>::clone() const
{
    return new MaskedComputeSymmetricNTensorTraceFunctional2D<T>(*this);
}

template<typename T>
void MaskedComputeSymmetricNTensorTraceFunctional2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT MaskedComputeSymmetricNTensorTraceFunctional2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T>
void BoxBulkNTensorVorticityFunctional2D<T>::process (
        Box2D domain, NTensorField2D<T>& vorticity,
                      NTensorField2D<T>& velocity )
{
    PLB_PRECONDITION( vorticity.getNdim()==1 );
    PLB_PRECONDITION( velocity.getNdim()==2 );
    Dot2D offset = computeRelativeDisplacement(vorticity, velocity);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            plint iX2 = iX+offset.x;
            plint iY2 = iY+offset.y;
            *vorticity.get(iX,iY) =
                fdNTensorField::bulkVorticity(velocity, iX2,iY2);
        }
    }
}

template<typename T>
BoxBulkNTensorVorticityFunctional2D<T>* BoxBulkNTensorVorticityFunctional2D<T>::clone() const
{
    return new BoxBulkNTensorVorticityFunctional2D<T>(*this);
}

template<typename T>
void BoxBulkNTensorVorticityFunctional2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT BoxBulkNTensorVorticityFunctional2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}



template<typename T>
void MaskedBoxBulkNTensorVorticityFunctional2D<T>::process (
        Box2D domain, NTensorField2D<T>& vorticity,
                      NTensorField2D<T>& velocity,
                      NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( vorticity.getNdim()==1 );
    PLB_PRECONDITION( velocity.getNdim()==2 );
    Dot2D maskOfs = computeRelativeDisplacement(vorticity, mask);
    Dot2D offset = computeRelativeDisplacement(vorticity, velocity);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                *vorticity.get(iX,iY) =
                    fdNTensorField::bulkVorticity(velocity, iX2,iY2);
            }
        }
    }
}

template<typename T>
MaskedBoxBulkNTensorVorticityFunctional2D<T>* MaskedBoxBulkNTensorVorticityFunctional2D<T>::clone() const
{
    return new MaskedBoxBulkNTensorVorticityFunctional2D<T>(*this);
}

template<typename T>
void MaskedBoxBulkNTensorVorticityFunctional2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT MaskedBoxBulkNTensorVorticityFunctional2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}



template<typename T>
void BoxNTensorVorticityFunctional2D<T>::processBulk (
        Box2D domain, NTensorField2D<T>& vorticity,
                      NTensorField2D<T>& velocity )
{
    PLB_PRECONDITION( vorticity.getNdim()==1 );
    PLB_PRECONDITION( velocity.getNdim()==2 );
    Dot2D offset = computeRelativeDisplacement(vorticity, velocity);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            plint iX2 = iX+offset.x;
            plint iY2 = iY+offset.y;
            *vorticity.get(iX,iY) =
                 fdNTensorField::bulkVorticity(velocity, iX2,iY2);
        }
    }
}

template<typename T>
void BoxNTensorVorticityFunctional2D<T>::processEdge (
        int direction, int orientation, Box2D domain,
        NTensorField2D<T>& vorticity, NTensorField2D<T>& velocity )
{
    PLB_PRECONDITION( vorticity.getNdim()==1 );
    PLB_PRECONDITION( velocity.getNdim()==2 );
    Dot2D offset = computeRelativeDisplacement(vorticity, velocity);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            plint iX2 = iX+offset.x;
            plint iY2 = iY+offset.y;
            *vorticity.get(iX,iY) =
                 fdNTensorField::edgeVorticity(velocity,direction,orientation, iX2,iY2);
        }
    }
}

template<typename T>
void BoxNTensorVorticityFunctional2D<T>::processCorner (
        int normalX, int normalY,  Box2D domain,
        NTensorField2D<T>& vorticity, NTensorField2D<T>& velocity )
{
    PLB_PRECONDITION( vorticity.getNdim()==1 );
    PLB_PRECONDITION( velocity.getNdim()==2 );
    Dot2D offset = computeRelativeDisplacement(vorticity, velocity);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            plint iX2 = iX+offset.x;
            plint iY2 = iY+offset.y;
            *vorticity.get(iX,iY) =
                 fdNTensorField::cornerVorticity(velocity,normalX,normalY, iX2,iY2);
        }
    }
}

template<typename T>
BoxNTensorVorticityFunctional2D<T>* BoxNTensorVorticityFunctional2D<T>::clone() const
{
    return new BoxNTensorVorticityFunctional2D<T>(*this);
}

template<typename T>
void BoxNTensorVorticityFunctional2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}


template<typename T>
BlockDomain::DomainT BoxNTensorVorticityFunctional2D<T>::appliesTo() const {
    // Don't apply to envelope, because nearest neighbors need to be accessed.
    return BlockDomain::bulk;
}


template<typename T>
void MaskedBoxNTensorVorticityFunctional2D<T>::processBulk (
        Box2D domain, NTensorField2D<T>& vorticity,
                      NTensorField2D<T>& velocity,
                      NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( vorticity.getNdim()==1 );
    PLB_PRECONDITION( velocity.getNdim()==2 );
    Dot2D maskOfs = computeRelativeDisplacement(vorticity, mask);
    Dot2D offset = computeRelativeDisplacement(vorticity, velocity);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                *vorticity.get(iX,iY) =
                     fdNTensorField::bulkVorticity(velocity, iX2,iY2);
            }
        }
    }
}

template<typename T>
void MaskedBoxNTensorVorticityFunctional2D<T>::processEdge (
        int direction, int orientation, Box2D domain,
        NTensorField2D<T>& vorticity, NTensorField2D<T>& velocity,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( vorticity.getNdim()==1 );
    PLB_PRECONDITION( velocity.getNdim()==2 );
    Dot2D maskOfs = computeRelativeDisplacement(vorticity, mask);
    Dot2D offset = computeRelativeDisplacement(vorticity, velocity);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                *vorticity.get(iX,iY) =
                     fdNTensorField::edgeVorticity(velocity,direction,orientation, iX2,iY2);
            }
        }
    }
}

template<typename T>
void MaskedBoxNTensorVorticityFunctional2D<T>::processCorner (
        int normalX, int normalY,  Box2D domain,
        NTensorField2D<T>& vorticity, NTensorField2D<T>& velocity,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( vorticity.getNdim()==1 );
    PLB_PRECONDITION( velocity.getNdim()==2 );
    Dot2D maskOfs = computeRelativeDisplacement(vorticity, mask);
    Dot2D offset = computeRelativeDisplacement(vorticity, velocity);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                *vorticity.get(iX,iY) =
                     fdNTensorField::cornerVorticity(velocity,normalX,normalY, iX2,iY2);
            }
        }
    }
}

template<typename T>
MaskedBoxNTensorVorticityFunctional2D<T>* MaskedBoxNTensorVorticityFunctional2D<T>::clone() const
{
    return new MaskedBoxNTensorVorticityFunctional2D<T>(*this);
}

template<typename T>
void MaskedBoxNTensorVorticityFunctional2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
}


template<typename T>
BlockDomain::DomainT MaskedBoxNTensorVorticityFunctional2D<T>::appliesTo() const {
    // Don't apply to envelope, because nearest neighbors need to be accessed.
    return BlockDomain::bulk;
}


template<typename T>
void BoxBulkNTensorStrainRateFunctional2D<T>::process (
        Box2D domain, NTensorField2D<T>& velocity,
                      NTensorField2D<T>& S )
{
    PLB_PRECONDITION( velocity.getNdim()==2 );
    PLB_PRECONDITION( S.getNdim()==3 );
    typedef SymmetricTensorImpl<T,2> tensor;
    Dot2D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            plint iX2 = iX+offset.x;
            plint iY2 = iY+offset.y;
            T* el = S.get(iX2,iY2);
            el[tensor::xx] = fdNTensorField::bulkXderiv(velocity, iX, iY, 0);
            el[tensor::yy] = fdNTensorField::bulkYderiv(velocity, iX, iY, 1);
            el[tensor::xy] = ( fdNTensorField::bulkXderiv(velocity, iX, iY, 1) +
                               fdNTensorField::bulkYderiv(velocity, iX, iY, 0) ) / (T)2;
        }
    }
}

template<typename T>
BoxBulkNTensorStrainRateFunctional2D<T>* BoxBulkNTensorStrainRateFunctional2D<T>::clone() const
{
    return new BoxBulkNTensorStrainRateFunctional2D<T>(*this);
}

template<typename T>
void BoxBulkNTensorStrainRateFunctional2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT BoxBulkNTensorStrainRateFunctional2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}



template<typename T>
void MaskedBoxBulkNTensorStrainRateFunctional2D<T>::process (
        Box2D domain, NTensorField2D<T>& velocity,
                      NTensorField2D<T>& S,
                      NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( velocity.getNdim()==2 );
    PLB_PRECONDITION( S.getNdim()==3 );
    Dot2D maskOfs = computeRelativeDisplacement(velocity, mask);
    typedef SymmetricTensorImpl<T,2> tensor;
    Dot2D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                T* el = S.get(iX2,iY2);
                el[tensor::xx] = fdNTensorField::bulkXderiv(velocity, iX, iY, 0);
                el[tensor::yy] = fdNTensorField::bulkYderiv(velocity, iX, iY, 1);
                el[tensor::xy] = ( fdNTensorField::bulkXderiv(velocity, iX, iY, 1) +
                                   fdNTensorField::bulkYderiv(velocity, iX, iY, 0) ) / (T)2;
            }
        }
    }
}

template<typename T>
MaskedBoxBulkNTensorStrainRateFunctional2D<T>* MaskedBoxBulkNTensorStrainRateFunctional2D<T>::clone() const
{
    return new MaskedBoxBulkNTensorStrainRateFunctional2D<T>(*this);
}

template<typename T>
void MaskedBoxBulkNTensorStrainRateFunctional2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT MaskedBoxBulkNTensorStrainRateFunctional2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T>
void BoxNTensorStrainRateFunctional2D<T>::processBulk (
        Box2D domain, NTensorField2D<T>& velocity,
                      NTensorField2D<T>& S )
{
    typedef SymmetricTensorImpl<T,2> tensor;
    PLB_PRECONDITION( velocity.getNdim()==2 );
    PLB_PRECONDITION( S.getNdim()==3 );
    Dot2D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            plint iX2 = iX+offset.x;
            plint iY2 = iY+offset.y;
            T* el = S.get(iX2,iY2);
            el[tensor::xx] = fdNTensorField::bulkXderiv(velocity, iX, iY, 0);
            el[tensor::yy] = fdNTensorField::bulkYderiv(velocity, iX, iY, 1);
            el[tensor::xy] = ( fdNTensorField::bulkXderiv(velocity, iX, iY, 1) +
                               fdNTensorField::bulkYderiv(velocity, iX, iY, 0) ) / (T)2;
        }
    }
}

template<typename T>
void BoxNTensorStrainRateFunctional2D<T>::processEdge (
        int direction, int orientation, Box2D domain,
        NTensorField2D<T>& velocity, NTensorField2D<T>& S )
{
    PLB_PRECONDITION( velocity.getNdim()==2 );
    PLB_PRECONDITION( S.getNdim()==3 );
    typedef SymmetricTensorImpl<T,2> tensor;
    Dot2D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            plint iX2 = iX+offset.x;
            plint iY2 = iY+offset.y;
            T* el = S.get(iX2,iY2);
            el[tensor::xx] = fdNTensorField::edgeXderiv(velocity, direction,orientation, iX, iY, 0);
            el[tensor::yy] = fdNTensorField::edgeYderiv(velocity, direction,orientation, iX, iY, 1);
            el[tensor::xy] = ( fdNTensorField::edgeXderiv(velocity, direction,orientation, iX, iY, 1) +
                               fdNTensorField::edgeYderiv(velocity, direction,orientation, iX, iY, 0) ) / (T)2;
        }
    }
}

template<typename T>
void BoxNTensorStrainRateFunctional2D<T>::processCorner (
        int normalX, int normalY,  Box2D domain,
        NTensorField2D<T>& velocity, NTensorField2D<T>& S )
{
    PLB_PRECONDITION( velocity.getNdim()==2 );
    PLB_PRECONDITION( S.getNdim()==3 );
    typedef SymmetricTensorImpl<T,2> tensor;
    Dot2D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            plint iX2 = iX+offset.x;
            plint iY2 = iY+offset.y;
            T* el = S.get(iX2,iY2);
            el[tensor::xx] = fdNTensorField::cornerXderiv(velocity, normalX,normalY, iX, iY, 0);
            el[tensor::yy] = fdNTensorField::cornerYderiv(velocity, normalX,normalY, iX, iY, 1);
            el[tensor::xy] = ( fdNTensorField::cornerXderiv(velocity, normalX,normalY, iX, iY, 1) +
                               fdNTensorField::cornerYderiv(velocity, normalX,normalY, iX, iY, 0) ) / (T)2;
        }
    }
}

template<typename T>
BoxNTensorStrainRateFunctional2D<T>* BoxNTensorStrainRateFunctional2D<T>::clone() const
{
    return new BoxNTensorStrainRateFunctional2D<T>(*this);
}

template<typename T>
void BoxNTensorStrainRateFunctional2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}


template<typename T>
BlockDomain::DomainT BoxNTensorStrainRateFunctional2D<T>::appliesTo() const {
    // Don't apply to envelope, because nearest neighbors need to be accessed.
    return BlockDomain::bulk;
}


template<typename T>
void MaskedBoxNTensorStrainRateFunctional2D<T>::processBulk (
        Box2D domain, NTensorField2D<T>& velocity,
                      NTensorField2D<T>& S,
                      NTensorField2D<int>& mask )
{
    typedef SymmetricTensorImpl<T,2> tensor;
    PLB_PRECONDITION( velocity.getNdim()==2 );
    PLB_PRECONDITION( S.getNdim()==3 );
    Dot2D maskOfs = computeRelativeDisplacement(velocity, mask);
    Dot2D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                T* el = S.get(iX2,iY2);
                el[tensor::xx] = fdNTensorField::bulkXderiv(velocity, iX, iY, 0);
                el[tensor::yy] = fdNTensorField::bulkYderiv(velocity, iX, iY, 1);
                el[tensor::xy] = ( fdNTensorField::bulkXderiv(velocity, iX, iY, 1) +
                                   fdNTensorField::bulkYderiv(velocity, iX, iY, 0) ) / (T)2;
            }
        }
    }
}

template<typename T>
void MaskedBoxNTensorStrainRateFunctional2D<T>::processEdge (
        int direction, int orientation, Box2D domain,
        NTensorField2D<T>& velocity, NTensorField2D<T>& S,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( velocity.getNdim()==2 );
    PLB_PRECONDITION( S.getNdim()==3 );
    Dot2D maskOfs = computeRelativeDisplacement(velocity, mask);
    typedef SymmetricTensorImpl<T,2> tensor;
    Dot2D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                T* el = S.get(iX2,iY2);
                el[tensor::xx] = fdNTensorField::edgeXderiv(velocity, direction,orientation, iX, iY, 0);
                el[tensor::yy] = fdNTensorField::edgeYderiv(velocity, direction,orientation, iX, iY, 1);
                el[tensor::xy] = ( fdNTensorField::edgeXderiv(velocity, direction,orientation, iX, iY, 1) +
                                   fdNTensorField::edgeYderiv(velocity, direction,orientation, iX, iY, 0) ) / (T)2;
            }
        }
    }
}

template<typename T>
void MaskedBoxNTensorStrainRateFunctional2D<T>::processCorner (
        int normalX, int normalY,  Box2D domain,
        NTensorField2D<T>& velocity, NTensorField2D<T>& S,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( velocity.getNdim()==2 );
    PLB_PRECONDITION( S.getNdim()==3 );
    Dot2D maskOfs = computeRelativeDisplacement(velocity, mask);
    typedef SymmetricTensorImpl<T,2> tensor;
    Dot2D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                T* el = S.get(iX2,iY2);
                el[tensor::xx] = fdNTensorField::cornerXderiv(velocity, normalX,normalY, iX, iY, 0);
                el[tensor::yy] = fdNTensorField::cornerYderiv(velocity, normalX,normalY, iX, iY, 1);
                el[tensor::xy] = ( fdNTensorField::cornerXderiv(velocity, normalX,normalY, iX, iY, 1) +
                                   fdNTensorField::cornerYderiv(velocity, normalX,normalY, iX, iY, 0) ) / (T)2;
            }
        }
    }
}

template<typename T>
MaskedBoxNTensorStrainRateFunctional2D<T>* MaskedBoxNTensorStrainRateFunctional2D<T>::clone() const
{
    return new MaskedBoxNTensorStrainRateFunctional2D<T>(*this);
}

template<typename T>
void MaskedBoxNTensorStrainRateFunctional2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}


template<typename T>
BlockDomain::DomainT MaskedBoxNTensorStrainRateFunctional2D<T>::appliesTo() const {
    // Don't apply to envelope, because nearest neighbors need to be accessed.
    return BlockDomain::bulk;
}


/* ******** A_plus_alpha_NTensor2D ************************************* */

template<typename T>
A_plus_alpha_NTensor2D<T>::A_plus_alpha_NTensor2D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_plus_alpha_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<T>& result )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot2D offset = computeRelativeDisplacement(A, result);

    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                *result.get(iX+offset.x,iY+offset.y)
                    = *A.get(iX,iY) + scalarAlpha;
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    result.get(iX+offset.x,iY+offset.y)[iDim]
                        = A.get(iX,iY)[iDim] + alpha[iDim];
                }
            }
        }
    }
}

template<typename T>
A_plus_alpha_NTensor2D<T>* A_plus_alpha_NTensor2D<T>::clone() const {
    return new A_plus_alpha_NTensor2D<T>(*this);
}

template<typename T>
void A_plus_alpha_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_plus_alpha_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_plus_alpha_NTensor2D ************************************* */

template<typename T>
Masked_A_plus_alpha_NTensor2D<T>::Masked_A_plus_alpha_NTensor2D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void Masked_A_plus_alpha_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<T>& result,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot2D maskOfs = computeRelativeDisplacement(A, mask);
    Dot2D offset = computeRelativeDisplacement(A, result);

    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    *result.get(iX+offset.x,iY+offset.y)
                        = *A.get(iX,iY) + scalarAlpha;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offset.x,iY+offset.y)[iDim]
                            = A.get(iX,iY)[iDim] + alpha[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_plus_alpha_NTensor2D<T>* Masked_A_plus_alpha_NTensor2D<T>::clone() const {
    return new Masked_A_plus_alpha_NTensor2D<T>(*this);
}

template<typename T>
void Masked_A_plus_alpha_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_plus_alpha_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** A_minus_alpha_NTensor2D ************************************** */

template<typename T>
A_minus_alpha_NTensor2D<T>::A_minus_alpha_NTensor2D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_minus_alpha_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<T>& result )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot2D offset = computeRelativeDisplacement(A, result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                *result.get(iX+offset.x,iY+offset.y)
                    = *A.get(iX,iY) - scalarAlpha;
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    result.get(iX+offset.x,iY+offset.y)[iDim]
                        = A.get(iX,iY)[iDim] - alpha[iDim];
                }
            }
        }
    }
}

template<typename T>
A_minus_alpha_NTensor2D<T>* A_minus_alpha_NTensor2D<T>::clone() const {
    return new A_minus_alpha_NTensor2D<T>(*this);
}

template<typename T>
void A_minus_alpha_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_minus_alpha_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}



/* ******** Masked_A_minus_alpha_NTensor2D ************************************** */

template<typename T>
Masked_A_minus_alpha_NTensor2D<T>::Masked_A_minus_alpha_NTensor2D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void Masked_A_minus_alpha_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<T>& result,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot2D maskOfs = computeRelativeDisplacement(A, mask);
    Dot2D offset = computeRelativeDisplacement(A, result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    *result.get(iX+offset.x,iY+offset.y)
                        = *A.get(iX,iY) - scalarAlpha;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offset.x,iY+offset.y)[iDim]
                            = A.get(iX,iY)[iDim] - alpha[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_minus_alpha_NTensor2D<T>* Masked_A_minus_alpha_NTensor2D<T>::clone() const {
    return new Masked_A_minus_alpha_NTensor2D<T>(*this);
}

template<typename T>
void Masked_A_minus_alpha_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_minus_alpha_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}



/* ******** Alpha_minus_A_NTensor2D ************************************* */

template<typename T>
Alpha_minus_A_NTensor2D<T>::Alpha_minus_A_NTensor2D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void Alpha_minus_A_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<T>& result )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot2D offset = computeRelativeDisplacement(A, result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                *result.get(iX+offset.x,iY+offset.y)
                    = scalarAlpha - *A.get(iX,iY);
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    result.get(iX+offset.x,iY+offset.y)[iDim]
                        = alpha[iDim] - A.get(iX,iY)[iDim];
                }
            }
        }
    }
}

template<typename T>
Alpha_minus_A_NTensor2D<T>* Alpha_minus_A_NTensor2D<T>::clone() const {
    return new Alpha_minus_A_NTensor2D<T>(*this);
}

template<typename T>
void Alpha_minus_A_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT Alpha_minus_A_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_Alpha_minus_A_NTensor2D ************************************* */

template<typename T>
Masked_Alpha_minus_A_NTensor2D<T>::Masked_Alpha_minus_A_NTensor2D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void Masked_Alpha_minus_A_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<T>& result,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot2D maskOfs = computeRelativeDisplacement(A, mask);
    Dot2D offset = computeRelativeDisplacement(A, result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    *result.get(iX+offset.x,iY+offset.y)
                        = scalarAlpha - *A.get(iX,iY);
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offset.x,iY+offset.y)[iDim]
                            = alpha[iDim] - A.get(iX,iY)[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_Alpha_minus_A_NTensor2D<T>* Masked_Alpha_minus_A_NTensor2D<T>::clone() const {
    return new Masked_Alpha_minus_A_NTensor2D<T>(*this);
}

template<typename T>
void Masked_Alpha_minus_A_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_Alpha_minus_A_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** A_times_alpha_NTensor2D ************************************* */

template<typename T>
A_times_alpha_NTensor2D<T>::A_times_alpha_NTensor2D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_times_alpha_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<T>& result )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot2D offset = computeRelativeDisplacement(A, result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                *result.get(iX+offset.x,iY+offset.y)
                    = *A.get(iX,iY) * scalarAlpha;
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    result.get(iX+offset.x,iY+offset.y)[iDim]
                        = A.get(iX,iY)[iDim] * alpha[iDim];
                }
            }
        }
    }
}

template<typename T>
A_times_alpha_NTensor2D<T>* A_times_alpha_NTensor2D<T>::clone() const {
    return new A_times_alpha_NTensor2D<T>(*this);
}

template<typename T>
void A_times_alpha_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_times_alpha_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_times_alpha_NTensor2D ************************************* */

template<typename T>
Masked_A_times_alpha_NTensor2D<T>::Masked_A_times_alpha_NTensor2D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void Masked_A_times_alpha_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<T>& result,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot2D maskOfs = computeRelativeDisplacement(A, mask);
    Dot2D offset = computeRelativeDisplacement(A, result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    *result.get(iX+offset.x,iY+offset.y)
                        = *A.get(iX,iY) * scalarAlpha;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offset.x,iY+offset.y)[iDim]
                            = A.get(iX,iY)[iDim] * alpha[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_times_alpha_NTensor2D<T>* Masked_A_times_alpha_NTensor2D<T>::clone() const {
    return new Masked_A_times_alpha_NTensor2D<T>(*this);
}

template<typename T>
void Masked_A_times_alpha_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_times_alpha_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** A_dividedBy_alpha_NTensor2D ************************************* */

template<typename T>
A_dividedBy_alpha_NTensor2D<T>::A_dividedBy_alpha_NTensor2D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_dividedBy_alpha_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<T>& result )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot2D offset = computeRelativeDisplacement(A, result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                *result.get(iX+offset.x,iY+offset.y)
                    = *A.get(iX,iY) / scalarAlpha;
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    result.get(iX+offset.x,iY+offset.y)[iDim]
                        = A.get(iX,iY)[iDim] / alpha[iDim];
                }
            }
        }
    }
}

template<typename T>
A_dividedBy_alpha_NTensor2D<T>* A_dividedBy_alpha_NTensor2D<T>::clone() const {
    return new A_dividedBy_alpha_NTensor2D<T>(*this);
}

template<typename T>
void A_dividedBy_alpha_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_dividedBy_alpha_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_dividedBy_alpha_NTensor2D ************************************* */

template<typename T>
Masked_A_dividedBy_alpha_NTensor2D<T>::Masked_A_dividedBy_alpha_NTensor2D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void Masked_A_dividedBy_alpha_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<T>& result,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot2D maskOfs = computeRelativeDisplacement(A, mask);
    Dot2D offset = computeRelativeDisplacement(A, result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    *result.get(iX+offset.x,iY+offset.y)
                        = *A.get(iX,iY) / scalarAlpha;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offset.x,iY+offset.y)[iDim]
                            = A.get(iX,iY)[iDim] / alpha[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_dividedBy_alpha_NTensor2D<T>* Masked_A_dividedBy_alpha_NTensor2D<T>::clone() const {
    return new Masked_A_dividedBy_alpha_NTensor2D<T>(*this);
}

template<typename T>
void Masked_A_dividedBy_alpha_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_dividedBy_alpha_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Alpha_dividedBy_A_NTensor2D ************************************* */

template<typename T>
Alpha_dividedBy_A_NTensor2D<T>::Alpha_dividedBy_A_NTensor2D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void Alpha_dividedBy_A_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<T>& result )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot2D offset = computeRelativeDisplacement(A, result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                *result.get(iX+offset.x,iY+offset.y)
                    = scalarAlpha / *A.get(iX,iY);
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    result.get(iX+offset.x,iY+offset.y)[iDim]
                        = alpha[iDim] / A.get(iX,iY)[iDim];
                }
            }
        }
    }
}

template<typename T>
Alpha_dividedBy_A_NTensor2D<T>* Alpha_dividedBy_A_NTensor2D<T>::clone() const {
    return new Alpha_dividedBy_A_NTensor2D<T>(*this);
}

template<typename T>
void Alpha_dividedBy_A_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT Alpha_dividedBy_A_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}



/* ******** Masked_Alpha_dividedBy_A_NTensor2D ************************************* */

template<typename T>
Masked_Alpha_dividedBy_A_NTensor2D<T>::Masked_Alpha_dividedBy_A_NTensor2D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void Masked_Alpha_dividedBy_A_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<T>& result,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot2D maskOfs = computeRelativeDisplacement(A, mask);
    Dot2D offset = computeRelativeDisplacement(A, result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    *result.get(iX+offset.x,iY+offset.y)
                        = scalarAlpha / *A.get(iX,iY);
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offset.x,iY+offset.y)[iDim]
                            = alpha[iDim] / A.get(iX,iY)[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_Alpha_dividedBy_A_NTensor2D<T>* Masked_Alpha_dividedBy_A_NTensor2D<T>::clone() const {
    return new Masked_Alpha_dividedBy_A_NTensor2D<T>(*this);
}

template<typename T>
void Masked_Alpha_dividedBy_A_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_Alpha_dividedBy_A_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}



/* ******** A_toThePower_alpha_NTensor2D ************************************* */

template<typename T>
A_toThePower_alpha_NTensor2D<T>::A_toThePower_alpha_NTensor2D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_toThePower_alpha_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<T>& result )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot2D offset = computeRelativeDisplacement(A, result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                *result.get(iX+offset.x,iY+offset.y)
                    = customPower(*A.get(iX,iY),scalarAlpha);
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    result.get(iX+offset.x,iY+offset.y)[iDim]
                        = customPower(A.get(iX,iY)[iDim],alpha[iDim]);
                }
            }
        }
    }
}

template<typename T>
A_toThePower_alpha_NTensor2D<T>* A_toThePower_alpha_NTensor2D<T>::clone() const {
    return new A_toThePower_alpha_NTensor2D<T>(*this);
}

template<typename T>
void A_toThePower_alpha_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_toThePower_alpha_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_toThePower_alpha_NTensor2D ************************************* */

template<typename T>
Masked_A_toThePower_alpha_NTensor2D<T>::Masked_A_toThePower_alpha_NTensor2D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void Masked_A_toThePower_alpha_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<T>& result,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot2D maskOfs = computeRelativeDisplacement(A, mask);
    Dot2D offset = computeRelativeDisplacement(A, result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    *result.get(iX+offset.x,iY+offset.y)
                        = customPower(*A.get(iX,iY),scalarAlpha);
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offset.x,iY+offset.y)[iDim]
                            = customPower(A.get(iX,iY)[iDim],alpha[iDim]);
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_toThePower_alpha_NTensor2D<T>* Masked_A_toThePower_alpha_NTensor2D<T>::clone() const {
    return new Masked_A_toThePower_alpha_NTensor2D<T>(*this);
}

template<typename T>
void Masked_A_toThePower_alpha_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_toThePower_alpha_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Alpha_toThePower_A_NTensor2D ************************************* */

template<typename T>
Alpha_toThePower_A_NTensor2D<T>::Alpha_toThePower_A_NTensor2D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void Alpha_toThePower_A_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<T>& result )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot2D offset = computeRelativeDisplacement(A, result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                *result.get(iX+offset.x,iY+offset.y)
                    = customPower(scalarAlpha,*A.get(iX,iY));
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    result.get(iX+offset.x,iY+offset.y)[iDim]
                        = customPower(alpha[iDim],A.get(iX,iY)[iDim]);
                }
            }
        }
    }
}

template<typename T>
Alpha_toThePower_A_NTensor2D<T>* Alpha_toThePower_A_NTensor2D<T>::clone() const {
    return new Alpha_toThePower_A_NTensor2D<T>(*this);
}

template<typename T>
void Alpha_toThePower_A_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT Alpha_toThePower_A_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}



/* ******** Masked_Alpha_toThePower_A_NTensor2D ************************************* */

template<typename T>
Masked_Alpha_toThePower_A_NTensor2D<T>::Masked_Alpha_toThePower_A_NTensor2D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void Masked_Alpha_toThePower_A_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<T>& result,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot2D maskOfs = computeRelativeDisplacement(A, mask);
    Dot2D offset = computeRelativeDisplacement(A, result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    *result.get(iX+offset.x,iY+offset.y)
                        = customPower(scalarAlpha,*A.get(iX,iY));
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offset.x,iY+offset.y)[iDim]
                            = customPower(alpha[iDim],A.get(iX,iY)[iDim]);
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_Alpha_toThePower_A_NTensor2D<T>* Masked_Alpha_toThePower_A_NTensor2D<T>::clone() const {
    return new Masked_Alpha_toThePower_A_NTensor2D<T>(*this);
}

template<typename T>
void Masked_Alpha_toThePower_A_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_Alpha_toThePower_A_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}



/* ******** A_plus_B_NTensor2D ************************************ */

template<typename T>
void A_plus_B_NTensor2D<T>::process (
        Box2D domain, std::vector<NTensorField2D<T>*> fields )
{
    PLB_PRECONDITION( fields.size()==3 );
    NTensorField2D<T>& A = *fields[0];
    NTensorField2D<T>& B = *fields[1];
    NTensorField2D<T>& result = *fields[2];
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot2D offsetB      = computeRelativeDisplacement(A,B);
    Dot2D offsetResult = computeRelativeDisplacement(A,result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                *result.get(iX+offsetResult.x,iY+offsetResult.y)
                    = *A.get(iX,iY) + *B.get(iX+offsetB.x,iY+offsetB.y);
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    result.get(iX+offsetResult.x,iY+offsetResult.y)[iDim]
                        = A.get(iX,iY)[iDim] + B.get(iX+offsetB.x,iY+offsetB.y)[iDim];
                }
            }
        }
    }
}

template<typename T>
A_plus_B_NTensor2D<T>* A_plus_B_NTensor2D<T>::clone() const {
    return new A_plus_B_NTensor2D<T>(*this);
}

template<typename T>
void A_plus_B_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_plus_B_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_plus_B_NTensor2D ************************************ */

template<typename T>
void Masked_A_plus_B_NTensor2D<T>::process (
        Box2D domain, std::vector<NTensorField2D<T>*> fields,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( fields.size()==3 );
    NTensorField2D<T>& A = *fields[0];
    NTensorField2D<T>& B = *fields[1];
    NTensorField2D<T>& result = *fields[2];
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot2D offsetB      = computeRelativeDisplacement(A,B);
    Dot2D offsetResult = computeRelativeDisplacement(A,result);
    Dot2D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    *result.get(iX+offsetResult.x,iY+offsetResult.y)
                        = *A.get(iX,iY) + *B.get(iX+offsetB.x,iY+offsetB.y);
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offsetResult.x,iY+offsetResult.y)[iDim]
                            = A.get(iX,iY)[iDim] + B.get(iX+offsetB.x,iY+offsetB.y)[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_plus_B_NTensor2D<T>* Masked_A_plus_B_NTensor2D<T>::clone() const {
    return new Masked_A_plus_B_NTensor2D<T>(*this);
}

template<typename T>
void Masked_A_plus_B_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
    modified[3] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_plus_B_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** A_minus_B_NTensor2D ************************************ */

template<typename T>
void A_minus_B_NTensor2D<T>::process (
        Box2D domain, std::vector<NTensorField2D<T>*> fields )
{
    PLB_PRECONDITION( fields.size()==3 );
    NTensorField2D<T>& A = *fields[0];
    NTensorField2D<T>& B = *fields[1];
    NTensorField2D<T>& result = *fields[2];
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot2D offsetB      = computeRelativeDisplacement(A,B);
    Dot2D offsetResult = computeRelativeDisplacement(A,result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                *result.get(iX+offsetResult.x,iY+offsetResult.y)
                    = *A.get(iX,iY) - *B.get(iX+offsetB.x,iY+offsetB.y);
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    result.get(iX+offsetResult.x,iY+offsetResult.y)[iDim]
                        = A.get(iX,iY)[iDim] - B.get(iX+offsetB.x,iY+offsetB.y)[iDim];
                }
            }
        }
    }
}

template<typename T>
A_minus_B_NTensor2D<T>* A_minus_B_NTensor2D<T>::clone() const {
    return new A_minus_B_NTensor2D<T>(*this);
}

template<typename T>
void A_minus_B_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_minus_B_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_minus_B_NTensor2D ************************************ */

template<typename T>
void Masked_A_minus_B_NTensor2D<T>::process (
        Box2D domain, std::vector<NTensorField2D<T>*> fields,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( fields.size()==3 );
    NTensorField2D<T>& A = *fields[0];
    NTensorField2D<T>& B = *fields[1];
    NTensorField2D<T>& result = *fields[2];
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot2D offsetB      = computeRelativeDisplacement(A,B);
    Dot2D offsetResult = computeRelativeDisplacement(A,result);
    Dot2D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    *result.get(iX+offsetResult.x,iY+offsetResult.y)
                        = *A.get(iX,iY) - *B.get(iX+offsetB.x,iY+offsetB.y);
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offsetResult.x,iY+offsetResult.y)[iDim]
                            = A.get(iX,iY)[iDim] - B.get(iX+offsetB.x,iY+offsetB.y)[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_minus_B_NTensor2D<T>* Masked_A_minus_B_NTensor2D<T>::clone() const {
    return new Masked_A_minus_B_NTensor2D<T>(*this);
}

template<typename T>
void Masked_A_minus_B_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
    modified[3] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_minus_B_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** A_times_B_NTensor2D ************************************ */

template<typename T>
void A_times_B_NTensor2D<T>::process (
        Box2D domain, std::vector<NTensorField2D<T>*> fields )
{
    PLB_PRECONDITION( fields.size()==3 );
    NTensorField2D<T>& A = *fields[0];
    NTensorField2D<T>& B = *fields[1];
    NTensorField2D<T>& result = *fields[2];
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot2D offsetB      = computeRelativeDisplacement(A,B);
    Dot2D offsetResult = computeRelativeDisplacement(A,result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                *result.get(iX+offsetResult.x,iY+offsetResult.y)
                    = *A.get(iX,iY) * *B.get(iX+offsetB.x,iY+offsetB.y);
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    result.get(iX+offsetResult.x,iY+offsetResult.y)[iDim]
                        = A.get(iX,iY)[iDim] * B.get(iX+offsetB.x,iY+offsetB.y)[iDim];
                }
            }
        }
    }
}

template<typename T>
A_times_B_NTensor2D<T>* A_times_B_NTensor2D<T>::clone() const {
    return new A_times_B_NTensor2D<T>(*this);
}

template<typename T>
void A_times_B_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_times_B_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_times_B_NTensor2D ************************************ */

template<typename T>
void Masked_A_times_B_NTensor2D<T>::process (
        Box2D domain, std::vector<NTensorField2D<T>*> fields,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( fields.size()==3 );
    NTensorField2D<T>& A = *fields[0];
    NTensorField2D<T>& B = *fields[1];
    NTensorField2D<T>& result = *fields[2];
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot2D offsetB      = computeRelativeDisplacement(A,B);
    Dot2D offsetResult = computeRelativeDisplacement(A,result);
    Dot2D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    *result.get(iX+offsetResult.x,iY+offsetResult.y)
                        = *A.get(iX,iY) * *B.get(iX+offsetB.x,iY+offsetB.y);
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offsetResult.x,iY+offsetResult.y)[iDim]
                            = A.get(iX,iY)[iDim] * B.get(iX+offsetB.x,iY+offsetB.y)[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_times_B_NTensor2D<T>* Masked_A_times_B_NTensor2D<T>::clone() const {
    return new Masked_A_times_B_NTensor2D<T>(*this);
}

template<typename T>
void Masked_A_times_B_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
    modified[3] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_times_B_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** A_dividedBy_B_NTensor2D ************************************ */

template<typename T>
void A_dividedBy_B_NTensor2D<T>::process (
        Box2D domain, std::vector<NTensorField2D<T>*> fields )
{
    PLB_PRECONDITION( fields.size()==3 );
    NTensorField2D<T>& A = *fields[0];
    NTensorField2D<T>& B = *fields[1];
    NTensorField2D<T>& result = *fields[2];
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot2D offsetB      = computeRelativeDisplacement(A,B);
    Dot2D offsetResult = computeRelativeDisplacement(A,result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                *result.get(iX+offsetResult.x,iY+offsetResult.y)
                    = *A.get(iX,iY) / *B.get(iX+offsetB.x,iY+offsetB.y);
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    result.get(iX+offsetResult.x,iY+offsetResult.y)[iDim]
                        = A.get(iX,iY)[iDim] / B.get(iX+offsetB.x,iY+offsetB.y)[iDim];
                }
            }
        }
    }
}

template<typename T>
A_dividedBy_B_NTensor2D<T>* A_dividedBy_B_NTensor2D<T>::clone() const {
    return new A_dividedBy_B_NTensor2D<T>(*this);
}

template<typename T>
void A_dividedBy_B_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_dividedBy_B_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_dividedBy_B_NTensor2D ************************************ */

template<typename T>
void Masked_A_dividedBy_B_NTensor2D<T>::process (
        Box2D domain, std::vector<NTensorField2D<T>*> fields,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( fields.size()==3 );
    NTensorField2D<T>& A = *fields[0];
    NTensorField2D<T>& B = *fields[1];
    NTensorField2D<T>& result = *fields[2];
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot2D offsetB      = computeRelativeDisplacement(A,B);
    Dot2D offsetResult = computeRelativeDisplacement(A,result);
    Dot2D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    *result.get(iX+offsetResult.x,iY+offsetResult.y)
                        = *A.get(iX,iY) / *B.get(iX+offsetB.x,iY+offsetB.y);
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offsetResult.x,iY+offsetResult.y)[iDim]
                            = A.get(iX,iY)[iDim] / B.get(iX+offsetB.x,iY+offsetB.y)[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_dividedBy_B_NTensor2D<T>* Masked_A_dividedBy_B_NTensor2D<T>::clone() const {
    return new Masked_A_dividedBy_B_NTensor2D<T>(*this);
}

template<typename T>
void Masked_A_dividedBy_B_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
    modified[3] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_dividedBy_B_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** A_toThePower_B_NTensor2D ************************************ */

template<typename T>
void A_toThePower_B_NTensor2D<T>::process (
        Box2D domain, std::vector<NTensorField2D<T>*> fields )
{
    PLB_PRECONDITION( fields.size()==3 );
    NTensorField2D<T>& A = *fields[0];
    NTensorField2D<T>& B = *fields[1];
    NTensorField2D<T>& result = *fields[2];
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot2D offsetB      = computeRelativeDisplacement(A,B);
    Dot2D offsetResult = computeRelativeDisplacement(A,result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                *result.get(iX+offsetResult.x,iY+offsetResult.y)
                    = customPower(*A.get(iX,iY),*B.get(iX+offsetB.x,iY+offsetB.y));
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    result.get(iX+offsetResult.x,iY+offsetResult.y)[iDim]
                        = customPower(A.get(iX,iY)[iDim],B.get(iX+offsetB.x,iY+offsetB.y)[iDim]);
                }
            }
        }
    }
}

template<typename T>
A_toThePower_B_NTensor2D<T>* A_toThePower_B_NTensor2D<T>::clone() const {
    return new A_toThePower_B_NTensor2D<T>(*this);
}

template<typename T>
void A_toThePower_B_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_toThePower_B_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_toThePower_B_NTensor2D ************************************ */

template<typename T>
void Masked_A_toThePower_B_NTensor2D<T>::process (
        Box2D domain, std::vector<NTensorField2D<T>*> fields,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( fields.size()==3 );
    NTensorField2D<T>& A = *fields[0];
    NTensorField2D<T>& B = *fields[1];
    NTensorField2D<T>& result = *fields[2];
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot2D offsetB      = computeRelativeDisplacement(A,B);
    Dot2D offsetResult = computeRelativeDisplacement(A,result);
    Dot2D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    *result.get(iX+offsetResult.x,iY+offsetResult.y)
                        = customPower(*A.get(iX,iY),*B.get(iX+offsetB.x,iY+offsetB.y));
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offsetResult.x,iY+offsetResult.y)[iDim]
                            = customPower(A.get(iX,iY)[iDim],B.get(iX+offsetB.x,iY+offsetB.y)[iDim]);
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_toThePower_B_NTensor2D<T>* Masked_A_toThePower_B_NTensor2D<T>::clone() const {
    return new Masked_A_toThePower_B_NTensor2D<T>(*this);
}

template<typename T>
void Masked_A_toThePower_B_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
    modified[3] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_toThePower_B_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** A_equals_B_NTensor2D ************************************ */

template<typename T>
void A_equals_B_NTensor2D<T>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==3 );
    NTensorField2D<T>& A = *dynamic_cast<NTensorField2D<T>*>(blocks[0]);
    NTensorField2D<T>& B = *dynamic_cast<NTensorField2D<T>*>(blocks[1]);
    NTensorField2D<int>& result = *dynamic_cast<NTensorField2D<int>*>(blocks[2]);
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot2D offsetB      = computeRelativeDisplacement(A,B);
    Dot2D offsetResult = computeRelativeDisplacement(A,result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                *result.get(iX+offsetResult.x,iY+offsetResult.y)
                    = *A.get(iX,iY) == *B.get(iX+offsetB.x,iY+offsetB.y) ? 1:0;
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    result.get(iX+offsetResult.x,iY+offsetResult.y)[iDim]
                        = A.get(iX,iY)[iDim] == B.get(iX+offsetB.x,iY+offsetB.y)[iDim] ? 1:0;
                }
            }
        }
    }
}

template<typename T>
A_equals_B_NTensor2D<T>* A_equals_B_NTensor2D<T>::clone() const {
    return new A_equals_B_NTensor2D<T>(*this);
}

template<typename T>
void A_equals_B_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_equals_B_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_equals_B_NTensor2D ************************************ */

template<typename T>
void Masked_A_equals_B_NTensor2D<T>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==3 );
    NTensorField2D<T>& A = *dynamic_cast<NTensorField2D<T>*>(blocks[0]);
    NTensorField2D<T>& B = *dynamic_cast<NTensorField2D<T>*>(blocks[1]);
    NTensorField2D<int>& result = *dynamic_cast<NTensorField2D<int>*>(blocks[2]);
    NTensorField2D<int>& mask = *dynamic_cast<NTensorField2D<int>*>(blocks[3]);
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot2D offsetB      = computeRelativeDisplacement(A,B);
    Dot2D offsetResult = computeRelativeDisplacement(A,result);
    Dot2D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    *result.get(iX+offsetResult.x,iY+offsetResult.y)
                        = *A.get(iX,iY) == *B.get(iX+offsetB.x,iY+offsetB.y) ? 1:0;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offsetResult.x,iY+offsetResult.y)[iDim]
                            = A.get(iX,iY)[iDim] == B.get(iX+offsetB.x,iY+offsetB.y)[iDim] ? 1:0;
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_equals_B_NTensor2D<T>* Masked_A_equals_B_NTensor2D<T>::clone() const {
    return new Masked_A_equals_B_NTensor2D<T>(*this);
}

template<typename T>
void Masked_A_equals_B_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
    modified[3] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_equals_B_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** A_equals_alpha_NTensor2D ************************************* */

template<typename T>
A_equals_alpha_NTensor2D<T>::A_equals_alpha_NTensor2D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_equals_alpha_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<int>& result )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot2D offset = computeRelativeDisplacement(A, result);

    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                *result.get(iX+offset.x,iY+offset.y)
                    = *A.get(iX,iY) == scalarAlpha ? 1 : 0;
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    result.get(iX+offset.x,iY+offset.y)[iDim]
                        = A.get(iX,iY)[iDim] == alpha[iDim] ? 1 : 0;
                }
            }
        }
    }
}

template<typename T>
A_equals_alpha_NTensor2D<T>* A_equals_alpha_NTensor2D<T>::clone() const {
    return new A_equals_alpha_NTensor2D<T>(*this);
}

template<typename T>
void A_equals_alpha_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_equals_alpha_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_equals_alpha_NTensor2D ************************************* */

template<typename T>
Masked_A_equals_alpha_NTensor2D<T>::Masked_A_equals_alpha_NTensor2D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void Masked_A_equals_alpha_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<int>& result,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot2D offset = computeRelativeDisplacement(A, result);
    Dot2D maskOfs = computeRelativeDisplacement(A, mask);

    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    *result.get(iX+offset.x,iY+offset.y)
                        = *A.get(iX,iY) == scalarAlpha ? 1 : 0;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offset.x,iY+offset.y)[iDim]
                            = A.get(iX,iY)[iDim] == alpha[iDim] ? 1 : 0;
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_equals_alpha_NTensor2D<T>* Masked_A_equals_alpha_NTensor2D<T>::clone() const {
    return new Masked_A_equals_alpha_NTensor2D<T>(*this);
}

template<typename T>
void Masked_A_equals_alpha_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_equals_alpha_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** A_lessThan_B_NTensor2D ************************************ */

template<typename T>
void A_lessThan_B_NTensor2D<T>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==3 );
    NTensorField2D<T>& A = *dynamic_cast<NTensorField2D<T>*>(blocks[0]);
    NTensorField2D<T>& B = *dynamic_cast<NTensorField2D<T>*>(blocks[1]);
    NTensorField2D<int>& result = *dynamic_cast<NTensorField2D<int>*>(blocks[2]);
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot2D offsetB      = computeRelativeDisplacement(A,B);
    Dot2D offsetResult = computeRelativeDisplacement(A,result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                *result.get(iX+offsetResult.x,iY+offsetResult.y)
                    = *A.get(iX,iY) < *B.get(iX+offsetB.x,iY+offsetB.y) ? 1:0;
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    result.get(iX+offsetResult.x,iY+offsetResult.y)[iDim]
                        = A.get(iX,iY)[iDim] < B.get(iX+offsetB.x,iY+offsetB.y)[iDim] ? 1:0;
                }
            }
        }
    }
}

template<typename T>
A_lessThan_B_NTensor2D<T>* A_lessThan_B_NTensor2D<T>::clone() const {
    return new A_lessThan_B_NTensor2D<T>(*this);
}

template<typename T>
void A_lessThan_B_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_lessThan_B_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}

/* ******** Masked_A_lessThan_B_NTensor2D ************************************ */

template<typename T>
void Masked_A_lessThan_B_NTensor2D<T>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==3 );
    NTensorField2D<T>& A = *dynamic_cast<NTensorField2D<T>*>(blocks[0]);
    NTensorField2D<T>& B = *dynamic_cast<NTensorField2D<T>*>(blocks[1]);
    NTensorField2D<int>& result = *dynamic_cast<NTensorField2D<int>*>(blocks[2]);
    NTensorField2D<int>& mask = *dynamic_cast<NTensorField2D<int>*>(blocks[3]);
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot2D offsetB      = computeRelativeDisplacement(A,B);
    Dot2D offsetResult = computeRelativeDisplacement(A,result);
    Dot2D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    *result.get(iX+offsetResult.x,iY+offsetResult.y)
                        = *A.get(iX,iY) < *B.get(iX+offsetB.x,iY+offsetB.y) ? 1:0;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offsetResult.x,iY+offsetResult.y)[iDim]
                            = A.get(iX,iY)[iDim] < B.get(iX+offsetB.x,iY+offsetB.y)[iDim] ? 1:0;
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_lessThan_B_NTensor2D<T>* Masked_A_lessThan_B_NTensor2D<T>::clone() const {
    return new Masked_A_lessThan_B_NTensor2D<T>(*this);
}

template<typename T>
void Masked_A_lessThan_B_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
    modified[3] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_lessThan_B_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** A_lessThan_alpha_NTensor2D ************************************* */

template<typename T>
A_lessThan_alpha_NTensor2D<T>::A_lessThan_alpha_NTensor2D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_lessThan_alpha_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<int>& result )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot2D offset = computeRelativeDisplacement(A, result);

    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                *result.get(iX+offset.x,iY+offset.y)
                    = *A.get(iX,iY) < scalarAlpha ? 1 : 0;
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    result.get(iX+offset.x,iY+offset.y)[iDim]
                        = A.get(iX,iY)[iDim] < alpha[iDim] ? 1 : 0;
                }
            }
        }
    }
}

template<typename T>
A_lessThan_alpha_NTensor2D<T>* A_lessThan_alpha_NTensor2D<T>::clone() const {
    return new A_lessThan_alpha_NTensor2D<T>(*this);
}

template<typename T>
void A_lessThan_alpha_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_lessThan_alpha_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_lessThan_alpha_NTensor2D ************************************* */

template<typename T>
Masked_A_lessThan_alpha_NTensor2D<T>::Masked_A_lessThan_alpha_NTensor2D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void Masked_A_lessThan_alpha_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<int>& result,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot2D offset = computeRelativeDisplacement(A, result);
    Dot2D maskOfs = computeRelativeDisplacement(A, mask);

    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    *result.get(iX+offset.x,iY+offset.y)
                        = *A.get(iX,iY) < scalarAlpha ? 1 : 0;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offset.x,iY+offset.y)[iDim]
                            = A.get(iX,iY)[iDim] < alpha[iDim] ? 1 : 0;
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_lessThan_alpha_NTensor2D<T>* Masked_A_lessThan_alpha_NTensor2D<T>::clone() const {
    return new Masked_A_lessThan_alpha_NTensor2D<T>(*this);
}

template<typename T>
void Masked_A_lessThan_alpha_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_lessThan_alpha_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** A_lessEqual_B_NTensor2D ************************************ */

template<typename T>
void A_lessEqual_B_NTensor2D<T>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==3 );
    NTensorField2D<T>& A = *dynamic_cast<NTensorField2D<T>*>(blocks[0]);
    NTensorField2D<T>& B = *dynamic_cast<NTensorField2D<T>*>(blocks[1]);
    NTensorField2D<int>& result = *dynamic_cast<NTensorField2D<int>*>(blocks[2]);
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot2D offsetB      = computeRelativeDisplacement(A,B);
    Dot2D offsetResult = computeRelativeDisplacement(A,result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                *result.get(iX+offsetResult.x,iY+offsetResult.y)
                    = *A.get(iX,iY) <= *B.get(iX+offsetB.x,iY+offsetB.y) ? 1:0;
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    result.get(iX+offsetResult.x,iY+offsetResult.y)[iDim]
                        = A.get(iX,iY)[iDim] <= B.get(iX+offsetB.x,iY+offsetB.y)[iDim] ? 1:0;
                }
            }
        }
    }
}

template<typename T>
A_lessEqual_B_NTensor2D<T>* A_lessEqual_B_NTensor2D<T>::clone() const {
    return new A_lessEqual_B_NTensor2D<T>(*this);
}

template<typename T>
void A_lessEqual_B_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_lessEqual_B_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_lessEqual_B_NTensor2D ************************************ */

template<typename T>
void Masked_A_lessEqual_B_NTensor2D<T>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==3 );
    NTensorField2D<T>& A = *dynamic_cast<NTensorField2D<T>*>(blocks[0]);
    NTensorField2D<T>& B = *dynamic_cast<NTensorField2D<T>*>(blocks[1]);
    NTensorField2D<int>& result = *dynamic_cast<NTensorField2D<int>*>(blocks[2]);
    NTensorField2D<int>& mask = *dynamic_cast<NTensorField2D<int>*>(blocks[3]);
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot2D offsetB      = computeRelativeDisplacement(A,B);
    Dot2D offsetResult = computeRelativeDisplacement(A,result);
    Dot2D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    *result.get(iX+offsetResult.x,iY+offsetResult.y)
                        = *A.get(iX,iY) <= *B.get(iX+offsetB.x,iY+offsetB.y) ? 1:0;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offsetResult.x,iY+offsetResult.y)[iDim]
                            = A.get(iX,iY)[iDim] <= B.get(iX+offsetB.x,iY+offsetB.y)[iDim] ? 1:0;
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_lessEqual_B_NTensor2D<T>* Masked_A_lessEqual_B_NTensor2D<T>::clone() const {
    return new Masked_A_lessEqual_B_NTensor2D<T>(*this);
}

template<typename T>
void Masked_A_lessEqual_B_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
    modified[3] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_lessEqual_B_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** A_lessEqual_alpha_NTensor2D ************************************* */

template<typename T>
A_lessEqual_alpha_NTensor2D<T>::A_lessEqual_alpha_NTensor2D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_lessEqual_alpha_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<int>& result )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot2D offset = computeRelativeDisplacement(A, result);

    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                *result.get(iX+offset.x,iY+offset.y)
                    = *A.get(iX,iY) <= scalarAlpha ? 1 : 0;
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    result.get(iX+offset.x,iY+offset.y)[iDim]
                        = A.get(iX,iY)[iDim] <= alpha[iDim] ? 1 : 0;
                }
            }
        }
    }
}

template<typename T>
A_lessEqual_alpha_NTensor2D<T>* A_lessEqual_alpha_NTensor2D<T>::clone() const {
    return new A_lessEqual_alpha_NTensor2D<T>(*this);
}

template<typename T>
void A_lessEqual_alpha_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_lessEqual_alpha_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_lessEqual_alpha_NTensor2D ************************************* */

template<typename T>
Masked_A_lessEqual_alpha_NTensor2D<T>::Masked_A_lessEqual_alpha_NTensor2D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void Masked_A_lessEqual_alpha_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<int>& result,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot2D offset = computeRelativeDisplacement(A, result);
    Dot2D maskOfs = computeRelativeDisplacement(A, mask);

    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    *result.get(iX+offset.x,iY+offset.y)
                        = *A.get(iX,iY) <= scalarAlpha ? 1 : 0;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offset.x,iY+offset.y)[iDim]
                            = A.get(iX,iY)[iDim] <= alpha[iDim] ? 1 : 0;
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_lessEqual_alpha_NTensor2D<T>* Masked_A_lessEqual_alpha_NTensor2D<T>::clone() const {
    return new Masked_A_lessEqual_alpha_NTensor2D<T>(*this);
}

template<typename T>
void Masked_A_lessEqual_alpha_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_lessEqual_alpha_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** A_greaterThan_B_NTensor2D ************************************ */

template<typename T>
void A_greaterThan_B_NTensor2D<T>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==3 );
    NTensorField2D<T>& A = *dynamic_cast<NTensorField2D<T>*>(blocks[0]);
    NTensorField2D<T>& B = *dynamic_cast<NTensorField2D<T>*>(blocks[1]);
    NTensorField2D<int>& result = *dynamic_cast<NTensorField2D<int>*>(blocks[2]);
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot2D offsetB      = computeRelativeDisplacement(A,B);
    Dot2D offsetResult = computeRelativeDisplacement(A,result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                *result.get(iX+offsetResult.x,iY+offsetResult.y)
                    = *A.get(iX,iY) > *B.get(iX+offsetB.x,iY+offsetB.y) ? 1:0;
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    result.get(iX+offsetResult.x,iY+offsetResult.y)[iDim]
                        = A.get(iX,iY)[iDim] > B.get(iX+offsetB.x,iY+offsetB.y)[iDim] ? 1:0;
                }
            }
        }
    }
}

template<typename T>
A_greaterThan_B_NTensor2D<T>* A_greaterThan_B_NTensor2D<T>::clone() const {
    return new A_greaterThan_B_NTensor2D<T>(*this);
}

template<typename T>
void A_greaterThan_B_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_greaterThan_B_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_greaterThan_B_NTensor2D ************************************ */

template<typename T>
void Masked_A_greaterThan_B_NTensor2D<T>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==3 );
    NTensorField2D<T>& A = *dynamic_cast<NTensorField2D<T>*>(blocks[0]);
    NTensorField2D<T>& B = *dynamic_cast<NTensorField2D<T>*>(blocks[1]);
    NTensorField2D<int>& result = *dynamic_cast<NTensorField2D<int>*>(blocks[2]);
    NTensorField2D<int>& mask = *dynamic_cast<NTensorField2D<int>*>(blocks[3]);
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot2D offsetB      = computeRelativeDisplacement(A,B);
    Dot2D offsetResult = computeRelativeDisplacement(A,result);
    Dot2D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    *result.get(iX+offsetResult.x,iY+offsetResult.y)
                        = *A.get(iX,iY) > *B.get(iX+offsetB.x,iY+offsetB.y) ? 1:0;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offsetResult.x,iY+offsetResult.y)[iDim]
                            = A.get(iX,iY)[iDim] > B.get(iX+offsetB.x,iY+offsetB.y)[iDim] ? 1:0;
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_greaterThan_B_NTensor2D<T>* Masked_A_greaterThan_B_NTensor2D<T>::clone() const {
    return new Masked_A_greaterThan_B_NTensor2D<T>(*this);
}

template<typename T>
void Masked_A_greaterThan_B_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
    modified[3] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_greaterThan_B_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** A_greaterThan_alpha_NTensor2D ************************************* */

template<typename T>
A_greaterThan_alpha_NTensor2D<T>::A_greaterThan_alpha_NTensor2D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_greaterThan_alpha_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<int>& result )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot2D offset = computeRelativeDisplacement(A, result);

    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                *result.get(iX+offset.x,iY+offset.y)
                    = *A.get(iX,iY) > scalarAlpha ? 1 : 0;
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    result.get(iX+offset.x,iY+offset.y)[iDim]
                        = A.get(iX,iY)[iDim] > alpha[iDim] ? 1 : 0;
                }
            }
        }
    }
}

template<typename T>
A_greaterThan_alpha_NTensor2D<T>* A_greaterThan_alpha_NTensor2D<T>::clone() const {
    return new A_greaterThan_alpha_NTensor2D<T>(*this);
}

template<typename T>
void A_greaterThan_alpha_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_greaterThan_alpha_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_greaterThan_alpha_NTensor2D ************************************* */

template<typename T>
Masked_A_greaterThan_alpha_NTensor2D<T>::Masked_A_greaterThan_alpha_NTensor2D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void Masked_A_greaterThan_alpha_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<int>& result,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot2D offset = computeRelativeDisplacement(A, result);
    Dot2D maskOfs = computeRelativeDisplacement(A, mask);

    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    *result.get(iX+offset.x,iY+offset.y)
                        = *A.get(iX,iY) > scalarAlpha ? 1 : 0;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offset.x,iY+offset.y)[iDim]
                            = A.get(iX,iY)[iDim] > alpha[iDim] ? 1 : 0;
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_greaterThan_alpha_NTensor2D<T>* Masked_A_greaterThan_alpha_NTensor2D<T>::clone() const {
    return new Masked_A_greaterThan_alpha_NTensor2D<T>(*this);
}

template<typename T>
void Masked_A_greaterThan_alpha_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_greaterThan_alpha_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** A_greaterEqual_B_NTensor2D ************************************ */

template<typename T>
void A_greaterEqual_B_NTensor2D<T>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==3 );
    NTensorField2D<T>& A = *dynamic_cast<NTensorField2D<T>*>(blocks[0]);
    NTensorField2D<T>& B = *dynamic_cast<NTensorField2D<T>*>(blocks[1]);
    NTensorField2D<int>& result = *dynamic_cast<NTensorField2D<int>*>(blocks[2]);
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot2D offsetB      = computeRelativeDisplacement(A,B);
    Dot2D offsetResult = computeRelativeDisplacement(A,result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                *result.get(iX+offsetResult.x,iY+offsetResult.y)
                    = *A.get(iX,iY) >= *B.get(iX+offsetB.x,iY+offsetB.y) ? 1:0;
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    result.get(iX+offsetResult.x,iY+offsetResult.y)[iDim]
                        = A.get(iX,iY)[iDim] >= B.get(iX+offsetB.x,iY+offsetB.y)[iDim] ? 1:0;
                }
            }
        }
    }
}

template<typename T>
A_greaterEqual_B_NTensor2D<T>* A_greaterEqual_B_NTensor2D<T>::clone() const {
    return new A_greaterEqual_B_NTensor2D<T>(*this);
}

template<typename T>
void A_greaterEqual_B_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_greaterEqual_B_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_greaterEqual_B_NTensor2D ************************************ */

template<typename T>
void Masked_A_greaterEqual_B_NTensor2D<T>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==3 );
    NTensorField2D<T>& A = *dynamic_cast<NTensorField2D<T>*>(blocks[0]);
    NTensorField2D<T>& B = *dynamic_cast<NTensorField2D<T>*>(blocks[1]);
    NTensorField2D<int>& result = *dynamic_cast<NTensorField2D<int>*>(blocks[2]);
    NTensorField2D<int>& mask = *dynamic_cast<NTensorField2D<int>*>(blocks[3]);
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot2D offsetB      = computeRelativeDisplacement(A,B);
    Dot2D offsetResult = computeRelativeDisplacement(A,result);
    Dot2D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    *result.get(iX+offsetResult.x,iY+offsetResult.y)
                        = *A.get(iX,iY) >= *B.get(iX+offsetB.x,iY+offsetB.y) ? 1:0;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offsetResult.x,iY+offsetResult.y)[iDim]
                            = A.get(iX,iY)[iDim] >= B.get(iX+offsetB.x,iY+offsetB.y)[iDim] ? 1:0;
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_greaterEqual_B_NTensor2D<T>* Masked_A_greaterEqual_B_NTensor2D<T>::clone() const {
    return new Masked_A_greaterEqual_B_NTensor2D<T>(*this);
}

template<typename T>
void Masked_A_greaterEqual_B_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
    modified[3] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_greaterEqual_B_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** A_greaterEqual_alpha_NTensor2D ************************************* */

template<typename T>
A_greaterEqual_alpha_NTensor2D<T>::A_greaterEqual_alpha_NTensor2D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_greaterEqual_alpha_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<int>& result )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot2D offset = computeRelativeDisplacement(A, result);

    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                *result.get(iX+offset.x,iY+offset.y)
                    = *A.get(iX,iY) >= scalarAlpha ? 1 : 0;
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    result.get(iX+offset.x,iY+offset.y)[iDim]
                        = A.get(iX,iY)[iDim] >= alpha[iDim] ? 1 : 0;
                }
            }
        }
    }
}

template<typename T>
A_greaterEqual_alpha_NTensor2D<T>* A_greaterEqual_alpha_NTensor2D<T>::clone() const {
    return new A_greaterEqual_alpha_NTensor2D<T>(*this);
}

template<typename T>
void A_greaterEqual_alpha_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_greaterEqual_alpha_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_greaterEqual_alpha_NTensor2D ************************************* */

template<typename T>
Masked_A_greaterEqual_alpha_NTensor2D<T>::Masked_A_greaterEqual_alpha_NTensor2D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void Masked_A_greaterEqual_alpha_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<int>& result,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot2D offset = computeRelativeDisplacement(A, result);
    Dot2D maskOfs = computeRelativeDisplacement(A, mask);

    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    *result.get(iX+offset.x,iY+offset.y)
                        = *A.get(iX,iY) >= scalarAlpha ? 1 : 0;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offset.x,iY+offset.y)[iDim]
                            = A.get(iX,iY)[iDim] >= alpha[iDim] ? 1 : 0;
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_greaterEqual_alpha_NTensor2D<T>* Masked_A_greaterEqual_alpha_NTensor2D<T>::clone() const {
    return new Masked_A_greaterEqual_alpha_NTensor2D<T>(*this);
}

template<typename T>
void Masked_A_greaterEqual_alpha_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_greaterEqual_alpha_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** A_and_B_NTensor2D ************************************ */

template<typename T>
void A_and_B_NTensor2D<T>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==3 );
    NTensorField2D<T>& A = *dynamic_cast<NTensorField2D<T>*>(blocks[0]);
    NTensorField2D<T>& B = *dynamic_cast<NTensorField2D<T>*>(blocks[1]);
    NTensorField2D<int>& result = *dynamic_cast<NTensorField2D<int>*>(blocks[2]);
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot2D offsetB      = computeRelativeDisplacement(A,B);
    Dot2D offsetResult = computeRelativeDisplacement(A,result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                *result.get(iX+offsetResult.x,iY+offsetResult.y)
                    = (bool)(*A.get(iX,iY)) && (bool)(*B.get(iX+offsetB.x,iY+offsetB.y)) ? 1:0;
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    result.get(iX+offsetResult.x,iY+offsetResult.y)[iDim]
                        = (bool)A.get(iX,iY)[iDim] && (bool)B.get(iX+offsetB.x,iY+offsetB.y)[iDim] ? 1:0;
                }
            }
        }
    }
}

template<typename T>
A_and_B_NTensor2D<T>* A_and_B_NTensor2D<T>::clone() const {
    return new A_and_B_NTensor2D<T>(*this);
}

template<typename T>
void A_and_B_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_and_B_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_and_B_NTensor2D ************************************ */

template<typename T>
void Masked_A_and_B_NTensor2D<T>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==3 );
    NTensorField2D<T>& A = *dynamic_cast<NTensorField2D<T>*>(blocks[0]);
    NTensorField2D<T>& B = *dynamic_cast<NTensorField2D<T>*>(blocks[1]);
    NTensorField2D<int>& result = *dynamic_cast<NTensorField2D<int>*>(blocks[2]);
    NTensorField2D<int>& mask = *dynamic_cast<NTensorField2D<int>*>(blocks[3]);
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot2D offsetB      = computeRelativeDisplacement(A,B);
    Dot2D offsetResult = computeRelativeDisplacement(A,result);
    Dot2D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    *result.get(iX+offsetResult.x,iY+offsetResult.y)
                        = (bool)(*A.get(iX,iY)) && (bool)(*B.get(iX+offsetB.x,iY+offsetB.y)) ? 1:0;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offsetResult.x,iY+offsetResult.y)[iDim]
                            = (bool)A.get(iX,iY)[iDim] && (bool)B.get(iX+offsetB.x,iY+offsetB.y)[iDim] ? 1:0;
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_and_B_NTensor2D<T>* Masked_A_and_B_NTensor2D<T>::clone() const {
    return new Masked_A_and_B_NTensor2D<T>(*this);
}

template<typename T>
void Masked_A_and_B_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
    modified[3] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_and_B_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}



/* ******** A_or_B_NTensor2D ************************************ */

template<typename T>
void A_or_B_NTensor2D<T>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==3 );
    NTensorField2D<T>& A = *dynamic_cast<NTensorField2D<T>*>(blocks[0]);
    NTensorField2D<T>& B = *dynamic_cast<NTensorField2D<T>*>(blocks[1]);
    NTensorField2D<int>& result = *dynamic_cast<NTensorField2D<int>*>(blocks[2]);
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot2D offsetB      = computeRelativeDisplacement(A,B);
    Dot2D offsetResult = computeRelativeDisplacement(A,result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                *result.get(iX+offsetResult.x,iY+offsetResult.y)
                    = (bool)(*A.get(iX,iY)) || (bool)(*B.get(iX+offsetB.x,iY+offsetB.y)) ? 1:0;
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    result.get(iX+offsetResult.x,iY+offsetResult.y)[iDim]
                        = (bool)A.get(iX,iY)[iDim] || (bool)B.get(iX+offsetB.x,iY+offsetB.y)[iDim] ? 1:0;
                }
            }
        }
    }
}

template<typename T>
A_or_B_NTensor2D<T>* A_or_B_NTensor2D<T>::clone() const {
    return new A_or_B_NTensor2D<T>(*this);
}

template<typename T>
void A_or_B_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_or_B_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_or_B_NTensor2D ************************************ */

template<typename T>
void Masked_A_or_B_NTensor2D<T>::processGenericBlocks (
        Box2D domain, std::vector<AtomicBlock2D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==3 );
    NTensorField2D<T>& A = *dynamic_cast<NTensorField2D<T>*>(blocks[0]);
    NTensorField2D<T>& B = *dynamic_cast<NTensorField2D<T>*>(blocks[1]);
    NTensorField2D<int>& result = *dynamic_cast<NTensorField2D<int>*>(blocks[2]);
    NTensorField2D<int>& mask = *dynamic_cast<NTensorField2D<int>*>(blocks[3]);
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot2D offsetB      = computeRelativeDisplacement(A,B);
    Dot2D offsetResult = computeRelativeDisplacement(A,result);
    Dot2D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    *result.get(iX+offsetResult.x,iY+offsetResult.y)
                        = (bool)(*A.get(iX,iY)) || (bool)(*B.get(iX+offsetB.x,iY+offsetB.y)) ? 1:0;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offsetResult.x,iY+offsetResult.y)[iDim]
                            = (bool)A.get(iX,iY)[iDim] || (bool)B.get(iX+offsetB.x,iY+offsetB.y)[iDim] ? 1:0;
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_or_B_NTensor2D<T>* Masked_A_or_B_NTensor2D<T>::clone() const {
    return new Masked_A_or_B_NTensor2D<T>(*this);
}

template<typename T>
void Masked_A_or_B_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
    modified[3] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_or_B_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Not_A_NTensor2D ************************************* */

template<typename T>
void Not_A_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<int>& result )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    Dot2D offset = computeRelativeDisplacement(A, result);

    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                *result.get(iX+offset.x,iY+offset.y)
                    = *A.get(iX,iY) == T() ? 1 : 0;
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    result.get(iX+offset.x,iY+offset.y)[iDim]
                        = A.get(iX,iY)[iDim] == T() ? 1 : 0;
                }
            }
        }
    }
}

template<typename T>
Not_A_NTensor2D<T>* Not_A_NTensor2D<T>::clone() const {
    return new Not_A_NTensor2D<T>(*this);
}

template<typename T>
void Not_A_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT Not_A_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_Not_A_NTensor2D ************************************* */

template<typename T>
void Masked_Not_A_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<int>& result,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    Dot2D offset = computeRelativeDisplacement(A, result);
    Dot2D maskOfs = computeRelativeDisplacement(A, mask);

    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    *result.get(iX+offset.x,iY+offset.y)
                        = *A.get(iX,iY) == T() ? 1 : 0;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offset.x,iY+offset.y)[iDim]
                            = A.get(iX,iY)[iDim] == T() ? 1 : 0;
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_Not_A_NTensor2D<T>* Masked_Not_A_NTensor2D<T>::clone() const {
    return new Masked_Not_A_NTensor2D<T>(*this);
}

template<typename T>
void Masked_Not_A_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_Not_A_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulk;
}



/* ******** A_plus_B_inplace_NTensor2D ************************************ */

template<typename T>
void A_plus_B_inplace_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<T>& B)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    plint ndim = A.getNdim();

    Dot2D offset = computeRelativeDisplacement(A,B);
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                *A.get(iX,iY) += *B.get(iX+offset.x,iY+offset.y);
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    A.get(iX,iY)[iDim] += B.get(iX+offset.x,iY+offset.y)[iDim];
                }
            }
        }
    }
}

template<typename T>
A_plus_B_inplace_NTensor2D<T>* A_plus_B_inplace_NTensor2D<T>::clone() const {
    return new A_plus_B_inplace_NTensor2D<T>(*this);
}

template<typename T>
void A_plus_B_inplace_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT A_plus_B_inplace_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Masked_A_plus_B_inplace_NTensor2D ************************************ */

template<typename T>
void Masked_A_plus_B_inplace_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<T>& B,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    plint ndim = A.getNdim();

    Dot2D offset = computeRelativeDisplacement(A,B);
    Dot2D maskOfs = computeRelativeDisplacement(A, mask);
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    *A.get(iX,iY) += *B.get(iX+offset.x,iY+offset.y);
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        A.get(iX,iY)[iDim] += B.get(iX+offset.x,iY+offset.y)[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_plus_B_inplace_NTensor2D<T>* Masked_A_plus_B_inplace_NTensor2D<T>::clone() const {
    return new Masked_A_plus_B_inplace_NTensor2D<T>(*this);
}

template<typename T>
void Masked_A_plus_B_inplace_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_plus_B_inplace_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_plus_alpha_inplace_NTensor2D ************************************ */

template<typename T>
A_plus_alpha_inplace_NTensor2D<T>::A_plus_alpha_inplace_NTensor2D (
        std::vector<T> const& alpha_ )
    : alpha(alpha_)
{ }

template<typename T>
void A_plus_alpha_inplace_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A)
{
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size());
    plint ndim = A.getNdim();

    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalar = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                *A.get(iX,iY) += scalar;
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    A.get(iX,iY)[iDim] += alpha[iDim];
                }
            }
        }
    }
}

template<typename T>
A_plus_alpha_inplace_NTensor2D<T>* A_plus_alpha_inplace_NTensor2D<T>::clone() const {
    return new A_plus_alpha_inplace_NTensor2D<T>(*this);
}

template<typename T>
void A_plus_alpha_inplace_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_plus_alpha_inplace_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Masked_A_plus_alpha_inplace_NTensor2D ************************************ */

template<typename T>
Masked_A_plus_alpha_inplace_NTensor2D<T>::Masked_A_plus_alpha_inplace_NTensor2D (
        std::vector<T> const& alpha_ )
    : alpha(alpha_)
{ }

template<typename T>
void Masked_A_plus_alpha_inplace_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<int>& mask)
{
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size());
    Dot2D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();

    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalar = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    *A.get(iX,iY) += scalar;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        A.get(iX,iY)[iDim] += alpha[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_plus_alpha_inplace_NTensor2D<T>* Masked_A_plus_alpha_inplace_NTensor2D<T>::clone() const {
    return new Masked_A_plus_alpha_inplace_NTensor2D<T>(*this);
}

template<typename T>
void Masked_A_plus_alpha_inplace_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_plus_alpha_inplace_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_minus_B_inplace_NTensor2D ************************************ */

template<typename T>
void A_minus_B_inplace_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<T>& B)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    plint ndim = A.getNdim();

    Dot2D offset = computeRelativeDisplacement(A,B);
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                *A.get(iX,iY) -= *B.get(iX+offset.x,iY+offset.y);
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    A.get(iX,iY)[iDim] -= B.get(iX+offset.x,iY+offset.y)[iDim];
                }
            }
        }
    }
}

template<typename T>
A_minus_B_inplace_NTensor2D<T>* A_minus_B_inplace_NTensor2D<T>::clone() const {
    return new A_minus_B_inplace_NTensor2D<T>(*this);
}

template<typename T>
void A_minus_B_inplace_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT A_minus_B_inplace_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

/* ******** Masked_A_minus_B_inplace_NTensor2D ************************************ */

template<typename T>
void Masked_A_minus_B_inplace_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<T>& B,
        NTensorField2D<int>& mask)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    Dot2D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();

    Dot2D offset = computeRelativeDisplacement(A,B);
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    *A.get(iX,iY) -= *B.get(iX+offset.x,iY+offset.y);
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        A.get(iX,iY)[iDim] -= B.get(iX+offset.x,iY+offset.y)[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_minus_B_inplace_NTensor2D<T>* Masked_A_minus_B_inplace_NTensor2D<T>::clone() const {
    return new Masked_A_minus_B_inplace_NTensor2D<T>(*this);
}

template<typename T>
void Masked_A_minus_B_inplace_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_minus_B_inplace_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

/* ******** A_minus_alpha_inplace_NTensor2D ************************************ */

template<typename T>
A_minus_alpha_inplace_NTensor2D<T>::A_minus_alpha_inplace_NTensor2D (
        std::vector<T> const& alpha_ )
    : alpha(alpha_)
{ }

template<typename T>
void A_minus_alpha_inplace_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A)
{
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size());
    plint ndim = A.getNdim();

    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalar = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                *A.get(iX,iY) -= scalar;
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    A.get(iX,iY)[iDim] -= alpha[iDim];
                }
            }
        }
    }
}

template<typename T>
A_minus_alpha_inplace_NTensor2D<T>* A_minus_alpha_inplace_NTensor2D<T>::clone() const {
    return new A_minus_alpha_inplace_NTensor2D<T>(*this);
}

template<typename T>
void A_minus_alpha_inplace_NTensor2D<T>::
    getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_minus_alpha_inplace_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}



/* ******** Masked_A_minus_alpha_inplace_NTensor2D ************************************ */

template<typename T>
Masked_A_minus_alpha_inplace_NTensor2D<T>::Masked_A_minus_alpha_inplace_NTensor2D (
        std::vector<T> const& alpha_ )
    : alpha(alpha_)
{ }

template<typename T>
void Masked_A_minus_alpha_inplace_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size());
    Dot2D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();

    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalar = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    *A.get(iX,iY) -= scalar;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        A.get(iX,iY)[iDim] -= alpha[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_minus_alpha_inplace_NTensor2D<T>* Masked_A_minus_alpha_inplace_NTensor2D<T>::clone() const {
    return new Masked_A_minus_alpha_inplace_NTensor2D<T>(*this);
}

template<typename T>
void Masked_A_minus_alpha_inplace_NTensor2D<T>::
    getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_minus_alpha_inplace_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}



/* ******** A_times_B_inplace_NTensor2D ************************************ */

template<typename T>
void A_times_B_inplace_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<T>& B)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    plint ndim = A.getNdim();

    Dot2D offset = computeRelativeDisplacement(A,B);
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                *A.get(iX,iY) *= *B.get(iX+offset.x,iY+offset.y);
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    A.get(iX,iY)[iDim] *= B.get(iX+offset.x,iY+offset.y)[iDim];
                }
            }
        }
    }
}

template<typename T>
A_times_B_inplace_NTensor2D<T>* A_times_B_inplace_NTensor2D<T>::clone() const {
    return new A_times_B_inplace_NTensor2D<T>(*this);
}

template<typename T>
void A_times_B_inplace_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT A_times_B_inplace_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Masked_A_times_B_inplace_NTensor2D ************************************ */

template<typename T>
void Masked_A_times_B_inplace_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<T>& B,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    plint ndim = A.getNdim();

    Dot2D offset = computeRelativeDisplacement(A,B);
    Dot2D maskOfs = computeRelativeDisplacement(A, mask);
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    *A.get(iX,iY) *= *B.get(iX+offset.x,iY+offset.y);
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        A.get(iX,iY)[iDim] *= B.get(iX+offset.x,iY+offset.y)[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_times_B_inplace_NTensor2D<T>* Masked_A_times_B_inplace_NTensor2D<T>::clone() const {
    return new Masked_A_times_B_inplace_NTensor2D<T>(*this);
}

template<typename T>
void Masked_A_times_B_inplace_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_times_B_inplace_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_times_alpha_inplace_NTensor2D ************************************ */

template<typename T>
A_times_alpha_inplace_NTensor2D<T>::A_times_alpha_inplace_NTensor2D (
        std::vector<T> const& alpha_ )
    : alpha(alpha_)
{ }

template<typename T>
void A_times_alpha_inplace_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A)
{
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size());
    plint ndim = A.getNdim();

    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalar = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                *A.get(iX,iY) *= scalar;
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    A.get(iX,iY)[iDim] *= alpha[iDim];
                }
            }
        }
    }
}

template<typename T>
A_times_alpha_inplace_NTensor2D<T>* A_times_alpha_inplace_NTensor2D<T>::clone() const {
    return new A_times_alpha_inplace_NTensor2D<T>(*this);
}

template<typename T>
void A_times_alpha_inplace_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_times_alpha_inplace_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}



/* ******** Masked_A_times_alpha_inplace_NTensor2D ************************************ */

template<typename T>
Masked_A_times_alpha_inplace_NTensor2D<T>::Masked_A_times_alpha_inplace_NTensor2D (
        std::vector<T> const& alpha_ )
    : alpha(alpha_)
{ }

template<typename T>
void Masked_A_times_alpha_inplace_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<int>& mask)
{
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size());
    Dot2D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();

    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalar = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    *A.get(iX,iY) *= scalar;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        A.get(iX,iY)[iDim] *= alpha[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_times_alpha_inplace_NTensor2D<T>* Masked_A_times_alpha_inplace_NTensor2D<T>::clone() const {
    return new Masked_A_times_alpha_inplace_NTensor2D<T>(*this);
}

template<typename T>
void Masked_A_times_alpha_inplace_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_times_alpha_inplace_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}



/* ******** A_dividedBy_B_inplace_NTensor2D ************************************ */

template<typename T>
void A_dividedBy_B_inplace_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<T>& B)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    plint ndim = A.getNdim();

    Dot2D offset = computeRelativeDisplacement(A,B);
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                *A.get(iX,iY) /= *B.get(iX+offset.x,iY+offset.y);
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    A.get(iX,iY)[iDim] /= B.get(iX+offset.x,iY+offset.y)[iDim];
                }
            }
        }
    }
}

template<typename T>
A_dividedBy_B_inplace_NTensor2D<T>* A_dividedBy_B_inplace_NTensor2D<T>::clone() const {
    return new A_dividedBy_B_inplace_NTensor2D<T>(*this);
}

template<typename T>
void A_dividedBy_B_inplace_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT A_dividedBy_B_inplace_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Masked_A_dividedBy_B_inplace_NTensor2D ************************************ */

template<typename T>
void Masked_A_dividedBy_B_inplace_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<T>& B,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    Dot2D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();

    Dot2D offset = computeRelativeDisplacement(A,B);
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    *A.get(iX,iY) /= *B.get(iX+offset.x,iY+offset.y);
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        A.get(iX,iY)[iDim] /= B.get(iX+offset.x,iY+offset.y)[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_dividedBy_B_inplace_NTensor2D<T>* Masked_A_dividedBy_B_inplace_NTensor2D<T>::clone() const {
    return new Masked_A_dividedBy_B_inplace_NTensor2D<T>(*this);
}

template<typename T>
void Masked_A_dividedBy_B_inplace_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_dividedBy_B_inplace_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_dividedBy_alpha_inplace_NTensor2D ************************************ */

template<typename T>
A_dividedBy_alpha_inplace_NTensor2D<T>::A_dividedBy_alpha_inplace_NTensor2D (
        std::vector<T> const& alpha_ )
    : alpha(alpha_)
{ }

template<typename T>
void A_dividedBy_alpha_inplace_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A)
{
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size());
    plint ndim = A.getNdim();

    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalar = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                *A.get(iX,iY) /= scalar;
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    A.get(iX,iY)[iDim] /= alpha[iDim];
                }
            }
        }
    }
}

template<typename T>
A_dividedBy_alpha_inplace_NTensor2D<T>* A_dividedBy_alpha_inplace_NTensor2D<T>::clone() const {
    return new A_dividedBy_alpha_inplace_NTensor2D<T>(*this);
}

template<typename T>
void A_dividedBy_alpha_inplace_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_dividedBy_alpha_inplace_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Masked_A_dividedBy_alpha_inplace_NTensor2D ************************************ */

template<typename T>
Masked_A_dividedBy_alpha_inplace_NTensor2D<T>::Masked_A_dividedBy_alpha_inplace_NTensor2D (
        std::vector<T> const& alpha_ )
    : alpha(alpha_)
{ }

template<typename T>
void Masked_A_dividedBy_alpha_inplace_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size());
    Dot2D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();

    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalar = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    *A.get(iX,iY) /= scalar;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        A.get(iX,iY)[iDim] /= alpha[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_dividedBy_alpha_inplace_NTensor2D<T>* Masked_A_dividedBy_alpha_inplace_NTensor2D<T>::clone() const {
    return new Masked_A_dividedBy_alpha_inplace_NTensor2D<T>(*this);
}

template<typename T>
void Masked_A_dividedBy_alpha_inplace_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_dividedBy_alpha_inplace_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_toThePower_B_inplace_NTensor2D ************************************ */

template<typename T>
void A_toThePower_B_inplace_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<T>& B)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    plint ndim = A.getNdim();

    Dot2D offset = computeRelativeDisplacement(A,B);
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                customInPlacePower(*A.get(iX,iY),*B.get(iX+offset.x,iY+offset.y));
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    customInPlacePower(A.get(iX,iY)[iDim],B.get(iX+offset.x,iY+offset.y)[iDim]);
                }
            }
        }
    }
}

template<typename T>
A_toThePower_B_inplace_NTensor2D<T>* A_toThePower_B_inplace_NTensor2D<T>::clone() const {
    return new A_toThePower_B_inplace_NTensor2D<T>(*this);
}

template<typename T>
void A_toThePower_B_inplace_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT A_toThePower_B_inplace_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Masked_A_toThePower_B_inplace_NTensor2D ************************************ */

template<typename T>
void Masked_A_toThePower_B_inplace_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<T>& B,
        NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    Dot2D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();

    Dot2D offset = computeRelativeDisplacement(A,B);
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    customInPlacePower(*A.get(iX,iY),*B.get(iX+offset.x,iY+offset.y));
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        customInPlacePower(A.get(iX,iY)[iDim],B.get(iX+offset.x,iY+offset.y)[iDim]);
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_toThePower_B_inplace_NTensor2D<T>* Masked_A_toThePower_B_inplace_NTensor2D<T>::clone() const {
    return new Masked_A_toThePower_B_inplace_NTensor2D<T>(*this);
}

template<typename T>
void Masked_A_toThePower_B_inplace_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_toThePower_B_inplace_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_toThePower_alpha_inplace_NTensor2D ************************************ */

template<typename T>
A_toThePower_alpha_inplace_NTensor2D<T>::A_toThePower_alpha_inplace_NTensor2D (
        std::vector<T> const& alpha_ )
    : alpha(alpha_)
{ }

template<typename T>
void A_toThePower_alpha_inplace_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A)
{
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size());
    plint ndim = A.getNdim();

    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalar = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                customInPlacePower(*A.get(iX,iY),scalar);
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    customInPlacePower(A.get(iX,iY)[iDim],alpha[iDim]);
                }
            }
        }
    }
}

template<typename T>
A_toThePower_alpha_inplace_NTensor2D<T>* A_toThePower_alpha_inplace_NTensor2D<T>::clone() const {
    return new A_toThePower_alpha_inplace_NTensor2D<T>(*this);
}

template<typename T>
void A_toThePower_alpha_inplace_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_toThePower_alpha_inplace_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

/* ******** Masked_A_toThePower_alpha_inplace_NTensor2D ************************************ */

template<typename T>
Masked_A_toThePower_alpha_inplace_NTensor2D<T>::Masked_A_toThePower_alpha_inplace_NTensor2D (
        std::vector<T> const& alpha_ )
    : alpha(alpha_)
{ }

template<typename T>
void Masked_A_toThePower_alpha_inplace_NTensor2D<T>::process (
        Box2D domain, NTensorField2D<T>& A, NTensorField2D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size());
    Dot2D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();

    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalar = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    customInPlacePower(*A.get(iX,iY),scalar);
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        customInPlacePower(A.get(iX,iY)[iDim],alpha[iDim]);
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_toThePower_alpha_inplace_NTensor2D<T>* Masked_A_toThePower_alpha_inplace_NTensor2D<T>::clone() const {
    return new Masked_A_toThePower_alpha_inplace_NTensor2D<T>(*this);
}

template<typename T>
void Masked_A_toThePower_alpha_inplace_NTensor2D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_toThePower_alpha_inplace_NTensor2D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

}  // namespace plb

#endif  // DATA_ANALYSIS_FUNCTIONAL_2D_HH
