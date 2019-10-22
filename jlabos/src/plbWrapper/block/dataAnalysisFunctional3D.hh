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
#ifndef DATA_ANALYSIS_FUNCTIONAL_3D_HH
#define DATA_ANALYSIS_FUNCTIONAL_3D_HH

#include "plbWrapper/block/dataAnalysisFunctional3D.h"
#include "plbWrapper/block/plbMath.h"
#include "core/blockStatistics.h"
#include "core/plbDebug.h"
#include "core/util.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include "finiteDifference/fdStencils1D.h"
#include "atomicBlock/atomicBlock3D.h"
#include "atomicBlock/blockLattice3D.h"
#include "atomicBlock/dataField3D.h"
#include <cmath>
#include <limits>

namespace plb {

namespace fdNTensorField {

template<typename T>
inline T bulkXderiv (
        NTensorField3D<T> const& field, plint iX, plint iY, plint iZ, int iD )
{
    T dxu = fd::ctl_diff( field.get(iX+1,iY,iZ)[iD],
                          field.get(iX-1,iY,iZ)[iD] );
    return dxu;
}

template<typename T>
inline T bulkYderiv (
        NTensorField3D<T> const& field, plint iX, plint iY, plint iZ, int iD )
{
    T dyu = fd::ctl_diff( field.get(iX,iY+1,iZ)[iD],
                          field.get(iX,iY-1,iZ)[iD] );
    return dyu;
}

template<typename T>
inline T bulkZderiv (
        NTensorField3D<T> const& field, plint iX, plint iY, plint iZ, int iD )
{
    T dzu = fd::ctl_diff( field.get(iX,iY,iZ+1)[iD],
                          field.get(iX,iY,iZ-1)[iD] );
    return dzu;
}

template<typename T>
inline T planeXderiv (
        NTensorField3D<T> const& field, int direction, int orientation,
        plint iX, plint iY, plint iZ, int iD )
{
    if (direction==0) {
        return -orientation *
            fd::o1_fwd_diff( field.get(iX              ,iY,iZ)[iD],
                             field.get(iX-1*orientation,iY,iZ)[iD] );
    }
    else {
        return bulkXderiv(field, iX,iY,iZ, iD);
    }
}

template<typename T>
inline T planeYderiv (
        NTensorField3D<T> const& field, int direction, int orientation,
        plint iX, plint iY, plint iZ, int iD )
{
    if (direction==1) {
        return -orientation *
            fd::o1_fwd_diff( field.get(iX,iY              ,iZ)[iD],
                             field.get(iX,iY-1*orientation,iZ)[iD] );
    }
    else {
        return bulkYderiv(field, iX,iY,iZ, iD);
    }
}

template<typename T>
inline T planeZderiv (
        NTensorField3D<T> const& field, int direction, int orientation,
        plint iX, plint iY, plint iZ, int iD )
{
    if (direction==2) {
        return -orientation *
            fd::o1_fwd_diff( field.get(iX,iY,iZ              )[iD],
                             field.get(iX,iY,iZ-1*orientation)[iD] );
    }
    else {
        return bulkZderiv(field, iX,iY,iZ, iD);
    }
}

template<typename T>
inline T edgeXderiv (
        NTensorField3D<T> const& field,
        int plane, int direction1, int direction2,
        plint iX, plint iY, plint iZ, int iD )
{
    if (plane==0) {
        return bulkXderiv(field, iX,iY,iZ, iD);
    }
    else {
        int orientation = plane==1 ? direction2 : direction1;
        return -orientation *
            fd::o1_fwd_diff( field.get(iX              ,iY,iZ)[iD],
                             field.get(iX-1*orientation,iY,iZ)[iD] );
    }
}

template<typename T>
inline T edgeYderiv (
        NTensorField3D<T> const& field,
        int plane, int direction1, int direction2,
        plint iX, plint iY, plint iZ, int iD )
{
    if (plane==1) {
        return bulkYderiv(field, iX,iY,iZ, iD);
    }
    else {
        int orientation = plane==0 ? direction1 : direction2;
        return -orientation *
            fd::o1_fwd_diff( field.get(iX,iY              ,iZ)[iD],
                             field.get(iX,iY-1*orientation,iZ)[iD] );
    }
}

template<typename T>
inline T edgeZderiv (
        NTensorField3D<T> const& field,
        int plane, int direction1, int direction2,
        plint iX, plint iY, plint iZ, int iD )
{
    if (plane==2) {
        return bulkZderiv(field, iX,iY,iZ, iD);
    }
    else {
        int orientation = plane==0 ? direction2 : direction1;
        return -orientation *
            fd::o1_fwd_diff( field.get(iX,iY,iZ              )[iD],
                             field.get(iX,iY,iZ-1*orientation)[iD] );
    }
}

template<typename T>
inline T cornerXderiv (
        NTensorField3D<T> const& field,
        int normalX, int normalY, int normalZ,
        plint iX, plint iY, plint iZ, int iD )
{
    int orientation = normalX;
    return -orientation *
        fd::o1_fwd_diff( field.get(iX              ,iY,iZ)[iD],
                         field.get(iX-1*orientation,iY,iZ)[iD] );
}

template<typename T>
inline T cornerYderiv (
        NTensorField3D<T> const& field,
        int normalX, int normalY, int normalZ,
        plint iX, plint iY, plint iZ, int iD )
{
    int orientation = normalY;
    return -orientation *
        fd::o1_fwd_diff( field.get(iX,iY              ,iZ)[iD],
                         field.get(iX,iY-1*orientation,iZ)[iD] );
}

template<typename T>
inline T cornerZderiv (
        NTensorField3D<T> const& field,
        int normalX, int normalY, int normalZ,
        plint iX, plint iY, plint iZ, int iD )
{
    int orientation = normalZ;
    return -orientation *
        fd::o1_fwd_diff( field.get(iX,iY,iZ              )[iD],
                         field.get(iX,iY,iZ-1*orientation)[iD] );
}

template<typename T>
inline T bulkVorticityX(NTensorField3D<T> const& velocity, plint iX, plint iY, plint iZ )
{
    T dyuz = bulkYderiv(velocity, iX,iY,iZ, 2);
    T dzuy = bulkZderiv(velocity, iX,iY,iZ, 1);

    return dyuz - dzuy;
}

template<typename T>
inline T bulkVorticityY(NTensorField3D<T> const& velocity, plint iX, plint iY, plint iZ )
{
    T dzux = bulkZderiv(velocity, iX,iY,iZ, 0);
    T dxuz = bulkXderiv(velocity, iX,iY,iZ, 2);

    return dzux - dxuz;
}

template<typename T>
inline T bulkVorticityZ(NTensorField3D<T> const& velocity, plint iX, plint iY, plint iZ )
{
    T dxuy = bulkXderiv(velocity, iX,iY,iZ, 1);
    T dyux = bulkYderiv(velocity, iX,iY,iZ, 0);

    return dxuy - dyux;
}

template<typename T>
inline T planeVorticityX( NTensorField3D<T> const& velocity, int direction, int orientation,
                          plint iX, plint iY, plint iZ )
{
    T dyuz = planeYderiv(velocity,direction,orientation, iX,iY,iZ, 2);
    T dzuy = planeZderiv(velocity,direction,orientation, iX,iY,iZ, 1);

    return dyuz - dzuy;
}

template<typename T>
inline T planeVorticityY( NTensorField3D<T> const& velocity, int direction, int orientation,
                          plint iX, plint iY, plint iZ )
{
    T dzux = planeZderiv(velocity,direction,orientation, iX,iY,iZ, 0);
    T dxuz = planeXderiv(velocity,direction,orientation, iX,iY,iZ, 2);

    return dzux - dxuz;
}

template<typename T>
inline T planeVorticityZ( NTensorField3D<T> const& velocity, int direction, int orientation,
                          plint iX, plint iY, plint iZ )
{
    T dxuy = planeXderiv(velocity,direction,orientation, iX,iY,iZ, 1);
    T dyux = planeYderiv(velocity,direction,orientation, iX,iY,iZ, 0);

    return dxuy - dyux;
}

template<typename T>
inline T edgeVorticityX( NTensorField3D<T> const& velocity, int plane, int normal1, int normal2,
                         plint iX, plint iY, plint iZ )
{
    T dyuz = edgeYderiv(velocity,plane,normal1,normal2, iX,iY,iZ, 2);
    T dzuy = edgeZderiv(velocity,plane,normal1,normal2, iX,iY,iZ, 1);

    return dyuz - dzuy;
}

template<typename T>
inline T edgeVorticityY( NTensorField3D<T> const& velocity, int plane, int normal1, int normal2,
                         plint iX, plint iY, plint iZ )
{
    T dzux = edgeZderiv(velocity,plane,normal1,normal2, iX,iY,iZ, 0);
    T dxuz = edgeXderiv(velocity,plane,normal1,normal2, iX,iY,iZ, 2);

    return dzux - dxuz;
}

template<typename T>
inline T edgeVorticityZ( NTensorField3D<T> const& velocity, int plane, int normal1, int normal2,
                         plint iX, plint iY, plint iZ )
{
    T dxuy = edgeXderiv(velocity,plane,normal1,normal2, iX,iY,iZ, 1);
    T dyux = edgeYderiv(velocity,plane,normal1,normal2, iX,iY,iZ, 0);

    return dxuy - dyux;
}

template<typename T>
inline T cornerVorticityX( NTensorField3D<T> const& velocity, int normalX, int normalY, int normalZ,
                           plint iX, plint iY, plint iZ )
{
    T dyuz = cornerYderiv(velocity,normalX,normalY,normalZ, iX,iY,iZ, 2);
    T dzuy = cornerZderiv(velocity,normalX,normalY,normalZ, iX,iY,iZ, 1);

    return dyuz - dzuy;
}

template<typename T>
inline T cornerVorticityY( NTensorField3D<T> const& velocity, int normalX, int normalY, int normalZ,
                           plint iX, plint iY, plint iZ )
{
    T dzux = cornerZderiv(velocity,normalX,normalY,normalZ, iX,iY,iZ, 0);
    T dxuz = cornerXderiv(velocity,normalX,normalY,normalZ, iX,iY,iZ, 2);

    return dzux - dxuz;
}

template<typename T>
inline T cornerVorticityZ( NTensorField3D<T> const& velocity, int normalX, int normalY, int normalZ,
                           plint iX, plint iY, plint iZ )
{
    T dxuy = cornerXderiv(velocity,normalX,normalY,normalZ, iX,iY,iZ, 1);
    T dyux = cornerYderiv(velocity,normalX,normalY,normalZ, iX,iY,iZ, 0);

    return dxuy - dyux;
}

}  // fdNTensorField

/* *************** Reductive Data Functionals for NTensorField ******** */

template<typename T>
BoxNTensorSumFunctional3D<T>::BoxNTensorSumFunctional3D(plint ndim)
    : sumVectorId(ndim)
{ 
    for (plint iDim=0; iDim<ndim; ++iDim) {
        sumVectorId[iDim] = this->getStatistics().subscribeSum();
    }
}

template<typename T>
void BoxNTensorSumFunctional3D<T>::process (
        Box3D domain, NTensorField3D<T>& vectorField )
{
    PLB_PRECONDITION( vectorField.getNdim() == (plint)sumVectorId.size() );
    plint ndim = vectorField.getNdim();
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    statistics.gatherSum(sumVectorId[iDim], (double)vectorField.get(iX,iY,iZ)[iDim]);
                }
            }
        }
    }
}

template<typename T>
BoxNTensorSumFunctional3D<T>* BoxNTensorSumFunctional3D<T>::clone() const
{
    return new BoxNTensorSumFunctional3D<T>(*this);
}

template<typename T>
void BoxNTensorSumFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;
}

template<typename T>
std::vector<T> BoxNTensorSumFunctional3D<T>::getSumVector() const {
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
MaskedBoxNTensorSumFunctional3D<T>::MaskedBoxNTensorSumFunctional3D(plint ndim)
    : sumVectorId(ndim)
{ 
    for (plint iDim=0; iDim<ndim; ++iDim) {
        sumVectorId[iDim] = this->getStatistics().subscribeSum();
    }
}

template<typename T>
void MaskedBoxNTensorSumFunctional3D<T>::process (
        Box3D domain, NTensorField3D<T>& vectorField,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( vectorField.getNdim() == (plint)sumVectorId.size() );
    Dot3D maskOfs = computeRelativeDisplacement(vectorField, mask);
    plint ndim = vectorField.getNdim();
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        statistics.gatherSum(sumVectorId[iDim], (double)vectorField.get(iX,iY,iZ)[iDim]);
                    }
                }
            }
        }
    }
}

template<typename T>
MaskedBoxNTensorSumFunctional3D<T>* MaskedBoxNTensorSumFunctional3D<T>::clone() const
{
    return new MaskedBoxNTensorSumFunctional3D<T>(*this);
}

template<typename T>
void MaskedBoxNTensorSumFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
}

template<typename T>
std::vector<T> MaskedBoxNTensorSumFunctional3D<T>::getSumVector() const {
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
BoxNTensorMinFunctional3D<T>::BoxNTensorMinFunctional3D(plint ndim)
    : maxVectorId(ndim)
{
    for (plint iDim=0; iDim<ndim; ++iDim) {
        maxVectorId[iDim] = this->getStatistics().subscribeMax();
    }
}

template<typename T>
void BoxNTensorMinFunctional3D<T>::process (
        Box3D domain, NTensorField3D<T>& vectorField )
{
    PLB_PRECONDITION( vectorField.getNdim() == (plint)maxVectorId.size() );
    plint ndim = vectorField.getNdim();
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                // BlockStatistics computes only maximum, no minimum. Therefore,
                //   the relation min(x) = -max(-x) is used.
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    statistics.gatherMax(maxVectorId[iDim], -(double)vectorField.get(iX,iY,iZ)[iDim]);
                }
            }
        }
    }
}

template<typename T>
BoxNTensorMinFunctional3D<T>* BoxNTensorMinFunctional3D<T>::clone() const
{
    return new BoxNTensorMinFunctional3D<T>(*this);
}

template<typename T>
void BoxNTensorMinFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;
}

template<typename T>
std::vector<T> BoxNTensorMinFunctional3D<T>::getMinVector() const {
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
MaskedBoxNTensorMinFunctional3D<T>::MaskedBoxNTensorMinFunctional3D(plint ndim)
    : maxVectorId(ndim)
{
    for (plint iDim=0; iDim<ndim; ++iDim) {
        maxVectorId[iDim] = this->getStatistics().subscribeMax();
    }
}

template<typename T>
void MaskedBoxNTensorMinFunctional3D<T>::process (
        Box3D domain, NTensorField3D<T>& vectorField,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( vectorField.getNdim() == (plint)maxVectorId.size() );
    Dot3D maskOfs = computeRelativeDisplacement(vectorField, mask);
    plint ndim = vectorField.getNdim();
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                    // BlockStatistics computes only maximum, no minimum. Therefore,
                    //   the relation min(x) = -max(-x) is used.
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        statistics.gatherMax(maxVectorId[iDim], -(double)vectorField.get(iX,iY,iZ)[iDim]);
                    }
                }
            }
        }
    }
}

template<typename T>
MaskedBoxNTensorMinFunctional3D<T>* MaskedBoxNTensorMinFunctional3D<T>::clone() const
{
    return new MaskedBoxNTensorMinFunctional3D<T>(*this);
}

template<typename T>
void MaskedBoxNTensorMinFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
}

template<typename T>
std::vector<T> MaskedBoxNTensorMinFunctional3D<T>::getMinVector() const {
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
BoxNTensorMaxFunctional3D<T>::BoxNTensorMaxFunctional3D(plint ndim)
    : maxVectorId(ndim)
{
    for (plint iDim=0; iDim<ndim; ++iDim) {
        maxVectorId[iDim] = this->getStatistics().subscribeMax();
    }
}

template<typename T>
void BoxNTensorMaxFunctional3D<T>::process (
        Box3D domain, NTensorField3D<T>& vectorField )
{
    PLB_PRECONDITION( vectorField.getNdim() == (plint)maxVectorId.size() );
    plint ndim = vectorField.getNdim();
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    statistics.gatherMax(maxVectorId[iDim], (double)vectorField.get(iX,iY,iZ)[iDim]);
                }
            }
        }
    }
}

template<typename T>
BoxNTensorMaxFunctional3D<T>* BoxNTensorMaxFunctional3D<T>::clone() const
{
    return new BoxNTensorMaxFunctional3D<T>(*this);
}

template<typename T>
void BoxNTensorMaxFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;
}

template<typename T>
std::vector<T> BoxNTensorMaxFunctional3D<T>::getMaxVector() const {
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
MaskedBoxNTensorMaxFunctional3D<T>::MaskedBoxNTensorMaxFunctional3D(plint ndim)
    : maxVectorId(ndim)
{
    for (plint iDim=0; iDim<ndim; ++iDim) {
        maxVectorId[iDim] = this->getStatistics().subscribeMax();
    }
}

template<typename T>
void MaskedBoxNTensorMaxFunctional3D<T>::process (
        Box3D domain, NTensorField3D<T>& vectorField,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( vectorField.getNdim() == (plint)maxVectorId.size() );
    Dot3D maskOfs = computeRelativeDisplacement(vectorField, mask);
    plint ndim = vectorField.getNdim();
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        statistics.gatherMax(maxVectorId[iDim], (double)vectorField.get(iX,iY,iZ)[iDim]);
                    }
                }
            }
        }
    }
}

template<typename T>
MaskedBoxNTensorMaxFunctional3D<T>* MaskedBoxNTensorMaxFunctional3D<T>::clone() const
{
    return new MaskedBoxNTensorMaxFunctional3D<T>(*this);
}

template<typename T>
void MaskedBoxNTensorMaxFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
}

template<typename T>
std::vector<T> MaskedBoxNTensorMaxFunctional3D<T>::getMaxVector() const {
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
BoundedBoxNTensorSumFunctional3D<T>::BoundedBoxNTensorSumFunctional3D(plint ndim)
    : sumVectorId(ndim)
{
    for (plint iDim=0; iDim<ndim; ++iDim) {
        sumVectorId[iDim] = this->getStatistics().subscribeSum();
    }
}

template<typename T>
void BoundedBoxNTensorSumFunctional3D<T>::processBulk (
        Box3D domain, NTensorField3D<T>& vectorField )
{
    PLB_PRECONDITION( vectorField.getNdim() == (plint)sumVectorId.size() );
    plint ndim = vectorField.getNdim();
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    statistics.gatherSum(sumVectorId[iDim], (double)vectorField.get(iX,iY,iZ)[iDim]);
                }
            }
        }
    }
}

template<typename T>
void BoundedBoxNTensorSumFunctional3D<T>::processPlane (
        int direction, int orientation,
        Box3D domain, NTensorField3D<T>& vectorField )
{
    PLB_PRECONDITION( vectorField.getNdim() == (plint)sumVectorId.size() );
    plint ndim = vectorField.getNdim();
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                // Plane nodes have a weight of 0.5, because only 50% of the
                //   cell centered at the node is inside the computational domain.
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    statistics.gatherSum( sumVectorId[iDim],
                                          (double)vectorField.get(iX,iY,iZ)[iDim] / 2.);
                }
            }
        }
    }
}

template<typename T>
void BoundedBoxNTensorSumFunctional3D<T>::processEdge (
        int plane, int normal1, int normal2,
        Box3D domain, NTensorField3D<T>& vectorField )
{
    PLB_PRECONDITION( vectorField.getNdim() == (plint)sumVectorId.size() );
    plint ndim = vectorField.getNdim();
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                // Edge nodes have a weight of 0.25, because only 25% of the
                //   cell centered at the node is inside the computational domain.
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    statistics.gatherSum( sumVectorId[iDim],
                                          (double)vectorField.get(iX,iY,iZ)[iDim] / 4.);
                }
            }
        }
    }
}

template<typename T>
void BoundedBoxNTensorSumFunctional3D<T>::processCorner (
        int normalX, int normalY, int normalZ,
        Box3D domain, NTensorField3D<T>& vectorField )
{
    PLB_PRECONDITION( vectorField.getNdim() == (plint)sumVectorId.size() );
    plint ndim = vectorField.getNdim();
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                // Corner nodes have a weight of 0.125, because only 1/8 of the
                //   cell centered at the node is inside the computational domain.
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    statistics.gatherSum( sumVectorId[iDim],
                                          (double)vectorField.get(iX,iY,iZ)[iDim] / 8.);
                }
            }
        }
    }
}

template<typename T>
BoundedBoxNTensorSumFunctional3D<T>* BoundedBoxNTensorSumFunctional3D<T>::clone() const
{
    return new BoundedBoxNTensorSumFunctional3D<T>(*this);
}

template<typename T>
void BoundedBoxNTensorSumFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;
}

template<typename T>
std::vector<T> BoundedBoxNTensorSumFunctional3D<T>::getSumVector() const {
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
BoundedMaskedBoxNTensorSumFunctional3D<T>::BoundedMaskedBoxNTensorSumFunctional3D(plint ndim)
    : sumVectorId(ndim)
{
    for (plint iDim=0; iDim<ndim; ++iDim) {
        sumVectorId[iDim] = this->getStatistics().subscribeSum();
    }
}

template<typename T>
void BoundedMaskedBoxNTensorSumFunctional3D<T>::processBulk (
        Box3D domain, NTensorField3D<T>& vectorField, NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( vectorField.getNdim() == (plint)sumVectorId.size() );
    Dot3D maskOfs = computeRelativeDisplacement(vectorField, mask);
    plint ndim = vectorField.getNdim();
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        statistics.gatherSum(sumVectorId[iDim], (double)vectorField.get(iX,iY,iZ)[iDim]);
                    }
                }
            }
        }
    }
}

template<typename T>
void BoundedMaskedBoxNTensorSumFunctional3D<T>::processPlane (
        int direction, int orientation,
        Box3D domain, NTensorField3D<T>& vectorField, NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( vectorField.getNdim() == (plint)sumVectorId.size() );
    Dot3D maskOfs = computeRelativeDisplacement(vectorField, mask);
    plint ndim = vectorField.getNdim();
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                    // Plane nodes have a weight of 0.5, because only 50% of the
                    //   cell centered at the node is inside the computational domain.
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        statistics.gatherSum( sumVectorId[iDim],
                                              (double)vectorField.get(iX,iY,iZ)[iDim] / 2.);
                    }
                }
            }
        }
    }
}

template<typename T>
void BoundedMaskedBoxNTensorSumFunctional3D<T>::processEdge (
        int plane, int normal1, int normal2,
        Box3D domain, NTensorField3D<T>& vectorField, NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( vectorField.getNdim() == (plint)sumVectorId.size() );
    Dot3D maskOfs = computeRelativeDisplacement(vectorField, mask);
    plint ndim = vectorField.getNdim();
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                    // Edge nodes have a weight of 0.25, because only 25% of the
                    //   cell centered at the node is inside the computational domain.
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        statistics.gatherSum( sumVectorId[iDim],
                                              (double)vectorField.get(iX,iY,iZ)[iDim] / 4.);
                    }
                }
            }
        }
    }
}

template<typename T>
void BoundedMaskedBoxNTensorSumFunctional3D<T>::processCorner (
        int normalX, int normalY, int normalZ,
        Box3D domain, NTensorField3D<T>& vectorField, NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( vectorField.getNdim() == (plint)sumVectorId.size() );
    Dot3D maskOfs = computeRelativeDisplacement(vectorField, mask);
    plint ndim = vectorField.getNdim();
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                    // Corner nodes have a weight of 0.125, because only 1/8 of the
                    //   cell centered at the node is inside the computational domain.
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        statistics.gatherSum( sumVectorId[iDim],
                                              (double)vectorField.get(iX,iY,iZ)[iDim] / 8.);
                    }
                }
            }
        }
    }
}

template<typename T>
BoundedMaskedBoxNTensorSumFunctional3D<T>* BoundedMaskedBoxNTensorSumFunctional3D<T>::clone() const
{
    return new BoundedMaskedBoxNTensorSumFunctional3D<T>(*this);
}

template<typename T>
void BoundedMaskedBoxNTensorSumFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
}

template<typename T>
std::vector<T> BoundedMaskedBoxNTensorSumFunctional3D<T>::getSumVector() const {
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
void MaskedCopyConvertNTensorFunctional3D<T1,T2>::process (
        Box3D domain, NTensorField3D<T1>& field1,
                      NTensorField3D<T2>& field2,
                      NTensorField3D<int>& mask )
{

    PLB_PRECONDITION( field1.getNdim() == field2.getNdim());
    Dot3D maskOfs = computeRelativeDisplacement(field1, mask);
    plint ndim = field1.getNdim();
    Dot3D offset = computeRelativeDisplacement(field1, field2);

    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (int iDim=0; iDim<ndim; ++iDim) {
                            *field2.get(iX+offset.x,iY+offset.y,iZ+offset.z) =
                                (T2) *field1.get(iX,iY,iZ);
                        }
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (int iDim=0; iDim<ndim; ++iDim) {
                            field2.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim] =
                                (T2) field1.get(iX,iY,iZ)[iDim];
                        }
                    }
                }
            }
        }
    }
}

template<typename T1, typename T2>
MaskedCopyConvertNTensorFunctional3D<T1,T2>* MaskedCopyConvertNTensorFunctional3D<T1,T2>::clone() const
{
    return new MaskedCopyConvertNTensorFunctional3D<T1,T2>(*this);
}

template<typename T1, typename T2>
void MaskedCopyConvertNTensorFunctional3D<T1,T2>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T1, typename T2>
BlockDomain::DomainT MaskedCopyConvertNTensorFunctional3D<T1,T2>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T>
ExtractNTensorComponentFunctional3D<T>::ExtractNTensorComponentFunctional3D(int iComponent_)
    : iComponent(iComponent_)
{ }

template<typename T>
void ExtractNTensorComponentFunctional3D<T>::process (
        Box3D domain, NTensorField3D<T>& scalarField,
                      NTensorField3D<T>& tensorField )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    PLB_PRECONDITION( tensorField.getNdim()>iComponent );
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                *scalarField.get(iX,iY,iZ) = tensorField.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iComponent];
            }
        }
    }
}

template<typename T>
ExtractNTensorComponentFunctional3D<T>* ExtractNTensorComponentFunctional3D<T>::clone() const
{
    return new ExtractNTensorComponentFunctional3D<T>(*this);
}

template<typename T>
void ExtractNTensorComponentFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT ExtractNTensorComponentFunctional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T>
MaskedExtractNTensorComponentFunctional3D<T>::MaskedExtractNTensorComponentFunctional3D(int iComponent_)
    : iComponent(iComponent_)
{ }

template<typename T>
void MaskedExtractNTensorComponentFunctional3D<T>::process (
        Box3D domain, NTensorField3D<T>& scalarField,
                      NTensorField3D<T>& tensorField,
                      NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    PLB_PRECONDITION( tensorField.getNdim()>iComponent );
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    Dot3D maskOfs = computeRelativeDisplacement(scalarField, mask);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                    *scalarField.get(iX,iY,iZ) = tensorField.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iComponent];
                }
            }
        }
    }
}

template<typename T>
MaskedExtractNTensorComponentFunctional3D<T>* MaskedExtractNTensorComponentFunctional3D<T>::clone() const
{
    return new MaskedExtractNTensorComponentFunctional3D<T>(*this);
}

template<typename T>
void MaskedExtractNTensorComponentFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT MaskedExtractNTensorComponentFunctional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T>
void ComputeNTensorNormFunctional3D<T>::process (
        Box3D domain, NTensorField3D<T>& scalarField,
                      NTensorField3D<T>& tensorField )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    plint ndim = tensorField.getNdim();
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    T normSqr = T();
                    T* vector = tensorField.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        normSqr += vector[iDim]*vector[iDim];
                    }
                    *scalarField.get(iX,iY,iZ) = std::sqrt(normSqr);
                }
            }
        }
    }
}

template<typename T>
ComputeNTensorNormFunctional3D<T>* ComputeNTensorNormFunctional3D<T>::clone() const
{
    return new ComputeNTensorNormFunctional3D<T>(*this);
}

template<typename T>
void ComputeNTensorNormFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT ComputeNTensorNormFunctional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T>
void MaskedComputeNTensorNormFunctional3D<T>::process (
        Box3D domain, NTensorField3D<T>& scalarField,
                      NTensorField3D<T>& tensorField,
                      NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    Dot3D maskOfs = computeRelativeDisplacement(scalarField, mask);
    plint ndim = tensorField.getNdim();
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        T normSqr = T();
                        T* vector = tensorField.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            normSqr += vector[iDim]*vector[iDim];
                        }
                        *scalarField.get(iX,iY,iZ) = std::sqrt(normSqr);
                    }
                }
            }
        }
    }
}

template<typename T>
MaskedComputeNTensorNormFunctional3D<T>* MaskedComputeNTensorNormFunctional3D<T>::clone() const
{
    return new MaskedComputeNTensorNormFunctional3D<T>(*this);
}

template<typename T>
void MaskedComputeNTensorNormFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT MaskedComputeNTensorNormFunctional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T>
void ComputeNTensorNormSqrFunctional3D<T>::process (
        Box3D domain, NTensorField3D<T>& scalarField,
                      NTensorField3D<T>& tensorField )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    plint ndim = tensorField.getNdim();
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    T normSqr = T();
                    T* vector = tensorField.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        normSqr += vector[iDim]*vector[iDim];
                    }
                    *scalarField.get(iX,iY,iZ) = normSqr;
                }
            }
        }
    }
}

template<typename T>
ComputeNTensorNormSqrFunctional3D<T>* ComputeNTensorNormSqrFunctional3D<T>::clone() const
{
    return new ComputeNTensorNormSqrFunctional3D<T>(*this);
}

template<typename T>
void ComputeNTensorNormSqrFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT ComputeNTensorNormSqrFunctional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T>
void MaskedComputeNTensorNormSqrFunctional3D<T>::process (
        Box3D domain, NTensorField3D<T>& scalarField,
                      NTensorField3D<T>& tensorField,
                      NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    Dot3D maskOfs = computeRelativeDisplacement(scalarField, mask);
    plint ndim = tensorField.getNdim();
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        T normSqr = T();
                        T* vector = tensorField.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            normSqr += vector[iDim]*vector[iDim];
                        }
                        *scalarField.get(iX,iY,iZ) = normSqr;
                    }
                }
            }
        }
    }
}

template<typename T>
MaskedComputeNTensorNormSqrFunctional3D<T>* MaskedComputeNTensorNormSqrFunctional3D<T>::clone() const
{
    return new MaskedComputeNTensorNormSqrFunctional3D<T>(*this);
}

template<typename T>
void MaskedComputeNTensorNormSqrFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT MaskedComputeNTensorNormSqrFunctional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T>
void ComputeSymmetricNTensorNormFunctional3D<T>::process (
        Box3D domain, NTensorField3D<T>& scalarField,
                      NTensorField3D<T>& tensorField )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    typedef SymmetricTensorImpl<T,3> tensor;
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                T* el = tensorField.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                *scalarField.get(iX,iY,iZ) = std::sqrt ( 
                        // Count diagonal components once ...
                                util::sqr(el[tensor::xx]) + util::sqr(el[tensor::yy]) + util::sqr(el[tensor::zz]) +
                        // .. and off-diagonal component twice, due to symmetry.
                        (T)2 * (util::sqr(el[tensor::xy]) + util::sqr(el[tensor::xz]) + util::sqr(el[tensor::yz]) ) );
            }
        }
    }
}

template<typename T>
ComputeSymmetricNTensorNormFunctional3D<T>* ComputeSymmetricNTensorNormFunctional3D<T>::clone() const
{
    return new ComputeSymmetricNTensorNormFunctional3D<T>(*this);
}

template<typename T>
void ComputeSymmetricNTensorNormFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT ComputeSymmetricNTensorNormFunctional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T>
void MaskedComputeSymmetricNTensorNormFunctional3D<T>::process (
        Box3D domain, NTensorField3D<T>& scalarField,
                      NTensorField3D<T>& tensorField,
                      NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    Dot3D maskOfs = computeRelativeDisplacement(scalarField, mask);
    typedef SymmetricTensorImpl<T,3> tensor;
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                    T* el = tensorField.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                    *scalarField.get(iX,iY,iZ) = std::sqrt ( 
                            // Count diagonal components once ...
                                    util::sqr(el[tensor::xx]) + util::sqr(el[tensor::yy]) + util::sqr(el[tensor::zz]) +
                            // .. and off-diagonal component twice, due to symmetry.
                            (T)2 * (util::sqr(el[tensor::xy]) + util::sqr(el[tensor::xz]) + util::sqr(el[tensor::yz]) ) );
                }
            }
        }
    }
}

template<typename T>
MaskedComputeSymmetricNTensorNormFunctional3D<T>* MaskedComputeSymmetricNTensorNormFunctional3D<T>::clone() const
{
    return new MaskedComputeSymmetricNTensorNormFunctional3D<T>(*this);
}

template<typename T>
void MaskedComputeSymmetricNTensorNormFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT MaskedComputeSymmetricNTensorNormFunctional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T>
void ComputeSymmetricNTensorNormSqrFunctional3D<T>::process (
        Box3D domain, NTensorField3D<T>& scalarField,
                      NTensorField3D<T>& tensorField )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    typedef SymmetricTensorImpl<T,3> tensor;
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                T* el = tensorField.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                *scalarField.get(iX,iY,iZ) =
                        // Count diagonal components once ...
                                util::sqr(el[tensor::xx]) + util::sqr(el[tensor::yy]) + util::sqr(el[tensor::zz]) +
                        // .. and off-diagonal component twice, due to symmetry.
                        (T)2 * (util::sqr(el[tensor::xy]) + util::sqr(el[tensor::xz]) + util::sqr(el[tensor::yz]) );
            }
        }
    }
}

template<typename T>
ComputeSymmetricNTensorNormSqrFunctional3D<T>* ComputeSymmetricNTensorNormSqrFunctional3D<T>::clone() const
{
    return new ComputeSymmetricNTensorNormSqrFunctional3D<T>(*this);
}

template<typename T>
void ComputeSymmetricNTensorNormSqrFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT ComputeSymmetricNTensorNormSqrFunctional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

template<typename T>
void MaskedComputeSymmetricNTensorNormSqrFunctional3D<T>::process (
        Box3D domain, NTensorField3D<T>& scalarField,
                      NTensorField3D<T>& tensorField,
                      NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    Dot3D maskOfs = computeRelativeDisplacement(scalarField, mask);
    typedef SymmetricTensorImpl<T,3> tensor;
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                    T* el = tensorField.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                    *scalarField.get(iX,iY,iZ) =
                            // Count diagonal components once ...
                                    util::sqr(el[tensor::xx]) + util::sqr(el[tensor::yy]) + util::sqr(el[tensor::zz]) +
                            // .. and off-diagonal component twice, due to symmetry.
                            (T)2 * (util::sqr(el[tensor::xy]) + util::sqr(el[tensor::xz]) + util::sqr(el[tensor::yz]) );
                }
            }
        }
    }
}

template<typename T>
MaskedComputeSymmetricNTensorNormSqrFunctional3D<T>* MaskedComputeSymmetricNTensorNormSqrFunctional3D<T>::clone() const
{
    return new MaskedComputeSymmetricNTensorNormSqrFunctional3D<T>(*this);
}

template<typename T>
void MaskedComputeSymmetricNTensorNormSqrFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT MaskedComputeSymmetricNTensorNormSqrFunctional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

template<typename T>
void ComputeSymmetricNTensorTraceFunctional3D<T>::process (
        Box3D domain, NTensorField3D<T>& scalarField,
                      NTensorField3D<T>& tensorField )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    typedef SymmetricTensorImpl<T,3> tensor;
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                T* el = tensorField.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                *scalarField.get(iX,iY,iZ) = el[tensor::xx] + el[tensor::yy] + el[tensor::zz];
            }
        }
    }
}

template<typename T>
ComputeSymmetricNTensorTraceFunctional3D<T>* ComputeSymmetricNTensorTraceFunctional3D<T>::clone() const
{
    return new ComputeSymmetricNTensorTraceFunctional3D<T>(*this);
}

template<typename T>
void ComputeSymmetricNTensorTraceFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT ComputeSymmetricNTensorTraceFunctional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T>
void MaskedComputeSymmetricNTensorTraceFunctional3D<T>::process (
        Box3D domain, NTensorField3D<T>& scalarField,
                      NTensorField3D<T>& tensorField,
                      NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( scalarField.getNdim()==1 );
    Dot3D maskOfs = computeRelativeDisplacement(scalarField, mask);
    typedef SymmetricTensorImpl<T,3> tensor;
    Dot3D offset = computeRelativeDisplacement(scalarField, tensorField);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                    T* el = tensorField.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                    *scalarField.get(iX,iY,iZ) = el[tensor::xx] + el[tensor::yy] + el[tensor::zz];
                }
            }
        }
    }
}

template<typename T>
MaskedComputeSymmetricNTensorTraceFunctional3D<T>* MaskedComputeSymmetricNTensorTraceFunctional3D<T>::clone() const
{
    return new MaskedComputeSymmetricNTensorTraceFunctional3D<T>(*this);
}

template<typename T>
void MaskedComputeSymmetricNTensorTraceFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT MaskedComputeSymmetricNTensorTraceFunctional3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


template<typename T>
void BoxBulkNTensorVorticityFunctional3D<T>::process (
        Box3D domain, NTensorField3D<T>& vorticity,
                      NTensorField3D<T>& velocity )
{
    PLB_PRECONDITION( vorticity.getNdim()==3 );
    PLB_PRECONDITION( velocity.getNdim()==3 );
    Dot3D offset = computeRelativeDisplacement(vorticity, velocity);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                plint iZ2 = iZ+offset.z;
                vorticity.get(iX,iY,iZ)[0] =
                    fdNTensorField::bulkVorticityX(velocity, iX2,iY2,iZ2);
                vorticity.get(iX,iY,iZ)[1] =
                    fdNTensorField::bulkVorticityY(velocity, iX2,iY2,iZ2);
                vorticity.get(iX,iY,iZ)[2] =
                    fdNTensorField::bulkVorticityZ(velocity, iX2,iY2,iZ2);
            }
        }
    }
}

template<typename T>
BoxBulkNTensorVorticityFunctional3D<T>* BoxBulkNTensorVorticityFunctional3D<T>::clone() const
{
    return new BoxBulkNTensorVorticityFunctional3D<T>(*this);
}

template<typename T>
void BoxBulkNTensorVorticityFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT BoxBulkNTensorVorticityFunctional3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}



template<typename T>
void MaskedBoxBulkNTensorVorticityFunctional3D<T>::process (
        Box3D domain, NTensorField3D<T>& vorticity,
                      NTensorField3D<T>& velocity,
                      NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( vorticity.getNdim()==1 );
    PLB_PRECONDITION( velocity.getNdim()==2 );
    Dot3D maskOfs = computeRelativeDisplacement(vorticity, mask);
    Dot3D offset = computeRelativeDisplacement(vorticity, velocity);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                    plint iX2 = iX+offset.x;
                    plint iY2 = iY+offset.y;
                    plint iZ2 = iZ+offset.z;
                    vorticity.get(iX,iY,iZ)[0] =
                        fdNTensorField::bulkVorticityX(velocity, iX2,iY2,iZ2);
                    vorticity.get(iX,iY,iZ)[1] =
                        fdNTensorField::bulkVorticityY(velocity, iX2,iY2,iZ2);
                    vorticity.get(iX,iY,iZ)[2] =
                        fdNTensorField::bulkVorticityZ(velocity, iX2,iY2,iZ2);
                }
            }
        }
    }
}

template<typename T>
MaskedBoxBulkNTensorVorticityFunctional3D<T>* MaskedBoxBulkNTensorVorticityFunctional3D<T>::clone() const
{
    return new MaskedBoxBulkNTensorVorticityFunctional3D<T>(*this);
}

template<typename T>
void MaskedBoxBulkNTensorVorticityFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT MaskedBoxBulkNTensorVorticityFunctional3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}



template<typename T>
void BoxNTensorVorticityFunctional3D<T>::processBulk (
        Box3D domain, NTensorField3D<T>& velocity, NTensorField3D<T>& vorticity )
{
    PLB_PRECONDITION( velocity.getNdim()==3 );
    PLB_PRECONDITION( vorticity.getNdim()==3 );

    Dot3D offset = computeRelativeDisplacement(velocity, vorticity);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                plint iZ2 = iZ+offset.z;
                vorticity.get(iX,iY,iZ)[0] =
                    fdNTensorField::bulkVorticityX(velocity, iX2,iY2,iZ2);
                vorticity.get(iX,iY,iZ)[1] =
                    fdNTensorField::bulkVorticityY(velocity, iX2,iY2,iZ2);
                vorticity.get(iX,iY,iZ)[2] =
                    fdNTensorField::bulkVorticityZ(velocity, iX2,iY2,iZ2);
            }
        }
    }
}

template<typename T>
void BoxNTensorVorticityFunctional3D<T>::processPlane (
        int direction, int orientation, Box3D domain,
        NTensorField3D<T>& velocity, NTensorField3D<T>& vorticity )
{
    PLB_PRECONDITION( velocity.getNdim()==3 );
    PLB_PRECONDITION( vorticity.getNdim()==3 );

    Dot3D offset = computeRelativeDisplacement(velocity, vorticity);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                plint iZ2 = iZ+offset.z;
                vorticity.get(iX,iY,iZ)[0] =
                    fdNTensorField::planeVorticityX(velocity,direction,orientation, iX2,iY2,iZ2);
                vorticity.get(iX,iY,iZ)[1] =
                    fdNTensorField::planeVorticityY(velocity,direction,orientation, iX2,iY2,iZ2);
                vorticity.get(iX,iY,iZ)[2] =
                    fdNTensorField::planeVorticityZ(velocity,direction,orientation, iX2,iY2,iZ2);
            }
        }
    }
}

template<typename T>
void BoxNTensorVorticityFunctional3D<T>::processEdge (
        int plane, int normal1, int normal2, Box3D domain,
        NTensorField3D<T>& velocity, NTensorField3D<T>& vorticity )
{
    PLB_PRECONDITION( velocity.getNdim()==3 );
    PLB_PRECONDITION( vorticity.getNdim()==3 );

    Dot3D offset = computeRelativeDisplacement(velocity, vorticity);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                plint iZ2 = iZ+offset.z;
                vorticity.get(iX,iY,iZ)[0] =
                    fdNTensorField::edgeVorticityX(velocity,plane,normal1,normal2, iX2,iY2,iZ2);
                vorticity.get(iX,iY,iZ)[1] =
                    fdNTensorField::edgeVorticityY(velocity,plane,normal1,normal2, iX2,iY2,iZ2);
                vorticity.get(iX,iY,iZ)[2] =
                    fdNTensorField::edgeVorticityZ(velocity,plane,normal1,normal2, iX2,iY2,iZ2);
            }
        }
    }
}

template<typename T>
void BoxNTensorVorticityFunctional3D<T>::processCorner (
        int normalX, int normalY, int normalZ, Box3D domain,
        NTensorField3D<T>& velocity, NTensorField3D<T>& vorticity )
{
    PLB_PRECONDITION( velocity.getNdim()==3 );
    PLB_PRECONDITION( vorticity.getNdim()==3 );

    Dot3D offset = computeRelativeDisplacement(velocity, vorticity);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                plint iZ2 = iZ+offset.z;
                vorticity.get(iX,iY,iZ)[0] =
                    fdNTensorField::cornerVorticityX(velocity,normalX,normalY,normalZ, iX2,iY2,iZ2);
                vorticity.get(iX,iY,iZ)[1] =
                    fdNTensorField::cornerVorticityY(velocity,normalX,normalY,normalZ, iX2,iY2,iZ2);
                vorticity.get(iX,iY,iZ)[2] =
                    fdNTensorField::cornerVorticityZ(velocity,normalX,normalY,normalZ, iX2,iY2,iZ2);
            }
        }
    }
}


template<typename T>
BoxNTensorVorticityFunctional3D<T>* BoxNTensorVorticityFunctional3D<T>::clone() const
{
    return new BoxNTensorVorticityFunctional3D<T>(*this);
}

template<typename T>
void BoxNTensorVorticityFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT BoxNTensorVorticityFunctional3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}




template<typename T>
void MaskedBoxNTensorVorticityFunctional3D<T>::processBulk (
        Box3D domain, NTensorField3D<T>& velocity, NTensorField3D<T>& vorticity,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( velocity.getNdim()==3 );
    PLB_PRECONDITION( vorticity.getNdim()==3 );

    Dot3D offset = computeRelativeDisplacement(velocity, vorticity);
    Dot3D maskOfs = computeRelativeDisplacement(velocity, mask);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                    plint iX2 = iX+offset.x;
                    plint iY2 = iY+offset.y;
                    plint iZ2 = iZ+offset.z;
                    vorticity.get(iX,iY,iZ)[0] =
                        fdNTensorField::bulkVorticityX(velocity, iX2,iY2,iZ2);
                    vorticity.get(iX,iY,iZ)[1] =
                        fdNTensorField::bulkVorticityY(velocity, iX2,iY2,iZ2);
                    vorticity.get(iX,iY,iZ)[2] =
                        fdNTensorField::bulkVorticityZ(velocity, iX2,iY2,iZ2);
                }
            }
        }
    }
}

template<typename T>
void MaskedBoxNTensorVorticityFunctional3D<T>::processPlane (
        int direction, int orientation, Box3D domain,
        NTensorField3D<T>& velocity, NTensorField3D<T>& vorticity,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( velocity.getNdim()==3 );
    PLB_PRECONDITION( vorticity.getNdim()==3 );

    Dot3D offset = computeRelativeDisplacement(velocity, vorticity);
    Dot3D maskOfs = computeRelativeDisplacement(velocity, mask);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                    plint iX2 = iX+offset.x;
                    plint iY2 = iY+offset.y;
                    plint iZ2 = iZ+offset.z;
                    vorticity.get(iX,iY,iZ)[0] =
                        fdNTensorField::planeVorticityX(velocity,direction,orientation, iX2,iY2,iZ2);
                    vorticity.get(iX,iY,iZ)[1] =
                        fdNTensorField::planeVorticityY(velocity,direction,orientation, iX2,iY2,iZ2);
                    vorticity.get(iX,iY,iZ)[2] =
                        fdNTensorField::planeVorticityZ(velocity,direction,orientation, iX2,iY2,iZ2);
                }
            }
        }
    }
}

template<typename T>
void MaskedBoxNTensorVorticityFunctional3D<T>::processEdge (
        int plane, int normal1, int normal2, Box3D domain,
        NTensorField3D<T>& velocity, NTensorField3D<T>& vorticity,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( velocity.getNdim()==3 );
    PLB_PRECONDITION( vorticity.getNdim()==3 );

    Dot3D offset = computeRelativeDisplacement(velocity, vorticity);
    Dot3D maskOfs = computeRelativeDisplacement(velocity, mask);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                    plint iX2 = iX+offset.x;
                    plint iY2 = iY+offset.y;
                    plint iZ2 = iZ+offset.z;
                    vorticity.get(iX,iY,iZ)[0] =
                        fdNTensorField::edgeVorticityX(velocity,plane,normal1,normal2, iX2,iY2,iZ2);
                    vorticity.get(iX,iY,iZ)[1] =
                        fdNTensorField::edgeVorticityY(velocity,plane,normal1,normal2, iX2,iY2,iZ2);
                    vorticity.get(iX,iY,iZ)[2] =
                        fdNTensorField::edgeVorticityZ(velocity,plane,normal1,normal2, iX2,iY2,iZ2);
                }
            }
        }
    }
}

template<typename T>
void MaskedBoxNTensorVorticityFunctional3D<T>::processCorner (
        int normalX, int normalY, int normalZ, Box3D domain,
        NTensorField3D<T>& velocity, NTensorField3D<T>& vorticity,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( velocity.getNdim()==3 );
    PLB_PRECONDITION( vorticity.getNdim()==3 );

    Dot3D offset = computeRelativeDisplacement(velocity, vorticity);
    Dot3D maskOfs = computeRelativeDisplacement(velocity, mask);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                    plint iX2 = iX+offset.x;
                    plint iY2 = iY+offset.y;
                    plint iZ2 = iZ+offset.z;
                    vorticity.get(iX,iY,iZ)[0] =
                        fdNTensorField::cornerVorticityX(velocity,normalX,normalY,normalZ, iX2,iY2,iZ2);
                    vorticity.get(iX,iY,iZ)[1] =
                        fdNTensorField::cornerVorticityY(velocity,normalX,normalY,normalZ, iX2,iY2,iZ2);
                    vorticity.get(iX,iY,iZ)[2] =
                        fdNTensorField::cornerVorticityZ(velocity,normalX,normalY,normalZ, iX2,iY2,iZ2);
                }
            }
        }
    }
}


template<typename T>
MaskedBoxNTensorVorticityFunctional3D<T>* MaskedBoxNTensorVorticityFunctional3D<T>::clone() const
{
    return new MaskedBoxNTensorVorticityFunctional3D<T>(*this);
}

template<typename T>
void MaskedBoxNTensorVorticityFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}


template<typename T>
BlockDomain::DomainT MaskedBoxNTensorVorticityFunctional3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}



template<typename T>
void BoxBulkNTensorStrainRateFunctional3D<T>::process (
        Box3D domain, NTensorField3D<T>& velocity,
                      NTensorField3D<T>& S )
{
    PLB_PRECONDITION( velocity.getNdim()==2 );
    PLB_PRECONDITION( S.getNdim()==3 );
    typedef SymmetricTensorImpl<T,3> tensor;
    Dot3D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                plint iZ2 = iZ+offset.z;
                T* el = S.get(iX2,iY2,iZ2);
                el[tensor::xx] = fdNTensorField::bulkXderiv(velocity, iX, iY, iZ, 0);
                el[tensor::yy] = fdNTensorField::bulkYderiv(velocity, iX, iY, iZ, 1);
                el[tensor::zz] = fdNTensorField::bulkZderiv(velocity, iX, iY, iZ, 2);
                el[tensor::xy] = ( fdNTensorField::bulkXderiv(velocity, iX, iY, iZ, 1) +
                                   fdNTensorField::bulkYderiv(velocity, iX, iY, iZ, 0) ) / (T)2;
                el[tensor::xz] = ( fdNTensorField::bulkXderiv(velocity, iX, iY, iZ, 2) +
                                   fdNTensorField::bulkZderiv(velocity, iX, iY, iZ, 0) ) / (T)2;
                el[tensor::yz] = ( fdNTensorField::bulkYderiv(velocity, iX, iY, iZ, 2) +
                                   fdNTensorField::bulkZderiv(velocity, iX, iY, iZ, 1) ) / (T)2;
            }
        }
    }
}

template<typename T>
BoxBulkNTensorStrainRateFunctional3D<T>* BoxBulkNTensorStrainRateFunctional3D<T>::clone() const
{
    return new BoxBulkNTensorStrainRateFunctional3D<T>(*this);
}

template<typename T>
void BoxBulkNTensorStrainRateFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT BoxBulkNTensorStrainRateFunctional3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}



template<typename T>
void MaskedBoxBulkNTensorStrainRateFunctional3D<T>::process (
        Box3D domain, NTensorField3D<T>& velocity,
                      NTensorField3D<T>& S,
                      NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( velocity.getNdim()==2 );
    PLB_PRECONDITION( S.getNdim()==3 );
    Dot3D maskOfs = computeRelativeDisplacement(velocity, mask);
    typedef SymmetricTensorImpl<T,3> tensor;
    Dot3D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                    plint iX2 = iX+offset.x;
                    plint iY2 = iY+offset.y;
                    plint iZ2 = iZ+offset.z;
                    T* el = S.get(iX2,iY2,iZ2);
                    el[tensor::xx] = fdNTensorField::bulkXderiv(velocity, iX, iY, iZ, 0);
                    el[tensor::yy] = fdNTensorField::bulkYderiv(velocity, iX, iY, iZ, 1);
                    el[tensor::zz] = fdNTensorField::bulkZderiv(velocity, iX, iY, iZ, 2);
                    el[tensor::xy] = ( fdNTensorField::bulkXderiv(velocity, iX, iY, iZ, 1) +
                                       fdNTensorField::bulkYderiv(velocity, iX, iY, iZ, 0) ) / (T)2;
                    el[tensor::xz] = ( fdNTensorField::bulkXderiv(velocity, iX, iY, iZ, 2) +
                                       fdNTensorField::bulkZderiv(velocity, iX, iY, iZ, 0) ) / (T)2;
                    el[tensor::yz] = ( fdNTensorField::bulkYderiv(velocity, iX, iY, iZ, 2) +
                                       fdNTensorField::bulkZderiv(velocity, iX, iY, iZ, 1) ) / (T)2;
                }
            }
        }
    }
}

template<typename T>
MaskedBoxBulkNTensorStrainRateFunctional3D<T>* MaskedBoxBulkNTensorStrainRateFunctional3D<T>::clone() const
{
    return new MaskedBoxBulkNTensorStrainRateFunctional3D<T>(*this);
}

template<typename T>
void MaskedBoxBulkNTensorStrainRateFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT MaskedBoxBulkNTensorStrainRateFunctional3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


template<typename T>
void BoxNTensorStrainRateFunctional3D<T>::processBulk (
        Box3D domain, NTensorField3D<T>& velocity, NTensorField3D<T>& S )
{
    typedef SymmetricTensorImpl<T,3> tensor;
    PLB_PRECONDITION( velocity.getNdim()==3 );
    PLB_PRECONDITION( S.getNdim()==6 );
    Dot3D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                plint iZ2 = iZ+offset.z;
                T* el = S.get(iX2,iY2,iZ2);
                el[tensor::xx] = fdNTensorField::bulkXderiv(velocity, iX, iY, iZ, 0);
                el[tensor::xy] = ( fdNTensorField::bulkXderiv(velocity, iX, iY, iZ, 1) +
                                   fdNTensorField::bulkYderiv(velocity, iX, iY, iZ, 0) ) / (T)2;
                el[tensor::xz] = ( fdNTensorField::bulkXderiv(velocity, iX, iY, iZ, 2) +
                                   fdNTensorField::bulkZderiv(velocity, iX, iY, iZ, 0) ) / (T)2;
                el[tensor::yy] = fdNTensorField::bulkYderiv(velocity, iX, iY, iZ, 1);
                el[tensor::yz] = ( fdNTensorField::bulkYderiv(velocity, iX, iY, iZ, 2) +
                                   fdNTensorField::bulkZderiv(velocity, iX, iY, iZ, 1) ) / (T)2;
                el[tensor::zz] = fdNTensorField::bulkZderiv(velocity, iX, iY, iZ, 2);
            }
        }
    }
}

template<typename T>
void BoxNTensorStrainRateFunctional3D<T>::processPlane (
        int direction, int orientation, Box3D domain,
        NTensorField3D<T>& velocity, NTensorField3D<T>& S )
{
    PLB_PRECONDITION( velocity.getNdim()==3 );
    PLB_PRECONDITION( S.getNdim()==6 );
    typedef SymmetricTensorImpl<T,3> tensor;
    Dot3D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                plint iZ2 = iZ+offset.z;
                T* el = S.get(iX2,iY2,iZ2);
                el[tensor::xx] = fdNTensorField::planeXderiv(velocity, direction,orientation, iX, iY, iZ, 0);
                el[tensor::xy] = ( fdNTensorField::planeXderiv(velocity, direction,orientation, iX, iY, iZ, 1) +
                               fdNTensorField::planeYderiv(velocity, direction,orientation, iX, iY, iZ, 0) ) / (T)2;
                el[tensor::xz] = ( fdNTensorField::planeXderiv(velocity, direction,orientation, iX, iY, iZ, 2) +
                               fdNTensorField::planeZderiv(velocity, direction,orientation, iX, iY, iZ, 0) ) / (T)2;
                el[tensor::yy] = fdNTensorField::planeYderiv(velocity, direction,orientation, iX, iY, iZ, 1);
                el[tensor::yz] = ( fdNTensorField::planeYderiv(velocity, direction,orientation, iX, iY, iZ, 2) +
                               fdNTensorField::planeZderiv(velocity, direction,orientation, iX, iY, iZ, 1) ) / (T)2;
                el[tensor::zz] = fdNTensorField::planeZderiv(velocity, direction,orientation, iX, iY, iZ, 2);
            }
        }
    }
}

template<typename T>
void BoxNTensorStrainRateFunctional3D<T>::processEdge (
        int plane, int normal1, int normal2, Box3D domain,
        NTensorField3D<T>& velocity, NTensorField3D<T>& S )
{
    PLB_PRECONDITION( velocity.getNdim()==3 );
    PLB_PRECONDITION( S.getNdim()==6 );
    typedef SymmetricTensorImpl<T,3> tensor;
    Dot3D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                plint iZ2 = iZ+offset.z;
                T* el = S.get(iX2,iY2,iZ2);
                el[tensor::xx] = fdNTensorField::edgeXderiv(velocity, plane,normal1,normal2, iX, iY, iZ, 0);
                el[tensor::xy] = ( fdNTensorField::edgeXderiv(velocity, plane,normal1,normal2, iX, iY, iZ, 1) +
                               fdNTensorField::edgeYderiv(velocity, plane,normal1,normal2, iX, iY, iZ, 0) ) / (T)2;
                el[tensor::xz] = ( fdNTensorField::edgeXderiv(velocity, plane,normal1,normal2, iX, iY, iZ, 2) +
                               fdNTensorField::edgeZderiv(velocity, plane,normal1,normal2, iX, iY, iZ, 0) ) / (T)2;
                el[tensor::yy] = fdNTensorField::edgeYderiv(velocity, plane,normal1,normal2, iX, iY, iZ, 1);
                el[tensor::yz] = ( fdNTensorField::edgeYderiv(velocity, plane,normal1,normal2, iX, iY, iZ, 2) +
                               fdNTensorField::edgeZderiv(velocity, plane,normal1,normal2, iX, iY, iZ, 1) ) / (T)2;
                el[tensor::zz] = fdNTensorField::edgeZderiv(velocity, plane,normal1,normal2, iX, iY, iZ, 2);
            }
        }
    }
}

template<typename T>
void BoxNTensorStrainRateFunctional3D<T>::processCorner (
        int normalX, int normalY, int normalZ, Box3D domain,
        NTensorField3D<T>& velocity, NTensorField3D<T>& S )
{
    PLB_PRECONDITION( velocity.getNdim()==3 );
    PLB_PRECONDITION( S.getNdim()==6 );
    typedef SymmetricTensorImpl<T,3> tensor;
    Dot3D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                plint iX2 = iX+offset.x;
                plint iY2 = iY+offset.y;
                plint iZ2 = iZ+offset.z;
                T* el = S.get(iX2,iY2,iZ2);
                el[tensor::xx] = fdNTensorField::cornerXderiv(velocity, normalX,normalY,normalZ, iX, iY, iZ, 0);
                el[tensor::xy] = ( fdNTensorField::cornerXderiv(velocity, normalX,normalY,normalZ, iX, iY, iZ, 1) +
                                   fdNTensorField::cornerYderiv(velocity, normalX,normalY,normalZ, iX, iY, iZ, 0) ) / (T)2;
                el[tensor::xz] = ( fdNTensorField::cornerXderiv(velocity, normalX,normalY,normalZ, iX, iY, iZ, 2) +
                                   fdNTensorField::cornerZderiv(velocity, normalX,normalY,normalZ, iX, iY, iZ, 0) ) / (T)2;
                el[tensor::yy] = fdNTensorField::cornerYderiv(velocity, normalX,normalY,normalZ, iX, iY, iZ, 1);
                el[tensor::yz] = ( fdNTensorField::cornerYderiv(velocity, normalX,normalY,normalZ, iX, iY, iZ, 2) +
                                   fdNTensorField::cornerZderiv(velocity, normalX,normalY,normalZ, iX, iY, iZ, 1) ) / (T)2;
                el[tensor::zz] = fdNTensorField::cornerZderiv(velocity, normalX,normalY,normalZ, iX, iY, iZ, 2);
            }
        }
    }
}


template<typename T>
BoxNTensorStrainRateFunctional3D<T>* BoxNTensorStrainRateFunctional3D<T>::clone() const
{
    return new BoxNTensorStrainRateFunctional3D<T>(*this);
}

template<typename T>
void BoxNTensorStrainRateFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}


template<typename T>
BlockDomain::DomainT BoxNTensorStrainRateFunctional3D<T>::appliesTo() const {
    // Don't apply to envelope, because nearest neighbors need to be accessed.
    return BlockDomain::bulk;
}


template<typename T>
void MaskedBoxNTensorStrainRateFunctional3D<T>::processBulk (
        Box3D domain, NTensorField3D<T>& velocity, NTensorField3D<T>& S,
        NTensorField3D<int>& mask)
{
    PLB_PRECONDITION( velocity.getNdim()==3 );
    PLB_PRECONDITION( S.getNdim()==6 );
    typedef SymmetricTensorImpl<T,3> tensor;
    Dot3D maskOfs = computeRelativeDisplacement(velocity, mask);
    Dot3D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                    plint iX2 = iX+offset.x;
                    plint iY2 = iY+offset.y;
                    plint iZ2 = iZ+offset.z;
                    T* el = S.get(iX2,iY2,iZ2);
                    el[tensor::xx] = fdNTensorField::bulkXderiv(velocity, iX, iY, iZ, 0);
                    el[tensor::xy] = ( fdNTensorField::bulkXderiv(velocity, iX, iY, iZ, 1) +
                                       fdNTensorField::bulkYderiv(velocity, iX, iY, iZ, 0) ) / (T)2;
                    el[tensor::xz] = ( fdNTensorField::bulkXderiv(velocity, iX, iY, iZ, 2) +
                                       fdNTensorField::bulkZderiv(velocity, iX, iY, iZ, 0) ) / (T)2;
                    el[tensor::yy] = fdNTensorField::bulkYderiv(velocity, iX, iY, iZ, 1);
                    el[tensor::yz] = ( fdNTensorField::bulkYderiv(velocity, iX, iY, iZ, 2) +
                                       fdNTensorField::bulkZderiv(velocity, iX, iY, iZ, 1) ) / (T)2;
                    el[tensor::zz] = fdNTensorField::bulkZderiv(velocity, iX, iY, iZ, 2);
                }
            }
        }
    }
}

template<typename T>
void MaskedBoxNTensorStrainRateFunctional3D<T>::processPlane (
        int direction, int orientation, Box3D domain,
        NTensorField3D<T>& velocity, NTensorField3D<T>& S,
        NTensorField3D<int>& mask)
{
    PLB_PRECONDITION( velocity.getNdim()==3 );
    PLB_PRECONDITION( S.getNdim()==6 );
    Dot3D maskOfs = computeRelativeDisplacement(velocity, mask);
    typedef SymmetricTensorImpl<T,3> tensor;
    Dot3D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                    plint iX2 = iX+offset.x;
                    plint iY2 = iY+offset.y;
                    plint iZ2 = iZ+offset.z;
                    T* el = S.get(iX2,iY2,iZ2);
                    el[tensor::xx] = fdNTensorField::planeXderiv(velocity, direction,orientation, iX, iY, iZ, 0);
                    el[tensor::xy] = ( fdNTensorField::planeXderiv(velocity, direction,orientation, iX, iY, iZ, 1) +
                                   fdNTensorField::planeYderiv(velocity, direction,orientation, iX, iY, iZ, 0) ) / (T)2;
                    el[tensor::xz] = ( fdNTensorField::planeXderiv(velocity, direction,orientation, iX, iY, iZ, 2) +
                                   fdNTensorField::planeZderiv(velocity, direction,orientation, iX, iY, iZ, 0) ) / (T)2;
                    el[tensor::yy] = fdNTensorField::planeYderiv(velocity, direction,orientation, iX, iY, iZ, 1);
                    el[tensor::yz] = ( fdNTensorField::planeYderiv(velocity, direction,orientation, iX, iY, iZ, 2) +
                                   fdNTensorField::planeZderiv(velocity, direction,orientation, iX, iY, iZ, 1) ) / (T)2;
                    el[tensor::zz] = fdNTensorField::planeZderiv(velocity, direction,orientation, iX, iY, iZ, 2);
                }
            }
        }
    }
}

template<typename T>
void MaskedBoxNTensorStrainRateFunctional3D<T>::processEdge (
        int plane, int normal1, int normal2, Box3D domain,
        NTensorField3D<T>& velocity, NTensorField3D<T>& S,
        NTensorField3D<int>& mask)
{
    PLB_PRECONDITION( velocity.getNdim()==3 );
    PLB_PRECONDITION( S.getNdim()==6 );
    Dot3D maskOfs = computeRelativeDisplacement(velocity, mask);
    typedef SymmetricTensorImpl<T,3> tensor;
    Dot3D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                    plint iX2 = iX+offset.x;
                    plint iY2 = iY+offset.y;
                    plint iZ2 = iZ+offset.z;
                    T* el = S.get(iX2,iY2,iZ2);
                    el[tensor::xx] = fdNTensorField::edgeXderiv(velocity, plane,normal1,normal2, iX, iY, iZ, 0);
                    el[tensor::xy] = ( fdNTensorField::edgeXderiv(velocity, plane,normal1,normal2, iX, iY, iZ, 1) +
                                   fdNTensorField::edgeYderiv(velocity, plane,normal1,normal2, iX, iY, iZ, 0) ) / (T)2;
                    el[tensor::xz] = ( fdNTensorField::edgeXderiv(velocity, plane,normal1,normal2, iX, iY, iZ, 2) +
                                   fdNTensorField::edgeZderiv(velocity, plane,normal1,normal2, iX, iY, iZ, 0) ) / (T)2;
                    el[tensor::yy] = fdNTensorField::edgeYderiv(velocity, plane,normal1,normal2, iX, iY, iZ, 1);
                    el[tensor::yz] = ( fdNTensorField::edgeYderiv(velocity, plane,normal1,normal2, iX, iY, iZ, 2) +
                                   fdNTensorField::edgeZderiv(velocity, plane,normal1,normal2, iX, iY, iZ, 1) ) / (T)2;
                    el[tensor::zz] = fdNTensorField::edgeZderiv(velocity, plane,normal1,normal2, iX, iY, iZ, 2);
                }
            }
        }
    }
}

template<typename T>
void MaskedBoxNTensorStrainRateFunctional3D<T>::processCorner (
        int normalX, int normalY, int normalZ, Box3D domain,
        NTensorField3D<T>& velocity, NTensorField3D<T>& S,
        NTensorField3D<int>& mask)
{
    PLB_PRECONDITION( velocity.getNdim()==3 );
    PLB_PRECONDITION( S.getNdim()==6 );
    Dot3D maskOfs = computeRelativeDisplacement(velocity, mask);
    typedef SymmetricTensorImpl<T,3> tensor;
    Dot3D offset = computeRelativeDisplacement(velocity, S);
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                    plint iX2 = iX+offset.x;
                    plint iY2 = iY+offset.y;
                    plint iZ2 = iZ+offset.z;
                    T* el = S.get(iX2,iY2,iZ2);
                    el[tensor::xx] = fdNTensorField::cornerXderiv(velocity, normalX,normalY,normalZ, iX, iY, iZ, 0);
                    el[tensor::xy] = ( fdNTensorField::cornerXderiv(velocity, normalX,normalY,normalZ, iX, iY, iZ, 1) +
                                       fdNTensorField::cornerYderiv(velocity, normalX,normalY,normalZ, iX, iY, iZ, 0) ) / (T)2;
                    el[tensor::xz] = ( fdNTensorField::cornerXderiv(velocity, normalX,normalY,normalZ, iX, iY, iZ, 2) +
                                       fdNTensorField::cornerZderiv(velocity, normalX,normalY,normalZ, iX, iY, iZ, 0) ) / (T)2;
                    el[tensor::yy] = fdNTensorField::cornerYderiv(velocity, normalX,normalY,normalZ, iX, iY, iZ, 1);
                    el[tensor::yz] = ( fdNTensorField::cornerYderiv(velocity, normalX,normalY,normalZ, iX, iY, iZ, 2) +
                                       fdNTensorField::cornerZderiv(velocity, normalX,normalY,normalZ, iX, iY, iZ, 1) ) / (T)2;
                    el[tensor::zz] = fdNTensorField::cornerZderiv(velocity, normalX,normalY,normalZ, iX, iY, iZ, 2);
                }
            }
        }
    }
}


template<typename T>
MaskedBoxNTensorStrainRateFunctional3D<T>* MaskedBoxNTensorStrainRateFunctional3D<T>::clone() const
{
    return new MaskedBoxNTensorStrainRateFunctional3D<T>(*this);
}

template<typename T>
void MaskedBoxNTensorStrainRateFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}


template<typename T>
BlockDomain::DomainT MaskedBoxNTensorStrainRateFunctional3D<T>::appliesTo() const {
    // Don't apply to envelope, because nearest neighbors need to be accessed.
    return BlockDomain::bulk;
}


/* ******** A_plus_alpha_NTensor3D ************************************* */

template<typename T>
A_plus_alpha_NTensor3D<T>::A_plus_alpha_NTensor3D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_plus_alpha_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<T>& result )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot3D offset = computeRelativeDisplacement(A, result);

    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    *result.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                        = *A.get(iX,iY,iX) + scalarAlpha;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim]
                            = A.get(iX,iY,iZ)[iDim] + alpha[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
A_plus_alpha_NTensor3D<T>* A_plus_alpha_NTensor3D<T>::clone() const {
    return new A_plus_alpha_NTensor3D<T>(*this);
}

template<typename T>
void A_plus_alpha_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_plus_alpha_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_plus_alpha_NTensor3D ************************************* */

template<typename T>
Masked_A_plus_alpha_NTensor3D<T>::Masked_A_plus_alpha_NTensor3D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void Masked_A_plus_alpha_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<T>& result,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot3D maskOfs = computeRelativeDisplacement(A, mask);
    Dot3D offset = computeRelativeDisplacement(A, result);

    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        *result.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                            = *A.get(iX,iY,iZ) + scalarAlpha;
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            result.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim]
                                = A.get(iX,iY,iZ)[iDim] + alpha[iDim];
                        }
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_plus_alpha_NTensor3D<T>* Masked_A_plus_alpha_NTensor3D<T>::clone() const {
    return new Masked_A_plus_alpha_NTensor3D<T>(*this);
}

template<typename T>
void Masked_A_plus_alpha_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_plus_alpha_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** A_minus_alpha_NTensor3D ************************************** */

template<typename T>
A_minus_alpha_NTensor3D<T>::A_minus_alpha_NTensor3D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_minus_alpha_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<T>& result )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot3D offset = computeRelativeDisplacement(A, result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    *result.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                        = *A.get(iX,iY,iZ) - scalarAlpha;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim]
                            = A.get(iX,iY,iZ)[iDim] - alpha[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
A_minus_alpha_NTensor3D<T>* A_minus_alpha_NTensor3D<T>::clone() const {
    return new A_minus_alpha_NTensor3D<T>(*this);
}

template<typename T>
void A_minus_alpha_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_minus_alpha_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}



/* ******** Masked_A_minus_alpha_NTensor3D ************************************** */

template<typename T>
Masked_A_minus_alpha_NTensor3D<T>::Masked_A_minus_alpha_NTensor3D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void Masked_A_minus_alpha_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<T>& result,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot3D maskOfs = computeRelativeDisplacement(A, mask);
    Dot3D offset = computeRelativeDisplacement(A, result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        *result.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                            = *A.get(iX,iY,iZ) - scalarAlpha;
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            result.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim]
                                = A.get(iX,iY,iZ)[iDim] - alpha[iDim];
                        }
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_minus_alpha_NTensor3D<T>* Masked_A_minus_alpha_NTensor3D<T>::clone() const {
    return new Masked_A_minus_alpha_NTensor3D<T>(*this);
}

template<typename T>
void Masked_A_minus_alpha_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_minus_alpha_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}



/* ******** Alpha_minus_A_NTensor3D ************************************* */

template<typename T>
Alpha_minus_A_NTensor3D<T>::Alpha_minus_A_NTensor3D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void Alpha_minus_A_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<T>& result )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot3D offset = computeRelativeDisplacement(A, result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    *result.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                        = scalarAlpha - *A.get(iX,iY,iZ);
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim]
                            = alpha[iDim] - A.get(iX,iY,iZ)[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
Alpha_minus_A_NTensor3D<T>* Alpha_minus_A_NTensor3D<T>::clone() const {
    return new Alpha_minus_A_NTensor3D<T>(*this);
}

template<typename T>
void Alpha_minus_A_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT Alpha_minus_A_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_Alpha_minus_A_NTensor3D ************************************* */

template<typename T>
Masked_Alpha_minus_A_NTensor3D<T>::Masked_Alpha_minus_A_NTensor3D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void Masked_Alpha_minus_A_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<T>& result,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot3D maskOfs = computeRelativeDisplacement(A, mask);
    Dot3D offset = computeRelativeDisplacement(A, result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        *result.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                            = scalarAlpha - *A.get(iX,iY,iZ);
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            result.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim]
                                = alpha[iDim] - A.get(iX,iY,iZ)[iDim];
                        }
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_Alpha_minus_A_NTensor3D<T>* Masked_Alpha_minus_A_NTensor3D<T>::clone() const {
    return new Masked_Alpha_minus_A_NTensor3D<T>(*this);
}

template<typename T>
void Masked_Alpha_minus_A_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_Alpha_minus_A_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** A_times_alpha_NTensor3D ************************************* */

template<typename T>
A_times_alpha_NTensor3D<T>::A_times_alpha_NTensor3D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_times_alpha_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<T>& result )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot3D offset = computeRelativeDisplacement(A, result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    *result.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                        = *A.get(iX,iY,iZ) * scalarAlpha;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim]
                            = A.get(iX,iY,iZ)[iDim] * alpha[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
A_times_alpha_NTensor3D<T>* A_times_alpha_NTensor3D<T>::clone() const {
    return new A_times_alpha_NTensor3D<T>(*this);
}

template<typename T>
void A_times_alpha_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_times_alpha_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_times_alpha_NTensor3D ************************************* */

template<typename T>
Masked_A_times_alpha_NTensor3D<T>::Masked_A_times_alpha_NTensor3D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void Masked_A_times_alpha_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<T>& result,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot3D maskOfs = computeRelativeDisplacement(A, mask);
    Dot3D offset = computeRelativeDisplacement(A, result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        *result.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                            = *A.get(iX,iY,iZ) * scalarAlpha;
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            result.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim]
                                = A.get(iX,iY,iZ)[iDim] * alpha[iDim];
                        }
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_times_alpha_NTensor3D<T>* Masked_A_times_alpha_NTensor3D<T>::clone() const {
    return new Masked_A_times_alpha_NTensor3D<T>(*this);
}

template<typename T>
void Masked_A_times_alpha_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_times_alpha_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** A_dividedBy_alpha_NTensor3D ************************************* */

template<typename T>
A_dividedBy_alpha_NTensor3D<T>::A_dividedBy_alpha_NTensor3D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_dividedBy_alpha_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<T>& result )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot3D offset = computeRelativeDisplacement(A, result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    *result.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                        = *A.get(iX,iY,iZ) / scalarAlpha;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim]
                            = A.get(iX,iY,iZ)[iDim] / alpha[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
A_dividedBy_alpha_NTensor3D<T>* A_dividedBy_alpha_NTensor3D<T>::clone() const {
    return new A_dividedBy_alpha_NTensor3D<T>(*this);
}

template<typename T>
void A_dividedBy_alpha_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_dividedBy_alpha_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_dividedBy_alpha_NTensor3D ************************************* */

template<typename T>
Masked_A_dividedBy_alpha_NTensor3D<T>::Masked_A_dividedBy_alpha_NTensor3D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void Masked_A_dividedBy_alpha_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<T>& result,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot3D maskOfs = computeRelativeDisplacement(A, mask);
    Dot3D offset = computeRelativeDisplacement(A, result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        *result.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                            = *A.get(iX,iY,iZ) / scalarAlpha;
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            result.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim]
                                = A.get(iX,iY,iZ)[iDim] / alpha[iDim];
                        }
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_dividedBy_alpha_NTensor3D<T>* Masked_A_dividedBy_alpha_NTensor3D<T>::clone() const {
    return new Masked_A_dividedBy_alpha_NTensor3D<T>(*this);
}

template<typename T>
void Masked_A_dividedBy_alpha_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_dividedBy_alpha_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Alpha_dividedBy_A_NTensor3D ************************************* */

template<typename T>
Alpha_dividedBy_A_NTensor3D<T>::Alpha_dividedBy_A_NTensor3D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void Alpha_dividedBy_A_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<T>& result )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot3D offset = computeRelativeDisplacement(A, result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    *result.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                        = scalarAlpha / *A.get(iX,iY,iZ);
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim]
                            = alpha[iDim] / A.get(iX,iY,iZ)[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
Alpha_dividedBy_A_NTensor3D<T>* Alpha_dividedBy_A_NTensor3D<T>::clone() const {
    return new Alpha_dividedBy_A_NTensor3D<T>(*this);
}

template<typename T>
void Alpha_dividedBy_A_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT Alpha_dividedBy_A_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}



/* ******** Masked_Alpha_dividedBy_A_NTensor3D ************************************* */

template<typename T>
Masked_Alpha_dividedBy_A_NTensor3D<T>::Masked_Alpha_dividedBy_A_NTensor3D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void Masked_Alpha_dividedBy_A_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<T>& result,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot3D maskOfs = computeRelativeDisplacement(A, mask);
    Dot3D offset = computeRelativeDisplacement(A, result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        *result.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                            = scalarAlpha / *A.get(iX,iY,iZ);
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                   if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            result.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim]
                                = alpha[iDim] / A.get(iX,iY,iZ)[iDim];
                        }
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_Alpha_dividedBy_A_NTensor3D<T>* Masked_Alpha_dividedBy_A_NTensor3D<T>::clone() const {
    return new Masked_Alpha_dividedBy_A_NTensor3D<T>(*this);
}

template<typename T>
void Masked_Alpha_dividedBy_A_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_Alpha_dividedBy_A_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}



/* ******** A_toThePower_alpha_NTensor3D ************************************* */

template<typename T>
A_toThePower_alpha_NTensor3D<T>::A_toThePower_alpha_NTensor3D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_toThePower_alpha_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<T>& result )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot3D offset = computeRelativeDisplacement(A, result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    *result.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                        = customPower(*A.get(iX,iY,iZ),scalarAlpha);
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim]
                            = customPower(A.get(iX,iY,iZ)[iDim],alpha[iDim]);
                    }
                }
            }
        }
    }
}

template<typename T>
A_toThePower_alpha_NTensor3D<T>* A_toThePower_alpha_NTensor3D<T>::clone() const {
    return new A_toThePower_alpha_NTensor3D<T>(*this);
}

template<typename T>
void A_toThePower_alpha_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_toThePower_alpha_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_toThePower_alpha_NTensor3D ************************************* */

template<typename T>
Masked_A_toThePower_alpha_NTensor3D<T>::Masked_A_toThePower_alpha_NTensor3D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void Masked_A_toThePower_alpha_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<T>& result,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot3D maskOfs = computeRelativeDisplacement(A, mask);
    Dot3D offset = computeRelativeDisplacement(A, result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        *result.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                            = customPower(*A.get(iX,iY,iZ),scalarAlpha);
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            result.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim]
                                = customPower(A.get(iX,iY,iZ)[iDim],alpha[iDim]);
                        }
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_toThePower_alpha_NTensor3D<T>* Masked_A_toThePower_alpha_NTensor3D<T>::clone() const {
    return new Masked_A_toThePower_alpha_NTensor3D<T>(*this);
}

template<typename T>
void Masked_A_toThePower_alpha_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_toThePower_alpha_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Alpha_toThePower_A_NTensor3D ************************************* */

template<typename T>
Alpha_toThePower_A_NTensor3D<T>::Alpha_toThePower_A_NTensor3D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void Alpha_toThePower_A_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<T>& result )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot3D offset = computeRelativeDisplacement(A, result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    *result.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                        = customPower(scalarAlpha,*A.get(iX,iY,iZ));
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim]
                            = customPower(alpha[iDim],A.get(iX,iY,iZ)[iDim]);
                    }
                }
            }
        }
    }
}

template<typename T>
Alpha_toThePower_A_NTensor3D<T>* Alpha_toThePower_A_NTensor3D<T>::clone() const {
    return new Alpha_toThePower_A_NTensor3D<T>(*this);
}

template<typename T>
void Alpha_toThePower_A_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT Alpha_toThePower_A_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}



/* ******** Masked_Alpha_toThePower_A_NTensor3D ************************************* */

template<typename T>
Masked_Alpha_toThePower_A_NTensor3D<T>::Masked_Alpha_toThePower_A_NTensor3D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void Masked_Alpha_toThePower_A_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<T>& result,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot3D maskOfs = computeRelativeDisplacement(A, mask);
    Dot3D offset = computeRelativeDisplacement(A, result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        *result.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                            = customPower(scalarAlpha,*A.get(iX,iY,iZ));
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            result.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim]
                                = customPower(alpha[iDim],A.get(iX,iY,iZ)[iDim]);
                        }
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_Alpha_toThePower_A_NTensor3D<T>* Masked_Alpha_toThePower_A_NTensor3D<T>::clone() const {
    return new Masked_Alpha_toThePower_A_NTensor3D<T>(*this);
}

template<typename T>
void Masked_Alpha_toThePower_A_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_Alpha_toThePower_A_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}



/* ******** A_plus_B_NTensor3D ************************************ */

template<typename T>
void A_plus_B_NTensor3D<T>::process (
        Box3D domain, std::vector<NTensorField3D<T>*> fields )
{
    PLB_PRECONDITION( fields.size()==3 );
    NTensorField3D<T>& A = *fields[0];
    NTensorField3D<T>& B = *fields[1];
    NTensorField3D<T>& result = *fields[2];
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    *result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                        = *A.get(iX,iY,iZ) + *B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z);
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)[iDim]
                            = A.get(iX,iY,iZ)[iDim] + B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z)[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
A_plus_B_NTensor3D<T>* A_plus_B_NTensor3D<T>::clone() const {
    return new A_plus_B_NTensor3D<T>(*this);
}

template<typename T>
void A_plus_B_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_plus_B_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_plus_B_NTensor3D ************************************ */

template<typename T>
void Masked_A_plus_B_NTensor3D<T>::process (
        Box3D domain, std::vector<NTensorField3D<T>*> fields,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( fields.size()==3 );
    NTensorField3D<T>& A = *fields[0];
    NTensorField3D<T>& B = *fields[1];
    NTensorField3D<T>& result = *fields[2];
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    Dot3D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        *result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                            = *A.get(iX,iY,iZ) + *B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z);
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)[iDim]
                                = A.get(iX,iY,iZ)[iDim] + B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z)[iDim];
                        }
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_plus_B_NTensor3D<T>* Masked_A_plus_B_NTensor3D<T>::clone() const {
    return new Masked_A_plus_B_NTensor3D<T>(*this);
}

template<typename T>
void Masked_A_plus_B_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
    modified[3] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_plus_B_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** A_minus_B_NTensor3D ************************************ */

template<typename T>
void A_minus_B_NTensor3D<T>::process (
        Box3D domain, std::vector<NTensorField3D<T>*> fields )
{
    PLB_PRECONDITION( fields.size()==3 );
    NTensorField3D<T>& A = *fields[0];
    NTensorField3D<T>& B = *fields[1];
    NTensorField3D<T>& result = *fields[2];
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    *result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                        = *A.get(iX,iY,iZ) - *B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z);
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)[iDim]
                            = A.get(iX,iY,iZ)[iDim] - B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z)[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
A_minus_B_NTensor3D<T>* A_minus_B_NTensor3D<T>::clone() const {
    return new A_minus_B_NTensor3D<T>(*this);
}

template<typename T>
void A_minus_B_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_minus_B_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_minus_B_NTensor3D ************************************ */

template<typename T>
void Masked_A_minus_B_NTensor3D<T>::process (
        Box3D domain, std::vector<NTensorField3D<T>*> fields,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( fields.size()==3 );
    NTensorField3D<T>& A = *fields[0];
    NTensorField3D<T>& B = *fields[1];
    NTensorField3D<T>& result = *fields[2];
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    Dot3D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        *result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                            = *A.get(iX,iY,iZ) - *B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z);
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)[iDim]
                                = A.get(iX,iY,iZ)[iDim] - B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z)[iDim];
                        }
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_minus_B_NTensor3D<T>* Masked_A_minus_B_NTensor3D<T>::clone() const {
    return new Masked_A_minus_B_NTensor3D<T>(*this);
}

template<typename T>
void Masked_A_minus_B_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
    modified[3] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_minus_B_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** A_times_B_NTensor3D ************************************ */

template<typename T>
void A_times_B_NTensor3D<T>::process (
        Box3D domain, std::vector<NTensorField3D<T>*> fields )
{
    PLB_PRECONDITION( fields.size()==3 );
    NTensorField3D<T>& A = *fields[0];
    NTensorField3D<T>& B = *fields[1];
    NTensorField3D<T>& result = *fields[2];
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    *result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                        = *A.get(iX,iY,iZ) * *B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z);
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)[iDim]
                            = A.get(iX,iY,iZ)[iDim] * B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z)[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
A_times_B_NTensor3D<T>* A_times_B_NTensor3D<T>::clone() const {
    return new A_times_B_NTensor3D<T>(*this);
}

template<typename T>
void A_times_B_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_times_B_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_times_B_NTensor3D ************************************ */

template<typename T>
void Masked_A_times_B_NTensor3D<T>::process (
        Box3D domain, std::vector<NTensorField3D<T>*> fields,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( fields.size()==3 );
    NTensorField3D<T>& A = *fields[0];
    NTensorField3D<T>& B = *fields[1];
    NTensorField3D<T>& result = *fields[2];
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    Dot3D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        *result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                            = *A.get(iX,iY,iZ) * *B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z);
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)[iDim]
                                = A.get(iX,iY,iZ)[iDim] * B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z)[iDim];
                        }
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_times_B_NTensor3D<T>* Masked_A_times_B_NTensor3D<T>::clone() const {
    return new Masked_A_times_B_NTensor3D<T>(*this);
}

template<typename T>
void Masked_A_times_B_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
    modified[3] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_times_B_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** A_dividedBy_B_NTensor3D ************************************ */

template<typename T>
void A_dividedBy_B_NTensor3D<T>::process (
        Box3D domain, std::vector<NTensorField3D<T>*> fields )
{
    PLB_PRECONDITION( fields.size()==3 );
    NTensorField3D<T>& A = *fields[0];
    NTensorField3D<T>& B = *fields[1];
    NTensorField3D<T>& result = *fields[2];
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    *result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                        = *A.get(iX,iY,iZ) / *B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z);
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)[iDim]
                            = A.get(iX,iY,iZ)[iDim] / B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z)[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
A_dividedBy_B_NTensor3D<T>* A_dividedBy_B_NTensor3D<T>::clone() const {
    return new A_dividedBy_B_NTensor3D<T>(*this);
}

template<typename T>
void A_dividedBy_B_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_dividedBy_B_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_dividedBy_B_NTensor3D ************************************ */

template<typename T>
void Masked_A_dividedBy_B_NTensor3D<T>::process (
        Box3D domain, std::vector<NTensorField3D<T>*> fields,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( fields.size()==3 );
    NTensorField3D<T>& A = *fields[0];
    NTensorField3D<T>& B = *fields[1];
    NTensorField3D<T>& result = *fields[2];
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    Dot3D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        *result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                            = *A.get(iX,iY,iZ) / *B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z);
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)[iDim]
                                = A.get(iX,iY,iZ)[iDim] / B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z)[iDim];
                        }
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_dividedBy_B_NTensor3D<T>* Masked_A_dividedBy_B_NTensor3D<T>::clone() const {
    return new Masked_A_dividedBy_B_NTensor3D<T>(*this);
}

template<typename T>
void Masked_A_dividedBy_B_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
    modified[3] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_dividedBy_B_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** A_toThePower_B_NTensor3D ************************************ */

template<typename T>
void A_toThePower_B_NTensor3D<T>::process (
        Box3D domain, std::vector<NTensorField3D<T>*> fields )
{
    PLB_PRECONDITION( fields.size()==3 );
    NTensorField3D<T>& A = *fields[0];
    NTensorField3D<T>& B = *fields[1];
    NTensorField3D<T>& result = *fields[2];
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    *result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                        = customPower(*A.get(iX,iY,iZ),*B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z));
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)[iDim]
                            = customPower(A.get(iX,iY,iZ)[iDim],B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z)[iDim]);
                    }
                }
            }
        }
    }
}

template<typename T>
A_toThePower_B_NTensor3D<T>* A_toThePower_B_NTensor3D<T>::clone() const {
    return new A_toThePower_B_NTensor3D<T>(*this);
}

template<typename T>
void A_toThePower_B_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_toThePower_B_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_toThePower_B_NTensor3D ************************************ */

template<typename T>
void Masked_A_toThePower_B_NTensor3D<T>::process (
        Box3D domain, std::vector<NTensorField3D<T>*> fields,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( fields.size()==3 );
    NTensorField3D<T>& A = *fields[0];
    NTensorField3D<T>& B = *fields[1];
    NTensorField3D<T>& result = *fields[2];
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    Dot3D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        *result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                            = customPower(*A.get(iX,iY,iZ),*B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z));
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)[iDim]
                                = customPower(A.get(iX,iY,iZ)[iDim],B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z)[iDim]);
                        }
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_toThePower_B_NTensor3D<T>* Masked_A_toThePower_B_NTensor3D<T>::clone() const {
    return new Masked_A_toThePower_B_NTensor3D<T>(*this);
}

template<typename T>
void Masked_A_toThePower_B_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
    modified[3] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_toThePower_B_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** A_equals_B_NTensor3D ************************************ */

template<typename T>
void A_equals_B_NTensor3D<T>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==3 );
    NTensorField3D<T>& A = *dynamic_cast<NTensorField3D<T>*>(blocks[0]);
    NTensorField3D<T>& B = *dynamic_cast<NTensorField3D<T>*>(blocks[1]);
    NTensorField3D<int>& result = *dynamic_cast<NTensorField3D<int>*>(blocks[2]);
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    *result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                        = *A.get(iX,iY,iZ) == *B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z) ? 1:0;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)[iDim]
                            = A.get(iX,iY,iZ)[iDim] == B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z)[iDim] ? 1:0;
                    }
                }
            }
        }
    }
}

template<typename T>
A_equals_B_NTensor3D<T>* A_equals_B_NTensor3D<T>::clone() const {
    return new A_equals_B_NTensor3D<T>(*this);
}

template<typename T>
void A_equals_B_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_equals_B_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_equals_B_NTensor3D ************************************ */

template<typename T>
void Masked_A_equals_B_NTensor3D<T>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==3 );
    NTensorField3D<T>& A = *dynamic_cast<NTensorField3D<T>*>(blocks[0]);
    NTensorField3D<T>& B = *dynamic_cast<NTensorField3D<T>*>(blocks[1]);
    NTensorField3D<int>& result = *dynamic_cast<NTensorField3D<int>*>(blocks[2]);
    NTensorField3D<int>& mask = *dynamic_cast<NTensorField3D<int>*>(blocks[3]);
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    Dot3D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        *result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                            = *A.get(iX,iY,iZ) == *B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z) ? 1:0;
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)[iDim]
                                = A.get(iX,iY,iZ)[iDim] == B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z)[iDim] ? 1:0;
                        }
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_equals_B_NTensor3D<T>* Masked_A_equals_B_NTensor3D<T>::clone() const {
    return new Masked_A_equals_B_NTensor3D<T>(*this);
}

template<typename T>
void Masked_A_equals_B_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
    modified[3] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_equals_B_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** A_equals_alpha_NTensor3D ************************************* */

template<typename T>
A_equals_alpha_NTensor3D<T>::A_equals_alpha_NTensor3D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_equals_alpha_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<int>& result )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot3D offset = computeRelativeDisplacement(A, result);

    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    *result.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                        = *A.get(iX,iY,iZ) == scalarAlpha ? 1 : 0;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim]
                            = A.get(iX,iY,iZ)[iDim] == alpha[iDim] ? 1 : 0;
                    }
                }
            }
        }
    }
}

template<typename T>
A_equals_alpha_NTensor3D<T>* A_equals_alpha_NTensor3D<T>::clone() const {
    return new A_equals_alpha_NTensor3D<T>(*this);
}

template<typename T>
void A_equals_alpha_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_equals_alpha_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_equals_alpha_NTensor3D ************************************* */

template<typename T>
Masked_A_equals_alpha_NTensor3D<T>::Masked_A_equals_alpha_NTensor3D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void Masked_A_equals_alpha_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<int>& result,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot3D offset = computeRelativeDisplacement(A, result);
    Dot3D maskOfs = computeRelativeDisplacement(A, mask);

    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        *result.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                            = *A.get(iX,iY,iZ) == scalarAlpha ? 1 : 0;
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            result.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim]
                                = A.get(iX,iY,iZ)[iDim] == alpha[iDim] ? 1 : 0;
                        }
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_equals_alpha_NTensor3D<T>* Masked_A_equals_alpha_NTensor3D<T>::clone() const {
    return new Masked_A_equals_alpha_NTensor3D<T>(*this);
}

template<typename T>
void Masked_A_equals_alpha_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_equals_alpha_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** A_lessThan_B_NTensor3D ************************************ */

template<typename T>
void A_lessThan_B_NTensor3D<T>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==3 );
    NTensorField3D<T>& A = *dynamic_cast<NTensorField3D<T>*>(blocks[0]);
    NTensorField3D<T>& B = *dynamic_cast<NTensorField3D<T>*>(blocks[1]);
    NTensorField3D<int>& result = *dynamic_cast<NTensorField3D<int>*>(blocks[2]);
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    *result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                        = *A.get(iX,iY,iZ) < *B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z) ? 1:0;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)[iDim]
                            = A.get(iX,iY,iZ)[iDim] < B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z)[iDim] ? 1:0;
                    }
                }
            }
        }
    }
}

template<typename T>
A_lessThan_B_NTensor3D<T>* A_lessThan_B_NTensor3D<T>::clone() const {
    return new A_lessThan_B_NTensor3D<T>(*this);
}

template<typename T>
void A_lessThan_B_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_lessThan_B_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}

/* ******** Masked_A_lessThan_B_NTensor3D ************************************ */

template<typename T>
void Masked_A_lessThan_B_NTensor3D<T>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==3 );
    NTensorField3D<T>& A = *dynamic_cast<NTensorField3D<T>*>(blocks[0]);
    NTensorField3D<T>& B = *dynamic_cast<NTensorField3D<T>*>(blocks[1]);
    NTensorField3D<int>& result = *dynamic_cast<NTensorField3D<int>*>(blocks[2]);
    NTensorField3D<int>& mask = *dynamic_cast<NTensorField3D<int>*>(blocks[3]);
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    Dot3D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        *result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                            = *A.get(iX,iY,iZ) < *B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z) ? 1:0;
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)[iDim]
                                = A.get(iX,iY,iZ)[iDim] < B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z)[iDim] ? 1:0;
                        }
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_lessThan_B_NTensor3D<T>* Masked_A_lessThan_B_NTensor3D<T>::clone() const {
    return new Masked_A_lessThan_B_NTensor3D<T>(*this);
}

template<typename T>
void Masked_A_lessThan_B_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
    modified[3] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_lessThan_B_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** A_lessThan_alpha_NTensor3D ************************************* */

template<typename T>
A_lessThan_alpha_NTensor3D<T>::A_lessThan_alpha_NTensor3D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_lessThan_alpha_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<int>& result )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot3D offset = computeRelativeDisplacement(A, result);

    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    *result.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                        = *A.get(iX,iY,iZ) < scalarAlpha ? 1 : 0;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim]
                            = A.get(iX,iY,iZ)[iDim] < alpha[iDim] ? 1 : 0;
                    }
                }
            }
        }
    }
}

template<typename T>
A_lessThan_alpha_NTensor3D<T>* A_lessThan_alpha_NTensor3D<T>::clone() const {
    return new A_lessThan_alpha_NTensor3D<T>(*this);
}

template<typename T>
void A_lessThan_alpha_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_lessThan_alpha_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_lessThan_alpha_NTensor3D ************************************* */

template<typename T>
Masked_A_lessThan_alpha_NTensor3D<T>::Masked_A_lessThan_alpha_NTensor3D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void Masked_A_lessThan_alpha_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<int>& result,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot3D offset = computeRelativeDisplacement(A, result);
    Dot3D maskOfs = computeRelativeDisplacement(A, mask);

    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        *result.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                            = *A.get(iX,iY,iZ) < scalarAlpha ? 1 : 0;
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            result.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim]
                                = A.get(iX,iY,iZ)[iDim] < alpha[iDim] ? 1 : 0;
                        }
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_lessThan_alpha_NTensor3D<T>* Masked_A_lessThan_alpha_NTensor3D<T>::clone() const {
    return new Masked_A_lessThan_alpha_NTensor3D<T>(*this);
}

template<typename T>
void Masked_A_lessThan_alpha_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_lessThan_alpha_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** A_lessEqual_B_NTensor3D ************************************ */

template<typename T>
void A_lessEqual_B_NTensor3D<T>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==3 );
    NTensorField3D<T>& A = *dynamic_cast<NTensorField3D<T>*>(blocks[0]);
    NTensorField3D<T>& B = *dynamic_cast<NTensorField3D<T>*>(blocks[1]);
    NTensorField3D<int>& result = *dynamic_cast<NTensorField3D<int>*>(blocks[2]);
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    *result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                        = *A.get(iX,iY,iZ) <= *B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z) ? 1:0;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)[iDim]
                            = A.get(iX,iY,iZ)[iDim] <= B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z)[iDim] ? 1:0;
                    }
                }
            }
        }
    }
}

template<typename T>
A_lessEqual_B_NTensor3D<T>* A_lessEqual_B_NTensor3D<T>::clone() const {
    return new A_lessEqual_B_NTensor3D<T>(*this);
}

template<typename T>
void A_lessEqual_B_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_lessEqual_B_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_lessEqual_B_NTensor3D ************************************ */

template<typename T>
void Masked_A_lessEqual_B_NTensor3D<T>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==3 );
    NTensorField3D<T>& A = *dynamic_cast<NTensorField3D<T>*>(blocks[0]);
    NTensorField3D<T>& B = *dynamic_cast<NTensorField3D<T>*>(blocks[1]);
    NTensorField3D<int>& result = *dynamic_cast<NTensorField3D<int>*>(blocks[2]);
    NTensorField3D<int>& mask = *dynamic_cast<NTensorField3D<int>*>(blocks[3]);
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    Dot3D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        *result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                            = *A.get(iX,iY,iZ) <= *B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z) ? 1:0;
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)[iDim]
                                = A.get(iX,iY,iZ)[iDim] <= B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z)[iDim] ? 1:0;
                        }
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_lessEqual_B_NTensor3D<T>* Masked_A_lessEqual_B_NTensor3D<T>::clone() const {
    return new Masked_A_lessEqual_B_NTensor3D<T>(*this);
}

template<typename T>
void Masked_A_lessEqual_B_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
    modified[3] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_lessEqual_B_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** A_lessEqual_alpha_NTensor3D ************************************* */

template<typename T>
A_lessEqual_alpha_NTensor3D<T>::A_lessEqual_alpha_NTensor3D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_lessEqual_alpha_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<int>& result )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot3D offset = computeRelativeDisplacement(A, result);

    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    *result.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                        = *A.get(iX,iY,iZ) <= scalarAlpha ? 1 : 0;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim]
                            = A.get(iX,iY,iZ)[iDim] <= alpha[iDim] ? 1 : 0;
                    }
                }
            }
        }
    }
}

template<typename T>
A_lessEqual_alpha_NTensor3D<T>* A_lessEqual_alpha_NTensor3D<T>::clone() const {
    return new A_lessEqual_alpha_NTensor3D<T>(*this);
}

template<typename T>
void A_lessEqual_alpha_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_lessEqual_alpha_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_lessEqual_alpha_NTensor3D ************************************* */

template<typename T>
Masked_A_lessEqual_alpha_NTensor3D<T>::Masked_A_lessEqual_alpha_NTensor3D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void Masked_A_lessEqual_alpha_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<int>& result,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot3D offset = computeRelativeDisplacement(A, result);
    Dot3D maskOfs = computeRelativeDisplacement(A, mask);

    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        *result.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                            = *A.get(iX,iY,iZ) <= scalarAlpha ? 1 : 0;
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            result.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim]
                                = A.get(iX,iY,iZ)[iDim] <= alpha[iDim] ? 1 : 0;
                        }
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_lessEqual_alpha_NTensor3D<T>* Masked_A_lessEqual_alpha_NTensor3D<T>::clone() const {
    return new Masked_A_lessEqual_alpha_NTensor3D<T>(*this);
}

template<typename T>
void Masked_A_lessEqual_alpha_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_lessEqual_alpha_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** A_greaterThan_B_NTensor3D ************************************ */

template<typename T>
void A_greaterThan_B_NTensor3D<T>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==3 );
    NTensorField3D<T>& A = *dynamic_cast<NTensorField3D<T>*>(blocks[0]);
    NTensorField3D<T>& B = *dynamic_cast<NTensorField3D<T>*>(blocks[1]);
    NTensorField3D<int>& result = *dynamic_cast<NTensorField3D<int>*>(blocks[2]);
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    *result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                        = *A.get(iX,iY,iZ) > *B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z) ? 1:0;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)[iDim]
                            = A.get(iX,iY,iZ)[iDim] > B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z)[iDim] ? 1:0;
                    }
                }
            }
        }
    }
}

template<typename T>
A_greaterThan_B_NTensor3D<T>* A_greaterThan_B_NTensor3D<T>::clone() const {
    return new A_greaterThan_B_NTensor3D<T>(*this);
}

template<typename T>
void A_greaterThan_B_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_greaterThan_B_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_greaterThan_B_NTensor3D ************************************ */

template<typename T>
void Masked_A_greaterThan_B_NTensor3D<T>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==3 );
    NTensorField3D<T>& A = *dynamic_cast<NTensorField3D<T>*>(blocks[0]);
    NTensorField3D<T>& B = *dynamic_cast<NTensorField3D<T>*>(blocks[1]);
    NTensorField3D<int>& result = *dynamic_cast<NTensorField3D<int>*>(blocks[2]);
    NTensorField3D<int>& mask = *dynamic_cast<NTensorField3D<int>*>(blocks[3]);
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    Dot3D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        *result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                            = *A.get(iX,iY,iZ) > *B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z) ? 1:0;
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)[iDim]
                                = A.get(iX,iY,iZ)[iDim] > B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z)[iDim] ? 1:0;
                        }
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_greaterThan_B_NTensor3D<T>* Masked_A_greaterThan_B_NTensor3D<T>::clone() const {
    return new Masked_A_greaterThan_B_NTensor3D<T>(*this);
}

template<typename T>
void Masked_A_greaterThan_B_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
    modified[3] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_greaterThan_B_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** A_greaterThan_alpha_NTensor3D ************************************* */

template<typename T>
A_greaterThan_alpha_NTensor3D<T>::A_greaterThan_alpha_NTensor3D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_greaterThan_alpha_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<int>& result )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot3D offset = computeRelativeDisplacement(A, result);

    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    *result.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                        = *A.get(iX,iY,iZ) > scalarAlpha ? 1 : 0;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim]
                            = A.get(iX,iY,iZ)[iDim] > alpha[iDim] ? 1 : 0;
                    }
                }
            }
        }
    }
}

template<typename T>
A_greaterThan_alpha_NTensor3D<T>* A_greaterThan_alpha_NTensor3D<T>::clone() const {
    return new A_greaterThan_alpha_NTensor3D<T>(*this);
}

template<typename T>
void A_greaterThan_alpha_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_greaterThan_alpha_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_greaterThan_alpha_NTensor3D ************************************* */

template<typename T>
Masked_A_greaterThan_alpha_NTensor3D<T>::Masked_A_greaterThan_alpha_NTensor3D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void Masked_A_greaterThan_alpha_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<int>& result,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot3D offset = computeRelativeDisplacement(A, result);
    Dot3D maskOfs = computeRelativeDisplacement(A, mask);

    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        *result.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                            = *A.get(iX,iY,iZ) > scalarAlpha ? 1 : 0;
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            result.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim]
                                = A.get(iX,iY,iZ)[iDim] > alpha[iDim] ? 1 : 0;
                        }
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_greaterThan_alpha_NTensor3D<T>* Masked_A_greaterThan_alpha_NTensor3D<T>::clone() const {
    return new Masked_A_greaterThan_alpha_NTensor3D<T>(*this);
}

template<typename T>
void Masked_A_greaterThan_alpha_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_greaterThan_alpha_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** A_greaterEqual_B_NTensor3D ************************************ */

template<typename T>
void A_greaterEqual_B_NTensor3D<T>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==3 );
    NTensorField3D<T>& A = *dynamic_cast<NTensorField3D<T>*>(blocks[0]);
    NTensorField3D<T>& B = *dynamic_cast<NTensorField3D<T>*>(blocks[1]);
    NTensorField3D<int>& result = *dynamic_cast<NTensorField3D<int>*>(blocks[2]);
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    *result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                        = *A.get(iX,iY,iZ) >= *B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z) ? 1:0;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)[iDim]
                            = A.get(iX,iY,iZ)[iDim] >= B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z)[iDim] ? 1:0;
                    }
                }
            }
        }
    }
}

template<typename T>
A_greaterEqual_B_NTensor3D<T>* A_greaterEqual_B_NTensor3D<T>::clone() const {
    return new A_greaterEqual_B_NTensor3D<T>(*this);
}

template<typename T>
void A_greaterEqual_B_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_greaterEqual_B_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_greaterEqual_B_NTensor3D ************************************ */

template<typename T>
void Masked_A_greaterEqual_B_NTensor3D<T>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==3 );
    NTensorField3D<T>& A = *dynamic_cast<NTensorField3D<T>*>(blocks[0]);
    NTensorField3D<T>& B = *dynamic_cast<NTensorField3D<T>*>(blocks[1]);
    NTensorField3D<int>& result = *dynamic_cast<NTensorField3D<int>*>(blocks[2]);
    NTensorField3D<int>& mask = *dynamic_cast<NTensorField3D<int>*>(blocks[3]);
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    Dot3D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        *result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                            = *A.get(iX,iY,iZ) >= *B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z) ? 1:0;
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)[iDim]
                                = A.get(iX,iY,iZ)[iDim] >= B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z)[iDim] ? 1:0;
                        }
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_greaterEqual_B_NTensor3D<T>* Masked_A_greaterEqual_B_NTensor3D<T>::clone() const {
    return new Masked_A_greaterEqual_B_NTensor3D<T>(*this);
}

template<typename T>
void Masked_A_greaterEqual_B_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
    modified[3] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_greaterEqual_B_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** A_greaterEqual_alpha_NTensor3D ************************************* */

template<typename T>
A_greaterEqual_alpha_NTensor3D<T>::A_greaterEqual_alpha_NTensor3D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void A_greaterEqual_alpha_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<int>& result )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot3D offset = computeRelativeDisplacement(A, result);

    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    *result.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                        = *A.get(iX,iY,iZ) >= scalarAlpha ? 1 : 0;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim]
                            = A.get(iX,iY,iZ)[iDim] >= alpha[iDim] ? 1 : 0;
                    }
                }
            }
        }
    }
}

template<typename T>
A_greaterEqual_alpha_NTensor3D<T>* A_greaterEqual_alpha_NTensor3D<T>::clone() const {
    return new A_greaterEqual_alpha_NTensor3D<T>(*this);
}

template<typename T>
void A_greaterEqual_alpha_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_greaterEqual_alpha_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_greaterEqual_alpha_NTensor3D ************************************* */

template<typename T>
Masked_A_greaterEqual_alpha_NTensor3D<T>::Masked_A_greaterEqual_alpha_NTensor3D (
        std::vector<T> const& alpha_)
    : alpha(alpha_)
{ }

template<typename T>
void Masked_A_greaterEqual_alpha_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<int>& result,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size() );
    Dot3D offset = computeRelativeDisplacement(A, result);
    Dot3D maskOfs = computeRelativeDisplacement(A, mask);

    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalarAlpha = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        *result.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                            = *A.get(iX,iY,iZ) >= scalarAlpha ? 1 : 0;
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            result.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim]
                                = A.get(iX,iY,iZ)[iDim] >= alpha[iDim] ? 1 : 0;
                        }
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_greaterEqual_alpha_NTensor3D<T>* Masked_A_greaterEqual_alpha_NTensor3D<T>::clone() const {
    return new Masked_A_greaterEqual_alpha_NTensor3D<T>(*this);
}

template<typename T>
void Masked_A_greaterEqual_alpha_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_greaterEqual_alpha_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** A_and_B_NTensor3D ************************************ */

template<typename T>
void A_and_B_NTensor3D<T>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==3 );
    NTensorField3D<T>& A = *dynamic_cast<NTensorField3D<T>*>(blocks[0]);
    NTensorField3D<T>& B = *dynamic_cast<NTensorField3D<T>*>(blocks[1]);
    NTensorField3D<int>& result = *dynamic_cast<NTensorField3D<int>*>(blocks[2]);
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    *result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                        = (bool)(*A.get(iX,iY,iZ)) && (bool)(*B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z)) ? 1:0;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)[iDim]
                            = (bool)A.get(iX,iY,iZ)[iDim] && (bool)B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z)[iDim] ? 1:0;
                    }
                }
            }
        }
    }
}

template<typename T>
A_and_B_NTensor3D<T>* A_and_B_NTensor3D<T>::clone() const {
    return new A_and_B_NTensor3D<T>(*this);
}

template<typename T>
void A_and_B_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_and_B_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_and_B_NTensor3D ************************************ */

template<typename T>
void Masked_A_and_B_NTensor3D<T>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==3 );
    NTensorField3D<T>& A = *dynamic_cast<NTensorField3D<T>*>(blocks[0]);
    NTensorField3D<T>& B = *dynamic_cast<NTensorField3D<T>*>(blocks[1]);
    NTensorField3D<int>& result = *dynamic_cast<NTensorField3D<int>*>(blocks[2]);
    NTensorField3D<int>& mask = *dynamic_cast<NTensorField3D<int>*>(blocks[3]);
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    Dot3D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        *result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                            = (bool)(*A.get(iX,iY,iZ)) && (bool)(*B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z)) ? 1:0;
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)[iDim]
                                = (bool)A.get(iX,iY,iZ)[iDim] && (bool)B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z)[iDim] ? 1:0;
                        }
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_and_B_NTensor3D<T>* Masked_A_and_B_NTensor3D<T>::clone() const {
    return new Masked_A_and_B_NTensor3D<T>(*this);
}

template<typename T>
void Masked_A_and_B_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
    modified[3] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_and_B_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}



/* ******** A_or_B_NTensor3D ************************************ */

template<typename T>
void A_or_B_NTensor3D<T>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==3 );
    NTensorField3D<T>& A = *dynamic_cast<NTensorField3D<T>*>(blocks[0]);
    NTensorField3D<T>& B = *dynamic_cast<NTensorField3D<T>*>(blocks[1]);
    NTensorField3D<int>& result = *dynamic_cast<NTensorField3D<int>*>(blocks[2]);
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    *result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                        = (bool)(*A.get(iX,iY,iZ)) || (bool)(*B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z)) ? 1:0;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)[iDim]
                            = (bool)A.get(iX,iY,iZ)[iDim] || (bool)B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z)[iDim] ? 1:0;
                    }
                }
            }
        }
    }
}

template<typename T>
A_or_B_NTensor3D<T>* A_or_B_NTensor3D<T>::clone() const {
    return new A_or_B_NTensor3D<T>(*this);
}

template<typename T>
void A_or_B_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_or_B_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_A_or_B_NTensor3D ************************************ */

template<typename T>
void Masked_A_or_B_NTensor3D<T>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> blocks )
{
    PLB_PRECONDITION( blocks.size()==3 );
    NTensorField3D<T>& A = *dynamic_cast<NTensorField3D<T>*>(blocks[0]);
    NTensorField3D<T>& B = *dynamic_cast<NTensorField3D<T>*>(blocks[1]);
    NTensorField3D<int>& result = *dynamic_cast<NTensorField3D<int>*>(blocks[2]);
    NTensorField3D<int>& mask = *dynamic_cast<NTensorField3D<int>*>(blocks[3]);
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );

    Dot3D offsetB      = computeRelativeDisplacement(A,B);
    Dot3D offsetResult = computeRelativeDisplacement(A,result);
    Dot3D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        *result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)
                            = (bool)(*A.get(iX,iY,iZ)) || (bool)(*B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z)) ? 1:0;
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            result.get(iX+offsetResult.x,iY+offsetResult.y,iZ+offsetResult.z)[iDim]
                                = (bool)A.get(iX,iY,iZ)[iDim] || (bool)B.get(iX+offsetB.x,iY+offsetB.y,iZ+offsetB.z)[iDim] ? 1:0;
                        }
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_or_B_NTensor3D<T>* Masked_A_or_B_NTensor3D<T>::clone() const {
    return new Masked_A_or_B_NTensor3D<T>(*this);
}

template<typename T>
void Masked_A_or_B_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::staticVariables;
    modified[3] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_or_B_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Not_A_NTensor3D ************************************* */

template<typename T>
void Not_A_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<int>& result )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    Dot3D offset = computeRelativeDisplacement(A, result);

    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    *result.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                        = *A.get(iX,iY,iZ) == T() ? 1 : 0;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        result.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim]
                            = A.get(iX,iY,iZ)[iDim] == T() ? 1 : 0;
                    }
                }
            }
        }
    }
}

template<typename T>
Not_A_NTensor3D<T>* Not_A_NTensor3D<T>::clone() const {
    return new Not_A_NTensor3D<T>(*this);
}

template<typename T>
void Not_A_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT Not_A_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}


/* ******** Masked_Not_A_NTensor3D ************************************* */

template<typename T>
void Masked_Not_A_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<int>& result,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == result.getNdim() );
    Dot3D offset = computeRelativeDisplacement(A, result);
    Dot3D maskOfs = computeRelativeDisplacement(A, mask);

    plint ndim = A.getNdim();
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        *result.get(iX+offset.x,iY+offset.y,iZ+offset.z)
                            = *A.get(iX,iY,iZ) == T() ? 1 : 0;
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            result.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim]
                                = A.get(iX,iY,iZ)[iDim] == T() ? 1 : 0;
                        }
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_Not_A_NTensor3D<T>* Masked_Not_A_NTensor3D<T>::clone() const {
    return new Masked_Not_A_NTensor3D<T>(*this);
}

template<typename T>
void Masked_Not_A_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_Not_A_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulk;
}



/* ******** A_plus_B_inplace_NTensor3D ************************************ */

template<typename T>
void A_plus_B_inplace_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<T>& B)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    plint ndim = A.getNdim();

    Dot3D offset = computeRelativeDisplacement(A,B);
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    *A.get(iX,iY,iZ) += *B.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        A.get(iX,iY,iZ)[iDim] += B.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
A_plus_B_inplace_NTensor3D<T>* A_plus_B_inplace_NTensor3D<T>::clone() const {
    return new A_plus_B_inplace_NTensor3D<T>(*this);
}

template<typename T>
void A_plus_B_inplace_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT A_plus_B_inplace_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Masked_A_plus_B_inplace_NTensor3D ************************************ */

template<typename T>
void Masked_A_plus_B_inplace_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<T>& B,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    plint ndim = A.getNdim();

    Dot3D offset = computeRelativeDisplacement(A,B);
    Dot3D maskOfs = computeRelativeDisplacement(A, mask);
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        *A.get(iX,iY,iZ) += *B.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            A.get(iX,iY,iZ)[iDim] += B.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim];
                        }
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_plus_B_inplace_NTensor3D<T>* Masked_A_plus_B_inplace_NTensor3D<T>::clone() const {
    return new Masked_A_plus_B_inplace_NTensor3D<T>(*this);
}

template<typename T>
void Masked_A_plus_B_inplace_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_plus_B_inplace_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_plus_alpha_inplace_NTensor3D ************************************ */

template<typename T>
A_plus_alpha_inplace_NTensor3D<T>::A_plus_alpha_inplace_NTensor3D (
        std::vector<T> const& alpha_ )
    : alpha(alpha_)
{ }

template<typename T>
void A_plus_alpha_inplace_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A)
{
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size());
    plint ndim = A.getNdim();

    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalar = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    *A.get(iX,iY,iZ) += scalar;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        A.get(iX,iY,iZ)[iDim] += alpha[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
A_plus_alpha_inplace_NTensor3D<T>* A_plus_alpha_inplace_NTensor3D<T>::clone() const {
    return new A_plus_alpha_inplace_NTensor3D<T>(*this);
}

template<typename T>
void A_plus_alpha_inplace_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_plus_alpha_inplace_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Masked_A_plus_alpha_inplace_NTensor3D ************************************ */

template<typename T>
Masked_A_plus_alpha_inplace_NTensor3D<T>::Masked_A_plus_alpha_inplace_NTensor3D (
        std::vector<T> const& alpha_ )
    : alpha(alpha_)
{ }

template<typename T>
void Masked_A_plus_alpha_inplace_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<int>& mask)
{
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size());
    Dot3D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();

    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalar = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        *A.get(iX,iY,iZ) += scalar;
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            A.get(iX,iY,iZ)[iDim] += alpha[iDim];
                        }
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_plus_alpha_inplace_NTensor3D<T>* Masked_A_plus_alpha_inplace_NTensor3D<T>::clone() const {
    return new Masked_A_plus_alpha_inplace_NTensor3D<T>(*this);
}

template<typename T>
void Masked_A_plus_alpha_inplace_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_plus_alpha_inplace_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_minus_B_inplace_NTensor3D ************************************ */

template<typename T>
void A_minus_B_inplace_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<T>& B)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    plint ndim = A.getNdim();

    Dot3D offset = computeRelativeDisplacement(A,B);
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    *A.get(iX,iY,iZ) -= *B.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        A.get(iX,iY,iZ)[iDim] -= B.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
A_minus_B_inplace_NTensor3D<T>* A_minus_B_inplace_NTensor3D<T>::clone() const {
    return new A_minus_B_inplace_NTensor3D<T>(*this);
}

template<typename T>
void A_minus_B_inplace_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT A_minus_B_inplace_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

/* ******** Masked_A_minus_B_inplace_NTensor3D ************************************ */

template<typename T>
void Masked_A_minus_B_inplace_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<T>& B,
        NTensorField3D<int>& mask)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    Dot3D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();

    Dot3D offset = computeRelativeDisplacement(A,B);
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        *A.get(iX,iY,iZ) -= *B.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            A.get(iX,iY,iZ)[iDim] -= B.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim];
                        }
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_minus_B_inplace_NTensor3D<T>* Masked_A_minus_B_inplace_NTensor3D<T>::clone() const {
    return new Masked_A_minus_B_inplace_NTensor3D<T>(*this);
}

template<typename T>
void Masked_A_minus_B_inplace_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_minus_B_inplace_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

/* ******** A_minus_alpha_inplace_NTensor3D ************************************ */

template<typename T>
A_minus_alpha_inplace_NTensor3D<T>::A_minus_alpha_inplace_NTensor3D (
        std::vector<T> const& alpha_ )
    : alpha(alpha_)
{ }

template<typename T>
void A_minus_alpha_inplace_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A)
{
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size());
    plint ndim = A.getNdim();

    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalar = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    *A.get(iX,iY,iZ) -= scalar;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        A.get(iX,iY,iZ)[iDim] -= alpha[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
A_minus_alpha_inplace_NTensor3D<T>* A_minus_alpha_inplace_NTensor3D<T>::clone() const {
    return new A_minus_alpha_inplace_NTensor3D<T>(*this);
}

template<typename T>
void A_minus_alpha_inplace_NTensor3D<T>::
    getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_minus_alpha_inplace_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}



/* ******** Masked_A_minus_alpha_inplace_NTensor3D ************************************ */

template<typename T>
Masked_A_minus_alpha_inplace_NTensor3D<T>::Masked_A_minus_alpha_inplace_NTensor3D (
        std::vector<T> const& alpha_ )
    : alpha(alpha_)
{ }

template<typename T>
void Masked_A_minus_alpha_inplace_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size());
    Dot3D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();

    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalar = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        *A.get(iX,iY,iZ) -= scalar;
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            A.get(iX,iY,iZ)[iDim] -= alpha[iDim];
                        }
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_minus_alpha_inplace_NTensor3D<T>* Masked_A_minus_alpha_inplace_NTensor3D<T>::clone() const {
    return new Masked_A_minus_alpha_inplace_NTensor3D<T>(*this);
}

template<typename T>
void Masked_A_minus_alpha_inplace_NTensor3D<T>::
    getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_minus_alpha_inplace_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}



/* ******** A_times_B_inplace_NTensor3D ************************************ */

template<typename T>
void A_times_B_inplace_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<T>& B)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    plint ndim = A.getNdim();

    Dot3D offset = computeRelativeDisplacement(A,B);
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    *A.get(iX,iY,iZ) *= *B.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        A.get(iX,iY,iZ)[iDim] *= B.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
A_times_B_inplace_NTensor3D<T>* A_times_B_inplace_NTensor3D<T>::clone() const {
    return new A_times_B_inplace_NTensor3D<T>(*this);
}

template<typename T>
void A_times_B_inplace_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT A_times_B_inplace_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Masked_A_times_B_inplace_NTensor3D ************************************ */

template<typename T>
void Masked_A_times_B_inplace_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<T>& B,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    plint ndim = A.getNdim();

    Dot3D offset = computeRelativeDisplacement(A,B);
    Dot3D maskOfs = computeRelativeDisplacement(A, mask);
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        *A.get(iX,iY,iZ) *= *B.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            A.get(iX,iY,iZ)[iDim] *= B.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim];
                        }
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_times_B_inplace_NTensor3D<T>* Masked_A_times_B_inplace_NTensor3D<T>::clone() const {
    return new Masked_A_times_B_inplace_NTensor3D<T>(*this);
}

template<typename T>
void Masked_A_times_B_inplace_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_times_B_inplace_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_times_alpha_inplace_NTensor3D ************************************ */

template<typename T>
A_times_alpha_inplace_NTensor3D<T>::A_times_alpha_inplace_NTensor3D (
        std::vector<T> const& alpha_ )
    : alpha(alpha_)
{ }

template<typename T>
void A_times_alpha_inplace_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A)
{
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size());
    plint ndim = A.getNdim();

    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalar = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    *A.get(iX,iY,iZ) *= scalar;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        A.get(iX,iY,iZ)[iDim] *= alpha[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
A_times_alpha_inplace_NTensor3D<T>* A_times_alpha_inplace_NTensor3D<T>::clone() const {
    return new A_times_alpha_inplace_NTensor3D<T>(*this);
}

template<typename T>
void A_times_alpha_inplace_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_times_alpha_inplace_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}



/* ******** Masked_A_times_alpha_inplace_NTensor3D ************************************ */

template<typename T>
Masked_A_times_alpha_inplace_NTensor3D<T>::Masked_A_times_alpha_inplace_NTensor3D (
        std::vector<T> const& alpha_ )
    : alpha(alpha_)
{ }

template<typename T>
void Masked_A_times_alpha_inplace_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<int>& mask)
{
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size());
    Dot3D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();

    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalar = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        *A.get(iX,iY,iZ) *= scalar;
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            A.get(iX,iY,iZ)[iDim] *= alpha[iDim];
                        }
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_times_alpha_inplace_NTensor3D<T>* Masked_A_times_alpha_inplace_NTensor3D<T>::clone() const {
    return new Masked_A_times_alpha_inplace_NTensor3D<T>(*this);
}

template<typename T>
void Masked_A_times_alpha_inplace_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_times_alpha_inplace_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}



/* ******** A_dividedBy_B_inplace_NTensor3D ************************************ */

template<typename T>
void A_dividedBy_B_inplace_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<T>& B)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    plint ndim = A.getNdim();

    Dot3D offset = computeRelativeDisplacement(A,B);
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    *A.get(iX,iY,iZ) /= *B.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        A.get(iX,iY,iZ)[iDim] /= B.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
A_dividedBy_B_inplace_NTensor3D<T>* A_dividedBy_B_inplace_NTensor3D<T>::clone() const {
    return new A_dividedBy_B_inplace_NTensor3D<T>(*this);
}

template<typename T>
void A_dividedBy_B_inplace_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT A_dividedBy_B_inplace_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Masked_A_dividedBy_B_inplace_NTensor3D ************************************ */

template<typename T>
void Masked_A_dividedBy_B_inplace_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<T>& B,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    Dot3D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();

    Dot3D offset = computeRelativeDisplacement(A,B);
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        *A.get(iX,iY,iZ) /= *B.get(iX+offset.x,iY+offset.y,iZ+offset.z);
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            A.get(iX,iY,iZ)[iDim] /= B.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim];
                        }
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_dividedBy_B_inplace_NTensor3D<T>* Masked_A_dividedBy_B_inplace_NTensor3D<T>::clone() const {
    return new Masked_A_dividedBy_B_inplace_NTensor3D<T>(*this);
}

template<typename T>
void Masked_A_dividedBy_B_inplace_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_dividedBy_B_inplace_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_dividedBy_alpha_inplace_NTensor3D ************************************ */

template<typename T>
A_dividedBy_alpha_inplace_NTensor3D<T>::A_dividedBy_alpha_inplace_NTensor3D (
        std::vector<T> const& alpha_ )
    : alpha(alpha_)
{ }

template<typename T>
void A_dividedBy_alpha_inplace_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A)
{
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size());
    plint ndim = A.getNdim();

    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalar = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    *A.get(iX,iY,iZ) /= scalar;
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        A.get(iX,iY,iZ)[iDim] /= alpha[iDim];
                    }
                }
            }
        }
    }
}

template<typename T>
A_dividedBy_alpha_inplace_NTensor3D<T>* A_dividedBy_alpha_inplace_NTensor3D<T>::clone() const {
    return new A_dividedBy_alpha_inplace_NTensor3D<T>(*this);
}

template<typename T>
void A_dividedBy_alpha_inplace_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_dividedBy_alpha_inplace_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Masked_A_dividedBy_alpha_inplace_NTensor3D ************************************ */

template<typename T>
Masked_A_dividedBy_alpha_inplace_NTensor3D<T>::Masked_A_dividedBy_alpha_inplace_NTensor3D (
        std::vector<T> const& alpha_ )
    : alpha(alpha_)
{ }

template<typename T>
void Masked_A_dividedBy_alpha_inplace_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size());
    Dot3D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();

    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalar = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        *A.get(iX,iY,iZ) /= scalar;
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            A.get(iX,iY,iZ)[iDim] /= alpha[iDim];
                        }
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_dividedBy_alpha_inplace_NTensor3D<T>* Masked_A_dividedBy_alpha_inplace_NTensor3D<T>::clone() const {
    return new Masked_A_dividedBy_alpha_inplace_NTensor3D<T>(*this);
}

template<typename T>
void Masked_A_dividedBy_alpha_inplace_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_dividedBy_alpha_inplace_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_toThePower_B_inplace_NTensor3D ************************************ */

template<typename T>
void A_toThePower_B_inplace_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<T>& B)
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    plint ndim = A.getNdim();

    Dot3D offset = computeRelativeDisplacement(A,B);
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    customInPlacePower(*A.get(iX,iY,iZ),*B.get(iX+offset.x,iY+offset.y,iZ+offset.z));
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        customInPlacePower(A.get(iX,iY,iZ)[iDim],B.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim]);
                    }
                }
            }
        }
    }
}

template<typename T>
A_toThePower_B_inplace_NTensor3D<T>* A_toThePower_B_inplace_NTensor3D<T>::clone() const {
    return new A_toThePower_B_inplace_NTensor3D<T>(*this);
}

template<typename T>
void A_toThePower_B_inplace_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT A_toThePower_B_inplace_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** Masked_A_toThePower_B_inplace_NTensor3D ************************************ */

template<typename T>
void Masked_A_toThePower_B_inplace_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<T>& B,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == B.getNdim() );
    Dot3D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();

    Dot3D offset = computeRelativeDisplacement(A,B);
    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        customInPlacePower(*A.get(iX,iY,iZ),*B.get(iX+offset.x,iY+offset.y,iZ+offset.z));
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            customInPlacePower(A.get(iX,iY,iZ)[iDim],B.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim]);
                        }
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_toThePower_B_inplace_NTensor3D<T>* Masked_A_toThePower_B_inplace_NTensor3D<T>::clone() const {
    return new Masked_A_toThePower_B_inplace_NTensor3D<T>(*this);
}

template<typename T>
void Masked_A_toThePower_B_inplace_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_toThePower_B_inplace_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* ******** A_toThePower_alpha_inplace_NTensor3D ************************************ */

template<typename T>
A_toThePower_alpha_inplace_NTensor3D<T>::A_toThePower_alpha_inplace_NTensor3D (
        std::vector<T> const& alpha_ )
    : alpha(alpha_)
{ }

template<typename T>
void A_toThePower_alpha_inplace_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A)
{
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size());
    plint ndim = A.getNdim();

    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalar = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    customInPlacePower(*A.get(iX,iY,iZ),scalar);
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        customInPlacePower(A.get(iX,iY,iZ)[iDim],alpha[iDim]);
                    }
                }
            }
        }
    }
}

template<typename T>
A_toThePower_alpha_inplace_NTensor3D<T>* A_toThePower_alpha_inplace_NTensor3D<T>::clone() const {
    return new A_toThePower_alpha_inplace_NTensor3D<T>(*this);
}

template<typename T>
void A_toThePower_alpha_inplace_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
}

template<typename T>
BlockDomain::DomainT A_toThePower_alpha_inplace_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}

/* ******** Masked_A_toThePower_alpha_inplace_NTensor3D ************************************ */

template<typename T>
Masked_A_toThePower_alpha_inplace_NTensor3D<T>::Masked_A_toThePower_alpha_inplace_NTensor3D (
        std::vector<T> const& alpha_ )
    : alpha(alpha_)
{ }

template<typename T>
void Masked_A_toThePower_alpha_inplace_NTensor3D<T>::process (
        Box3D domain, NTensorField3D<T>& A, NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( A.getNdim() == (plint)alpha.size());
    Dot3D maskOfs = computeRelativeDisplacement(A, mask);
    plint ndim = A.getNdim();

    // Improve computational speed if the field is scalar-valued.
    if (ndim==1) {
        T scalar = alpha[0];
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        customInPlacePower(*A.get(iX,iY,iZ),scalar);
                    }
                }
            }
        }
    }
    // Generic implementation for any number of dimensions.
    else {
        for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
            for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
                for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                    if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                        for (plint iDim=0; iDim<ndim; ++iDim) {
                            customInPlacePower(A.get(iX,iY,iZ)[iDim],alpha[iDim]);
                        }
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_A_toThePower_alpha_inplace_NTensor3D<T>* Masked_A_toThePower_alpha_inplace_NTensor3D<T>::clone() const {
    return new Masked_A_toThePower_alpha_inplace_NTensor3D<T>(*this);
}

template<typename T>
void Masked_A_toThePower_alpha_inplace_NTensor3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const {
    modified[0] = modif::staticVariables;
    modified[1] = modif::nothing;
}

template<typename T>
BlockDomain::DomainT Masked_A_toThePower_alpha_inplace_NTensor3D<T>::appliesTo() const {
    return BlockDomain::bulkAndEnvelope;
}


/* *************** UPO  ******** */

template<typename T>
UPO_ScalarProductFunctional3D<T>::UPO_ScalarProductFunctional3D()
    : sumId(this->getStatistics().subscribeSum())
{ }

template<typename T>
void UPO_ScalarProductFunctional3D<T>::process (
        Box3D domain,
        NTensorField3D<T>& a,
        NTensorField3D<T>& b )
{
    PLB_PRECONDITION( a.getNdim()==b.getNdim() );
    Dot3D offset = computeRelativeDisplacement(a, b);
    plint ndim = a.getNdim();
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                for (plint iDim=0; iDim<ndim; ++iDim) {
                    statistics.gatherSum( sumId,
                            (double)( a.get(iX,iY,iZ)[iDim] *
                                      b.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim] ) );
                }
            }
        }
    }
}

template<typename T>
UPO_ScalarProductFunctional3D<T>* UPO_ScalarProductFunctional3D<T>::clone() const
{
    return new UPO_ScalarProductFunctional3D<T>(*this);
}

template<typename T>
void UPO_ScalarProductFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
}

template<typename T>
T UPO_ScalarProductFunctional3D<T>::getSum() const {
    T sum;
    double doubleSum = this->getStatistics().getSum(sumId);
    // The sum is internally computed on floating-point values. If T is
    //   integer, the value must be rounded at the end.
    if (std::numeric_limits<T>::is_integer) {
        sum = (T) util::roundToInt(doubleSum);
    }
    else {
        sum = (T) doubleSum;
    }
    return sum;
}


template<typename T>
Masked_UPO_ScalarProductFunctional3D<T>::Masked_UPO_ScalarProductFunctional3D()
    : sumId(this->getStatistics().subscribeSum())
{ }

template<typename T>
void Masked_UPO_ScalarProductFunctional3D<T>::process (
        Box3D domain,
        NTensorField3D<T>& a,
        NTensorField3D<T>& b,
        NTensorField3D<int>& mask )
{
    PLB_PRECONDITION( a.getNdim()==b.getNdim() );
    Dot3D offset = computeRelativeDisplacement(a, b);
    Dot3D maskOfs = computeRelativeDisplacement(a, mask);
    plint ndim = a.getNdim();
    BlockStatistics& statistics = this->getStatistics();
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                if (*mask.get(iX+maskOfs.x, iY+maskOfs.y, iZ+maskOfs.z)) {
                    for (plint iDim=0; iDim<ndim; ++iDim) {
                        statistics.gatherSum( sumId,
                                (double)( a.get(iX,iY,iZ)[iDim] *
                                          b.get(iX+offset.x,iY+offset.y,iZ+offset.z)[iDim] ) );
                    }
                }
            }
        }
    }
}

template<typename T>
Masked_UPO_ScalarProductFunctional3D<T>* Masked_UPO_ScalarProductFunctional3D<T>::clone() const
{
    return new Masked_UPO_ScalarProductFunctional3D<T>(*this);
}

template<typename T>
void Masked_UPO_ScalarProductFunctional3D<T>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::nothing;
    modified[1] = modif::nothing;
    modified[2] = modif::nothing;
}

template<typename T>
T Masked_UPO_ScalarProductFunctional3D<T>::getSum() const {
    T sum;
    double doubleSum = this->getStatistics().getSum(sumId);
    // The sum is internally computed on floating-point values. If T is
    //   integer, the value must be rounded at the end.
    if (std::numeric_limits<T>::is_integer) {
        sum = (T) util::roundToInt(doubleSum);
    }
    else {
        sum = (T) doubleSum;
    }
    return sum;
}

}  // namespace plb

#endif  // DATA_ANALYSIS_FUNCTIONAL_3D_HH
