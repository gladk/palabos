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

/** \file
 * Serial implementation of scalar, vector and tensor fields for 3D data analysis.
 * -- header file
 */

#ifndef DATA_FIELD_3D_H
#define DATA_FIELD_3D_H

#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include "atomicBlock/dataField2D.h"
#include "core/dataFieldBase2D.h"
#include "core/dataFieldBase3D.h"
#include "atomicBlock/atomicBlock3D.h"

namespace plb {

template<typename T> class ScalarField3D;

template<typename T>
class ScalarFieldDataTransfer3D : public BlockDataTransfer3D {
public:
    ScalarFieldDataTransfer3D(ScalarField3D<T>& field_);
    virtual plint staticCellSize() const;
    /// Send data from the block into a byte-stream.
    virtual void send(Box3D domain, std::vector<char>& buffer, modif::ModifT kind) const;
    /// Receive data from a byte-stream into the block.
    virtual void receive(Box3D domain, std::vector<char> const& buffer, modif::ModifT kind );
    virtual void receive(Box3D domain, std::vector<char> const& buffer, modif::ModifT kind, Dot3D absoluteOffset) {
        receive(domain, buffer, kind);
    }
    virtual void receive(Box3D domain, std::vector<char> const& buffer,
                         modif::ModifT kind, std::map<int,std::string> const& foreignIds )
    {
        receive(domain, buffer, kind);
    }
    /// Attribute data between two blocks.
    virtual void attribute(Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
                           AtomicBlock3D const& from, modif::ModifT kind);
    virtual void attribute(Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
                           AtomicBlock3D const& from, modif::ModifT kind, Dot3D absoluteOffset)
    {
        attribute(toDomain, deltaX, deltaY, deltaZ, from, kind);
    }
private:
    ScalarField3D<T>& field;
};


template<typename T>
class ScalarField3D : public ScalarFieldBase3D<T>, public AtomicBlock3D {
public:
    ScalarField3D(plint nx_, plint ny_, plint nz_, T iniVal=T());
    ~ScalarField3D();
    ScalarField3D(ScalarField3D<T> const& rhs);
    ScalarField3D<T>& operator=(ScalarField3D<T> const& rhs);
    void swap(ScalarField3D<T>& rhs);
public:
    virtual void reset();
    virtual pluint getSize() const { return (pluint)this->getNx()*(pluint)this->getNy()*(pluint)this->getNz(); }
    virtual T& get(plint iX, plint iY, plint iZ) {
        PLB_PRECONDITION(iX>=0 && iX<this->getNx());
        PLB_PRECONDITION(iY>=0 && iY<this->getNy());
        PLB_PRECONDITION(iZ>=0 && iZ<this->getNz());
        return field[iX][iY][iZ];
    }
    virtual T const& get(plint iX, plint iY, plint iZ) const {
        PLB_PRECONDITION(iX>=0 && iX<this->getNx());
        PLB_PRECONDITION(iY>=0 && iY<this->getNy());
        PLB_PRECONDITION(iZ>=0 && iZ<this->getNz());
        return field[iX][iY][iZ];
    }
    T& operator[] (plint ind) {
        PLB_PRECONDITION(ind>=0 && ind<this->getNx()*this->getNy()*this->getNz());
        return rawData[ind];
    }
    T const& operator[] (plint ind) const {
        PLB_PRECONDITION(ind>=0 && ind<this->getNx()*this->getNy()*this->getNz());
        return rawData[ind];
    }
    /// Get access to data transfer between blocks
    virtual ScalarFieldDataTransfer3D<T>& getDataTransfer();
    /// Get access to data transfer between blocks (const version)
    virtual ScalarFieldDataTransfer3D<T> const& getDataTransfer() const;
private:
    void allocateMemory();
    void releaseMemory();
private:
    T   *rawData;
    T   ***field;
    ScalarFieldDataTransfer3D<T> dataTransfer;
};

template<typename T, int nDim> class TensorField3D;

template<typename T, int nDim>
class TensorFieldDataTransfer3D : public BlockDataTransfer3D {
public:
    TensorFieldDataTransfer3D(TensorField3D<T,nDim>& field_);
    virtual plint staticCellSize() const;
    /// Send data from the block into a byte-stream.
    virtual void send(Box3D domain, std::vector<char>& buffer, modif::ModifT kind) const;
    /// Receive data from a byte-stream into the block.
    virtual void receive(Box3D domain, std::vector<char> const& buffer, modif::ModifT kind);
    virtual void receive(Box3D domain, std::vector<char> const& buffer, modif::ModifT kind, Dot3D absoluteOffset) {
        receive(domain, buffer, kind);
    }
    virtual void receive(Box3D domain, std::vector<char> const& buffer,
                         modif::ModifT kind, std::map<int,std::string> const& foreignIds )
    {
        receive(domain, buffer, kind);
    }
    /// Attribute data between two blocks.
    virtual void attribute(Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
                           AtomicBlock3D const& from, modif::ModifT kind);
    virtual void attribute(Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
                           AtomicBlock3D const& from, modif::ModifT kind, Dot3D absoluteOffset)
    {
        attribute(toDomain, deltaX, deltaY, deltaZ, from, kind);
    }
private:
    TensorField3D<T,nDim>& field;
};


template<typename T, int nDim>
class TensorField3D : public TensorFieldBase3D<T,nDim>, public AtomicBlock3D {
public:
    TensorField3D(plint nx_, plint ny_, plint nz_);
    TensorField3D(plint nx_, plint ny_, plint nz_, Array<T,nDim> const& iniVal);
    ~TensorField3D();
    TensorField3D(TensorField3D<T,nDim> const& rhs);
    TensorField3D<T,nDim>& operator=(TensorField3D<T,nDim> const& rhs);
    void swap(TensorField3D<T,nDim>& rhs);
public:
    virtual void reset();
    virtual Array<T,nDim>& get(plint iX, plint iY, plint iZ) {
        PLB_PRECONDITION(iX>=0 && iX<this->getNx());
        PLB_PRECONDITION(iY>=0 && iY<this->getNy());
        PLB_PRECONDITION(iZ>=0 && iZ<this->getNz());
        return field[iX][iY][iZ];
    }
    virtual Array<T,nDim> const& get(plint iX, plint iY, plint iZ) const {
        PLB_PRECONDITION(iX>=0 && iX<this->getNx());
        PLB_PRECONDITION(iY>=0 && iY<this->getNy());
        PLB_PRECONDITION(iZ>=0 && iZ<this->getNz());
        return field[iX][iY][iZ];
    }
    Array<T,nDim>& operator[] (plint ind) {
        PLB_PRECONDITION(ind>=0 && ind<this->getNx()*this->getNy()*this->getNz());
        return rawData[ind];
    }
    Array<T,nDim> const& operator[] (plint ind) const {
        PLB_PRECONDITION(ind>=0 && ind<this->getNx()*this->getNy()*this->getNz());
        return rawData[ind];
    }
    /// Get access to data transfer between blocks
    virtual TensorFieldDataTransfer3D<T,nDim>& getDataTransfer();
    /// Get access to data transfer between blocks (const version)
    virtual TensorFieldDataTransfer3D<T,nDim> const& getDataTransfer() const;
private:
    void allocateMemory();
    void releaseMemory();
private:
    Array<T,nDim> *rawData;
    Array<T,nDim> ***field;
    TensorFieldDataTransfer3D<T,nDim> dataTransfer;
};


template<typename T> class NTensorField3D;

template<typename T>
class NTensorFieldDataTransfer3D : public BlockDataTransfer3D {
public:
    NTensorFieldDataTransfer3D(NTensorField3D<T>& field_);
    virtual plint staticCellSize() const;
    /// Send data from the block into a byte-stream.
    virtual void send(Box3D domain, std::vector<char>& buffer, modif::ModifT kind) const;
    /// Receive data from a byte-stream into the block.
    virtual void receive(Box3D domain, std::vector<char> const& buffer, modif::ModifT kind);
    virtual void receive(Box3D domain, std::vector<char> const& buffer, modif::ModifT kind, Dot3D absoluteOffset) {
        receive(domain, buffer, kind);
    }
    virtual void receive(Box3D domain, std::vector<char> const& buffer,
                         modif::ModifT kind, std::map<int,std::string> const& foreignIds )
    {
        receive(domain, buffer, kind);
    }
    /// Attribute data between two blocks.
    virtual void attribute(Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
                           AtomicBlock3D const& from, modif::ModifT kind);
    virtual void attribute(Box3D toDomain, plint deltaX, plint deltaY, plint deltaZ,
                           AtomicBlock3D const& from, modif::ModifT kind, Dot3D absoluteOffset)
    {
        attribute(toDomain, deltaX, deltaY, deltaZ, from, kind);
    }
private:
    NTensorField3D<T>& field;
};

template<typename T>
class NTensorField3D : public NTensorFieldBase3D<T>, public AtomicBlock3D {
public:
    NTensorField3D(plint nx_, plint ny_, plint nz_, plint ndim_);
    NTensorField3D(plint nx_, plint ny_, plint nz_, plint ndim_, T const* iniVal);
    ~NTensorField3D();
    NTensorField3D(NTensorField3D<T> const& rhs);
    NTensorField3D<T>& operator=(NTensorField3D<T> const& rhs);
    void swap(NTensorField3D<T>& rhs);
public:
    virtual void reset();
    virtual T* get(plint iX, plint iY, plint iZ) {
        PLB_PRECONDITION(iX>=0 && iX<this->getNx());
        PLB_PRECONDITION(iY>=0 && iY<this->getNy());
        PLB_PRECONDITION(iZ>=0 && iZ<this->getNz());
        return field[iX][iY][iZ];
    }
    virtual T const* get(plint iX, plint iY, plint iZ) const {
        PLB_PRECONDITION(iX>=0 && iX<this->getNx());
        PLB_PRECONDITION(iY>=0 && iY<this->getNy());
        PLB_PRECONDITION(iZ>=0 && iZ<this->getNz());
        return field[iX][iY][iZ];
    }
    T& operator[] (plint ind) {
        PLB_PRECONDITION(ind>=0 && ind<this->getNx()*this->getNy()*this->getNz()*this->getNdim());
        return rawData[ind];
    }
    T const& operator[] (plint ind) const {
        PLB_PRECONDITION(ind>=0 && ind<this->getNx()*this->getNy()*this->getNz()*this->getNdim());
        return rawData[ind];
    }
    /// Get access to data transfer between blocks
    virtual NTensorFieldDataTransfer3D<T>& getDataTransfer();
    /// Get access to data transfer between blocks (const version)
    virtual NTensorFieldDataTransfer3D<T> const& getDataTransfer() const;
private:
    void allocateMemory();
    void releaseMemory();
private:
    T *rawData;
    T ****field;
    NTensorFieldDataTransfer3D<T> dataTransfer;
};

}  // namespace plb


#endif  // DATA_FIELD_3D_H
