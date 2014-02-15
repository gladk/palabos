/*
 * multiBlock/multiDataField.i
 */

%{
#include "PALABOS_ROOT/src/multiBlock/multiDataField2D.h"
#include "JLABOS_ROOT/plbWrapper/block/multiBlockInfo2D.h"
#include "PALABOS_ROOT/src/multiBlock/multiDataField3D.h"
#include "JLABOS_ROOT/plbWrapper/block/multiBlockInfo3D.h"
%}

namespace plb {
        typedef long plint;
            typedef size_t pluint;
template<typename T>
class MultiNTensorField2D {
public:
    Box2D getBoundingBox() const;
    plint getNdim() const;
    void initialize();
private:
    MultiNTensorField2D();
};

template<typename T>
class MultiNTensorInfo2D {
public:
    MultiNTensorInfo2D(MultiNTensorField2D<T> const& multiBlock);
    plint getNx() const;
    plint getNy() const;
    plint getNumBlocks() const;
    Box2D getSmallestBlock() const;
    Box2D getLargestBlock() const;
    plint getNumAllocatedCells() const;
};

template<typename T>
class MultiNTensorField3D {
public:
    Box3D getBoundingBox() const;
    plint getNdim() const;
    void initialize();
private:
    MultiNTensorField3D();
};

template<typename T>
class MultiNTensorInfo3D {
public:
    MultiNTensorInfo3D(MultiNTensorField3D<T> const& multiBlock);
    plint getNx() const;
    plint getNy() const;
    plint getNz() const;
    plint getNumBlocks() const;
    Box3D getSmallestBlock() const;
    Box3D getLargestBlock() const;
    plint getNumAllocatedCells() const;
};

}  // namespace plb

%template(PlbMultiNTensorField2DDouble) plb::MultiNTensorField2D<double>;
%template(PlbMultiBlockInfo2DDouble) plb::MultiNTensorInfo2D<double>;

%template(PlbMultiNTensorField2DInt) plb::MultiNTensorField2D<int>;
%template(PlbMultiBlockInfo2DInt) plb::MultiNTensorInfo2D<int>;

%template(PlbMultiNTensorField2DFloat) plb::MultiNTensorField2D<float>;
%template(PlbMultiBlockInfo2DFloat) plb::MultiNTensorInfo2D<float>;



%template(PlbMultiNTensorField3DDouble) plb::MultiNTensorField3D<double>;
%template(PlbMultiBlockInfo3DDouble) plb::MultiNTensorInfo3D<double>;

%template(PlbMultiNTensorField3DInt) plb::MultiNTensorField3D<int>;
%template(PlbMultiBlockInfo3DInt) plb::MultiNTensorInfo3D<int>;

%template(PlbMultiNTensorField3DFloat) plb::MultiNTensorField3D<float>;
%template(PlbMultiBlockInfo3DFloat) plb::MultiNTensorInfo3D<float>;
