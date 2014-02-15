namespace plb {
    template<typename T> class MultiNTensorField3D;
}
/*
 * core/geometry3D.i
 */

%{
#include "PALABOS_ROOT/src/core/geometry3D.h"
%}

namespace plb {

struct Box3D {
    Box3D();
    Box3D(plint x0_, plint x1_, plint y0_, plint y1_, plint z0_, plint z1_);
    Box3D shift(plint deltaX, plint deltaY, plint deltaZ) const;
    plint getNx() const;
    plint getNy() const;
    plint getNz() const;
    plint nCells() const;

    plint x0, x1, y0, y1, z0, z1;
};

bool intersect(Box3D const& box1, Box3D const& box2, Box3D& inters);

}  // namespace plb

