namespace plb {
    template<typename T> class MultiNTensorField3D;
}
/*
 * core/geometry2D.i
 */

%{
#include "PALABOS_ROOT/src/core/geometry2D.h"
%}

namespace plb {

struct Box2D {
    Box2D();
    Box2D(plint x0_, plint x1_, plint y0_, plint y1_);
    Box2D shift(plint deltaX, plint deltaY) const;
    plint getNx() const;
    plint getNy() const;
    plint nCells() const;

    plint x0, x1, y0, y1;
};

bool intersect(Box2D const& box1, Box2D const& box2, Box2D& inters);

}  // namespace plb

