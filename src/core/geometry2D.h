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


#ifndef GEOMETRY_2D_H
#define GEOMETRY_2D_H

#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include "core/util.h"
#include "core/array.h"
#include <algorithm>
#include <iterator>
#include <vector>

namespace plb {

/// Coordinates of a 2D Box
struct Box2D {
    Box2D() : x0(), x1(), y0(), y1() { }
    Box2D(plint x0_, plint x1_, plint y0_, plint y1_)
        : x0(x0_), x1(x1_), y0(y0_), y1(y1_)
    { }
    /// Return same box, shifted by (deltaX,deltaY)
    Box2D shift(plint deltaX, plint deltaY) const {
        return Box2D(x0+deltaX, x1+deltaX, y0+deltaY, y1+deltaY);
    }
    /// Return same box, rescaled by a factor scaling
    Box2D multiply(plint scaling) const {
        return Box2D(scaling*x0, scaling*x1, scaling*y0, scaling*y1);
    }
    /// Return same box, rescaled by a factor 1/scaling
    Box2D divide(plint scaling) const {
        return Box2D(x0/scaling, x1/scaling, y0/scaling, y1/scaling);
    }
    /// Rescale by 1/scaling and make sure to fit into higher-level box by trimming excess.
    /** This function should be used with grid refinement if the 
     *  coordinates of a fine grid are rescaled to fit on a coarse grid.
     *  The lower bounds of the resulting coarse box are increased by one
     *  if their fine-grid original version was not divisible by scaling;
     *  and the upper bounds are decreased in the same way.
     *  This makes sure that the coarse grid result does not exceed the
     *  bounds of the fine grid.
     */
    Box2D divideAndFitSmaller(plint scaling) const {
        return Box2D (util::roundUp(x0,scaling)/scaling, util::roundDown(x1,scaling)/scaling,
                      util::roundUp(y0,scaling)/scaling, util::roundDown(y1,scaling)/scaling);
    }
    /// Rescale by 1/scaling and make sure to fit into higher-level box by snapping to larger size.
    /** This function should be used with grid refinement if the 
     *  coordinates of a fine grid are rescaled to fit on a coarse grid.
     *  The lower bounds of the resulting coarse box are decreased by one
     *  if their fine-grid original version was not divisible by scaling;
     *  and the upper bounds are increased in the same way.
     *  This makes sure that the coarse grid result contains the
     *  bounds of the fine grid.
     */
    Box2D divideAndFitLarger(plint scaling) const {
        return Box2D (util::roundDown(x0,scaling)/scaling, util::roundUp(x1,scaling)/scaling,
                      util::roundDown(y0,scaling)/scaling, util::roundUp(y1,scaling)/scaling);
    }
    /// Add a border of nCells cells to the box
    Box2D enlarge(plint nCells) const {
        return Box2D(x0-nCells, x1+nCells, y0-nCells, y1+nCells);
    }
    /// Number of cells in x-direction
    plint getNx()  const { return (x1-x0+1); }
    /// Number of cells in y-direction
    plint getNy()  const { return (y1-y0+1); }
    /// Total number of cells in the box
    plint nCells() const { return getNx()*getNy(); }
    /// Return the maximum of getNx() and getNy()
    plint getMaxWidth() const { return std::max(getNx(), getNy()); }
    /// Copy the data into a 4-element array.
    Array<plint,4> to_plbArray() const {
        Array<plint,4> array;
        array[0] = x0; array[1] = x1;
        array[2] = y0; array[3] = y1;
        return array;
    }
    /// Initialize the data from a 4-element array.
    void from_plbArray(Array<plint,4> const& array) {
        x0 = array[0]; x1 = array[1];
        y0 = array[2]; y1 = array[3];
    }

    bool operator==(Box2D const& rhs) const {
        return x0 == rhs.x0 && y0 == rhs.y0 &&
               x1 == rhs.x1 && y1 == rhs.y1;
    }

    plint x0, x1, y0, y1;
};

/// Coordinates of a 2D point
struct Dot2D {
    Dot2D() : x(), y() { }
    Dot2D(plint x_, plint y_) : x(x_), y(y_) { }
    Dot2D& operator+=(Dot2D const& rhs) {
        x += rhs.x;
        y += rhs.y;
        return *this;
    }
    Dot2D& operator-=(Dot2D const& rhs) {
        x -= rhs.x;
        y -= rhs.y;
        return *this;
    }
    
    template<typename T>
    Array<T,2> to_plbArray() const { return Array<T,2>((T)x,(T)y); }
    
    plint x, y;
};

inline Dot2D operator+(Dot2D const& dot1, Dot2D const& dot2) {
    return Dot2D(dot1) += dot2;
}

inline Dot2D operator-(Dot2D const& dot1, Dot2D const& dot2) {
    return Dot2D(dot1) -= dot2;
}

/// Comparison operator based on Palabos' x-major ordering.
inline bool operator<(Dot2D const& dot1, Dot2D const& dot2) {
    return (dot1.x < dot2.x) || (
               (dot1.x == dot2.x) && (dot1.y < dot2.y) );
}

inline bool operator==(Dot2D const& dot1, Dot2D const& dot2) {
    return (dot1.x == dot2.x) && (dot1.y == dot2.y);
}

/// The global ordering on boxes has no specific interpretation, but
///   is useful to hold boxes in ordered data structures like sets.
inline bool operator<(Box2D const& box1, Box2D const& box2) {
    return ( Dot2D(box1.x0,box1.y0) < Dot2D(box2.x0,box2.y0) ) || (
                ( Dot2D(box1.x0,box1.y0) == Dot2D(box2.x0,box2.y0) ) &&
                ( Dot2D(box1.x1,box1.y1) <  Dot2D(box2.x1,box2.y1) ) );
}

/// List of 2D points, used to describe a subdomain
struct DotList2D {
    DotList2D() { }
    DotList2D(std::vector<Dot2D> const& dots_) : dots(dots_) { }
    /// Add one more point to the list
    void addDot(Dot2D dot) {
        dots.push_back(dot);
    }
    /// Add more points to the list
    void addDots(std::vector<Dot2D> const& dots_) {
        dots.insert(dots.end(), dots_.begin(), dots_.end());
    }
    /// Get const reference to one of the dots
    Dot2D const& getDot(plint whichDot) const {
        PLB_PRECONDITION( whichDot<getN() );
        return dots[whichDot];
    }
    /// Get non-const reference to one of the dots
    Dot2D& getDot(plint whichDot) {
        PLB_PRECONDITION( whichDot<getN() );
        return dots[whichDot];
    }
    /// Return same point list, shifted by (deltaX,deltaY)
    DotList2D shift(plint deltaX, plint deltaY) const {
        std::vector<Dot2D> newDots;
        std::vector<Dot2D>::const_iterator it = dots.begin();
        for (; it!=dots.end(); ++it) {
            newDots.push_back( Dot2D(it->x+deltaX, it->y+deltaY) );
        }
        return DotList2D(newDots);
    }
    /// Return same box, rescaled by a factor scaling
    DotList2D multiply(plint scaling) const {
        std::vector<Dot2D> newDots;
        std::vector<Dot2D>::const_iterator it = dots.begin();
        for (; it!=dots.end(); ++it) {
            newDots.push_back( Dot2D(it->x*scaling, it->y*scaling) );
        }
        return DotList2D(newDots);

    }
    /// Return same box, rescaled by a factor 1/scaling
    DotList2D divide(plint scaling) const {
        std::vector<Dot2D> newDots;
        std::vector<Dot2D>::const_iterator it = dots.begin();
        for (; it!=dots.end(); ++it) {
            newDots.push_back( Dot2D(it->x/scaling, it->y/scaling) );
        }
        return DotList2D(newDots);

    }
    /// Get total number of points
    plint getN() const {
        return dots.size();
    }

    std::vector<Dot2D> dots;
};

/// Decide if lattice point is contained in 2D box, boundaries inclusive
inline bool contained(plint x, plint y, Box2D const& box) {
    return x>=box.x0 && x<=box.x1 &&
           y>=box.y0 && y<=box.y1;
}

/// Decide if a Lagrangian point is contained in 2D box, boundaries exclusive
template<typename T>
inline bool contained(Array<T,2> const& x, Box2D const& box) {
    return x[0]>box.x0 && x[0]<box.x1 &&
           x[1]>box.y0 && x[1]<box.y1;
}

template<>
inline bool contained<plint>(Array<plint,2> const& x, Box2D const& box) {
    //IMPORTANT: the behavior of this function (for int) has changed in Palabos!
    //use contained(plint, plint, Box2D) instead.
    //TODO: In the future, this function will do the following instead of an assert:
    //return contained(x[0], x[1], box);
    PLB_ASSERT( false );
    return false;
}

/// Decide if a Lagrangian point is contained in 2D box, boundaries exclusive at the epsilon level of accuracy.
template<typename T>
inline bool contained(Array<T,2> const& x, Box2D const& box, T epsilon) {
    return x[0]-epsilon>box.x0 && x[0]+epsilon<box.x1 &&
           x[1]-epsilon>box.y0 && x[1]+epsilon<box.y1;
}

/// Decide if a Lagrangian point is contained in 2D box, boundaries inclusive at the epsilon level of accuracy.
template<typename T>
inline bool containedInclusive(Array<T,2> const& x, Box2D const& box, T epsilon)
{
    return util::greaterEqual_abs(x[0], (T) box.x0, epsilon) && util::lessEqual_abs(x[0], (T) box.x1, epsilon) &&
           util::greaterEqual_abs(x[1], (T) box.y0, epsilon) && util::lessEqual_abs(x[1], (T) box.y1, epsilon);
}

/// Decide if a Lagrangian point is contained in 2D box, boundaries inclusive at the floating point precision level of accuracy.
template<typename T>
inline bool containedInclusive(Array<T,2> const& x, Box2D const& box)
{
    return util::greaterEqual_abs(x[0], (T) box.x0) && util::lessEqual_abs(x[0], (T) box.x1) &&
           util::greaterEqual_abs(x[1], (T) box.y0) && util::lessEqual_abs(x[1], (T) box.y1);
}

/// Decide if lattice point is contained in 2D box, boundaries inclusive
inline bool contained(Dot2D dot, Box2D const& box) {
    return contained(dot.x, dot.y, box);
}

/// Decide if 2D box1 is contained in 2D box2, boundaries inclusive
inline bool contained(Box2D const& box1, Box2D const& box2) {
    return box1.x0>=box2.x0 && box1.x1<=box2.x1 &&
           box1.y0>=box2.y0 && box1.y1<=box2.y1;
}

/// Compute intersection between two 2D boxes
/** \return false if boxes don't intersect
 */
inline bool intersect(Box2D const& box1, Box2D const& box2, Box2D& inters) {
    inters.x0 = std::max(box1.x0,box2.x0);
    inters.y0 = std::max(box1.y0,box2.y0);

    inters.x1 = std::min(box1.x1,box2.x1);
    inters.y1 = std::min(box1.y1,box2.y1);

    return inters.x1>=inters.x0 && inters.y1>=inters.y0;
}

/// Determine if the intersection between box1 and box2 is non-empty
inline bool doesIntersect(Box2D const& box1, Box2D const& box2) {
    return (std::min(box1.x1,box2.x1) >= std::max(box1.x0,box2.x0)) &&
           (std::min(box1.y1,box2.y1) >= std::max(box1.y0,box2.y0));
}

/// If the two boxes can be merged into one, do it, and write the result
///   into the first box. Return value states if merging was successful.
inline bool merge(Box2D& box1, Box2D const& box2) {
    if (box1.x0==box2.x0 && box1.x1==box2.x1) {
        if (box2.y0==box1.y1+1) {
            box1.y1=box2.y1;
            return true;
        }
        if (box2.y1+1==box1.y0) {
            box1.y0=box2.y0;
            return true;
        }
    }
    if (box1.y0==box2.y0 && box1.y1==box2.y1) {
        if (box2.x0==box1.x1+1) {
            box1.x1=box2.x1;
            return true;
        }
        if (box2.x1+1==box1.x0) {
            box1.x0=box2.x0;
            return true;
        }
    }
    return false;
}


/// Compute intersection between a 2D box and a 2D DotList
/** \return false if the two don't intersect
 */
inline bool intersect(Box2D const& box, DotList2D const& dotlist, DotList2D& inters) {
    inters = DotList2D();
    std::vector<Dot2D>::const_iterator it = dotlist.dots.begin();
    for (; it != dotlist.dots.end(); ++it) {
        if ( contained(it->x,it->y, box) ) {
            inters.addDot(*it);
        }
    }
    return !inters.dots.empty();
}

/// Except the domain of box "toExcept" form the domain of box "originalBox".
/** The result consists of at most three boxes, which are added to the vector "result"
 */
inline void except(Box2D const& originalBox, Box2D const& toExcept, std::vector<Box2D>& result)
{
    Box2D intersection;
    bool doesIntersect = intersect(originalBox, toExcept, intersection);
    if (!doesIntersect) {
        result.push_back(originalBox);
    }
    else {
        if (intersection.x0 != originalBox.x0) {
            result.push_back(Box2D(originalBox.x0, intersection.x0-1, originalBox.y0, originalBox.y1));
        }
        if (intersection.x1 != originalBox.x1) {
            result.push_back(Box2D(intersection.x1+1, originalBox.x1, originalBox.y0, originalBox.y1));
        }
        if (intersection.y0 != originalBox.y0) {
            result.push_back(Box2D(intersection.x0, intersection.x1, originalBox.y0, intersection.y0-1));
        }
        if (intersection.y1 != originalBox.y1) {
            result.push_back(Box2D(intersection.x0, intersection.x1, intersection.y1+1, originalBox.y1));
        }
    }
}

/// Return a minimal box which contains both box1 and box2, boundaries inclusive.
inline Box2D bound(Box2D const& box1, Box2D const& box2) {
    return Box2D (
            std::min(box1.x0, box2.x0), std::max(box1.x1, box2.x1),
            std::min(box1.y0, box2.y0), std::max(box1.y1, box2.y1) );
}

/// Adjust two domains so they are of equal size (but not necessarily overlapping).
inline void adjustEqualSize(Box2D& fromDomain, Box2D& toDomain)
{
    plint deltaX = fromDomain.x0 - toDomain.x0;
    plint deltaY = fromDomain.y0 - toDomain.y0;
    Box2D intersection;
    intersect(fromDomain, toDomain.shift(deltaX,deltaY), intersection);
    fromDomain = intersection;
    toDomain = intersection.shift(-deltaX,-deltaY);
}

} // namespace plb

#endif  // GEOMETRY_2D_H
