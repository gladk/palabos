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


#ifndef GEOMETRY_3D_H
#define GEOMETRY_3D_H

#include "core/globalDefs.h"
#include "core/plbDebug.h"
#include "core/util.h"
#include "core/array.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"
#include <algorithm>
#include <iterator>
#include <vector>
#include <cmath>

namespace plb {

/// Coordinates of a 3D Box
struct Box3D {
    Box3D() : x0(), x1(), y0(), y1(), z0(), z1() { }
    Box3D(plint x0_, plint x1_, plint y0_, plint y1_, plint z0_, plint z1_)
        : x0(x0_), x1(x1_), y0(y0_), y1(y1_), z0(z0_), z1(z1_)
    { }
    /// Return same box, shifted by (deltaX,deltaY,deltaZ)
    Box3D shift(plint deltaX, plint deltaY, plint deltaZ) const {
        return Box3D(x0+deltaX, x1+deltaX, y0+deltaY, y1+deltaY, z0+deltaZ, z1+deltaZ);
    }
    /// Return same box, rescaled by a factor scaling
    Box3D multiply(plint scaling) const {
        return Box3D(scaling*x0, scaling*x1, scaling*y0,
                     scaling*y1, scaling*z0, scaling*z1);
    }
    /// Return same box, rescaled by a factor 1/scaling
    Box3D divide(plint scaling) const {
        return Box3D(x0/scaling, x1/scaling, y0/scaling,
                     y1/scaling, z0/scaling, z1/scaling);
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
    Box3D divideAndFitSmaller(plint scaling) const {
        return Box3D (util::roundUp(x0,scaling)/scaling, util::roundDown(x1,scaling)/scaling,
                      util::roundUp(y0,scaling)/scaling, util::roundDown(y1,scaling)/scaling,
                      util::roundUp(z0,scaling)/scaling, util::roundDown(z1,scaling)/scaling);
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
    Box3D divideAndFitLarger(plint scaling) const {
        return Box3D (util::roundDown(x0,scaling)/scaling, util::roundUp(x1,scaling)/scaling,
                      util::roundDown(y0,scaling)/scaling, util::roundUp(y1,scaling)/scaling,
                      util::roundDown(z0,scaling)/scaling, util::roundUp(z1,scaling)/scaling);
    }

    /// Add a border of nCells cells to the box
    Box3D enlarge(plint nCells) const {
        return Box3D(x0-nCells, x1+nCells, y0-nCells, y1+nCells, z0-nCells, z1+nCells);
    }
    /// Number of cells in x-direction
    plint getNx()  const { return (x1-x0+1); }
    /// Number of cells in y-direction
    plint getNy()  const { return (y1-y0+1); }
    /// Number of cells in z-direction
    plint getNz()  const { return (z1-z0+1); }
    /// Total number of cells in the box
    plint nCells() const { return getNx()*getNy()*getNz(); }
    /// Return the maximum of getNx(), getNy(), and getNz()
    plint getMaxWidth() const { return std::max(std::max(getNx(), getNy()), getNz()); }
    /// Copy the data into a 6-element array.
    Array<plint,6> to_plbArray() const {
        Array<plint,6> array;
        array[0] = x0; array[1] = x1;
        array[2] = y0; array[3] = y1;
        array[4] = z0; array[5] = z1;
        return array;
    }
    /// Initialize the data from a 6-element array.
    void from_plbArray(Array<plint,6> const& array) {
        x0 = array[0]; x1 = array[1];
        y0 = array[2]; y1 = array[3];
        z0 = array[4]; z1 = array[5];
    }

    bool operator==(Box3D const& rhs) {
        return x0 == rhs.x0 && y0 == rhs.y0 && z0 == rhs.z0 &&
               x1 == rhs.x1 && y1 == rhs.y1 && z1 == rhs.z1;
    }

    plint x0, x1, y0, y1, z0, z1;
};

/// Coordinates of a 3D point
struct Dot3D {
    Dot3D() : x(), y(), z() { }
    Dot3D(plint x_, plint y_, plint z_) : x(x_), y(y_), z(z_) { }
    Dot3D& operator+=(Dot3D const& rhs) {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }
    Dot3D& operator-=(Dot3D const& rhs) {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    }
    
    template<typename T>
    Array<T,3> to_plbArray() const { return Array<T,3>((T)x,(T)y,(T)z); }
    
    plint x, y, z;
};

inline Dot3D operator+(Dot3D const& dot1, Dot3D const& dot2) {
    return Dot3D(dot1) += dot2;
}

inline Dot3D operator-(Dot3D const& dot1, Dot3D const& dot2) {
    return Dot3D(dot1) -= dot2;
}

/// Comparison operator based on Palabos' x-major ordering.
inline bool operator<(Dot3D const& dot1, Dot3D const& dot2) {
    return (dot1.x < dot2.x) || (
               (dot1.x == dot2.x) && (
                   (dot1.y < dot2.y) || (
                       (dot1.y == dot2.y) && (dot1.z < dot2.z) ) ) );
}

inline bool operator==(Dot3D const& dot1, Dot3D const& dot2) {
    return (dot1.x == dot2.x) && (dot1.y == dot2.y) &&
           (dot1.z == dot2.z);
}

/// The global ordering on boxes has no specific interpretation, but
///   is useful to hold boxes in ordered data structures like sets.
inline bool operator<(Box3D const& box1, Box3D const& box2) {
    return ( Dot3D(box1.x0,box1.y0,box1.z0) < Dot3D(box2.x0,box2.y0,box2.z0) ) || (
                ( Dot3D(box1.x0,box1.y0,box1.z0) == Dot3D(box2.x0,box2.y0,box2.z0) ) &&
                ( Dot3D(box1.x1,box1.y1,box1.z1) <  Dot3D(box2.x1,box2.y1,box2.z1) ) );
}

/// List of 3D points, used to describe a subdomain
struct DotList3D {
    DotList3D() { }
    DotList3D(std::vector<Dot3D> const& dots_) : dots(dots_) { }
    /// Add one more point to the list
    void addDot(Dot3D dot) {
        dots.push_back(dot);
    }
    /// Add more points to the list
    void addDots(std::vector<Dot3D> const& dots_) {
        dots.insert(dots.end(), dots_.begin(), dots_.end());
    }
    /// Get const reference to one of the dots
    Dot3D const& getDot(plint whichDot) const {
        PLB_PRECONDITION( whichDot<getN() );
        return dots[whichDot];
    }
    /// Get non-const reference to one of the dots
    Dot3D& getDot(plint whichDot) {
        PLB_PRECONDITION( whichDot<getN() );
        return dots[whichDot];
    }
    /// Return same point list, shifted by (deltaX,deltaY,deltaZ)
    DotList3D shift(plint deltaX, plint deltaY, plint deltaZ) const {
        std::vector<Dot3D> newDots;
        std::vector<Dot3D>::const_iterator it = dots.begin();
        for (; it!=dots.end(); ++it) {
            newDots.push_back( Dot3D(it->x+deltaX, it->y+deltaY, it->z+deltaZ) );
        }
        return DotList3D(newDots);
    }
    /// Return same box, rescaled by a factor scaling
    DotList3D multiply(plint scaling) const {
        std::vector<Dot3D> newDots;
        std::vector<Dot3D>::const_iterator it = dots.begin();
        for (; it!=dots.end(); ++it) {
            newDots.push_back( Dot3D(it->x*scaling, it->y*scaling, it->z*scaling) );
        }
        return DotList3D(newDots);
    }
    /// Return same box, rescaled by a factor 1/scaling
    DotList3D divide(plint scaling) const {
        std::vector<Dot3D> newDots;
        std::vector<Dot3D>::const_iterator it = dots.begin();
        for (; it!=dots.end(); ++it) {
            newDots.push_back( Dot3D(it->x/scaling, it->y/scaling, it->z/scaling) );
        }
        return DotList3D(newDots);

    }
    /// Get total number of points
    plint getN() const {
        return dots.size();
    }

    std::vector<Dot3D> dots;
};

/// Decide if lattice point is contained in 3D box, boundaries inclusive
inline bool contained(plint x, plint y, plint z, Box3D const& box) {
    return x>=box.x0 && x<=box.x1 &&
           y>=box.y0 && y<=box.y1 &&
           z>=box.z0 && z<=box.z1;
}

/// Decide if a Lagrangian point is contained in 3D box, boundaries exclusive
template<typename T>
inline bool contained(Array<T,3> const& x, Box3D const& box) {
    return x[0]>box.x0 && x[0]<box.x1 &&
           x[1]>box.y0 && x[1]<box.y1 &&
           x[2]>box.z0 && x[2]<box.z1;
}

template<>
inline bool contained<plint>(Array<plint,3> const& x, Box3D const& box) {
    //IMPORTANT: the behavior of this function (for int) has changed in Palabos!
    //use contained(plint, plint, plint, Box3D) instead.
    //TODO: In the future, this function will do the following instead of an assert:
    //return contained(x[0], x[1], x[2], box);
    PLB_ASSERT( false );
    return false;
}

/// Decide if a Lagrangian point is contained in 3D box, boundaries exclusive at the epsilon level of accuracy.
template<typename T>
inline bool contained(Array<T,3> const& x, Box3D const& box, T epsilon) {
    return x[0]-epsilon>box.x0 && x[0]+epsilon<box.x1 &&
           x[1]-epsilon>box.y0 && x[1]+epsilon<box.y1 &&
           x[2]-epsilon>box.z0 && x[2]+epsilon<box.z1;
}

/// Decide if a Lagrangian point is contained in 3D box, boundaries inclusive at the epsilon level of accuracy.
template<typename T>
inline bool containedInclusive(Array<T,3> const& x, Box3D const& box, T epsilon)
{
    return util::greaterEqual_abs(x[0], (T) box.x0, epsilon) && util::lessEqual_abs(x[0], (T) box.x1, epsilon) &&
           util::greaterEqual_abs(x[1], (T) box.y0, epsilon) && util::lessEqual_abs(x[1], (T) box.y1, epsilon) &&
           util::greaterEqual_abs(x[2], (T) box.z0, epsilon) && util::lessEqual_abs(x[2], (T) box.z1, epsilon);
}

/// Decide if a Lagrangian point is contained in 3D box, boundaries inclusive at the floating point precision level of accuracy.
template<typename T>
inline bool containedInclusive(Array<T,3> const& x, Box3D const& box)
{
    return util::greaterEqual_abs(x[0], (T) box.x0) && util::lessEqual_abs(x[0], (T) box.x1) &&
           util::greaterEqual_abs(x[1], (T) box.y0) && util::lessEqual_abs(x[1], (T) box.y1) &&
           util::greaterEqual_abs(x[2], (T) box.z0) && util::lessEqual_abs(x[2], (T) box.z1);
}

/// Decide if lattice point is contained in 3D box, boundaries inclusive
inline bool contained(Dot3D dot, Box3D const& box) {
    return contained(dot.x, dot.y, dot.z, box);
}

/// Decide if 3D box1 is contained in 3D box2, boundaries inclusive
inline bool contained(Box3D const& box1, Box3D const& box2) {
    return box1.x0>=box2.x0 && box1.x1<=box2.x1 &&
           box1.y0>=box2.y0 && box1.y1<=box2.y1 &&
           box1.z0>=box2.z0 && box1.z1<=box2.z1;
}

/// Compute intersection between two 3D boxes
/** \return false if boxes don't intersect
 */
inline bool intersect(Box3D const& box1, Box3D const& box2, Box3D& inters) {
    inters.x0 = std::max(box1.x0,box2.x0);
    inters.y0 = std::max(box1.y0,box2.y0);
    inters.z0 = std::max(box1.z0,box2.z0);

    inters.x1 = std::min(box1.x1,box2.x1);
    inters.y1 = std::min(box1.y1,box2.y1);
    inters.z1 = std::min(box1.z1,box2.z1);

    return inters.x1>=inters.x0 && inters.y1>=inters.y0 && inters.z1>=inters.z0;
}

/// Determine if the intersection between box1 and box2 is non-empty
inline bool doesIntersect(Box3D const& box1, Box3D const& box2) {
    return (std::min(box1.x1,box2.x1) >= std::max(box1.x0,box2.x0)) &&
           (std::min(box1.y1,box2.y1) >= std::max(box1.y0,box2.y0)) &&
           (std::min(box1.z1,box2.z1) >= std::max(box1.z0,box2.z0));
}

/// If the two boxes can be merged into one, do it, and write the result
///   into the first box. Return value states if merging was successful.
inline bool merge(Box3D& box1, Box3D const& box2) {
    if (box1.x0==box2.x0 && box1.x1==box2.x1 &&
        box1.y0==box2.y0 && box1.y1==box2.y1) {
        if (box2.z0==box1.z1+1) {
            box1.z1=box2.z1;
            return true;
        }
        if (box2.z1+1==box1.z0) {
            box1.z0=box2.z0;
            return true;
        }
    }
    if (box1.x0==box2.x0 && box1.x1==box2.x1 &&
        box1.z0==box2.z0 && box1.z1==box2.z1) {
        if (box2.y0==box1.y1+1) {
            box1.y1=box2.y1;
            return true;
        }
        if (box2.y1+1==box1.y0) {
            box1.y0=box2.y0;
            return true;
        }
    }
    if (box1.y0==box2.y0 && box1.y1==box2.y1 &&
        box1.z0==box2.z0 && box1.z1==box2.z1) {
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

/// Compute intersection between a 3D box and a 3D DotList
/** \return false if the two don't intersect
 */
inline bool intersect(Box3D const& box, DotList3D const& dotlist, DotList3D& inters) {
    inters = DotList3D();
    std::vector<Dot3D>::const_iterator it = dotlist.dots.begin();
    for (; it != dotlist.dots.end(); ++it) {
        if ( contained(it->x,it->y,it->z, box) ) {
            inters.addDot(*it);
        }
    }
    return !inters.dots.empty();
}

/// Except the domain of box "toExcept" form the domain of box "originalBox"
/** The result consists of three boxes, which are added to the vector "result"
 */
inline void except(Box3D const& originalBox, Box3D const& toExcept, std::vector<Box3D>& result)
{
    Box3D intersection;
    bool doesIntersect = intersect(originalBox, toExcept, intersection);
    if (!doesIntersect) {
        result.push_back(originalBox);
    }
    else {
        if (intersection.x0 != originalBox.x0) {
            result.push_back(Box3D(originalBox.x0, intersection.x0-1, originalBox.y0, originalBox.y1, originalBox.z0, originalBox.z1));
        }
        if (intersection.x1 != originalBox.x1) {
            result.push_back(Box3D(intersection.x1+1, originalBox.x1, originalBox.y0, originalBox.y1, originalBox.z0, originalBox.z1));
        }
        if (intersection.y0 != originalBox.y0) {
            result.push_back(Box3D(intersection.x0, intersection.x1, originalBox.y0, intersection.y0-1, originalBox.z0, originalBox.z1));
        }
        if (intersection.y1 != originalBox.y1) {
            result.push_back(Box3D(intersection.x0, intersection.x1, intersection.y1+1, originalBox.y1, originalBox.z0, originalBox.z1));
        }
        if (intersection.z0 != originalBox.z0) {
            result.push_back(Box3D(intersection.x0, intersection.x1, intersection.y0, intersection.y1, originalBox.z0, intersection.z0-1));
        }
        if (intersection.z1 != originalBox.z1) {
            result.push_back(Box3D(intersection.x0, intersection.x1, intersection.y0, intersection.y1, intersection.z1+1, originalBox.z1));
        }
    }
}


/// Return a minimal box which contains both box1 and box2, boundaries inclusive.
inline Box3D bound(Box3D const& box1, Box3D const& box2) {
    return Box3D (
            std::min(box1.x0, box2.x0), std::max(box1.x1, box2.x1),
            std::min(box1.y0, box2.y0), std::max(box1.y1, box2.y1),
            std::min(box1.z0, box2.z0), std::max(box1.z1, box2.z1) );
}

/// Adjust two domains so they are of equal size (but not necessarily overlapping).
inline void adjustEqualSize(Box3D& fromDomain, Box3D& toDomain)
{
    plint deltaX = fromDomain.x0 - toDomain.x0;
    plint deltaY = fromDomain.y0 - toDomain.y0;
    plint deltaZ = fromDomain.z0 - toDomain.z0;
    Box3D intersection;
    intersect(fromDomain, toDomain.shift(deltaX,deltaY,deltaZ), intersection);
    fromDomain = intersection;
    toDomain = intersection.shift(-deltaX,-deltaY,-deltaZ);
}

/// Definition of a 3D unbounded plane.
template<typename T>
struct Plane {
    Plane(Array<T,3> point_, Array<T,3> normal_)
        : point(point_),
          normal(normal_)
    { }
    Plane() {
        point.resetToZero();
        normal.resetToZero();
    }
    Array<T,3> point;
    Array<T,3> normal;
};

/// Definition of a cuboid (rectangular parallelepiped).
template<typename T>
struct Cuboid {
    Cuboid(Array<T,3> lowerLeftCorner_, Array<T,3> upperRightCorner_)
        : lowerLeftCorner(lowerLeftCorner_),
          upperRightCorner(upperRightCorner_)
    { }
    Cuboid() {
        lowerLeftCorner.resetToZero();
        upperRightCorner.resetToZero();
    }
    Array<T,3> lowerLeftCorner;
    Array<T,3> upperRightCorner;
};

/// Decide if a Lagrangian point is contained in a cuboid, boundaries inclusive
template<typename T>
inline bool contained(Array<T,3> const& x, Cuboid<T> const& c)
{
    Array<T,3> llc = c.lowerLeftCorner;
    Array<T,3> urc = c.upperRightCorner;

    return x[0]>=llc[0] && x[0]<=urc[0] &&
           x[1]>=llc[1] && x[1]<=urc[1] &&
           x[2]>=llc[2] && x[2]<=urc[2];
}

/// Function to compute the intersection between a plane "plane" which is defined by
///   a point and a unit normal, and a line segment between points "point1"
///   and "point2". "intersection" is the object whose state is changed by this function.
///   This state is undefined, and cannot be used by the caller function, when the
///   return value of this function is not 1. This function returns 1 if a unique
///   intersection is found, -1 if the line belongs to the plane or does not intersect
///   it, and 0 otherwise.
template<typename T>
inline int lineIntersectionWithPlane (
        Plane<T> const& plane, Array<T,3> const& point1,
        Array<T,3> const& point2, Precision precision, Array<T,3>& intersection )
{
    T epsilon = getEpsilon<T>(precision);

    Array<T,3> direction = point2 - point1;

    T num = dot(plane.point, plane.normal) - dot(point1, plane.normal);
    T denom = dot(direction, plane.normal);

    if (std::fabs(denom) <=  epsilon) {
        return -1; // Line belongs to the plane or does not intersect it.
    }

    T t = num / denom;

    if ((t < 0.0 && !(std::fabs(t) <= epsilon)) ||
        (t > 1.0 && !(std::fabs(t - 1.0) <= epsilon))) {
        return 0;
    }

    intersection = point1 + direction * t; // Intersection point with the plane.

    return 1;
}

/* This functions checks if the point "point" belongs to the half-space which is
 * defined by the plane "plane" and pointed-into by the plane's normal.
 */
template<typename T>
bool isInHalfSpace(Array<T,3> const& point, Plane<T> const& plane, T& signedDistanceFromPlane)
{
    Array<T,3> dx = point - plane.point;
    signedDistanceFromPlane = dot<T,3>(dx, plane.normal);
    if (signedDistanceFromPlane >= (T) 0) {
        return true;
    }
    return false;
}

} // namespace plb

#endif  // GEOMETRY_3D_H
