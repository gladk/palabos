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
 * Global Definitions -- header file.
 */
#ifndef GLOBAL_DEFS_H
#define GLOBAL_DEFS_H

#ifdef PLB_MPI_PARALLEL
#include "mpi.h"
#endif

#include "core/plbDebug.h"

#include <limits>
#include <string>
#include <cstddef>
#include <algorithm>

namespace plb {

// We take the definition that the base type of T is just T, except
// in special cases. For example, for complex numbers, the base type
// of Complex<U> is U. PlbTraits is correspondingly overloaded in
// plbComplex.h.
template<typename T>
struct PlbTraits {
    typedef T BaseType;
};

/// Integer type for Palabos
/** On some architectures, this type is larger
 *  than int. Using plint instead of int ensures
 *  64-bit compatibility of the code.
 */
#ifdef PLB_BGP
typedef long long int plint;
#else
typedef ptrdiff_t plint;
#endif

/// Unsigned integer type for Palabos
/** On some architectures, this type is larger
 *  than int. Using fluplint instead of unsigned plint 
 *  ensures 64-bit compatibility of the code.
 **/
#ifdef PLB_BGP
typedef unsigned long long int pluint;
#else
typedef size_t pluint;
#endif

/// Define the endianness
/** The endianess is little-endian by default.
 *  On the Blue Gene/P, the endianness is
 *  big-endian.
 */
#ifdef PLB_BGP
#define PLB_BIG_ENDIAN
#endif

/// Unsigned integer type for tracking global ids of individual objects.
typedef unsigned long id_t;

/// Enumeration type that sets the precision for further use.
/// Single precision (float): FLT
/// Double precision (double): DBL
/// Extended precision (long double): LDBL
/// "Infinite" precision: INF
enum Precision { FLT, DBL, LDBL, INF };

template<typename T>
inline Precision floatingPointPrecision()
{
    if (sizeof(T) == sizeof(float)) {
        return FLT;
    } else if (sizeof(T) == sizeof(long double)) {
        return LDBL;
    }
    return DBL;
}

template<typename T>
inline T getEpsilon(Precision precision, T scale = T(1))
{
    PLB_ASSERT(scale > (T) 0);

    T epsilon;
    switch (precision) {
    case FLT:
        epsilon = std::numeric_limits<float>::epsilon();
        break;
    case DBL:
        epsilon = std::numeric_limits<double>::epsilon();
        break;
    case LDBL:
        epsilon = std::numeric_limits<long double>::epsilon();
        break;
    case INF: default:
        epsilon = (T) 0;
        break;
    }

    return scale * epsilon;
}

// Version that works also for integral types, and always refers to the
// type defined by the template argument. This means that if T is a floating point
// type, then "getEpsilon<T>(scale)" is the same as "getEplsilon<T>(floatingPointPrecision<T>(), scale)".
// The function "getEpsilon<T>(Precision, scale)" works for floating point types only, and can
// be used with two different precisions (e.g. T is double, but Precision is float). This
// version works with one type (T) and T can be also integral (in such a case, epsilon is zero).
template<typename T>
inline T getEpsilon(T scale = T(1))
{
    PLB_ASSERT(scale > (T) 0);
    static T epsilon = std::numeric_limits<T>::epsilon();
    return scale * epsilon;
}

/// Enumeration type that sets the file format for the triangular surface meshes.
/// Stereolithography format: STL
/// Object file format: OFF
enum SurfaceGeometryFileFormat { STL, OFF };

/// Ordering of indices when a BlockXD is converted into a serial data stream.
/** Signification of constants:
 *    - forward:  Right-most index (y in 2D and z in 3D) is contiguous in memory.
 *                For non-allocated parts of the Block, the value 0 is produced for
 *                output, and values are ignored during input.
 *    - backward: Left-most index (x) is contiguous in memory.
 *                For non-allocated parts of the Block, the value 0 is produced for
 *                output, and values are ignored during input.
 *    - memorySaving: Ordering is forward (this respects the natural ordering in
 *                    Palabos). Non-allocated parts of the Block are neither written
 *                    or read: memory savings in the program are reflected by memory
 *                    savings on the disk.
 **/
namespace IndexOrdering {
    enum OrderingT {forward, backward, memorySaving};
}

/// Sub-domain of an atomic-block, on which for example a data processor is executed.
/** Signification of constants:
 *      - bulk: Refers to bulk-nodes, without envelope.
 *      - envelope: Refers to nodes on the envelope, without the bulk.
 *      - bulkAndEnvelope: Refers to the full domain.
 **/
namespace BlockDomain {

    enum DomainT {bulk, envelope, bulkAndEnvelope};
    inline bool usesEnvelope(DomainT domain) {
        return domain==envelope || domain==bulkAndEnvelope;
    }

}

namespace modif {

    /// Indicates what kind of cell content was modified and must
    ///   be updated in a multi-block structure.
    enum ModifT {
        nothing          =0,  //< No modification.
        staticVariables  =1,  //< Static cell content (populations+externals).
        dynamicVariables =2,  //< Only content of dynamics objects, but no static content.
        allVariables     =3,  //< Both the static and dynamic cell content.
        dataStructure    =4,  //< Recreate dynamics and copy both static and dynamic content.
        undefined        =5   
    };

    enum { numConstants=5 };

    /// If two data processors act on a block, combine their ModifT modification to
    ///   determine a resulting (worst-case) ModifT modification.
    inline ModifT combine(ModifT type1, ModifT type2) {
        ModifT result = std::max(type1, type2);
        // Note: static+dynamic = all.
        if ( result==dynamicVariables &&
             (type1==staticVariables || type2==staticVariables) )
        {
            result = allVariables;
        }
        return result;
    }
}

namespace global {

class IOpolicyClass {
public:
    void setIndexOrderingForStreams(IndexOrdering::OrderingT streamOrdering_);
    IndexOrdering::OrderingT getIndexOrderingForStreams() const;

    void setEndianSwitchOnBase64out(bool doSwitch);
    bool getEndianSwitchOnBase64out();

    void setEndianSwitchOnBase64in(bool doSwitch);
    bool getEndianSwitchOnBase64in();

    void setStlFilesHaveLowerBound(bool flag);
    bool stlFilesHaveLowerBound() const;

    void setLowerBoundForStlFiles(double stlLowerBound_);
    double getLowerBoundForStlFiles() const;

    void activateParallelIO(bool activate);
    bool useParallelIO() const;
private:
    IOpolicyClass();
private:
    IndexOrdering::OrderingT streamOrdering;
    bool endianSwitchOnBase64out;
    bool endianSwitchOnBase64in;
    bool stlLowerBoundFlag;
    double stlLowerBound;
    bool parallelIOflag;
    friend IOpolicyClass& IOpolicy();
};
    
class Directories {
public:
    /// Collectively set all output directories to the value specified by outputDir.
    void setOutputDir(std::string outputDir_);
    /// Specify output directory for logfile.
    void setLogOutDir(std::string logOutDir_);
    /// Specify output directory for the ImageWriter.
    void setImageOutDir(std::string imageOutDir_);
    /// Specify output directory for the VtkDataWriter.
    void setVtkOutDir(std::string inputDir_);
    /// Specify directory for input files.
    void setInputDir(std::string inputDir);

    /// Get output directory for logfile.
    std::string getLogOutDir() const;
    /// Get output directory for the ImageWriter.
    std::string getImageOutDir() const;
    /// Get output directory for the VtkDataWriter.
    std::string getVtkOutDir() const;
    /// Get directory for input files.
    std::string getInputDir() const;
    /// Get generic output directory.
    std::string getOutputDir() const;
private:
    Directories();
private:
    std::string logOutDir;
    std::string imageOutDir;
    std::string vtkOutDir;
    std::string inputDir;
    std::string outputDir;
friend Directories& directories();
};

inline Directories& directories() {
    static Directories singleton;
    return singleton;
}

inline IOpolicyClass& IOpolicy() {
    static IOpolicyClass singleton;
    return singleton;
}

}  // namespace global

}  // namespace plb

#endif
