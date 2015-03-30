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

#ifndef FREE_SURFACE_MODEL_3D_H
#define FREE_SURFACE_MODEL_3D_H

#include <algorithm>
#include "core/globalDefs.h"
#include "atomicBlock/dataProcessingFunctional3D.h"
#include "multiBlock/defaultMultiBlockPolicy3D.h"
#include "multiPhysics/freeSurfaceUtil3D.h"
#include "multiPhysics/freeSurfaceInitializer3D.h"
#include "dataProcessors/dataInitializerWrapper3D.h"
#include "basicDynamics/dynamicsProcessor3D.h"

namespace plb {

template<typename T, template<typename U> class Descriptor>
class TwoPhaseComputeNormals3D : public BoxProcessingFunctional3D {
public:
    TwoPhaseComputeNormals3D()
    {
        precision = floatingPointPrecision<T>();
    }
    virtual TwoPhaseComputeNormals3D<T,Descriptor>* clone() const {
        return new TwoPhaseComputeNormals3D<T,Descriptor>(*this);
    }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks);
    virtual void getTypeOfModification (std::vector<modif::ModifT>& modified) const {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;         // Fluid.
        modified[1] = modif::nothing;         // rhoBar.
        modified[2] = modif::nothing;         // j.
        modified[3] = modif::nothing;         // Mass.
        modified[4] = modif::nothing;         // Volume fraction.
        modified[5] = modif::nothing;         // Flag-status.
        modified[6] = modif::staticVariables; // Normal.
        modified[7] = modif::nothing;         // Interface-lists.
        modified[8] = modif::nothing;         // Curvature.
        modified[9] = modif::nothing;         // Outside density.
    }
private:
    Precision precision;
};

template<typename T, template<typename U> class Descriptor>
class FreeSurfaceGeometry3D : public BoxProcessingFunctional3D {
public:
    FreeSurfaceGeometry3D(T contactAngle_)
        : contactAngle(contactAngle_)
    {
        Precision precision = floatingPointPrecision<T>();
        eps = getEpsilon<T>(precision);

        // The contact angle must take values between 0 and 180 degrees. If it is negative,
        // this means that contact angle effects will not be modeled.
        PLB_ASSERT(contactAngle < (T) 180.0 || std::fabs(contactAngle - (T) 180.0) <= eps);

        if (contactAngle < (T) 0.0 && std::fabs(contactAngle) > eps) {
            useContactAngle = 0;
        } else {
            useContactAngle = 1;
        }
    }
    virtual FreeSurfaceGeometry3D<T,Descriptor>* clone() const {
        return new FreeSurfaceGeometry3D<T,Descriptor>(*this);
    }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks);
    virtual void getTypeOfModification (std::vector<modif::ModifT>& modified) const {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;             // Fluid.
        modified[1] = modif::nothing;             // rhoBar.
        modified[2] = modif::nothing;             // j.
        modified[3] = modif::nothing;             // Mass.
        modified[4] = modif::nothing;             // Volume fraction.
        modified[5] = modif::nothing;             // Flag-status.
        modified[6] = modif::staticVariables;     // Normal.
        modified[7] = modif::nothing;             // Interface-lists.
        modified[8] = modif::staticVariables;     // Curvature.
        modified[9] = modif::nothing;             // Outside density.
    }
private:
    ScalarField3D<int> *getInterfaceFlags(Box3D domain, FreeSurfaceProcessorParam3D<T,Descriptor>& param);
    void computeHeights3D(FreeSurfaceProcessorParam3D<T,Descriptor>& param, int integrationDirection,
            plint iX, plint iY, plint iZ, T h[3][3]);
    void computeHeights2D(FreeSurfaceProcessorParam3D<T,Descriptor>& param, Array<int,3>& wallTangent0,
            Array<int,3>& wallTangent1, int integrationDirection, plint iX, plint iY, plint iZ, T h[3]);
private:
    enum { unTagged = 0, notInterface = 1, regular = 2, contactLine = 4, adjacent = 10 };
private:
    T contactAngle;
    int useContactAngle;
    T eps;
};

template<typename T, template<typename U> class Descriptor>
class TwoPhaseComputeCurvature3D : public BoxProcessingFunctional3D {
public:
    TwoPhaseComputeCurvature3D(T contactAngle_, Box3D globalBoundingBox_)
        : contactAngle(contactAngle_),
          globalBoundingBox(globalBoundingBox_)
    {
        precision = floatingPointPrecision<T>();
        T eps = getEpsilon<T>(precision);

        // The contact angle must take values between 0 and 180 degrees. If it is negative,
        // this means that contact angle effects will not be modeled.
        PLB_ASSERT(contactAngle < (T) 180.0 || std::fabs(contactAngle - (T) 180.0) <= eps);

        if (contactAngle < (T) 0.0 && std::fabs(contactAngle) > eps) {
            useContactAngle = 0;
        } else {
            useContactAngle = 1;
        }

        if (useContactAngle) {
            T pi = 3.14159265358979323844;
            contactAngle *= pi / (T) 180.0;
        }
    }
    virtual TwoPhaseComputeCurvature3D<T,Descriptor>* clone() const {
        return new TwoPhaseComputeCurvature3D<T,Descriptor>(*this);
    }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks);
    virtual void getTypeOfModification (std::vector<modif::ModifT>& modified) const {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;             // Fluid.
        modified[1] = modif::nothing;             // rhoBar.
        modified[2] = modif::nothing;             // j.
        modified[3] = modif::nothing;             // Mass.
        modified[4] = modif::nothing;             // Volume fraction.
        modified[5] = modif::nothing;             // Flag-status.
        modified[6] = modif::nothing;             // Normal.
        modified[7] = modif::nothing;             // Interface-lists.
        modified[8] = modif::staticVariables;     // Curvature.
        modified[9] = modif::nothing;             // Outside density.
    }
private:
    T contactAngle;
    int useContactAngle;
    Box3D globalBoundingBox;
    Precision precision;
};

/// Compute the mass balance on every node in the domain, and store in mass matrix.
/** Input:
  *   - Flag-status:   needed in bulk+1
  *   - Mass:          needed in bulk
  *   - Volume fraction: needed in bulk
  *   - Populations:   needed in bulk+1
  * Output:
  *   - mass.
  **/
template< typename T,template<typename U> class Descriptor>
class FreeSurfaceMassChange3D : public BoxProcessingFunctional3D {
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks);
    virtual FreeSurfaceMassChange3D<T,Descriptor>* clone() const {
        return new FreeSurfaceMassChange3D<T,Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;         // Fluid.
        modified[1] = modif::nothing;         // rhoBar.
        modified[2] = modif::nothing;         // j.
        modified[3] = modif::staticVariables; // Mass.
        modified[4] = modif::nothing;         // Volume fraction.
        modified[5] = modif::nothing;         // Flag-status.
        modified[6] = modif::nothing;         // Normal.
        modified[7] = modif::nothing;         // Interface-lists.
        modified[8] = modif::nothing;         // Curvature.
        modified[9] = modif::nothing;         // Outside density.
   }
};

/// Completion scheme on the post-collide populations on interface cells.
/** Input:
  *   - Flag-status:   needed in bulk+1
  *   - Volume fraction: needed in bulk+1
  *   - Populations:   needed in bulk+1
  *   - Momentum:      needed in bulk+1
  *   - Density:       needed in bulk+1
  * Output:
  *   - Populations.
  **/
// ASK: This data processor loops over the whole volume. Is this really
//      necessary, or could one of the lists be used instead?
template<typename T, template<typename U> class Descriptor>
class FreeSurfaceCompletion3D : public BoxProcessingFunctional3D {
public:
    virtual FreeSurfaceCompletion3D<T,Descriptor>* clone() const {
        return new FreeSurfaceCompletion3D<T,Descriptor>(*this);
    }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks);
    virtual void getTypeOfModification (std::vector<modif::ModifT>& modified) const {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;         // Fluid. Should be: staticVariables.
        modified[1] = modif::nothing;         // rhoBar.
        modified[2] = modif::nothing;         // j.
        modified[3] = modif::nothing;         // Mass.
        modified[4] = modif::nothing;         // Volume fraction.
        modified[5] = modif::nothing;         // Flag-status.
        modified[6] = modif::nothing;         // Normal.
        modified[7] = modif::nothing;         // Interface-lists.
        modified[8] = modif::nothing;         // Curvature.
        modified[9] = modif::nothing;         // Outside density.
    }
};

/// Compute and store mass-fraction and macroscopic variables.
/** Input:
  *   - Flag-status:   needed in bulk
  *   - Mass:          needed in bulk
  *   - Populations:   needed in bulk
  * Output:
  *   - mass-fraction, density, momentum, flag (because setting bounce-back).
  **/
template<typename T, template<typename U> class Descriptor>
class FreeSurfaceMacroscopic3D : public BoxProcessingFunctional3D {
public:
    FreeSurfaceMacroscopic3D(T rhoDefault_)
        : rhoDefault(rhoDefault_)
    { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks);
    virtual FreeSurfaceMacroscopic3D<T,Descriptor>* clone() const {
        return new FreeSurfaceMacroscopic3D<T,Descriptor>(*this);
    }
    virtual void getTypeOfModification (std::vector<modif::ModifT>& modified) const {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;         // Fluid. Should be: staticVariables.
        modified[1] = modif::staticVariables; // rhoBar.
        modified[2] = modif::staticVariables; // j.
        modified[3] = modif::staticVariables; // Mass. Should be: staticVariables.
        modified[4] = modif::staticVariables; // Volume fraction.
        modified[5] = modif::staticVariables; // Flag-status.
        modified[6] = modif::nothing;         // Normal.
        modified[7] = modif::nothing;         // Interface-lists.
        modified[8] = modif::nothing;         // Curvature.
        modified[9] = modif::nothing;         // Outside density.
    }
private:
    T rhoDefault;
};

/// Add the surface tension contribution.
/** Input:
  *   - Flag-status:   needed in bulk
  *   - Mass:          needed in bulk
  *   - Populations:   needed in bulk
  * Output:
  *   - mass-fraction, density, momentum, flag (because setting bounce-back).
  **/
template<typename T, template<typename U> class Descriptor>
class TwoPhaseAddSurfaceTension3D : public BoxProcessingFunctional3D {
public:
    TwoPhaseAddSurfaceTension3D(T surfaceTension_, T rhoDefault_)
        : surfaceTension(surfaceTension_),
          rhoDefault(rhoDefault_)
    { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks);
    virtual TwoPhaseAddSurfaceTension3D<T,Descriptor>* clone() const {
        return new TwoPhaseAddSurfaceTension3D<T,Descriptor>(*this);
    }
    virtual void getTypeOfModification (std::vector<modif::ModifT>& modified) const {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;         // Fluid. Should be: staticVariables.
        modified[1] = modif::staticVariables; // rhoBar.
        modified[2] = modif::staticVariables; // j.
        modified[3] = modif::staticVariables; // Mass.
        modified[4] = modif::staticVariables; // Volume fraction.
        modified[5] = modif::staticVariables; // Flag-status.
        modified[6] = modif::nothing;         // Normal.
        modified[7] = modif::nothing;         // Interface-lists.
        modified[8] = modif::nothing;         // Curvature.
        modified[9] = modif::nothing;         // Outside density.
    }
private:
    T surfaceTension;
    T rhoDefault;
};

/// Based on the current flag status, decide, upon the value of mass fraction, which nodes shall
///   switch state.
/** Input:
  *   - Volume fraction: needed in bulk+2
  *   - Flag-status:   needed in bulk+2
  * Output:
  *   - interface-to-fluid list: defined in bulk+2
  *   - interface-to-empty list: defined in bulk+1
  *   - empty-to-interface list: defined in bulk+1
  **/
template<typename T,template<typename U> class Descriptor>
class FreeSurfaceComputeInterfaceLists3D : public BoxProcessingFunctional3D
{
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks);
    virtual FreeSurfaceComputeInterfaceLists3D<T,Descriptor>* clone() const {
        return new FreeSurfaceComputeInterfaceLists3D<T,Descriptor>(*this);
    }
    virtual void getTypeOfModification (std::vector<modif::ModifT>& modified) const {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;         // Fluid (not used in this processor).
        modified[1] = modif::nothing;         // rhoBar.
        modified[2] = modif::nothing;         // j.
        modified[3] = modif::nothing;         // Mass (not used in this processor).
        modified[4] = modif::nothing;         // Volume fraction, read-only.
        modified[5] = modif::nothing;         // Flag-status, read-only.
        modified[6] = modif::nothing;         // Normal.
        modified[7] = modif::staticVariables; // Interface-lists; all lists are created here.
        modified[8] = modif::nothing;         // Curvature.
        modified[9] = modif::nothing;         // Outside density.
    }
private:
    static T kappa; // Safety threshold for state-change, to prevent back-and-forth oscillations.
};

/** Input:
  *   - interface-to-fluid list: needed in bulk+1
  *   - interface-to-empty list: needed in bulk+1
  *   - density: needed in bulk+1
  *   - mass:    needed in bulk+1
  *   - flag:    needed in bulk+1
  * Output:
  *   - flag, dynamics, mass, volumeFraction, density, force, momentum
  *   - mass-excess-list: defined in bulk+1
  **/
template<typename T,template<typename U> class Descriptor>
class FreeSurfaceIniInterfaceToAnyNodes3D : public BoxProcessingFunctional3D {
public:
    FreeSurfaceIniInterfaceToAnyNodes3D(T rhoDefault_);
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks);

    virtual FreeSurfaceIniInterfaceToAnyNodes3D<T,Descriptor>* clone() const {
        return new FreeSurfaceIniInterfaceToAnyNodes3D<T,Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;           // Fluid. Gets assigned new dynamics. Should be: dataStructure
        modified[1] = modif::staticVariables;   // rhoBar.
        modified[2] = modif::nothing;           // j. Should be: staticVariables.
        modified[3] = modif::staticVariables;   // Mass. Is redistributed and initialized from neighborying density.
        modified[4] = modif::nothing;           // Volume fraction. Is default-initialized. Should be: staticVariables.
        modified[5] = modif::staticVariables;   // Flag-status. Is adapted according to cell-change lists.
        modified[6] = modif::nothing;           // Normal.
        modified[7] = modif::nothing;           // Interface-lists. Read-only.
        modified[8] = modif::nothing;           // Curvature.
        modified[9] = modif::nothing;           // Outside density.
    }
private:
    T rhoDefault; 
};

/// Based on the previously computed empty->interface list, initialize flow variables for
///   new interface cells.
/** Input:
  *   - Populations: needed in bulk+0
  *   - Momentum:    needed in bulk+1
  *   - Density:     needed in bulk+1
  *   - Flag-status: needed in bulk+0
  * Output:
  *   - flag-status:   initialized to "interface" on corresponding cells.
  *   - lattice:       initialized from neighbor averages on new interface cells.
  *   - mass:          initialized to zero on new interface cells.
  *   - mass-fraction: initialized to zero on new interface cells.
  *   - momentum
  **/
template<typename T,template<typename U> class Descriptor>
class FreeSurfaceIniEmptyToInterfaceNodes3D: public BoxProcessingFunctional3D {
public:
    FreeSurfaceIniEmptyToInterfaceNodes3D(Dynamics<T,Descriptor>* dynamicsTemplate_, Array<T,Descriptor<T>::d> force_)
        : dynamicsTemplate(dynamicsTemplate_), force(force_)
    { }
    FreeSurfaceIniEmptyToInterfaceNodes3D(FreeSurfaceIniEmptyToInterfaceNodes3D<T,Descriptor> const& rhs)
        : dynamicsTemplate(rhs.dynamicsTemplate->clone()),
          force(rhs.force)
    { }
    FreeSurfaceIniEmptyToInterfaceNodes3D<T,Descriptor>* operator=(
            FreeSurfaceIniEmptyToInterfaceNodes3D<T,Descriptor> const& rhs)
    { 
        FreeSurfaceIniEmptyToInterfaceNodes3D<T,Descriptor>(rhs).swap(*this);
        return *this;
    }
    void swap(FreeSurfaceIniEmptyToInterfaceNodes3D<T,Descriptor>& rhs) {
        std::swap(dynamicsTemplate, rhs.dynamicsTemplate);
        std::swap(force, rhs.force);
    }
    ~FreeSurfaceIniEmptyToInterfaceNodes3D() {
        delete dynamicsTemplate;
    }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks);
    virtual FreeSurfaceIniEmptyToInterfaceNodes3D<T,Descriptor>* clone() const {
        return new FreeSurfaceIniEmptyToInterfaceNodes3D<T,Descriptor>(*this);
    }
    virtual void getTypeOfModification (std::vector<modif::ModifT>& modified) const {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;         // Fluid. Should be: dataStructure
        modified[1] = modif::staticVariables; // rhoBar.
        modified[2] = modif::nothing;         // j. Should be: staticVariables.
        modified[3] = modif::staticVariables; // Mass.
        modified[4] = modif::nothing;         // Volume fraction, read-only. Should be: staticVariables
        modified[5] = modif::staticVariables; // Flag-status, read-only.
        modified[6] = modif::nothing;         // Normal.
        modified[7] = modif::nothing;         // Interface-lists. Read access to gasCellToInitializeData.
        modified[8] = modif::nothing;         // Curvature.
        modified[9] = modif::nothing;         // Outside density.
    }
private:
    Dynamics<T,Descriptor>* dynamicsTemplate;
    Array<T,Descriptor<T>::d> force; // Body force, for initialization of the new interface cell.
};

/// Isolated cells cannot be part of the interface. This data processor spots and
/// removes them.
/** Input:
  *   - Flag-status: needed in bulk+2
  *   - mass:        needed in bulk+1
  *   - density:     needed in bulk+1
  * Output:
  *   - interfaceToFluidNodes:   initialized in bulk+1
  *   - interfaceToEmptyNodes:   initialized in bulk+1
  *   - massExcess list:         initialized in bulk+1
  *   - mass, density, mass-fraction, dynamics, force, momentum, flag: in bulk+1
  **/
template<typename T,template<typename U> class Descriptor>
class FreeSurfaceRemoveFalseInterfaceCells3D : public BoxProcessingFunctional3D {
public:
    FreeSurfaceRemoveFalseInterfaceCells3D(T rhoDefault_)
        : rhoDefault(rhoDefault_)
    { }
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks);
    virtual FreeSurfaceRemoveFalseInterfaceCells3D<T,Descriptor>* clone() const {
        return new FreeSurfaceRemoveFalseInterfaceCells3D<T,Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;          // Fluid: Gets NoDynamics when node changes to empty. Should be: dataStructure.
        modified[1] = modif::staticVariables;  // rhoBar.
        modified[2] = modif::nothing;          // j. Should be: staticVariables.
        modified[3] = modif::staticVariables;  // Mass.
        modified[4] = modif::nothing;          // Volume fraction. Should be: staticVariables.
        modified[5] = modif::staticVariables;  // Flag-status.
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface lists.
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
    }
private:
    T rhoDefault;  
};

/// Enforce exact mass balance when interface cells become fluid or empty.
/** Input:
  *   - mass-excess list: needed in bulk+1
  *   - Flag-status: needed in bulk+2
  *   - mass:        needed in bulk+2
  *   - density:     needed in bulk+2
  * Output:
  *   - mass, mass-fraction
  **/
template<typename T,template<typename U> class Descriptor>
class FreeSurfaceEqualMassExcessReDistribution3D : public BoxProcessingFunctional3D {
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks);
    virtual void redistribute( Box3D const& domain, Box3D const& originalDomain,
                               FreeSurfaceProcessorParam3D<T,Descriptor>& param,
                               plint iX, plint iY, plint iZ, T mass );
    virtual FreeSurfaceEqualMassExcessReDistribution3D<T,Descriptor>* clone() const {
        return new FreeSurfaceEqualMassExcessReDistribution3D(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::dataStructure;    // Fluid.
        modified[1] = modif::staticVariables;  // rhoBar.
        modified[2] = modif::staticVariables;  // j.
        modified[3] = modif::staticVariables;  // Mass.
        modified[4] = modif::staticVariables;  // Volume fraction.
        modified[5] = modif::nothing;          // Flag-status.
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface lists.
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
    }
};

/// Enforce exact mass balance when interface cells become fluid or empty, redistribute along interface normal.
/** Input:
  *   - mass-excess list: needed in bulk+1
  *   - interface normal: needed in bulk+1
  *   - Flag-status: needed in bulk+2
  *   - mass:        needed in bulk+2
  *   - density:     needed in bulk+2
  * Output:
  *   - mass, mass-fraction
  **/
template<typename T,template<typename U> class Descriptor>
class FreeSurfaceWeightedMassExcessReDistribution3D : public BoxProcessingFunctional3D {
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks);
    virtual void redistribute( Box3D const& domain, Box3D const& originalDomain,
                               FreeSurfaceProcessorParam3D<T,Descriptor>& param,
                               plint iX, plint iY, plint iZ, T mass, T sign );
    virtual FreeSurfaceWeightedMassExcessReDistribution3D<T,Descriptor>* clone() const {
        return new FreeSurfaceWeightedMassExcessReDistribution3D(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::dataStructure;    // Fluid.
        modified[1] = modif::staticVariables;  // rhoBar.
        modified[2] = modif::staticVariables;  // j.
        modified[3] = modif::staticVariables;  // Mass.
        modified[4] = modif::staticVariables;  // Volume fraction.
        modified[5] = modif::nothing;          // Flag-status.
        modified[6] = modif::nothing;          // Normal.
        modified[7] = modif::nothing;          // Interface lists.
        modified[8] = modif::nothing;          // Curvature.
        modified[9] = modif::nothing;          // Outside density.
    }
};

template<typename T,template<typename U> class Descriptor>
class TwoPhaseComputeStatistics3D : public BoxProcessingFunctional3D {
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks);
    virtual TwoPhaseComputeStatistics3D<T,Descriptor>* clone() const {
        return new TwoPhaseComputeStatistics3D(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;         // Fluid.
        modified[1] = modif::nothing;         // rhoBar.
        modified[2] = modif::nothing;         // j.
        modified[3] = modif::nothing;         // Mass.
        modified[4] = modif::nothing;         // Volume fraction.
        modified[5] = modif::nothing;         // Flag-status.
        modified[6] = modif::nothing;         // Normal.
        modified[7] = modif::nothing;         // Interface lists.
        modified[8] = modif::nothing;         // Curvature.
        modified[9] = modif::nothing;         // Outside density.
    }
};

template< typename T,template<typename U> class Descriptor>
class FreeSurfaceInterfaceFilter : public BoxProcessingFunctional3D {
public:
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks);
    virtual FreeSurfaceInterfaceFilter<T,Descriptor>* clone() const {
        return new FreeSurfaceInterfaceFilter<T,Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::nothing;         // Fluid.
        modified[1] = modif::staticVariables; // rhoBar.
        modified[2] = modif::staticVariables; // j.
        modified[3] = modif::nothing;         // Mass.
        modified[4] = modif::nothing;         // Volume fraction.
        modified[5] = modif::nothing;         // Flag-status.
        modified[6] = modif::nothing;         // Normal.
        modified[7] = modif::nothing;         // Interface-lists.
        modified[8] = modif::nothing;         // Curvature.
        modified[9] = modif::nothing;         // Outside density.
   }
};

template<typename T,template<typename U> class Descriptor>
class InitializeInterfaceLists3D : public BoxProcessingFunctional3D {
    virtual void processGenericBlocks(Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks)
    {
        PLB_ASSERT(atomicBlocks.size()==1);
        
        AtomicContainerBlock3D* containerInterfaceLists = dynamic_cast<AtomicContainerBlock3D*>(atomicBlocks[0]);
        PLB_ASSERT(containerInterfaceLists);
        InterfaceLists<T,Descriptor>* interfaceLists = new InterfaceLists<T,Descriptor>;
        containerInterfaceLists->setData(interfaceLists);

    }
    virtual InitializeInterfaceLists3D<T,Descriptor>* clone() const {
        return new InitializeInterfaceLists3D<T,Descriptor>(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const {
        // Default-assign potential other parameters present in a multi-fluid system.
        std::fill(modified.begin(), modified.end(), modif::nothing);
        modified[0] = modif::staticVariables;
    }
};

/// Wrapper for execution of InitializeInterfaceLists3D.
template<typename T,template<typename U> class Descriptor>
void initializeInterfaceLists3D(MultiContainerBlock3D& interfaceListBlock) {
    std::vector<MultiBlock3D*> arg;
    arg.push_back(&interfaceListBlock);
    applyProcessingFunctional (
            new InitializeInterfaceLists3D<T,Descriptor>,
            interfaceListBlock.getBoundingBox(), arg );
}

template<typename T, template<typename U> class Descriptor>
struct FreeSurfaceFields3D {
    static const int envelopeWidth;
    static const int smallEnvelopeWidth;
    static const int envelopeWidthForImmersedWalls;

    FreeSurfaceFields3D(SparseBlockStructure3D const& blockStructure,
                        Dynamics<T,Descriptor>* dynamics_,
                        T rhoDefault_, T surfaceTension_, T contactAngle_, Array<T,3> force_,
                        bool useImmersedWalls = false)
        : dynamics(dynamics_),
          rhoDefault(rhoDefault_), surfaceTension(surfaceTension_), contactAngle(contactAngle_), force(force_),
          lattice (
                MultiBlockManagement3D (
                        blockStructure, defaultMultiBlockPolicy3D().getThreadAttribution(),
                        smallEnvelopeWidth ),
                defaultMultiBlockPolicy3D().getBlockCommunicator(),
                defaultMultiBlockPolicy3D().getCombinedStatistics(),
                defaultMultiBlockPolicy3D().getMultiCellAccess<T,Descriptor>(), dynamics->clone() ),
          helperLists(lattice),
          mass(lattice),
          flag (
                MultiBlockManagement3D (
                        blockStructure,
                        defaultMultiBlockPolicy3D().getThreadAttribution(),
                        useImmersedWalls ? envelopeWidthForImmersedWalls : envelopeWidth ),
                defaultMultiBlockPolicy3D().getBlockCommunicator(),
                defaultMultiBlockPolicy3D().getCombinedStatistics(),
                defaultMultiBlockPolicy3D().getMultiScalarAccess<int>() ),
          volumeFraction((MultiBlock3D&) flag),
          curvature (
                MultiBlockManagement3D (
                        blockStructure,
                        defaultMultiBlockPolicy3D().getThreadAttribution(),
                        envelopeWidth ),
                defaultMultiBlockPolicy3D().getBlockCommunicator(),
                defaultMultiBlockPolicy3D().getCombinedStatistics(),
                defaultMultiBlockPolicy3D().getMultiScalarAccess<T>() ),
          outsideDensity((MultiBlock3D&) curvature),
          rhoBar (
                MultiBlockManagement3D (
                        blockStructure,
                        defaultMultiBlockPolicy3D().getThreadAttribution(),
                        useImmersedWalls ? envelopeWidthForImmersedWalls : smallEnvelopeWidth ),
                defaultMultiBlockPolicy3D().getBlockCommunicator(),
                defaultMultiBlockPolicy3D().getCombinedStatistics(),
                defaultMultiBlockPolicy3D().getMultiScalarAccess<T>() ),
          j (
                MultiBlockManagement3D (
                        blockStructure,
                        defaultMultiBlockPolicy3D().getThreadAttribution(),
                        useImmersedWalls ? envelopeWidthForImmersedWalls : smallEnvelopeWidth ),
                defaultMultiBlockPolicy3D().getBlockCommunicator(),
                defaultMultiBlockPolicy3D().getCombinedStatistics(),
                defaultMultiBlockPolicy3D().getMultiTensorAccess<T,3>() ),
          normal((MultiBlock3D&) curvature)
    {
        Precision precision = floatingPointPrecision<T>();
        T eps = getEpsilon<T>(precision);
        // The contact angle must take values between 0 and 180 degrees. If it is negative,
        // this means that contact angle effects will not be modeled.
        PLB_ASSERT(contactAngle < (T) 180.0 || std::fabs(contactAngle - (T) 180.0) <= eps);

        if (std::fabs(surfaceTension) <= eps) {
            useSurfaceTension = 0;
        } else {
            useSurfaceTension = 1;
        }

        twoPhaseArgs = aggregateFreeSurfaceParams(lattice, rhoBar, j, mass, volumeFraction,
                    flag, normal, helperLists, curvature, outsideDensity);

        initializeInterfaceLists3D<T,Descriptor>(helperLists);
        lattice.periodicity().toggleAll(true);
        mass.periodicity().toggleAll(true);
        flag.periodicity().toggleAll(true);
        volumeFraction.periodicity().toggleAll(true);
        curvature.periodicity().toggleAll(true);
        outsideDensity.periodicity().toggleAll(true);
        rhoBar.periodicity().toggleAll(true);
        j.periodicity().toggleAll(true);
        normal.periodicity().toggleAll(true);
        setToConstant(flag, flag.getBoundingBox(), (int)twoPhaseFlag::empty);
        setToConstant(outsideDensity, outsideDensity.getBoundingBox(), rhoDefault);
        rhoBarJparam.push_back(&lattice);
        rhoBarJparam.push_back(&rhoBar);
        rhoBarJparam.push_back(&j);

        lattice.internalStatSubscription().subscribeSum();     // Total mass.
        lattice.internalStatSubscription().subscribeSum();     // Lost mass.
        lattice.internalStatSubscription().subscribeIntSum();  // Num interface cells.

        freeSurfaceDataProcessors(rhoDefault, force, *dynamics);
        setExternalVector(lattice, lattice.getBoundingBox(), Descriptor<T>::ExternalField::forceBeginsAt, force);
    }
    FreeSurfaceFields3D(MultiBlockManagement3D const& blockManagement,
                        Dynamics<T,Descriptor>* dynamics_,
                        T rhoDefault_, T surfaceTension_, T contactAngle_, Array<T,3> force_,
                        bool useImmersedWalls = false)
        : dynamics(dynamics_),
          rhoDefault(rhoDefault_), surfaceTension(surfaceTension_), contactAngle(contactAngle_), force(force_),
          lattice (
                MultiBlockManagement3D (
                        blockManagement.getSparseBlockStructure(), blockManagement.getThreadAttribution().clone(),
                        smallEnvelopeWidth ),
                defaultMultiBlockPolicy3D().getBlockCommunicator(),
                defaultMultiBlockPolicy3D().getCombinedStatistics(),
                defaultMultiBlockPolicy3D().getMultiCellAccess<T,Descriptor>(), dynamics->clone() ),
          helperLists(lattice),
          mass(lattice),
          flag (
                MultiBlockManagement3D (
                        blockManagement.getSparseBlockStructure(),
                        blockManagement.getThreadAttribution().clone(),
                        useImmersedWalls ? envelopeWidthForImmersedWalls : envelopeWidth ),
                defaultMultiBlockPolicy3D().getBlockCommunicator(),
                defaultMultiBlockPolicy3D().getCombinedStatistics(),
                defaultMultiBlockPolicy3D().getMultiScalarAccess<int>() ),
          volumeFraction((MultiBlock3D&) flag),
          curvature (
                MultiBlockManagement3D (
                        blockManagement.getSparseBlockStructure(),
                        blockManagement.getThreadAttribution().clone(),
                        envelopeWidth ),
                defaultMultiBlockPolicy3D().getBlockCommunicator(),
                defaultMultiBlockPolicy3D().getCombinedStatistics(),
                defaultMultiBlockPolicy3D().getMultiScalarAccess<T>() ),
          outsideDensity((MultiBlock3D&) curvature),
          rhoBar (
                MultiBlockManagement3D (
                        blockManagement.getSparseBlockStructure(),
                        blockManagement.getThreadAttribution().clone(),
                        useImmersedWalls ? envelopeWidthForImmersedWalls : smallEnvelopeWidth ),
                defaultMultiBlockPolicy3D().getBlockCommunicator(),
                defaultMultiBlockPolicy3D().getCombinedStatistics(),
                defaultMultiBlockPolicy3D().getMultiScalarAccess<T>() ),
          j (
                MultiBlockManagement3D (
                        blockManagement.getSparseBlockStructure(),
                        blockManagement.getThreadAttribution().clone(),
                        useImmersedWalls ? envelopeWidthForImmersedWalls : smallEnvelopeWidth ),
                defaultMultiBlockPolicy3D().getBlockCommunicator(),
                defaultMultiBlockPolicy3D().getCombinedStatistics(),
                defaultMultiBlockPolicy3D().getMultiTensorAccess<T,3>() ),
          normal((MultiBlock3D&) curvature)
    {
        Precision precision = floatingPointPrecision<T>();
        T eps = getEpsilon<T>(precision);
        // The contact angle must take values between 0 and 180 degrees. If it is negative,
        // this means that contact angle effects will not be modeled.
        PLB_ASSERT(contactAngle < (T) 180.0 || std::fabs(contactAngle - (T) 180.0) <= eps);

        if (std::fabs(surfaceTension) <= eps) {
            useSurfaceTension = 0;
        } else {
            useSurfaceTension = 1;
        }

        twoPhaseArgs = aggregateFreeSurfaceParams(lattice, rhoBar, j, mass, volumeFraction,
                    flag, normal, helperLists, curvature, outsideDensity);

        initializeInterfaceLists3D<T,Descriptor>(helperLists);
        lattice.periodicity().toggleAll(true);
        mass.periodicity().toggleAll(true);
        flag.periodicity().toggleAll(true);
        volumeFraction.periodicity().toggleAll(true);
        curvature.periodicity().toggleAll(true);
        outsideDensity.periodicity().toggleAll(true);
        rhoBar.periodicity().toggleAll(true);
        j.periodicity().toggleAll(true);
        normal.periodicity().toggleAll(true);
        setToConstant(flag, flag.getBoundingBox(), (int)twoPhaseFlag::empty);
        setToConstant(outsideDensity, outsideDensity.getBoundingBox(), rhoDefault);
        rhoBarJparam.push_back(&lattice);
        rhoBarJparam.push_back(&rhoBar);
        rhoBarJparam.push_back(&j);

        lattice.internalStatSubscription().subscribeSum();     // Total mass.
        lattice.internalStatSubscription().subscribeSum();     // Lost mass.
        lattice.internalStatSubscription().subscribeIntSum();  // Num interface cells.

        freeSurfaceDataProcessors(rhoDefault, force, *dynamics);
        setExternalVector(lattice, lattice.getBoundingBox(), Descriptor<T>::ExternalField::forceBeginsAt, force);
    }

    FreeSurfaceFields3D(FreeSurfaceFields3D<T,Descriptor> const& rhs)
        : dynamics(rhs.dynamics->clone()),
          rhoDefault(rhs.rhoDefault),
          surfaceTension(rhs.surfaceTension),
          contactAngle(rhs.contactAngle),
          useSurfaceTension(rhs.useSurfaceTension),
          force(rhs.force),
          lattice(rhs.lattice),
          helperLists(rhs.helperLists),
          mass(rhs.mass),
          flag(rhs.flag),
          volumeFraction(rhs.volumeFraction),
          curvature(rhs.curvature),
          outsideDensity(rhs.outsideDensity),
          rhoBar(rhs.rhoBar),
          j(rhs.j),
          normal(rhs.normal),
          rhoBarJparam(rhs.rhoBarJparam),
          twoPhaseArgs(rhs.twoPhaseArgs)
    { }

    void swap(FreeSurfaceFields3D<T,Descriptor>& rhs)
    {
        std::swap(dynamics, rhs.dynamics);
        std::swap(rhoDefault, rhs.rhoDefault);
        std::swap(surfaceTension, rhs.surfaceTension);
        std::swap(contactAngle, rhs.contactAngle);
        std::swap(useSurfaceTension, rhs.useSurfaceTension);
        std::swap(force, rhs.force);
        std::swap(lattice, rhs.lattice);
        std::swap(helperLists, rhs.helperLists);
        std::swap(mass, rhs.mass);
        std::swap(flag, rhs.flag);
        std::swap(volumeFraction, rhs.volumeFraction);
        std::swap(curvature, rhs.curvature);
        std::swap(outsideDensity, rhs.outsideDensity);
        std::swap(rhoBar, rhs.rhoBar);
        std::swap(j, rhs.j);
        std::swap(normal, rhs.normal);
        std::swap(rhoBarJparam, rhs.rhoBarJparam);
        std::swap(twoPhaseArgs, rhs.twoPhaseArgs);
    }

    FreeSurfaceFields3D<T,Descriptor>& operator=(FreeSurfaceFields3D<T,Descriptor> const& rhs)
    {
        FreeSurfaceFields3D<T,Descriptor>(rhs).swap(*this);
        return *this;
    }

    FreeSurfaceFields3D<T,Descriptor>* clone() const
    {
        return new FreeSurfaceFields3D<T,Descriptor>(*this);
    }

    ~FreeSurfaceFields3D() {
        delete dynamics;
    }

    void periodicityToggle(plint direction, bool periodic)
    {
        PLB_ASSERT(direction == 0 || direction == 1 || direction == 2);

        lattice.periodicity().toggle(direction, periodic);
        mass.periodicity().toggle(direction, periodic);
        flag.periodicity().toggle(direction, periodic);
        volumeFraction.periodicity().toggle(direction, periodic);
        curvature.periodicity().toggle(direction, periodic);
        outsideDensity.periodicity().toggle(direction, periodic);
        rhoBar.periodicity().toggle(direction, periodic);
        j.periodicity().toggle(direction, periodic);
        normal.periodicity().toggle(direction, periodic);
    }

    void periodicityToggleAll(bool periodic)
    {
        lattice.periodicity().toggleAll(periodic);
        mass.periodicity().toggleAll(periodic);
        flag.periodicity().toggleAll(periodic);
        volumeFraction.periodicity().toggleAll(periodic);
        curvature.periodicity().toggleAll(periodic);
        outsideDensity.periodicity().toggleAll(periodic);
        rhoBar.periodicity().toggleAll(periodic);
        j.periodicity().toggleAll(periodic);
        normal.periodicity().toggleAll(periodic);
    }

    void defaultInitialize(bool useConstRho = true) {
        applyProcessingFunctional (
           new DefaultInitializeFreeSurface3D<T,Descriptor>( dynamics->clone(), force,
                                                             rhoDefault, useConstRho ),
                   lattice.getBoundingBox(), twoPhaseArgs );
    }

    void partiallyDefaultInitialize() {
        applyProcessingFunctional (
           new PartiallyDefaultInitializeFreeSurface3D<T,Descriptor>(dynamics->clone(), force, rhoDefault),
                   lattice.getBoundingBox(), twoPhaseArgs );
    }
    void freeSurfaceDataProcessors(T rhoDefault, Array<T,3> force, Dynamics<T,Descriptor>& dynamics)
    {
        plint pl; // Processor level.

        /***** Initial level ******/
        pl = 0;

        integrateProcessingFunctional (
                new ExternalRhoJcollideAndStream3D<T,Descriptor>,
                lattice.getBoundingBox(), rhoBarJparam, pl );

        integrateProcessingFunctional (
                new TwoPhaseComputeNormals3D<T,Descriptor>,
                lattice.getBoundingBox(), twoPhaseArgs, pl );

        /***** New level ******/
        pl++;

        if (useSurfaceTension) {
            integrateProcessingFunctional (
                    new TwoPhaseComputeCurvature3D<T,Descriptor>(contactAngle, lattice.getBoundingBox()),
                    lattice.getBoundingBox(), twoPhaseArgs, pl );
            
            // To change to the curvature calculation with height functions, uncomment the next data processor and
            // comment out the two previous ones. If only the next data processor is used and there is no
            // surface tension, the normals are not computed at all. Be careful if you intent to use
            // the normals and do not have the surface tension algorithm enabled.
            //integrateProcessingFunctional (
            //        new FreeSurfaceGeometry3D<T,Descriptor>(contactAngle),
            //        lattice.getBoundingBox(), twoPhaseArgs, pl );
        }

        integrateProcessingFunctional (
            new FreeSurfaceMassChange3D<T,Descriptor>, lattice.getBoundingBox(),
            twoPhaseArgs, pl );
       
        integrateProcessingFunctional (
            new FreeSurfaceCompletion3D<T,Descriptor>,
            lattice.getBoundingBox(), twoPhaseArgs, pl );
                                    
        integrateProcessingFunctional (
            new FreeSurfaceMacroscopic3D<T,Descriptor>(rhoDefault), 
            lattice.getBoundingBox(), twoPhaseArgs, pl );
        /***** New level ******/
        //pl++;

        //integrateProcessingFunctional (
        //        new FreeSurfaceInterfaceFilter<T,Descriptor>(),
        //        lattice.getBoundingBox(), twoPhaseArgs, pl );

        /***** New level ******/
        //pl++;

        //integrateProcessingFunctional (
        //        new FreeSurfaceInterfaceFilter<T,Descriptor>(),
        //        lattice.getBoundingBox(), twoPhaseArgs, pl );

        if (useSurfaceTension) {
            integrateProcessingFunctional (
                new TwoPhaseAddSurfaceTension3D<T,Descriptor>(surfaceTension, rhoDefault), 
                lattice.getBoundingBox(), twoPhaseArgs, pl );
        }

        /***** New level ******/
        pl++;

        integrateProcessingFunctional (
            new FreeSurfaceComputeInterfaceLists3D<T,Descriptor>(),
            lattice.getBoundingBox(), twoPhaseArgs, pl );

        integrateProcessingFunctional (
            new FreeSurfaceIniInterfaceToAnyNodes3D<T,Descriptor>(rhoDefault),
            lattice.getBoundingBox(), twoPhaseArgs, pl );
            
        integrateProcessingFunctional (
            new FreeSurfaceIniEmptyToInterfaceNodes3D<T,Descriptor>(dynamics.clone(), force),
                                    lattice.getBoundingBox(),
                                    twoPhaseArgs, pl ); 

        /***** New level ******/
        pl++;

        integrateProcessingFunctional (
            new FreeSurfaceRemoveFalseInterfaceCells3D<T,Descriptor>(rhoDefault),
            lattice.getBoundingBox(), twoPhaseArgs, pl);

        /***** New level ******/
        pl++;

        integrateProcessingFunctional (
            new FreeSurfaceEqualMassExcessReDistribution3D<T,Descriptor>(),
            //new FreeSurfaceWeightedMassExcessReDistribution3D<T,Descriptor>(),
            lattice.getBoundingBox(), twoPhaseArgs, pl );

        integrateProcessingFunctional (
            new TwoPhaseComputeStatistics3D<T,Descriptor>,
            lattice.getBoundingBox(), twoPhaseArgs, pl );
    }

    void appendBlocksToCheckpointVector(std::vector<MultiBlock3D*>& checkpointBlocks)
    {
        checkpointBlocks.push_back(&lattice);
        checkpointBlocks.push_back(&mass);
        checkpointBlocks.push_back(&flag);
        checkpointBlocks.push_back(&volumeFraction);
        checkpointBlocks.push_back(&outsideDensity);
        checkpointBlocks.push_back(&rhoBar);
        checkpointBlocks.push_back(&j);
    }

    Dynamics<T,Descriptor>* dynamics;
    T rhoDefault;
    T surfaceTension;
    T contactAngle;
    int useSurfaceTension;
    Array<T,3> force;
    MultiBlockLattice3D<T, Descriptor> lattice;
    MultiContainerBlock3D helperLists;
    MultiScalarField3D<T> mass;
    MultiScalarField3D<int> flag;
    MultiScalarField3D<T> volumeFraction;
    MultiScalarField3D<T> curvature;
    MultiScalarField3D<T> outsideDensity;
    MultiScalarField3D<T> rhoBar;
    MultiTensorField3D<T,3> j;
    MultiTensorField3D<T,3> normal;
    std::vector<MultiBlock3D*> rhoBarJparam;
    std::vector<MultiBlock3D*> twoPhaseArgs;
};

template<typename T, template<typename U> class Descriptor>
const int FreeSurfaceFields3D<T,Descriptor>::envelopeWidth = 3; // Necessary when we use height functions to compute the curvature,
                                                                // or when double smoothing is used at the data processor that
                                                                // computes the normals from the volume fraction.
//template<typename T, template<typename U> class Descriptor>
//const int FreeSurfaceFields3D<T,Descriptor>::envelopeWidth = 4; // Necessary when we use height functions to compute the curvature and
                                                                  // use the old contact angle algorithm.
template<typename T, template<typename U> class Descriptor>
const int FreeSurfaceFields3D<T,Descriptor>::smallEnvelopeWidth = 1;

template<typename T, template<typename U> class Descriptor>
const int FreeSurfaceFields3D<T,Descriptor>::envelopeWidthForImmersedWalls = 4;

}  // namespace plb

#endif  // FREE_SURFACE_MODEL_3D_H

