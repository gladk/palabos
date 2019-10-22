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
 *  This example illustrates a simple multi-phase dynamic mixer.
 *  The impeller is a small rectangle which rotates with constant
 *  angular velocity inside a mixer domain. This example serves to
 *  demonstrate the combination of the Shan-Chen multi-component
 *  model, with the immersed boundary method.
 **/

#include "palabos3D.h"
#include "palabos3D.hh"

#include <cstdlib>
#include <iostream>

using namespace plb;
using namespace std;

// Use double-precision arithmetics
typedef double T;
// Use a grid which additionally to the f's stores three variables for
//   the external force term.
#define DESCRIPTOR descriptors::ForcedD3Q19Descriptor

// Global variables that relate to the placement and movement of the
// immersed surface (the impleller of this simple mixer).
static Array<T,3> rectangleCenter;
static T angVel;

// Directory for output of the results
static std::string outDir("./tmp/");

/// Initial condition: fluid1 on top, fluid2 on bottom.
/** This functional is going to be used as an argument to the function "applyIndexed",
 *  to setup the initial condition. For efficiency reasons, this approach should
 *  always be preferred over explicit space loops in end-user codes.
 */
template<typename T, template<typename U> class Descriptor>
class TwoLayerInitializer : public OneCellIndexedWithRandFunctional3D<T,Descriptor> {
public:
    TwoLayerInitializer(plint h_, bool topLayer_)
        : h(h_),
          topLayer(topLayer_)
    { }
    TwoLayerInitializer<T,Descriptor>* clone() const {
        return new TwoLayerInitializer<T,Descriptor>(*this);
    }
    virtual void execute(plint iX, plint iY, plint iZ, T rand_val, Cell<T,Descriptor>& cell) const {
        T densityFluctuations = 0.0; //1.e-2;
        T almostNoFluid       = 1.e-4;
        Array<T,3> zeroVelocity ((T) 0., (T) 0., (T) 0.);

        T rho = (T)1;
        // Add a random perturbation to the initial condition (not always necessary).
        if ( (topLayer && iY>h) || (!topLayer && iY <= h) ) {
            rho += rand_val * densityFluctuations;
        }
        else {
            rho = almostNoFluid;
        }

        iniCellAtEquilibrium(cell, rho, zeroVelocity);
    }
private:
    plint h;
    bool topLayer;
};

// Geometry of the mixer domain (not the impeller).
template<typename T>
class MixerShapeDomain3D : public DomainFunctional3D {
public:
    MixerShapeDomain3D(plint cx_, plint cy_, plint radius)
        : cx(cx_),
          cy(cy_),
          radiusSqr(util::sqr(radius))
    { }
    virtual bool operator() (plint iX, plint iY, plint iZ) const
    {
        return iX <= cx && util::sqr(iX-cx) + util::sqr(iY-cy) > radiusSqr;
    }
    virtual MixerShapeDomain3D<T>* clone() const {
        return new MixerShapeDomain3D<T>(*this);
    }
private:
    plint cx;
    plint cy;
    plint radiusSqr;
};

void mixerSetup( MultiBlockLattice3D<T, DESCRIPTOR>& fluid1,
                 MultiBlockLattice3D<T, DESCRIPTOR>& fluid2,
                 T force )
{
    plint nx = fluid1.getNx();
    plint ny = fluid1.getNy();
    plint nz = fluid1.getNz();
   
    // Initialize top layer.
    applyIndexed(fluid1, Box3D(0, nx-1, 0, ny-1, 0, nz-1),
                 new TwoLayerInitializer<T,DESCRIPTOR>(3*ny/4, true) );
    // Initialize bottom layer.
    applyIndexed(fluid2, Box3D(0, nx-1, 0, ny-1, 0, nz-1),
                 new TwoLayerInitializer<T,DESCRIPTOR>(3*ny/4, false) );

    // Specify the mixer domain around the impeller.
    plint cx = util::roundToInt((T) 0.5 * (nx - 1));
    plint cy = util::roundToInt((T) 0.5 * (ny - 1));
    plint radius = util::roundToInt((T) 0.5 * (ny - 1));
    defineDynamics(fluid1, fluid1.getBoundingBox(), new MixerShapeDomain3D<T>(cx,cy,radius), new BounceBack<T, DESCRIPTOR>);
    defineDynamics(fluid2, fluid2.getBoundingBox(), new MixerShapeDomain3D<T>(cx,cy,radius), new BounceBack<T, DESCRIPTOR>);

    // Let's have gravity acting on fluid1 only. This represents a situation
    //   where the molecular mass of fluid2 is very small, and thus the
    //   action of gravity on this species is negligible.
    setExternalVector(fluid1, fluid1.getBoundingBox(),
                      DESCRIPTOR<T>::ExternalField::forceBeginsAt, Array<T,3>((T)0.,-force,(T)0.));
    setExternalVector(fluid2, fluid2.getBoundingBox(),
                      DESCRIPTOR<T>::ExternalField::forceBeginsAt, Array<T,3>((T)0.,(T)0.,(T)0.));
}

// The angle of rotation of the immersed rectangle is given as a linear function of time.
T getAngle(T t)
{
    return angVel * t;
}

// The angular velocity of the immersed rectangle is the time derivative of the angle function.
T getAngularVelocity(T t)
{
    return angVel;
}

// This class computes the velocity of the immersed surface at every
// time step. The immersed surface has a defined position given by 
// an angle of the form:
//   phi = angularVelocity * t
// So its velocity is given by the cross product between its 
// angular velocity (time derivative of the angle) and its position.
// The axis of rotation is the z-axis.
class SurfaceVelocity {
public:
    SurfaceVelocity(T t_)
        : t(t_)
    { }
    Array<T,3> operator()(Array<T,3> const& pos)
    {
        T xp = pos[0] - rectangleCenter[0];
        T yp = pos[1] - rectangleCenter[1];
        T angularVelocity = getAngularVelocity(t);
        Array<T,3> velocity;
        velocity[0] = -angularVelocity * yp;
        velocity[1] =  angularVelocity * xp;
        velocity[2] =  0.0;
        return velocity;
    }
private:
    T t;
};

void writePPM(MultiBlockLattice3D<T, DESCRIPTOR>& fluid1,
              MultiBlockLattice3D<T, DESCRIPTOR>& fluid2, plint iT)
{
    const plint nx = fluid1.getNx();
    const plint ny = fluid1.getNy();
    const plint nz = fluid1.getNz();
    Box3D slice(0, nx-1, 0, ny-1, nz/2, nz/2);

    ImageWriter<T> imageWriter("leeloo");
    imageWriter.writeScaledPpm(createFileName("rho1_", iT, 6), *computeDensity(fluid1, slice));
    imageWriter.writeScaledPpm(createFileName("rho2_", iT, 6), *computeDensity(fluid2, slice));
}

void writeVTK(MultiBlockLattice3D<T, DESCRIPTOR>& fluid1,
              MultiBlockLattice3D<T, DESCRIPTOR>& fluid2, plint iT)
{
    const plint nx = fluid1.getNx();
    const plint ny = fluid1.getNy();
    const plint nz = fluid1.getNz();

    plint cx = nx/2;
    plint cy = ny/2;
    plint cz = nz/2;

    Box3D box_x(cx - 1, cx + 1, 0, ny - 1, 0, nz - 1);
    Box3D box_y(0, nx - 1, cy - 1, cy + 1, 0, nz - 1);
    Box3D box_z(0, nx - 1, 0, ny - 1, cz - 1, cz + 1);

    {
        VtkImageOutput3D<T> vtkOut(createFileName("slice_x_", iT, 6), 1.0);
        std::auto_ptr<MultiScalarField3D<T> > rho1 = computeDensity(fluid1, box_x);
        std::auto_ptr<MultiScalarField3D<T> > rho2 = computeDensity(fluid2, box_x);
        vtkOut.writeData<float>(*rho1, "rho1", 1.0);
        vtkOut.writeData<float>(*rho2, "rho2", 1.0);
    }
    {
        VtkImageOutput3D<T> vtkOut(createFileName("slice_y_", iT, 6), 1.0);
        std::auto_ptr<MultiScalarField3D<T> > rho1 = computeDensity(fluid1, box_y);
        std::auto_ptr<MultiScalarField3D<T> > rho2 = computeDensity(fluid2, box_y);
        vtkOut.writeData<float>(*rho1, "rho1", 1.0);
        vtkOut.writeData<float>(*rho2, "rho2", 1.0);
    }
    {
        VtkImageOutput3D<T> vtkOut(createFileName("slice_z_", iT, 6), 1.0);
        std::auto_ptr<MultiScalarField3D<T> > rho1 = computeDensity(fluid1, box_z);
        std::auto_ptr<MultiScalarField3D<T> > rho2 = computeDensity(fluid2, box_z);
        vtkOut.writeData<float>(*rho1, "rho1", 1.0);
        vtkOut.writeData<float>(*rho2, "rho2", 1.0);
    }

    {
        VtkImageOutput3D<T> vtkOut(createFileName("full_domain_", iT, 6), 1.0);
        std::auto_ptr<MultiScalarField3D<T> > rho1 = computeDensity(fluid1);
        std::auto_ptr<MultiScalarField3D<T> > rho2 = computeDensity(fluid2);
        vtkOut.writeData<float>(*rho1, "rho1", 1.0);
        vtkOut.writeData<float>(*rho2, "rho2", 1.0);
    }
}

// Write an STL file of the instantaneous geometry of the
// moving immersed surface.
void writeSTL(TriangleSet<T> triangleSet, plint iT)
{
    static Array<T,3> normedAxis((T) 0, (T) 0, (T) 1);

    T initialAngle = 0.0;
    T angle = getAngle(iT) - initialAngle;
    triangleSet.translate(-rectangleCenter);
    triangleSet.rotateAtOrigin(normedAxis, angle);
    triangleSet.translate( rectangleCenter);
    triangleSet.writeBinarySTL(createFileName(outDir + "surface_", iT, 6) + ".stl");
}

int main(int argc, char *argv[])
{
    plbInit(&argc, &argv);
    global::directories().setOutputDir(outDir);
    srand(global::mpi().getRank());
    
    const T omega1 = 1.0;
    const T omega2 = 1.0;
    const plint nx   = 75;
    const plint ny   = 75;
    const plint nz   = 75;
    const T G      = 2.0;
    T force        = 0.0; // We neglect gravity.
    const plint maxIter  = 40000;
    const plint saveIter = 100;
    const plint statIter = 50;

    T pi   = std::acos((T) -1);
    angVel = 2.0 * pi / 10000.0;

    plint wallThickness      = 2;
    plint largeEnvelopeWidth = 4;

    // Use regularized BGK dynamics to improve numerical stability (but note that
    //   BGK dynamics works well too).
    MultiBlockLattice3D<T, DESCRIPTOR> fluid1(nx,ny,nz, new RegularizedBGKdynamics<T, DESCRIPTOR>(omega1));
    bool incompressibleModel1 = false;
    defineDynamics(fluid1, fluid1.getBoundingBox(), new BounceBack<T, DESCRIPTOR>());
    defineDynamics(fluid1, fluid1.getBoundingBox().enlarge(-wallThickness), new RegularizedBGKdynamics<T, DESCRIPTOR>(omega1));
    MultiScalarField3D<T>* rhoBar1 = generateMultiScalarField<T>((MultiBlock3D&) fluid1, largeEnvelopeWidth).release();
    MultiTensorField3D<T,3>* j1 = generateMultiTensorField<T,3>((MultiBlock3D&) fluid1, largeEnvelopeWidth).release();

    MultiBlockLattice3D<T, DESCRIPTOR> fluid2(nx,ny,nz, new RegularizedBGKdynamics<T, DESCRIPTOR>(omega2));
    bool incompressibleModel2 = false;
    defineDynamics(fluid2, fluid2.getBoundingBox(), new BounceBack<T, DESCRIPTOR>());
    defineDynamics(fluid2, fluid2.getBoundingBox().enlarge(-wallThickness), new RegularizedBGKdynamics<T, DESCRIPTOR>(omega2));
    MultiScalarField3D<T>* rhoBar2 = generateMultiScalarField<T>((MultiBlock3D&) fluid2, largeEnvelopeWidth).release();
    MultiTensorField3D<T,3>* j2 = generateMultiTensorField<T,3>((MultiBlock3D&) fluid2, largeEnvelopeWidth).release();

    fluid1.periodicity().toggleAll(false);
    rhoBar1->periodicity().toggleAll(false);
    j1->periodicity().toggleAll(false);

    fluid2.periodicity().toggleAll(false);
    rhoBar2->periodicity().toggleAll(false);
    j2->periodicity().toggleAll(false);

    vector<MultiBlock3D*> blocksFluid1;
    blocksFluid1.push_back(&fluid1);
    blocksFluid1.push_back(rhoBar1);
    blocksFluid1.push_back(j1);
    integrateProcessingFunctional(new ExternalRhoJcollideAndStream3D<T,DESCRIPTOR>(), fluid1.getBoundingBox(), blocksFluid1, 0);

    vector<MultiBlock3D*> blocksFluid2;
    blocksFluid2.push_back(&fluid2);
    blocksFluid2.push_back(rhoBar2);
    blocksFluid2.push_back(j2);
    integrateProcessingFunctional(new ExternalRhoJcollideAndStream3D<T,DESCRIPTOR>(), fluid2.getBoundingBox(), blocksFluid2, 0);
    
    // Store a pointer to all blocks (six in the present application) in a vector to
    //   create the Shan/Chen coupling term. fluid1 being at the first place
    //   in the vector, the coupling term is going to be executed at the end of the call
    //   to executeInternalProcessors() for fluid1.
    vector<MultiBlock3D*> blocks;
    blocks.push_back(&fluid1);
    blocks.push_back(rhoBar1);
    blocks.push_back(j1);
    blocks.push_back(&fluid2);
    blocks.push_back(rhoBar2);
    blocks.push_back(j2);
    
    // The argument "constOmegaValues" to the Shan/Chen processor is optional,
    //   and is used for efficiency reasons only. It tells the data processor
    //   that the relaxation times are constant, and that their inverse must be
    //   computed only once.
    std::vector<T> constOmegaValues;
    constOmegaValues.push_back(omega1);
    constOmegaValues.push_back(omega2);
    plint processorLevel = 1;
    integrateProcessingFunctional (
            new ShanChenExternalMultiComponentProcessor3D<T,DESCRIPTOR>(G,constOmegaValues),
            fluid1.getBoundingBox().enlarge(-wallThickness),
            blocks, processorLevel );

    mixerSetup(fluid1, fluid2, force);

    // Immersed surface.

    pcout << "Creating the immersed rectangle surface." << std::endl;

    // The immersed boundary method needs a set of vertices and a set of areas
    // that correspond to each vertex. These vertices and areas describe the
    // time dependent geometry of the surface at each time step.
    std::vector<Array<T,3> > vertices;
    std::vector<T> areas;

    // The initial geometry of the immersed surface can be created analyticaly or, alternatively,
    // it can be loaded by reading a user provided STL file.
    T rectangleLx = (T) 2 * (util::roundToInt((T) 0.5 * (T) (ny - 1)) - (T) (wallThickness + 1));
    T rectangleLy = rectangleLx;
    plint rectangleNx = util::roundToInt((T) 2 * rectangleLx);
    plint rectangleNy = util::roundToInt((T) 2 * rectangleLy);
    rectangleCenter = Array<T,3>((T) 0.5 * (T) (nx - 1), (T) 0.5 * (T) (ny - 1), (T) 0.5 * (T) (nz - 1));
    Array<T,3> rectangleNormal((T) 0, (T) 1, (T) 0);
    TriangleSet<T> rectangleTriangleSet = constructGenericRectangle<T>(rectangleLx, rectangleLy,
            rectangleNx, rectangleNy, rectangleCenter, rectangleNormal);

    DEFscaledMesh<T> *rectangleDef = new DEFscaledMesh<T>(rectangleTriangleSet, 0, 0, 0, Dot3D(0, 0, 0));
    pcout << "The rectangle has " << rectangleDef->getMesh().getNumVertices() << " vertices and " <<
        rectangleDef->getMesh().getNumTriangles() << " triangles." << std::endl;
    for (pluint iVertex = 0; iVertex < (pluint) rectangleDef->getMesh().getNumVertices(); iVertex++) {
        vertices.push_back(rectangleDef->getMesh().getVertex(iVertex));
        areas.push_back(rectangleDef->getMesh().computeVertexArea(iVertex));
    }
    delete rectangleDef;

    // The next container block is necessary for the immersed-wall algorithm.
    MultiContainerBlock3D container(*rhoBar1);

    // Initialization.

    applyProcessingFunctional (
            new ShanChenExternalMultiComponentProcessor3D<T,DESCRIPTOR>(G,constOmegaValues),
            fluid1.getBoundingBox().enlarge(-wallThickness), blocks );

    // Start of iterations.
	
    pcout << "Starting simulation" << endl;
    // Main loop over time iterations.
    for (plint iT=0; iT<maxIter; ++iT) {
        if (iT%saveIter==0) {
            writePPM(fluid1, fluid2, iT);
            writeVTK(fluid1, fluid2, iT);
            writeSTL(rectangleTriangleSet, iT);
        }

        // Time iteration for fluid2.
        fluid2.executeInternalProcessors();
        // Time iteration for fluid1 must come after fluid2, because the coupling is
        //   executed here. You should understand this as follows.
        //   The effect of the coupling is to compute the interaction force between
        //   species, and to precompute density and momentum for each species. This must
        //   be executed *before* collide-and-streaming the fluids, because the collision
        //   step needs to access all these values. In the present case, it is done after
        //   both collide-and-stream step, which means, before the collide-and-stream of
        //   the next iteration (it's the same if you are before or after; the important
        //   point is not to be between the two collide-and-streams of fluid2 and fluid1
        //   As for the initial condition, the coupling is initially performed once
        //   before the start of the main time loop.
        fluid1.executeInternalProcessors();

        if (iT%statIter==0) {
            pcout << "At iteration " << iT << ": Average energy fluid one = "
                  << computeAverageEnergy<T>(fluid1);
            pcout << ", average energy fluid two = "
                  << computeAverageEnergy<T>(fluid2) << endl;
        }

        // Immersed walls algorithm.

        T currentTime = iT;
        T nextTime = iT + 1.0;

        // New position of the immersed surface at the next time iteration.
        // The angle of rotation is given by a function of the form:
        //   phi = angularVelocity * t
        // The rotation axis is the z-axis.
        for (plint iVertex = 0; iVertex < (plint) vertices.size(); iVertex++) {
            T x = vertices[iVertex][0] - rectangleCenter[0];
            T y = vertices[iVertex][1] - rectangleCenter[1];
            T dphi = getAngle(nextTime) - getAngle(currentTime);
            T c = std::cos(dphi);
            T s = std::sin(dphi);
            vertices[iVertex][0] = x * c - y * s + rectangleCenter[0];
            vertices[iVertex][1] = x * s + y * c + rectangleCenter[1];
        }

        // Instantiate the immersed wall data and performed the immersed boundary iterations.
        instantiateImmersedWallData(vertices, areas, container);
        plint ibIter = 6;
        for (plint i = 0; i < ibIter; i++) {
            inamuroIteration(SurfaceVelocity(nextTime),
                    *rhoBar1, *j1, container, (T) 1.0 / omega1, incompressibleModel1);
        }
        for (plint i = 0; i < ibIter; i++) {
            inamuroIteration(SurfaceVelocity(nextTime),
                    *rhoBar2, *j2, container, (T) 1.0 / omega2, incompressibleModel2);
        }
    }

    delete j2;
    delete rhoBar2;
    delete j1;
    delete rhoBar1;

    return 0;
}
