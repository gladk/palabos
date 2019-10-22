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
  * Flow in a lid-driven 3D cavity. The cavity is square and has no-slip walls,
  * except for the top wall which is diagonally driven with a constant
  * velocity. The benchmark is challenging because of the velocity
  * discontinuities on corner nodes.
  **/

#include "palabos3D.h"
#include "palabos3D.hh"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace plb;
using namespace std;

typedef double T;
#define DESCRIPTOR descriptors::D3Q19Descriptor

void cavitySetup( MultiBlockLattice3D<T,DESCRIPTOR>& lattice,
                  IncomprFlowParam<T> const& parameters,
                  OnLatticeBoundaryCondition3D<T,DESCRIPTOR>& boundaryCondition )
{
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz = parameters.getNz();
    Box3D topLid = Box3D(0, nx-1, ny-1, ny-1, 0, nz-1);
    Box3D everythingButTopLid = Box3D(0, nx-1, 0, ny-2, 0, nz-1);

    // All walls implement a Dirichlet velocity condition.
    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice);

    T u = std::sqrt((T)2)/(T)2 * parameters.getLatticeU();
    initializeAtEquilibrium(lattice, everythingButTopLid, (T)1., Array<T,3>((T)0.,(T)0.,(T)0.) );
    initializeAtEquilibrium(lattice, topLid, (T)1., Array<T,3>(u,(T)0.,u) );
    setBoundaryVelocity(lattice, topLid, Array<T,3>(u,(T)0.,u) );

    lattice.initialize();
}

template<class BlockLatticeT>
void writeGifs(BlockLatticeT& lattice,
               IncomprFlowParam<T> const& parameters, plint iter)
{
    const plint imSize = 600;
    const plint nx = parameters.getNx();
    const plint ny = parameters.getNy();
    const plint nz = parameters.getNz();
    const plint zComponent = 2;

    Box3D slice(0, nx-1, 0, ny-1, nz/2, nz/2);
    ImageWriter<T> imageWriter("leeloo");

    imageWriter.writeScaledGif( createFileName("uz", iter, 6),
                                *computeVelocityComponent (lattice, slice, zComponent),
                                imSize, imSize );

    imageWriter.writeScaledGif( createFileName("uNorm", iter, 6),
                                *computeVelocityNorm (lattice, slice),
                                imSize, imSize );
    imageWriter.writeScaledGif( createFileName("omega", iter, 6),
                                *computeNorm(*computeVorticity (
                                        *computeVelocity(lattice) ), slice ),
                                imSize, imSize );
}

template<class BlockLatticeT>
void writeVTK(BlockLatticeT& lattice,
              IncomprFlowParam<T> const& parameters, plint iter)
{
    T dx = parameters.getDeltaX();
    T dt = parameters.getDeltaT();
    VtkImageOutput3D<T> vtkOut(createFileName("vtk", iter, 6), dx);
    vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", dx/dt);
    vtkOut.writeData<3,float>(*computeVelocity(lattice), "velocity", dx/dt);
    vtkOut.writeData<3,float>(*computeVorticity(*computeVelocity(lattice)), "vorticity", 1./dt);
}


SparseBlockStructure3D createRegularDistribution3D (
        std::vector<plint> const& xVal, std::vector<plint> const& yVal, std::vector<plint> const& zVal )
{
    PLB_ASSERT(xVal.size()>=2);
    PLB_ASSERT(yVal.size()>=2);
    PLB_ASSERT(zVal.size()>=2);
    SparseBlockStructure3D dataGeometry (
            Box3D(xVal[0], xVal.back()-1, yVal[0], yVal.back()-1, zVal[0], zVal.back()-1) );
    for (plint iX=0; iX<(plint)xVal.size()-1; ++iX) {
        for (plint iY=0; iY<(plint)yVal.size()-1; ++iY) {
            for (plint iZ=0; iZ<(plint)zVal.size()-1; ++iZ) {
                plint nextID = dataGeometry.nextIncrementalId();
                Box3D domain( xVal[iX], xVal[iX+1]-1, yVal[iY],
                              yVal[iY+1]-1, zVal[iZ], zVal[iZ+1]-1 );
                dataGeometry.addBlock(domain, nextID);
                pcout << "Adding block with ID=" << nextID << ": ["
                      << domain.x0 << "," << domain.x1 << " | "
                      << domain.y0 << "," << domain.y1 << " | "
                      << domain.z0 << "," << domain.z1 << "]" << std::endl;
            }
        }
    }
    return dataGeometry;
}

SparseBlockStructure3D createCavityDistribution3D(plint nx, plint ny, plint nz)
{
    static const plint numval=4;
    plint x[numval] = {0, nx/4, 3*nx/4, nx};
    plint y[numval] = {0, ny/4, 3*ny/4, ny};
    plint z[numval] = {0, nz/4, 3*nz/4, nz};
    std::vector<plint> xVal(x, x+numval);
    std::vector<plint> yVal(y, y+numval);
    std::vector<plint> zVal(z, z+numval);
    return createRegularDistribution3D(xVal, yVal, zVal);
}


void saveAtomicBlock(MultiBlockLattice3D<T,DESCRIPTOR>& lattice, plint blockId) {
    BlockLattice3D<T,DESCRIPTOR>& atomicBlock = lattice.getComponent(blockId);
    plb_ofstream ofile("block13.dat");
    ofile << atomicBlock;
}

int main(int argc, char* argv[]) {

    plbInit(&argc, &argv);
    global::directories().setOutputDir("./tmp/");

    IncomprFlowParam<T> parameters(
            (T) 1e-4,  // uMax
            (T) 10.,   // Re
            30,        // N
            1.,        // lx
            1.,        // ly
            1.         // lz
    );
    const T logT     = (T)1/(T)1000;
    const T imSave   = (T)1/(T)40;
    const T vtkSave  = (T)1;
    const T maxT     = (T)10.1;

    pcout << "omega= " << parameters.getOmega() << std::endl;
    writeLogFile(parameters, "3D diagonal cavity");

    // @Tomasz: STEP 1
    // Instead of simply creating a MultiBlockLattice3D as usual, the internal
    // structure and parallelization of the block are created manually. The block
    // is covered by 5x5x5 sub-domains. Each of these 125 sub-domains, if it
    // contains no boundary area (i.e. only pure BGK), will then be off-loaded to
    // the co-processor.
    
    // Here the 5x5x5 cover-up is instantiated.
    plint numBlocksX = 3;
    plint numBlocksY = 3;
    plint numBlocksZ = 3;
    plint numBlocks = numBlocksX*numBlocksY*numBlocksZ;
    plint envelopeWidth = 1;
    SparseBlockStructure3D blockStructure (
            createCavityDistribution3D (
                parameters.getNx(), parameters.getNy(), parameters.getNz() ) );

    // In case of MPI parallelism, the blocks are explicitly assigned to processors,
    // with equal load.
    ExplicitThreadAttribution* threadAttribution = new ExplicitThreadAttribution;
    std::vector<std::pair<plint,plint> > ranges;
    plint numRanges = std::min(numBlocks, (plint)global::mpi().getSize());
    util::linearRepartition(0, numBlocks-1, numRanges, ranges);
    for (pluint iProc=0; iProc<ranges.size(); ++iProc) {
        for (plint blockId=ranges[iProc].first; blockId<=ranges[iProc].second; ++blockId) {
            threadAttribution -> addBlock(blockId, iProc);
        }
    }


    // Create a lattice with the above specified internal structure.
    MultiBlockLattice3D<T, DESCRIPTOR> lattice (
        MultiBlockManagement3D (
                        blockStructure, threadAttribution, envelopeWidth ),
        defaultMultiBlockPolicy3D().getBlockCommunicator(),
                    defaultMultiBlockPolicy3D().getCombinedStatistics(),
                    defaultMultiBlockPolicy3D().getMultiCellAccess<T,DESCRIPTOR>(),
                    new BGKdynamics<T,DESCRIPTOR>(parameters.getOmega()) );

    saveAtomicBlock(lattice, 13);

    OnLatticeBoundaryCondition3D<T,DESCRIPTOR>* boundaryCondition
        = createInterpBoundaryCondition3D<T,DESCRIPTOR>();
    cavitySetup(lattice, parameters, *boundaryCondition);

    // @Tomasz: STEP 2
    // Here the co-processor is informed about the domains for which it will need to do computations.
    bool printInfo=true;
    initiateCoProcessors(lattice, BGKdynamics<T,DESCRIPTOR>(parameters.getOmega()).getId(), printInfo);

    T previousIterationTime = T();
    // Loop over main time iteration.
    for (plint iT=0; iT<parameters.nStep(maxT); ++iT) {
        global::timer("mainLoop").restart();

        if (iT%parameters.nStep(imSave)==0) {
            pcout << "Writing Gif ..." << endl;
            writeGifs(lattice, parameters, iT);
        }

        if (iT%parameters.nStep(vtkSave)==0 && iT>0) {
            pcout << "Saving VTK file ..." << endl;
            writeVTK(lattice, parameters, iT);
        }

        if (iT%parameters.nStep(logT)==0) {
            pcout << "step " << iT
                  << "; t=" << iT*parameters.getDeltaT();
        }

        // @Tomasz: STEP 3
        // For now, all the data assigned to co-processors is transferred from CPU memory
        // to co-processor memory before the collision-streaming step, and back to CPU
        // memory after the collision-streaming step. This is not efficient at all: we
        // will replace this by communicating outer boundary layers only. It is however
        // simpler for now, for our first proof of concept.

        // Execute a time iteration.
        transferToCoProcessors(lattice);
        lattice.collideAndStream();
        transferFromCoProcessors(lattice);
        // After transferring back from co-processor to CPU memory, data in the envelopes
        // of the CPU blocks must be synchronized.
        lattice.duplicateOverlaps(modif::staticVariables);


        // Access averages from internal statistics ( their value is defined
        //   only after the call to lattice.collideAndStream() )
        if (iT%parameters.nStep(logT)==0) {
            pcout << "; av energy="
                  << setprecision(10) << computeAverageEnergy<T>(lattice)
                  << "; av rho="
                  << setprecision(10) << computeAverageDensity<T>(lattice) << endl;
            pcout << "Time spent during previous iteration: "
                  << previousIterationTime << endl;
        }

        previousIterationTime = global::timer("mainLoop").stop();
    }

    delete boundaryCondition;
}
