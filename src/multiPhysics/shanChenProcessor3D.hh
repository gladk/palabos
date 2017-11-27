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

/* The original version of this file was written by Orestis Malaspinas
 * and Andrea Parmigiani.
 */

#ifndef SHAN_CHEN_PROCESSOR_3D_HH
#define SHAN_CHEN_PROCESSOR_3D_HH

#include "multiPhysics/shanChenProcessor3D.h"
#include "core/util.h"
#include "finiteDifference/finiteDifference3D.h"
#include "latticeBoltzmann/momentTemplates.h"
#include "latticeBoltzmann/externalFieldAccess.h"
#include "multiPhysics/multiPhaseTemplates3D.h"

namespace plb {

/* *************** ShanChenMultiComponentProcessor3D ***************** */

template<typename T, template<typename U> class Descriptor>
ShanChenMultiComponentProcessor3D <T,Descriptor>::ShanChenMultiComponentProcessor3D(T G_)
    : G(G_)
{ }

template<typename T, template<typename U> class Descriptor>
ShanChenMultiComponentProcessor3D <T,Descriptor>::ShanChenMultiComponentProcessor3D(
        std::vector<std::vector<T> > const& speciesG_)
    : G((T) 0)
{
    pluint numSpecies = speciesG_.size();
    // Although speciesG_ has a 2D "matrix structure", speciesG has a 1D "array structure".
    speciesG.resize(numSpecies * numSpecies);
    for (pluint iSpecies = 0; iSpecies < numSpecies; iSpecies++) {
        PLB_ASSERT(speciesG_[iSpecies].size() == numSpecies);
        for (pluint jSpecies = 0; jSpecies < numSpecies; jSpecies++) {
            speciesG[iSpecies * numSpecies + jSpecies] = speciesG_[iSpecies][jSpecies];
        }
    }
}

template<typename T, template<typename U> class Descriptor>
ShanChenMultiComponentProcessor3D <T,Descriptor>::ShanChenMultiComponentProcessor3D (
        T G_, std::vector<T> const& imposedOmega_)
    : G(G_),
      imposedOmega(imposedOmega_)
{ }

template<typename T, template<typename U> class Descriptor>
ShanChenMultiComponentProcessor3D <T,Descriptor>::ShanChenMultiComponentProcessor3D (
        std::vector<std::vector<T> > const& speciesG_, std::vector<T> const& imposedOmega_)
    : G((T) 0),
      imposedOmega(imposedOmega_)
{
    pluint numSpecies = speciesG_.size();
    // Although speciesG_ has a 2D "matrix structure", speciesG has a 1D "array structure".
    speciesG.resize(numSpecies * numSpecies);
    for (pluint iSpecies = 0; iSpecies < numSpecies; iSpecies++) {
        PLB_ASSERT(speciesG_[iSpecies].size() == numSpecies);
        for (pluint jSpecies = 0; jSpecies < numSpecies; jSpecies++) {
            speciesG[iSpecies * numSpecies + jSpecies] = speciesG_[iSpecies][jSpecies];
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void ShanChenMultiComponentProcessor3D<T,Descriptor>::process (
        Box3D domain,
        std::vector<BlockLattice3D<T,Descriptor>*> lattices )
{
    // Number of species (or components) which are coupled in this Shan/Chen multi-component fluid.
    plint numSpecies = (plint) lattices.size();
    // Short-hand notation for the lattice descriptor
    typedef Descriptor<T> D;
    // Handle to external scalars
    enum {
        densityOffset  = D::ExternalField::densityBeginsAt,
        momentumOffset = D::ExternalField::momentumBeginsAt
    };
    
    // Compute per-lattice density  and momentum on every site and on each
    //   lattice, and store result in external scalars;  envelope cells are included,
    //   because they are needed to compute the interaction potential in the following.
    //   Note that the per-lattice value of the momentum is stored temporarily only, as
    //   it is corrected later on, based on the common fluid velocity.
    for (plint iSpecies=0; iSpecies<numSpecies; ++iSpecies) {
        for (plint iX=domain.x0-1; iX<=domain.x1+1; ++iX) {
            for (plint iY=domain.y0-1; iY<=domain.y1+1; ++iY) {
                for (plint iZ=domain.z0-1; iZ<=domain.z1+1; ++iZ) {
                    // Get "intelligent" value of density through cell object, to account
                    //   for the fact that the density value can be user-defined, for example
                    //   on boundaries.
                    Cell<T,Descriptor>& cell = lattices[iSpecies]->get(iX,iY,iZ);
                    Array<T,Descriptor<T>::d> j;
                    T rhoBar = cell.getDynamics().computeRhoBar(cell);
                    momentTemplates<T,Descriptor>::get_j(cell,j);
                    *cell.getExternal(densityOffset) = Descriptor<T>::fullRho(rhoBar);
                    j.to_cArray(cell.getExternal(momentumOffset));
                }
            }
        }
    }

    // Temporary variable for the relaxation parameters omega.
    std::vector<T> omega(numSpecies), invOmega(numSpecies);
    // Temporary variable for total velocity.
    Array<T,Descriptor<T>::d> uTot;
    // Temporary variable for interaction potential.
    std::vector<Array<T,D::d> > rhoContribution(numSpecies);

    // If omega is constant and imposed by the user, copy its value to
    //   the vector "omega", and compute the inverse.
    if (!imposedOmega.empty()) {
        PLB_ASSERT((plint) imposedOmega.size() == numSpecies);
        omega = imposedOmega;
        for (pluint iOmega=0; iOmega<omega.size(); ++iOmega) {
            invOmega[iOmega] = (T)1 / omega[iOmega];
        }
    }

    if (speciesG.empty()) {
        speciesG.resize(numSpecies * numSpecies, G);
    }

    // Compute the interaction force between the species, and store it by
    //   means of a velocity correction in the external velocity field.
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                // Computation of the common density over all populations, weighted by
                //   the relaxation parameters omega.
                T weightedDensity = T();
                for (plint iSpecies=0; iSpecies<numSpecies; ++iSpecies) {
                    Cell<T,Descriptor> const& cell = lattices[iSpecies]->get(iX,iY,iZ);
                    // Take this opportunity to read omega from the cell, unless the value
                    //   of omega is constant and imposed by the user.
                    if (imposedOmega.empty()) {
                        omega[iSpecies] = cell.getDynamics().getOmega();
                        invOmega[iSpecies] = (T)1/omega[iSpecies];
                    }
                    weightedDensity += omega[iSpecies] * (*cell.getExternal(densityOffset));
                }
                // Computation of the common velocity, shared among all populations.
                for (int iD = 0; iD < Descriptor<T>::d; ++iD) {
                    uTot[iD] = T();
                    for (plint iSpecies=0; iSpecies<numSpecies; ++iSpecies) {
                        T *momentum = lattices[iSpecies]->get(iX,iY,iZ).getExternal(momentumOffset);
                        uTot[iD] += momentum[iD] * omega[iSpecies];
                    }
                    uTot[iD] /= weightedDensity;
                }

                // Computation of the interaction potential.
                for (plint iSpecies=0; iSpecies<numSpecies; ++iSpecies) {
                    multiPhaseTemplates3D<T,Descriptor>::shanChenInteraction (
                            *lattices[iSpecies],rhoContribution[iSpecies],iX,iY,iZ );
                }

                // Computation and storage of the final velocity, consisting
                //   of uTot plus the momentum difference due to interaction
                //   potential and external force
                for (plint iSpecies=0; iSpecies<numSpecies; ++iSpecies) {
                    Cell<T,Descriptor>& cell = lattices[iSpecies]->get(iX,iY,iZ);
                    T *momentum = cell.getExternal(momentumOffset);
                    for (int iD = 0; iD < D::d; ++iD) {
                        momentum[iD] = uTot[iD];
                        // Initialize force contribution with force from external fields if there
                        //   is any, or with zero otherwise.
                        T forceContribution = getExternalForceComponent(cell, iD);
                        // Then, add a contribution from the potential of all other species.
                        for (plint iPartnerSpecies=0; iPartnerSpecies<numSpecies; ++iPartnerSpecies) {
                            if (iPartnerSpecies != iSpecies) {
                                forceContribution -= speciesG[iSpecies * numSpecies + iPartnerSpecies] *
                                    rhoContribution[iPartnerSpecies][iD];
                            }
                        }
                        momentum[iD] += invOmega[iSpecies]*forceContribution;
                        // Multiply by rho to convert from velocity to momentum.
                        momentum[iD] *= *cell.getExternal(densityOffset);
                    }
                }
            }
        }
    }
}


template<typename T, template<typename U> class Descriptor>
ShanChenMultiComponentProcessor3D<T,Descriptor>*
    ShanChenMultiComponentProcessor3D<T,Descriptor>::clone() const
{
    return new ShanChenMultiComponentProcessor3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void ShanChenMultiComponentProcessor3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    // All blocks are modified by the Shan/Chen processor.
    for (pluint iBlock=0; iBlock<modified.size(); ++iBlock) {
        modified[iBlock] = modif::staticVariables;
    }
}


/* *************** ShanChenSingleComponentProcessor3D ***************** */

template<typename T, template<typename U> class Descriptor>
ShanChenSingleComponentProcessor3D<T,Descriptor>::ShanChenSingleComponentProcessor3D (
        T G_, interparticlePotential::PsiFunction<T>* Psi_ )
    : G(G_),
      Psi(Psi_)
{ }

template<typename T, template<typename U> class Descriptor>
ShanChenSingleComponentProcessor3D<T,Descriptor>::~ShanChenSingleComponentProcessor3D() {
    // Pointer to Psi function is owned; delete it in the destructor.
    delete Psi;
}

template<typename T, template<typename U> class Descriptor>
ShanChenSingleComponentProcessor3D<T,Descriptor>::ShanChenSingleComponentProcessor3D (
        ShanChenSingleComponentProcessor3D<T,Descriptor> const& rhs )
    : G(rhs.G),
      Psi(rhs.Psi->clone())
{ }

template<typename T, template<typename U> class Descriptor>
ShanChenSingleComponentProcessor3D<T,Descriptor>&
    ShanChenSingleComponentProcessor3D<T,Descriptor>::operator= (
        ShanChenSingleComponentProcessor3D<T,Descriptor> const& rhs )
{
    G = rhs.G;
    delete Psi; Psi = rhs.Psi->clone();
    return *this;
}

template<typename T, template<typename U> class Descriptor>
ShanChenSingleComponentProcessor3D<T,Descriptor>*
    ShanChenSingleComponentProcessor3D<T,Descriptor>::clone() const
{
    return new ShanChenSingleComponentProcessor3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void ShanChenSingleComponentProcessor3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    modified[0] = modif::staticVariables;
}

template<typename T, template<typename U> class Descriptor>
void ShanChenSingleComponentProcessor3D<T,Descriptor>::process (
        Box3D domain, BlockLattice3D<T,Descriptor>& lattice )
{
    // Short-hand notation for the lattice descriptor
    typedef Descriptor<T> D;
    // Handle to external scalars
    enum {
        densityOffset  = D::ExternalField::densityBeginsAt,
        momentumOffset = D::ExternalField::momentumBeginsAt
    };

    plint nx = domain.getNx() + 2;  // Include a one-cell boundary
    plint ny = domain.getNy() + 2;  // Include a one-cell boundary
    plint nz = domain.getNz() + 2;  // Include a one-cell boundary
    plint offsetX = domain.x0-1;
    plint offsetY = domain.y0-1;
    plint offsetZ = domain.z0-1;
    ScalarField3D<T> psiField(nx,ny,nz);
    
    // Compute density and momentum on every site and store result in external scalars;
    //   furthermore, evaluate the interaction potential Psi and store it into a ScalarField.
    //   Envelope cells are included, because they are needed to compute the interaction potential
    //   in the following. Note that the value of the momentum is stored temporarily only, as
    //   it is corrected later on to include corrections due to the interaction potential.
    for (plint iX=domain.x0-1; iX<=domain.x1+1; ++iX) {
        for (plint iY=domain.y0-1; iY<=domain.y1+1; ++iY) {
            for (plint iZ=domain.z0-1; iZ<=domain.z1+1; ++iZ) {
                // Get "intelligent" value of density through cell object, to account
                //   for the fact that the density value can be user-defined, for example
                //   on boundaries.
                Cell<T,Descriptor>& cell = lattice.get(iX,iY,iZ);
                T rho = cell.computeDensity();
                // Evaluate potential function psi.
                psiField.get(iX-offsetX, iY-offsetY, iZ-offsetZ) = Psi->compute(rho);
                // Store density into the corresponding external scalar.
                *cell.getExternal(densityOffset) = rho;
                // Compute momentum through direct access to particle populations, and store
                //   result in corresponding external scalars. Note that Cell::computeVelocity
                //   cannot be used, because it returns the velocity of the external scalars,
                //   not the velocity computed from the particle populations.
                Array<T,Descriptor<T>::d> j;
                momentTemplates<T,Descriptor>::get_j(cell,j);
                for (int iD=0; iD<Descriptor<T>::d; ++iD) {
                    *(cell.getExternal(momentumOffset)+iD) = j[iD];
                }
            }
        }
    }

    // Compute the interparticle forces, and store they by means of a 
    //   velocity correction in the external velocity field.
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                Array<T,D::d> rhoContribution;
                rhoContribution.resetToZero();
                // Compute the term \sum_i ( t_i psi(x+c_i,t) c_i )
                for (plint iPop = 0; iPop < D::q; ++iPop) {
                    plint nextX = iX + D::c[iPop][0];
                    plint nextY = iY + D::c[iPop][1];
                    plint nextZ = iZ + D::c[iPop][2];
                    T psi = psiField.get(nextX-offsetX, nextY-offsetY, nextZ-offsetZ);
                    for (int iD = 0; iD < D::d; ++iD) {
                        rhoContribution[iD] += D::t[iPop] * psi * D::c[iPop][iD];
                    }
                }

                // Computation and storage of the final momentum, including tho momentum
                //   difference due to interaction potential and the external force.
                Cell<T,Descriptor>& cell = lattice.get(iX,iY,iZ);
                T *momentum = cell.getExternal(momentumOffset);
                for (int iD = 0; iD < D::d; ++iD) {
                    // Initialize force contribution with force from external fields if there
                    //   is any, or with zero otherwise.
                    T forceContribution = getExternalForceComponent(cell, iD);
                    // Add interaction term.
                    T psi = psiField.get(iX-offsetX, iY-offsetY, iZ-offsetZ);
                    forceContribution -= G * psi * rhoContribution[iD];
                    // Include into total momentum.
                    momentum[iD] += (T)1/cell.getDynamics().getOmega()*forceContribution;
                }
            }
        }
    }
}

/* *************** ShanChenExternalMultiComponentProcessor3D ***************** */

template<typename T, template<typename U> class Descriptor>
ShanChenExternalMultiComponentProcessor3D <T,Descriptor>::ShanChenExternalMultiComponentProcessor3D(T G_)
    : G(G_)
{ }

template<typename T, template<typename U> class Descriptor>
ShanChenExternalMultiComponentProcessor3D <T,Descriptor>::ShanChenExternalMultiComponentProcessor3D(
        std::vector<std::vector<T> > const& speciesG_)
    : G((T) 0)
{
    pluint numSpecies = speciesG_.size();
    // Although speciesG_ has a 2D "matrix structure", speciesG has a 1D "array structure".
    speciesG.resize(numSpecies * numSpecies);
    for (pluint iSpecies = 0; iSpecies < numSpecies; iSpecies++) {
        PLB_ASSERT(speciesG_[iSpecies].size() == numSpecies);
        for (pluint jSpecies = 0; jSpecies < numSpecies; jSpecies++) {
            speciesG[iSpecies * numSpecies + jSpecies] = speciesG_[iSpecies][jSpecies];
        }
    }
}

template<typename T, template<typename U> class Descriptor>
ShanChenExternalMultiComponentProcessor3D <T,Descriptor>::ShanChenExternalMultiComponentProcessor3D (
        T G_, std::vector<T> const& imposedOmega_)
    : G(G_),
      imposedOmega(imposedOmega_)
{ }

template<typename T, template<typename U> class Descriptor>
ShanChenExternalMultiComponentProcessor3D <T,Descriptor>::ShanChenExternalMultiComponentProcessor3D (
        std::vector<std::vector<T> > const& speciesG_, std::vector<T> const& imposedOmega_)
    : G((T) 0),
      imposedOmega(imposedOmega_)
{
    pluint numSpecies = speciesG_.size();
    // Although speciesG_ has a 2D "matrix structure", speciesG has a 1D "array structure".
    speciesG.resize(numSpecies * numSpecies);
    for (pluint iSpecies = 0; iSpecies < numSpecies; iSpecies++) {
        PLB_ASSERT(speciesG_[iSpecies].size() == numSpecies);
        for (pluint jSpecies = 0; jSpecies < numSpecies; jSpecies++) {
            speciesG[iSpecies * numSpecies + jSpecies] = speciesG_[iSpecies][jSpecies];
        }
    }
}

template<typename T, template<typename U> class Descriptor>
void ShanChenExternalMultiComponentProcessor3D<T,Descriptor>::processGenericBlocks (
        Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks )
{
    // Short-hand notation for the lattice descriptor
    typedef Descriptor<T> D;

    PLB_ASSERT( atomicBlocks.size() % 3 == 0 );
    std::vector<BlockLattice3D<T,Descriptor>*> lattices;
    std::vector<ScalarField3D<T>*> rhoBar;
    std::vector<TensorField3D<T,D::d>*> j;
    std::vector<Dot3D> ofsLatt;
    std::vector<Dot3D> ofsRhoBar;
    std::vector<Dot3D> ofsJ;
    for (pluint i=0; i<atomicBlocks.size(); i+=3) {
        BlockLattice3D<T,Descriptor>* lattice_ = dynamic_cast<BlockLattice3D<T,Descriptor>*>(atomicBlocks[i]);
        PLB_ASSERT( lattice_ );
        lattices.push_back(lattice_);
        ofsLatt.push_back(computeRelativeDisplacement(*lattices[0], *lattice_));

        ScalarField3D<T>* rhoBar_ = dynamic_cast<ScalarField3D<T>*>(atomicBlocks[i+1]);
        PLB_ASSERT( rhoBar_ );
        rhoBar.push_back(rhoBar_);
        ofsRhoBar.push_back(computeRelativeDisplacement(*lattices[0], *rhoBar_));

        TensorField3D<T,D::d>* j_ = dynamic_cast<TensorField3D<T,D::d>*>(atomicBlocks[i+2]);
        PLB_ASSERT( j_ );
        j.push_back(j_);
        ofsJ.push_back(computeRelativeDisplacement(*lattices[0], *j_));

    }
    // Number of species (or components) which are coupled in this Shan/Chen multi-component fluid.
    plint numSpecies = (plint) lattices.size();
    
    // Compute per-lattice density  and momentum on every site and on each
    //   lattice, and store result in external blocks; envelope cells are included,
    //   because they are needed to compute the interaction potential in the following.
    //   Note that the per-lattice value of the momentum is stored temporarily only, as
    //   it is corrected later on, based on the common fluid velocity.
    for (plint iSpecies=0; iSpecies<numSpecies; ++iSpecies) {
        for (plint iX=domain.x0-1; iX<=domain.x1+1; ++iX) {
            for (plint iY=domain.y0-1; iY<=domain.y1+1; ++iY) {
                for (plint iZ=domain.z0-1; iZ<=domain.z1+1; ++iZ) {
                    // Get "intelligent" value of density through cell object, to account
                    //   for the fact that the density value can be user-defined, for example
                    //   on boundaries.
                    Cell<T,Descriptor>& cell =
                        lattices[iSpecies]->get (
                                iX+ofsLatt[iSpecies].x,iY+ofsLatt[iSpecies].y,iZ+ofsLatt[iSpecies].z );
                    //if (!cell.getDynamics().hasMoments()) {
                    //    break;
                    //}
                    cell.getDynamics().computeRhoBarJ (
                            cell, rhoBar[iSpecies]->get (
                                      iX+ofsRhoBar[iSpecies].x,iY+ofsRhoBar[iSpecies].y,iZ+ofsRhoBar[iSpecies].z ),
                                  j[iSpecies]->get (
                                      iX+ofsJ[iSpecies].x,iY+ofsJ[iSpecies].y,iZ+ofsJ[iSpecies].z) );
                    //momentTemplates<T,Descriptor>::get_j (
                    //        cell, j[iSpecies]->get (
                    //                  iX+ofsJ[iSpecies].x,iY+ofsJ[iSpecies].y,iZ+ofsJ[iSpecies].z) );
                }
            }
        }
    }

    // Temporary variable for the relaxation parameters omega.
    std::vector<T> omega(numSpecies), invOmega(numSpecies);
    // Temporary variable for total velocity.
    Array<T,D::d> uTot;
    // Temporary variable for interaction potential.
    std::vector<Array<T,D::d> > rhoContribution(numSpecies);

    // If omega is constant and imposed by the user, copy its value to
    //   the vector "omega", and compute the inverse.
    if (!imposedOmega.empty()) {
        PLB_ASSERT((plint) imposedOmega.size() == numSpecies);
        omega = imposedOmega;
        for (pluint iOmega=0; iOmega<omega.size(); ++iOmega) {
            invOmega[iOmega] = (T)1 / omega[iOmega];
        }
    }

    if (speciesG.empty()) {
        speciesG.resize(numSpecies * numSpecies, G);
    }

    // Compute the interaction force between the species, and store it by
    //   means of a velocity correction in the external velocity field.
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
            for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
                // Computation of the common density over all populations, weighted by
                //   the relaxation parameters omega.
                T weightedDensity = T();
                // A cell is only invalid if all species have bounce-back.
                bool validCell = false;
                for (plint iSpecies=0; iSpecies<numSpecies; ++iSpecies) {
                    Cell<T,Descriptor> const& cell =
                        lattices[iSpecies]->get(iX+ofsLatt[iSpecies].x,iY+ofsLatt[iSpecies].y,iZ+ofsLatt[iSpecies].z);
                    if (!cell.getDynamics().hasMoments()) {
                        break;
                    }
                    else {
                        validCell = true;
                    }
                    // Take this opportunity to read omega from the cell, unless the value
                    //   of omega is constant and imposed by the user.
                    if (imposedOmega.empty()) {
                        omega[iSpecies] = cell.getDynamics().getOmega();
                        invOmega[iSpecies] = (T)1/omega[iSpecies];
                    }
                    weightedDensity += omega[iSpecies] * D::fullRho(rhoBar[iSpecies]->get(
                                iX+ofsRhoBar[iSpecies].x,iY+ofsRhoBar[iSpecies].y,iZ+ofsRhoBar[iSpecies].z));
                }
                // Computation of the common velocity, shared among all populations.
                uTot.resetToZero();
                if (validCell) {
                    for (int iD = 0; iD < D::d; ++iD) {
                        for (plint iSpecies=0; iSpecies<numSpecies; ++iSpecies) {
                            Array<T,D::d> momentum = j[iSpecies]->get(iX+ofsJ[iSpecies].x,iY+ofsJ[iSpecies].y,iZ+ofsJ[iSpecies].z);
                            uTot[iD] += momentum[iD] * omega[iSpecies];
                        }
                        uTot[iD] /= weightedDensity;
                    }
                }

                // Computation of the interaction potential.
                for (plint iSpecies=0; iSpecies<numSpecies; ++iSpecies) {
                    Cell<T,Descriptor> const& cell =
                        lattices[iSpecies]->get(iX+ofsLatt[iSpecies].x,iY+ofsLatt[iSpecies].y,iZ+ofsLatt[iSpecies].z);
                    if (!cell.getDynamics().hasMoments()) {
                        break;
                    }
                    multiPhaseTemplates3D<T,Descriptor>::shanChenInteraction (
                            *lattices[iSpecies], *rhoBar[iSpecies], rhoContribution[iSpecies],
                            iX+ofsLatt[iSpecies].x,iY+ofsLatt[iSpecies].y,iZ+ofsLatt[iSpecies].z );
                }

                // Computation and storage of the final velocity, consisting
                //   of uTot plus the momentum difference due to interaction
                //   potential and external force
                for (plint iSpecies=0; iSpecies<numSpecies; ++iSpecies) {
                    Cell<T,Descriptor>& cell =
                        lattices[iSpecies]->get (
                            iX+ofsLatt[iSpecies].x,iY+ofsLatt[iSpecies].y,iZ+ofsLatt[iSpecies].z );
                    if (!cell.getDynamics().hasMoments()) {
                        break;
                    }
                    Array<T,D::d>& momentum = j[iSpecies]->get(iX+ofsJ[iSpecies].x,iY+ofsJ[iSpecies].y,iZ+ofsJ[iSpecies].z);
                    for (int iD = 0; iD < D::d; ++iD) {
                        momentum[iD] = uTot[iD];
                        // Initialize force contribution with force from external fields if there
                        //   is any, or with zero otherwise.
                        T forceContribution = getExternalForceComponent(cell, iD);
                        // Then, add a contribution from the potential of all other species.
                        for (plint iPartnerSpecies=0; iPartnerSpecies<numSpecies; ++iPartnerSpecies) {
                            if (iPartnerSpecies != iSpecies) {
                                forceContribution -= speciesG[iSpecies * numSpecies + iPartnerSpecies] *
                                    rhoContribution[iPartnerSpecies][iD];
                            }
                        }
                        momentum[iD] += invOmega[iSpecies]*forceContribution;
                        // Multiply by rho to convert from velocity to momentum.
                        momentum[iD] *= D::fullRho(rhoBar[iSpecies]->get (
                                            iX+ofsRhoBar[iSpecies].x,iY+ofsRhoBar[iSpecies].y,iZ+ofsRhoBar[iSpecies].z) );
                    }
                }
            }
        }
    }
}


template<typename T, template<typename U> class Descriptor>
ShanChenExternalMultiComponentProcessor3D<T,Descriptor>*
    ShanChenExternalMultiComponentProcessor3D<T,Descriptor>::clone() const
{
    return new ShanChenExternalMultiComponentProcessor3D<T,Descriptor>(*this);
}

template<typename T, template<typename U> class Descriptor>
void ShanChenExternalMultiComponentProcessor3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
{
    // All blocks are modified by the Shan/Chen processor.
    for (pluint iBlock=0; iBlock<modified.size(); ++iBlock) {
        modified[iBlock] = modif::staticVariables;
    }
}


}  // namespace plb

#endif  // SHAN_CHEN_PROCESSOR_3D_HH
