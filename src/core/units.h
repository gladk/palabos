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

#ifndef UNITS_H
#define UNITS_H

#include "core/globalDefs.h"
#include "core/globalDefs.h"
#include "parallelism/mpiManager.h"
#include "io/parallelIO.h"
#include <string>
#include <fstream>

namespace plb {

/// Numeric parameters for isothermal, incompressible flow.
template<typename T>
class IncomprFlowParam {
public:
    /// Constructor
    /** \param latticeU_  Characteristic velocity in lattice units (proportional to Mach number).
     *  \param Re_ Reynolds number.
     *  \param N_  Resolution (a lattice of size 1 has N_+1 cells).
     *  \param lx_ x-length in dimensionless units (e.g. 1).
     *  \param ly_ y-length in dimensionless units (e.g. 1).
     *  \param lz_ z-length in dimensionless units (e.g. 1).
     */
    IncomprFlowParam(T physicalU_, T latticeU_, T Re_, T physicalLength_, plint resolution_, T lx_, T ly_, T lz_=T() )
        : physicalU(physicalU_), latticeU(latticeU_), physicalLength(physicalLength_), Re(Re_),
          resolution(resolution_), lx(lx_), ly(ly_), lz(lz_)
    { }
    
    IncomprFlowParam(T latticeU_, T Re_, plint resolution_, T lx_, T ly_, T lz_=T() )
    : latticeU(latticeU_), Re(Re_), resolution(resolution_), lx(lx_), ly(ly_), lz(lz_)
    { 
        physicalU      = (T)1;
        physicalLength = (T)1; 
    }
    /// velocity in lattice units (proportional to Mach number)
    T getLatticeU() const { return latticeU; }
    /// velocity in physical units
    T getPhysicalU() const { return physicalU; }
    /// Reynolds number
    T getRe() const      { return Re; }
    /// physical resolution
    T getPhysicalLength() const { return physicalLength; }
    /// resolution
    plint getResolution() const { return resolution; }
    /// x-length in dimensionless units
    T getLx() const      { return getPhysicalLength()*lx; }
    /// y-length in dimensionless units
    T getLy() const      { return getPhysicalLength()*ly; }
    /// z-length in dimensionless units
    T getLz() const      { return getPhysicalLength()*lz; }
    /// lattice spacing in dimensionless units
    T getDeltaX() const  { return (T)getPhysicalLength()/(T)getResolution(); }
    /// time step in dimensionless units
    T getDeltaT() const  { return getDeltaX()*getLatticeU()/getPhysicalU(); }
    /// conversion from dimensionless to lattice units for space coordinate
    plint nCell(T l) const { return (int)(l/getDeltaX()+(T)0.5); }
    /// conversion from dimensionless to lattice units for time coordinate
    plint nStep(T t) const { return (int)(t/getDeltaT()+(T)0.5); }
    /// number of lattice cells in x-direction
    plint getNx(bool offLattice=false) const { return nCell(getLx())+1+(int)offLattice; }
    /// number of lattice cells in y-direction
    plint getNy(bool offLattice=false) const { return nCell(getLy())+1+(int)offLattice; }
    /// number of lattice cells in z-direction
    plint getNz(bool offLattice=false) const { return nCell(getLz())+1+(int)offLattice; }
    /// viscosity in lattice units
    T getLatticeNu() const { return getLatticeU()*(T)getResolution()/Re; }
    /// relaxation time
    T getTau() const       { return (T)3*getLatticeNu()+(T)0.5; }
    /// relaxation frequency
    T getOmega() const     { return (T)1 / getTau(); }
private:
    T physicalU, latticeU, physicalLength, Re;
    plint resolution;
    T lx, ly, lz;
};

template<typename T>
void writeLogFile(IncomprFlowParam<T> const& parameters,
                  std::string const& title)
{
    std::string fullName = global::directories().getLogOutDir() + "plbLog.dat";
    plb_ofstream ofile(fullName.c_str());
    ofile << title << "\n\n";
    ofile << "Velocity in lattice units: u=" << parameters.getLatticeU() << "\n";
    ofile << "Reynolds number:           Re=" << parameters.getRe() << "\n";
    ofile << "Lattice resolution:        N=" << parameters.getResolution() << "\n";
    ofile << "Relaxation frequency:      omega=" << parameters.getOmega() << "\n";
    ofile << "Extent of the system:      lx=" << parameters.getLx() << "\n";
    ofile << "Extent of the system:      ly=" << parameters.getLy() << "\n";
    ofile << "Extent of the system:      lz=" << parameters.getLz() << "\n";
    ofile << "Grid spacing deltaX:       dx=" << parameters.getDeltaX() << "\n";
    ofile << "Time step deltaT:          dt=" << parameters.getDeltaT() << "\n";
}

/// Numeric parameters for isothermal, incompressible flow.
template<typename T>
class ComprFlowParam {
public:
    /// Constructor
    /** \param latticeU_  Characteristic velocity in lattice units (proportional to Mach number).
     *  \param Re_ Reynolds number.
     *  \param Pe_ Peclet number.
     *  \param N_  Resolution (a lattice of size 1 has N_+1 cells).
     *  \param lx_ x-length in dimensionless units (e.g. 1).
     *  \param ly_ y-length in dimensionless units (e.g. 1).
     *  \param lz_ z-length in dimensionless units (e.g. 1).
     */
    ComprFlowParam(T latticeU_, T latticeRho_, T latticeTemp_, T physU_, T physRho_, T physTemp_, 
                   T Re_, T Pe_, plint resolution_, T lx_, T ly_, T lz_=T() )
        : latticeU(latticeU_), latticeRho(latticeRho_), latticeTemp(latticeTemp_), 
          physRho(physRho_), physU(physU_), physTemp(physTemp_), 
          Re(Re_), Pe(Pe_), resolution(resolution_), lx(lx_), ly(ly_), lz(lz_)
    { }
    /// velocity in lattice units 
    T getLatticeU() const { return latticeU; }
    /// density in lattice units 
    T getLatticeRho() const { return latticeRho; }
    /// temperature in lattice units 
    T getLatticeTemp() const { return latticeTemp; }
    /// velocity in physical units 
    T getPhysicalU() const { return physU; }
    /// density in physical units 
    T getPhysicalRho() const { return physRho; }
    /// temperature in physical units 
    T getPhysicalTemp() const { return physTemp; }
    /// Reynolds number
    T getRe() const      { return Re; }
    /// Peclet number (= Re*Pr)
    T getPe() const      { return Pe; }
    /// resolution
    plint getResolution() const { return resolution; }
    /// x-length in dimensionless units
    T getLx() const      { return lx; }
    /// y-length in dimensionless units
    T getLy() const      { return ly; }
    /// z-length in dimensionless units
    T getLz() const      { return lz; }
    /// lattice spacing in dimensionless units
    T getDeltaX() const  { return (T)1/(T)getResolution(); }
    /// time step in dimensionless units
    T getDeltaT() const  { return getDeltaX() * getLatticeU() / getPhysicalU(); }
    /// density in dimensionless units
    T getDeltaRho() const  { return getPhysicalRho() / getLatticeRho(); }
    /// density in dimensionless units
    T getDeltaTemp() const  { return getPhysicalTemp() / getLatticeTemp(); }
    /// conversion from dimensionless to lattice units for space coordinate
    plint nCell(T l) const { return (int)(l/getDeltaX()+(T)0.5); }
    /// conversion from dimensionless to lattice units for time coordinate
    plint nStep(T t) const { return (int)(t/getDeltaT()+(T)0.5); }
    /// number of lattice cells in x-direction
    plint getNx(bool offLattice=false) const { return nCell(lx)+1+(int)offLattice; }
    /// number of lattice cells in y-direction
    plint getNy(bool offLattice=false) const { return nCell(ly)+1+(int)offLattice; }
    /// number of lattice cells in z-direction
    plint getNz(bool offLattice=false) const { return nCell(lz)+1+(int)offLattice; }
    /// viscosity in lattice units
    T getLatticeMu() const { return getLatticeU()*getResolution()*getLatticeRho() / Re; }
    /// viscosity in lattice units 
    T getLatticeKappa() const { return getLatticeU()*getResolution()*getLatticeRho() / Pe; }
    /// relaxation time
    T getTau() const       { return getLatticeMu() / getLatticeRho() / getLatticeTemp()+(T)0.5; }
    // TODO for the moment only Pr = 1 (Pe = Re) fluids are simulable....
    /// relaxation frequency
    T getOmega() const     { return (T)1 / getTau(); }
private:
    T latticeU, latticeRho, latticeTemp, physRho, physU, physTemp, Re, Pe;
    plint resolution;
    T lx, ly, lz;
};

template<typename T>
void writeLogFile(ComprFlowParam<T> const& parameters,
                  std::string const& title)
{
    std::string fullName = global::directories().getLogOutDir() + "plbLog.dat";
    plb_ofstream ofile(fullName.c_str());
    ofile << title << "\n\n";
    ofile << "Velocity in lattice units:    u=" << parameters.getLatticeU() << "\n";
    ofile << "Density in lattice units:     rho=" << parameters.getLatticeRho() << "\n";
    ofile << "Temperature in lattice units: T=" << parameters.getLatticeTemp() << "\n";
    ofile << "Reynolds number:              Re=" << parameters.getRe() << "\n";
    ofile << "Peclet number:                Pe=" << parameters.getPe() << "\n";
    ofile << "Lattice resolution:           N=" << parameters.getResolution() << "\n";
    ofile << "Extent of the system:         lx=" << parameters.getLx() << "\n";
    ofile << "Extent of the system:         ly=" << parameters.getLy() << "\n";
    ofile << "Extent of the system:         lz=" << parameters.getLz() << "\n";
    ofile << "Grid spacing deltaX:          dx=" << parameters.getDeltaX() << "\n";
    ofile << "Time step deltaT:             dt=" << parameters.getDeltaT() << "\n";
    ofile << "Density deltaRho:             dRho=" << parameters.getDeltaRho() << "\n";
    ofile << "Temp deltaTemp:               dTemp=" << parameters.getDeltaTemp() << "\n";
    ofile << "Relaxation time:              Tau=" << parameters.getTau() << "\n";
}

}  // namespace plb

#endif

