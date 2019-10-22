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

#ifndef UNITS_H
#define UNITS_H

#include "core/globalDefs.h"
#include "core/array.h"
#include "core/geometry2D.h"
#include "core/geometry3D.h"
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
        : physicalU(physicalU_), latticeU(latticeU_), Re(Re_), physicalLength(physicalLength_),
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
    T physicalU, latticeU, Re, physicalLength;
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

class Units3D {
public:
    Units3D();
    /// This constructor explicitly sets the physical domain. The
    /// resolution isn't specified here. You should specify it later on,
    /// for example with set_dx() or set_xResolution().
    Units3D(Cuboid<double> const& physDomain);
    /// This constructor explicitly sets the physical domain. The
    /// resolution isn't specified here. You should specify it later on,
    /// for example with set_dx() or set_xResolution().
    Units3D(Array<double,3> const& lowerLeftCorner, Array<double,3> const& upperRightCorner);
    /// This constructor sets the physical domain automatically, using
    /// the origin as lower-left corner.
    Units3D(plint nx, plint ny, plint nz, double dx_);
    /// Sets the discrete spacing dx, adjusts the physical domain
    /// so that all lengths are multiples of dx, and recomputes
    /// the LB domain.
    void set_dx(double dx__);
    void set_dt(double dt__);
    /// This is an indirect way to fix the discrete time step,
    /// by imposing a lattice velocity.
    void set_uLB(double uPhys, double uLB);
    /// This is an indirect way to fix the discrete time step,
    /// by imposing the relaxation time.
    void set_nuLB(double nuPhys, double tau);
    /// Recomputes the LB domain according to the indicated
    /// resolution. There will be nx+1 nodes in x-direction. The
    /// other directions follow proportionally.
    void set_xResolution(plint nx);
    /// Recomputes the LB domain according to the indicated
    /// resolution. There will be ny+1 nodes in y-direction. The
    /// other directions follow proportionally.
    void set_yResolution(plint ny);
    /// Recomputes the LB domain according to the indicated
    /// resolution. There will be nz+1 nodes in z-direction. The
    /// other directions follow proportionally.
    void set_zResolution(plint nz);
    /// Recomputes the LB domain according to the indicated
    /// resolution. There will be n+1 nodes in the indicated
    /// direction. The other directions follow proportionally.
    /// The direction must be 0, 1, or 2 for x, y, and z.
    void setResolution(plint n, plint direction);
    /// Define is a given space direction is periodic or not.
    /// By default, no direction is periodic.
    /// In space directions which switch from non-periodic to periodic,
    /// the computational domain is reduced by one node at the rightmost
    /// bound (both the physical and the LB domain). In space directions
    /// which switch from periodic to non-periodic, the computational
    /// domain is enlarged by one node.
    void set_periodic(bool periodicX, bool periodicY, bool periodicZ);
    /// Enlarge the domain (both the physical and the LB one) by "delta"
    /// lattice nodes in each (positive and negative) direction. As a
    /// result, the total domain is 2*delta bigger than initially in
    /// each direction.
    void extend(plint delta);
    /// Enlarge the domain (both the physical and the LB one) by the number
    /// of indicated nodes in each (positive and negative) direction. In 
    /// x-direction for example, the domain gets bigger by dx1+dx0 nodes.
    void extend(plint dx0, plint dx1, plint dy0, plint dy1, plint dz0, plint dz1);
public:
    /// Computational domain in physical units.
    Cuboid<double> physDomain() const;
    /// Value of lower-left corner, corresponding to the origin in lattice units.
    Array<double,3> physOffset() const;
    /// Computational domain in lattice units (lower-left corner is always (0,0,0)).
    Box3D lbDomain() const;
    /// Returns the resolution in x-direction. If the x-direction is periodic,
    /// this is equal to the number of nodes of the computational domain in
    /// x-direction. If it is non-periodic, it is equal to the number of nodes
    /// minus 1.
    plint nx() const;
    /// Returns the resolution in x-direction. If the x-direction is periodic,
    /// this is equal to the number of nodes of the computational domain in
    /// x-direction. If it is non-periodic, it is equal to the number of nodes
    /// minus 1.
    plint ny() const;
    /// Returns the resolution in x-direction. If the x-direction is periodic,
    /// this is equal to the number of nodes of the computational domain in
    /// x-direction. If it is non-periodic, it is equal to the number of nodes
    /// minus 1.
    plint nz() const;
    /// Returns the discrete grid spacing.
    double dx() const;
    /// Returns the discrete time step.
    double dt() const;
    /// Converts a physical time to LB iteration (rounds to the nearest integer).
    plint lbIter(double physTime) const;
    /// Converts a LB iteration to physical time.
    double physTime(plint lbIter) const;
    /// Converts a physical length to lattice units (rounds to the nearest integer).
    plint numCells(double physLength) const;
    /// Converts a length from lattice to physical units.
    double physLength(plint numCells) const;
    /// Converts a velocity from physical to lattice units.
    double lbVel(double physVel) const;
    /// Converts a velocity from lattice to physical units.
    double physVel(double lbVel) const;
    /// Converts a velocity from physical to lattice units.
    Array<double,3> lbVel(Array<double,3> const& physVel) const;
    /// Converts a velocity from lattice to physical units.
    Array<double,3> physVel(Array<double,3> const& lbVel) const;
    /// Converts an acceleration from physical to lattice units.
    double lbAcceleration(double physAcceleration) const;
    /// Converts an acceleration from lattice to physical units.
    double physAcceleration(double lbAcceleration) const;
    /// Converts an acceleration from physical to lattice units.
    Array<double,3> lbAcceleration(Array<double,3> const& physAcceleration) const;
    /// Converts an acceleration from lattice to physical units.
    Array<double,3> physAcceleration(Array<double,3> const& lbAcceleration) const;
    /// Converts a kinematic viscosity from physical to lattice units.
    double lbVisc(double physVisc) const;
    /// Converts a kinematic viscosity from lattice to physical units.
    double physVisc(double lbVisc) const;
    /// Converts a physical-units viscosity to relaxation time.
    double tau(double physVisc, double cs2 = 1./3.) const;
    /// Converts a physical-units viscosity to relaxation frequency.
    double omega(double physVisc, double cs2 = 1./3.) const;
    /// Converts a pressure from physical to lattice units. The parameter
    /// rho0 is the fluid density in physical units.
    double lbPressure(double physPressure, double rho0=1.0) const;
    /// Converts a pressure from lattice to physical units. The parameter
    /// rho0 is the fluid density in physical units.
    double physPressure(double lbPressure, double rho0=1.0) const;
    /// Converts a scalar force from physical to lattice units. The parameter
    /// rho0 is the fluid density in physical units.
    double lbForce(double physForce, double rho0=1.0) const;
    /// Converts a scalar force from lattice to physical units. The parameter
    /// rho0 is the fluid density in physical units.
    double physForce(double lbForce, double rho0=1.0) const;
    /// Converts a vector force from physical to lattice units. The parameter
    /// rho0 is the fluid density in physical units.
    Array<double,3> lbForce(Array<double,3> const& physForce, double rho0=1.0) const;
    /// Converts a vector force from lattice to physical units. The parameter
    /// rho0 is the fluid density in physical units.
    Array<double,3> physForce(Array<double,3> const& lbForce, double rho0=1.0) const;
    /// Converts a scalar torque from physical to lattice units. The parameter
    /// rho0 is the fluid density in physical units.
    double lbTorque(double physTorque, double rho0=1.0) const;
    /// Converts a scalar torque from lattice to physical units. The parameter
    /// rho0 is the fluid density in physical units.
    double physTorque(double lbTorque, double rho0=1.0) const;
    /// Converts a vector torque from physical to lattice units. The parameter
    /// rho0 is the fluid density in physical units.
    Array<double,3> lbTorque(Array<double,3> const& physTorque, double rho0=1.0) const;
    /// Converts a vector torque from lattice to physical units. The parameter
    /// rho0 is the fluid density in physical units.
    Array<double,3> physTorque(Array<double,3> const& lbTorque, double rho0=1.0) const;
    /// Converts a vector coordinate from physical to lattice units. No rounding
    /// involved.
    Array<double,3> lbCoord(Array<double,3> const& physCoord) const;
    /// Converts a vector coordinate from lattice to physical units. No rounding
    /// involved.
    Array<double,3> physCoord(Array<double,3> const& lbCoord) const;
private:
    void computePhysDomain();
    plint adjustLength(double x0, double& x1);
private:
    Box3D lbDomain_;
    Cuboid<double> physDomain_;
    double dx_, dt_;
    bool periodicX_, periodicY_, periodicZ_;
};

}  // namespace plb

#endif  // UNITS_H

