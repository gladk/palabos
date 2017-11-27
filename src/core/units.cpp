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
 * Units3D -- source code.
 */
#include "core/units.h"

namespace plb {

Units3D::Units3D()
    : lbDomain_(0,100, 0,100, 0,100),
      dx_(1.), dt_(1.),
      periodicX_(false),
      periodicY_(false),
      periodicZ_(false)
{
    computePhysDomain();
}

Units3D::Units3D(Cuboid<double> const& physDomain)
    : lbDomain_(0,0, 0,0, 0,0),
      physDomain_(physDomain),
      dx_(1.),
      dt_(1.),
      periodicX_(false),
      periodicY_(false),
      periodicZ_(false)
{ }

Units3D::Units3D(Array<double,3> const& lowerLeftCorner, Array<double,3> const& upperRightCorner)
    : lbDomain_(0,0, 0,0, 0,0),
      physDomain_(lowerLeftCorner, upperRightCorner),
      dx_(1.),
      dt_(1.),
      periodicX_(false),
      periodicY_(false),
      periodicZ_(false)
{ }

Units3D::Units3D(plint nx, plint ny, plint nz, double dx)
    : lbDomain_(0, nx, 0, ny, 0, nz),
      dx_(dx),
      dt_(1.),
      periodicX_(false),
      periodicY_(false),
      periodicZ_(false)
{
    computePhysDomain();
}

Cuboid<double> Units3D::physDomain() const {
    return physDomain_;
}

Array<double,3> Units3D::physOffset() const {
    return physDomain_.lowerLeftCorner;
}

Box3D Units3D::lbDomain() const {
    return lbDomain_;
}

plint Units3D::nx() const {
    return periodicX_ ? lbDomain_.getNx() : lbDomain_.getNx()-1;
}

plint Units3D::ny() const {
    return periodicY_ ? lbDomain_.getNy() : lbDomain_.getNy()-1;
}

plint Units3D::nz() const {
    return periodicZ_ ? lbDomain_.getNz() : lbDomain_.getNz()-1;
}


double Units3D::dx() const {
    return dx_;
}

double Units3D::dt() const {
    return dt_;
}

plint Units3D::lbIter(double physTime) const {
    return util::roundToInt(physTime/dt_);
}

double Units3D::physTime(plint lbIter) const {
    return lbIter*dt_;
}

plint Units3D::numCells(double physLength) const {
    return util::roundToInt(physLength/dx_);
}

double Units3D::physLength(plint numCells) const {
    return numCells*dx_;
}

double Units3D::lbVel(double physVel) const {
    return physVel * dt_/dx_;
}

double Units3D::physVel(double lbVel) const {
    return lbVel * dx_/dt_;
}

Array<double,3> Units3D::lbVel(Array<double,3> const& physVel) const {
    return physVel * dt_/dx_;
}

Array<double,3> Units3D::physVel(Array<double,3> const& lbVel) const {
    return lbVel * dx_/dt_;
}

double Units3D::lbAcceleration(double physAcceleration) const {
    return physAcceleration * dt_*dt_/dx_;
}

double Units3D::physAcceleration(double lbAcceleration) const {
    return lbAcceleration * dx_/(dt_*dt_);
}

Array<double,3> Units3D::lbAcceleration(Array<double,3> const& physAcceleration) const {
    return physAcceleration * dt_*dt_/dx_;
}

Array<double,3> Units3D::physAcceleration(Array<double,3> const& lbAcceleration) const {
    return lbAcceleration * dx_/(dt_*dt_);
}

double Units3D::lbVisc(double physVisc) const {
    return physVisc * dt_/(dx_*dx_);
}

double Units3D::physVisc(double lbVisc) const {
    return lbVisc * dx_*dx_/dt_;
}

double Units3D::tau(double physVisc, double cs2) const {
    return 0.5 + lbVisc(physVisc)/cs2;
}

double Units3D::omega(double physVisc, double cs2) const {
    return 1. / tau(physVisc, cs2);
}

double Units3D::lbPressure(double physPressure, double rho0) const {
    return physPressure / rho0 * util::sqr(dt_) / util::sqr(dx_);
}

double Units3D::physPressure(double lbPressure, double rho0) const {
    return rho0 * util::sqr(dx_)/util::sqr(dt_) * lbPressure;
}

double Units3D::lbForce(double physForce, double rho0) const {
    return physForce / rho0 * util::sqr(dt_) / util::sqr(dx_*dx_);
}

double Units3D::physForce(double lbForce, double rho0) const {
    return rho0 * util::sqr(dx_*dx_)/util::sqr(dt_) * lbForce;
}

Array<double,3> Units3D::lbForce(Array<double,3> const& physForce, double rho0) const {
    return Array<double,3> (
            lbForce(physForce[0], rho0), 
            lbForce(physForce[1], rho0), 
            lbForce(physForce[2], rho0) );
}

Array<double,3> Units3D::physForce(Array<double,3> const& lbForce, double rho0) const {
    return Array<double,3> (
            physForce(lbForce[0], rho0), 
            physForce(lbForce[1], rho0), 
            physForce(lbForce[2], rho0) );
}

double Units3D::lbTorque(double physTorque, double rho0) const {
    return physTorque / rho0 * util::sqr(dt_) / (util::sqr(dx_*dx_)*dx_);
}

double Units3D::physTorque(double lbTorque, double rho0) const {
    return rho0 * (util::sqr(dx_*dx_)*dx_)/util::sqr(dt_) * lbTorque;
}

Array<double,3> Units3D::lbTorque(Array<double,3> const& physTorque, double rho0) const {
    return Array<double,3> (
            lbTorque(physTorque[0], rho0), 
            lbTorque(physTorque[1], rho0), 
            lbTorque(physTorque[2], rho0) );
}

Array<double,3> Units3D::physTorque(Array<double,3> const& lbTorque, double rho0) const {
    return Array<double,3> (
            physTorque(lbTorque[0], rho0), 
            physTorque(lbTorque[1], rho0), 
            physTorque(lbTorque[2], rho0) );
}

Array<double,3> Units3D::lbCoord(Array<double,3> const& physCoord) const {
    return (physCoord-physDomain_.lowerLeftCorner)/dx_;
}

Array<double,3> Units3D::physCoord(Array<double,3> const& lbCoord) const {
    return physDomain_.lowerLeftCorner + lbCoord*dx_;
}


plint Units3D::adjustLength(double x0, double& x1)
{
    double l = x1-x0;
    plint nAve = (plint)(l/dx_+0.5);
    x1 = x0 + nAve*dx_;
    return nAve;
}

void Units3D::set_dx(double dx__) {
    dx_ = dx__;
    plint nx = adjustLength( physDomain_.lowerLeftCorner[0],
                             physDomain_.upperRightCorner[0] );
    plint ny = adjustLength( physDomain_.lowerLeftCorner[1],
                             physDomain_.upperRightCorner[1] );
    plint nz = adjustLength( physDomain_.lowerLeftCorner[2],
                             physDomain_.upperRightCorner[2] );
    if (periodicX_ && nx>0) {
        --nx;
    }
    if (periodicY_ && ny>0) {
        --ny;
    }
    if (periodicZ_ && nz>0) {
        --nz;
    }
    lbDomain_ = Box3D(0, nx, 0, ny, 0, nz);
}

void Units3D::set_dt(double dt__) {
    dt_ = dt__;
}

void Units3D::set_uLB(double uPhys, double uLB) {
    dt_ = dx_ * uLB/uPhys;
}

void Units3D::set_nuLB(double nuPhys, double nuLB) {
    dt_ = dx_*dx_ * nuLB/nuPhys;
}

void Units3D::set_xResolution(plint nx) {
    setResolution(nx, 0);
}

void Units3D::set_yResolution(plint ny) {
    setResolution(ny, 1);
}

void Units3D::set_zResolution(plint nz) {
    setResolution(nz, 2);
}

void Units3D::setResolution(plint n, plint direction) {
    PLB_ASSERT( direction==0 || direction==1 || direction==2 );
    PLB_ASSERT( n>0 );
    double l = physDomain_.upperRightCorner[direction] -
               physDomain_.lowerLeftCorner[direction];
    set_dx(l/n);
}

void Units3D::set_periodic(bool periodicX, bool periodicY, bool periodicZ) {
    periodicX_ = periodicX;
    periodicY_ = periodicY;
    periodicZ_ = periodicZ;
}

void Units3D::extend(plint delta) {
    extend(delta, delta, delta, delta, delta, delta);
}

void Units3D::extend(plint dx0, plint dx1, plint dy0, plint dy1, plint dz0, plint dz1)
{
    // Makef sure the new domain is at least one node large.
    PLB_ASSERT( physDomain_.lowerLeftCorner[0]-dx0*dx_ + dx_ < physDomain_.upperRightCorner[0]+dx1*dx_ );
    PLB_ASSERT( physDomain_.lowerLeftCorner[1]-dy0*dx_ + dx_ < physDomain_.upperRightCorner[1]+dy1*dx_ );
    PLB_ASSERT( physDomain_.lowerLeftCorner[2]-dz0*dx_ + dx_ < physDomain_.upperRightCorner[2]+dz1*dx_ );

    physDomain_.lowerLeftCorner[0] -= dx0*dx_;
    physDomain_.lowerLeftCorner[1] -= dy0*dx_;
    physDomain_.lowerLeftCorner[2] -= dz0*dx_;

    physDomain_.upperRightCorner[0] += dx1*dx_;
    physDomain_.upperRightCorner[1] += dy1*dx_;
    physDomain_.upperRightCorner[2] += dz1*dx_;

    lbDomain_ = Box3D (
            0, lbDomain_.getNx()-1 + dx0 + dx1,
            0, lbDomain_.getNy()-1 + dy0 + dy1,
            0, lbDomain_.getNz()-1 + dz0 + dz1 );
}

void Units3D::computePhysDomain() {
    Array<double,3> llCorner(lbDomain_.x0*dx_, lbDomain_.y0*dx_, lbDomain_.z0*dx_);
    Array<double,3> urCorner(lbDomain_.x1*dx_, lbDomain_.y1*dx_, lbDomain_.z1*dx_);
    physDomain_ = Cuboid<double>(llCorner, urCorner);
}

}  // namespace plb
