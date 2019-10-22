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

/*
 * Main author: Yann Sagon <yann.sagon@unige.ch>
 */

public class Cylinder2d{
   public static void main(String[] argv){
      jlabos.JlabosBase base = jlabos.JlabosBase.getSingletonObject();

      int lx,ly;                     // Dimensions of the system in lattice units.
      lx = 800;
      ly = 120;       
      double dx     = 1./ly;         // Discrete space step (ly is the reference length).
      int  nx,ny;                    // Lattice resolution.
      nx = lx+1;
      ny = ly+1;    
      double uMax   = 0.02;          // Maximum velocity in lattice units.
      double Re     = 200.;          // Reynolds number (wrt ly).
      double nu     = uMax*ly / Re;  // Viscosity in lattice units.
      double tau    = 3*nu+0.5;      // Relaxation time.
      long maxIter  = 4000;          // Maximum number of iterations.
      long statIter = 1000;          // Intervals at which some statistics is printed.

      jlabos.core.plbInit();

      jlabos.Box2D myDomain = new jlabos.Box2D(0, lx, 0, ly);

      // Allocate memory for the block-lattice, and assign the dynamics.
      jlabos.Box2D domain = new jlabos.Box2D(0, lx, 0, ly);
      jlabos.MultiBlockLattice2D lattice = new jlabos.MultiBlockLattice2D(domain, jlabos.Descriptor.D2Q9, new jlabos.BGK(1./tau));

      jlabos.Utils.pout(jlabos.MultiBlockUtils.info(lattice));

      // Here, x and y are (parallelized) matrices. It's similar to Matlab's meshgrid.
      jlabos.MultiNTensorField2DInt x = lattice.meshGrid(0);
      jlabos.MultiNTensorField2DInt y = lattice.meshGrid(1);

      jlabos.MultiNTensorField2DDouble xdouble = x.toDouble();
      jlabos.MultiNTensorField2DDouble ydouble = y.toDouble();

      // Compute pressure and velocity according to the analytical Poiseuille solution.
      jlabos.MultiNTensorField2DDouble pois_pressure = xdouble.subtract(new double[] {lx/2.0}).multiply(new double[] { 8*nu*uMax/(ly*ly)}); 

      jlabos.MultiNTensorField2DDouble pois_velocity = new jlabos.MultiNTensorField2DDouble(domain, 2);

      pois_velocity.setToConstant(domain, new double[] {0.,0.});

      jlabos.MultiNTensorField2DDouble pois_res = (ydouble.subtract(ydouble.multiply(ydouble).multiply(new double[] {dx}))).multiply(new double[] {4.*uMax*dx});

      pois_velocity.setComponent(0, pois_res, domain);

      // Use a regularized boundary condition. Inlet and outlet are Poiseuille profile.
      jlabos.OuterBoxBC.setVelocityCondition(lattice, domain);

      lattice.setBoundaryVelocity(pois_velocity, domain);

      lattice.initializeAtEquilibrium(pois_pressure.add(new double[] {1.}), pois_velocity, domain);

      // Define the shape of a cylinder and associate bounce-back nodes to its domain.
      int cx,cy;
      cx = nx/4;
      cy = ny/2+2;
      int r  = 10;

      jlabos.MultiNTensorField2DInt cylinder = ((xdouble.subtract(new double[] {cx})).toThePower(new double[] {2.})).
         add((ydouble.subtract(new double[] {cy}).toThePower(new double[] {2.}))).
         lessThan(new double[] {Math.pow(r,2)});


      lattice.defineDynamics(cylinder, new jlabos.BounceBack(1.));

      jlabos.Utils.pout("Iterating");

      java.util.Vector<jlabos.MultiNTensorField2DDouble> snapshots = new java.util.Vector<jlabos.MultiNTensorField2DDouble>();

      long startTime = System.currentTimeMillis()/1000;
      for (int i=0;i<maxIter;i++){
         if (i%statIter==0 && i>0){
            // Compute a velocity snapshot of the system, and keep it safe.
            snapshots.add(lattice.velocityNorm(domain));

            // Apply a reduction operation to compute the average kinetic energy.
            jlabos.Utils.pout("The energy is " + lattice.averageEnergy(domain));

            long elapsed = System.currentTimeMillis()/1000-startTime+1; //we add 1 to avoid division by zero
            double Msups = i*nx*ny / elapsed / 1e6;
            jlabos.Utils.pout("Efficiency: " + Msups + " Mega site updates per second.");
         }
         // This is the real computation: evolving through a collision-streaming cycle.
         lattice.collideAndStream();
      }

      int index = 0;
      for(jlabos.MultiNTensorField2DDouble item : snapshots){
         jlabos.Utils.saveAsJpg(item, "res"+ index);
         item.writeToFile(new String("res"+ index++ + ".txt"));
      }
   }
}
