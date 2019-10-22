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

public class Cavity3d{
   public static void main(String[] argv){
      jlabos.JlabosBase base = jlabos.JlabosBase.getSingletonObject();

      int lx,ly,lz;                 // Dimensions of the system in lattice units.
      lx = ly = lz  = 100;
      double dx     = 1./lx;        // Discrete space step (ly is the reference length).
      int nx,ny,nz;                 // Lattice resolution.
      nx  = lx+1;
      ny  = ly+1;
      nz  = lz+1;
      double uMax   = 0.02;         // Maximum velocity in lattice units.
      double Re     = 100.;         // Reynolds number (wrt ly).
      double nu     = uMax*ly / Re; // Viscosity in lattice units.
      double tau    = 3*nu+0.5;     // Relaxation time.
      long maxIter  = 40;           // Maximum number of iterations.
      long statIter = 10;           // Intervals at which some statistics is printed.

      jlabos.core.plbInit(); 

      // Allocate memory for the block-lattice, and assign the dynamics.
      jlabos.Box3D domain = new jlabos.Box3D(0, lx, 0, ly, 0, lz);
      jlabos.MultiBlockLattice3D lattice = new jlabos.MultiBlockLattice3D(domain, jlabos.Descriptor.D3Q19, new jlabos.BGK(1./tau));

      jlabos.Utils.pout(jlabos.MultiBlockUtils.info(lattice));  

      // Use a regularized boundary condition.
      jlabos.OuterBoxBC.setVelocityCondition(lattice, domain);

      // Default all external walls to zero velocity.
      lattice.setBoundaryVelocity(new double[] {0.,0.,0.}, domain);

      // Top lid has constant in-plane velocity.
      lattice.setBoundaryVelocity(new double[] {uMax,0., 0.}, new jlabos.Box3D(1, nx-2, 1, ny-2, nz-1, nz-1));

      lattice.initializeAtEquilibrium(1., new double[] {0.,0.,0.}, domain);

      jlabos.MultiNTensorField3DDouble populations = lattice.getPopulations(new jlabos.Box3D(0, 10, 0, 10, 0, 10));

      populations.writeToFile("populations.txt");

      jlabos.Utils.pout("Iterating");

      java.util.Vector<jlabos.MultiNTensorField3DDouble> snapshots = new java.util.Vector<jlabos.MultiNTensorField3DDouble>();

      long startTime = System.currentTimeMillis()/1000;
      for(int i=0;i<maxIter+1; i++){
         if (i%statIter==0 && i>0){
            // Compute a velocity snapshot of the system, and keep it safe.
            snapshots.add(lattice.velocityNorm(new jlabos.Box3D(0, lx, ly/2, ly/2, 0, lz)));

            // Apply a reduction operation to compute the average kinetic energy.
            jlabos.Utils.pout("The energy is " + lattice.averageEnergy());

            long elapsed = System.currentTimeMillis()/1000-startTime;

            double Msups = i*nx*ny*nz / elapsed / 1e6;
            jlabos.Utils.pout("Efficiency: " + Msups + " Mega site updates per second.");
         }

         // This is the real computation: evolving through a collision-streaming cycle.
         lattice.collideAndStream();
      }

      long elapsed = System.currentTimeMillis()/1000-startTime;       
      double Msups = (maxIter+1)*nx*ny*nz / elapsed / 1e6;

      jlabos.Utils.pout("Efficiency: " + Msups + " Mega site updates per second.");

      jlabos.Utils.pout("Plotting results");

      int index = 0;
      for(jlabos.MultiNTensorField3DDouble item : snapshots){
         jlabos.Utils.saveAsJpg(item, "res" + index);
         item.writeToFile(new String("res"+ index++ + ".txt"));
      } 
   }
}

