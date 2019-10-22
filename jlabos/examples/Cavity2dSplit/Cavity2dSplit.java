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

public class Cavity2dSplit{
   public static void main(String[] argv){
      jlabos.JlabosBase base = jlabos.JlabosBase.getSingletonObject();

      int lx,ly;                       // Dimensions of the system in lattice units.
      lx = ly = 300;
      double dx = 1./lx;               // Discrete space step (ly is the reference length).
      int nx,ny;                       // Lattice resolution.
      double uMax     = 0.02;          // Maximum velocity in lattice units.
      double Re       = 100.;          // Reynolds number (wrt ly).
      double nu       = uMax*ly / Re;  // Viscosity in lattice units.
      double tau      = 3*nu+0.5;      // Relaxation time.
      long maxIter  = 4000;            // Maximum number of iterations.
      long statIter = 1000;            // Intervals at which some statistics is printed.

      nx = lx/2+1;
      ny = ly+1;

      jlabos.core.plbInit();

      jlabos.Utils.pout("dx : " + dx);
      jlabos.Utils.pout("nu : " + nu);
      jlabos.Utils.pout("tau : " + tau);

      // Allocate memory for the block-lattice, and assign the dynamics.

      jlabos.Box2D domain = new jlabos.Box2D(0, lx/2, 0, ly);

      jlabos.MultiBlockLattice2D latticeLeft = new jlabos.MultiBlockLattice2D(domain, jlabos.Descriptor.D2Q9, new jlabos.BGK(1./tau));
      jlabos.MultiBlockLattice2D latticeRight = new jlabos.MultiBlockLattice2D(domain, jlabos.Descriptor.D2Q9, new jlabos.BGK(1./tau));

      jlabos.Utils.pout(jlabos.MultiBlockUtils.info(latticeLeft));
      jlabos.Utils.pout(jlabos.MultiBlockUtils.info(latticeRight));

      
      jlabos.OuterBoxBC.setVelocityCondition(latticeLeft, domain);
      jlabos.OuterBoxBC.setVelocityCondition(latticeRight, domain);

      latticeLeft.setBoundaryVelocity(new double[] {0.,0.}, domain);
      latticeRight.setBoundaryVelocity(new double[] {0.,0.}, domain);

      // Top lid has constant in-plane velocity.
      latticeLeft.setBoundaryVelocity(new double[] {uMax,0.}, new jlabos.Box2D(1, nx-2, ny-1, ny-1));
      latticeRight.setBoundaryVelocity(new double[] {uMax,0.}, new jlabos.Box2D(1, nx-2, ny-1, ny-1));

      latticeLeft.initializeAtEquilibrium(1., new double[] {0.,0.}, domain);
      latticeRight.initializeAtEquilibrium(1., new double[] {0.,0.}, domain);

      java.util.Vector<jlabos.MultiNTensorField2DDouble> snapshotLeft = new java.util.Vector<jlabos.MultiNTensorField2DDouble>();
      java.util.Vector<jlabos.MultiNTensorField2DDouble> snapshotRight = new java.util.Vector<jlabos.MultiNTensorField2DDouble>();

      long startTime = System.currentTimeMillis()/1000;
      for(int i=0;i<maxIter+1; i++){
         if (i%statIter==0 && i>0){

            // Compute a velocity snapshot of the system, and keep it safe.
            snapshotLeft.add(latticeLeft.velocityNorm(domain));
            snapshotRight.add(latticeRight.velocityNorm(domain));

            // Apply a reduction operation to compute the average kinetic energy.
            jlabos.Utils.pout("Average energy : " + latticeLeft.averageEnergy());
            jlabos.Utils.pout("Average energy : " + latticeRight.averageEnergy());

            // Use Java's timer to compute the efficiency of the computations.
            long elapsed = System.currentTimeMillis()/1000-startTime;
            double Msups = i*nx*ny / elapsed / 1e6;
            jlabos.Utils.pout("Efficiency: " + Msups + " Mega site updates per second.");
         }

         // This is the real computation: evolving through a collision-streaming cycle.
         latticeLeft.collideAndStream();
         double[] rowRight = latticeRight.toArray(new jlabos.Box2D(2,2,1, ny-1));
         latticeLeft.fromArray(rowRight, new jlabos.Box2D(nx-1, nx-1, 1, ny-1));
         
         latticeRight.collideAndStream();
         double[] rowLeft = latticeLeft.toArray(new jlabos.Box2D(nx-2, nx-2, 1, ny-1));
         latticeRight.fromArray(rowLeft, new jlabos.Box2D(1,1,1, ny-1));

      }

      jlabos.Utils.pout("Plotting results");
      int index = 0;
      for(jlabos.MultiNTensorField2DDouble item : snapshotLeft){
         jlabos.Utils.saveAsJpg(item, "resLeft" + index, 800);
         item.writeToFile(new String("resLeft"+ index++ + ".txt"));
      }
      index = 0;
      for(jlabos.MultiNTensorField2DDouble item : snapshotRight){
         jlabos.Utils.saveAsJpg(item, "resRight" + index, 800);
         item.writeToFile(new String("resRight"+ index++ + ".txt"));
      }
   }
}
