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

public class Cylinder3d {
   public static void main(String[] argv){
      jlabos.JlabosBase base = jlabos.JlabosBase.getSingletonObject();
      int lx,ly,lz;              // Dimensions of the system in lattice units.
      lx = 400;
      ly = 50;
      lz = 50;
      double dx     = 1./lx;     // Discrete space step (ly is the reference length).
      int  nx,ny,nz;             // Lattice resolution.
      nx = lx+1;
      ny = ly+1;
      nz = lz+1;
      double rhoIn, rhoOut;
      rhoIn  = 1.02;
      rhoOut = 1.;
      double tau    = 0.55;      // Relaxation time.
      long maxIter  = 120;       // Maximum number of iterations.
      long statIter = 20;        // Intervals at which some statistics is printed.
      jlabos.core.plbInit();

      // Allocate memory for the block-lattice, and assign the dynamics.
      jlabos.Box3D domain = new jlabos.Box3D(0, lx, 0, ly, 0, lz);
      jlabos.MultiBlockLattice3D lattice = new jlabos.MultiBlockLattice3D(domain, jlabos.Descriptor.D3Q19, new jlabos.BGK(1./tau));

      jlabos.Utils.pout(jlabos.MultiBlockUtils.info(lattice));

      // Use a regularized boundary condition.
      jlabos.OuterBoxBC.setPressureCondition(lattice, domain);

      lattice.setBoundaryDensity(rhoIn, new jlabos.Box3D(0, 0, 1, ny-2, 1, nz-2));

      lattice.setBoundaryDensity(rhoOut, new jlabos.Box3D(nx-1, nx-1, 1, ny-2, 1, nz-2));

      // Here, x , y and z are (parallelized) matrices. It's similar to Matlab's meshgrid.
      jlabos.MultiNTensorField3DInt x = lattice.meshGrid(0);

      jlabos.MultiNTensorField3DInt y = lattice.meshGrid(1);

      jlabos.MultiNTensorField3DInt z = lattice.meshGrid(2);

      int cy = ny/2;
      int cz = nz/2;
      int rTube = 23;
   
      jlabos.MultiNTensorField3DInt ysubcy    = y.subtract(cy);

      jlabos.MultiNTensorField3DInt ysubcypow = ysubcy.toThePower(2);
      
      jlabos.MultiNTensorField3DInt zsubczpow = z.subtract(cz).toThePower(2);

      jlabos.MultiNTensorField3DInt tube      = ysubcypow.add(zsubczpow).lessThan((int)Math.pow(rTube,2));

      lattice.defineDynamics(tube.negate(), new jlabos.BounceBack(1.));
      int cx = nx/4;
      int rCyl = 4;

      jlabos.MultiNTensorField3DInt xsubcxpow     = x.subtract(cx).toThePower(2);

      jlabos.MultiNTensorField3DInt ysubcysub2pow = ysubcy.subtract(2).toThePower(2);

      jlabos.MultiNTensorField3DInt resCylinder   = xsubcxpow.add(ysubcysub2pow);

      jlabos.MultiNTensorField3DInt cylinder      = resCylinder.lessThan((int)Math.pow(rCyl,2));

      lattice.defineDynamics(cylinder, new jlabos.BounceBack(1.));


      lattice.initializeAtEquilibrium(1., new double[] {0., 0., 0.}, domain);

      jlabos.Utils.pout("Iterating from 0 to " + maxIter);

      java.util.Vector<jlabos.MultiNTensorField3DDouble> snapshots = new java.util.Vector<jlabos.MultiNTensorField3DDouble>();

      long startTime = System.currentTimeMillis()/1000;
      for (int i=0;i<maxIter+1;i++){
         if (i%statIter==0 && i>0){
            // Compute a velocity snapshot of the system, and keep it safe.
            snapshots.add(lattice.velocityNorm(new jlabos.Box3D(0, lx, 0, ly, cz, cz)));

            // Apply a reduction operation to compute the average kinetic energy.
            jlabos.Utils.pout("At iteration step " + i + " the energy is " + lattice.averageEnergy());
            long elapsed = System.currentTimeMillis()/1000-startTime+1; //we add 1 to avoid division by zero
            double Msups = i*nx*ny*nz / elapsed / 1e6;
            jlabos.Utils.pout("Efficiency: " + Msups + " Mega site updates per second.");
         }
         // This is the real computation: evolving through a collision-streaming cycle.
         lattice.collideAndStream();
      }

      jlabos.Utils.pout("Plotting results");
      int index = 0;
      for(jlabos.MultiNTensorField3DDouble item : snapshots){
         jlabos.Utils.saveAsJpg(item, "res" + index);
         item.writeToFile(new String("res"+ index++ + ".txt"));
      }
   }
}
