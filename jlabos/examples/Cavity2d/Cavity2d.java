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
import java.io.*;
import java.util.*;
public class Cavity2d { 
   public static void main(String[] args) throws FileNotFoundException{

      // load all the needed libraries and initialize Palabos
      jlabos.JlabosBase base = jlabos.JlabosBase.getSingletonObject();

      // parse a parameter file in yaml format
      jlabos.Param param = new jlabos.Param("param.yml");

      // Dimensions of the system in lattice units.
      int lx = param.getIntValue(Arrays.asList("geometry", "lx"));
      int ly = param.getIntValue(Arrays.asList("geometry", "ly"));
      int nx = lx+1;
      int ny = ly+1;

      // Discrete space step (ly is the reference length).
      double dx = 1./lx; 
      // Maximum velocity in lattice units.
      double uMax = param.getDblValue(Arrays.asList("simulation", "uMax"));
      // Reynolds number (wrt ly)
      double Re   = param.getDblValue(Arrays.asList("simulation", "Re")); 
      double nu = uMax*ly / Re;  // Viscosity in lattice units.
      double tau = 3*nu+0.5;     // Relaxation time.
      // Maximum number of iterations.
      int maxIter = param.getIntValue(Arrays.asList("iteration", "maxIter")); 
      // Intervals at which some statistics is printed.
      int statIter = param.getIntValue(Arrays.asList("iteration", "statIter")); 

      // filename for writing output (without ext)
      String fName = param.getStrValue(Arrays.asList("general", "files", "out"));       


      jlabos.Utils.pout("dx : " + dx);
      jlabos.Utils.pout("nu : " + nu);
      jlabos.Utils.pout("tau : " + tau);

      // Allocate memory for the block-lattice, and assign the dynamics.

      jlabos.Box2D domain = new jlabos.Box2D(0, lx, 0, ly);

      jlabos.MultiBlockLattice2D lattice = new jlabos.MultiBlockLattice2D(domain, jlabos.Descriptor.D2Q9, new jlabos.BGK(1./tau));

      jlabos.Utils.pout(jlabos.MultiBlockUtils.info(lattice));

      jlabos.OuterBoxBC.setVelocityCondition(lattice, domain);

      lattice.setBoundaryVelocity(new double[] {0.,0.}, domain);

      // Top lid has constant in-plane velocity.
      lattice.setBoundaryVelocity(new double[] {uMax,0.}, new jlabos.Box2D(1, nx-2, ny-1, ny-1));

      lattice.initializeAtEquilibrium(1., new double[] {0.0,0.}, domain);

      java.util.Vector<jlabos.MultiNTensorField2DDouble> snapshot = new java.util.Vector<jlabos.MultiNTensorField2DDouble>();

      long startTime = System.currentTimeMillis()/1000;
      for(int i=0;i<maxIter+1; i++){
         if (i%statIter==0 && i>0){

            // Compute a velocity snapshot of the system, and keep it safe.
            snapshot.add(lattice.velocityNorm(domain));

            // Apply a reduction operation to compute the average kinetic energy.
            jlabos.Utils.pout("Average energy : " + lattice.averageEnergy());

            // Use Java's timer to compute the efficiency of the computations.
            long elapsed = System.currentTimeMillis()/1000-startTime;
            double Msups = i*nx*ny / elapsed / 1e6;
            jlabos.Utils.pout("Efficiency: " + Msups + " Mega site updates per second.");
         }
         // This is the real computation: evolving through a collision-streaming cycle.
         lattice.collideAndStream();
      }

      jlabos.Utils.pout("Plotting results");
      int index = 0;
      for(jlabos.MultiNTensorField2DDouble item : snapshot){
         jlabos.Utils.saveAsJpg(item, fName + index, 800);
         // write as svg if needed (large files!)
         //jlabos.Utils.saveAsSvg(item, "res" + index, 800);
         // write as text file if needed (to be used with octave for example)
         item.writeToFile(new String(fName + index++ + ".txt"));
      }

   }

}
