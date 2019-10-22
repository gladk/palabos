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

/*
 * Palabos example program: Flow through a porous media, driven by a pressure-
 * gradient. This simple 100-liner illustrates how to formulate relatively
 * advanced programs.  It includes the creation of the geometry of a porous
 * media, definition of boundary conditions, and evaluation of the result. And,
 * as usual, all steps (including the creation of the porous media) are fully
 * parallelized.
 */

public class Porous {
   // This convenience function generates a list with "num_val" random numbers
   // comprized between 0 and max_val-1.
   public static java.util.Vector<Integer> randInt(java.util.Random handle, int numVal, int maxVal){
      java.util.Vector<Integer> randVect = new java.util.Vector<Integer>();
      for(int i=0; i<numVal; i++){
         randVect.add(handle.nextInt(maxVal));
      }
      return randVect;
   }

   public static void main(String[] argv){
      jlabos.JlabosBase base = jlabos.JlabosBase.getSingletonObject(); 

      // Initialize the random number generator. 
      java.util.Random random = new java.util.Random(10);

      double tau      = 0.6;      // Relaxation time.

      jlabos.core.plbInit();

      jlabos.Utils.pout("Creating the lattice");

      // This example uses a relatively small lattice, in order to get a result
      // quickly on a single processor. On a parallel machine, don't hesitate to
      // go for substantially larger systems.
      int lx, ly, lz;
      int nx, ny, nz;
      lx = 200;
      ly = 100;
      lz = 100;
      nx = lx+1;
      ny = ly+1;
      nz = lz+1;

      // The integer-valued block "media" shall hold a description of the porous media.
      jlabos.Box3D domain = new jlabos.Box3D(0, lx, 0, ly, 0, lz);

      jlabos.MultiNTensorField3DInt media = new jlabos.MultiNTensorField3DInt(domain, 1);

      jlabos.MultiNTensorField3DInt x = media.meshGrid(0);

      jlabos.MultiNTensorField3DInt y = media.meshGrid(1);

      jlabos.MultiNTensorField3DInt z = media.meshGrid(2);

      jlabos.Utils.pout("Creating the porous media");

      // The media consists of 50 spheres with radius 14 around a random position.
      int nobst = 50;
      int r     = 14;

      java.util.Vector<Integer> vectCx = randInt(random, nobst, nx);
      java.util.Vector<Integer> vectCy = randInt(random, nobst, ny);
      java.util.Vector<Integer> vectCz = randInt(random, nobst, nz);
      for(int i=0;i<nobst;i++){
         int cx = vectCx.get(i);
         int cy = vectCy.get(i);
         int cz = vectCz.get(i);
         jlabos.Utils.pout("  Creating a sphere centered at [" + cx + " " + cy + " " + cz +"]"); 

         jlabos.MultiNTensorField3DInt ysubcysquare = y.subtract(cy).toThePower(2);

         jlabos.MultiNTensorField3DInt zsubczsquare = z.subtract(cz).toThePower(2);

         jlabos.MultiNTensorField3DInt res          = x.subtract(cx).toThePower(2).add(ysubcysquare).add(zsubczsquare);

         jlabos.MultiNTensorField3DInt sphere       = res.lessThan((int)Math.pow(r,2));

         media = media.add(sphere);
      }

      // to see the media matrix under matlab:
      // jlabos.Utils.pout("Displaying the porous media");
      // jlabos.MultiNTensorField3DDouble mediaDouble = media.toDouble();
      // mediaDouble.writeToFile("media");
      // Matlab code
      // fp=load('media');
      // mat = reshape(fp, 101, 101, 201);
      // isosurface(mat);

      // A boolean operation is used to convert the integer values in "media" to
      // a zero-or-one representation.
      jlabos.MultiNTensorField3DInt isoSurface = media.greaterThan(0);

      // A relaxation parameter tau=0.6 is used here for illustration purposes. In
      // a real simulation, this should preferrably be computed from the fluid
      // viscosity, after proper unit version.
      jlabos.MultiBlockLattice3D lattice = new jlabos.MultiBlockLattice3D(domain, jlabos.Descriptor.D3Q19, new jlabos.BGK(1./tau));

      lattice.initializeAtEquilibrium(1., new double[] {0., 0., 0.}, domain);

      // Set a pressure boundary condition at inlet and outlet, with a larger
      // pressure for the inlet (remember that density is proportional to pressure).

      jlabos.OuterBoxBC.setPressureCondition(lattice, new jlabos.Box3D(0, 0, 1,ny-2, 1, nz-2));     

      jlabos.OuterBoxBC.setPressureCondition(lattice, new jlabos.Box3D(nx-1, nx-1, 1,ny-2, 1,nz-2));     

      lattice.setBoundaryDensity(1.02, new jlabos.Box3D(0,    0,    0, ny-1, 0, nz-1));

      lattice.setBoundaryDensity(1., new jlabos.Box3D(nx-1, nx-1, 0, ny-1, 0, nz-1));

      // Set bounce-back dynamics on all nodes on which "media" is non-zero.
      lattice.defineDynamics(media, new jlabos.BounceBack(1.));

      // Simulate 200 steps, and print statistics to the terminal every 20th step.
      int nsteps, nvisual;
      nsteps  = 200;
      nvisual = 20;

      jlabos.Utils.pout("Starting simulation (" +  nsteps +  " time steps)");
      long startTime = System.currentTimeMillis()/1000; 
      for(int i=0;i<nsteps; i++){
         if (i%nvisual==0 && i>0){
            long elapsed = System.currentTimeMillis()/1000-startTime;
            double Msups = i*nx*ny*nz / elapsed / 1e6;
            jlabos.Utils.pout(String.format("%d%% performed. Efficiency:",100*i/nsteps));
            jlabos.Utils.pout(String.format("%f Mega site updates per second.", Msups));
            jlabos.Utils.pout(String.format("The energy is %f", lattice.averageEnergy()));
         }
         lattice.collideAndStream();
      }


      // Compute the flow through the media on two selected y-z planes. The flow
      // must be computed on fluid-nodes only (excluding the spherical obstacles),
      // and it is based on the x-component of the velocity, normal to the plane.
      jlabos.MultiNTensorField3DInt flow1 = media.equals(0, new jlabos.Box3D(nx/3, nx/3, 0, ny, 0, nz));
      double flow1Average = lattice.velocityComponent(flow1, 0).average();

      jlabos.MultiNTensorField3DInt flow2 = media.equals(0, new jlabos.Box3D(2*nx/3, 2*nx/3, 0, ny, 0, nz));
      double flow2Average = lattice.velocityComponent(flow2, 0).average();

      // Note that the two values should be the same (up to the accuracy of the
      // numerical scheme), because the flow is incompressible. The difference
      // between the two values can indicate that the flow is still in an initial
      // transient state, or is under-resolved, or has a too large Mach number,
      // leading to compressibility effects.
      jlabos.Utils.pout("The average flow at 1/2 distance is " +  flow1Average +
            " and at 3/4 distance: " + flow2Average );

      // Display two snapshots of the (norm of the) strain-rate. On the first image,
      // the strain rate is evaluated with second-order accuracy from the velocity-
      // gradients, using a finite difference scheme. On the second image, it is
      // evaluated from the stress tensor, using the particle populations. The second
      // approach is first-order accurate, and the image obtained by this approach
      // has visibly lower quality.
      jlabos.MultiNTensorField3DDouble velocity = lattice.velocity().computeStrainRate().computeSymmetricTensorNorm(new jlabos.Box3D(0,nx, 0,ny, 20, 20));
      jlabos.Utils.saveAsJpg(velocity, "velocity");

      jlabos.MultiNTensorField3DDouble strainRateFromStress = lattice.strainRateFromStress().computeSymmetricTensorNorm(new jlabos.Box3D(0,nx, 0,ny, 20, 20));
      jlabos.Utils.saveAsJpg(strainRateFromStress, "strainRateFromStress");
   }
}
