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

package jlabos;

public class MultiBlockUtils {
   public static String info(jlabos.MultiBlockLattice2D lattice){
      char nl = '\n';
      StringBuffer buf = new java.lang.StringBuffer();
      jlabos.double_D2Q9Descriptor_PlbMultiBlockLatticeInfo2D info = new jlabos.double_D2Q9Descriptor_PlbMultiBlockLatticeInfo2D(lattice.getObject());
      long nx = info.getNx();
      long ny = info.getNy();
      buf.append(nx).append(nl);
      buf.append(ny).append(nl);
      jlabos.Box2D smallest = info.getSmallestBlock();
      jlabos.Box2D largest = info.getLargestBlock();
      double millionAllocCells = info.getNumAllocatedCells()/1.e6;
      float cellPercentage    = (float)info.getNumAllocatedCells()/(float)(nx*ny)*100;
      buf.append("Internal structure of the " + nx + "-by-" + ny + " multi-block:" + nl);
      buf.append("Number of blocks in multi-block:" + info.getNumBlocks() + nl);
      buf.append("Smallest atomic-block: [" + smallest.getX0() + "," + smallest.getX1() + ", ");
      buf.append(smallest.getY0() + "," + smallest.getY1() + "]" + nl );
      buf.append("Largest atomic-block: [" + largest.getX0() + "," + largest.getX1() + ", ");
      buf.append(largest.getY0() + "," + largest.getY1() + "]").append(nl);
      buf.append("Number of allocated cells: " + millionAllocCells + "million").append(nl);
      buf.append("Percentage of allocated cells in multi-block: " + cellPercentage).append(nl);
      return buf.toString();
   }

   public static String info(jlabos.MultiBlockLattice3D lattice){
      char nl = '\n';
      StringBuffer buf = new java.lang.StringBuffer();
      jlabos.double_D3Q19Descriptor_PlbMultiBlockLatticeInfo3D info = new jlabos.double_D3Q19Descriptor_PlbMultiBlockLatticeInfo3D(lattice.getObject());
      long nx = info.getNx();
      long ny = info.getNy();
      long nz = info.getNz();
      buf.append(nx).append(nl);
      buf.append(ny).append(nl);
      buf.append(nz).append(nl);
      jlabos.Box3D smallest = info.getSmallestBlock();
      jlabos.Box3D largest = info.getLargestBlock();
      double millionAllocCells = info.getNumAllocatedCells()/1.e6;
      float cellPercentage    = (float)info.getNumAllocatedCells()/(float)(nx*ny*nz)*100;
      buf.append("Internal structure of the " + nx + "-by-" + ny + "-by-" + nz + " multi-block:").append(nl);
      buf.append("Number of blocks in multi-block:" + info.getNumBlocks()).append(nl);
      buf.append("Smallest atomic-block: [" + smallest.getX0() + "," + smallest.getX1() + ", ");
      buf.append(smallest.getY0() + "," + smallest.getY1() + ",");
      buf.append(smallest.getZ0() + "," + smallest.getZ1() + "]").append(nl);
      buf.append("Largest atomic-block: [" + largest.getX0() + "," + largest.getX1() + ", ");
      buf.append(largest.getY0() + "," + largest.getY1() + ",");
      buf.append(largest.getZ0() + "," + largest.getZ1() + "]").append(nl);
      buf.append("Number of allocated cells: " + millionAllocCells + "million").append(nl);
      buf.append("Percentage of allocated cells in multi-block: " + cellPercentage).append(nl);
      return buf.toString();
   }
}
