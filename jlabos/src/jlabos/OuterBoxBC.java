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

public class OuterBoxBC {
   private void OuterBoxBC(){

   }

   public static void setVelocityCondition(MultiBlockLattice2D lattice){
   }

   public static void setVelocityCondition(MultiBlockLattice2D lattice, Box2D domain){
      lattice.getOuterBoxBC().setVelocityCondition(lattice.getObject(), domain);
   }

   public static void setVelocityCondition(MultiBlockLattice3D lattice, Box3D domain){
      lattice.getOuterBoxBC().setVelocityCondition(lattice.getObject(), domain);
   }

   public static void setPressureCondition(MultiBlockLattice2D lattice, Box2D domain){
      lattice.getOuterBoxBC().setPressureCondition(lattice.getObject(), domain);
   }

   public static void setPressureCondition(MultiBlockLattice3D lattice, Box3D domain){
      lattice.getOuterBoxBC().setPressureCondition(lattice.getObject(), domain);
   }

}
