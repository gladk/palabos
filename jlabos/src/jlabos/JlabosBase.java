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

public class JlabosBase {
   private static JlabosBase singletonObject;
   private jlabos.Mpi mpi;
   private JlabosBase() {
      System.loadLibrary("_core");
      System.loadLibrary("_int_block");
      System.loadLibrary("_float_d3q19");
      System.loadLibrary("_float_d2q9");
      System.loadLibrary("_float_block");
      System.loadLibrary("_double_d3q19");
      System.loadLibrary("_double_d2q9");
      System.loadLibrary("_double_block");
      mpi = new jlabos.Mpi();
      mpi.loadMpi();
      jlabos.core.plbInit();
      Runtime.getRuntime().addShutdownHook(new Thread(new Runnable() {
               public void run() {
               finalizeMpi();
               }
               }));
   }

   public static synchronized JlabosBase getSingletonObject() {
      if (singletonObject == null) {
         singletonObject = new JlabosBase();
      }
      return singletonObject;
   }

   private int finalizeMpi(){
      return mpi.finalizeMpi();
   }

   public int getRankMpi(){
      return mpi.getRankMpi();
   }

   public Object clone() throws CloneNotSupportedException {
      throw new CloneNotSupportedException();
   }

   protected void finalize() throws Throwable {
      try {
         mpi.finalizeMpi();
      } finally {
         super.finalize();
      }
   }

}
