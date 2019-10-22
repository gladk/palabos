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

#ifndef EMPIRICAL_DATA_H
#define EMPIRICAL_DATA_H

#include "core/globalDefs.h"

namespace plb {

/// Returns the drag coefficient.
/// CD = 2 Fd/(rho u^2 r^2 pi)
/// Input: Re = 2 r u / nu
double empirical_sphere_drag(double Re);

double computeTerminalVelocity(double densityRatio, double vSphere, double r, double kinematicViscosity, double gravity=9.8, bool doOutput=true);

double computeOptimalRadius(double targetRe, double densityRatio, double kinematicViscosity, double gravity=9.8, bool doOutput=true);

} // namespace plb

#endif  // EMPIRICAL_DATA_H
