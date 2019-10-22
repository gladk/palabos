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

/* Main author: Orestis Malaspinas
 */

/** \file
 * Global Definitions -- header file.
 */
#ifndef CARREAU_GLOBAL_DEFS_H
#define CARREAU_GLOBAL_DEFS_H

#include "core/globalDefs.h"

namespace plb {

namespace global {

class CarreauParametersClass
{
public:
    void setNu0(double nu0_);
	void setNuInf(double nuInf_);
    void setLambda(double lambda_);
    void setExponent(double n_);
    
    double getNu0() const;
	double getNuInf() const;
    double getLambda() const;
    double getExponent() const;
private:
    CarreauParametersClass() {   };
private:
    double nu0, nuInf, lambda, n;
    friend CarreauParametersClass& CarreauParameters();
};
    
inline CarreauParametersClass& CarreauParameters() {
    static CarreauParametersClass singleton;
    return singleton;
}

}  // namespace global

}  // namespace plb

#endif  // CARREAU_GLOBAL_DEFS_H
