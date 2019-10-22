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

/** \file
 * Helper functions for domain initialization -- header file.
 */
#ifndef NTENSOR_ANALYSIS_FUNCTIONAL_2D_H
#define NTENSOR_ANALYSIS_FUNCTIONAL_2D_H

#include "core/globalDefs.h"
#include "atomicBlock/dataProcessingFunctional2D.h"
#include "atomicBlock/reductiveDataProcessingFunctional2D.h"
#include "latticeBoltzmann/geometricOperationTemplates.h"

namespace plb {

template<typename T1, typename T2>
class CopyConvertNTensorFunctional2D : public BoxProcessingFunctional2D_NN<T1,T2>
{
public:
    virtual void process(Box2D domain, NTensorField2D<T1>& field1, NTensorField2D<T2>& field2);
    virtual CopyConvertNTensorFunctional2D<T1,T2>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    virtual BlockDomain::DomainT appliesTo() const;
};

}  // namespace plb

#endif  // NTENSOR_ANALYSIS_FUNCTIONAL_2D_H

