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

#ifndef CREATE_BUBBLES_3D_HH
#define CREATE_BUBBLES_3D_HH

#include "multiPhysics/createBubbles3D.h"
#include "offLattice/makeSparse3D.h"
#include "atomicBlock/dataProcessingFunctional3D.h"
#include "multiPhysics/freeSurfaceUtil3D.h"
#include <limits>


namespace plb {

template<typename T, template<typename U> class Descriptor>
void punchSphere(FreeSurfaceFields3D<T,Descriptor>& fields, Array<T,3> const& center, T radius,
                 T rhoEmpty, T rho0, Dynamics<T,Descriptor>& dynamics)
{
    applyProcessingFunctional (
            new PunchSphere3D<T,Descriptor>(center, radius, rho0),
            fields.lattice.getBoundingBox(), fields.freeSurfaceArgs );

    applyProcessingFunctional (
        new FreeSurfaceComputeInterfaceLists3D<T,Descriptor>(),
        fields.lattice.getBoundingBox(), fields.freeSurfaceArgs );

    applyProcessingFunctional (
        new FreeSurfaceIniInterfaceToAnyNodes3D<T,Descriptor>(rhoEmpty),
        fields.lattice.getBoundingBox(), fields.freeSurfaceArgs );

    applyProcessingFunctional (
        new FreeSurfaceIniEmptyToInterfaceNodes3D<T,Descriptor>(dynamics.clone(), Array<T,3>((T)0.,(T)0.,(T)0.)),
                                fields.lattice.getBoundingBox(),
                                fields.freeSurfaceArgs );

    applyProcessingFunctional (
        new FreeSurfaceRemoveFalseInterfaceCells3D<T,Descriptor>(rhoEmpty),
        fields.lattice.getBoundingBox(), fields.freeSurfaceArgs );

    applyProcessingFunctional (
        new FreeSurfaceEqualMassExcessReDistribution3D<T,Descriptor>(),
        fields.lattice.getBoundingBox(), fields.freeSurfaceArgs );

    applyProcessingFunctional (
        new FreeSurfaceComputeStatistics3D<T,Descriptor>,
        fields.lattice.getBoundingBox(), fields.freeSurfaceArgs );
}

template<typename T, template<typename U> class Descriptor>
void analyticalPunchSphere(FreeSurfaceFields3D<T,Descriptor>& fields, Array<T,3> const& center, T radius,
                           T rhoEmpty, T rho0, plint subDivision, Dynamics<T,Descriptor>& dynamics)
{
    applyProcessingFunctional (
            new AnalyticalPunchSphere3D<T,Descriptor>(center, radius, rho0, subDivision),
            fields.lattice.getBoundingBox(), fields.freeSurfaceArgs );

    applyProcessingFunctional (
        new FreeSurfaceComputeInterfaceLists3D<T,Descriptor>(),
        fields.lattice.getBoundingBox(), fields.freeSurfaceArgs );

    applyProcessingFunctional (
        new FreeSurfaceIniInterfaceToAnyNodes3D<T,Descriptor>(rhoEmpty),
        fields.lattice.getBoundingBox(), fields.freeSurfaceArgs );

    applyProcessingFunctional (
        new FreeSurfaceIniEmptyToInterfaceNodes3D<T,Descriptor>(dynamics.clone(), Array<T,3>((T)0.,(T)0.,(T)0.)),
                                fields.lattice.getBoundingBox(),
                                fields.freeSurfaceArgs );

    applyProcessingFunctional (
        new FreeSurfaceRemoveFalseInterfaceCells3D<T,Descriptor>(rhoEmpty),
        fields.lattice.getBoundingBox(), fields.freeSurfaceArgs );

    applyProcessingFunctional (
        new FreeSurfaceEqualMassExcessReDistribution3D<T,Descriptor>(),
        fields.lattice.getBoundingBox(), fields.freeSurfaceArgs );

    applyProcessingFunctional (
        new FreeSurfaceComputeStatistics3D<T,Descriptor>,
        fields.lattice.getBoundingBox(), fields.freeSurfaceArgs );
}

template<typename T, template<typename U> class Descriptor>
T computeAverageSphereDensity(FreeSurfaceFields3D<T,Descriptor>& fields, Array<T,3> const& center, T radius)
{
    CalculateAverageSphereDensity3D<T,Descriptor> functional(center, radius);
    applyProcessingFunctional (
            functional, fields.lattice.getBoundingBox(), fields.freeSurfaceArgs );
    return functional.getAverageDensity();
}

template<typename T, template<typename U> class Descriptor>
void punchSphere(FreeSurfaceSetup<T,Descriptor>& setup, Array<T,3> const& center, T radius,
                 T rhoEmpty, T rho0, Dynamics<T,Descriptor>& dynamics)
{
    applyProcessingFunctional (
            new PunchSphere3D<T,Descriptor>(center, radius, rho0),
            setup.getGroup().getBoundingBox(), setup.freeSurfaceArgs );

    applyProcessingFunctional (
        new FreeSurfaceComputeInterfaceLists3D<T,Descriptor>(),
        setup.getGroup().getBoundingBox(), setup.freeSurfaceArgs );

    applyProcessingFunctional (
        new FreeSurfaceIniInterfaceToAnyNodes3D<T,Descriptor>(rhoEmpty),
        setup.getGroup().getBoundingBox(), setup.freeSurfaceArgs );

    applyProcessingFunctional (
        new FreeSurfaceIniEmptyToInterfaceNodes3D<T,Descriptor>(dynamics.clone(), Array<T,3>((T)0.,(T)0.,(T)0.)),
                                setup.getGroup().getBoundingBox(),
                                setup.freeSurfaceArgs );

    applyProcessingFunctional (
        new FreeSurfaceRemoveFalseInterfaceCells3D<T,Descriptor>(rhoEmpty),
        setup.getGroup().getBoundingBox(), setup.freeSurfaceArgs );

    applyProcessingFunctional (
        new FreeSurfaceEqualMassExcessReDistribution3D<T,Descriptor>(),
        setup.getGroup().getBoundingBox(), setup.freeSurfaceArgs );

    applyProcessingFunctional (
        new FreeSurfaceComputeStatistics3D<T,Descriptor>,
        setup.getGroup().getBoundingBox(), setup.freeSurfaceArgs );
}

template<typename T, template<typename U> class Descriptor>
void analyticalPunchSphere(FreeSurfaceSetup<T,Descriptor>& setup, Array<T,3> const& center, T radius,
                           T rhoEmpty, T rho0, plint subDivision, Dynamics<T,Descriptor>& dynamics)
{
    applyProcessingFunctional (
            new AnalyticalPunchSphere3D<T,Descriptor>(center, radius, rho0, subDivision),
            setup.getGroup().getBoundingBox(), setup.freeSurfaceArgs );

    applyProcessingFunctional (
        new FreeSurfaceComputeInterfaceLists3D<T,Descriptor>(),
        setup.getGroup().getBoundingBox(), setup.freeSurfaceArgs );

    applyProcessingFunctional (
        new FreeSurfaceIniInterfaceToAnyNodes3D<T,Descriptor>(rhoEmpty),
        setup.getGroup().getBoundingBox(), setup.freeSurfaceArgs );

    applyProcessingFunctional (
        new FreeSurfaceIniEmptyToInterfaceNodes3D<T,Descriptor>(dynamics.clone(), Array<T,3>((T)0.,(T)0.,(T)0.)),
                                setup.getGroup().getBoundingBox(),
                                setup.freeSurfaceArgs );

    applyProcessingFunctional (
        new FreeSurfaceRemoveFalseInterfaceCells3D<T,Descriptor>(rhoEmpty),
        setup.getGroup().getBoundingBox(), setup.freeSurfaceArgs );

    applyProcessingFunctional (
        new FreeSurfaceEqualMassExcessReDistribution3D<T,Descriptor>(),
        setup.getGroup().getBoundingBox(), setup.freeSurfaceArgs );

    applyProcessingFunctional (
        new FreeSurfaceComputeStatistics3D<T,Descriptor>,
        setup.getGroup().getBoundingBox(), setup.freeSurfaceArgs );
}

template<typename T, template<typename U> class Descriptor>
T computeAverageSphereDensity(FreeSurfaceSetup<T,Descriptor>& setup, Array<T,3> const& center, T radius)
{
    CalculateAverageSphereDensity3D<T,Descriptor> functional(center, radius);
    applyProcessingFunctional (
            functional, setup.getGroup().getBoundingBox(), setup.freeSurfaceArgs );
    return functional.getAverageDensity();
}

}  // namespace plb

#endif  // CREATE_BUBBLES_3D_HH

