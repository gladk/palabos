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

/* Main author: Daniel Lagrava
 **/

#ifndef MULTI_GRID_DATA_FIELD_2D_H
#define MULTI_GRID_DATA_FIELD_2D_H

#include "core/globalDefs.h"

namespace plb {

template<typename T>
class MultiGridScalarField2D : public ScalarFieldBase2D<T>, public MultiGrid2D {
    public:
        MultiGridScalarField2D (
                        MultiGridManagement2D management_,
                        std::vector<BlockCommunicator2D* > communicators_,
                        std::vector<CombinedStatistics*> combinedStatistics_, 
                        plint behaviorLevel_=0 );

        MultiGridScalarField2D (
                        MultiGridManagement2D management_,
                        plint behaviorLevel_=0 );
                        
        
                        
        MultiGridScalarField2D(MultiGridScalarField2D<T> const& rhs);
        MultiGridScalarField2D(MultiGrid2D const& rhs);
        MultiGridScalarField2D(MultiGrid2D const& rhs, Box2D subDomain, bool crop=true);
        
        ~MultiGridScalarField2D();
        
    /* ************ ScalarFieldBase2D ************* */
    public:
        void reset();
        T& get(plint iX, plint iY);
        T const& get(plint iX, plint iY) const;
        
    /* ************ MultiGrid2D ************* */    
    public:
        MultiScalarField2D<T>& getComponent(plint level);
        MultiScalarField2D<T> const& getComponent(plint level) const;
        int getBlockId() const;
        
        std::auto_ptr<MultiScalarField2D<T> > convertToCoarsest(plint dimDx, plint dimDt);
        std::auto_ptr<MultiScalarField2D<T> > convertToFinest(plint dimDx, plint dimDt);
   
    private:
        /// Create the multiScalarFields for each level
        void allocateFields();
        void allocateFields( std::vector<BlockCommunicator2D* > communicators,
                             std::vector<CombinedStatistics*> combinedStatistics );
    private:
        std::vector<MultiScalarField2D<T> *> fields;
};


template<typename T, int nDim>
class MultiGridTensorField2D : public TensorFieldBase2D<T,nDim>, public MultiGrid2D {
    public:
        MultiGridTensorField2D (
                        MultiGridManagement2D management_,
                        std::vector<BlockCommunicator2D* > communicators_,
                        std::vector<CombinedStatistics*> combinedStatistics_, 
                        plint behaviorLevel_=0 );

        MultiGridTensorField2D (
                        MultiGridManagement2D management_,
                        plint behaviorLevel_=0 );
                        
        MultiGridTensorField2D(MultiGridTensorField2D<T,nDim> const& rhs);
        MultiGridTensorField2D(MultiGrid2D const& rhs);
        MultiGridTensorField2D(MultiGrid2D const& rhs, Box2D subDomain, bool crop=true);

        ~MultiGridTensorField2D();
        
    /* ************ TensorFieldBase2D ************* */
    public:
        void reset();
        virtual Array<T,nDim>& get(plint iX, plint iY);
        virtual Array<T,nDim> const& get(plint iX, plint iY) const;
        
    /* ************ MultiGrid2D ************* */    
    public:
        MultiTensorField2D<T,nDim>& getComponent(plint level);
        MultiTensorField2D<T,nDim> const& getComponent(plint level) const;
        int getBlockId() const;
        
        std::auto_ptr<MultiTensorField2D<T,nDim> > convertToCoarsest(plint dimDx, plint dimDt);
        std::auto_ptr<MultiTensorField2D<T,nDim> > convertToFinest(plint dimDx, plint dimDt);
        
    private:
        /// Create the multiTensorFields for each level
        void allocateFields();
        void allocateFields( std::vector<BlockCommunicator2D* > communicators,
                             std::vector<CombinedStatistics*> combinedStatistics );

    private:
        std::vector<MultiTensorField2D<T,nDim> *> fields;
};



} // namespace plb

#endif  // MULTI_GRID_DATA_FIELD_2D_H

