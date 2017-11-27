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

#ifndef MULTI_GRID_DATA_FIELD_2D_HH
#define MULTI_GRID_DATA_FIELD_2D_HH

#include "core/globalDefs.h"
#include "multiGrid/multiGridDataField2D.h"
#include "multiGrid/multiGridGenerator2D.h"


namespace plb {

template<typename T>
MultiGridScalarField2D<T>::MultiGridScalarField2D (
                        MultiGridManagement2D management_,
                        std::vector<BlockCommunicator2D* > communicators_,
                        std::vector<CombinedStatistics*> combinedStatistics_, 
                        plint behaviorLevel_ )
            : MultiGrid2D(management_, behaviorLevel_)
{
    allocateFields(communicators_,combinedStatistics_);
}

template<typename T>
MultiGridScalarField2D<T>::MultiGridScalarField2D (
                        MultiGridManagement2D management_,
                        plint behaviorLevel_ )
    : MultiGrid2D(management_, behaviorLevel_)
{
    allocateFields();
}

template<typename T>
MultiGridScalarField2D<T>::MultiGridScalarField2D(MultiGridScalarField2D<T> const& rhs)
    : MultiGrid2D(rhs)
{
    allocateFields();
}

template<typename T>
MultiGridScalarField2D<T>::MultiGridScalarField2D(MultiGrid2D const& rhs, Box2D subDomain, bool crop)
    : MultiGrid2D(extractManagement(rhs.getMultiGridManagement(),subDomain,crop),rhs.getBehaviorLevel())
{
    allocateFields();
}


template<typename T>
MultiGridScalarField2D<T>::~MultiGridScalarField2D(){
    for (plint iLevel=0; iLevel<(plint)fields.size(); ++iLevel){
        delete fields[iLevel];
    }
}

template<typename T>
void MultiGridScalarField2D<T>::reset(){
    for (plint iLevel=0; iLevel<(plint)fields.size(); ++iLevel){
        fields[iLevel]->reset();
    }
}

template<typename T>
T& MultiGridScalarField2D<T>::get(plint iX, plint iY){
    return fields[this->getBehaviorLevel()]->get(iX,iY);
}

template<typename T>
T const& MultiGridScalarField2D<T>::get(plint iX, plint iY) const {
    return fields[this->getBehaviorLevel()]->get(iX,iY);
}

template<typename T>
MultiScalarField2D<T>& MultiGridScalarField2D<T>::getComponent(plint level){
    PLB_PRECONDITION(level<(plint)fields.size());
    return *fields[level];
}

template<typename T>
MultiScalarField2D<T> const& MultiGridScalarField2D<T>::getComponent(plint level) const{
    PLB_PRECONDITION(level<(plint)fields.size());
    return *fields[level];
}

/// Using the management object, create the corresponding fields for each level
template<typename T>
void MultiGridScalarField2D<T>::allocateFields(){
    fields = generateScalarFields<T>( this->getMultiGridManagement(),
                                      defaultMultiGridPolicy2D().getBlockCommunicator<T>(this->getNumLevels()),
                                      defaultMultiGridPolicy2D().getCombinedStatistics(this->getNumLevels()) );
}

template<typename T>
void MultiGridScalarField2D<T>::allocateFields( std::vector<BlockCommunicator2D* > communicators,
                                                std::vector<CombinedStatistics*> combinedStatistics ){
    fields = generateScalarFields<T>( this->getMultiGridManagement(),
                                      communicators, combinedStatistics );
}



template<typename T>
int MultiGridScalarField2D<T>::getBlockId() const {
    return fields[0]->getStaticId();
}

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > MultiGridScalarField2D<T>::convertToCoarsest(plint dimDx, plint dimDt){
    
//TODO: modify avec new functions
//     plint levels = this->getNumLevels();
// 
//     MultiScalarField2D<T> *copy, *tmp;
//     copy = joinMultiScalarInCoarsest(
//         *fields[levels-2],*fields[levels-1], dimDx, dimDt );
// 
//     tmp=copy;
// 
//     for (plint iLevel=levels-2; iLevel>0; --iLevel){
//         copy = joinMultiScalarInCoarsest(
//                     *fields[iLevel-1],*tmp, dimDx, dimDt );
//         delete tmp; // erase the old value of copy
//         tmp = copy; // keep always a pointer over copy
//     }
//     
//     return std::auto_ptr<MultiScalarField2D<T> >(copy);
}

template<typename T>
std::auto_ptr<MultiScalarField2D<T> > MultiGridScalarField2D<T>::convertToFinest(plint dimDx, plint dimDt){
    
//TODO: modify avec new functions
//     MultiScalarField2D<T> *copy, *tmp;
//     copy = joinMultiScalarInFinest(
//         *fields[0],*fields[1], dimDx, dimDt );
//         
//     tmp=copy;
//                     
//     for (plint iLevel=2; iLevel<(plint)fields.size(); ++iLevel){
//         copy = joinMultiScalarInFinest(
//                     *tmp,*fields[iLevel], dimDx, dimDt );
//         delete tmp; // erase the old value of copy
//         tmp = copy; // keep always a pointer over copy
//     }
//         
//     return std::auto_ptr<MultiScalarField2D<T> >(copy);
}


/* ************** MultiGridTensorField2D ******************* */

template<typename T, int nDim>
MultiGridTensorField2D<T,nDim>::MultiGridTensorField2D ( MultiGridManagement2D management_,
                                                    std::vector<BlockCommunicator2D* > communicators_,
                                                    std::vector<CombinedStatistics*> combinedStatistics_, 
                                                    plint behaviorLevel_ )
{
    allocateFields(communicators_,combinedStatistics_);
}

template<typename T, int nDim>
MultiGridTensorField2D<T,nDim>::MultiGridTensorField2D (
                MultiGridManagement2D management_,
                plint behaviorLevel_ )
{
    allocateFields( defaultMultiGridPolicy2D().getBlockCommunicator<T>(this->getNumLevels()),
                    defaultMultiGridPolicy2D().getCombinedStatistics(this->getNumLevels()) );
}
                
template<typename T, int nDim>
MultiGridTensorField2D<T,nDim>::MultiGridTensorField2D(MultiGridTensorField2D<T,nDim> const& rhs)
{
    allocateFields();
}

template<typename T, int nDim>
MultiGridTensorField2D<T,nDim>::MultiGridTensorField2D(MultiGrid2D const& rhs)
        : MultiGrid2D(rhs.getMultiGridManagement(),rhs.getBehaviorLevel())
{
    allocateFields();
}

template<typename T, int nDim>
MultiGridTensorField2D<T,nDim>::MultiGridTensorField2D(
                        MultiGrid2D const& rhs, 
                        Box2D subDomain, bool crop)
        : MultiGrid2D(extractManagement(rhs.getMultiGridManagement(),subDomain,crop),rhs.getBehaviorLevel())
{
    allocateFields();
}


template<typename T, int nDim>
MultiGridTensorField2D<T,nDim>::~MultiGridTensorField2D(){
    for (pluint iLevel=0; iLevel<fields.size(); ++iLevel){
        delete fields[iLevel];
    }
}


template<typename T, int nDim>
void MultiGridTensorField2D<T,nDim>::reset(){
    for (pluint iLevel=0; iLevel<fields.size(); ++iLevel){
        fields[iLevel]->reset();
    }
}

template<typename T, int nDim>
Array<T,nDim>& MultiGridTensorField2D<T,nDim>::get(plint iX, plint iY){
    return fields[this->getBehaviorLevel()]->get(iX,iY);
}

template<typename T, int nDim>
Array<T,nDim> const& MultiGridTensorField2D<T,nDim>::get(plint iX, plint iY) const{
    return fields[this->getBehaviorLevel()]->get(iX,iY);
}

template<typename T, int nDim>
MultiTensorField2D<T,nDim>& MultiGridTensorField2D<T,nDim>::getComponent(plint level){
    PLB_PRECONDITION(level<(plint)fields.size());
    return *fields[level];
}

template<typename T, int nDim>
MultiTensorField2D<T,nDim> const& MultiGridTensorField2D<T,nDim>::getComponent(plint level) const{
    PLB_PRECONDITION(level<(plint)fields.size());
    return *fields[level];
}

template<typename T, int nDim>
int MultiGridTensorField2D<T,nDim>::getBlockId() const{
    return fields[this->getBehaviorLevel()]->getStaticId();
}

/// Using the management object, create the corresponding fields for each level
template<typename T, int nDim>
void MultiGridTensorField2D<T,nDim>::allocateFields(){
    fields = generateTensorFields<T,nDim>( this->getMultiGridManagement(),
                                      defaultMultiGridPolicy2D().getBlockCommunicator<T>(this->getNumLevels()),
                                      defaultMultiGridPolicy2D().getCombinedStatistics(this->getNumLevels()) );
}

template<typename T, int nDim>
void MultiGridTensorField2D<T,nDim>::allocateFields( std::vector<BlockCommunicator2D* > communicators,
                                                std::vector<CombinedStatistics*> combinedStatistics ){
    fields = generateTensorFields<T,nDim>(  this->getMultiGridManagement(),
                                            communicators, combinedStatistics );
}


template<typename T, int nDim>
std::auto_ptr<MultiTensorField2D<T,nDim> > MultiGridTensorField2D<T,nDim>::convertToCoarsest(plint dimDx, plint dimDt){
//TODO: modify avec new functions
//     plint levels = this->getNumLevels();
// 
//     MultiTensorField2D<T,nDim> *copy, *tmp;
//     copy = joinMultiTensorInCoarsest(
//         *fields[levels-2],*fields[levels-1], dimDx, dimDt );
// 
//     tmp=copy;
// 
//     for (plint iLevel=levels-2; iLevel>0; --iLevel){
//         copy = joinMultiTensorInCoarsest(
//                     *fields[iLevel-1],*tmp, dimDx, dimDt );
//         delete tmp; // erase the old value of copy
//         tmp = copy; // keep always a pointer over copy
//     }
//     
//     return std::auto_ptr<MultiTensorField2D<T,nDim> >(copy);
}

template<typename T, int nDim>
std::auto_ptr<MultiTensorField2D<T,nDim> > MultiGridTensorField2D<T,nDim>::convertToFinest(plint dimDx, plint dimDt){
//TODO: modify avec new functions
//     MultiTensorField2D<T,nDim> *copy, *tmp;
//     copy = joinMultiTensorInFinest(
//         *fields[0],*fields[1], dimDx, dimDt );
//         
//     tmp=copy;
//                     
//     for (plint iLevel=2; iLevel<(plint)fields.size(); ++iLevel){
//         copy = joinMultiTensorInFinest(
//                     *tmp,*fields[iLevel], dimDx, dimDt );
//         delete tmp; // erase the old value of copy
//         tmp = copy; // keep always a pointer over copy
//     }
//         
//     return std::auto_ptr<MultiTensorField2D<T,nDim> >(copy);
}




} // namespace plb

#endif  // MULTI_GRID_DATA_FIELD_2D_HH

