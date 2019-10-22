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
 * Helper functions for domain initialization -- precompilation.
 */

#include "plbWrapper/block/dataAnalysisFunctional2D.h"
#include "plbWrapper/block/dataAnalysisFunctional2D.hh"
#include "dataProcessors/ntensorAnalysisFunctional2D.h"
#include "dataProcessors/ntensorAnalysisFunctional2D.hh"

namespace plb {

template class BoxNTensorSumFunctional2D<PRECOMP_T>;
template class BoxNTensorMinFunctional2D<PRECOMP_T>;
template class BoxNTensorMaxFunctional2D<PRECOMP_T>;
template class BoundedBoxNTensorSumFunctional2D<PRECOMP_T>;
template class CopyConvertNTensorFunctional2D<PRECOMP_T,int>;
template class CopyConvertNTensorFunctional2D<PRECOMP_T,float>;
template class CopyConvertNTensorFunctional2D<PRECOMP_T,double>;
template class A_plus_alpha_NTensor2D<PRECOMP_T>;
template class A_minus_alpha_NTensor2D<PRECOMP_T>;
template class Alpha_minus_A_NTensor2D<PRECOMP_T>;
template class A_times_alpha_NTensor2D<PRECOMP_T>;
template class A_dividedBy_alpha_NTensor2D<PRECOMP_T>;
template class Alpha_dividedBy_A_NTensor2D<PRECOMP_T>;
template class A_toThePower_alpha_NTensor2D<PRECOMP_T>;
template class Alpha_toThePower_A_NTensor2D<PRECOMP_T>;
template class A_plus_alpha_inplace_NTensor2D<PRECOMP_T>;
template class A_minus_alpha_inplace_NTensor2D<PRECOMP_T>;
template class A_times_alpha_inplace_NTensor2D<PRECOMP_T>;
template class A_dividedBy_alpha_inplace_NTensor2D<PRECOMP_T>;
template class A_toThePower_alpha_inplace_NTensor2D<PRECOMP_T>;
template class A_plus_B_NTensor2D<PRECOMP_T>;
template class A_minus_B_NTensor2D<PRECOMP_T>;
template class A_times_B_NTensor2D<PRECOMP_T>;
template class A_dividedBy_B_NTensor2D<PRECOMP_T>;
template class A_toThePower_B_NTensor2D<PRECOMP_T>;
template class A_equals_B_NTensor2D<PRECOMP_T>;
template class A_equals_alpha_NTensor2D<PRECOMP_T>;
template class A_lessThan_B_NTensor2D<PRECOMP_T>;
template class A_lessThan_alpha_NTensor2D<PRECOMP_T>;
template class A_lessEqual_B_NTensor2D<PRECOMP_T>;
template class A_lessEqual_alpha_NTensor2D<PRECOMP_T>;
template class A_greaterThan_B_NTensor2D<PRECOMP_T>;
template class A_greaterThan_alpha_NTensor2D<PRECOMP_T>;
template class A_greaterEqual_B_NTensor2D<PRECOMP_T>;
template class A_greaterEqual_alpha_NTensor2D<PRECOMP_T>;
template class A_and_B_NTensor2D<PRECOMP_T>;
template class A_or_B_NTensor2D<PRECOMP_T>;
template class Not_A_NTensor2D<PRECOMP_T>;
template class A_plus_B_inplace_NTensor2D<PRECOMP_T>;
template class A_minus_B_inplace_NTensor2D<PRECOMP_T>;
template class A_times_B_inplace_NTensor2D<PRECOMP_T>;
template class A_dividedBy_B_inplace_NTensor2D<PRECOMP_T>;
template class A_toThePower_B_inplace_NTensor2D<PRECOMP_T>;
template class ExtractNTensorComponentFunctional2D<PRECOMP_T>;
template class ComputeNTensorNormFunctional2D<PRECOMP_T>;
template class ComputeNTensorNormSqrFunctional2D<PRECOMP_T>;
template class ComputeSymmetricNTensorNormFunctional2D<PRECOMP_T>;
template class ComputeSymmetricNTensorNormSqrFunctional2D<PRECOMP_T>;
template class ComputeSymmetricNTensorTraceFunctional2D<PRECOMP_T>;
template class BoxBulkNTensorVorticityFunctional2D<PRECOMP_T>;
template class BoxNTensorVorticityFunctional2D<PRECOMP_T>;
template class BoxBulkNTensorStrainRateFunctional2D<PRECOMP_T>;
template class BoxNTensorStrainRateFunctional2D<PRECOMP_T>;

template class MaskedBoxNTensorSumFunctional2D<PRECOMP_T>;
template class MaskedBoxNTensorMinFunctional2D<PRECOMP_T>;
template class MaskedBoxNTensorMaxFunctional2D<PRECOMP_T>;
template class BoundedMaskedBoxNTensorSumFunctional2D<PRECOMP_T>;
template class MaskedCopyConvertNTensorFunctional2D<PRECOMP_T,int>;
template class MaskedCopyConvertNTensorFunctional2D<PRECOMP_T,float>;
template class MaskedCopyConvertNTensorFunctional2D<PRECOMP_T,double>;
template class Masked_A_plus_alpha_NTensor2D<PRECOMP_T>;
template class Masked_A_minus_alpha_NTensor2D<PRECOMP_T>;
template class Masked_Alpha_minus_A_NTensor2D<PRECOMP_T>;
template class Masked_A_times_alpha_NTensor2D<PRECOMP_T>;
template class Masked_A_dividedBy_alpha_NTensor2D<PRECOMP_T>;
template class Masked_Alpha_dividedBy_A_NTensor2D<PRECOMP_T>;
template class Masked_A_toThePower_alpha_NTensor2D<PRECOMP_T>;
template class Masked_Alpha_toThePower_A_NTensor2D<PRECOMP_T>;
template class Masked_A_plus_alpha_inplace_NTensor2D<PRECOMP_T>;
template class Masked_A_minus_alpha_inplace_NTensor2D<PRECOMP_T>;
template class Masked_A_times_alpha_inplace_NTensor2D<PRECOMP_T>;
template class Masked_A_dividedBy_alpha_inplace_NTensor2D<PRECOMP_T>;
template class Masked_A_toThePower_alpha_inplace_NTensor2D<PRECOMP_T>;
template class Masked_A_plus_B_NTensor2D<PRECOMP_T>;
template class Masked_A_minus_B_NTensor2D<PRECOMP_T>;
template class Masked_A_times_B_NTensor2D<PRECOMP_T>;
template class Masked_A_dividedBy_B_NTensor2D<PRECOMP_T>;
template class Masked_A_toThePower_B_NTensor2D<PRECOMP_T>;
template class Masked_A_equals_B_NTensor2D<PRECOMP_T>;
template class Masked_A_equals_alpha_NTensor2D<PRECOMP_T>;
template class Masked_A_lessThan_B_NTensor2D<PRECOMP_T>;
template class Masked_A_lessThan_alpha_NTensor2D<PRECOMP_T>;
template class Masked_A_lessEqual_B_NTensor2D<PRECOMP_T>;
template class Masked_A_lessEqual_alpha_NTensor2D<PRECOMP_T>;
template class Masked_A_greaterThan_B_NTensor2D<PRECOMP_T>;
template class Masked_A_greaterThan_alpha_NTensor2D<PRECOMP_T>;
template class Masked_A_greaterEqual_B_NTensor2D<PRECOMP_T>;
template class Masked_A_greaterEqual_alpha_NTensor2D<PRECOMP_T>;
template class Masked_A_and_B_NTensor2D<PRECOMP_T>;
template class Masked_A_or_B_NTensor2D<PRECOMP_T>;
template class Masked_Not_A_NTensor2D<PRECOMP_T>;
template class Masked_A_plus_B_inplace_NTensor2D<PRECOMP_T>;
template class Masked_A_minus_B_inplace_NTensor2D<PRECOMP_T>;
template class Masked_A_times_B_inplace_NTensor2D<PRECOMP_T>;
template class Masked_A_dividedBy_B_inplace_NTensor2D<PRECOMP_T>;
template class Masked_A_toThePower_B_inplace_NTensor2D<PRECOMP_T>;
template class MaskedExtractNTensorComponentFunctional2D<PRECOMP_T>;
template class MaskedComputeNTensorNormFunctional2D<PRECOMP_T>;
template class MaskedComputeNTensorNormSqrFunctional2D<PRECOMP_T>;
template class MaskedComputeSymmetricNTensorNormFunctional2D<PRECOMP_T>;
template class MaskedComputeSymmetricNTensorNormSqrFunctional2D<PRECOMP_T>;
template class MaskedComputeSymmetricNTensorTraceFunctional2D<PRECOMP_T>;
template class MaskedBoxBulkNTensorVorticityFunctional2D<PRECOMP_T>;
template class MaskedBoxNTensorVorticityFunctional2D<PRECOMP_T>;
template class MaskedBoxBulkNTensorStrainRateFunctional2D<PRECOMP_T>;
template class MaskedBoxNTensorStrainRateFunctional2D<PRECOMP_T>;

}  // namespace plb
