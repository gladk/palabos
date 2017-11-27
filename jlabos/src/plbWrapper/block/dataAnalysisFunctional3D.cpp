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

#include "plbWrapper/block/dataAnalysisFunctional3D.h"
#include "plbWrapper/block/dataAnalysisFunctional3D.hh"
#include "dataProcessors/ntensorAnalysisFunctional3D.h"
#include "dataProcessors/ntensorAnalysisFunctional3D.hh"

namespace plb {

template class BoxNTensorSumFunctional3D<PRECOMP_T>;
template class BoxNTensorMinFunctional3D<PRECOMP_T>;
template class BoxNTensorMaxFunctional3D<PRECOMP_T>;
template class BoundedBoxNTensorSumFunctional3D<PRECOMP_T>;
template class CopyConvertNTensorFunctional3D<PRECOMP_T,int>;
template class CopyConvertNTensorFunctional3D<PRECOMP_T,float>;
template class CopyConvertNTensorFunctional3D<PRECOMP_T,double>;
template class A_plus_alpha_NTensor3D<PRECOMP_T>;
template class A_minus_alpha_NTensor3D<PRECOMP_T>;
template class Alpha_minus_A_NTensor3D<PRECOMP_T>;
template class A_times_alpha_NTensor3D<PRECOMP_T>;
template class A_dividedBy_alpha_NTensor3D<PRECOMP_T>;
template class Alpha_dividedBy_A_NTensor3D<PRECOMP_T>;
template class A_toThePower_alpha_NTensor3D<PRECOMP_T>;
template class Alpha_toThePower_A_NTensor3D<PRECOMP_T>;
template class A_plus_alpha_inplace_NTensor3D<PRECOMP_T>;
template class A_minus_alpha_inplace_NTensor3D<PRECOMP_T>;
template class A_times_alpha_inplace_NTensor3D<PRECOMP_T>;
template class A_dividedBy_alpha_inplace_NTensor3D<PRECOMP_T>;
template class A_toThePower_alpha_inplace_NTensor3D<PRECOMP_T>;
template class A_plus_B_NTensor3D<PRECOMP_T>;
template class A_minus_B_NTensor3D<PRECOMP_T>;
template class A_times_B_NTensor3D<PRECOMP_T>;
template class A_dividedBy_B_NTensor3D<PRECOMP_T>;
template class A_toThePower_B_NTensor3D<PRECOMP_T>;
template class A_equals_B_NTensor3D<PRECOMP_T>;
template class A_equals_alpha_NTensor3D<PRECOMP_T>;
template class A_lessThan_B_NTensor3D<PRECOMP_T>;
template class A_lessThan_alpha_NTensor3D<PRECOMP_T>;
template class A_lessEqual_B_NTensor3D<PRECOMP_T>;
template class A_lessEqual_alpha_NTensor3D<PRECOMP_T>;
template class A_greaterThan_B_NTensor3D<PRECOMP_T>;
template class A_greaterThan_alpha_NTensor3D<PRECOMP_T>;
template class A_greaterEqual_B_NTensor3D<PRECOMP_T>;
template class A_greaterEqual_alpha_NTensor3D<PRECOMP_T>;
template class A_and_B_NTensor3D<PRECOMP_T>;
template class A_or_B_NTensor3D<PRECOMP_T>;
template class Not_A_NTensor3D<PRECOMP_T>;
template class A_plus_B_inplace_NTensor3D<PRECOMP_T>;
template class A_minus_B_inplace_NTensor3D<PRECOMP_T>;
template class A_times_B_inplace_NTensor3D<PRECOMP_T>;
template class A_dividedBy_B_inplace_NTensor3D<PRECOMP_T>;
template class A_toThePower_B_inplace_NTensor3D<PRECOMP_T>;
template class ExtractNTensorComponentFunctional3D<PRECOMP_T>;
template class ComputeNTensorNormFunctional3D<PRECOMP_T>;
template class ComputeNTensorNormSqrFunctional3D<PRECOMP_T>;
template class ComputeSymmetricNTensorNormFunctional3D<PRECOMP_T>;
template class ComputeSymmetricNTensorNormSqrFunctional3D<PRECOMP_T>;
template class ComputeSymmetricNTensorTraceFunctional3D<PRECOMP_T>;
template class BoxBulkNTensorVorticityFunctional3D<PRECOMP_T>;
template class BoxNTensorVorticityFunctional3D<PRECOMP_T>;
template class BoxBulkNTensorStrainRateFunctional3D<PRECOMP_T>;
template class BoxNTensorStrainRateFunctional3D<PRECOMP_T>;
template class UPO_ScalarProductFunctional3D<PRECOMP_T>;

template class MaskedBoxNTensorSumFunctional3D<PRECOMP_T>;
template class MaskedBoxNTensorMinFunctional3D<PRECOMP_T>;
template class MaskedBoxNTensorMaxFunctional3D<PRECOMP_T>;
template class BoundedMaskedBoxNTensorSumFunctional3D<PRECOMP_T>;
template class MaskedCopyConvertNTensorFunctional3D<PRECOMP_T,int>;
template class MaskedCopyConvertNTensorFunctional3D<PRECOMP_T,float>;
template class MaskedCopyConvertNTensorFunctional3D<PRECOMP_T,double>;
template class Masked_A_plus_alpha_NTensor3D<PRECOMP_T>;
template class Masked_A_minus_alpha_NTensor3D<PRECOMP_T>;
template class Masked_Alpha_minus_A_NTensor3D<PRECOMP_T>;
template class Masked_A_times_alpha_NTensor3D<PRECOMP_T>;
template class Masked_A_dividedBy_alpha_NTensor3D<PRECOMP_T>;
template class Masked_Alpha_dividedBy_A_NTensor3D<PRECOMP_T>;
template class Masked_A_toThePower_alpha_NTensor3D<PRECOMP_T>;
template class Masked_Alpha_toThePower_A_NTensor3D<PRECOMP_T>;
template class Masked_A_plus_alpha_inplace_NTensor3D<PRECOMP_T>;
template class Masked_A_minus_alpha_inplace_NTensor3D<PRECOMP_T>;
template class Masked_A_times_alpha_inplace_NTensor3D<PRECOMP_T>;
template class Masked_A_dividedBy_alpha_inplace_NTensor3D<PRECOMP_T>;
template class Masked_A_toThePower_alpha_inplace_NTensor3D<PRECOMP_T>;
template class Masked_A_plus_B_NTensor3D<PRECOMP_T>;
template class Masked_A_minus_B_NTensor3D<PRECOMP_T>;
template class Masked_A_times_B_NTensor3D<PRECOMP_T>;
template class Masked_A_dividedBy_B_NTensor3D<PRECOMP_T>;
template class Masked_A_toThePower_B_NTensor3D<PRECOMP_T>;
template class Masked_A_equals_B_NTensor3D<PRECOMP_T>;
template class Masked_A_equals_alpha_NTensor3D<PRECOMP_T>;
template class Masked_A_lessThan_B_NTensor3D<PRECOMP_T>;
template class Masked_A_lessThan_alpha_NTensor3D<PRECOMP_T>;
template class Masked_A_lessEqual_B_NTensor3D<PRECOMP_T>;
template class Masked_A_lessEqual_alpha_NTensor3D<PRECOMP_T>;
template class Masked_A_greaterThan_B_NTensor3D<PRECOMP_T>;
template class Masked_A_greaterThan_alpha_NTensor3D<PRECOMP_T>;
template class Masked_A_greaterEqual_B_NTensor3D<PRECOMP_T>;
template class Masked_A_greaterEqual_alpha_NTensor3D<PRECOMP_T>;
template class Masked_A_and_B_NTensor3D<PRECOMP_T>;
template class Masked_A_or_B_NTensor3D<PRECOMP_T>;
template class Masked_Not_A_NTensor3D<PRECOMP_T>;
template class Masked_A_plus_B_inplace_NTensor3D<PRECOMP_T>;
template class Masked_A_minus_B_inplace_NTensor3D<PRECOMP_T>;
template class Masked_A_times_B_inplace_NTensor3D<PRECOMP_T>;
template class Masked_A_dividedBy_B_inplace_NTensor3D<PRECOMP_T>;
template class Masked_A_toThePower_B_inplace_NTensor3D<PRECOMP_T>;
template class MaskedExtractNTensorComponentFunctional3D<PRECOMP_T>;
template class MaskedComputeNTensorNormFunctional3D<PRECOMP_T>;
template class MaskedComputeNTensorNormSqrFunctional3D<PRECOMP_T>;
template class MaskedComputeSymmetricNTensorNormFunctional3D<PRECOMP_T>;
template class MaskedComputeSymmetricNTensorNormSqrFunctional3D<PRECOMP_T>;
template class MaskedComputeSymmetricNTensorTraceFunctional3D<PRECOMP_T>;
template class MaskedBoxBulkNTensorVorticityFunctional3D<PRECOMP_T>;
template class MaskedBoxNTensorVorticityFunctional3D<PRECOMP_T>;
template class MaskedBoxBulkNTensorStrainRateFunctional3D<PRECOMP_T>;
template class MaskedBoxNTensorStrainRateFunctional3D<PRECOMP_T>;
template class Masked_UPO_ScalarProductFunctional3D<PRECOMP_T>;

}  // namespace plb
