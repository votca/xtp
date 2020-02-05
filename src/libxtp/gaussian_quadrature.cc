/*
 *            Copyright 2009-2019 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <votca/tools/constants.h>
#include <votca/xtp/gauss_hermite_quadrature_constants.h>
#include <votca/xtp/gauss_laguerre_quadrature_constants.h>
#include <votca/xtp/gaussian_quadrature.h>
#include <votca/xtp/threecenter.h>

namespace votca {
namespace xtp {

// Constructor
GaussianQuadrature::GaussianQuadrature(const Eigen::VectorXd& energies,
                                       const TCMatrix_gwbse& Mmn)
    : _energies(energies), _Mmn(Mmn) {}

void GaussianQuadrature::configure(options opt, const RPA& rpa) {
  _opt = opt;
  //Gauss_Laguerre_Quadrature_Constants glqc;
  Gauss_Hermite_Quadrature_Constants glqc;
  _quadpoints = glqc.getPoints(_opt.order);
  _quadadaptedweights = glqc.getAdaptedWeights(_opt.order);
  CalcDielInvVector(rpa);
}

// This function calculates and stores inverses of the microscopic dielectric
// matrix in a matrix vector
void GaussianQuadrature::CalcDielInvVector(const RPA& rpa) {

  _dielinv_matrices_r.resize(_opt.order);
#pragma openmp parallel schedule (guided)
  for (Index j = 0; j < _opt.order; j++) {
    Eigen::MatrixXcd eps_inv_j =
        rpa.calculate_epsilon_complex(0.0, _quadpoints(j)).inverse();
      eps_inv_j.diagonal().array() -= 1.0;
    _dielinv_matrices_r[j] = -eps_inv_j;  // (eps_inv_zero_alpha - eps_inv_j);
  }
}

double GaussianQuadrature::SigmaGQDiag(double frequency, Index gw_level) const {

  const Eigen::MatrixXd& Imx = _Mmn[gw_level];
  Eigen::ArrayXd DeltaE = frequency - _energies.array();
  std::complex<double> result = std::complex<double>(0.0, 0.0);
  for (Index j = 0; j < _opt.order; ++j) {
    Eigen::VectorXcd coeffs1 =
        (DeltaE) / (DeltaE.square() + std::pow(_quadpoints(j), 2));
    Eigen::MatrixXcd Amx = coeffs1.asDiagonal() * Imx;
    Eigen::MatrixXcd Cmx = Imx * (_dielinv_matrices_r[j]);  
    std::complex<double> value = (Cmx.cwiseProduct(Amx)).sum();
    result += _quadadaptedweights(j) * value;
  }

  return result.real() / (tools::conv::Pi);
}

double GaussianQuadrature::SigmaGQHDiag(double frequency, Index gw_level,
                                        double eta) const {
  Index homo = _opt.homo - _opt.rpamin;
  Index lumo = homo + 1;
  
  const Eigen::MatrixXd& Imx = _Mmn[gw_level];
  Eigen::ArrayXcd DeltaE = frequency - _energies.array();
  std::complex<double> result = std::complex<double>(0.0, 0.0);

  
  for ( Index k = 0; k < lumo; ++k){
    DeltaE(k) += std::complex<double>(0.0,1.0*eta); 
  }
  for ( Index k = lumo; k < _energies.size() ; ++k){
    DeltaE(k) += std::complex<double>(0.0,-1.0*eta); 
  }

  for (Index j = 0; j < _opt.order; ++j) {

    Eigen::VectorXcd denominator =
        (1.0) / (DeltaE + std::complex<double>(0.0, _quadpoints(j)));
    Eigen::MatrixXcd Amx = denominator.asDiagonal() * Imx;
    Eigen::MatrixXcd Cmx = Imx * (_dielinv_matrices_r[j]);  // C = I * B
    std::complex<double> value = (Cmx.cwiseProduct(Amx)).sum();
    result += _quadadaptedweights(j)  * value;
  }

  return result.real() / (2.0*tools::conv::Pi);
}

}  // namespace xtp

}  // namespace votca
