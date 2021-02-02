/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

// VOTCA includes
#include "votca/xtp/erdiabatization.h"
#include "votca/xtp/ERIs.h"
#include <votca/tools/constants.h>

using std::flush;

namespace votca {
namespace xtp {

void ERDiabatization::setUpMatrices() {

  XTP_LOG(Log::debug, *_pLog) << "Setting up basis" << flush;
  this->_dftbasis = _orbitals.SetupDftBasis();
  this->_auxbasis = _orbitals.SetupAuxBasis();
  this->_bse_cmax = _orbitals.getBSEcmax();
  this->_bse_cmin = _orbitals.getBSEcmin();
  this->_bse_vmax = _orbitals.getBSEvmax();
  this->_bse_vmin = _orbitals.getBSEvmin();
  this->_bse_vtotal = _bse_vmax - _bse_vmin + 1;
  this->_bse_ctotal = _bse_cmax - _bse_cmin + 1;
  this->_basis = _orbitals.getBasisSetSize();
  this->_bse_size_ao = _basis * _basis;
  this->_occlevels = _orbitals.MOs().eigenvectors().block(
      0, _bse_vmin, _orbitals.MOs().eigenvectors().rows(), _bse_vtotal);
  this->_virtlevels = _orbitals.MOs().eigenvectors().block(
      0, _bse_cmin, _orbitals.MOs().eigenvectors().rows(), _bse_ctotal);

  _eris.Initialize_4c(_dftbasis);
}

void ERDiabatization::configure(const options_erdiabatization& opt) {
  _opt = opt;
}

double ERDiabatization::CalculateR(const Eigen::MatrixXd& D_LM,
                                   const Eigen::MatrixXd& D_JK) const {

  // Here I want to do \sum_{kl} (ij|kl) D^{LM}_{jk}. Is it right?
  XTP_LOG(Log::debug, *_pLog) << "Calculating 4c ERIs" << flush;
  Eigen::MatrixXd contracted = _eris.CalculateERIs_4c(D_LM, 1e-12);

  return D_JK.cwiseProduct(contracted).sum();
}

Eigen::MatrixXd ERDiabatization::CalculateU(const double phi) const {
  Eigen::MatrixXd U(2, 2);
  U(0, 0) = std::cos(phi);
  U(0, 1) = -1.0 * std::sin(phi);
  U(1, 0) = std::sin(phi);
  U(1, 1) = std::cos(phi);
  return U;
}

Eigen::MatrixXd ERDiabatization::Calculate_diabatic_H(
    const double E1, const double E2, const double angle) const {
  Eigen::VectorXd ad_energies(2);
  ad_energies << E1, E2;
  XTP_LOG(Log::error, *_pLog)
      << "Adiabatic energies in eV "
      << "E1: " << E1 * votca::tools::conv::hrt2ev
      << " E2: " << E2 * votca::tools::conv::hrt2ev << flush;
  XTP_LOG(Log::error, *_pLog)
      << "Rotation angle (degrees) " << angle * 57.2958 << flush;
  Eigen::MatrixXd U = CalculateU(angle);
  return U.transpose() * ad_energies.asDiagonal() * U;
}

Eigen::Tensor<double, 4> ERDiabatization::CalculateRtensor(
    const Orbitals& orb, QMStateType type) const {
  XTP_LOG(Log::debug, *_pLog) << "Computing R tensor" << flush;
  Eigen::Tensor<double, 4> r_tensor(2, 2, 2, 2);
  for (Index J = 0; J < 2; J++) {
    for (Index K = 0; K < 2; K++) {
      Eigen::MatrixXd D_JK = CalculateD(orb, type, J, K);
      for (Index L = 0; L < 2; L++) {
        for (Index M = 0; M < 2; M++) {
          Eigen::MatrixXd D_LM = CalculateD(orb, type, L, M);
          r_tensor(J, K, L, M) = CalculateR(D_LM, D_JK);
        }
      }
    }
  }
  return r_tensor;
}
Eigen::VectorXd ERDiabatization::CalculateER(const Orbitals& orb,
                                             QMStateType type) const {

  Eigen::Tensor<double, 4> R_JKLM = CalculateRtensor(orb, type);

  const double pi = votca::tools::conv::Pi;
  // Scanning through angles
  Eigen::VectorXd results = Eigen::VectorXd::Zero(360);
  // Initial mixing angle
  double phi_in = 0.;
  // Final mixing angle
  double phi_fin = 2. * pi;
  // We divide the interval into equal bits
  double step = (phi_fin - phi_in) / results.size();
  // Define angle we are iterating
  double phi;
  for (Index n = 0; n < results.size(); n++) {
    // Update angle
    phi = phi_in + n * step;
    Eigen::MatrixXd U = CalculateU(phi);
    // Complicated loop to handle. Can we make it a bit better?
    double result = 0.;
    for (Index I = 0; I < 2; I++) {
      for (Index J = 0; J < 2; J++) {
        for (Index K = 0; K < 2; K++) {
          for (Index L = 0; L < 2; L++) {
            for (Index M = 0; M < 2; M++) {
              result +=
                  U(J, I) * U(K, I) * U(L, I) * U(M, I) * R_JKLM(J, K, L, M);
            }
          }
        }
      }
    }
    results(n) = result;
  }
  return results;
}

Eigen::MatrixXd ERDiabatization::CalculateD(const Orbitals& orb,
                                            QMStateType type, Index stateindex1,
                                            Index stateindex2) const {

  Index index1;
  Index index2;

  if (stateindex1 == 0) {
    index1 = _opt.state_idx_1;
  }
  if (stateindex1 == 1) {
    index1 = _opt.state_idx_2;
  }
  if (stateindex2 == 0) {
    index2 = _opt.state_idx_1;
  }
  if (stateindex2 == 1) {
    index2 = _opt.state_idx_2;
  }

  Eigen::VectorXd exciton1;
  Eigen::VectorXd exciton2;
  if (type == QMStateType::Singlet) {
    exciton1 = orb.BSESinglets().eigenvectors().col(index1 - 1);
    exciton2 = orb.BSESinglets().eigenvectors().col(index2 - 1);
  } else {
    exciton1 = orb.BSETriplets().eigenvectors().col(index1 - 1);
    exciton2 = orb.BSETriplets().eigenvectors().col(index2 - 1);
  }
  Eigen::Map<const Eigen::MatrixXd> mat1(exciton1.data(), _bse_ctotal,
                                         _bse_vtotal);
  Eigen::Map<const Eigen::MatrixXd> mat2(exciton2.data(), _bse_ctotal,
                                         _bse_vtotal);

  // Here I ignored the diagonal term related to the stationary unexcited
  // electrons. It seems it doesn't play a huge role in the overall diabatization
  // process.
  Eigen::MatrixXd AuxMat_vv = mat1.transpose() * mat2;
  // This is the same as in the paper.
  Eigen::MatrixXd AuxMat_cc = mat1 * mat2.transpose();
  // This defines D = X + Y where X = occupied and Y = unoccupied contribution
  return _virtlevels * AuxMat_cc * _virtlevels.transpose() -
         _occlevels * AuxMat_vv * _occlevels.transpose();
}

}  // namespace xtp
}  // namespace votca
