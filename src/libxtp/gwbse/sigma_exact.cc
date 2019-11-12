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

#include "votca/xtp/vc2index.h"
#include <votca/xtp/customtools.h>
#include <votca/xtp/rpa.h>
#include <votca/xtp/sigma_exact.h>
#include <votca/xtp/threecenter.h>

namespace votca {
namespace xtp {

void Sigma_Exact::PrepareScreening() {
  // Solve eigenvalue problem
  _EigenSol = _rpa.diagonalize_H2p();
  // Cache residues
  _residues = CalcResidues();
  return;
}

Eigen::VectorXd Sigma_Exact::CalcCorrelationDiag(
    const Eigen::VectorXd& frequencies) const {
  const int number_eigenvectors = _EigenSol._omega.size();
  Eigen::VectorXd result = Eigen::VectorXd::Zero(_qptotal);

#pragma omp parallel for
  for (int m = 0; m < _qptotal; m++) {
    double res = 0.0;
    const Eigen::MatrixXd& rm = _residues[m];
    for (int s = 0; s < number_eigenvectors; s++) {
      const Eigen::VectorXd rm_x_rm = rm.col(s).cwiseAbs2();
      double eigenvalue = _EigenSol._omega(s);
      res += Equation47(rm_x_rm, eigenvalue, frequencies(m));
    }
    // Multiply with factor 2.0 to sum over both (identical) spin states
    result(m) = 2.0 * res;
  }  // State m

  return result;
}

Eigen::MatrixXd Sigma_Exact::CalcCorrelationOffDiag(
    const Eigen::VectorXd& frequencies) const {
  const int rpasize = _EigenSol._omega.size();
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(_qptotal, _qptotal);

  if (CustomOpts::SigmaCNoOffdiags()) {
    return result;
  }

#pragma omp parallel for
  for (int m = 0; m < _qptotal; m++) {
    const Eigen::MatrixXd& rm = _residues[m];
    for (int n = m + 1; n < _qptotal; n++) {
      double res = 0.0;
      const Eigen::MatrixXd& rn = _residues[n];
      for (int s = 0; s < rpasize; s++) {
        Eigen::VectorXd rm_x_rn = rm.col(s).cwiseProduct(rn.col(s));
        double eigenvalue = _EigenSol._omega(s);
        double res_m = Equation47(rm_x_rn, eigenvalue, frequencies(m));
        double res_n = Equation47(rm_x_rn, eigenvalue, frequencies(n));
        res += res_m + res_n;
      }  // Eigenvalue s
      // Multiply with factor 2.0 to sum over both (identical) spin states
      result(m, n) = 0.5 * 2.0 * res;
      result(n, m) = 0.5 * 2.0 * res;
    }  // State n
  }    // State m

  return result;
}

std::vector<Eigen::MatrixXd> Sigma_Exact::CalcResidues() const {
  const int lumo = _opt.homo + 1;
  const int n_occup = lumo - _opt.rpamin;
  const int n_unocc = _opt.rpamax - _opt.homo;
  const int rpasize = n_occup * n_unocc;
  const int qpoffset = _opt.qpmin - _opt.rpamin;
  const int auxsize = _Mmn.auxsize();
  vc2index vc = vc2index(0, 0, n_unocc);
  std::vector<Eigen::MatrixXd> residues(_qptotal);

  // To do the 4c integrals (mn|vc) efficiently, loop over m, v first
#pragma omp parallel for
  for (int m = 0; m < _qptotal; m++) {
    const Eigen::MatrixXd Mmn_mT = _Mmn[m + qpoffset].transpose();
    Eigen::MatrixXd res = Eigen::MatrixXd::Zero(_rpatotal, rpasize);
    for (int v = 0; v < n_occup; v++) {  // Sum over v
      const Eigen::MatrixXd fc =
          _Mmn[v].block(n_occup, 0, n_unocc, auxsize) * Mmn_mT;  // Sum over chi
      res += fc.transpose() * _EigenSol._XpY.block(vc.I(v, 0), 0, n_unocc,
                                                   rpasize);  // Sum over c
    }
    residues[m] = res;
  }

  return residues;
}

double Sigma_Exact::Equation47(const Eigen::VectorXd& A12, double eigenvalue,
                               double freq) const {
  const double eta = _opt.eta;
  const int lumo = _opt.homo + 1;
  const int n_occup = lumo - _opt.rpamin;
  const int n_unocc = _opt.rpamax - _opt.homo;
  Eigen::ArrayXd B12 = -_rpa.getRPAInputEnergies().array() + freq;
  B12.segment(0, n_occup) += eigenvalue;
  B12.segment(n_occup, n_unocc) -= eigenvalue;
  const Eigen::ArrayXd numer = A12.array() * B12;
  const Eigen::ArrayXd denom = B12.abs2() + eta * eta;
  return (numer / denom).sum();
}

double Sigma_Exact::Equation48(const Eigen::VectorXd& A12,
                               double eigenvalue) const {
  const double eta = _opt.eta;
  const int lumo = _opt.homo + 1;
  const int n_occup = lumo - _opt.rpamin;
  const int n_unocc = _opt.rpamax - _opt.homo;
  Eigen::ArrayXd B12 =
      Eigen::VectorXd::Zero(_rpatotal);  // eigenvalue >> |freq - energy|
  B12.segment(0, n_occup) += eigenvalue;
  B12.segment(n_occup, n_unocc) -= eigenvalue;
  const Eigen::ArrayXd numer = A12.array() * B12;
  const Eigen::ArrayXd denom = B12.abs2() + eta * eta;
  return (numer / denom).sum();
}

}  // namespace xtp
};  // namespace votca
