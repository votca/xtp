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
#include <votca/xtp/sigma_spectral.h>
#include <votca/xtp/threecenter.h>

namespace votca {
namespace xtp {

void Sigma_Spectral::PrepareScreening() {
  // Solve eigenvalue problem
  _EigenSol = _rpa.calculate_eigenvalues();
  // Cache residues
  _residues = CalcResidues();
  // Set Options
  _COHSEX = CustomOpts::COHSEX();
  return;
}

Eigen::VectorXd Sigma_Spectral::CalcCorrelationDiag(
    const Eigen::VectorXd& frequencies) const {
  const int rpasize = _EigenSol._Omega.size();
  Eigen::VectorXd result = Eigen::VectorXd::Zero(_qptotal);

  if (_COHSEX) {

#pragma omp parallel for
    for (int m = 0; m < _qptotal; m++) {
      double res = 0.0;
      const Eigen::MatrixXd& rm = _residues[m];
      for (int s = 0; s < rpasize; s++) {
        const Eigen::VectorXd rm_x_rm = rm.col(s).cwiseAbs2();
        double omega = _EigenSol._Omega(s);
        res += Equation48(rm_x_rm, omega);
      }  // Eigenvalues s
      result(m) = res;
    }  // State m

  } else {

#pragma omp parallel for
    for (int m = 0; m < _qptotal; m++) {
      double res = 0.0;
      const Eigen::MatrixXd& rm = _residues[m];
      for (int s = 0; s < rpasize; s++) {
        const Eigen::VectorXd rm_x_rm = rm.col(s).cwiseAbs2();
        double omega = _EigenSol._Omega(s);
        res += Equation47(rm_x_rm, omega, frequencies(m));
      }  // Eigenvalue s
      result(m) = res;
    }  // State m
  }

  return result;
}

Eigen::MatrixXd Sigma_Spectral::CalcCorrelationOffDiag(
    const Eigen::VectorXd& frequencies) const {
  const int rpasize = _EigenSol._Omega.size();
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(_qptotal, _qptotal);

  if (_COHSEX) {

#pragma omp parallel for
    for (int m = 0; m < _qptotal; m++) {
      const Eigen::MatrixXd& rm = _residues[m];
      for (int n = m + 1; n < _qptotal; n++) {
        double res = 0.0;
        const Eigen::MatrixXd& rn = _residues[n];
        for (int s = 0; s < rpasize; s++) {
          Eigen::VectorXd rm_x_rn =
              rm.col(s).cwiseProduct(rn.col(s));
          double omega = _EigenSol._Omega(s);
          res += Equation48(rm_x_rn, omega);
        }  // Eigenvalue s
        result(m, n) = res;
        result(n, m) = res;
      }  // State n
    }    // State m

  } else {

#pragma omp parallel for
    for (int m = 0; m < _qptotal; m++) {
      const Eigen::MatrixXd& rm = _residues[m];
      for (int n = m + 1; n < _qptotal; n++) {
        double res = 0.0;
        const Eigen::MatrixXd& rn = _residues[n];
        for (int s = 0; s < rpasize; s++) {
          Eigen::VectorXd rm_x_rn =
              rm.col(s).cwiseProduct(rn.col(s));
          double omega = _EigenSol._Omega(s);
          double res_m = Equation47(rm_x_rn, omega, frequencies(m));
          double res_n = Equation47(rm_x_rn, omega, frequencies(n));
          res += 0.5 * (res_m + res_n);
        }  // Eigenvalue s
        result(m, n) = res;
        result(n, m) = res;
      }  // State n
    }    // State m
  }

  return result;
}

std::vector<Eigen::MatrixXd> Sigma_Spectral::CalcResidues() const {
  const int lumo = _opt.homo + 1;
  const int n_occup = lumo - _opt.rpamin;
  const int n_unocc = _opt.rpamax - _opt.homo;
  const int rpasize = n_occup * n_unocc;
  const int qpoffset = _opt.qpmin - _opt.rpamin;
  const int auxsize = _Mmn.auxsize();
  vc2index vc = vc2index(0, 0, n_unocc);

  // Initialize residues object m*n*s
  std::vector<Eigen::MatrixXd> residues;
  residues.resize(_qptotal);
  for (int m = 0; m < _qptotal; m++ ) {
    residues[m] = Eigen::MatrixXd::Zero(_rpatotal, rpasize);
  }  // State m
  
  // To do the 4c integrals (mn|vc) efficiently, loop over m, v first
#pragma omp parallel for
  for (int m = 0; m < _qptotal; m++ ) {
    const Eigen::MatrixXd Mmn_mT =
        _Mmn[m + qpoffset].transpose();
    Eigen::MatrixXd res = Eigen::MatrixXd::Zero(_rpatotal, rpasize);
    for (int v = 0; v < n_occup; v++ ) { // Sum over v
      const Eigen::MatrixXd fcT = _Mmn[v].block(n_occup, 0, n_unocc, auxsize) * Mmn_mT; // Sum over chi
      res += fcT.transpose() * _EigenSol._XpY.block(vc.I(v, 0), 0, n_unocc, rpasize); // Sum over c
    }
    residues[m] += res;
  }
  
  return residues;
}

double Sigma_Spectral::Equation47(const Eigen::VectorXd& A12, double omega,
                                  double freq) const {
  const double eta = CustomOpts::SigmaSpectralEta();
  const int lumo = _opt.homo + 1;
  const int n_occup = lumo - _opt.rpamin;
  const int n_unocc = _opt.rpamax - _opt.homo;
  Eigen::ArrayXd B12 = _rpa.getRPAInputEnergies().array() + freq;
  B12.segment(0, n_occup) += omega;
  B12.segment(n_occup, n_unocc) -= omega;
  const Eigen::ArrayXd numer = A12.array() * B12;
  const Eigen::ArrayXd denom = B12.abs2() + eta;
  return (numer / denom).sum();
}

double Sigma_Spectral::Equation48(const Eigen::VectorXd& A12,
                                  double omega) const {
  const int lumo = _opt.homo + 1;
  const int n_occup = lumo - _opt.rpamin;
  double s1 = A12.head(n_occup).sum();
  double s2 = A12.sum();
  return 2 * (s1 - s2) / omega;
}

}  // namespace xtp
};  // namespace votca
