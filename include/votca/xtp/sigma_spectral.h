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

#ifndef _VOTCA_XTP_SIGMA_SPECTRAL_H
#define _VOTCA_XTP_SIGMA_SPECTRAL_H

#include <votca/xtp/rpa.h>
#include <votca/xtp/sigma_base.h>

namespace votca {
namespace xtp {

class TCMatrix_gwbse;
class RPA;

class Sigma_Spectral : public Sigma_base {

 public:
  Sigma_Spectral(TCMatrix_gwbse& Mmn, RPA& rpa) : Sigma_base(Mmn, rpa){};

  bool get_HedinApprox() { return _HedinApprox; }
  void set_HedinApprox(bool value) { _HedinApprox = value; }

  // Sets up the screening parametrisation
  void PrepareScreening();
  // Calculates Sigma_c diag elements
  Eigen::VectorXd CalcCorrelationDiag(const Eigen::VectorXd& frequencies) const;
  // Calculates Sigma_c offdiag elements
  Eigen::MatrixXd CalcCorrelationOffDiag(
      const Eigen::VectorXd& frequencies) const;

 private:
  bool _HedinApprox = false;  // Hedin's static approximation

  rpa_eigensolution _EigenSol;             // Eigenvalues, eigenvectors from RPA
  std::vector<Eigen::MatrixXd> _residues;  // Residues

  // Bruneval, F. et al. molgw 1: Many-body perturbation theory software for
  // atoms, molecules, and clusters. Computer Physics Communications 208,
  // 149–161 (2016).
  // Eq. 45, 47, 48
  std::vector<Eigen::MatrixXd> CalcResidues() const;
  double Equation47(const Eigen::VectorXd& A12, double omega, double freq) const;
  double Equation48(const Eigen::VectorXd& A12, double omega) const;
  // A12 = residues[m, :] .* residues[n, :]
};
}  // namespace xtp
}  // namespace votca

#endif /* _VOTCA_XTP_SIGMA_SPECTRAL_H */
