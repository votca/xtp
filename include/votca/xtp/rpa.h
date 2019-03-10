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

#ifndef _VOTCA_XTP_RPA_H
#define _VOTCA_XTP_RPA_H
#include <vector>
#include <votca/xtp/eigen.h>
#include <votca/ctp/logger.h>

namespace votca {
namespace xtp {
class TCMatrix_gwbse;

struct rpa_eigensolution {
  Eigen::VectorXd _Omega; // Eigenvalues
  Eigen::MatrixXd _XpY; // Eigenvector components (X + Y)
};

class RPA {
 public:
  RPA(ctp::Logger &log, const TCMatrix_gwbse& Mmn) : _log(log), _Mmn(Mmn){};

  void configure(int homo, int rpamin, int rpamax) {
    _homo = homo;
    _rpamin = rpamin;
    _rpamax = rpamax;
  }
  
  rpa_eigensolution calculate_eigenvalues() const;

  Eigen::MatrixXd calculate_epsilon_i(double frequency) const {
    return calculate_epsilon<true>(frequency);
  }

  Eigen::MatrixXd calculate_epsilon_r(double frequency) const {
    return calculate_epsilon<false>(frequency);
  }

  const Eigen::VectorXd& getRPAInputEnergies() const { return _energies; }

  void setRPAInputEnergies(const Eigen::VectorXd& rpaenergies) {
    _energies = rpaenergies;
  }

  // calculates full RPA vector of energies from gwa and dftenergies and qpmin
  // RPA energies have three parts, lower than qpmin: dftenergies,between qpmin
  // and qpmax:gwa_energies,above:dftenergies+homo-lumo shift
  void UpdateRPAInputEnergies(const Eigen::VectorXd& dftenergies,
                              const Eigen::VectorXd& gwaenergies, int qpmin);

 private:
  int _homo;  // HOMO index with respect to dft energies
  int _rpamin;
  int _rpamax;
  const double _eta = 0.0001;

  Eigen::VectorXd _energies;

  const TCMatrix_gwbse& _Mmn;
  
  ctp::Logger &_log;
  
  // Bruneval, F. et al. molgw 1: Many-body perturbation theory software for
  // atoms, molecules, and clusters. Computer Physics Communications 208,
  // 149–161 (2016).
  // Eqs. 36-41
  Eigen::VectorXd calculate_spectral_AmB() const;
  Eigen::MatrixXd calculate_spectral_ApB() const;
  Eigen::MatrixXd calculate_spectral_C(const Eigen::VectorXd& AmB, const Eigen::MatrixXd& ApB) const;
  rpa_eigensolution diag_C(const Eigen::VectorXd& AmB, const Eigen::MatrixXd& C) const;

  template <bool imag>
  Eigen::MatrixXd calculate_epsilon(double frequency) const;
};
}  // namespace xtp
}  // namespace votca

#endif /* _VOTCA_RPA_RPA_H */
