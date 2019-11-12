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

#pragma once
#ifndef _VOTCA_XTP_GW_H
#define _VOTCA_XTP_GW_H

#include "votca/xtp/logger.h"
#include <votca/xtp/orbitals.h>
#include <votca/xtp/rpa.h>
#include <votca/xtp/sigma_base.h>
#include <votca/xtp/threecenter.h>
namespace votca {
namespace xtp {

class GW {
 public:
  GW(Logger& log, TCMatrix_gwbse& Mmn, const Eigen::MatrixXd& vxc,
     const Eigen::VectorXd& dft_energies)
      : _log(log),
        _Mmn(Mmn),
        _vxc(vxc),
        _dft_energies(dft_energies),
        _rpa(log, Mmn){};

  struct options {
    Index homo;
    Index qpmin;
    Index qpmax;
    Index rpamin;
    Index rpamax;
    double g_sc_limit = 1e-5;  // default 1e-5
    Index g_sc_max_iterations = 50;
    double gw_sc_limit = 1e-5;
    Index gw_sc_max_iterations = 50;
    double shift = 0;
    double ScaHFX = 0.0;
    std::string sigma_integration = "ppm";
    double sigma_eta = 1e-3;
    Index reset_3c = 5;  // how often the 3c integrals in iterate should be
                         // rebuild
    // QP equation root finder
    int    gw_sc_root_finder_method = 0;
    double gw_sc_root_finder_range  = 1.0;
    Index  gw_sc_root_finder_steps  = 101;
    double gw_sc_root_finder_refine = 1.0; // Grid refinement factor for each gw iter
  };

  void configure(const options& opt);

  const Eigen::VectorXd& getGWAResults() const { return _gwa_energies; }
  // Calculates the diagonal elements up to self consistency
  void CalculateGWPerturbation();

  // Calculated offdiagonal elements as well
  void CalculateHQP();

  Eigen::MatrixXd getHQP() const;

  // Diagonalize QP particle Hamiltonian
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> DiagonalizeQPHamiltonian()
      const;

  void ExportSigmaC(Eigen::VectorXd frequencies) const;

 private:
  Index _qptotal;

  Eigen::VectorXd _gwa_energies;

  Eigen::MatrixXd _Sigma_x;
  Eigen::MatrixXd _Sigma_c;

  options _opt;

  std::unique_ptr<Sigma_base> _sigma = nullptr;
  Logger& _log;
  TCMatrix_gwbse& _Mmn;
  const Eigen::MatrixXd& _vxc;
  const Eigen::VectorXd& _dft_energies;

  RPA _rpa;

  Eigen::VectorXd CalculateExcitationFreq(Eigen::VectorXd frequencies);
  double CalcHomoLumoShift() const;
  Eigen::VectorXd ScissorShift_DFTlevel(
      const Eigen::VectorXd& dft_energies) const;
  void PrintQP_Energies(const Eigen::VectorXd& qp_diag_energies) const;
  void PrintGWA_Energies() const;
  Eigen::VectorXd CalcDiagonalEnergies() const;
  bool Converged(const Eigen::VectorXd& e1, const Eigen::VectorXd& e2,
                 double epsilon) const;
  bool IterConverged(int i_freq, const Eigen::MatrixXd& frequencies) const;
  void ExportCorrelationDiags(const Eigen::VectorXd& frequencies) const;
};
}  // namespace xtp
}  // namespace votca

#endif /* _VOTCA_XTP_BSE_H */
