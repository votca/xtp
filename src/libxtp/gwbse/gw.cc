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

#include "votca/xtp/rpa.h"
#include "votca/xtp/sigma_ppm.h"
#include <votca/xtp/customtools.h>
#include <votca/xtp/gw.h>
#include <votca/xtp/sigma_spectral.h>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <eigen3/Eigen/src/Core/util/Constants.h>

namespace votca {
namespace xtp {

void GW::configure(const options& opt) {
  _opt = opt;
  _qptotal = _opt.qpmax - _opt.qpmin + 1;
  _rpa.configure(_opt.homo, _opt.rpamin, _opt.rpamax);
  if (_opt.sigma_integration == "ppm") {
    _sigma = std::unique_ptr<Sigma_base>(new Sigma_PPM(_Mmn, _rpa));
  } else if (_opt.sigma_integration == "exact") {
    _sigma = std::unique_ptr<Sigma_base>(new Sigma_Spectral(_Mmn, _rpa));
  }
  Sigma_base::options sigma_opt;
  sigma_opt.homo = _opt.homo;
  sigma_opt.qpmax = _opt.qpmax;
  sigma_opt.qpmin = _opt.qpmin;
  sigma_opt.rpamin = _opt.rpamin;
  sigma_opt.rpamax = _opt.rpamax;
  sigma_opt.eta = _opt.sigma_eta;
  _sigma->configure(sigma_opt);
  _Sigma_x = Eigen::MatrixXd::Zero(_qptotal, _qptotal);
  _Sigma_c = Eigen::MatrixXd::Zero(_qptotal, _qptotal);
  if (CustomOpts::GWSCExport()) {
    GWSelfConsistencyLogger::Initialize(_qptotal, _opt.g_sc_max_iterations);
  }
}

double GW::CalcHomoLumoShift() const {
  double DFTgap = _dft_energies(_opt.homo + 1) - _dft_energies(_opt.homo);
  double QPgap = _gwa_energies(_opt.homo + 1 - _opt.qpmin) -
                 _gwa_energies(_opt.homo - _opt.qpmin);
  return QPgap - DFTgap;
}

Eigen::VectorXd GW::CalcDiagonalEnergies() const {
  return _Sigma_x.diagonal() + _Sigma_c.diagonal() - _vxc.diagonal() +
         _dft_energies.segment(_opt.qpmin, _qptotal);
}

Eigen::MatrixXd GW::getHQP() const {
  return _Sigma_x + _Sigma_c - _vxc +
         Eigen::MatrixXd(
             _dft_energies.segment(_opt.qpmin, _qptotal).asDiagonal());
}

Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> GW::DiagonalizeQPHamiltonian()
    const {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> qpdiag(getHQP());
  PrintQP_Energies(qpdiag.eigenvalues());
  return qpdiag;
}

Eigen::MatrixXd GW::getGWAResults() const {
  Eigen::MatrixXd qp_energies_store = Eigen::MatrixXd::Zero(_qptotal, 5);
  qp_energies_store.col(0) = _dft_energies.segment(_opt.qpmin, _qptotal);
  qp_energies_store.col(1) = _Sigma_x.diagonal();
  qp_energies_store.col(2) = _Sigma_c.diagonal();
  qp_energies_store.col(3) = _vxc.diagonal();
  qp_energies_store.col(4) = _gwa_energies;
  return qp_energies_store;
}
void GW::PrintGWA_Energies() const {

  double shift = CalcHomoLumoShift();

  CTP_LOG(ctp::logINFO, _log)
      << (boost::format(
              "  ====== Perturbative quasiparticle energies (Hartree) ====== "))
             .str()
      << std::flush;
  CTP_LOG(ctp::logINFO, _log)
      << (boost::format("   DeltaHLGap = %1$+1.6f Hartree") % shift).str()
      << std::flush;

  for (int i = 0; i < _qptotal; i++) {
    std::string level = "  Level";
    if ((i + _opt.qpmin) == _opt.homo) {
      level = "  HOMO ";
    } else if ((i + _opt.qpmin) == _opt.homo + 1) {
      level = "  LUMO ";
    }

    CTP_LOG(ctp::logINFO, _log)
        << level
        << (boost::format(" = %1$4d DFT = %2$+1.4f VXC = %3$+1.4f S-X = "
                          "%4$+1.4f S-C = %5$+1.4f GWA = %6$+1.4f") %
            (i + _opt.qpmin) % _dft_energies(i + _opt.qpmin) % _vxc(i, i) %
            _Sigma_x(i, i) % _Sigma_c(i, i) % _gwa_energies(i))
               .str()
        << std::flush;
  }
  return;
}

void GW::PrintQP_Energies(const Eigen::VectorXd& qp_diag_energies) const {
  CTP_LOG(ctp::logDEBUG, _log)
      << ctp::TimeStamp() << " Full quasiparticle Hamiltonian  " << std::flush;
  CTP_LOG(ctp::logINFO, _log)
      << (boost::format(
              "  ====== Diagonalized quasiparticle energies (Hartree) "
              "====== "))
             .str()
      << std::flush;
  for (int i = 0; i < _qptotal; i++) {
    std::string level = "  Level";
    if ((i + _opt.qpmin) == _opt.homo) {
      level = "  HOMO ";
    } else if ((i + _opt.qpmin) == _opt.homo + 1) {
      level = "  LUMO ";
    }
    CTP_LOG(ctp::logINFO, _log)
        << level
        << (boost::format(" = %1$4d PQP = %2$+1.4f DQP = %3$+1.4f ") %
            (i + _opt.qpmin) % _gwa_energies(i) % qp_diag_energies(i))
               .str()
        << std::flush;
  }
  return;
}

Eigen::VectorXd GW::ScissorShift_DFTlevel(
    const Eigen::VectorXd& dft_energies) const {
  CTP_LOG(ctp::logDEBUG, _log)
      << ctp::TimeStamp() << " Scissor shifting DFT energies by: " << _opt.shift
      << " Hrt" << std::flush;
  Eigen::VectorXd shifted_energies = dft_energies;
  shifted_energies.segment(_opt.homo + 1, dft_energies.size() - _opt.homo - 1)
      .array() += _opt.shift;
  return shifted_energies;
}

bool GW::Converged(const Eigen::VectorXd& e1, const Eigen::VectorXd& e2,
                   double epsilon) const {
  int state = 0;
  bool energies_converged = true;
  double diff_max = (e1 - e2).cwiseAbs().maxCoeff(&state);
  if (diff_max > epsilon) {
    energies_converged = false;
  }
  if (tools::globals::verbose) {
    CTP_LOG(ctp::logDEBUG, _log)
        << ctp::TimeStamp() << " E_diff max=" << diff_max
        << " StateNo:" << state << std::flush;
  }
  return energies_converged;
}

Eigen::VectorXd GW::CalculateExcitationFreq(Eigen::VectorXd frequencies) {
  const double alpha = 0.0; // TODO: Mixing parameter
  // TODO: Make "Update" function that updates members variables: _Sigma_c,
  // _gwa_energies after each iteration.
  
  if (_opt.gw_sc_root_finder_method == 0) {
    // Fixed Point Method
    
    for (int i_freq = 0; i_freq < _opt.g_sc_max_iterations; ++i_freq) {
      _Sigma_c.diagonal() = _sigma->CalcCorrelationDiag(frequencies);
      _gwa_energies = CalcDiagonalEnergies();
      if (IterConverged(i_freq, frequencies)) {
        break;
      } else {
        frequencies = (1 - alpha) * _gwa_energies + alpha * frequencies;
      }
    }
    
  } else if (_opt.gw_sc_root_finder_method == 1 || _opt.gw_sc_root_finder_method == 2) {
    // Bisection Method, Regula Falsi Method

    // Define constants
    const bool regulaFalsi = _opt.gw_sc_root_finder_method == 2;
    const Eigen::VectorXd c = _Sigma_x.diagonal() - _vxc.diagonal() + _dft_energies.segment(_opt.qpmin, _qptotal);
    // First guess, two points required
    Eigen::VectorXd freq_1 = frequencies;
    Eigen::VectorXd freq_2 = _sigma->CalcCorrelationDiag(frequencies) + c; // Second point found by fixed point method
    // Compute left, right function values
    Eigen::VectorXd func_1 = _sigma->CalcCorrelationDiag(freq_1) + c - freq_1;
    Eigen::VectorXd func_2 = _sigma->CalcCorrelationDiag(freq_2) + c - freq_2;
    // Check whether a root is bounded
    Eigen::Array<bool, Eigen::Dynamic, 1> bounded = (func_1.cwiseProduct(func_2).array() <= 0);
    for (int i_freq = 0; i_freq < _opt.g_sc_max_iterations; ++i_freq) {
      // Compute next guess
      Eigen::VectorXd freq_3;
      if (regulaFalsi) {
        freq_3 = freq_1 - func_1.cwiseProduct(freq_2 - freq_1).cwiseQuotient(func_2 - func_1);
      } else {
        freq_3 = (freq_1 + freq_2) / 2;
      }
      // Without a bounded root, use fixed point method
      freq_3 = bounded.select(freq_3, _sigma->CalcCorrelationDiag(frequencies) + c);
      // Compute new function values
      Eigen::VectorXd sigc_3 = _sigma->CalcCorrelationDiag(freq_3);
      Eigen::VectorXd func_3 = sigc_3 + c - freq_3;
      // Update member variables
      _Sigma_c.diagonal() = sigc_3;
      _gwa_energies = freq_3;
      // Check converged
      if (IterConverged(i_freq, frequencies)) {
        break;
      } else {
        // Compute sign change
        Eigen::Array<bool, Eigen::Dynamic, 1> sign_1 = (func_1.cwiseProduct(func_3).array() <= 0);
        Eigen::Array<bool, Eigen::Dynamic, 1> sign_2 = (func_2.cwiseProduct(func_3).array() <= 0);
        // Update left, right frequencies
        freq_1 = sign_1.select(freq_1, freq_3); // TODO: Mixing
        freq_2 = sign_2.select(freq_2, freq_3);
        // Update left, right function values
        func_1 = sign_1.select(func_1, func_3);
        func_2 = sign_2.select(func_2, func_3);
        // Update current guess
        frequencies = freq_3; // TODO: Mixing
      }
    }
    
  } else if (_opt.gw_sc_root_finder_method == 3) {
    // Grid Method
    // TODO: Grid refinement

    // Options
    // TODO: Make grid root-finder base class? Or use function delegates?
    const int root_value_method = 1;
    const int root_score_method = 1;
    // Define constants
    const double rx0 = _opt.gw_sc_root_finder_range; // Range
    const int    nx  = _opt.gw_sc_root_finder_steps; // Steps
    const double inf = std::numeric_limits<double>::infinity();
    // Grid refinement
    const double f_refine     = _opt.gw_sc_root_finder_refine;
    const bool   b_refine     = f_refine > 0.0 && f_refine < 1.0;
    const double dx_min       = _opt.g_sc_limit * (nx - 1.0) / (2.0 * rx0);
    const int    n_refine_max = std::ceil(std::log(dx_min) / std::log(f_refine));

    int n_refine = 0; // Refinement number
    for (int i_freq = 0; i_freq < _opt.g_sc_max_iterations; ++i_freq) {
      // Grid refinement
      const double rx = rx0 * std::pow(f_refine, n_refine);//std::min(n_refine_max, n_refine));
      const double dx = (2.0 * rx) / (nx - 1.0);
      // Prepare vectors
      const Eigen::VectorXd xx_off = Eigen::VectorXd::LinSpaced(nx, -rx, +rx);
      const Eigen::VectorXd sx_vxc = _Sigma_x.diagonal() - _vxc.diagonal();
      // Evaluate sigma_c on all grid points
      Eigen::MatrixXd xx = Eigen::MatrixXd::Zero(nx, _qptotal); // TODO: Do not cache this
      Eigen::MatrixXd fx = Eigen::MatrixXd::Zero(nx, _qptotal);
      if (tools::globals::verbose) {
        CTP_LOG(ctp::logDEBUG, _log)
            << ctp::TimeStamp() << " Grid:"
            << " nx: " << nx << " rx: " << rx << " dx: " << dx << std::flush;
      }
      // TODO: Is it faster to first fill the columns and then transpose the entire matrix?
      for (int ix = 0; ix < nx; ix++) {
        Eigen::VectorXd xx_cur = frequencies.array() + xx_off[ix];
        xx.row(ix) = xx_cur;
        fx.row(ix) = _sigma->CalcCorrelationDiag(xx_cur);
      } // Grid point ix
      // For each state, find best root
      Eigen::VectorXd root_values = Eigen::VectorXd::Zero(_qptotal);
      Eigen::VectorXd root_scores = Eigen::VectorXd::Zero(_qptotal);
      // TODO: Multi-thread?
      for (int i_qp = 0; i_qp < _qptotal; i_qp++) {
        //    e_GW = sigma_c(e_GW) + sigma_x - v_xc + e_DFT
        // =>    0 = sigma_c(e_GW) + sigma_x - v_xc + e_DFT - e_GW = f(e_GW)
        const double c = sx_vxc[i_qp] + _dft_energies[_opt.qpmin + i_qp];
        Eigen::VectorXd xx_cur = xx.col(i_qp);             // lhs
        Eigen::VectorXd fx_cur = fx.col(i_qp).array() + c; // rhs
        Eigen::VectorXd gx_cur = fx_cur - xx_cur;          // target
        // Find best root
        double root_value_max =  0.0;
        double root_score_max = -1.0;
        int root_idx = 0;
        for (int ix = 0; ix < nx - 1; ix++) { // TODO: Loop only over sign-changes
          if (gx_cur[ix] * gx_cur[ix + 1] < 0.0) { // We have a sign change
            // Estimate the root
            double root_value_cur;
            if (root_value_method == 0 ) { // Average
              root_value_cur = (xx_cur[ix] + xx_cur[ix + 1]) / 2.0;
            } else if (root_value_method == 1) { // Fixed-point
              double dgdx = (gx_cur[ix + 1] - gx_cur[ix]) / (xx_cur[ix + 1] - xx_cur[ix]);
              root_value_cur = xx_cur[ix] - gx_cur[ix] / dgdx;
            } else {
              throw std::runtime_error("Grid root-finder: Invalid value method");
            }
            // Score the root
            double root_score_cur;
            if (root_score_method == 0 ) { // Distance
              root_score_cur = rx - std::abs(root_value_cur - frequencies[i_qp]);
            } else if (root_score_method == 1) { // Pole/spectral/QP weight
              double dfdx = (fx_cur[ix + 1] - fx_cur[ix]) / (xx_cur[ix + 1] - xx_cur[ix]);
              root_score_cur = 1.0 / (1.0 - dfdx); // Should be in (0, 1)
              if (root_score_cur < 1e-5) { continue; } // Invalid root
            } else {
              throw std::runtime_error("Grid root-finder: Invalid score method");
            }
            // Check if the root is better
            if (root_score_cur > root_score_max) { // We found a closer root
              root_value_max = root_value_cur;
              root_score_max = root_score_cur;
            }
            // Display the root
            if (tools::globals::verbose && i_qp <= _opt.homo) {
              CTP_LOG(ctp::logINFO, _log)
                  << boost::format(
                      "Level = %1$4d Index = %2$4d Value = %3$+1.6f Ha Score = %4$+1.6f") %
                      i_qp % root_idx % root_value_cur % root_score_cur
                  << std::flush;
            }
            root_idx++;
          }
        } // Grid point ix
        root_values[i_qp] = root_value_max;
        root_scores[i_qp] = root_score_max;
      } // State i_qp
      // Display all roots
      if (tools::globals::verbose) {
        CTP_LOG(ctp::logDEBUG, _log)
            << ctp::TimeStamp() << " QP roots " << std::flush;
        for (int i_qp = 0; i_qp < _qptotal; i_qp++) {
          CTP_LOG(ctp::logINFO, _log)
              << boost::format(
                  "Level = %1$4d E_0 = %2$+1.6f Ha E_GW = %3$+1.6f Ha Score = %4$+1.6f") %
                  i_qp % frequencies[i_qp] % root_values[i_qp] % root_scores[i_qp]
              << std::flush;
        }
      }
      // Update member variables
      _gwa_energies = (root_scores.array() >= 0).select(root_values, frequencies);
      _Sigma_c.diagonal() = _sigma->CalcCorrelationDiag(_gwa_energies);
      if (IterConverged(i_freq, frequencies)) {
        break;
      } else {
        frequencies = (1 - alpha) * _gwa_energies + alpha * frequencies;
        if (b_refine) {
          n_refine++;
        } else {
          break;
        }
      }
    }

  } else {
    throw std::runtime_error("Invalid GW SC root finder");
  }
  
  return frequencies;
}

bool GW::IterConverged(int i_freq, const Eigen::MatrixXd& frequencies) const {
  if (tools::globals::verbose) {
    CTP_LOG(ctp::logDEBUG, _log)
        << ctp::TimeStamp() << " G_Iteration:" << i_freq
        << " Shift[Hrt]:" << CalcHomoLumoShift() << std::flush;
  }
  if (CustomOpts::GWSCExport()) {
    GWSelfConsistencyLogger::LogFrequencies(frequencies);
  }
  if (Converged(_gwa_energies, frequencies, _opt.g_sc_limit)) {
    CTP_LOG(ctp::logDEBUG, _log)
        << ctp::TimeStamp() << " Converged after " << i_freq + 1
        << " G iterations." << std::flush;
    if (CustomOpts::GWSCExport()) {
      GWSelfConsistencyLogger::WriteGWIter(true);
    }
    return true;
  } else if (i_freq == _opt.g_sc_max_iterations - 1 &&
             _opt.g_sc_max_iterations > 1) {
    CTP_LOG(ctp::logDEBUG, _log)
        << ctp::TimeStamp()
        << " G-self-consistency cycle not converged after "
        << _opt.g_sc_max_iterations << " iterations." << std::flush;
    if (CustomOpts::GWSCExport()) {
      GWSelfConsistencyLogger::WriteGWIter(false);
    }
    return true;
  } else {
    return false;
  }
}

void GW::CalculateGWPerturbation() {

  _Sigma_x = (1 - _opt.ScaHFX) * _sigma->CalcExchange();
  CTP_LOG(ctp::logDEBUG, _log)
      << ctp::TimeStamp() << " Calculated Hartree exchange contribution  "
      << std::flush;
  // dft energies has size aobasissize
  // rpaenergies has siye rpatotal so has Mmn
  // gwaenergies/frequencies has qpmin,qpmax
  // homo index is relative to dftenergies
  Eigen::VectorXd dft_shifted_energies = ScissorShift_DFTlevel(_dft_energies);
  Eigen::VectorXd rpa_energies =
      dft_shifted_energies.segment(_opt.rpamin, _opt.rpamax - _opt.rpamin + 1);
  _rpa.setRPAInputEnergies(rpa_energies);
  Eigen::VectorXd frequencies =
      dft_shifted_energies.segment(_opt.qpmin, _qptotal);
  
  if (CustomOpts::GWEnergiesImport()) {
    CTP_LOG(ctp::logDEBUG, _log)
        << ctp::TimeStamp() << " Importing GW energies (/frequencies) " << std::flush;
    Eigen::VectorXd gw_energies_vec;
    Eigen::MatrixXd gw_energies_mat = CustomTools::ImportMatBinary("gw_energies.bin");
    gw_energies_vec.resize(gw_energies_mat.rows());
    gw_energies_vec << gw_energies_mat;
    frequencies = gw_energies_vec;
  }
  if (CustomOpts::SigmaExportRange() > 0 && !CustomOpts::SigmaExportConverged()) {
    ExportCorrelationDiags(frequencies);
  }
  
  for (int i_gw = 0; i_gw < _opt.gw_sc_max_iterations; ++i_gw) {

    if (i_gw % _opt.reset_3c == 0 && i_gw != 0) {
      _Mmn.Rebuild();
      CTP_LOG(ctp::logDEBUG, _log)
          << ctp::TimeStamp() << " Rebuilding 3c integrals" << std::flush;
    }
    _sigma->PrepareScreening();
    CTP_LOG(ctp::logDEBUG, _log)
        << ctp::TimeStamp() << " Calculated screening via RPA  " << std::flush;
    frequencies = CalculateExcitationFreq(frequencies);
    CTP_LOG(ctp::logDEBUG, _log)
        << ctp::TimeStamp() << " Calculated diagonal part of Sigma  "
        << std::flush;
    Eigen::VectorXd rpa_energies_old = _rpa.getRPAInputEnergies();
    _rpa.UpdateRPAInputEnergies(_dft_energies, frequencies, _opt.qpmin);

    CTP_LOG(ctp::logDEBUG, _log)
        << ctp::TimeStamp() << " GW_Iteration:" << i_gw
        << " Shift[Hrt]:" << CalcHomoLumoShift() << std::flush;
    if (Converged(_rpa.getRPAInputEnergies(), rpa_energies_old,
                  _opt.gw_sc_limit)) {
      CTP_LOG(ctp::logDEBUG, _log)
          << ctp::TimeStamp() << " Converged after " << i_gw + 1
          << " GW iterations." << std::flush;
      if (CustomOpts::GWSCExport()) {
        GWSelfConsistencyLogger::WriteCount(true);
      }
      break;
    } else if (i_gw == _opt.gw_sc_max_iterations - 1 &&
               _opt.gw_sc_max_iterations > 1) {
      CTP_LOG(ctp::logDEBUG, _log)
          << ctp::TimeStamp()
          << " WARNING! GW-self-consistency cycle not converged after "
          << _opt.gw_sc_max_iterations << " iterations." << std::flush;
      CTP_LOG(ctp::logDEBUG, _log)
          << ctp::TimeStamp()
          << "      Run continues. Inspect results carefully!" << std::flush;
      if (CustomOpts::GWSCExport()) {
        GWSelfConsistencyLogger::WriteCount(false);
      }
      break;
    } else {
      const double alpha = 0.0; // TODO: Mixing parameter
      rpa_energies = (1 - alpha) * rpa_energies + alpha * rpa_energies_old;
    }
  }

  PrintGWA_Energies();

  if (CustomOpts::SigmaExportRange() > 0 && CustomOpts::SigmaExportConverged()) {
    ExportCorrelationDiags(_gwa_energies);
  }
  if (CustomOpts::GWEnergiesExport() && !CustomOpts::GWEnergiesImport()) {
    CTP_LOG(ctp::logDEBUG, _log)
        << ctp::TimeStamp() << " Exporting gw energies (/frequencies) " << std::flush;
    Eigen::VectorXd gw_energies_vec = _gwa_energies;
    Eigen::MatrixXd gw_energies_mat;
    gw_energies_mat.resize(gw_energies_vec.size(), 1);
    gw_energies_mat << gw_energies_vec;
    CustomTools::ExportMatBinary("gw_energies.bin", gw_energies_mat);
  }
}

void GW::CalculateHQP() {
  _rpa.UpdateRPAInputEnergies(_dft_energies, _gwa_energies, _opt.qpmin);
  Eigen::VectorXd diag_backup = _Sigma_c.diagonal();
  _Sigma_c = _sigma->CalcCorrelationOffDiag(_gwa_energies);
  _Sigma_c.diagonal() = diag_backup;
  if (CustomOpts::SigmaMatrixExport()) {
    CTP_LOG(ctp::logDEBUG, _log)
        << ctp::TimeStamp() << " Exporting SigmaC matrix " << std::flush;
    if (CustomOpts::ExportBinary()) {
      CustomTools::ExportMatBinary("sigma_c_matrix.bin", _Sigma_c);
    } else {
      CustomTools::ExportMat("sigma_c_matrix.txt", _Sigma_c);
    }
  }
}

void GW::ExportCorrelationDiags(const Eigen::VectorXd& frequencies) const {
  const int range = CustomOpts::SigmaExportRange();
  const double delta = CustomOpts::SigmaExportDelta();
  const int size = 2 * range + 1;
  CTP_LOG(ctp::logDEBUG, _log)
      << ctp::TimeStamp() << " Exporting SigmaC diagonals " << std::flush;
  _sigma->PrepareScreening();
  Eigen::VectorXd offsets =
      Eigen::VectorXd::LinSpaced(size, -range * delta, range * delta);
  CTP_LOG(ctp::logDEBUG, _log)
      << ctp::TimeStamp() << " Calculating SigmaC diagonals " << std::flush;
  Eigen::MatrixXd results = Eigen::MatrixXd::Zero(_qptotal, size);
  for (int i = 0; i < size; i++) {
    results.col(i) = _sigma->CalcCorrelationDiag(
        frequencies + Eigen::VectorXd::Constant(_qptotal, offsets[i]));
    if ((i % ((size - 1) / 100)) == 0) {
      CTP_LOG(ctp::logDEBUG, _log)
          << ctp::TimeStamp() << " Progress: "
          << i << "/" << size << std::flush;
    }
  }
  CTP_LOG(ctp::logDEBUG, _log)
      << ctp::TimeStamp() << " Formatting SigmaC diagonals " << std::flush;
  Eigen::MatrixXd table = Eigen::MatrixXd::Zero(size + 1, _qptotal + 1);
  table.block(0, 1, 1, _qptotal) = frequencies.transpose();
  table.block(1, 0, size, 1) = offsets;
  table.block(1, 1, size, _qptotal) = results.transpose();
  CTP_LOG(ctp::logDEBUG, _log)
      << ctp::TimeStamp() << " Writing SigmaC diagonals " << std::flush;
  if (CustomOpts::ExportBinary()) {
    CustomTools::ExportMatBinary("sigma_c.bin", table);
  } else {
    CustomTools::ExportMat("sigma_c.txt", table);
  }
}

}  // namespace xtp
};  // namespace votca
