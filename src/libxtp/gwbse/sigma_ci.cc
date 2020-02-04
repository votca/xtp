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
#include <votca/xtp/gw.h>
#include <votca/xtp/sigma_ci.h>

namespace votca {
namespace xtp {

void Sigma_CI::PrepareScreening() {
  GaussianQuadrature::options opt;
  opt.homo = _opt.homo;
  opt.order = _opt.order;
  opt.qptotal = _qptotal;
  opt.qpmin = _opt.qpmin;
  opt.rpamin = _opt.rpamin;
  opt.alpha = _opt.alpha;
  _gq.configure(opt, _rpa);
}

// This function is used in the calculation of the residues
// This calculates eps^-1 (inverse of the dielectric function) for complex
// frequencies of the kind omega = delta + i*eta
double Sigma_CI::CalcDiagContribution(Eigen::RowVectorXd Imx_row, double delta,
                                      double eta) const {

  // Eigen::MatrixXd DielMxInv = _rpa.calculate_real_epsilon_inverse(delta,
  // eta);
  Eigen::MatrixXcd DielMxInv =
      _rpa.calculate_epsilon_complex(delta, eta).inverse();
  DielMxInv.diagonal().array() -= 1.0;
  std::complex<double> lambda =
      ((Imx_row * DielMxInv).cwiseProduct(Imx_row)).sum();
  return lambda.real();
}

// This function is used in the calculation of the residues
// This calculates eps^-1 (inverse of the dielectric function) for complex
// frequency omega = 0 + i* 0 (origin)
double Sigma_CI::CalcDiagContributionValue_alpha(Eigen::RowVectorXd Imx_row,
                                                 double delta,
                                                 double alpha) const {
  Eigen::MatrixXcd R = _rpa.calculate_epsilon_complex(0.0, 0.0).inverse();
  // Eigen::MatrixXd R = _rpa.calculate_real_epsilon_inverse(0.0, 0.0);
  R.diagonal().array() -= 1.0;
  // We add this alpha term to have a nice behaviour around omega = 0. Please be
  // careful with the delta sign delta here is a generic number but in practice
  // it should be e^qp - e^ks
  double erfc_factor = -0.5 * std::copysign(1.0, delta) *
                       std::exp(std::pow(alpha * delta, 2)) *
                       std::erfc(std::abs(alpha * delta));

  double value = ((Imx_row * R.real()).cwiseProduct(Imx_row)).sum();

  return value * erfc_factor;
}

double Sigma_CI::CalcResiduePrefactor(double e_f, double e_m,
                                      double frequency) const {
  double factor = 0.0;

  if (e_f < e_m && e_m < frequency) {
    factor = 1.0;
  } else if (e_f > e_m && e_m > frequency) {
    factor = -1.0;
  } else if (e_m - frequency == 0.0 && e_f > e_m) {
    factor = -0.5;
  } else if (e_m - frequency == 0.0 && e_f < e_m) {
    factor = 0.5;
  }
  return factor;
}
double Sigma_CI::CalcResidueContribution(Eigen::VectorXd rpa_energies,
                                         double frequency,
                                         Index gw_level) const {
  Index rpatotal = rpa_energies.size();
  double sigma_c = 0.0;
  double result_alpha = 0.0;
  Index homo = _opt.homo - _opt.rpamin;
  Index lumo = homo + 1;
  double fermi_rpa = (rpa_energies(lumo) + rpa_energies(homo)) / 2.0;

  const Eigen::MatrixXd& Imx = _Mmn[gw_level];
  for (Index i = 0; i < rpatotal; ++i) {
    double delta = std::abs(rpa_energies(i) - frequency);
    double factor = CalcResiduePrefactor(fermi_rpa, rpa_energies(i), frequency);

    if (factor != 0.0) {
      sigma_c += factor * CalcDiagContribution(Imx.row(i), delta, _eta);
    }

    // This is what is left from the alpha correction in the Quadrature term
    // result_alpha += CalcDiagContributionValue_alpha(Imx.row(i), delta,
    // _opt.alpha);
  }
  return sigma_c;  //+ result_alpha;
}

double Sigma_CI::CalcCorrelationDiagElement(Index gw_level,
                                            double frequency) const {

  const Eigen::VectorXd& RPAenergies = _rpa.getRPAInputEnergies();

  double sigma_c_residue =
      CalcResidueContribution(RPAenergies, frequency, gw_level);

  double sigma_c_integral = _gq.SigmaGQDiag(frequency, gw_level);

  return sigma_c_residue + sigma_c_integral;
}

}  // namespace xtp
}  // namespace votca
