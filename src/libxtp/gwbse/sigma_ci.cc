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

double Sigma_CI::CalcDiagContribution(const Eigen::RowVectorXd& Imx_row,
                                      double delta, double eta) const {

  // This function is used in the calculation of the residues
  // This calculates eps^-1 (inverse of the dielectric function) for complex
  // frequencies of the kind omega = delta + i*eta
  Eigen::MatrixXd DielMxInv = _rpa.calculate_real_epsilon_inverse(delta, eta);
  // I subract the Identity to obtain what is usually called the susceptibility
  DielMxInv.diagonal().array() -= 1.0;
  // This evaluate the diagonal contribution for a specific dielectric function
  // at a generic complex frequency omega (What in the documentation is called
  // Lambda
  double value = ((Imx_row * DielMxInv).cwiseProduct(Imx_row)).sum();

  return value;
}

double Sigma_CI::CalcDiagContributionValue_alpha(
    const Eigen::RowVectorXd& Imx_row, double delta, double alpha) const {
  // This function is used in the calculation of the residues
  // This calculates eps^-1 (inverse of the dielectric function) for complex
  // frequency omega = 0 + i* 0 (origin)
  Eigen::MatrixXd R = _rpa.calculate_real_epsilon_inverse(0.0, 0.0);
  // I subract the Identity to obtain what is usually called the susceptibility
  R.diagonal().array() -= 1.0;
  // We add this alpha term to have a nice behaviour around omega = 0. Please be
  // careful with the delta sign delta here is a generic number but in practice
  // it should be e^qp - e^ks

  double erfc_factor = -0.5 * std::copysign(1.0, delta) *
                       std::exp(std::pow(alpha * delta, 2)) *
                       std::erfc(std::abs(alpha * delta));

  R *= erfc_factor;

  double value = ((Imx_row * R).cwiseProduct(Imx_row)).sum();

  return value;
}

double Sigma_CI::CalcCorrelationDiagElement(Index gw_level,
                                            double frequency) const {

  Index homo = _opt.homo - _opt.rpamin;
  Index lumo = homo + 1;

  const Eigen::VectorXd& RPAenergies = _rpa.getRPAInputEnergies();
  double fermi_rpa = (RPAenergies(lumo) + RPAenergies(homo)) / 2;

  Index rpatotal = RPAenergies.size();

  double sigma_c = 0.0;

  double result_alpha = 0.0;
  double result_occ = 0.0;
  double result_unocc = 0.0;

  const Eigen::MatrixXd& Imx = _Mmn[gw_level];

  for (Index i = 0; i < rpatotal; ++i) {
    // delta_ji = E^qp(j) - E^{KS}(i)
    double delta = -RPAenergies(i) + frequency;
    double value = 0;
    double value_occ = 0;
    double value_unocc = 0;
    double discriminant = fermi_rpa - RPAenergies(i);
    if (delta > 0 && discriminant < 0) {
      value = CalcDiagContribution(Imx.row(i), delta, _eta);
    } else if (delta < 0 && discriminant > 0) {
      value = -CalcDiagContribution(Imx.row(i), delta, -_eta);
    } else if (delta == 0 && discriminant > 0) {
      value_occ = -0.5 * CalcDiagContribution(Imx.row(i), 0, -_eta);
    } else if (delta == 0 && discriminant < 0) {
      value_unocc = 0.5 * CalcDiagContribution(Imx.row(i), 0, _eta);
    }
    sigma_c += value;
    result_occ += value_occ;
    result_unocc += value_unocc;
    // This is what is left from the alpha correction in the Quadrature term
    result_alpha +=
        CalcDiagContributionValue_alpha(Imx.row(i), delta, _opt.alpha);
  }

  double GAUSSQ = _gq.SigmaGQDiag(frequency, gw_level);

  sigma_c += GAUSSQ + result_occ + result_unocc + result_alpha;

  return sigma_c;
}

}  // namespace xtp
}  // namespace votca
