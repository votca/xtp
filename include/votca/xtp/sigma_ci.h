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

#ifndef _VOTCA_XTP_SIGMA_CI_H
#define _VOTCA_XTP_SIGMA_CI_H
#include <complex>
#include <votca/xtp/eigen.h>
#include <votca/xtp/gaussian_quadrature.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/rpa.h>
#include <votca/xtp/sigma_base.h>

// This computes the whole expectation matrix for the correlational part of the
// self-energy: so, both the residual and the Gauss-Hermite quadrature
// contribution
namespace votca {
namespace xtp {

class TCMatrix_gwbse;
class RPA;

class Sigma_CI : public Sigma_base {

 public:
  Sigma_CI(TCMatrix_gwbse& Mmn, RPA& rpa, Eigen::MatrixXd vxc)
      : Sigma_base(Mmn, rpa),
        _gq(rpa.getRPAInputEnergies(), Mmn),
        _eta(rpa.getEta()),
        _vxc(vxc){};

  ~Sigma_CI(){};

  void PrepareScreening();

  Eigen::VectorXd CalcCorrelationDiag(const Eigen::VectorXd& frequencies) const;
  Eigen::VectorXd CalcCorrelationDiag_imag(
      const Eigen::VectorXd& frequencies) const;

  Eigen::VectorXd ExactCorrelationDiag(
      const Eigen::VectorXd& frequencies) const;

  Eigen::MatrixXd CalcCorrelationOffDiag(
      const Eigen::VectorXd& frequencies) const;

  Eigen::MatrixXd ExactCorrelationOffDiag(
      const Eigen::VectorXd& frequencies) const;

  void SetSigmaX(Eigen::MatrixXd sigmax) { _sigmaX = sigmax; };

 private:
  double CalcDiagContributionValue(const Eigen::RowVectorXd& Imx_row,
                                   double delta, double eta) const;
  double CalcDiagContributionValue_i(const Eigen::RowVectorXd& Imx_row,
                                     double delta, double eta) const;

  /*double CalcDiagContributionValue(const Eigen::MatrixXd& IMatrix, double eta,
                                   double delta) const;*/

  double CalcOffDiagContributionValue(const Eigen::RowVectorXd& Imx_row1,
                                      const Eigen::RowVectorXd& Imx_row2,
                                      double eta, double delta) const;

  double CalcDiagContributionValue_alpha(const Eigen::RowVectorXd& Imx_row,
                                         double delta, double alpha) const;
  double CalcDiagContributionValue_alpha_i(const Eigen::RowVectorXd& Imx_row,
                                           double delta, double alpha) const;
  GaussianQuadrature _gq;

  double _eta;

  // double _alpha;

  Eigen::MatrixXd _vxc;
  Eigen::MatrixXd _sigmaX;
};

}  // namespace xtp

}  // namespace votca

#endif /* _VOTCA_XTP_SIGMA_CI_H */
