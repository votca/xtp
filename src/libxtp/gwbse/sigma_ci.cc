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

#include "votca/xtp/sigma_ppm.h"
#include <c++/5/complex>
#include <cmath>
#include <time.h>
#include <votca/tools/constants.h>
#include <votca/xtp/customtools.h>
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
  opt.alpha = _opt.alphaa;
  _gq.configure(opt, _rpa);
}

double Sigma_CI::CalcDiagContributionValue(const Eigen::RowVectorXd& Imx_row,
                                           double delta, double eta) const {

  // This function is used in the calculation of the residues
  // This calculates eps^-1 (inverse of the dielectric function) for complex
  // frequencies of the kind omega = delta + i*eta
  Eigen::MatrixXd DielMxInv = _rpa.calculate_real_epsilon_inverse(delta, eta);
  // I subract the Identity to obtain what is usually called the susceptibility
  DielMxInv -= Eigen::MatrixXd::Identity(DielMxInv.rows(), DielMxInv.cols());
  // This evaluate the diagonal contribution for a specific dielectric function
  // at a generic complex frequency omega (What in the documentation is called
  // Lambda
  double value = ((Imx_row * DielMxInv).cwiseProduct(Imx_row)).sum();

  return value;
}

double Sigma_CI::CalcDiagContributionValue_i(const Eigen::RowVectorXd& Imx_row,
                                             double delta, double eta) const {

  // This function is used in the calculation of the residues
  // This calculates eps^-1 (inverse of the dielectric function) for complex
  // frequencies of the kind omega = delta + i*eta
  Eigen::MatrixXd DielMxInv = _rpa.calculate_imag_epsilon_inverse(delta, eta);
  // I subract the Identity to obtain what is usually called the susceptibility
  DielMxInv -= Eigen::MatrixXd::Identity(DielMxInv.rows(), DielMxInv.cols());
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
  R -= Eigen::MatrixXd::Identity(R.rows(), R.cols());
  // We add this alpha term to have a nice behaviour around omega = 0. Please be
  // careful with the delta sign delta here is a generic number but in practice
  // it should be e^qp - e^ks

  // double erfc_factor = 0.5 * std::exp(std::pow(_opt.alpha*delta,2))*
  // std::erfc(_opt.alpha*delta);

  double erfc_factor = -0.5 * std::copysign(1.0, delta) *
                       std::exp(std::pow(alpha * delta, 2)) *
                       std::erfc(std::abs(alpha * delta));

  R *= erfc_factor;

  double value = ((Imx_row * R).cwiseProduct(Imx_row)).sum();

  return value;
}

double Sigma_CI::CalcDiagContributionValue_alpha_i(
    const Eigen::RowVectorXd& Imx_row, double delta, double alpha) const {
  // This function is used in the calculation of the residues
  // This calculates eps^-1 (inverse of the dielectric function) for complex
  // frequency omega = 0 + i* 0 (origin)
  Eigen::MatrixXd R = _rpa.calculate_imag_epsilon_inverse(0.0, 0.0);
  // I subract the Identity to obtain what is usually called the susceptibility
  R -= Eigen::MatrixXd::Identity(R.rows(), R.cols());
  // We add this alpha term to have a nice behaviour around omega = 0. Please be
  // careful with the delta sign delta here is a generic number but in practice
  // it should be e^qp - e^ks

  // double erfc_factor = 0.5 * std::exp(std::pow(_opt.alpha*delta,2))*
  // std::erfc(_opt.alpha*delta);

  double erfc_factor = -0.5 * std::copysign(1.0, delta) *
                       std::exp(std::pow(alpha * delta, 2)) *
                       std::erfc(std::abs(alpha * delta));

  R *= erfc_factor;

  double value = ((Imx_row * R).cwiseProduct(Imx_row)).sum();

  return value;
}

double Sigma_CI::CalcOffDiagContributionValue(
    const Eigen::RowVectorXd& Imx_row1, const Eigen::RowVectorXd& Imx_row2,
    double eta, double delta) const {
  std::complex<double> omega(delta, eta);
  Eigen::MatrixXcd DielMxInvC = (_rpa.calculate_epsilon(omega)).inverse();
  return ((Imx_row1 * DielMxInvC.real() - Imx_row1).cwiseProduct(Imx_row2))
      .sum();
}

// Eigen::VectorXd Sigma_CI::CalcCorrelationDiag(const Eigen::VectorXd& E) const
// {

//   int homo = _opt.homo - _opt.rpamin;
//   int lumo = homo + 1;

//   const Eigen::VectorXd& RPAenergies = _rpa.getRPAInputEnergies();
//   double fermi_rpa = (RPAenergies(lumo) + RPAenergies(homo)) / 2;

//   int rpatotal = RPAenergies.size();

//   Eigen::VectorXd result = Eigen::VectorXd::Zero(_qptotal);
//   Eigen::VectorXd E_arr =  E.array() ;

//   Eigen::VectorXd result_alpha = Eigen::VectorXd::Zero(_qptotal);
//   Eigen::VectorXd result_occ = Eigen::VectorXd::Zero(_qptotal);
//   Eigen::VectorXd result_unocc = Eigen::VectorXd::Zero(_qptotal);

//  std::cout << " INTO RESIDUE Calculations" << std::endl;

// #pragma omp parallel for
//                 for (int j = 0; j < _qptotal; ++j) {
//                     if ( j%10 == 0 ){
//                         std::cout << " Level \t " << j << "Correction " <<
//                         std::endl;
//                     }
//                     const Eigen::MatrixXd& Imx = _Mmn[j];
//                     //This loop stands for the contribution of the residual
//                     part for (int i = 0; i < rpatotal; ++i) {
//                         // delta_ji = E^qp(j) - E^{KS}(i)
//                         double delta = - RPAenergies(i) + E_arr(j);
//                         double value = 0;
//                         double value_occ = 0;
//                         double value_unocc = 0;
//                         double value_alpha = 0;
//                         double discriminant = fermi_rpa-RPAenergies(i);
//                         if (delta > 0 && discriminant < 0 ) {
//                             value =   CalcDiagContributionValue(Imx.row(i),
//                             delta, _eta);
//                         }
//                         else if (delta < 0 && discriminant > 0 ) {
//                             value =  - CalcDiagContributionValue(Imx.row(i),
//                             delta, -_eta);
//                         }
//                         else if (delta == 0 && discriminant > 0) {
//                             value_occ =
//                             -0.5*CalcDiagContributionValue(Imx.row(i), 0,
//                             -_eta);
//                         }
//                         else if (delta == 0 && discriminant < 0 ) {
//                             value_unocc =
//                             0.5*CalcDiagContributionValue(Imx.row(i), 0,
//                             _eta);
//                         }
//                         result(j) += value ;
//                         result_occ(j) +=value_occ;
//                         result_unocc(j) +=value_unocc;
//                         // This is what is left from the alpha correction in
//                         the Quadrature term result_alpha(j) +=
//                         CalcDiagContributionValue_alpha(Imx.row(i), delta,
//                         _opt.alphaa );
//                     }
//                 }

//  std::cout << " INTO GAUSSQ Calculation" << std::endl;
// Eigen::VectorXd GAUSSQ = _gq.SigmaGQDiag(E_arr, _rpa);

// Eigen::MatrixXd niceresults = Eigen::MatrixXd::Zero(result.rows(),5);
// niceresults.col(0) = result;
// niceresults.col(1) = result_occ;
// niceresults.col(2) = result_unocc;
// niceresults.col(3) = result_alpha;
// niceresults.col(4) = GAUSSQ;

// std::cout << "\n" << "Res Res_0_occ Res_0_unocc Res_alpha QD " << "\n" <<
// std::endl; std::cout << niceresults << std::endl;

// // update E_arr
// result += GAUSSQ  + result_occ + result_unocc + result_alpha;

// return result;
// }

Eigen::VectorXd Sigma_CI::CalcCorrelationDiag(const Eigen::VectorXd& E) const {

  int homo = _opt.homo - _opt.rpamin;
  int lumo = homo + 1;

  const Eigen::VectorXd& RPAenergies = _rpa.getRPAInputEnergies();
  double fermi_rpa = (RPAenergies(lumo) + RPAenergies(homo)) / 2;

  int rpatotal = RPAenergies.size();

  Eigen::VectorXd result = Eigen::VectorXd::Zero(_qptotal);
  Eigen::VectorXd E_arr = E.array();

  Eigen::VectorXd result_alpha = Eigen::VectorXd::Zero(_qptotal);
  Eigen::VectorXd result_occ = Eigen::VectorXd::Zero(_qptotal);
  Eigen::VectorXd result_unocc = Eigen::VectorXd::Zero(_qptotal);

  std::cout << " INTO RESIDUE Calculations" << std::endl;

#pragma omp parallel for
  for (int j = 0; j < _qptotal; ++j) {
    if (j % 10 == 0) {
      std::cout << " Level \t " << j << "Correction " << std::endl;
    }
    const Eigen::MatrixXd& Imx = _Mmn[j];
    // This loop stands for the contribution of the residual part
    for (int i = 0; i < rpatotal; ++i) {
      // delta_ji = E^qp(j) - E^{KS}(i)
      double delta = -RPAenergies(i) + E_arr(j);
      double value = 0;
      double value_occ = 0;
      double value_unocc = 0;
      double value_alpha = 0;
      double discriminant = fermi_rpa - RPAenergies(i);
      if (delta > 0 && discriminant < 0) {
        value = CalcDiagContributionValue(Imx.row(i), delta, _eta);
      } else if (delta < 0 && discriminant > 0) {
        value = -CalcDiagContributionValue(Imx.row(i), delta, -_eta);
      } else if (delta == 0 && discriminant > 0) {
        value_occ = -0.5 * CalcDiagContributionValue(Imx.row(i), 0, -_eta);
      } else if (delta == 0 && discriminant < 0) {
        value_unocc = 0.5 * CalcDiagContributionValue(Imx.row(i), 0, _eta);
      }
      result(j) += value;
      result_occ(j) += value_occ;
      result_unocc(j) += value_unocc;
      // This is what is left from the alpha correction in the Quadrature term
      result_alpha(j) +=
          CalcDiagContributionValue_alpha(Imx.row(i), delta, _opt.alphaa);
    }
  }

  std::cout << " INTO GAUSSQ Calculation" << std::endl;
  Eigen::VectorXd GAUSSQ = _gq.SigmaGQDiag(E_arr, _rpa);

  Eigen::MatrixXd niceresults = Eigen::MatrixXd::Zero(result.rows(), 5);
  niceresults.col(0) = result;
  niceresults.col(1) = result_occ;
  niceresults.col(2) = result_unocc;
  niceresults.col(3) = result_alpha;
  niceresults.col(4) = GAUSSQ;

  std::cout << "\n"
            << "Res Res_0_occ Res_0_unocc Res_alpha QD "
            << "\n"
            << std::endl;
  std::cout << niceresults << std::endl;

  // update E_arr
  result += GAUSSQ + result_occ + result_unocc + result_alpha;

  return result;
}

Eigen::VectorXd Sigma_CI::CalcCorrelationDiag_imag(
    const Eigen::VectorXd& E) const {

  int homo = _opt.homo - _opt.rpamin;
  int lumo = homo + 1;

  const Eigen::VectorXd& RPAenergies = _rpa.getRPAInputEnergies();
  double fermi_rpa = (RPAenergies(lumo) + RPAenergies(homo)) / 2;

  int rpatotal = RPAenergies.size();

  Eigen::VectorXd result = Eigen::VectorXd::Zero(_qptotal);
  Eigen::VectorXd E_arr = E.array();

  Eigen::VectorXd result_alpha = Eigen::VectorXd::Zero(_qptotal);
  Eigen::VectorXd result_occ = Eigen::VectorXd::Zero(_qptotal);
  Eigen::VectorXd result_unocc = Eigen::VectorXd::Zero(_qptotal);

  std::cout << " INTO RESIDUE Calculations" << std::endl;

#pragma omp parallel for
  for (int j = 0; j < _qptotal; ++j) {
    if (j % 10 == 0) {
      std::cout << " Level \t " << j << "Correction " << std::endl;
    }
    const Eigen::MatrixXd& Imx = _Mmn[j];
    // This loop stands for the contribution of the residual part
    for (int i = 0; i < rpatotal; ++i) {
      // delta_ji = E^qp(j) - E^{KS}(i)
      double delta = -RPAenergies(i) + E_arr(j);
      double value = 0;
      double value_occ = 0;
      double value_unocc = 0;
      double value_alpha = 0;
      double discriminant = fermi_rpa - RPAenergies(i);
      if (delta > 0 && discriminant < 0) {
        value = CalcDiagContributionValue_i(Imx.row(i), delta, _eta);
      } else if (delta < 0 && discriminant > 0) {
        value = -CalcDiagContributionValue_i(Imx.row(i), delta, -_eta);
      } else if (delta == 0 && discriminant > 0) {
        value_occ = -0.5 * CalcDiagContributionValue_i(Imx.row(i), 0, -_eta);
      } else if (delta == 0 && discriminant < 0) {
        value_unocc = 0.5 * CalcDiagContributionValue_i(Imx.row(i), 0, _eta);
      }
      result(j) += value;
      result_occ(j) += value_occ;
      result_unocc(j) += value_unocc;
      // This is what is left from the alpha correction in the Quadrature term
      result_alpha(j) +=
          CalcDiagContributionValue_alpha_i(Imx.row(i), delta, _opt.alphaa);
    }
  }

  std::cout << " INTO GAUSSQ Calculation" << std::endl;
  Eigen::VectorXd GAUSSQ = _gq.SigmaGQDiag_i(E_arr, _rpa);

  Eigen::MatrixXd niceresults = Eigen::MatrixXd::Zero(result.rows(), 5);
  niceresults.col(0) = result;
  niceresults.col(1) = result_occ;
  niceresults.col(2) = result_unocc;
  niceresults.col(3) = result_alpha;
  niceresults.col(4) = GAUSSQ;

  std::cout << "\n"
            << "Res Res_0_occ Res_0_unocc Res_alpha QD "
            << "\n"
            << std::endl;
  std::cout << niceresults << std::endl;

  // update E_arr
  result += GAUSSQ + result_occ + result_unocc + result_alpha;

  return result;
}

Eigen::MatrixXd Sigma_CI::CalcCorrelationOffDiag(
    const Eigen::VectorXd& frequencies) const {
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(_qptotal, _qptotal);
  int homo = _opt.homo - _opt.rpamin;
  int lumo = homo + 1;
  const Eigen::VectorXd& energies = _rpa.getRPAInputEnergies();
  double middle_of_gap_energy = (energies(lumo) + energies(homo)) / 2;
  const Eigen::VectorXd shiftedenergies =
      energies.array() - middle_of_gap_energy;
  const Eigen::VectorXd shiftedfrequencies =
      frequencies.array() - middle_of_gap_energy;
  int rpatotal = energies.size();
  int auxsize = _Mmn.auxsize();
#pragma omp parallel for schedule(dynamic)
  for (int m = 0; m < _qptotal; ++m) {
    Eigen::MatrixXd Rmxm = Eigen::MatrixXd::Zero(rpatotal, auxsize);
    const Eigen::MatrixXd& Imxm = _Mmn[m];
    for (int n = m + 1; n < _qptotal; ++n) {
      Eigen::MatrixXd Rmxn = Eigen::MatrixXd::Zero(rpatotal, auxsize);
      const Eigen::MatrixXd& Imxn = _Mmn[n];
      double resultmn = 0.0;
      for (int i = 0; i < rpatotal; ++i) {
        double deltam = shiftedenergies(i) - shiftedfrequencies(m);
        double deltan = shiftedenergies(i) - shiftedfrequencies(n);
        double value = 0;
        if (i < lumo) {
          if (deltam > 0) {
            value += CalcOffDiagContributionValue(Imxm.row(i), Imxn.row(i),
                                                  _eta, deltam);
          }
          if (deltan > 0) {
            value += CalcOffDiagContributionValue(Imxn.row(i), Imxm.row(i),
                                                  _eta, deltan);
          }
        } else if (i > homo) {
          if (deltam < 0) {
            value += CalcOffDiagContributionValue(Imxm.row(i), Imxn.row(i),
                                                  -_eta, deltam);
          }
          if (deltan < 0) {
            value += CalcOffDiagContributionValue(Imxn.row(i), Imxm.row(i),
                                                  -_eta, deltan);
          }
        }
        resultmn -= value / 2;
      }
      result(m, n) = resultmn;
      result(n, m) = resultmn;
    }
  }
  result += _gq.SigmaGQ(frequencies, _rpa);
  return result;
}

Eigen::VectorXd Sigma_CI::ExactCorrelationDiag(
    const Eigen::VectorXd& frequencies) const {
  Eigen::VectorXd result = Eigen::VectorXd::Zero(_qptotal, _qptotal);
  Eigen::VectorXcd complexresult = Eigen::VectorXcd::Zero(_qptotal, _qptotal);
  const Eigen::VectorXd& energies = _rpa.getRPAInputEnergies();
  const Eigen::VectorXd& shiftedenergies =
      energies.array() - (energies(_opt.homo - _opt.rpamin) +
                          energies(_opt.homo - _opt.rpamin + 1)) /
                             2;
  const Eigen::VectorXd& shiftedfrequencies =
      frequencies.array() - (energies(_opt.homo - _opt.rpamin) +
                             energies(_opt.homo - _opt.rpamin + 1)) /
                                2;
  int rpatotal = energies.size();
  int auxsize = _Mmn.auxsize();
  Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(auxsize, auxsize);
  for (int m = 0; m < _qptotal; ++m) {
    const Eigen::MatrixXd Imxm = _Mmn[m].cast<double>();
    for (int i = 0; i < _opt.homo - _opt.rpamin + 1; ++i) {
      std::complex<double> omegam(shiftedenergies(i) - shiftedfrequencies(m),
                                  _eta);
      Eigen::MatrixXcd DielMxInvm = Eigen::MatrixXcd::Zero(auxsize, auxsize);
      DielMxInvm = _rpa.calculate_epsilon(omegam).inverse();
      if (shiftedfrequencies(m) < shiftedenergies(i)) {
        for (int mu = 0; mu < auxsize; ++mu) {
          for (int nu = 0; nu < auxsize; ++nu) {
            complexresult(m) -=
                Imxm(i, mu) * Imxm(i, nu) * (DielMxInvm(mu, nu) - Id(mu, nu));
          }
        }
      }
    }
    for (int i = _opt.homo - _opt.rpamin + 1; i < rpatotal; ++i) {
      std::complex<double> omegam(shiftedenergies(i) - shiftedfrequencies(m),
                                  -_eta);
      Eigen::MatrixXcd DielMxInvm = Eigen::MatrixXcd::Zero(auxsize, auxsize);
      DielMxInvm = _rpa.calculate_epsilon(omegam).inverse();
      if (shiftedfrequencies(m) > shiftedenergies(i)) {
        for (int mu = 0; mu < auxsize; ++mu) {
          for (int nu = 0; nu < auxsize; ++nu) {
            complexresult(m) -=
                Imxm(i, mu) * Imxm(i, nu) * (DielMxInvm(mu, nu) - Id(mu, nu));
          }
        }
      }
    }
  }
  result = complexresult.real();
  result += _gq.ExactSigmaGQDiag(frequencies, _rpa);
  return result;
}

Eigen::MatrixXd Sigma_CI::ExactCorrelationOffDiag(
    const Eigen::VectorXd& frequencies) const {
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(_qptotal, _qptotal);
  Eigen::MatrixXcd complexresult = Eigen::MatrixXcd::Zero(_qptotal, _qptotal);
  const Eigen::VectorXd& energies = _rpa.getRPAInputEnergies();
  const Eigen::VectorXd& shiftedenergies =
      energies.array() - (energies(_opt.homo - _opt.rpamin) +
                          energies(_opt.homo - _opt.rpamin + 1)) /
                             2;
  const Eigen::VectorXd& shiftedfrequencies =
      frequencies.array() - (energies(_opt.homo - _opt.rpamin) +
                             energies(_opt.homo - _opt.rpamin + 1)) /
                                2;
  int rpatotal = energies.size();
  int auxsize = _Mmn.auxsize();
  Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(auxsize, auxsize);
  for (int m = 0; m < _qptotal; ++m) {
    const Eigen::MatrixXd Imxm = _Mmn[m].cast<double>();
    for (int n = 0; n < _qptotal; ++n) {
      const Eigen::MatrixXd Imxn = _Mmn[n].cast<double>();
      for (int i = 0; i < _opt.homo - _opt.rpamin + 1; ++i) {
        std::complex<double> omegam(shiftedenergies(i) - shiftedfrequencies(m),
                                    _eta);
        std::complex<double> omegan(shiftedenergies(i) - shiftedfrequencies(n),
                                    _eta);
        Eigen::MatrixXcd DielMxInvm = Eigen::MatrixXcd::Zero(auxsize, auxsize);
        Eigen::MatrixXcd DielMxInvn = Eigen::MatrixXcd::Zero(auxsize, auxsize);
        DielMxInvm = _rpa.calculate_epsilon(omegam).inverse();
        DielMxInvn = _rpa.calculate_epsilon(omegan).inverse();
        if (shiftedfrequencies(m) < shiftedenergies(i)) {
          for (int mu = 0; mu < auxsize; ++mu) {
            for (int nu = 0; nu < auxsize; ++nu) {
              complexresult(m, n) -=
                  Imxm(i, mu) * Imxn(i, nu) * (DielMxInvm(mu, nu) - Id(mu, nu));
            }
          }
        }
        if (shiftedfrequencies(n) < shiftedenergies(i)) {
          for (int mu = 0; mu < auxsize; ++mu) {
            for (int nu = 0; nu < auxsize; ++nu) {
              complexresult(m, n) -=
                  Imxm(i, mu) * Imxn(i, nu) * (DielMxInvn(mu, nu) - Id(mu, nu));
            }
          }
        }
      }
      for (int i = _opt.homo - _opt.rpamin + 1; i < rpatotal; ++i) {
        std::complex<double> omegam(shiftedenergies(i) - shiftedfrequencies(m),
                                    -_eta);
        std::complex<double> omegan(shiftedenergies(i) - shiftedfrequencies(n),
                                    -_eta);
        Eigen::MatrixXcd DielMxInvm = Eigen::MatrixXcd::Zero(auxsize, auxsize);
        Eigen::MatrixXcd DielMxInvn = Eigen::MatrixXcd::Zero(auxsize, auxsize);
        DielMxInvm = _rpa.calculate_epsilon(omegam).inverse();
        DielMxInvn = _rpa.calculate_epsilon(omegan).inverse();
        if (shiftedfrequencies(m) > shiftedenergies(i)) {
          for (int mu = 0; mu < auxsize; ++mu) {
            for (int nu = 0; nu < auxsize; ++nu) {
              complexresult(m, n) -=
                  Imxm(i, mu) * Imxn(i, nu) * (DielMxInvm(mu, nu) - Id(mu, nu));
            }
          }
        }
        if (shiftedfrequencies(n) > shiftedenergies(i)) {
          for (int mu = 0; mu < auxsize; ++mu) {
            for (int nu = 0; nu < auxsize; ++nu) {
              complexresult(m, n) -=
                  Imxm(i, mu) * Imxn(i, nu) * (DielMxInvn(mu, nu) - Id(mu, nu));
            }
          }
        }
      }
    }
  }
  result = complexresult.real();
  result /= 2;
  result += _gq.ExactSigmaGQOffDiag(frequencies, _rpa);
  return result;
}

}  // namespace xtp
}  // namespace votca
