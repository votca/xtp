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
#include <cmath>
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
  _gq.configure(opt, _rpa);
}

double Sigma_CI::CalcDiagContributionValue(const Eigen::RowVectorXd& Imx_row,
                                           double eta, double delta) const {
  std::complex<double> omega(delta, eta);
  Eigen::MatrixXcd DielMxInvC = (_rpa.calculate_epsilon(omega)).inverse();
  return ((Imx_row * DielMxInvC.real() - Imx_row).cwiseProduct(Imx_row)).sum();
}

double Sigma_CI::CalcOffDiagContributionValue(
    const Eigen::RowVectorXd& Imx_row1, const Eigen::RowVectorXd& Imx_row2,
    double eta, double delta) const {
  std::complex<double> omega(delta, eta);
  Eigen::MatrixXcd DielMxInvC = (_rpa.calculate_epsilon(omega)).inverse();
  return ((Imx_row1 * DielMxInvC.real() - Imx_row1).cwiseProduct(Imx_row2))
      .sum();
}

Eigen::VectorXd Sigma_CI::CalcCorrelationDiag(
    const Eigen::VectorXd& frequencies) const {
  Eigen::VectorXd result = Eigen::VectorXd::Zero(_qptotal);
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
#pragma omp parallel for
  for (int m = 0; m < _qptotal; ++m) {
    const Eigen::MatrixXd& Imx = _Mmn[m];
    for (int i = 0; i < rpatotal; ++i) {
      double delta = shiftedenergies(i) - shiftedfrequencies(m);
      double value = 0;
      /* std::cout << " " << std::endl;
      std::cout << "energies" << std::endl;
      std::cout << " " << std::endl;
      std::cout << shiftedenergies << std::endl;
      std::cout << " " << std::endl;
      std::cout << "frequencies" << std::endl;
      std::cout << " " << std::endl;
      std::cout << shiftedfrequencies << std::endl;
      std::cout << " " << std::endl;*/

      if (delta > 0 && i < lumo) {
        value = CalcDiagContributionValue(Imx.row(i), _eta, delta);
      } else if (delta < 0 && i > homo) {
        value = - CalcDiagContributionValue(Imx.row(i), -_eta, delta);
      }
      result(m) -= value;
    }
  }
   result += _gq.SigmaGQDiag(frequencies, _rpa);
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
