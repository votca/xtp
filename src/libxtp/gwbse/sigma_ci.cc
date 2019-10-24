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
#include <time.h>

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
  
  
  Eigen::MatrixXcd EPSQ = _gq.EpsilonGQ(omega, _rpa);
  
 
  
  
  Eigen::MatrixXcd DielMxInvC = (_rpa.calculate_epsilon(omega)).inverse();
  //// Wouter's implementation
  double value = (( Imx_row * DielMxInvC.real() - Imx_row).cwiseProduct(Imx_row)).sum();
  
  /// Gianluca's implementation
  //Eigen::MatrixXcd M = DielMxInvC - Eigen::MatrixXd::Identity(DielMxInvC.rows(),DielMxInvC.cols()); 
  //double value = (Imx_row * M.real() * Imx_row.transpose()).trace();
    
  return value;
}

double Sigma_CI::CalcOffDiagContributionValue(
    const Eigen::RowVectorXd& Imx_row1, const Eigen::RowVectorXd& Imx_row2,
    double eta, double delta) const {
  std::complex<double> omega(delta, eta);
  Eigen::MatrixXcd DielMxInvC = (_rpa.calculate_epsilon(omega)).inverse();
  return ((Imx_row1 * DielMxInvC.real() - Imx_row1).cwiseProduct(Imx_row2)).sum();
}


Eigen::VectorXd Sigma_CI::CalcCorrelationDiag(
    const Eigen::VectorXd& E) const {
 
  int homo = _opt.homo - _opt.rpamin;
  int lumo = homo + 1;
  std::cout << " Using eta " << std::setprecision(10) << _eta << std::endl;
  
  const Eigen::VectorXd& RPAenergies = _rpa.getRPAInputEnergies();
  double fermi = (RPAenergies(lumo) + RPAenergies(homo)) / 2;
  
  int rpatotal = RPAenergies.size();
  
  // const Eigen::VectorXd shiftedenergies = energies.array(); // - middle_of_gap_energy;
  Eigen::VectorXd shiftedE =  E.array() ; 
  
  // Open the gap for iteration 0
  
 
  Eigen::VectorXd result = Eigen::VectorXd::Zero(_qptotal);
                result = Eigen::VectorXd::Zero(_qptotal);
                Eigen::VectorXd res0p = Eigen::VectorXd::Zero(_qptotal);
                Eigen::VectorXd res0m = Eigen::VectorXd::Zero(_qptotal);

std::cout << " INTO RESIDUE Calculations" << std::endl;
#pragma omp parallel for
                for (int m = 0; m < _qptotal; ++m) {
                    if ( m%10 == 0 ){
                        std::cout << " Level \t " << m << "Correction " << std::endl; 
                    }
                    const Eigen::MatrixXd& Imx = _Mmn[m];
                    for (int i = 0; i < rpatotal; ++i) {
                        double delta = RPAenergies(i) - shiftedE(m);
			//double abs_delta = std::abs(delta);
                        double value = 0;
                        double value_res0p = 0;
                        double value_res0m = 0;

                        if (delta > 0 && RPAenergies(i) < fermi) {
                            // This is for occupied levels
                            value = CalcDiagContributionValue(Imx.row(i), _eta, delta);
                        } else if (delta == 0 && RPAenergies(i) < fermi) {
                            // This is for occupied levels
                            value_res0p = 0.5 * CalcDiagContributionValue(Imx.row(i), _eta, delta);
                            //std::cout << " I should not be in this p0 case!" << std::endl;
                        } else if (delta < 0 && RPAenergies(i) > fermi) {
                            value = -CalcDiagContributionValue(Imx.row(i), -_eta, delta);
                        } else if (delta == 0 && RPAenergies(i) > fermi) {
                            value_res0m = 0.5 * CalcDiagContributionValue(Imx.row(i), -_eta, delta);
                            //std::cout << " I should not be in this m0 case!" << std::endl;
                        }
                        res0p(m) -= value_res0p;
                        res0m(m) -= value_res0m;
                        result(m) -= value;
                    }
                }

                std::cout << " INTO GAUSSQ Calculation" << std::endl;
                Eigen::VectorXd GAUSSQ = _gq.SigmaGQDiag(shiftedE, _rpa);
                
                result += GAUSSQ + res0p + res0m;

                // update shiftedE 
                shiftedE = E + result;


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
