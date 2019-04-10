/*
 *            Copyright 2009-2018 The VOTCA Development Team
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

#include <boost/math/constants/constants.hpp>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>
#include <votca/tools/constants.h>
#include <votca/xtp/eigen.h>
#include <votca/xtp/sigma_ci.h>
#include <votca/xtp/threecenter.h>

namespace votca {
namespace xtp {

void Sigma_CI::PrepareScreening() {
  GaussianQuadrature::options opt;
  opt.homo = _opt.homo;
  opt.order = _opt.order;
  opt.qptotal = _qptotal;
  opt.qpmin = _opt.qpmin;
  opt.rpamin = _opt.rpamin;
  _gq.configure(opt);
}

Eigen::VectorXd Sigma_CI::CalcCorrelationDiag(
    const Eigen::VectorXd& frequencies) const {
  Eigen::VectorXd result = Eigen::VectorXd::Zero(_qptotal);
  const Eigen::VectorXd& energies = _rpa.getRPAInputEnergies();
  const Eigen::VectorXd& shiftedenergies =
      energies.array() - (energies(_opt.homo - _opt.rpamin) +
                          energies(_opt.homo - _opt.rpamin + 1)) /
                             2;
  int rpatotal = energies.size();
  int auxsize = _Mmn.auxsize();
  Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(auxsize, auxsize);
  for (int m = 0; m < _qptotal; ++m) {
    Eigen::MatrixXd Rmx = Eigen::MatrixXd::Zero(rpatotal, auxsize);
#if (GWBSE_DOUBLE)
    const Eigen::MatrixXd& Imx = _Mmn[m];
#else
    const Eigen::MatrixXd Imx = _Mmn[m].cast<double>();
#endif
    for (int i = 0; i < _opt.homo - _opt.rpamin + 1; ++i) {
      std::complex<double> omega(shiftedenergies(i) - frequencies(m), -_eta);
      Eigen::MatrixXd DielMxInv = Eigen::MatrixXd::Zero(auxsize, auxsize);
      Eigen::MatrixXcd DielMxInvC = _rpa.calculate_epsilon(omega).inverse();
      DielMxInv = DielMxInvC.real();
      Eigen::MatrixXd Jmx = Eigen::MatrixXd::Zero(rpatotal, auxsize);
      for (int mu = 0; mu < auxsize; ++mu) {
        for (int nu = 0; nu < auxsize; ++nu) {
          if (frequencies(m) < shiftedenergies(i)) {
            result(m) -=
                Imx(i, mu) * Imx(i, nu) * (DielMxInv(mu, nu) - Id(mu, nu));
          }
        }
      }
      /*if ( frequencies(m) < shiftedenergies(i) ){
          for ( int mu = 0; mu < auxsize; mu++ ){
              Jmx(i,mu) = Imx(i,mu);
         }
      }
      Eigen::MatrixXd JmxTimesDielMxInvmx =
      Eigen::MatrixXd::Zero(rpatotal,auxsize); JmxTimesDielMxInvmx = Jmx *
      DielMxInv; for ( int mu = 0; mu < auxsize; mu++ ){ Rmx(i,mu) =
      JmxTimesDielMxInvmx(i,mu)-Jmx(i,mu); result(m)-= Rmx(i,mu)*Imx(i,mu);
      }*/
    }
    for (int i = _opt.homo - _opt.rpamin + 1; i < rpatotal; ++i) {
      std::complex<double> omega(shiftedenergies(i) - frequencies(m), _eta);
      Eigen::MatrixXd DielMxInv = Eigen::MatrixXd::Zero(auxsize, auxsize);
      Eigen::MatrixXcd DielMxInvC = _rpa.calculate_epsilon(omega).inverse();
      DielMxInv = DielMxInvC.real();
      Eigen::MatrixXd Jmx = Eigen::MatrixXd::Zero(rpatotal, auxsize);
      for (int mu = 0; mu < auxsize; ++mu) {
        for (int nu = 0; nu < auxsize; ++nu) {
          if (frequencies(m) > shiftedenergies(i)) {
            result(m) -=
                Imx(i, mu) * Imx(i, nu) * (DielMxInv(mu, nu) - Id(mu, nu));
          }
        }
      }
      /*if ( frequencies(m) > shiftedenergies(i) ){
          for ( int mu = 0; mu < auxsize; mu++ ){
              Jmx(i,mu) = Imx(i,mu);
         }
      }
      Eigen::MatrixXd JmxTimesDielMxInvmx =
      Eigen::MatrixXd::Zero(rpatotal,auxsize); JmxTimesDielMxInvmx = Jmx *
      DielMxInv; for ( int mu = 0; mu < auxsize; mu++ ){ Rmx(i,mu) =
      JmxTimesDielMxInvmx(i,mu)-Jmx(i,mu); result(m)-=Rmx(i,mu)*Imx(i,mu);
      }*/
    }
    // result(m) = -(Rmx.cwiseProduct(Imx)).sum();
  }
  result += _gq.SigmaGQDiag(frequencies, _rpa);
  return result;
}

/*void Sigma_CI::ExportCorrelationDiags(const Eigen::VectorXd& frequencies)
  const { const int    range = CustomOpts::SigmaExportRange(); const double
  delta = CustomOpts::SigmaExportDelta(); const int    size  = 2 * range + 1;
      CTP_LOG(ctp::logDEBUG, _log)
      << ctp::TimeStamp() << " Writing SigmaC log "
      << "(" << size << ", " << _qptotal << ")"
      << std::flush;
      Eigen::VectorXd offsets = Eigen::VectorXd::LinSpaced(size, -range * delta,
  range * delta); Eigen::MatrixXd results = Eigen::MatrixXd::Zero(size,
  _qptotal); _sigma->PrepareScreening(); for (int i = 0; i < size; i++) {
          results.row(i) =
          _sigma->CalcCorrelationDiag(frequencies +
  Eigen::VectorXd::Constant(_qptotal, offsets[i]));
      }
      Eigen::MatrixXd table = Eigen::MatrixXd::Zero(size + 1, _qptotal + 1);
      table.block(0, 1,    1, _qptotal) = frequencies.transpose();
      table.block(1, 0, size,        1) = offsets;
      table.block(1, 1, size, _qptotal) = results;
      CustomTools::ExportMat("sigma_c.txt", table);
  }*/

Eigen::MatrixXd Sigma_CI::CalcCorrelationOffDiag(
    const Eigen::VectorXd& frequencies) const {
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(_qptotal, _qptotal);
  const Eigen::VectorXd& energies = _rpa.getRPAInputEnergies();
  const Eigen::VectorXd& shiftedenergies =
      energies.array() - (energies(_opt.homo - _opt.rpamin) +
                          energies(_opt.homo - _opt.rpamin + 1)) /
                             2;
  int rpatotal = energies.size();
  int auxsize = _Mmn.auxsize();
  for (int m = 0; m < _qptotal; ++m) {
    Eigen::MatrixXd Rmxm = Eigen::MatrixXd::Zero(rpatotal, auxsize);
#if (GWBSE_DOUBLE)
    const Eigen::MatrixXd& Imxm = _Mmn[m];
#else
    const Eigen::MatrixXd Imxm = _Mmn[m].cast<double>();
#endif
    for (int n = 0; n < m; ++n) {
      Eigen::MatrixXd Rmxn = Eigen::MatrixXd::Zero(rpatotal, auxsize);
#if (GWBSE_DOUBLE)
      const Eigen::MatrixXd& Imxn = _Mmn[n];
#else
      const Eigen::MatrixXd Imxn = _Mmn[n].cast<double>();
#endif
      for (int i = 0; i < _opt.homo - _opt.rpamin + 1; ++i) {
        std::complex<double> omegam(shiftedenergies(i) - frequencies(m), -_eta);
        std::complex<double> omegan(shiftedenergies(i) - frequencies(n), -_eta);
        Eigen::MatrixXd DielMxInvm = Eigen::MatrixXd::Zero(auxsize, auxsize);
        Eigen::MatrixXd DielMxInvn = Eigen::MatrixXd::Zero(auxsize, auxsize);
        Eigen::MatrixXcd DielMxInvmC = _rpa.calculate_epsilon(omegam).inverse();
        DielMxInvm = DielMxInvmC.real();
        Eigen::MatrixXcd DielMxInvnC = _rpa.calculate_epsilon(omegan).inverse();
        DielMxInvn = DielMxInvnC.real();
        Eigen::MatrixXd Jmxm = Eigen::MatrixXd::Zero(rpatotal, auxsize);
        Eigen::MatrixXd Jmxn = Eigen::MatrixXd::Zero(rpatotal, auxsize);
        if (frequencies(m) < shiftedenergies(i)) {
          for (int mu = 0; mu < auxsize; mu++) {
            Jmxm(i, mu) = Imxm(i, mu);
          }
        }
        if (frequencies(n) < shiftedenergies(i)) {
          for (int mu = 0; mu < auxsize; mu++) {
            Jmxn(i, mu) = Imxn(i, mu);
          }
        }
        for (int mu = 0; mu < auxsize; mu++) {
          Rmxm(i, mu) = DielMxInvm(i, mu) - Jmxm(i, mu);
          Rmxn(i, mu) = DielMxInvn(i, mu) - Jmxn(i, mu);
        }
      }
      for (int i = _opt.homo - _opt.rpamin + 1; i < rpatotal; ++i) {
        std::complex<double> omegam(shiftedenergies(i) - frequencies(m), _eta);
        std::complex<double> omegan(shiftedenergies(i) - frequencies(n), _eta);
        Eigen::MatrixXd DielMxInvm = Eigen::MatrixXd::Zero(auxsize, auxsize);
        Eigen::MatrixXd DielMxInvn = Eigen::MatrixXd::Zero(auxsize, auxsize);
        Eigen::MatrixXcd DielMxInvmC = _rpa.calculate_epsilon(omegam).inverse();
        DielMxInvm = DielMxInvmC.real();
        Eigen::MatrixXcd DielMxInvnC = _rpa.calculate_epsilon(omegan).inverse();
        DielMxInvn = DielMxInvnC.real();
        Eigen::MatrixXd Jmxm = Eigen::MatrixXd::Zero(rpatotal, auxsize);
        Eigen::MatrixXd Jmxn = Eigen::MatrixXd::Zero(rpatotal, auxsize);
        if (frequencies(m) > shiftedenergies(i)) {
          for (int mu = 0; mu < auxsize; mu++) {
            Jmxm(i, mu) = Imxm(i, mu);
          }
        }
        if (frequencies(n) > shiftedenergies(i)) {
          for (int mu = 0; mu < auxsize; mu++) {
            Jmxn(i, mu) = Imxn(i, mu);
          }
        }
        for (int mu = 0; mu < auxsize; mu++) {
          Rmxm(i, mu) = DielMxInvm(i, mu) - Jmxm(i, mu);
          Rmxn(i, mu) = DielMxInvn(i, mu) - Jmxn(i, mu);
        }
      }
      result(m, n) =
          (Rmxm.cwiseProduct(Imxn) + Rmxn.cwiseProduct(Imxm)).sum() / (-2);
    }
  }
  result += result.transpose();
  result += _gq.SigmaGQ(frequencies, _rpa);
  result.diagonal() = CalcCorrelationDiag(frequencies);
  return result;
}

}  // namespace xtp
}  // namespace votca
