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

#include <cmath>
#include <complex>
#include <iostream>
#include <vector>
#include <votca/xtp/eigen.h>
#include <votca/xtp/gauss_hermite_quadrature_constants.h>
#include <votca/xtp/gaussian_quadrature.h>
#include <votca/xtp/rpa.h>
#include <votca/xtp/threecenter.h>

namespace votca {
namespace xtp {

// Constructor
GaussianQuadrature::GaussianQuadrature(const Eigen::VectorXd& energies,
                                       const TCMatrix_gwbse& Mmn)
    : _energies(energies), _Mmn(Mmn) {}

void GaussianQuadrature::configure(options opt) {
  _opt = opt;
  Gauss_Hermite_Quadrature_Constants ghqc;
  _quadpoints = ghqc.getPoints(_opt.order);
  _quadweights = ghqc.getWeights(_opt.order);
}

Eigen::VectorXd GaussianQuadrature::AdaptedWeights() const {
  // We temporarily move to arrays to enable component-wise operations
  Eigen::ArrayXd quadpoints_sq = (_quadpoints.array()).square();
  Eigen::ArrayXd quadpoints_sq_exp = (quadpoints_sq).exp();
  return (_quadweights.array() * (quadpoints_sq_exp)).matrix();
}

// This function calculates and stores inverses of the microscopic dielectric
// matrix in a matrix vector
std::vector<Eigen::MatrixXcd> GaussianQuadrature::CalcDielInvVector(
    const RPA& rpa) const {
  std::vector<Eigen::MatrixXcd> result;
  for (int j = 0; j < _opt.order; j++) {
    std::complex<double> omega(0, _quadpoints(j));
    result.push_back(rpa.calculate_epsilon(omega).inverse());
  }
  return result;
}

Eigen::VectorXd GaussianQuadrature::SigmaGQDiag(
    const Eigen::VectorXd& frequencies, const RPA& rpa) const {
  Eigen::VectorXd result = Eigen::VectorXd::Zero(_opt.qptotal);
  const std::vector<Eigen::MatrixXcd> DielInvVector = CalcDielInvVector(rpa);
  Eigen::VectorXd shiftedenergies =
      _energies.array() - (_energies(_opt.homo - _opt.rpamin) +
                           _energies(_opt.homo - _opt.rpamin + 1)) /
                              2;
  Eigen::VectorXd AdapWeights = AdaptedWeights();
  double eta = rpa.getEta();  //=1e-4
  int rpatotal = _energies.size();
  int auxsize = _Mmn.auxsize();
  for (int m = 0; m < _opt.qptotal; ++m) {
#if (GWBSE_DOUBLE)
    const Eigen::MatrixXd& Imx = _Mmn[m];
#else
    const Eigen::MatrixXd Imx = _Mmn[m].cast<double>();
#endif
    Eigen::VectorXd DeltaE = frequencies(m) - shiftedenergies.array();
    Eigen::VectorXd DeltaESq = (DeltaE.array()).square();
    for (int j = 0; j < _opt.order; ++j) {
      Eigen::MatrixXd Amx = Eigen::MatrixXd::Zero(rpatotal, auxsize);
      Eigen::MatrixXd Bmx = Eigen::MatrixXd::Zero(rpatotal, auxsize);
      Eigen::MatrixXd Pmx = Eigen::MatrixXd::Zero(rpatotal, auxsize);
      Eigen::MatrixXd Qmx = Eigen::MatrixXd::Zero(rpatotal, auxsize);
      Eigen::MatrixXd RealDielMxInv = Eigen::MatrixXd::Zero(auxsize, auxsize);
      RealDielMxInv = (DielInvVector[j]).real();
      Eigen::MatrixXd ImagDielMxInv = Eigen::MatrixXd::Zero(auxsize, auxsize);
      ImagDielMxInv = (DielInvVector[j]).imag();
      Pmx = Imx * RealDielMxInv - Imx;
      Qmx = Imx * ImagDielMxInv;
      for (int mu = 0; mu < auxsize; ++mu) {
        for (int i = 0; i < _opt.homo - _opt.rpamin + 1; ++i) {
          Amx(i, mu) = DeltaE(i) * Imx(i, mu) /
                       (DeltaESq(i) + std::pow(_quadpoints(j) - eta, 2));
          Bmx(i, mu) = (eta - _quadpoints(j)) * Imx(i, mu) /
                       (DeltaESq(i) + std::pow(_quadpoints(j) - eta, 2));
        }
        for (int i = _opt.homo - _opt.rpamin + 1; i < rpatotal; ++i) {
          Amx(i, mu) = DeltaE(i) * Imx(i, mu) /
                       (DeltaESq(i) + std::pow(_quadpoints(j) + eta, 2));
          Bmx(i, mu) = (-1) * (eta + _quadpoints(j)) * Imx(i, mu) /
                       (DeltaESq(i) + std::pow(_quadpoints(j) + eta, 2));
        }
      }
      result(m) += AdapWeights(j) *
                   ((Pmx.cwiseProduct(Amx) + Qmx.cwiseProduct(Bmx)).sum());
    }
  }
  result /= (-2 * M_PI);
  return result;
}

Eigen::MatrixXd GaussianQuadrature::SigmaGQ(const Eigen::VectorXd& frequencies,
                                            const RPA& rpa) const {
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(_opt.qptotal, _opt.qptotal);
  const std::vector<Eigen::MatrixXcd> DielInvVector = CalcDielInvVector(rpa);
  Eigen::VectorXd shiftedenergies =
      _energies.array() - (_energies(_opt.homo - _opt.rpamin) +
                           _energies(_opt.homo - _opt.rpamin + 1)) /
                              2;
  Eigen::VectorXd AdapWeights = AdaptedWeights();
  double eta = rpa.getEta();  //=1e-4
  int rpatotal = _energies.size();
  int auxsize = _Mmn.auxsize();
  for (int m = 0; m < _opt.qptotal; ++m) {
#if (GWBSE_DOUBLE)
    const Eigen::MatrixXd& Imxm = _Mmn[m];
#else
    const Eigen::MatrixXd Imxm = _Mmn[m].cast<double>();
#endif
    for (int n = 0; n < m; ++n) {
#if (GWBSE_DOUBLE)
      const Eigen::MatrixXd& Imxn = _Mmn[n];
#else
      const Eigen::MatrixXd Imxn = _Mmn[n].cast<double>();
#endif
      Eigen::VectorXd DeltaEm = frequencies(m) - shiftedenergies.array();
      Eigen::VectorXd DeltaEn = frequencies(n) - shiftedenergies.array();
      Eigen::VectorXd DeltaEmSq = (DeltaEm.array()).square();
      Eigen::VectorXd DeltaEnSq = (DeltaEn.array()).square();
      for (int j = 0; j < _opt.order; ++j) {
        Eigen::MatrixXd Amxm = Eigen::MatrixXd::Zero(rpatotal, auxsize);
        Eigen::MatrixXd Amxn = Eigen::MatrixXd::Zero(rpatotal, auxsize);
        Eigen::MatrixXd Bmxm = Eigen::MatrixXd::Zero(rpatotal, auxsize);
        Eigen::MatrixXd Bmxn = Eigen::MatrixXd::Zero(rpatotal, auxsize);
        Eigen::MatrixXd Pmxm = Eigen::MatrixXd::Zero(rpatotal, auxsize);
        Eigen::MatrixXd Pmxn = Eigen::MatrixXd::Zero(rpatotal, auxsize);
        Eigen::MatrixXd Qmxm = Eigen::MatrixXd::Zero(rpatotal, auxsize);
        Eigen::MatrixXd Qmxn = Eigen::MatrixXd::Zero(rpatotal, auxsize);
        Eigen::MatrixXd RealDielMxInv = Eigen::MatrixXd::Zero(auxsize, auxsize);
        RealDielMxInv = (DielInvVector[j]).real();
        Eigen::MatrixXd ImagDielMxInv = Eigen::MatrixXd::Zero(auxsize, auxsize);
        ImagDielMxInv = (DielInvVector[j]).imag();
        Pmxm = Imxm * RealDielMxInv - Imxm;
        Pmxn = Imxn * RealDielMxInv - Imxn;
        Qmxm = Imxm * ImagDielMxInv;
        Qmxn = Imxn * ImagDielMxInv;
        for (int mu = 0; mu < auxsize; ++mu) {
          for (int i = 0; i < _opt.homo - _opt.rpamin + 1; ++i) {
            Amxm(i, mu) = DeltaEm(i) * Imxm(i, mu) /
                          (DeltaEmSq(i) + std::pow(_quadpoints(j) - eta, 2));
            Amxn(i, mu) = DeltaEn(i) * Imxn(i, mu) /
                          (DeltaEnSq(i) + std::pow(_quadpoints(j) - eta, 2));
            Bmxm(i, mu) = (eta - _quadpoints(j)) * Imxm(i, mu) /
                          (DeltaEmSq(i) + std::pow(_quadpoints(j) - eta, 2));
            Bmxn(i, mu) = (eta - _quadpoints(j)) * Imxn(i, mu) /
                          (DeltaEnSq(i) + std::pow(_quadpoints(j) - eta, 2));
          }
          for (int i = _opt.homo - _opt.rpamin + 1; i < rpatotal - 1; ++i) {
            Amxm(i, mu) = DeltaEm(i) * Imxm(i, mu) /
                          (DeltaEmSq(i) + std::pow(_quadpoints(j) + eta, 2));
            Amxn(i, mu) = DeltaEn(i) * Imxn(i, mu) /
                          (DeltaEnSq(i) + std::pow(_quadpoints(j) + eta, 2));
            Bmxm(i, mu) = (-1) * (eta + _quadpoints(j)) * Imxm(i, mu) /
                          (DeltaEmSq(i) + std::pow(_quadpoints(j) + eta, 2));
            Bmxn(i, mu) = (-1) * (eta + _quadpoints(j)) * Imxn(i, mu) /
                          (DeltaEnSq(i) + std::pow(_quadpoints(j) + eta, 2));
          }
        }
        result(m, n) += AdapWeights(j) *
                        ((Pmxn.cwiseProduct(Amxm) + Pmxm.cwiseProduct(Amxn) +
                          Qmxn.cwiseProduct(Bmxm) + Qmxm.cwiseProduct(Bmxn))
                             .sum());
      }
    }
  }
  result += result.transpose();
  result /= (-4 * M_PI);
  result.diagonal() = SigmaGQDiag(frequencies, rpa);
  return result;
}

}  // namespace xtp
}  // namespace votca
