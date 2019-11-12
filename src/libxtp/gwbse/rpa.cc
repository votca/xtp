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

#include "votca/xtp/customtools.h"
#include "votca/xtp/threecenter.h"
#include "votca/xtp/vc2index.h"
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/rpa.h>

using std::flush;

namespace votca {
namespace xtp {

void RPA::UpdateRPAInputEnergies(const Eigen::VectorXd& dftenergies,
                                 const Eigen::VectorXd& gwaenergies,
                                 Index qpmin) {
  Index rpatotal = _rpamax - _rpamin + 1;
  _energies = dftenergies.segment(_rpamin, rpatotal);
  Index gwsize = Index(gwaenergies.size());
  Index lumo = _homo + 1;

  Index qpmax = qpmin + gwsize - 1;
  _energies.segment(qpmin - _rpamin, gwsize) = gwaenergies;
  double DFTgap = dftenergies(lumo) - dftenergies(_homo);
  double QPgap = gwaenergies(lumo - qpmin) - gwaenergies(_homo - qpmin);
  double shift = QPgap - DFTgap;
  Index levelaboveqpmax = _rpamax - qpmax;
  _energies.segment(qpmax + 1 - _rpamin, levelaboveqpmax).array() += shift;
}

template <bool imag>
Eigen::MatrixXd RPA::calculate_epsilon(double frequency) const {
  const Index size = _Mmn.auxsize();
  std::vector<Eigen::MatrixXd> thread_result = std::vector<Eigen::MatrixXd>(
      OPENMP::getMaxThreads(), Eigen::MatrixXd::Zero(size, size));
  const Index lumo = _homo + 1;
  const Index n_occ = lumo - _rpamin;
  const Index n_unocc = _rpamax - lumo + 1;
  const double freq2 = frequency * frequency;
  const double eta2 = _eta * _eta;
#pragma omp parallel for
  for (Index m_level = 0; m_level < n_occ; m_level++) {
    const double qp_energy_m = _energies(m_level);

    const Eigen::MatrixXd Mmn_RPA = _Mmn[m_level].bottomRows(n_unocc);

    const Eigen::ArrayXd deltaE = _energies.tail(n_unocc).array() - qp_energy_m;
    Eigen::VectorXd denom;
    if (imag) {
      denom = 4 * deltaE / (deltaE.square() + freq2);
    } else {
      Eigen::ArrayXd deltEf = deltaE - frequency;
      Eigen::ArrayXd sum = deltEf / (deltEf.square() + eta2);
      deltEf = deltaE + frequency;
      sum += deltEf / (deltEf.square() + eta2);
      denom = 2 * sum;
    }
    thread_result[OPENMP::getThreadId()] +=
        Mmn_RPA.transpose() * denom.asDiagonal() * Mmn_RPA;
  }
  Eigen::MatrixXd result = Eigen::MatrixXd::Identity(size, size);
  for (const auto& mat : thread_result) {
    result += mat;
  }
  return result;
}

template Eigen::MatrixXd RPA::calculate_epsilon<true>(double frequency) const;
template Eigen::MatrixXd RPA::calculate_epsilon<false>(double frequency) const;

RPA::rpa_eigensolution RPA::calculate_eigenvalues() const {
  Eigen::VectorXd AmB = calculate_spectral_AmB();
  return diag_C(AmB, calculate_spectral_C(AmB, calculate_spectral_ApB()));
}

Eigen::VectorXd RPA::calculate_spectral_AmB() const {
  const int lumo = _homo + 1;
  const int n_occup = lumo - _rpamin;
  const int n_unocc = _rpamax - _homo;
  const int rpasize = n_occup * n_unocc;
  vc2index vc = vc2index(0, 0, n_unocc);
  Eigen::VectorXd AmB = Eigen::VectorXd::Zero(rpasize);

  for (int v = 0; v < n_occup; v++) {
    int i = vc.I(v, 0);
    AmB.segment(i, n_unocc) =
        _energies.segment(n_occup, n_unocc).array() - _energies(v);
  }

  return AmB;
}

Eigen::MatrixXd RPA::calculate_spectral_ApB() const {
  const int lumo = _homo + 1;
  const int n_occup = lumo - _rpamin;
  const int n_unocc = _rpamax - _homo;
  const int rpasize = n_occup * n_unocc;
  const int auxsize = _Mmn.auxsize();
  vc2index vc = vc2index(0, 0, n_unocc);
  Eigen::MatrixXd ApB = Eigen::MatrixXd::Zero(rpasize, rpasize);

#pragma omp parallel for schedule(guided)
  for (int v2 = 0; v2 < n_occup; v2++) {
    int i2 = vc.I(v2, 0);
    const Eigen::MatrixXd Mmn_v2T =
        _Mmn[v2].block(n_occup, 0, n_unocc, auxsize).transpose();
    for (int v1 = v2; v1 < n_occup; v1++) {
      int i1 = vc.I(v1, 0);
      // Multiply with factor 2 to sum over both (identical) spin states
      ApB.block(i1, i2, n_unocc, n_unocc) =
          2 * 2 * _Mmn[v1].block(n_occup, 0, n_unocc, auxsize) * Mmn_v2T;
    }
  }

  ApB.diagonal() += calculate_spectral_AmB();

  return ApB;
}

Eigen::MatrixXd RPA::calculate_spectral_C(const Eigen::VectorXd& AmB,
                                          const Eigen::MatrixXd& ApB) const {
  return AmB.cwiseSqrt().asDiagonal() * ApB * AmB.cwiseSqrt().asDiagonal();
}

RPA::rpa_eigensolution RPA::diag_C(const Eigen::VectorXd& AmB,
                                   const Eigen::MatrixXd& C) const {
  const int lumo = _homo + 1;
  const int n_occup = lumo - _rpamin;
  const int n_unocc = _rpamax - _homo;
  const int rpasize = n_occup * n_unocc;

  XTP_LOG_SAVE(logDEBUG, _log)
      << TimeStamp() << " Diagonalizing RPA Hamiltonian " << flush;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(C);  // Uses lower triangle
  XTP_LOG_SAVE(logDEBUG, _log)
      << TimeStamp() << " Diagonalization done " << flush;

  double mc = es.eigenvalues().minCoeff();
  if (mc <= 0.0) {
    std::string msg =
        (boost::format("Detected non-positive eigenvalue: %s") % mc).str();
    XTP_LOG_SAVE(logDEBUG, _log) << TimeStamp() << " " << msg << flush;
    throw std::runtime_error(msg);
  }
  rpa_eigensolution sol;
  // Omega has to have correct size otherwise MKL does not rescale for Sqrt
  sol._Omega = Eigen::VectorXd::Zero(rpasize);
  sol._Omega = es.eigenvalues().cwiseSqrt();

  XTP_LOG_SAVE(logDEBUG, _log)
      << TimeStamp() << " Lowest neutral excitation energy (eV): "
      << tools::conv::hrt2ev * sol._Omega.minCoeff() << flush;

  // X     = 0.5 * [Omega^(-1/2) * (A-B)^(+1/2) + Omega^(+1/2) * (A-B)^(-1/2)] *
  // Z Y     = 0.5 * [Omega^(-1/2) * (A-B)^(+1/2) - Omega^(+1/2) * (A-B)^(-1/2)]
  // * Z X + Y =       [Omega^(-1/2) * (A-B)^(+1/2) ] * Z
  sol._XpY = Eigen::MatrixXd(rpasize, rpasize);
  Eigen::VectorXd AmB_sqrt = AmB.cwiseSqrt();
  Eigen::VectorXd Omega_sqrt_inv = sol._Omega.cwiseSqrt().cwiseInverse();
  for (int s = 0; s < rpasize; s++) {
    sol._XpY.col(s) =
        Omega_sqrt_inv(s) * AmB_sqrt.cwiseProduct(es.eigenvectors().col(s));
  }

  return sol;
}

}  // namespace xtp
}  // namespace votca
