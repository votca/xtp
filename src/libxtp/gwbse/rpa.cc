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

#include "votca/xtp/threecenter.h"
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/rpa.h>
#include "votca/xtp/vc2index.h"

using std::flush;

namespace votca {
namespace xtp {

void RPA::UpdateRPAInputEnergies(const Eigen::VectorXd& dftenergies,
                                 const Eigen::VectorXd& gwaenergies,
                                 int qpmin) {
  int rpatotal = _rpamax - _rpamin + 1;
  _energies = dftenergies.segment(_rpamin, rpatotal);
  int gwsize = gwaenergies.size();
  int lumo = _homo + 1;

  int qpmax = qpmin + gwsize - 1;
  _energies.segment(qpmin - _rpamin, gwsize) = gwaenergies;
  double DFTgap = dftenergies(lumo) - dftenergies(_homo);
  double QPgap = gwaenergies(lumo - qpmin) - gwaenergies(_homo - qpmin);
  double shift = QPgap - DFTgap;
  int levelaboveqpmax = _rpamax - qpmax;
  _energies.segment(qpmax + 1 - _rpamin, levelaboveqpmax).array() += shift;
}

template <bool imag>
Eigen::MatrixXd RPA::calculate_epsilon(double frequency) const {
  const int size = _Mmn.auxsize();
  Eigen::MatrixXd result = Eigen::MatrixXd::Identity(size, size);
  const int lumo = _homo + 1;
  const int n_occ = lumo - _rpamin;
  const int n_unocc = _rpamax - lumo + 1;
  const double freq2 = frequency * frequency;
  const double eta2 = _eta * _eta;
#pragma omp parallel for
  for (int m_level = 0; m_level < n_occ; m_level++) {
    const double qp_energy_m = _energies(m_level);

#if (GWBSE_DOUBLE)
    const Eigen::MatrixXd Mmn_RPA =
        _Mmn[m_level].block(n_occ, 0, n_unocc, size);
#else
    const Eigen::MatrixXd Mmn_RPA =
        _Mmn[m_level].block(n_occ, 0, n_unocc, size).cast<double>();
#endif
    const Eigen::ArrayXd deltaE =
        _energies.segment(n_occ, n_unocc).array() - qp_energy_m;
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
    auto temp = Mmn_RPA.transpose() * denom.asDiagonal();
    Eigen::MatrixXd tempresult = temp * Mmn_RPA;

#pragma omp critical
    { result += tempresult; }
  }
  return result;
}

  template Eigen::MatrixXd RPA::calculate_epsilon<true>(double frequency) const;
  template Eigen::MatrixXd RPA::calculate_epsilon<false>(double frequency) const;

    rpa_eigensolution RPA::calculate_eigenvalues() const {
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
      } // Occupied MO v
      
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
      
      ApB.diagonal() = calculate_spectral_AmB();
      
      // TODO: cache (A+B)? It's independent of the energy
      for (int v2 = 0; v2 < n_occup; v2++) {
        int i2 = vc.I(v2, 0);
        const Eigen::MatrixXd Mmn_v2T =
                _Mmn[v2].block(n_occup, 0, n_unocc, auxsize).transpose();
        for (int v1 = 0; v1 <= v2; v1++) {
          int i1 = vc.I(v1, 0);
          ApB.block(i2, i1, n_unocc, n_unocc) -=
                  2 * _Mmn[v1].block(n_occup, 0, n_unocc, auxsize) * Mmn_v2T;
        } // Occupied MO v1
      } // Occupied MO v2
      
      return ApB;
    }

    Eigen::MatrixXd RPA::calculate_spectral_C(const Eigen::VectorXd& AmB, const Eigen::MatrixXd& ApB) const {
      return AmB.cwiseSqrt().asDiagonal() * ApB * AmB.cwiseSqrt().asDiagonal();
    }

    rpa_eigensolution RPA::diag_C(const Eigen::VectorXd& AmB, const Eigen::MatrixXd& C) const {
      const int lumo = _homo + 1;
      const int n_occup = lumo - _rpamin;
      const int n_unocc = _rpamax - _homo;
      const int rpasize = n_occup * n_unocc;
      rpa_eigensolution sol;
      
      CTP_LOG(ctp::logDEBUG, _log) << ctp::TimeStamp()
              << " Solving for RPA eigenvalues " << flush;
      
      // Note: Eigen's SelfAdjointEigenSolver only uses the lower triangular part of C
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(C);
      
      double minEigenvalue = es.eigenvalues().minCoeff();
      if (minEigenvalue <= 0.0) {
        throw std::runtime_error(
          (boost::format("Detected non-positive eigenvalue: %s") % minEigenvalue).str());
      }

      // Note: Omega has to have correct size otherwise MKL does not rescale for Sqrt
      sol._Omega = Eigen::VectorXd::Zero(rpasize);
      sol._Omega = es.eigenvalues().cwiseSqrt();
      sol._XpY = Eigen::MatrixXd(rpasize, rpasize);

      Eigen::VectorXd AmB_sqrt = AmB.cwiseSqrt();
      Eigen::VectorXd AmB_sqrt_inv = AmB_sqrt.cwiseInverse();
      Eigen::VectorXd Omega_sqrt = sol._Omega.cwiseSqrt();

      for (int s = 0; s < rpasize; s++) {
        Eigen::VectorXd lhs = (1 / Omega_sqrt(s)) * AmB_sqrt;
        Eigen::VectorXd rhs = (1 * Omega_sqrt(s)) * AmB_sqrt_inv;
        const Eigen::VectorXd& z = es.eigenvectors().col(s);
        sol._XpY.col(s) = 0.50 * ((lhs + rhs).cwiseProduct(z) + (lhs - rhs).cwiseProduct(z));
      }

      return sol;
    }

}  // namespace xtp
}  // namespace votca
