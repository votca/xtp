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

#include <votca/xtp/sigma_spectral.h>
#include <votca/xtp/rpa.h>
#include <votca/xtp/threecenter.h>
#include "votca/xtp/vc2index.h"

namespace votca {
  namespace xtp {

    void Sigma_Spectral::PrepareScreening() {
      _HedinApprox = false;
      _EigenSol = _rpa.calculate_eigenvalues();
      return;
    }

    Eigen::VectorXd Sigma_Spectral::CalcCorrelationDiag(const Eigen::VectorXd& frequencies) const {
      const Eigen::VectorXd& RPAEnergies = _rpa.getRPAInputEnergies();
      const int numeigenvalues = _EigenSol._Omega.size();
      Eigen::VectorXd result = Eigen::VectorXd::Zero(_qptotal);

      if (_HedinApprox) {
        
        for (int s = 0; s < numeigenvalues; s++) {
          double omega = _EigenSol._Omega(s);
          Eigen::MatrixXd residues = CalcResidues(s);
          for (int m = 0; m < _qptotal; m++) {
            Eigen::VectorXd rm_x_rm = residues.row(m).cwiseProduct(residues.row(m));
            result(m) += Equation48(rm_x_rm, omega);
          } // Energy level m
        } // Eigenvalues/poles s

      } else {

        for (int s = 0; s < numeigenvalues; s++) {
          double omega = _EigenSol._Omega(s);
          Eigen::MatrixXd residues = CalcResidues(s);
          for (int m = 0; m < _qptotal; m++) {
            Eigen::VectorXd rm_x_rm = residues.row(m).cwiseProduct(residues.row(m));
            result(m) += Equation47(rm_x_rm, omega, RPAEnergies(m)); // TODO: Pass frequency or energy?
          } // Energy level m
        } // Eigenvalues/poles s

      }
      
      return result;
    }

    Eigen::MatrixXd Sigma_Spectral::CalcCorrelationOffDiag(const Eigen::VectorXd& frequencies) const {
      const Eigen::VectorXd& energies = _rpa.getRPAInputEnergies();
      const int numeigenvalues = _EigenSol._Omega.size();
      Eigen::MatrixXd result = Eigen::MatrixXd::Zero(_qptotal, _qptotal);
      
      if (_HedinApprox) {

        for (int s = 0; s < numeigenvalues; s++) {
          double omega = _EigenSol._Omega(s);
          Eigen::MatrixXd residues = CalcResidues(s);
          for (int m = 0; m < _qptotal; m++) {
            for (int n = m + 1; n < _qptotal; n++) {
              Eigen::VectorXd rm_x_rn = residues.row(m).cwiseProduct(residues.row(n));
              double res = Equation48(rm_x_rn, omega);
              result(m, n) += res;
              result(n, m) += res;
            } // Energy level n
          } // Energy level m
        } // Eigenvalues/poles s

      } else {

        for (int s = 0; s < numeigenvalues; s++) {
          double omega = _EigenSol._Omega(s);
          Eigen::MatrixXd residues = CalcResidues(s);
          for (int m = 0; m < _qptotal; m++) {
            for (int n = m + 1; n < _qptotal; n++) {
              Eigen::VectorXd rm_x_rn = residues.row(m).cwiseProduct(residues.row(n));
              double result_m = Equation47(rm_x_rn, omega, energies(m)); // TODO: Pass frequency or energy?
              double result_n = Equation47(rm_x_rn, omega, energies(n));
              // (m|S(w)|n) = 0.5 * (m|S(e_m)|n) + 0.5 * (m|S(e_n)|n)
              double res = 0.5 * (result_m + result_n);
              result(m, n) += res;
              result(n, m) += res;
            } // Energy level n
          } // Energy level m
        } // Eigenvalues/poles s

      }

      return result;
    }

    Eigen::MatrixXd Sigma_Spectral::CalcResidues(int s) const {
      const int lumo = _opt.homo + 1;
      const int n_occup = lumo - _opt.qpmin;
      const int n_unocc = _opt.qpmax - _opt.homo;
      vc2index vc = vc2index(_opt.qpmin, lumo, n_unocc);
      const int auxsize = _Mmn.auxsize(); // Size of gwbasis
      Eigen::MatrixXd residues = Eigen::MatrixXd::Zero(_qptotal, _qptotal);
      Eigen::VectorXd xpy = _EigenSol._XpY.col(s);
      
      // TODO: Loop over m first, if we want to paralellize
      for (int v = _opt.qpmin; v <= _opt.homo; v++) {
        const Eigen::MatrixXd temp = _Mmn[v - _opt.qpmin].block(lumo - _opt.qpmin, 0, n_unocc, auxsize).transpose();
        const Eigen::VectorXd xpyv = xpy.segment(vc.I(v, lumo), n_unocc);
        for (int m = 0; m < _qptotal; m++) {
          const Eigen::MatrixXd fc = _Mmn[m].block(0, 0, _qptotal, auxsize) * temp; // Sum over aux. basis functions chi
          residues.row(m) += (fc * xpyv).transpose(); // Sum over unoccupied MOs c
        } // MO m
      } // Occupied MO v

      return residues;
    }

    // TODO: Name
    double Sigma_Spectral::Equation47(const Eigen::VectorXd& rm_x_rn, double omega, double w) const {
      const double eta = 1e-6;
      const Eigen::VectorXd& energies = _rpa.getRPAInputEnergies();
      
      // TODO: Combine loops, use Eigen for summations?
      double s1_Real = 0.0;
      double s1_Imag = 0.0;
      double s2_Real = 0.0;
      double s2_Imag = 0.0;

      // Eq. 47, part 1
      for (int v = 0; v <= _opt.homo; v++) {
        double A = rm_x_rn[v];
        double B = w - energies(v) + omega;
        double C_Real = A * B / (B * B + eta * eta);
        double C_Imag = A * eta / (B * B + eta * eta);
        s1_Real += C_Real;
        s1_Imag += C_Imag;
      } // Occupied MO v

      // Eq. 47, part 2
      for (int c = _opt.homo + 1; c < _opt.qpmax; c++) {
        double A = rm_x_rn[c];
        double B = w - energies(c) - omega;
        double C_Real = A * B / (B * B + eta * eta);
        double C_Imag = A * eta / (B * B + eta * eta);
        s2_Real += C_Real;
        s2_Imag += C_Imag;
      } // Unoccupied MO c

      // TODO: Return complex result
      return s1_Real + s2_Real;
    }

    // TODO: Name
    double Sigma_Spectral::Equation48(const Eigen::VectorXd& rm_x_rn, double omega) const {
      const int lumo = _opt.homo + 1;
      const int n_occup = lumo - _opt.qpmin;
      
      double s1 = rm_x_rn.head(n_occup).sum(); // Occupied MOs v
      double s2 = rm_x_rn.sum(); // All MOs k
      return (2 * s1 - s2) / omega;
    }

  }
};
