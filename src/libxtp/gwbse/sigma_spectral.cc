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
          const Eigen::MatrixXd residues = CalcResidues(s);
          for (int m = 0; m < _qptotal; m++) {
            Eigen::VectorXd rm_x_rm = residues.row(m).cwiseAbs2();
            result(m) += Equation48(rm_x_rm, omega);
          } // Energy level m
        } // Eigenvalues/poles s

      } else {

        for (int s = 0; s < numeigenvalues; s++) {
          double omega = _EigenSol._Omega(s);
          const Eigen::MatrixXd residues = CalcResidues(s);
          for (int m = 0; m < _qptotal; m++) {
            Eigen::VectorXd rm_x_rm = residues.row(m).cwiseAbs2();
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
          const Eigen::MatrixXd residues = CalcResidues(s);
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
          const Eigen::MatrixXd residues = CalcResidues(s);
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
      const int n_unocc = _opt.qpmax - _opt.homo;
      vc2index vc = vc2index(_opt.qpmin, lumo, n_unocc);
      const int auxsize = _Mmn.auxsize(); // Size of gwbasis
      const Eigen::VectorXd xpy = _EigenSol._XpY.col(s);
      Eigen::MatrixXd residues = Eigen::MatrixXd::Zero(_qptotal, _qptotal);
      
      for (int v = _opt.qpmin; v <= _opt.homo; v++) {
        const Eigen::MatrixXd temp = _Mmn[v - _opt.rpamin].block(lumo - _opt.rpamin, 0, n_unocc, auxsize).transpose();
        const Eigen::VectorXd xpyv = xpy.segment(vc.I(v, lumo), n_unocc);
        for (int m = 0; m < _qptotal; m++) {
          // Sum over aux. basis functions chi
          const Eigen::MatrixXd fc = _Mmn[m].block(0, 0, _qptotal, auxsize) * temp;
          // Sum over unoccupied MOs c
          residues.row(m) += (fc * xpyv).transpose();
        } // MO m
      } // Occupied MO v

      return residues;
    }

    double Sigma_Spectral::Equation47(const Eigen::VectorXd& rm_x_rn, double omega, double w) const {
      const double eta = 1e-6;
      const int lumo = _opt.homo + 1;
      const int n_occup = lumo - _opt.qpmin;
      const int n_unocc = _opt.qpmax - _opt.homo;
      
      Eigen::VectorXd signs(_qptotal);
      signs << Eigen::VectorXd::Ones(n_occup), Eigen::VectorXd::Ones(n_unocc) * -1.0;
      
      const Eigen::VectorXd ones = Eigen::VectorXd::Ones(_qptotal);
      const Eigen::VectorXd temp = _rpa.getRPAInputEnergies().segment(_opt.qpmin - _opt.rpamin, _qptotal) * -1.0 + signs * omega + ones * w;
      const Eigen::VectorXd numer = rm_x_rn.cwiseProduct(temp);
      const Eigen::VectorXd denom = temp.cwiseAbs2() + ones * eta * eta;
      
      return numer.cwiseQuotient(denom).sum();
    }

    double Sigma_Spectral::Equation48(const Eigen::VectorXd& rm_x_rn, double omega) const {
      const int lumo = _opt.homo + 1;
      const int n_occup = lumo - _opt.qpmin;
      
      double s1 = rm_x_rn.head(n_occup).sum();
      double s2 = rm_x_rn.sum();

      return 2 * (s1 - s2) / omega;
    }

  }
};
