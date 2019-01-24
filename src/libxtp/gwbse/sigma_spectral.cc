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
            result(m) += Equation48(m, m, omega, residues);
          } // Energy level m
        } // Eigenvalues/poles s

      } else {

        for (int s = 0; s < numeigenvalues; s++) {
          double omega = _EigenSol._Omega(s);
          Eigen::MatrixXd residues = CalcResidues(s);
          for (int m = 0; m < _qptotal; m++) {
            result(m) += Equation47(m, m, RPAEnergies, RPAEnergies(m), omega, residues); // TODO: Pass frequency or energy?
          } // Energy level m
        } // Eigenvalues/poles s

      }
      
      return result;
    }

    Eigen::MatrixXd Sigma_Spectral::CalcCorrelationOffDiag(const Eigen::VectorXd& frequencies) const {
      const Eigen::VectorXd& RPAEnergies = _rpa.getRPAInputEnergies();
      const int numeigenvalues = _EigenSol._Omega.size();
      Eigen::MatrixXd result = Eigen::MatrixXd::Zero(_qptotal, _qptotal);
      
      if (_HedinApprox) {

        for (int s = 0; s < numeigenvalues; s++) {
          double omega = _EigenSol._Omega(s);
          Eigen::MatrixXd residues = CalcResidues(s);
          for (int m = 0; m < _qptotal; m++) {
            for (int n = m + 1; n < _qptotal; n++) {
              double res = Equation48(m, n, omega, residues);
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
              double result_m = Equation47(m, n, RPAEnergies, RPAEnergies(m), omega, residues); // TODO: Pass frequency or energy?
              double result_n = Equation47(m, n, RPAEnergies, RPAEnergies(n), omega, residues);
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
      const int lumo = _homo + 1;
      const int n_occup = lumo - _qpmin;
      const int n_unocc = _qpmax - _homo;
      vc2index vc = vc2index(_qpmin, lumo, n_unocc);
      const int auxsize = _Mmn.auxsize(); // Size of gwbasis
      Eigen::MatrixXd residues = Eigen::MatrixXd::Zero(_qptotal, _qptotal);
      Eigen::VectorXd xpy = _EigenSol._XpY.col(s);

      for (int m = 0; m < _qptotal; m++) {
        for (int n = 0; n < _qptotal; n++) {
          for (int v = _qpmin; v <= _homo; v++) {
            for (int c = lumo; c <= _qpmax; c++) {
              int i = vc.I(v, c); // Composite index i
              double fc = 0.0;
              for (int i_aux = 0; i_aux < auxsize; i_aux++) {
                fc += _Mmn[m].col(i_aux)[n] * _Mmn[v].col(i_aux)[c];
              } // Auxiliary basis function
              residues(m, n) += fc * xpy(i); // Eq. 45
            } // Unoccupied MO c
          } // Occupied MO v
        } // Energy level n
      } // Energy level m

      return residues;
    }

    // TODO: Name, input args
    double Sigma_Spectral::Equation47(int m, int n, const Eigen::VectorXd& energies, double w, double omega, Eigen::MatrixXd& residues) const {
      const double eta = 1e-6;
      
      double s1_Real = 0.0;
      double s1_Imag = 0.0;
      double s2_Real = 0.0;
      double s2_Imag = 0.0;

      // Eq. 47, part 1
      for (int v = 0; v <= _homo; v++) {
        double A = residues(m, v) * residues(n, v);
        double B = w - energies(v) + omega;
        double C_Real = A * B / (B * B + eta * eta);
        double C_Imag = A * eta / (B * B + eta * eta);
        s1_Real += C_Real;
        s1_Imag += C_Imag;
      } // Occupied MO v

      // Eq. 47, part 2
      for (int c = _homo + 1; c < _qpmax; c++) {
        double A = residues(m, c) * residues(n, c);
        double B = w - energies(c) - omega;
        double C_Real = A * B / (B * B + eta * eta);
        double C_Imag = A * eta / (B * B + eta * eta);
        s2_Real += C_Real;
        s2_Imag += C_Imag;
      } // Occupied MO c

      // TODO: Return complex result
      return s1_Real + s2_Real;
    }

    // TODO: Name, input args
    double Sigma_Spectral::Equation48(int m, int n, double omega, Eigen::MatrixXd& residues) const {
      
      double s1 = 0.0;
      double s2 = 0.0;

      // TODO: Use Eigen for summations
      for (int v = 0; v <= _homo; v++) {
        s1 += residues(m, v) * residues(n, v);
      } // Occupied energy levels v
      for (int k = 0; k < _qpmax; k++) {
        s2 += residues(m, k) * residues(n, k);
      } // All energy levels m

      return (2 * s1 - s2) / omega;
    }

  }
};
