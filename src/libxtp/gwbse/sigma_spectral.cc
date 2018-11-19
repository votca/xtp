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

#include <votca/xtp/sigma_spectral.h>
#include <votca/xtp/threecenter.h>

namespace votca {
    namespace xtp {

        void Sigma_Spectral::prepare_decomp(const TCMatrix_gwbse& Mmn) {
            
            Eigen::MatrixXd C = Eigen::MatrixXd::Zero(_bse_size, _bse_size); // Fill with (A + B)
            Eigen::VectorXd AmB_sqrt = Eigen::VectorXd::Zero(_bse_size); // Fill with (A - B)
            
            // Fill (A + B), (A - B)
            Fill_AB(Mmn, C, AmB_sqrt);
            // Compute (A - B)^(1 / 2)
            AmB_sqrt = AmB_sqrt.cwiseSqrt();
            // Compute C
            C = AmB_sqrt.asDiagonal() * C * AmB_sqrt.asDiagonal();
            // Diagonalize C
            Diag_C(AmB_sqrt, C);
            
            std::cout << _Omega << std::endl;

            return;
        }
        
        void Sigma_Spectral::Fill_AB(const TCMatrix_gwbse& Mmn, Eigen::MatrixXd& ApB, Eigen::VectorXd& AmB) {

            // TODO: Make use of the fact that C is symmetric!

            for (int i_1 = 0; i_1 < _bse_size; i_1++) {

                // eps_{c_1} - eps_{v_1}
                double denergy = _gwa_energies(_index2c[i_1]) - _gwa_energies(_index2v[i_1]);

                ApB(i_1, i_1) = denergy; // Fill (A + B)
                AmB(i_1) = denergy; // Fill (A - B)

            } // Composite index 1

            for (int i_1 = 0; i_1 < _bse_size; i_1++) {

                int v_1 = _index2v[i_1];
                int c_1 = _index2c[i_1];

                for (int i_aux = 0; i_aux < Mmn.getAuxDimension(); ++i_aux) {

                    // Get three-center column for index 1
                    VectorXfd tc_1 = Mmn[v_1].col(i_aux);

                    // TODO: Can we use Eigen sums instead of this for-loop?
                    for (int i_2 = 0; i_2 < _bse_size; i_2++) {

                        int v_2 = _index2v[i_2];
                        int c_2 = _index2c[i_2];

                        // Get three-center column for index 2
                        VectorXfd tc_2 = Mmn[v_2].col(i_aux);
                        
                        // -2 * (v1c1|v2c2)
                        ApB(i_1, i_2) -= 2 * tc_1(c_1) * tc_2(c_2); // Fill (A + B)
                        
                    } // Composite index 2
                } // Auxiliary basis functions
            } // Composite index 1

            return;

        }

        void Sigma_Spectral::Diag_C(Eigen::VectorXd& AmB_sqrt, Eigen::MatrixXd& C) {

            // Solve eigenvalue problem (eq. 41)
            // Diagonalize C to find the eigenvalues Sigma.
            // Using SelfAdjointEigenSolver since C is symmetric!
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(C);
            
            // Eigenvalues
            _Omega = es.eigenvalues().cwiseSqrt();
            
            // Eigenvectors
            Eigen::MatrixXd ZS = es.eigenvectors();

            // Setup matrices XS and YS (eq. 42)
            
            _X = Eigen::MatrixXd(_bse_size, _bse_size);
            _Y = Eigen::MatrixXd(_bse_size, _bse_size);

            Eigen::VectorXd AmB_sqrt_inv = AmB_sqrt.cwiseInverse();
            Eigen::VectorXd Omega_sqrt = _Omega.cwiseSqrt();

            for (int s = 0; s < _bse_size; s++) {
                
                Eigen::VectorXd lhs = (1 / Omega_sqrt(s)) * AmB_sqrt;
                Eigen::VectorXd rhs = (1 * Omega_sqrt(s)) * AmB_sqrt_inv;

                _X.col(s) = (lhs + rhs).cwiseProduct(ZS.col(s)) / 2.0;
                _Y.col(s) = (lhs - rhs).cwiseProduct(ZS.col(s)) / 2.0;
                
            }

            // TODO: More efficient way to solve eq. 36 without sqrt is
            // described in section 6.2.

            return;

        }
        
        void Sigma_Spectral::compute_sigma(const TCMatrix_gwbse& Mmn, double freq) {

            const double eta = 1e-6;

            _Sigma = Eigen::MatrixXd::Zero(_qp_size, _qp_size);
            
            for (int s = 0; s < _bse_size; s++ ) {
                
                Eigen::MatrixXd res = Eigen::MatrixXd::Zero(_bse_size, _bse_size);

                Calculate_Residues(Mmn, s, res);

                for (int m = 0; m < _qp_size; m++) {
                
                    for (int n = 0; n < _qp_size; n++) {

                        for (int i = 0; i < _bse_size; i++ ) {

                            int v = _index2v[i];
                            int c = _index2c[i];
                            
                            // Eq. 47
                            //
                            // w_{m, v}^s * w_{n, v}^s      A
                            // ----------------------- = ------- = A * B
                            //   w - ε_v + Ω_s + iη      (1 / B)
                            //
                            // C = w - ε_v + Ω_s
                            //
                            // B = 1 / (C - iη) = (C + iη) / ((C + iη) * (C - iη))
                            //   = ...
                            //   = (C + iη) / (C² - η²)

                            double A1 = res(m, v) * res(n, v);
                            double A2 = res(m, c) * res(n, c);

                            double C1 = freq - _gwa_energies(v) + _Omega(s);
                            double C2 = freq - _gwa_energies(c) + _Omega(s);

                            double B1_Real = C1 / (C1 * C1 - eta * eta);
                            double B2_Real = C2 / (C2 * C2 - eta * eta);

                            double B1_Imag = eta / (C1 * C1 - eta * eta);
                            double B2_Imag = eta / (C2 * C2 - eta * eta);

                            _Sigma(m, n) += A1 * B1_Real + A2 * B2_Real; // Eq. 47

                        } // Composite index i
                    } // Energy level n
                } // Energy level m
            } // Eigenvalues s

            return;

        }
        
        void Sigma_Spectral::Calculate_Residues(const TCMatrix_gwbse& Mmn, int s, Eigen::MatrixXd& res) {

            Eigen::VectorXd x = _X.col(s);
            Eigen::VectorXd y = _Y.col(s);

            for (int m = 0; m < _qp_size; m++) {

                for (int n = 0; n < _qp_size; n++) {

                    for (int i_aux = 0; i_aux < Mmn.getAuxDimension(); ++i_aux) {

                        // Get three-center column for index (m, n)
                        VectorXfd tc_mn = Mmn[m].col(i_aux);

                        for (int i = 0; i < _bse_size; i++) {

                            int v = _index2v[i];
                            int c = _index2v[i];

                            // Get three-center column for index (v, c)
                            VectorXfd tc_vc = Mmn[v].col(i_aux);

                            // Fill residues vector
                            res(m, n) += (tc_mn(n) * tc_vc(c)) * (x(i) + y(i)); // Eq. 45

                        } // Composite index i
                    } // i_aux
                } // Energy level n
            } // Energy level m
            
            return;
            
        }

    }
};
