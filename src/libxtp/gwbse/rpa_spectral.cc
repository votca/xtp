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

#include <votca/xtp/rpa_spectral.h>
#include <votca/xtp/threecenter.h>

namespace votca {
    namespace xtp {

        void RPA_Spectral::prepare_decomp(const TCMatrix_gwbse& Mmn) {
            
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
        
        void RPA_Spectral::Fill_AB(const TCMatrix_gwbse& Mmn, Eigen::MatrixXd& ApB, Eigen::VectorXd& AmB) {

            // TODO: Make use of the fact that C is symmetric!

            for (int i = 0; i < _bse_size; i++) {

                double denergy = _gwa_energies(_index2c[i]) - _gwa_energies(_index2v[i]);

                ApB(i, i) = denergy; // Fill (A + B)
                AmB(i) = denergy; // Fill (A - B)

            } // Composite index i

            for (int i_1 = 0; i_1 < _bse_size; i_1++) {

                int v_1 = _index2v[i_1];
                int c_1 = _index2c[i_1];
                
                for (int i_2 = 0; i_2 < _bse_size; i_2++) {
                    
                    int v_2 = _index2v[i_2];
                    int c_2 = _index2c[i_2];
                    
                    double fourcenter = 0.0;
                    
                    for (int i_aux = 0; i_aux < Mmn.getAuxDimension(); ++i_aux) {
                        
                        VectorXfd tc_vc_1 = Mmn[v_1].col(i_aux);
                        VectorXfd tc_vc_2 = Mmn[v_2].col(i_aux);
                        
                        fourcenter += tc_vc_1(c_1) * tc_vc_2(c_2);
                        
                    } // Auxiliary basis function

                    ApB(i_1, i_2) -= 2 * fourcenter; // Fill (A + B)
                    
                } // Composite index i_2
            } // Composite index i_1

            return;

        }

        void RPA_Spectral::Diag_C(Eigen::VectorXd& AmB_sqrt, Eigen::MatrixXd& C) {

            // Solve eigenvalue problem (eq. 41)
            // Diagonalize C to find the eigenvalues Sigma.
            // Using SelfAdjointEigenSolver since C is symmetric!
            // (*) In MolGW, they diagonalize the block matrix (A, B; -B, -A) directly!
            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(C);
            
            // Eigenvalues
            // (*) In MolGW (.../src/linear_response.f90, line 279), they take
            // the absolute values of the "neutral excitation energies"
            _Omega = es.eigenvalues().cwiseAbs().cwiseSqrt();
            
            // Eigenvectors
            Eigen::MatrixXd Z = es.eigenvectors();

            // Setup matrix (X + Y) (eq. 42)
            _XpY = Eigen::MatrixXd(_bse_size, _bse_size);

            Eigen::VectorXd AmB_sqrt_inv = AmB_sqrt.cwiseInverse();
            Eigen::VectorXd Omega_sqrt = _Omega.cwiseSqrt();

            for (int s = 0; s < _bse_size; s++) {

                Eigen::VectorXd lhs = (1 / Omega_sqrt(s)) * AmB_sqrt;
                Eigen::VectorXd rhs = (1 * Omega_sqrt(s)) * AmB_sqrt_inv;

                _XpY.col(s) = 0.50 * (
                        (lhs + rhs).cwiseProduct(Z.col(s)) + // X
                        (lhs - rhs).cwiseProduct(Z.col(s))); // Y

            }

            // TODO: More efficient way to solve eq. 36 without sqrt is
            // described in section 6.2.

            return;

        }

        void RPA_Spectral::compute_residues(const TCMatrix_gwbse& Mmn, int s, Eigen::MatrixXd& res) const {

            Eigen::VectorXd xpy = _XpY.col(s);

            for (int m = 0; m < _qp_size; m++) {

                for (int n = 0; n < _qp_size; n++) {

                    for (int i_aux = 0; i_aux < Mmn.getAuxDimension(); i_aux++) {

                        // Get three-center column for index (m, n)
                        VectorXfd tc_mn = Mmn[m].col(i_aux);

                        for (int i = 0; i < _bse_size; i++) {

                            int v = _index2v[i];
                            int c = _index2c[i];

                            // Get three-center column for index (v, c)
                            VectorXfd tc_vc = Mmn[v].col(i_aux);

                            // Fill residues vector
                            res(m, n) += tc_mn(n) * tc_vc(c) * xpy(i); // Eq. 45

                        } // Composite index i
                    } // Auxiliary basis function
                } // Energy level n
            } // Energy level m

            return;

        }

    }
};
