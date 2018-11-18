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
#include <votca/xtp/orbitals.h>
#include <votca/xtp/threecenter.h>

namespace votca {
    namespace xtp {

        void Sigma_Spectral::prepare_decomp(const Orbitals& orbitals, const TCMatrix_gwbse& Mmn) {
            
            Eigen::MatrixXd C = Eigen::MatrixXd::Zero(_bse_size, _bse_size); // Fill with (A + B)
            Eigen::VectorXd AmB_sqrt = Eigen::VectorXd::Zero(_bse_size); // Fill with (A - B)
            
            // Fill (A + B), (A - B)
            Fill_AB(orbitals, Mmn, C, AmB_sqrt);
            // Compute (A - B)^(1 / 2)
            AmB_sqrt = AmB_sqrt.cwiseSqrt();
            // Compute C
            C = AmB_sqrt.asDiagonal() * C * AmB_sqrt.asDiagonal();
            // Diagonalize C
            Diag_C(AmB_sqrt, C);
            
            std::cout << _Omega << std::endl;

            return;
        }
        
        void Sigma_Spectral::Fill_AB(const Orbitals& orbitals, const TCMatrix_gwbse& Mmn, Eigen::MatrixXd& ApB, Eigen::VectorXd& AmB) {

            // TODO: Make use of the fact that C is symmetric!

            for (int i_1 = 0; i_1 < _bse_size; i_1++) {

                // eps_{c_1} - eps_{v_1}
                double denergy = orbitals.MOEnergies()(_index2c[i_1]) - orbitals.MOEnergies()(_index2v[i_1]);

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

    }
};
