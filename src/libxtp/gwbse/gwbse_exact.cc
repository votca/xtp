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

#include <votca/xtp/gwbse_exact.h>
#include <eigen3/Eigen/src/Core/PlainObjectBase.h>

#include "votca/xtp/threecenter.h"

namespace votca {
    namespace xtp {

        class GWBSE_Exact {

            void Initialize() {

                return;

            }

            void Fill_C() {

                // C = (A - B)^(1 / 2) * (A + B) * (A - B)^(1 / 2)

                // (*)
                // Multiply left with diagonal: row scaling (slow)
                // Multiply right with diagonal: column scaling (fast)

                // Other idea:
                // Loop v1
                //   Cache denergy
                //   Loop basis function
                //     Cache tc1
                //     Loop v2
                //       Cache tc2
                //       Loop c1
                //         Loop c2

                // TODO: Make use of the fact that C is symmetric!

                // First, we compute (A - B)^(1 / 2) * (A + B).
                // Next, we multiply right with (A - B)^(1 / 2) to get C.
                _C = Eigen::MatrixXd::Zero(_bse_size, _bse_size);

                // (A - B)^(1 / 2)
                _AmB_sqrt = Eigen::VectorXd::Zero(_bse_size);

                for (int i_1 = 0; i_1 < _bse_size; i_1++) {

                    int v_1 = _index2v[i_1];
                    int c_1 = _index2v[i_1];

                    // e_c1 - e_v1
                    double denergy = _orbitals.MOEnergies()(c_1) - _orbitals.MOEnergies()(v_1);
                    double denergy_sqrt = std::sqrt(denergy);

                    _C(i_1, i_1) = denergy_sqrt * denergy; // Fill (A - B)^(1 / 2) * (A + B)
                    _AmB_sqrt(i_1) = denergy_sqrt; // Fill (A - B)

                    // TODO
                    for (int i_aux = 0; i_aux < _Mmn.getAuxDimension(); ++i_aux) {

                        // Get three-center column for index 1
                        Eigen::VectorXd tc_1 = _Mmn[v_1].col(i_f);

                        for (int i_2 = 0; i_2 < _bse_size; i_2++) {

                            int v_2 = _index2v[i_2];
                            int c_2 = _index2v[i_2];

                            // Get three-center column for index 2
                            Eigen::VectorXd tc_2 = _Mmn[v_2].col(i_f);

                            // -2 * (v1c1|v2c2)
                            _C(i_1, i_2) -= denergy_sqrt * 2 * tc_1(c_1) * tc_2(c_2);

                        } // Composite index 2

                    } // Auxiliary basis functions

                } // Composite index 1

                // Now, (A - B)^(1 / 2) * (A + B) is filled.
                // Next, we multiply right with (A - B)^(1 / 2) to get C.

                _C *= _AmB_sqrt.asDiagonal();

                // Now, C is filled.

                return;

            }

            void Diag_C() {

                // Diagonalize C to find the eigenvalues Sigma.
                // For this, use Eigen's EigenSolver:

                // https://eigen.tuxfamily.org/dox/classEigen_1_1EigenSolver.html

                // TODO: Use SelfAdjointEigenSolver
                // Because C is symmetric!

                // Solve eigenvalue problem (eq. 41)
                
                Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(_C);
                _Sigma = es.eigenvalues().sqrt();
                Eigen::MatrixXd ZS = es.eigenvectors();
                
                // Setup matrices XS and YS (eq. 42)

                Eigen::DiagonalMatrix<double> lhs, rhs;

                lhs = _Sigma.sqrt().pow(-1).asDiagonal() * _AmB_sqrt.asDiagonal();
                rhs = _Sigma.sqrt().asDiagonal() * _AmB_sqrt.pow(-1).asDiagonal();
                
                _XS = (1 / 2) * (lhs + rhs) * ZS;
                _YS = (1 / 2) * (lhs - rhs) * ZS;
                
                // TODO: More efficient way to solve eq. 36 without sqrt is
                // described in section 6.2.

            }

        };

    }
}