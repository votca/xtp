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

#include <votca/ctp/logger.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/aobasis.h>
#include <votca/xtp/threecenter.h>
#include <votca/xtp/gwbse_exact.h>

namespace votca {
    namespace xtp {

        void GWBSE_Exact::Fill_C() {

            // TODO: Make use of the fact that C is symmetric!

            // Fill with (A + B)
            _C = Eigen::MatrixXd::Zero(_bse_size, _bse_size);
            
            // Fill with (A - B)
            _AmB_sqrt = Eigen::VectorXd::Zero(_bse_size);

            for (int i_1 = 0; i_1 < _bse_size; i_1++) {

                int v_1 = _index2v[i_1];
                int c_1 = _index2c[i_1];
                
                // ε_{c_1} - ε_{v_1}
                double denergy = _orbitals.MOEnergies()(c_1) - _orbitals.MOEnergies()(v_1);

                for (int i_aux = 0; i_aux < _Mmn.getAuxDimension(); ++i_aux) {

                    // Get three-center column for index 1
                    VectorXfd tc_1 = _Mmn[v_1].col(i_aux);

                    for (int i_2 = 0; i_2 < _bse_size; i_2++) {

                        int v_2 = _index2v[i_2];
                        int c_2 = _index2c[i_2];
                        
                        // (ε_{c_1} - ε_{v_1}) * 𝛿_{v_1, v_2} * 𝛿_{c_1, c_2}
                        if (v_1 == v_2 && c_1 == c_2) {

                            _C(i_1, i_2) = denergy; // Fill (A + B)
                            _AmB_sqrt(i_1) = denergy; // Fill (A - B)
                        }
                        
                        // Get three-center column for index 2
                        VectorXfd tc_2 = _Mmn[v_2].col(i_aux);
                        
                        // -2 * (v1c1|v2c2)
                        _C(i_1, i_2) -= 2 * tc_1(c_1) * tc_2(c_2);
                        
                    } // Composite index 2
                    
                } // Auxiliary basis functions
            
            } // Composite index 1

            _AmB_sqrt = _AmB_sqrt.cwiseSqrt();
            
            _C = _AmB_sqrt.asDiagonal() * _C * _AmB_sqrt.asDiagonal();

            return;

        }

        void GWBSE_Exact::Diag_C() {

            // Diagonalize C to find the eigenvalues Sigma.
            // For this, use Eigen's EigenSolver:

            // https://eigen.tuxfamily.org/dox/classEigen_1_1EigenSolver.html

            // TODO: Use SelfAdjointEigenSolver
            // Because C is symmetric!

            // Solve eigenvalue problem (eq. 41)

            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(_C);
            
            // Eigenvalues
            _Omega = es.eigenvalues();
            _Omega(0) = _Omega(1); // Cheat: Get rid of negative eigenvalue
            _Omega = _Omega.cwiseSqrt();
            
            // Eigenvectors
            Eigen::MatrixXd ZS = es.eigenvectors();

            // Setup matrices XS and YS (eq. 42)
            
            _X = Eigen::MatrixXd(_bse_size, _bse_size);
            _Y = Eigen::MatrixXd(_bse_size, _bse_size);

            Eigen::VectorXd AmB_sqrt_inv = _AmB_sqrt.cwiseInverse();
            Eigen::VectorXd Omega_sqrt = _Omega.cwiseSqrt();

            for (int s = 0; s < _bse_size; s++) {
                
                Eigen::VectorXd lhs = (1 / Omega_sqrt(s)) * _AmB_sqrt;
                Eigen::VectorXd rhs = (1 * Omega_sqrt(s)) * AmB_sqrt_inv;

                _X.col(s) = (lhs + rhs).cwiseProduct(ZS.col(s)) / 2.0;
                _Y.col(s) = (lhs - rhs).cwiseProduct(ZS.col(s)) / 2.0;
                
            }

            // TODO: More efficient way to solve eq. 36 without sqrt is
            // described in section 6.2.

            return;

        }
        
        void GWBSE_Exact::Calculate_Residues(int s, Eigen::MatrixXd& res) {

            Eigen::VectorXd x = _X.col(s);
            Eigen::VectorXd y = _Y.col(s);

            for (int m = 0; m < _qp_total; m++) {

                for (int n = 0; n < _qp_total; n++) {

                    for (int i_aux = 0; i_aux < _Mmn.getAuxDimension(); ++i_aux) {

                        // Get three-center column for index (m, n)
                        VectorXfd tc_mn = _Mmn[m].col(i_aux);

                        for (int i = 0; i < _bse_size; i++) {

                            int v = _index2v[i];
                            int c = _index2v[i];

                            // Get three-center column for index (v, c)
                            VectorXfd tc_vc = _Mmn[v].col(i_aux);

                            // Fill residues vector
                            res(m, n) += (tc_mn(n) * tc_vc(c)) * (x(i) + y(i)); // Eq. 45

                        }
                        
                    }

                }
                
            }
            
            return;
            
        }

        void GWBSE_Exact::Fill_Sigma_Diagonal(double freq) {
            
            double eta = 1e-6;

            _Sigma_Diagonal = Eigen::VectorXd::Zero(_qp_total);

            for (int s = 0; s < _bse_size; s++ ) {

                Eigen::MatrixXd res = Eigen::MatrixXd::Zero(_bse_size, _bse_size);

                Calculate_Residues(s, res);

                for (int n = 0; n < _qp_total; n++) {
                    
                    for (int i = 0; i < _bse_size; i++ ) {

                        int v = _index2v[i];
                        int c = _index2c[i];

                        double numer1 = res(n, v) * res(n, v);
                        double denom1 = freq - _orbitals.MOEnergies()(v) + _Omega(s) - v * eta;
                        
                        double numer2 = res(n, c) * res(n, c);
                        double denom2 = freq - _orbitals.MOEnergies()(c) + _Omega(s) - v * eta;

                        _Sigma_Diagonal(n) += (numer1 / denom1) + (numer2 / denom2); // Eq. 47
                        
                    }
                    
                }
                
            }

            return;

        }
        
        void GWBSE_Exact::Fill_Sigma(double freq) {
            
            double eta = 1e-6;

            _Sigma = Eigen::MatrixXd::Zero(_qp_total, _qp_total);
            
            for (int s = 0; s < _bse_size; s++ ) {
                
                Eigen::MatrixXd res = Eigen::MatrixXd::Zero(_bse_size, _bse_size);

                Calculate_Residues(s, res);

                for (int m = 0; m < _qp_total; m++) {
                
                    for (int n = 0; n < _qp_total; n++) {

                        for (int i = 0; i < _bse_size; i++ ) {

                            int v = _index2v[i];
                            int c = _index2c[i];

                            double numer1 = res(m, v) * res(n, v);
                            double denom1 = freq - _orbitals.MOEnergies()(v) + _Omega(s) - v * eta;

                            double numer2 = res(m, c) * res(n, c);
                            double denom2 = freq - _orbitals.MOEnergies()(c) + _Omega(s) - v * eta;

                            _Sigma(m, n) += (numer1 / denom1) + (numer2 / denom2); // Eq. 47

                        }

                    }
                
                }
                
            }

            return;
            
        }

        Eigen::MatrixXd GWBSE_Exact::SetupFullQPHamiltonian() {

            // Construct full QP Hamiltonian
            Eigen::MatrixXd Hqp = _Sigma - (*_vxc);
            
            // Diagonal elements are given by _qp_energies
            for (int m = 0; m < Hqp.rows(); m++) {
                Hqp(m, m) = _gwa_energies(m + _qpmin);
            }
            
            return Hqp;
        }
        
        bool GWBSE_Exact::Evaluate() {

            std::cout << ctp::TimeStamp() << " Molecule Coordinates [A] " << std::endl;

            for (QMAtom* atom : _orbitals.QMAtoms()) {

                std::string output = (boost::format("  %1$s"
                        "   %2$+1.4f %3$+1.4f %4$+1.4f")
                        % (atom->getType())
                        % (atom->getPos().getX() * tools::conv::bohr2ang)
                        % (atom->getPos().getY() * tools::conv::bohr2ang)
                        % (atom->getPos().getZ() * tools::conv::bohr2ang)).str();
                
                std::cout << output << std::endl;
            }

            // Check which QC program was used for the DFT run
            // -> Implicit info about MO coefficient storage order
            std::string dft_package = _orbitals.getQMpackage();
            std::cout << ctp::TimeStamp() << " DFT data was created by " << dft_package << std::endl;

            // Store information in _orbitals for later use
            _orbitals.setRPAindices(_rpamin, _rpamax);
            _orbitals.setGWAindices(_qpmin, _qpmax);
            _orbitals.setBSEindices(_bse_vmin, _bse_cmax, _bse_maxeigenvectors);

            if (_do_full_BSE)
                _orbitals.setTDAApprox(false);
            else {
                _orbitals.setTDAApprox(true);
            }

            // Load DFT basis
            BasisSet dftbs;
            dftbs.LoadBasisSet(_dftbasis_name);
            _orbitals.setDFTbasis(_dftbasis_name);
            std::cout << ctp::TimeStamp() << " Loaded DFT Basis Set " << _dftbasis_name << std::flush;
            
            // Fill DFT AO basis
            
            // Load aux. basis
            
            // Fill Exchange-Correlation Potential (VXC) matrix
            
            // Fill overlap matrix
            
            // Fill Coulomb matrix
            
            // Prepare 3-center integral object
            
            // Initialize QP energies
            
            // Start iterations

            return false;
        }

    }
}