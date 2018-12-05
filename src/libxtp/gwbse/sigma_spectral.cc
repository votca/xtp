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
#include <votca/xtp/rpa_spectral.h>

namespace votca {
    namespace xtp {
        
        void Sigma_Spectral::compute_sigma(TCMatrix_gwbse& Mmn, RPA_Spectral& rpa, double scaHFX) {
            
            compute_sigma_x(Mmn, scaHFX);
            compute_sigma_c(Mmn, rpa);
            
            return;
            
        }

        void Sigma_Spectral::refine_energies(TCMatrix_gwbse& Mmn, RPA_Spectral& rpa,
                double scaHFX, const Eigen::VectorXd& dft_energies, const Eigen::MatrixXd& vxc) { // TODO: Pass object containing DFT data instead?

            if (_gwa_energies.size() < 1) {
                throw std::runtime_error("Sigma gwa_energies not set!");
            }
            
            std::cout << "_gwa_energies:" << std::endl << _gwa_energies << std::endl;
            
            // TODO: Only compute diagonal of sigma_x
            compute_sigma_x(Mmn, scaHFX);

            // Previous energies
            Eigen::VectorXd gwa_energies_prev = _gwa_energies; // Creates a copy

            // TODO: Only diagonal elements except for in final iteration
            for (int g_iter = 0; g_iter < _g_sc_max_iterations; g_iter++) {
                
                // TODO: Only re-compute diagonal of Sigma_C
                compute_sigma_c(Mmn, rpa);
                
                // Update energies
                for (int m = 0; m < _qp_size; m++) {
                    _gwa_energies(m + _qp_min) = dft_energies(m + _qp_min) + _Sigma_c(m, m) + _Sigma_x(m, m) - vxc(m, m);
                }

                Eigen::VectorXd diff = gwa_energies_prev - _gwa_energies;

                int index = 0;
                double diff_max = diff.cwiseAbs().maxCoeff(&index);
                double qp_gap = _gwa_energies(_qp_homo + 1) - _gwa_energies(_qp_homo);
                double dft_gap = dft_energies(_qp_homo + 1) - dft_energies(_qp_homo);
                
                std::cout
                        << "g_iter: " << g_iter + 1
                        << ", shift = " << qp_gap - dft_gap
                        << ", diff max = " << diff_max
                        << ", index = " << index
                        << std::endl;

                const double alpha = 0.0;
                _gwa_energies = (1 - alpha) * _gwa_energies + alpha * gwa_energies_prev;
                
                bool conv = (diff_max <= _g_sc_limit);

                if (conv) {
                    std::cout << "G self consistency cycle converged after " << g_iter + 1 << " iterations." << std::endl; break;
                } else if (g_iter == _g_sc_max_iterations - 1) {
                    std::cout << "G self consistency cycle not converged after " << _g_sc_max_iterations << " iterations." << std::endl; break;
                } else {
                    gwa_energies_prev = _gwa_energies; // Creates a copy
                }
                
            } // g_iter

            std::cout << "_gwa_energies:" << std::endl << _gwa_energies << std::endl;

            return;

        }
        
        void Sigma_Spectral::compute_sigma_x(const TCMatrix_gwbse& Mmn, double scaHFX) {

            // TODO: Make use of the fact that sigma_x is symmetric
            
            _Sigma_x = Eigen::MatrixXd::Zero(_qp_size, _qp_size);

            for (int m = 0; m < _qp_size; m++) {

                const MatrixXfd& Mmn_m = Mmn[ m + _qp_min ];

                for (int n = 0; n < _qp_size; n++) {

                    const MatrixXfd & Mmn_n = Mmn[ n + _qp_min ];
                    
                    double sigma_x = 0.0;

                    for (int i_aux = 0; i_aux < Mmn.getAuxDimension(); ++i_aux) {

                        // Loop over all occupied bands used in screening
                        for (int i_occ = 0; i_occ <= _qp_homo; i_occ++) {

                            sigma_x -= Mmn_m(i_occ, i_aux) * Mmn_n(i_occ, i_aux);

                        } // Occupied bands
                    } // Auxiliary basis functions
                    
                    _Sigma_x(m, n) = (1.0 - scaHFX) * sigma_x;
                    
                } // Energy level n
            } // Energy level m
            
            return;
            
        }
        
        void Sigma_Spectral::compute_sigma_c(const TCMatrix_gwbse& Mmn, const RPA_Spectral& rpa) {

            const double eta = 1e-6;

            _Sigma_c = Eigen::MatrixXd::Zero(_qp_size, _qp_size);
            
            if (_Hedin) {
                
                // TODO: Perhaps use Hedin's approximation for certain values of sigma?
                for (int s = 0; s < _bse_size; s++ ) {
                    
                    Eigen::MatrixXd res = Eigen::MatrixXd::Zero(_bse_size, _bse_size);
                    rpa.compute_residues(Mmn, s, res);
                    
                    double omega = rpa.getOmega()(s);
                    
                    for (int m = 0; m < _qp_size; m++) {
                        
                        for (int n = 0; n < _qp_size; n++) {

                            double s1 = 0.0;
                            double s2 = 0.0;

                            // TODO: Use Eigen for summations
                            for (int v = 0; v <= _qp_homo; v++) { s1 += res(m, v) * res(n, v); } // Occupied energy levels v
                            for (int k = 0; k <  _qp_size; k++) { s2 += res(m, k) * res(n, k); } // All energy levels m

                            _Sigma_c(m, n) += (2 * s1 - s2) / omega; // Eq. 48

                        } // Energy level n
                    } // Energy level m
                } // Eigenvalues/poles s
                
            } else {
                
                for (int s = 0; s < _bse_size; s++ ) {

                    Eigen::MatrixXd res = Eigen::MatrixXd::Zero(_bse_size, _bse_size);
                    rpa.compute_residues(Mmn, s, res);
                    
                    double omega = rpa.getOmega()(s);

                    for (int m = 0; m < _qp_size; m++) {
                        
                        for (int n = 0; n < _qp_size; n++) {
                            
                            // TODO: Use both e_m, e_n to correctly compute sigma:
                            // (m|S(w)|n) = 0.5 * (m|S(e_m)|n) + 0.5 * (m|S(e_n)|n)
                            double w = _gwa_energies(m); // Frequency/energy
                            
                            double s1 = 0.0;
                            double s2 = 0.0;
                            
                            for (int v = 0; v <= _qp_homo; v++) {
                                
                                double A = res(m, v) * res(n, v);
                                double B = w - _gwa_energies(v) + omega;
                                double C_Real = A * B / (B * B + eta * eta);
                                
                                s1 += C_Real; // Eq. 47, part 1
                                
                            } // Occupied energy levels v
                            
                            for (int c = _qp_homo + 1; c < _qp_size; c++) {
                                
                                double A = res(m, c) * res(n, c);
                                double B = w - _gwa_energies(c) - omega;
                                double C_Real = A * B / (B * B + eta * eta);
                                
                                s2 += C_Real; // Eq. 47, part 2
                                
                            } // Occupied energy levels c

                            _Sigma_c(m, n) += s1 + s2; // Eq. 47
                            
                        } // Energy level n
                    } // Energy level m
                } // Eigenvalues/poles s
                
            }

            return;

        }

    }
};
