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
        
        void Sigma_Spectral::updateEnergies(const TCMatrix_gwbse& Mmn, const RPA_Spectral& rpa, double freq) {
            
            compute_sigma_x(Mmn, freq);
            compute_sigma_c(Mmn, rpa, freq);
            
            std::cout << "Sigma_x:" << std::endl << _Sigma_x << std::endl;
            std::cout << "Sigma_c:" << std::endl << _Sigma_c << std::endl;
            
            return;
            
        }
        
        void Sigma_Spectral::compute_sigma_x(const TCMatrix_gwbse& Mmn, double freq) {

            // TODO: Make use of the fact that sigma_x is symmetric
            
            _Sigma_x = Eigen::MatrixXd::Zero(_qp_size, _qp_size);

            for (int m = 0; m < _qp_size; m++) {

                const MatrixXfd& Mmn_m = Mmn[ m + _qp_min ];

                for (int n = 0; n < _qp_size; n++) {

                    const MatrixXfd & Mmn_n = Mmn[ n + _qp_min ];

                    for (int i_aux = 0; i_aux < Mmn.getAuxDimension(); ++i_aux) {

                        // Loop over all occupied bands used in screening
                        for (int i_occ = 0; i_occ <= _qp_homo; i_occ++) {

                            _Sigma_x(m, n) -= Mmn_m(i_occ, i_aux) * Mmn_n(i_occ, i_aux);

                        } // Occupied bands
                    } // Auxiliary basis functions
                } // Energy level n
            } // Energy level m
            
            return;
            
        }
        
        void Sigma_Spectral::compute_sigma_c(const TCMatrix_gwbse& Mmn, const RPA_Spectral& rpa, double freq) {

            const double eta = 1e-6;

            _Sigma_c = Eigen::MatrixXd::Zero(_qp_size, _qp_size);
            
            for (int s = 0; s < _bse_size; s++ ) {
                
                Eigen::MatrixXd res = Eigen::MatrixXd::Zero(_bse_size, _bse_size);
                rpa.compute_residues(Mmn, s, res);

                for (int m = 0; m < _qp_size; m++) {
                
                    for (int n = 0; n < _qp_size; n++) {

                        for (int i = 0; i < _bse_size; i++ ) {

                            int v = _index2v[i];
                            int c = _index2c[i];
                            
                            // Eq. 47
                            //
                            // w_{m, v}^s * w_{n, v}^s      A
                            // ----------------------- = ------- = A * B
                            // w - e_v + e_s + i * etaη  (1 / B)
                            //
                            // C = w - e_v + e_s
                            //
                            // B = 1 / (C - i * eta) = (C + i * eta) / ((C + i * eta) * (C - i * eta))
                            //   = ...
                            //   = (C + i * eta) / (C² - eta²)

                            double A1 = res(m, v) * res(n, v);
                            double A2 = res(m, c) * res(n, c);

                            double C1 = freq - _gwa_energies(v) + rpa.getOmega()(s);
                            double C2 = freq - _gwa_energies(c) + rpa.getOmega()(s);

                            double B1_Real = C1 / (C1 * C1 - eta * eta);
                            double B2_Real = C2 / (C2 * C2 - eta * eta);

                            double B1_Imag = eta / (C1 * C1 - eta * eta);
                            double B2_Imag = eta / (C2 * C2 - eta * eta);

                            _Sigma_c(m, n) += A1 * B1_Real + A2 * B2_Real; // Eq. 47

                        } // Composite index i
                    } // Energy level n
                } // Energy level m
            } // Eigenvalues s

            return;

        }

    }
};
