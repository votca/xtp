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

#ifndef _VOTCA_XTP_SIGMA_SPECTRAL_H
#define _VOTCA_XTP_SIGMA_SPECTRAL_H

#include <votca/xtp/eigen.h>
#include <votca/xtp/rpa_spectral.h>

namespace votca {
namespace xtp {

class TCMatrix_gwbse;

class Sigma_Spectral {
    
private:
    
    // BSE Options
    int _bse_size; // Total number of BSE energy levels

    // Composite Indexing
    std::vector<int> _index2v;
    std::vector<int> _index2c;

    // QP Options
    int _qp_homo; // Index of homo energy level
    int _qp_min; // Minimal energy level index
    int _qp_size; // Total number of QP energy levels <= _bse_size

    // Convergence criteria for g iteration (Hartree)
    double _g_sc_limit;
    int _g_sc_max_iterations;

    // Energies
    Eigen::VectorXd _gwa_energies;

    // Sigma matrix
    Eigen::MatrixXd _Sigma_x; // Exchange part of sigma
    Eigen::MatrixXd _Sigma_c; // Correlation part of sigma

    // Compute sigma (eq. 47)
    void compute_sigma_x(double freq, const TCMatrix_gwbse& Mmn, double scaHFX);
    void compute_sigma_c(double freq, const TCMatrix_gwbse& Mmn, const RPA_Spectral& rpa);
    
public:

    Sigma_Spectral() {
        _gwa_energies.resize(0);
    }

    void configure_bse(int homo, int vmin, int cmax, int nmax) {

        // TODO: Which of these values should we keep as member?
        int bse_homo = homo;
        int bse_vmin = vmin;
        int bse_vmax = homo;
        int bse_cmin = homo + 1;
        int bse_cmax = cmax;
        int bse_nmax = nmax;
        int bse_vtotal = bse_vmax - bse_vmin + 1;
        int bse_ctotal = bse_cmax - bse_cmin + 1;
        _bse_size = bse_vtotal * bse_ctotal;

        for (int v = 0; v < bse_vtotal; v++) {

            for (int c = 0; c < bse_ctotal; c++) {

                _index2v.push_back(bse_vmin + v);
                _index2c.push_back(bse_cmin + c);
            } // Unoccupied (/coulomb) energy level index c
        } // Occupied (/valance) energy level index v
        
        return;

    }

    void configure_qp(int homo, int min, int max) {

        // TODO: Which of these values should we keep as member?
        _qp_homo = homo;
        _qp_min = min;
        int qp_max = max;
        _qp_size = qp_max - _qp_min + 1;
        
        return;
        
    }

    void configure_g_iter(int g_sc_max_iterations, double g_sc_limit) {
        
        _g_sc_limit = g_sc_limit;
        _g_sc_max_iterations = g_sc_max_iterations;
        
        if (_g_sc_max_iterations < 1) {
            _g_sc_max_iterations = 1;
        }
        
        return;
        
    }
    
    const Eigen::VectorXd& getGWAEnergies() const {
        return _gwa_energies;
    }
    
    void setGWAEnergies(Eigen::VectorXd& gwa_energies) {
        _gwa_energies = gwa_energies; // Creates a copy
    }

    void refine_energies(double freq, TCMatrix_gwbse& Mmn, RPA_Spectral& rpa, double scaHFX, const Eigen::VectorXd& dft_energies, const Eigen::MatrixXd& vxc);

    void FreeMatrices() {
        
        // TODO

        return;
    }
    
};
}
}

#endif /* _VOTCA_XTP_SIGMA_SPECTRAL_H */
