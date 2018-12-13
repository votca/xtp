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

#include <complex>
#include <cmath>
#include <votca/xtp/eigen.h>
#include <votca/xtp/rpa_spectral.h>

namespace votca {
namespace xtp {

class TCMatrix_gwbse;

class Sigma_Spectral {
    
private:
    
    // BSE Options
    int _bse_size; // Total number of BSE energy levels

    // QP Options
    int _qp_homo; // Index of homo energy level
    int _qp_min; // Minimal energy level index
    int _qp_size; // Total number of QP energy levels <= _bse_size
    
    // Internal Options
    bool _Hedin; // Heddin's Static Approximation

    // Convergence criteria for g iteration (Hartree)
    double _g_sc_limit;
    int _g_sc_max_iterations;

    // Energies
    Eigen::VectorXd _gwa_energies;

    // Sigma matrix
    Eigen::MatrixXd _Sigma_x; // Exchange part of sigma
    Eigen::MatrixXcd _Sigma_c; // Correlation part of sigma

    void compute_sigma_x(const TCMatrix_gwbse& Mmn, double scaHFX);
    void compute_sigma_c(const TCMatrix_gwbse& Mmn, const RPA_Spectral& rpa);
    
    // Compute sigma_c (eq. 47, 48)
    void clear_sigma_c();
    void fill_sigma_c_diag(const TCMatrix_gwbse& Mmn, const RPA_Spectral& rpa);
    void fill_sigma_c_offdiag(const TCMatrix_gwbse& Mmn, const RPA_Spectral& rpa);
    std::complex<double> Equation47(int m, int n, double w, Eigen::MatrixXd& res, double omega);
    double Equation48(int m, int n, Eigen::MatrixXd& res, double omega);
    
public:

    Sigma_Spectral() {
        _gwa_energies.resize(0);
    }

    void configure_bse(int homo, int vmin, int cmax) {

        // TODO: Which of these values should we keep as member?
        int bse_homo = homo;
        int bse_vmin = vmin;
        int bse_vmax = homo;
        int bse_cmin = homo + 1;
        int bse_cmax = cmax;
        int bse_vtotal = bse_vmax - bse_vmin + 1;
        int bse_ctotal = bse_cmax - bse_cmin + 1;
        _bse_size = bse_vtotal * bse_ctotal;

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
    
    bool getHedin() {
        return _Hedin;
    }
    
    void setHedin(bool value) {
        _Hedin = value;
    }

    void configure_g_iter(int g_sc_max_iterations, double g_sc_limit) {
        
        _g_sc_limit = g_sc_limit;
        _g_sc_max_iterations = g_sc_max_iterations;
        
        if (_g_sc_max_iterations < 1) {
            _g_sc_max_iterations = 1;
        }
        
        return;
        
    }
    
    const Eigen::VectorXd& get_GWAEnergies() const {
        return _gwa_energies;
    }
    
    void set_GWAEnergies(Eigen::VectorXd& gwa_energies) {
        _gwa_energies = gwa_energies; // Creates a copy
    }

    const Eigen::MatrixXd& get_sigma_x() const {
        return _Sigma_x;
    }
    
    const Eigen::MatrixXcd& get_sigma_c() const {
        return _Sigma_c;
    }

    void compute_sigma(TCMatrix_gwbse& Mmn, RPA_Spectral& rpa, double scaHFX);
    Eigen::MatrixXd SetupFullQPHamiltonian(Eigen::MatrixXd& vxc);
    //void refine_energies(TCMatrix_gwbse& Mmn, RPA_Spectral& rpa, double scaHFX, const Eigen::VectorXd& dft_energies, const Eigen::MatrixXd& vxc);

    void FreeMatrices() {
        
        // TODO

        return;
    }
    
};
}
}

#endif /* _VOTCA_XTP_SIGMA_SPECTRAL_H */
