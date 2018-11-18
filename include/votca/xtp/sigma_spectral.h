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

#ifndef _VOTCA_XTP_RPA_SPECTRAL_H
#define _VOTCA_XTP_RPA_SPECTRAL_H

#include <votca/xtp/eigen.h>
#include <votca/xtp/orbitals.h>

namespace votca {
namespace xtp {

class TCMatrix_gwbse;

class Sigma_Spectral {
    
private:
    
    int _bse_size; // Total number of energy levels

    // Composite indexing
    std::vector<int> _index2v;
    std::vector<int> _index2c;
    
    // Spectral decomposition
    Eigen::VectorXd _Omega; // Eigenvalues
    Eigen::MatrixXd _X, _Y; // Eigenvector components
    
    // Compute (A + B) and (A - B)
    void Fill_AB(const Orbitals& orbitals, const TCMatrix_gwbse& Mmn, Eigen::MatrixXd& ApB, Eigen::VectorXd& AmB);
    
    // Do the spectral decomposition
    void Diag_C(Eigen::VectorXd& AmB_sqrt, Eigen::MatrixXd& C);
    
public:

    void configure(int homo, int vmin, int cmax, int nmax) {

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

    // Prepare the spectral decomposition
    void prepare_decomp(const Orbitals& orbitals, const TCMatrix_gwbse& Mmn);
    
    // TODO: Calc sigma diagonal
    
    // TODO: Calc sigma off diagonal

    void FreeMatrices() {
        
        // TODO

        return;
    }
    
};
}
}

#endif /* _VOTCA_XTP_RPA_SPECTRAL_H */
