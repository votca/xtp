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

#ifndef _VOTCA_XTP_GWBSE_EXACT_H
#define _VOTCA_XTP_GWBSE_EXACT_H

#include <votca/xtp/orbitals.h>
#include <votca/xtp/aobasis.h>
#include <votca/xtp/threecenter.h>
#include <unsupported/Eigen/MatrixFunctions>
#include <eigen3/Eigen/src/Core/PlainObjectBase.h>

namespace votca
{
namespace xtp
{

class GWBSE_Exact {
    
private:
    
    // ****** Previous results ******
    
    Orbitals &_orbitals;
    AOBasis &_gwbasis;
    TCMatrix_gwbse &_Mmn;
    
    // ****** BSE ******
    
    int _homo; // Index of homo energy level
    int _bse_vmin; // v: valance, i.e. occupied
    int _bse_vmax;
    int _bse_cmin; // c: coulomb, i.e. unoccupied
    int _bse_cmax;
    int _bse_nmax; // Number of eigenvectors
    int _bse_vtotal;
    int _bse_ctotal;
    int _bse_size; // = _bse_vtotal * _bse_ctotal
    
    std::vector<int> _index2v;
    std::vector<int> _index2c;
    
    // ****** Intermediate results ******

    Eigen::VectorXd _AmB_sqrt; // (A - B)⁽1 / 2)
    Eigen::MatrixXd _C;

    Eigen::VectorXd _Sigma; // Eigenvalues

    Eigen::MatrixXd _XS; // Eigenvectors
    Eigen::MatrixXd _YS; // Eigenvectors

    // ****** RPA ******
    
    /*
    // (*) Use _home from BSE
    int _rpamin;
    int _rpamax;
    int _rpa_total;
    */
    
    // ****** QP ******
    
    // (*) Use _home from BSE
    int _qpmin;
    int _qpmax;
    int _qp_total;
    
    // ****** Screening ******

    // Container for the epsilon matrix
    std::vector<Eigen::MatrixXd > _epsilon_r;
    std::vector<Eigen::MatrixXd > _epsilon_i;
    
    // ****** Final results ******

    Eigen::VectorXd _sigma_c;

public:
    
    // ****** Constructor ******
    
    GWBSE_Exact(Orbitals& orbitals, AOBasis& gwbasis, TCMatrix_gwbse& Mmn) :
        _orbitals(orbitals),
        _gwbasis(gwbasis),
        _Mmn(Mmn) {
        
    }
    
    // ****** Init ******
    
    void Initialize() {
        
    }
    
    // ****** BSE ******

    void Configure_BSE(int homo, int vmin, int cmax, int nmax) {

        _homo = homo;
        _bse_vmin = vmin;
        _bse_vmax = homo;
        _bse_cmin = homo + 1;
        _bse_cmax = cmax;
        _bse_nmax = nmax;
        _bse_vtotal = _bse_vmax - _bse_vmin + 1;
        _bse_ctotal = _bse_cmax - _bse_cmin + 1;
        _bse_size = _bse_vtotal * _bse_ctotal;

        for (int v = 0; v < _bse_vtotal; v++) {

            for (int c = 0; c < _bse_ctotal; c++) {

                _index2v.push_back(_bse_vmin + v);
                _index2c.push_back(_bse_cmin + c);
            }
        }

        return;
        
    }
    
    // ****** Intermediate results ******
    
    void Fill_C();
    
    void Diag_C();
    
    Eigen::MatrixXd& Get_C() {
        return _C;
    }
    
    // ****** RPA ******
    
    /*
    void Configure_RPA(int homo, int rpamin, int rpamax) {
        
        _homo = homo;
        _rpamin = rpamin;
        _rpamax = rpamax;
        _rpa_total = _rpamax - _rpamin + 1;
        
        return;
        
    }
    */
    
    // ****** QP ******

    void Configure_QP(int homo, int qpmin, int qpmax) {
        
        _homo = homo;
        _qpmin = qpmin;
        _qpmax = qpmax;
        _qp_total = _qpmax - _qpmin + 1;
        
        return;
        
    }
    
    void Calculate_Residues(int s, Eigen::MatrixXd& res);
    
    // ****** Final results ******
    
    void Calculate_Sigma_Diagonal(double freq);
    
    // ****** Clean-up ******

    void FreeMatrices() {
        
        for (Eigen::MatrixXd& matrix : _epsilon_r) {
            matrix.resize(0, 0);
        }
        
        for (Eigen::MatrixXd& matrix : _epsilon_i) {
            matrix.resize(0, 0);
        }
        
        return;
        
    }

};

}
}

#endif /* _VOTCA_XTP_GWBSE_EXACT_H */