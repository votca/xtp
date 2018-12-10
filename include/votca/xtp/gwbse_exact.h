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

#include <votca/ctp/logger.h>
#include <votca/xtp/eigen.h>

namespace votca
{
namespace xtp
{

class GWBSE_Exact {
    
private:
    
    // ****** Misc. ******
    
    ctp::Logger* _pLog;
    
    // ****** Orbitals ******
    
    Orbitals& _orbitals;

    int _homo; // HOMO index
    
    // RPA
    int _rpamin;
    int _rpamax;
    
    // QP
    int _qpmin;
    int _qpmax;
    int _qptotal;

    // BSE
    int _bse_vmin;
    int _bse_vmax;
    int _bse_cmin;
    int _bse_cmax;
    int _bse_maxeigenvectors;
    
    // ****** Basis sets ******
    
    std::string _auxbasis_name;
    std::string _dftbasis_name;
    
public:

    GWBSE_Exact(Orbitals& orbitals)
        : _orbitals(orbitals) {
    };
    
    void set_Logger(ctp::Logger* pLog) {
        _pLog = pLog;
    }
    
    void set_homo(int homo) {
        _homo = homo;
    }
    
    void configure_RPA(int rpamin, int rpamax) {
        _rpamin = rpamin;
        _rpamax = rpamax;
    }
    
    void configure_QP(int qpmin, int qpmax, int qptotal) {
        _qpmin = qpmin;
        _qpmax = qpmax;
        _qptotal = qptotal;
    }
    
    void configure_BSE(int bse_vmin, int bse_vmax, int bse_cmin, int bse_cmax, int bse_maxeigenvectors ) {
        _bse_vmin = bse_vmin;
        _bse_vmax = bse_vmax;
        _bse_cmin = bse_cmin;
        _bse_cmax = bse_cmax;
        _bse_maxeigenvectors = bse_maxeigenvectors;
    }
    
    void set_basis(std::string auxbasis_name, std::string dftbasis_name) {
        _auxbasis_name = auxbasis_name;
        _dftbasis_name = dftbasis_name;
    }

    Eigen::MatrixXd CalculateVXC(const AOBasis& dftbasis);
    bool Evaluate();

};

}
}

#endif /* _VOTCA_XTP_GWBSE_EXACT_H */